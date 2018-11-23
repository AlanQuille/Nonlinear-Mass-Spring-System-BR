#define _USE_MATH_DEFINES
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <chrono>
#include <iomanip>
#include "Simulation.cpp"
#include "Eigen/Dense"
#include "Eigen/QR"


using namespace std;
using namespace Eigen;


unsigned long long rdtsc()
{
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}


int main(int argc, char** argv)
{

   //Start the clock to measure how long the code takes to execute
   auto begin = std::chrono::high_resolution_clock::now();


    //data structure with the input values
    InitialDataValues data;
    //send from processor
    srand(rdtsc());

    //The target signal, Volterra change to target signal
    vector<double> Volterra;
    //The input signal
    vector<double> Input;
    
    
    //This takes the input in raw string from from the CSV file.
    vector<string> Input_Lines;
    
    //This takes the target signal in raw string form from the CSV file.
    vector<string> Volterra_Lines;  

    //Cue to start programme
    cout << "-- Start ---------------------------------------- " << endl;


    //open up ifstream to access csv file which contains input signal
    ifstream file_Input ( "Data/inputsignal.csv" ); file_Input.precision(15);
    string tmp;

    //This reads in each line frmo the input signal CSV file and it is added to the string vector
    while (getline(file_Input, tmp,'\n'))
    {
        Input_Lines.push_back(tmp); //Get each line of the file as a string
    }

     
    //Length of the vector containing the input signal in string form
    unsigned long s = Input_Lines.size();
    double x;
    
    //Convert the string to a double and add to the Input vector, which is the vecotr containg the input signal
    for (unsigned int i=1; i<s; ++i)
    {
        size_t pos = Input_Lines[i].find(",");      // position of the end of the name of each one in the respective string
        x = stod(Input_Lines[i].substr(pos+1,Input_Lines[i].size()));
        x = x*1;
        Input.push_back(x); // convert string age to a double
    }




    //open up ifstream to access csv file which contains input signal
     ifstream file_Volterra ("Data/volterra.csv");

    while (getline(file_Volterra, tmp,'\n'))
    {
        Volterra_Lines.push_back(tmp); //Get each line of the file as a string
    }

    //Length of the vector containing the target signal in string form
    s = Volterra_Lines.size();
    for (unsigned int i=1; i<s; ++i)
    {

        size_t pos = Volterra_Lines[i].find(",");      // position of the end of the name of each one in the respective string
        x = stod(Volterra_Lines[i].substr(pos+1,Volterra_Lines[i].size()));
        x = x*1;
        Volterra.push_back(x); // convert string age to a doubl
    }


    //The size of the input and target signal
    cout <<"Size of input signal is: "<< Input.size() << endl;
    cout <<"Size of target signal is: "<< Volterra.size() << endl;

    
	//The part of the target signal that is eliminated. This is determined via an impulse response, values before the impulse response settles down to a cycle or a steady state value are discarded
	double wash_out_time = 20000;
   //The part of the target signal that is used to determine the optimal weight vector
    double learning_time = 200000;
    //The part of the target signal that is used to test the optimal weight vector
    double learning_time_test = 15000;


    //These values are specifically for the random dynamical reservoir RNN (recurrent neural network), NOT THE SPIDER WEB
    //The number of nodes in a random RNN
    data.N = 25;
    //The horizontal input force on the input nodes
    data.ux=0;
    //The vertical input force on the input nodes
    data.uy= 0;

    //the proportion of nodes that take in the input signal
    data.input_connectivity = 0.05;

    //This is the minimum value for the rand_in_range function for the input weights
    data.min_input_weight = -0.0001 * 1;
    //This is the maximum value for the rand_in_range function for the input weights
    data.max_input_weight = 0.0001 * 1;

    //This is the minimum value for the rand_in_range functio for the input weights
    data.min_input_weight = -0.005 * 1;
    //This is the maximum value for the rand_in_range functio for the input weights
    data.max_input_weight = -0.005 * 1;

    //These are the minimum and maximum x andy value sof the nodes for the random dynamical reservoir RNN
    data.min_x_position = 0;
    data.max_x_position = 10;
    data.min_y_position  = 0;
    data.max_y_position = 10;


    //These are the minimum and maximum values for the spring and damping coefficients
    data.min_k3 = 50;
    data.max_k3  = 50;

    data.min_d3 = 50;
    data.max_d3  =50;

    data.min_k1 = 100;
    data.max_k1  = 100;

    data.min_d1 = 100;
    data.max_d1  = 100;
    
    //For nodes with identical masses, this is the value for the mass of each and every node
    data.mass_of_nodes = 1;
    //This is the scaling factor for the r
    data.scaling_factor = 0.1;


    //The time step in this case it is 0.001, 1 ms.
    data.dt = 0.001;
    //The initial time t0, the clock starts when the wash out time ends
    data.t0 = wash_out_time*data.dt;
    //The maximum time, this is the final time of the target signal 
    data.tmax = (wash_out_time+learning_time+learning_time_test)*data.dt;

    //Sine wave, simple example for testing.
    vector<double> Sine_Wave;

    //The first value of the input signal
    cout << "Initial Input is: "<< Input[0] << endl;
    
    
 //////////////////////////////////////////////////////////////////////////////////////////////////
 //THIS SIMULATES THE  RANDOM DYNAMICAL RESERVOIR RNN (recurrent neural network), NOT THE SPIDER WEB
 //////////////////////////////////////////////////////////////////////////////////////////////////
 Simulation sim1(data, Input, Volterra, wash_out_time, learning_time, learning_time_test);

//radius of the threads in the spider web
double radius = 10.0;
//number of spirals in the web
int rounds = 2;
//  int rounds = 24;
//Number of nodes (connection points between threads in the web) per spiral
int no_of_points_per_round= 5;

//This is to name the spider web
string str = "1";
  
vector<double> range_d1_d3_list;
double best_range_d1_d3;

double a = 100;


//Minimum value of spring coefficient for rand_in_range function
data.min_k1 = 1;
//Minimum value of spring coefficient for rand_in_range function
data.max_k1  = 1;

//Minimum value of spring coefficient for rand_in_range function
data.min_k3 = 1;
//Minimum value of spring coefficient for rand_in_range function
data.max_k3  = 1;

//Minimum value of spring coefficient for rand_in_range function
data.min_d1 = 1;
//Minimum value of spring coefficient for rand_in_range function
data.max_d1  = 1;

//Minimum value of spring coefficient for rand_in_range function
data.min_d3 = 1;
//Minimum value of spring coefficient for rand_in_range function
data.max_d3  = 1;
  

//This determines if the springs have identical spring and damping coefficients
    bool springs_identical = false;
//This determines if the bias column (all 1's) for linear regression is appended to teh learning matrix
    bool bias_learning = true;
//This determines if the input signal is a unit impulse at time = or the input signal Input
	bool impulse_response_or_input_signal = false;
//This determines if there is a random factor (in this case a random Gaussian variable) added to the node positions in the spider web so the web is slightly randomly shaken
	bool random_node_positions = false;
	
//This is the mean and standard deviation of the random Gaussian variable added to the node positions
	double mean =0;
	double stdev = 0;

	springs_identical = false;
    bias_learning = true;
	impulse_response_or_input_signal = true;
	
////////////////////////////////////////////////////////////////////////
//THIS SIMULATES THE SPIDER WEB
//////////////////////////////////////////////////////////////////////

//This is the overloaded constructor for the simulation class
//This initialises the spider web and loads the starting variables, the input and target signals into the Simulation object sim
 Simulation sim(radius, rounds, no_of_points_per_round, data, Input, Volterra, wash_out_time, learning_time, learning_time_test, springs_identical, random_node_positions, mean, stdev);

 

//This function carries out the simulation in time. 
//The parameters bias_learning and impulse_response_or_input_signal determines whether bias learning is on and whether the input signal is an impulse respnse.
 sim.update(bias_learning, impulse_response_or_input_signal);
 
 
//the return_thread_Number returns the number of the thread that has two particular nodes with certain node numbers
cout << sim.return_thread_Number(0, 1) << endl;
cout << sim.return_thread_Number(1, 2) << endl;
cout << sim.return_thread_Number(2, 3) << endl;
cout << sim.return_thread_Number(3, 4) << endl;
cout << sim.return_thread_Number(4, 0) << endl;
cout << sim.return_thread_Number(0, 10) << endl;
cout << sim.return_thread_Number(1, 10) << endl;
cout << sim.return_thread_Number(2, 10) << endl;
cout << sim.return_thread_Number(3, 10) << endl;
cout << sim.return_thread_Number(4, 10) << endl;
 


  
//This resets the simulation so that all positions of threads and nodes are what they were at time t=0.
sim.Reset_Simulation();
impulse_response_or_input_signal = true;
 //  sim.update(bias_learning, impulse_response_or_input_signal);
   
//This outputs the mean squared error and loads the Output vector and Target vector in the Simulatino Object RENAME THIS IT DOES NOT OUTPUT LEARNING MATRIX
sim.output_LearningMatrix_and_MeanSquaredError();

//This outputs the output signal, target signal (to check if loaded in correctly) and the learning matrix.
sim.output_Output_Signal(str);
  
  //This is to measure how long the code has taken in total to execute  
auto end = std::chrono::high_resolution_clock::now();
cout << "The time it took for the programme to run in total in milliseconds: ";
std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms";

return 0;
}

