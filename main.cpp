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

   auto begin = std::chrono::high_resolution_clock::now();


    InitialDataValues data;
    //send from processor
    srand(rdtsc());

    vector<double> Volterra;
    vector<double> Input;
    vector<string> Input_Lines;
    vector<string> Volterre_Lines;  // Todo: ???

    cout << "-- Start ---------------------------------------- " << endl;


    // Todo: I think it would be useful to define class for reading in files from csv, maybe best even encapsulated in a class
    // there are plenty of these out there
  //  ifstream file_Input ( "/Users/hh14185/Leverhulme_Trust_Proposal_spider_web/Xcode/Nonlinear-Mass-Spring-System-BR/input.csv" );
    ifstream file_Input ( "Data/inputsignal.csv" ); file_Input.precision(15);
    string tmp;

    while (getline(file_Input, tmp,'\n'))
    {
        Input_Lines.push_back(tmp); //Get each line of the file as a string
    }

    unsigned long s = Input_Lines.size();
    double x;
    for (unsigned int i=1; i<s; ++i)
    {
        size_t pos = Input_Lines[i].find(",");      // position of the end of the name of each one in the respective string
        x = stod(Input_Lines[i].substr(pos+1,Input_Lines[i].size()));
        x = x*1;
        Input.push_back(x); // convert string age to a double
    }



    // TODO:  again best to have that encapsulated in a class - makes the code much cleaner and we will use this code a lot
    //  ifstream file2 ("inputsignal.csv"); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/

  //  ifstream file_Volterra ("/Users/hh14185/Leverhulme_Trust_Proposal_spider_web/Xcode/Nonlinear-Mass-Spring-System-BR/volterra.csv");
  
  //This takes in Volterra.CSV. will change this to narma and it should be the same.
   //   ifstream file_Volterra ("Data/narma.csv");
     ifstream file_Volterra ("Data/volterra.csv");

    while (getline(file_Volterra, tmp,'\n'))
    {
        Volterre_Lines.push_back(tmp); //Get each line of the file as a string
    }

    s = Volterre_Lines.size();
    for (unsigned int i=1; i<s; ++i)
    {

        size_t pos = Volterre_Lines[i].find(",");      // position of the end of the name of each one in the respective string
        x = stod(Volterre_Lines[i].substr(pos+1,Volterre_Lines[i].size()));
        x = x*1;
        Volterra.push_back(x); // convert string age to a doubl
    }


    cout <<"Size of input signal is: "<< Input.size() << endl;
    cout <<"Size of target signal is: "<< Volterra.size() << endl;

    
	//Original Volterra wash out time etc.
	
	double wash_out_time = 20000;
  //  wash_out_time = 20000;
    double learning_time = 200000;
    double learning_time_test = 15000;


    //Wash out time, learning_time and learning_time_test
 //   wash_out_time = 3000;
 //   learning_time = 400000;
 //   learning_time_test = 15000;
    
    // setting parameters for simulation
    // This should be possible to read in from a text file
    data.N = 25;
    data.ux=0;
    data.uy= 0;

    data.input_connectivity = 0.05;
    //data.w_in_initial = -1;

//    data.min_input_weight = -1;
//    data.max_input_weight = 1;

    data.min_input_weight = -0.0001 * 1;
    data.max_input_weight = 0.0001 * 1;

    data.min_input_weight = -0.005 * 1;
    data.max_input_weight = -0.005 * 1;

    data.min_x_position = 0;
    data.max_x_position = 10;
    data.min_y_position  = 0;
    data.max_y_position = 10;
    /*

    data.min_k3 = 1;
    data.max_k3  = 100;

    //Increase d3 and d1.
    double range_d1_d3 = 0;
    data.min_d3 = 1+range_d1_d3;
    data.max_d3  = 100+range_d1_d3;

    data.min_k1 = 1;
    data.max_k1  = 200;

    //Increase d3 and d1.
    data.min_d1 = 1+range_d1_d3;
    data.max_d1  = 200+range_d1_d3;
    */

    data.min_k3 = 50;
    data.max_k3  = 50;

    data.min_d3 = 50;
    data.max_d3  =50;

    data.min_k1 = 100;
    data.max_k1  = 100;

    //Increase d3 and d1.
    data.min_d1 = 100;
    data.max_d1  = 100;
    
    data.mass_of_nodes = 1;
    data.scaling_factor = 0.1;


    data.dt = 0.001;
    data.t0 = wash_out_time*data.dt;
    data.tmax = (wash_out_time+learning_time+learning_time_test)*data.dt;

    vector<double> Sine_Wave;

    cout << "Initial Input is: "<< Input[0] << endl;

    double new_MSE = 0;
    double old_MSE = 100;
    vector<double> MSE_list;
    string st;




  double radius = 10.0;
  int rounds = 2;
//  int rounds = 24;
  int no_of_points_per_round= 5;

  string str = "1";
  new_MSE = 0;
  old_MSE = 0;

  vector<double> range_d1_d3_list;
  double best_range_d1_d3;

  double a = 100;


 data.min_k1 = 1;
data.max_k1  = 1;

data.min_k3 = 1;
data.max_k3  = 1;


  data.min_d1 = 1;
  data.max_d1  = 1;

  data.min_d3 = 1;
  data.max_d3  = 1;
  

    bool springs_identical = false;
    bool bias_learning = true;
	bool impulse_response_or_input_signal = false;
	bool random_node_positions = false;
	
	double mean =0;
	double stdev = 0;

	springs_identical = false;
    bias_learning = true;
	impulse_response_or_input_signal = true;

 Simulation sim(radius, rounds, no_of_points_per_round, data, Input, Volterra, wash_out_time, learning_time, learning_time_test, springs_identical, random_node_positions, mean, stdev);

 //Simulation sim(data, Input, Volterra, wash_out_time, learning_time, learning_time_test);
 

 sim.update(bias_learning, impulse_response_or_input_signal);
 
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
 
 

  
auto end = std::chrono::high_resolution_clock::now();
   cout << "The time it took for the programme to run in total in milliseconds: ";
   std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms";
  
  

  //sim.output_Output_Signal(str);
  
 // sim.Reset_Simulation();
  // impulse_response_or_input_signal = true;
 //  sim.update(bias_learning, impulse_response_or_input_signal);
   
  // sim.output_LearningMatrix_and_MeanSquaredError();

  

  
  
 // data.scaling_factor= 1/(pow(10,18));

  // data.scaling_factor = 1;
  
    data.min_k3 = 50;
    data.max_k3  = 50;

    data.min_d3 = 50;
    data.max_d3  =50;

    data.min_k1 = 100;
    data.max_k1  = 100;

    //Increase d3 and d1.
    data.min_d1 = 100;
    data.max_d1  = 100;
    
    //output k1, d1, k3, d3 etc.
    
    
    
    
    
    
    
    
    
    //ofstream k1_d1_k3_d3_stab("k1_d1_k3_d3_stab_2.csv"); k1_d1_k3_d3_stab.precision(15);
    
    
    
    
    
    
    
    
    
    
    
    
    
    //Spring and damping coefficietns and stability vector of vectors.
    vector<vector<double>> sdc_and_stability;
    


	
//	Simulation sim(data, Input, Volterra, wash_out_time, learning_time, learning_time_test);
	
	
	




  

//  cout <<"The number of nodes is: " << data.N << endl;
//  cout <<"The number of springs is: " << sim.Spring_List() << endl;

//     auto end = std::chrono::high_resolution_clock::now();
 //  cout << "The time it took for the programme to run in total in milliseconds: ";
//   std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms";




    return 0;
}

