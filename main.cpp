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

  // auto begin = std::chrono::high_resolution_clock::now();


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
      ifstream file_Volterra ("Data/narma.csv");
    //  ifstream file_Volterra ("Data/volterra.csv");

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
    wash_out_time = 3000;
    learning_time = 400000;
    learning_time_test = 15000;
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

  //  Simulation sim(data, Input, Volterra, wash_out_time, learning_time, learning_time_test);
/*
  for(int i=0; i<10; i++)
  {

  //Simulation sim(data, Input, Volterra, wash_out_time, learning_time, learning_time_test);
  new_MSE = sim.return_MSE();
  MSE_list.push_back(new_MSE);
  st = to_string(new_MSE);

  if(old_MSE>new_MSE)
    {
  sim.Output_For_Plot();
  sim.output_Output_Signal(st);
  old_MSE = new_MSE;
  cout << "New one" << endl;
  cout << "New one" << endl;
  cout << "New one" << endl;
  cout << "New one" << endl;
    }

  }

  cout <<"The best MSE is: " << old_MSE << endl;
  cout <<"The average MSE is: " << accumulate( MSE_list.begin(), MSE_list.end(), 0.0)/ MSE_list.size() << endl;
*/
  // sim.Reset_Simulation();
//   sim.execute();
//   sim.output_LearningMatrix_and_MeanSquaredError();


/*
  cout <<"The number of nodes is: " << data.N << endl;
  cout <<"The number of springs is: " << sim.Spring_List() << endl;

  data.min_d3 = 1+range_d1_d3;
  data.max_d3  = 100+range_d1_d3;
*/




  double radius = 10.0;
  int rounds = 2;
//  int rounds = 24;
  int no_of_points_per_round= 5;

  string str = "1";
  new_MSE = 0;
  old_MSE = 0;

  vector<double> range_d1_d3_list;
  double best_range_d1_d3;

//  range_d1_d3 = 5000000;
  //best_range_d1_d3 = 5000;
  //ange_d1_d3 = 0;
  double a = 100;

//  range_d1_d3 += a;

//data.k_lim = [10 50;0.1 1];

/*
data.k_lim = [10 50;0.1 1];
data.d_lim = [10 500;100 500];
*/
/*
data.min_k1 = 10;
data.max_k1  = 50;

data.min_k3 = 0.1;
data.max_k3  = 1;
*/
/*
data.min_k3 = 10;
data.max_k3  = 50;

data.min_k1 = 0.1;
data.max_k1  = 1;

  data.min_d1 = 100;
  data.max_d1  = 500;

  data.min_d3 = 10;
  data.max_d3  = 500;
  */

  //data.min_d1 = 100;
//  data.max_d1  = 20000;

/*
  data.min_d3 = 100000;
  data.max_d3  =1000000;

  data.min_d1 = 100000;
  data.max_d1  = 1000000;
  */

//Input, Volterra.
//Volterra, Input

// Simulation sim(radius, rounds, no_of_points_per_round, data, Input, Volterra, wash_out_time, learning_time, learning_time_test);

  //sim.output_LearningMatrix_and_MeanSquaredError();
 // sim.output_Output_Signal(str);
  

   
   ofstream MSE_and_scaling_factor("MSE_and_scaling_factor_1.csv"); MSE_and_scaling_factor.precision(15);
  data.scaling_factor = 1;
  
  
  //Simulation for plotting scaling factor with MSE;
  
  /*
  
  for(int i=0; i<20; i++)
  {  
  data.scaling_factor= 0.0000000000000000001*pow(10, i);
  cout << data.scaling_factor << endl;
  Simulation sim(radius, rounds, no_of_points_per_round, data, Input, Volterra, wash_out_time, learning_time, learning_time_test);  
  MSE_and_scaling_factor << data.scaling_factor << "," << sim.return_MSE() << endl;
 // delete sim;
  }
  
  ofstream MSE_and_scaling_factor2("MSE_and_scaling_factor_0.5.csv"); MSE_and_scaling_factor.precision(15);
  data.scaling_factor = 0.5;
  
  for(int i=0; i<20; i++)
  {  
  data.scaling_factor= 0.0000000000000000001*pow(10, i);
  cout << data.scaling_factor << endl;
  Simulation sim(radius, rounds, no_of_points_per_round, data, Input, Volterra, wash_out_time, learning_time, learning_time_test);  
  MSE_and_scaling_factor2 << data.scaling_factor << "," << sim.return_MSE() << endl;
 // delete sim;
  }
  
  ofstream MSE_and_scaling_factor3("MSE_and_scaling_factor_0.12.csv"); MSE_and_scaling_factor.precision(15);
  data.scaling_factor = 0.12;
  
  for(int i=0; i<20; i++)
  {  
  data.scaling_factor= 0.0000000000000000001*pow(10, i);
  cout << data.scaling_factor << endl;
  Simulation sim(radius, rounds, no_of_points_per_round, data, Input, Volterra, wash_out_time, learning_time, learning_time_test);  
  MSE_and_scaling_factor3 << data.scaling_factor << "," << sim.return_MSE() << endl;
 // delete sim;
  }
  
  */
  
  
  data.scaling_factor= 1/(pow(10,18));

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
    
    ofstream k1_d1_k3_d3_stab("k1_d1_k3_d3_stab.csv"); k1_d1_k3_d3_stab.precision(15);
    
    //Spring and damping coefficietns and stability vector of vectors.
    vector<vector<double>> sdc_and_stability;
    
	
	//This adds in spring and damping coefficients and stabilityu 
	
	
	bool stab = true;
	int returnk = 0;
	
	auto begin = std::chrono::high_resolution_clock::now();
	
	
    
    for(int i=0; i<10; i++)
    {
    	for(int j=0; j<10; j++)
    	{
    		for(int k=0; k<10; k++)
    		{
    			for(int l=0; l<10; l++)
    			{
    				   //auto begin = std::chrono::high_resolution_clock::now();
    				
    				data.min_k1 = 0.0000001 * pow(10, i);
    				data.max_k1 = 0.0000001 *  pow(10, i);
    				
                    data.min_k3 = 0.0000001 * pow(10, j);
    				data.max_k3 = 0.0000001 * pow(10, j);
    				
    				data.min_d1 = 0.0000001 * pow(10, k);
    				data.max_d1 =  0.0000001 * pow(10, k);
    				
    			    data.min_d3 = 0.0000001* pow(10, l);
    				data.max_d3 = 0.0000001* pow(10, l);
    				
    				Simulation sim(radius, rounds, no_of_points_per_round, data, Input, Volterra, wash_out_time, learning_time, learning_time_test); 
    				
    				stab = sim.Stability_return();
    				
    				k1_d1_k3_d3_stab <<data.min_k1 <<"," << data.min_k3 <<"," <<data.min_d1 <<"," << data.min_d3 <<","<< stab << endl;
    				

					 

				}
			}
		}
	}
	
	   auto end = std::chrono::high_resolution_clock::now();
	   cout << "The time it took for the programme to run in total in milliseconds: ";
       std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms";
	

	
	/*
	data.min_k3 =0.0000001;
    data.max_k3  = 0.0000001;

    data.min_d3 = 0.0000001;
    data.max_d3  =0.0000001;

    data.min_k1 = 0.0000001;
    data.max_k1  = 0.0000001;

    //Increase d3 and d1.
    data.min_d1 = 0.0000001;
    data.max_d1  = 0.0000001;

  */
  
  
 //  Simulation sim(radius, rounds, no_of_points_per_round, data, Input, Volterra, wash_out_time, learning_time, learning_time_test); 
   
 //  cout <<"Is this structure stable: " << sim.Stability_return() << endl;
   
 //   sim.output_LearningMatrix_and_MeanSquaredError();
  //  sim.output_Output_Signal(str);
   





  

//  cout <<"The number of nodes is: " << data.N << endl;
//  cout <<"The number of springs is: " << sim.Spring_List() << endl;

 //  auto end = std::chrono::high_resolution_clock::now();
 //  cout << "The time it took for the programme to run in total in milliseconds: ";
 //  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms";


    return 0;
}
