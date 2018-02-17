#define _USE_MATH_DEFINES
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <chrono>
#include "Simulation.cpp"
#include "Eigen/Dense"
#include "Eigen/QR"


using namespace std;
using namespace Eigen;




/*
<<<<<<< HEAD
int main(int argc, char** argv)
{
  clock_t start_time,stop_time;
  start_time = clock();

  InitialDataValues data;
  //send from processor
  srand(rdtsc());

  vector<double> Volterra;
  vector<double> Input_Signal;
  vector<string> classData;
  vector<string> classData2;  // Todo: ???

  cout << "Test " << endl;


  // Todo: I think it would be useful to define class for reading in files from csv
  // there are plenty of these out there, e.g. https://github.com/iofish/CSV/blob/master/csv.h
  ifstream file ( "volterra.csv" ); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
  string value;

  while (getline(file, value,'\n'))
  {
    classData.push_back(value); //Get each line of the file as a string
  }

  unsigned long s = classData.size();
  double x;
  for (unsigned int i=1; i<s; ++i)
  {

    size_t pos = classData[i].find(";");      // position of the end of the name of each one in the respective string
    x = stod(classData[i].substr(pos+1,classData[i].size()));
    Volterra.push_back(x); // convert string age to a double
  //  cout <<Volterra.at(i) << endl;
  }


again best to have that encapsulated in a class - makes the code much cleaner and we will use this code a lot
  ifstream file2 ("inputsignal.csv"); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
    string value2;  // Todo: naming value2 and file2 is quite cryptic.

  while (getline(file2, value2,'\n'))
  {
    classData2.push_back(value2); //Get each line of the file as a string
  }

  s = classData2.size();
  for (unsigned int i=1; i<s; ++i)
  {

    size_t pos = classData2[i].find(";");      // position of the end of the name of each one in the respective string
    x = stod(classData2[i].substr(pos+1,classData2[i].size()));
    Input_Signal.push_back(x); // convert string age to a doubl
    //cout <<Input_Signal.at(i) << endl;
  }
  cout <<"Size of new signal is: "<< Input_Signal.size();


  // setting parameters for simulation
  data.N = 50;

  data.ux=0;
  data.uy= 0;

  data.input_connectivity_percentage = 0.2;
  //data.w_in_initial = -1;
  data.input_weight_smallest_value = -1;
  data.input_weight_largest_value = 1;

  //range-doubled, was fooling around with something, didn't change  back. Will change later.
  data.smallest_x_position = 0;
  data.largest_x_position = 10;
  data.smallest_y_position  = 0;
  data.largest_y_position = 10;

  data.log_uniform_smallest_value = 1;
  data.log_uniform_largest_value  = 10;
  data.uniform_smallest_value = 100;
  data.uniform_largest_value= 200;

  data.t0 = 0;
  data.tmax =50;
  //  data.tmax = 1;
  data.dt = 0.001;
=======
>>>>>>> 2a3190bd33ae3eb03a7a23043ecc19180087e5bf

*/

//This takes the pseudo-cycles from the processor for srand()

unsigned long long rdtsc()
{
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}


int main(int argc, char** argv)
{
// ways to get current directory
//    ostringstream fullpath;
//    fullpath << argv[0]; // add current working path
//    cout <<"pwd: " << argv[0] << endl;


//    clock_t start_time,stop_time;
//    start_time = clock();

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
    ifstream file_Input ( "inputsignal.csv" );
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
        Volterra.push_back(x); // convert string age to a double
    }



    // TODO:  again best to have that encapsulated in a class - makes the code much cleaner and we will use this code a lot
    //  ifstream file2 ("inputsignal.csv"); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/

  //  ifstream file_Volterra ("/Users/hh14185/Leverhulme_Trust_Proposal_spider_web/Xcode/Nonlinear-Mass-Spring-System-BR/volterra.csv");
      ifstream file_Volterra ("volterra.csv");

    while (getline(file_Volterra, tmp,'\n'))
    {
        Volterre_Lines.push_back(tmp); //Get each line of the file as a string
    }

    s = Volterre_Lines.size();
    for (unsigned int i=1; i<s; ++i)
    {

        size_t pos = Volterre_Lines[i].find(",");      // position of the end of the name of each one in the respective string
        x = stod(Volterre_Lines[i].substr(pos+1,Volterre_Lines[i].size()));
        Input.push_back(x); // convert string age to a doubl
    }


    cout <<"Size of input signal is: "<< Input.size() << endl;
    cout <<"Size of target signal is: "<< Volterra.size() << endl;


    // setting parameters for simulation
    // This should be possible to read in from a text file
    data.N = 50;
    data.ux=0;
    data.uy= 0;

    data.input_connectivity_percentage = 0.2;
    //data.w_in_initial = -1;
    data.min_input_weight = -1;
    data.max_input_weight = 1;

    //range-doubled, was fooling around with something, didn't change  back. Will change later.
    // Todo: Please change that!
    data.min_x_position = 0;
    data.max_x_position = 10;
    data.min_y_position  = 0;
    data.max_y_position = 10;

    data.min_log_uniform = 1;
    data.max_log_uniform  = 10;
    data.min_uniform = 100;
    data.max_uniform= 200;

    data.t0 = 0;
    data.tmax = 10;
    data.dt = 0.001;


    vector<double> Sine_Wave;

    Simulation sim(data, Input, Volterra);
    cout <<"The number of nodes is: " << data.N << endl;
    cout <<"The number of springs is: " << sim.Spring_List() << endl;

    sim.Output_Signal_And_MSE();   // Todo: what is that doing? Naming is confusing.



    //stop_time = clock();
  //  double difference = (1000)*((stop_time - start_time)/CLOCKS_PER_SEC);

   auto end = std::chrono::high_resolution_clock::now();
   cout << "The time it took for the programme to run in total in milliseconds: ";
   std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms";


    //cout << difference << endl;



    return 0;
}
