#define _USE_MATH_DEFINES
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "Simulation.cpp"
#include "Eigen/Dense"
#include "Eigen/QR"

using namespace std;
using namespace Eigen;


//This takes the pseudo-cycles from the processor for srand()
unsigned long long rdtsc()
{
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}


int main(int argc, char** argv)
{
  clock_t start_time,stop_time;
  start_time = clock();

  InitialDataValues data;
  //send from processor

  srand(rdtsc());


  vector<double> Volterra;
  vector<string> classData;

  ifstream file ( "volterra.csv" ); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
  string value;

  while (getline(file, value,'\n'))
  {
    classData.push_back(value); //Get each line of the file as a string
  }

  int s = classData.size();
  double x;
  for (unsigned int i=1; i<s; ++i)
  {

    size_t pos = classData[i].find(";");      // position of the end of the name of each one in the respective string
    x = stod(classData[i].substr(pos+1,classData[i].size()));
    Volterra.push_back(x); // convert string age to a double
  }




  data.N =15;
  data.ux=1;
  data.uy= 0;
  data.input_connectivity = 0.2;
  //data.w_in_initial = -1;
  data.w_in_initial = -1;
  data.w_in_final = 1;
  data.w_out_initial = -1;
  data.w_out_final = 1;
  //range-doubled, was fooling around with something, didn't change  back. Will change later.
  data.range0x = 0;
  data.range1x = 10;
  data.range0y = 0;
  data.range1y = 10;
  data.initial_log_uniform = 1;
  data.final_log_uniform = 10;
  data.initial_uniform = 100;
  data.final_uniform = 200;
  data.t0 = 0;
  data.tmax = 1;
//  data.tmax = 1;
  data.dt = 0.001;


  //  sys1.LotkaVolterra(LotkaX, LotkaY);
  //  sys1.SineWave(Sine_Wave)

    //Proportion of input signal that is used to generate learning weights.
    x = (int)(0.666666*((data.tmax - data.t0)/data.dt));

    int maxtimesteps = (int)((data.tmax-data.t0)/data.dt);

    std::vector<double> Volterra2(Volterra.begin(), Volterra.end() - 0.999*Volterra.size());

    cout <<"x is: " << x << endl;

    cout << Volterra2.size() << endl;




    Simulation sim(data, Volterra2, x);
    sim.Output_Signal_And_MSE();

    stop_time = clock();
    double difference = (1000)*((stop_time - start_time)/CLOCKS_PER_SEC);

    cout << "The time it took for the programme to run in total in milliseconds: ";
    cout << difference << endl;
//  int rounds =1;
//  double radius =1.0;
//  int no_of_points_per_round = 5;


  //send from processor
  /*
  srand(rdtsc());


  vector<double> Volterra;
  vector<string> classData;

  ifstream file ( "volterra.csv" ); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
  string value;

  while (getline(file, value,'\n'))
  {
    classData.push_back(value); //Get each line of the file as a string
  }

  int s = classData.size();
  double x;
  for (unsigned int i=1; i<s; ++i)
  {

    size_t pos = classData[i].find(";");      // position of the end of the name of each one in the respective string
    x = stod(classData[i].substr(pos+1,classData[i].size()));
    Volterra.push_back(x); // convert string age to a double
  }
*/

/*
  cout <<"Size of Volterra is: " << Volterra.size() << endl;

  data.initial_log_uniform = 1;
  data.final_log_uniform = 10;
  data.initial_uniform = 100;
  data.final_uniform = 200;
  //Washout time = 200ms
  data.t0 = 0.2;
//  data.tmax = 2*M_PI;
  data.ux=0.0;
  data.uy= 0;
  data.input_connectivity = 0.2;
//data.w_in_initial = -1;
  data.w_in_initial = -1;
  data.w_in_final = 1;
  data.w_out_initial = -1;
  data.w_out_final = 1;

  data.tmax = 33.33;
  data.dt = 0.001;



  vector<double> LW;
  vector<double> Output_Signal;

  //Lotka Volterra vectors
  vector<double> LotkaX;
  vector<double> LotkaY;

  vector<double> LeWe;


  DynamicalSystems sys1(data.t0, data.tmax, data.dt);
  sys1.LotkaVolterra(LotkaX, LotkaY);

  cout<<"Lotka initial is: "<<LotkaX[50] << endl;
  MatrixXd LM;
  Simulation sim(radius, rounds, no_of_points_per_round, data, Volterra);

  sim.Output_Signal_And_MSE();

  data.t0 = 33.33;
  data.tmax = 50.199;
  cout <<"Tmax is: " <<data.tmax << endl;


  Simulation sim2(radius, rounds, no_of_points_per_round, data, Volterra);

  LeWe = sim.Return_Learning_Weights();

  sim2.Output_Signal_And_MSE(LeWe);
  cout <<"Size of old learning weights" <<  sim.Return_Learning_Weights().size() << endl;
  cout  <<"Size of cols" << sim2.Return_Learning_Matrix().cols() << endl;


*/

//  stop_time = clock();
//  double difference = (1000)*((stop_time - start_time)/CLOCKS_PER_SEC);

//  cout << "The time it took for the programme to run in total in milliseconds: ";
//  cout << difference << endl;


  return 0;
}
