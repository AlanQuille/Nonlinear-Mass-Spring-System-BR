#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>

#include "Simulation.hpp"

using namespace std;


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


  data.N =10;
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
  data.range1x = 5;
  data.range0y = 0;
  data.range1y = 5;
  data.initial_log_uniform = 1;
  data.final_log_uniform = 10;
  data.initial_uniform = 100;
  data.final_uniform = 200;
  data.t0 = 0;
  data.tmax = 2*M_PI;
  data.dt = 0.001;

  //run same simulation x timeb_s

   //Do the simulation, get learning weights and than use it to get learning matrix.
  vector<double> LW;
  vector<double> Output_Signal;
  MatrixXd LM;
  Simulation sim(data);
  LW = sim.Return_Learning_Weights();
//  Simulation sim2(data);
  LM = sim.Return_Learning_Matrix();

  ofstream output("outputsignal.csv");

  int maxtimesteps = ((data.tmax - data.t0)/(data.dt));
  double outputsignal = 0;
  double currenttime = 0;
  for(int i=0; i<maxtimesteps; i++)
  {
  for(int j=0; j<LM.cols(); j++)
  {
    outputsignal += LW.at(j) * LM(i, j);
  }
  Output_Signal.push_back(outputsignal);
  cout <<"The output signal is "<<Output_Signal[i] << endl;
  currenttime = data.t0 + i*data.dt;
  output << currenttime<<"," << Output_Signal[i] << endl;
  outputsignal = 0;
  }



  cout << LW[0];
  cout <<  endl;


  stop_time = clock();
  double difference = (1000)*((stop_time - start_time)/CLOCKS_PER_SEC);

  cout << "The time it took for the programme to run in total in milliseconds: ";
  cout << difference << endl;


  return 0;
}
