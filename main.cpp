#define _USE_MATH_DEFINES
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




void LotkaVolterra(vector<double> &LVx, vector<double> &LVy, int maxtimesteps, double dt)
{
  //USe Euler's to get quick Lotka Volterra. Parameters = 1, 1 just for speed.
  double x0;
  double xnext;
  double y0;
  double ynext;

  //Initial conditions
  x0 = 1.2;
  y0 = 1.2;

  LVx.push_back(x0);
  LVy.push_back(y0);

  for(int i=1; i<maxtimesteps; i++)
  {

    xnext = x0 + dt*(x0 - x0*y0);
    ynext = y0 + dt*(x0*y0 - y0);
    x0 = xnext;
    y0 = ynext;
    LVx.push_back(xnext);
    LVy.push_back(ynext);
  }
}


int main(int argc, char** argv)
{
  clock_t start_time,stop_time;
  start_time = clock();

  InitialDataValues data;
  //send from processor
  srand(rdtsc());


  data.N =70;
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
//  data.tmax = 2*M_PI;
  data.tmax = 0.002;
  data.dt = 0.001;

  //run same simulation x timeb_s

   //Do the simulation, get learning weights and than use it to get learning matrix.
   /*
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
  */

  int rounds =3;
  double radius =1.0;
  int no_of_points_per_round = 8;
  int maxtimesteps = ((data.tmax - data.t0)/(data.dt));

  vector<double> LW;
  vector<double> Output_Signal;

  //Lotka Volterra vectors
  vector<double> LotkaX;
  vector<double> LotkaY;

  LotkaVolterra(LotkaX, LotkaY, maxtimesteps, data.dt);

  cout<<"Lotka initial is: "<<LotkaX[50] << endl;
  MatrixXd LM;
  Simulation sim(radius, rounds, no_of_points_per_round, data, LotkaX, LotkaY);

  cout <<"Here?" <<endl;

  LW = sim.Return_Learning_Weights();
//  Simulation sim2(data);
  LM = sim.Return_Learning_Matrix();


  ofstream output("outputsignal.csv");
  ofstream output2("learningweights.csv");

  double outputsignal = 0;
  double currenttime = 0;
  for(int i=0; i<maxtimesteps; i++)
  {
  for(int j=0; j<LM.cols(); j++)
  {
    outputsignal += LW.at(j) * LM(i, j);
    if(i==0) output2 << LW[j] << endl;
  }
  Output_Signal.push_back(outputsignal);
  currenttime = data.t0 + i*data.dt;
  cout << outputsignal;
  cout << endl;
  output << currenttime <<"," << Output_Signal.at(i);
  output << endl;
  outputsignal = 0;

  }


  stop_time = clock();
  double difference = (1000)*((stop_time - start_time)/CLOCKS_PER_SEC);

  cout << "The time it took for the programme to run in total in milliseconds: ";
  cout << difference << endl;


  return 0;
}
