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
  data.tmax = 10;
//  data.tmax = 1;
  data.dt = 0.001;



    vector<double> LotkaX;
    vector<double> LotkaY;
    vector<double> Sine_Wave;

    DynamicalSystems sys1(data.t0, data.tmax, data.dt);

    sys1.LotkaVolterra(LotkaX, LotkaY);
  //  sys1.SineWave(Sine_Wave);

    cout <<"WHat?" << endl;



    Simulation sim(data, LotkaX);
    sim.Output_Signal_And_MSE();

    stop_time = clock();
    double difference = (1000)*((stop_time - start_time)/CLOCKS_PER_SEC);

    cout << "The time it took for the programme to run in total in milliseconds: ";
    cout << difference << endl;

/*
  int rounds =3;
  double radius =1.0;
  int no_of_points_per_round = 8;

  data.initial_log_uniform = 1;
  data.final_log_uniform = 10;
  data.initial_uniform = 100;
  data.final_uniform = 200;
  data.t0 = 0;
//  data.tmax = 2*M_PI;
  data.ux=0;
  data.uy= 0;
  data.input_connectivity = 0.2;
//data.w_in_initial = -1;
  data.w_in_initial = -1;
  data.w_in_final = 1;
  data.w_out_initial = -1;
  data.w_out_final = 1;

  data.tmax = 1;
  data.dt = 0.001;


  vector<double> LW;
  vector<double> Output_Signal;

  //Lotka Volterra vectors
  vector<double> LotkaX;
  vector<double> LotkaY;

  DynamicalSystems sys1(data.t0, data.tmax, data.dt);
  sys1.LotkaVolterra(LotkaX, LotkaY);

  cout<<"Lotka initial is: "<<LotkaX[50] << endl;
  MatrixXd LM;
  //Simulation sim(radius, rounds, no_of_points_per_round, data, LotkaX, LotkaY);
//  sim.Output_Signal_And_MSE();

  stop_time = clock();
  double difference = (1000)*((stop_time - start_time)/CLOCKS_PER_SEC);

  cout << "The time it took for the programme to run in total in milliseconds: ";
  cout << difference << endl;
*/

  return 0;
}
