#include <iostream>
#include <cstdlib>

#include "Simulation.hpp"

using namespace std;



unsigned long long rdtsc()
{
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}

int main(int argc, char** argv)
{

  InitialDataValues data;
  srand(rdtsc());

  data.N =3;
  data.ux=1;
  data.uy= 0;
  data.input_connectivity = 0.2;
  data.w_in_initial = -1;
  data.w_in_final = 1;
  data.w_out_initial = -1;
  data.w_out_final = 1;
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
  data.dt = 0.001;

  Simulation sim(data);


  return 0;
}
