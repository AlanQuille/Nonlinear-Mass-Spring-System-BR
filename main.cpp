#include <iostream>
#include "Simulation.hpp"
using namespace std;

int main(int argc, char** argv)
{
  double N =3;
	double ux =1;
	double uy = 0;
	double input_connectivity = 0.2;
	double w_in_initial = -1;
	double w_in_final = 1;
	double w_out_initial = -1;
	double w_out_final = 1;
	double range0x = 0;
	double range1x = 1;
	double range0y = 0;
	double range1y = 1;
	double initial_log_uniform = 0;
	double final_log_uniform = 1;
	double initial_uniform = 0;
	double final_uniform = 1;
  double t0 = 0;
	double tmax = 0.1;
	double dt = 0.001;

  Simulation s = new Simulation(N, ux, uy, input_connectivity, w_in_initial, w_in_final, w_out_initial, w_out_final, range0x, range1x, range0y, range1y, initial_log_uniform, final_log_uniform, t0, tmax, dt)
  return 0;
}
