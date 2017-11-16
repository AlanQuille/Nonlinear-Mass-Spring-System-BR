#include <iostream>
#include "Simulation.hpp"
#include "Eigen/QR"
using namespace std;
using namespace Eigen;


int main(int argc, char** argv)
{
  InitialDataValues data;

  data.N =3;
	data.ux=1;
	data.uy= 0;
	data.input_connectivity = 0.2;
	data.w_in_initial = -1;
	data.w_in_final = 1;
	data.w_out_initial = -1;
	data.w_out_final = 1;
	data.range0x = 0;
	data.range1x = 1;
	data.range0y = 0;
	data.range1y = 1;
	data.initial_log_uniform = 1;
	data.final_log_uniform = 100;
	data.initial_uniform = 100;
	data.final_uniform = 200;
  data.t0 = 0;
	data.tmax = 1;
	data.dt = 0.001;

  Simulation sim(data);

  cout << endl;

  for(int i =0; i<sim.Spring_List(); i++)
  {
    cout <<"Weight for spring" <<" " << sim.Spring_Return(i).Output_Weight() << endl;
  }

//  sim.Springs

  double test = 0;

  int lessthan10 =0;
  int lessthan100 =0;
  int lessthan1000 = 0;

  sim.Log_10_Uniform(1, 100);
  for(int i=0; i<100000; i++)
  {
    test = sim.Log_10_Uniform(0, 1000);
    if(test <=10) lessthan10++;
    if(test <=100 && test>10) lessthan100++;
    if(test <=1000 && test>100) lessthan1000++;
  }
  cout <<"If the log to the base 10 uniform distribution is functioning correctly, these numbers should be approx. equal.";
  cout << endl;
  cout << lessthan10;
  cout << endl;
  cout << lessthan100;
  cout <<endl;
  cout << lessthan1000;
  cout <<endl;
 cout <<"Test to see if the Moore-Penrose pseudoinverse is functioning correctly.";
 MatrixXd m(2,2);
 m(0,0) = 1;
 m(1,0) = 0;
 m(0,1) = 1;
 m(1,1) = 0;
 cout << endl;
 cout << m << endl;
 m = m.completeOrthogonalDecomposition().pseudoInverse();
 cout << m << endl;

  return 0;
}
