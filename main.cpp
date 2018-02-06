#define _USE_MATH_DEFINES
//#define EIGEN_USE_BLAS
//#define EIGEN_USE_LAPACKE
//#define EIGEN_USE_LAPACKE_STRICT
#include <iostream>
#include <complex>
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



//Volterra target signal and input signal.
  vector<double> Volterra;
  vector<double> Input_Signal;
  vector<string> classData;
  vector<string> classData2;

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
  //  cout <<Volterra.at(i) << endl;
  }

  ifstream file2 ("inputsignal.csv"); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
  string value2;

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




  data.N =40;
  data.ux=0;
  data.uy= 0;
  data.input_connectivity = 0.2;
  //data.w_in_initial = -1;
  data.w_in_initial = -3;
  data.w_in_final = 3;
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
  data.tmax = 0.1;
//  data.tmax = 1;
  data.dt = 0.001;


  //  sys1.LotkaVolterra(LotkaX, LotkaY);
  //  sys1.SineWave(Sine_Wave)

    //Proportion of input signal that is used to generate learning weights.
    //x = (int)(0.666666*((data.tmax - data.t0)/data.dt));

    int maxtimesteps = (int)((data.tmax-data.t0)/data.dt);

    double sum2 = std::accumulate(Input_Signal.begin(), Input_Signal.end(), 0.0);
    double mean2 = sum2/ Input_Signal.size();

    double sq_sum2 = std::inner_product(Input_Signal.begin(), Input_Signal.end(), Input_Signal.begin(), 0.0);
    double stdev2 = std::sqrt(sq_sum2 /Input_Signal.size() - mean2*mean2);

    double sum = std::accumulate(Volterra.begin(), Volterra.end(), 0.0);
    double mean = sum / Volterra.size();

    double sq_sum = std::inner_product(Volterra.begin(), Volterra.end(), Volterra.begin(), 0.0);
    double stdev = std::sqrt(sq_sum /Volterra.size() - mean * mean);



  //  for( int i =0; i<Input_Signal.size(); i++)
  //  {
    // Input_Signal.at(i) = (Input_Signal.at(i) - mean2)/(stdev2);
    //}




    //std::vector<double> Volterra2(Volterra.begin()+0.5*Volterra.size(), Volterra.end() - 0.49*Volterra.size());
    std::vector<double> Volterra2(Volterra.begin()+0, Volterra.begin()+235000);
    std::vector<double> Input_Signal2(Input_Signal.begin()+0, Input_Signal.begin()+235000);



    cout << mean << endl;
    cout << stdev << endl;

    cout <<"is it normalised? " << endl;




    sum = std::accumulate(Volterra2.begin(), Volterra2.end(), 0.0);
    mean = sum / Volterra2.size();
    sq_sum = std::inner_product(Volterra2.begin(), Volterra2.end(), Volterra2.begin(), 0.0);
    stdev = std::sqrt(sq_sum / Volterra2.size() - mean * mean);
    cout <<"Volterra mean and sd is:" << mean << endl;
    cout <<"Volterra mean and sd is:" << stdev << endl;


  //  std::vector<double> Volterra3(Volterra.begin()+x*(1-0.3333333)*Volterra.Size(), x*Volterra.end());

  cout <<"Size of input vector is: " << Volterra2.size() << endl;

    cout <<"x is: " << x << endl;

    //cout << Volterra2.size() << endl;

    double twothirdsprotocol = 0.930232;


   //Testing eigen lapacke
    Simulation sim(data, Volterra2, twothirdsprotocol, Input_Signal2);
    sim.Output_Signal_And_MSE();

    stop_time = clock();
    double difference = (1000)*((stop_time - start_time)/CLOCKS_PER_SEC);

    cout << "The time it took for the programme to run in total in milliseconds: ";
    cout << difference << endl;
    cout <<Volterra2.size() << endl;
    cout <<Volterra.size() << endl;


  MatrixXd m(3,2);

  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(2,1) = 50;
  m(2,0) = 40;
  m(1,1) = m(1,0) + m(0,1) +5;

  std::cout << m << std::endl;
  start_time = clock();
  cout << m.completeOrthogonalDecomposition().pseudoInverse() << endl;
  stop_time = clock();
  difference = (1000)*((stop_time - start_time)/CLOCKS_PER_SEC);
  cout <<"The time is took for the original pseduoinverse is: " << difference << endl;
  //Eigen::P =  m.colPivHouseholderQr();
  //cout <<

  HouseholderQR<MatrixXd> qr(m);
  cout << endl;
  MatrixXd Q = qr.householderQ();
  cout << Q << endl;
  MatrixXd R = qr.matrixQR().triangularView<Upper>();
  cout <<R;

  cout << (R.transpose()*R) * R.transpose();
  cout << endl;
  cout <<Q.transpose();
  cout << endl;
  //cout << ((R.transpose()*R) * R.transpose())*Q.transpose();
//  .inverse()*R.tranpose();
 start_time = clock();
  cout <<Q*R << endl;

  //m.transposeInPlace();

  cout << m;
  m.transposeInPlace();
  cout << endl;
  cout << m;
  MatrixXd n(3,2);
  n(0,0) = 3;
  n(1,0) = 2.5;
  n(0,1) = -1;
  n(2,1) = 50;
  n(2,0) = 40;
  n(1,1) = n(1,0) + n(0,1) +5;

  cout << endl << ((m*n).inverse())*m;


  //cout << m*((m * n).inverse() );
//  cout << endl;
  //cout << m.transpose() * ((m.transpose() * m).inverse());
//  cout << endl;
  stop_time = clock();
  difference = (1000)*((stop_time - start_time)/CLOCKS_PER_SEC);
  cout <<"The time is took for the new method is: " << difference << endl;
  //cout <<R.inverse();


















    // Initialize the random generator
  //  arma::arma_rng::set_seed_random();

      // Create a 4x4 random matrix and print it on the screen
  //  arma::Mat<double> A = arma::randu(4,4);
  //  std::cout << "A:\n" << inv(A) << "\n";
    //arma::Mat<double> B = A * A;

    //arma::Mat<double> B = arma::pinv(A);        // use default tolerance




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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
