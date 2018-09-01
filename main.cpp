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
//Matlab
#include "mat.h"
#include "matrix.h"
#include "tmwtypes.h"



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

  //Matlab inputs
  MATFile * pmat; char ** you ; int ndir; int ndir2; mxArray * pa; mxArray * qa; mxArray *ra; mxArray *sa; mxArray *ta; mxArray *ua; mxArray *va; mxArray *xa;  mxArray *ya; mxArray *za; mxArray *zb; mxArray *in; mxArray *in2; const char * name; char ** dir;
  char ** dir2;
  pmat = matOpen ("init_net.mat", "r" ) ;
  //dir = matGetDir ( pmat, &ndir ) ;
  //double * paData;
  //Double and bool vectors for .mat input.
  double * qaData;
  double * raData;
  double * saData;
  double * taData;
  double * uaData;
  double * xaData;
  double * yaData;

  double * zaData;
  double * zbData;

  double * inData;

  double * inData2;

  //Get positions of nodes
pa = matGetVariable(pmat, "P");
qa = mxGetField(pa, 0, "states");
in = mxGetField(pa, 0, "fixed");

//For input nodes
in2 = matGetVariable(pmat, "W_in");

//Get spring variables.
ra = matGetVariable(pmat, "W");
sa = mxGetField(ra, 0, "k1");
ta = mxGetField(ra, 0, "k3");
ua = mxGetField(ra, 0, "d1");
xa = mxGetField(ra, 0, "d3");
ya = mxGetField(ra, 0, "l0");

za = mxGetField(ra, 0, "from");
zb = mxGetField(ra, 0, "to");

//in = mxGetField(ra, 0, "fixed");

qaData = ( double * ) mxGetData ( qa ) ;

raData = ( double * ) mxGetData ( ra ) ;
saData = ( double * ) mxGetData ( sa ) ;
taData = ( double * ) mxGetData ( ta ) ;
uaData = ( double * ) mxGetData ( ua ) ;
xaData = ( double * ) mxGetData ( xa ) ;
yaData = ( double * ) mxGetData ( ya ) ;

zaData = ( double * ) mxGetData ( za ) ;
zbData = ( double * ) mxGetData ( zb ) ;

inData = ( double * ) mxGetData ( in ) ;
inData2 = ( double * ) mxGetData ( in2 ) ;
//X positions of nodes
vector<double> x_nodes;
vector<double> y_nodes;
//input nodes.
//vector<bool> input_nodes;

//Get k1, k3, d1, d3, l0, from, to
vector<double> k1;
vector<double> k3;
vector<double> d1;
vector<double> d3;
vector<double> l0;
vector<int> node1;
vector<int> node2;

vector<bool> fixed_nodes;
vector<double> W_in;

for(int i=0; i<78; i++)
{
  if(i<30)
  {
  x_nodes.push_back(qaData[i]);
//  cout << x_nodes[i] << endl;
  }
  if(i>=30 && i<60)
  {
  y_nodes.push_back(qaData[i]);
//  cout << y_nodes[i] << endl;
  }

  k1.push_back(saData[i]);
//  cout << k1[i] << endl;
  k3.push_back(taData[i]);
//  cout << k3[i] << endl;
  d1.push_back(uaData[i]);
//  cout << d1[i] << endl;
  d3.push_back(xaData[i]);
//  cout << d3[i] << endl;
  l0.push_back(yaData[i]);
//  cout << l0[i] << endl;;

  node1.push_back((int)zaData[i]);
  //cout << node1[i] << endl;
  node2.push_back((int)zbData[i]);
//  cout << node2[i] << endl;

  if(i<30)
   {
  fixed_nodes.push_back((bool)inData[i]);
  cout <<fixed_nodes[i] << endl;
  //cout << inData[i] << endl;
   }

   //Number of fixed nodes, manual load in not present yet.
   if(i<30)
   {
   W_in.push_back(inData2[i]);
   cout <<"W_in is" << W_in[i] << endl;
   }

}




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
      ifstream file_Volterra ("Data/volterra.csv");

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

    double wash_out_time = 20000;
    double learning_time = 200000;
    double learning_time_test = 15000;


    // setting parameters for simulation
    // This should be possible to read in from a text file
    data.N = 27;
    data.ux=0;
    data.uy= 0;

    data.input_connectivity_percentage = 20;
    //data.w_in_initial = -1;
    data.min_input_weight = -1;
    data.max_input_weight = 1;
    data.min_x_position = 0;
    data.max_x_position = 10;
    data.min_y_position  = 0;
    data.max_y_position = 10;

    data.min_k3 = 1;
    data.max_k3  = 100;
    data.min_d3 = 1;
    data.max_d3  = 100;

    data.min_k1 = 1;
    data.max_k1  = 200;
    data.min_d1 = 1;
    data.max_d1  = 200;


    data.dt = 0.001;
    data.t0 = wash_out_time*data.dt;
    data.tmax = (wash_out_time+learning_time+learning_time_test)*data.dt;

    vector<double> Sine_Wave;

    cout << "Initial Input is: "<< Input[0] << endl;

    Simulation sim(Input, Volterra, wash_out_time, learning_time, learning_time_test, data.min_input_weight, data.max_input_weight, x_nodes, y_nodes, fixed_nodes, W_in, k1, k3, d1, d3, l0, node1, node2, data.dt);

   //Simulation sim(data, Input, Volterra, wash_out_time, learning_time, learning_time_test);
   //sim.Reset_Simulation();
   //sim.execute();
   //sim.output_LearningMatrix_and_MeanSquaredError();
/*
   sim.Reset_Simulation();
   sim.execute();
   sim.output_LearningMatrix_and_MeanSquaredError();

   sim.Reset_Simulation();
   sim.execute();
   sim.output_LearningMatrix_and_MeanSquaredError();

   sim.Reset_Simulation();
   sim.execute();
   sim.output_LearningMatrix_and_MeanSquaredError();
   */

  //cout <<"The number of nodes is: " << data.N << endl;
  //cout <<"The number of springs is: " << sim.Spring_List() << endl;

   auto end = std::chrono::high_resolution_clock::now();
   cout << "The time it took for the programme to run in total in milliseconds: ";
   std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms";

   sim.output_LearningMatrix_and_MeanSquaredError();


    return 0;
}
