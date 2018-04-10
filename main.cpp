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
#include <cstdio>
//For .mat file, maybe only works with C
//#include <string.h> /* For strcmp() */
//#include <stdlib.h> /* For EXIT_FAILURE, EXIT_SUCCESS */
#include <cstring> /* For strcmp() */
#include <vector> /* For STL */
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
  //This code block is for .Mat input.
  /*
  *
 *   matcreat
 *
 * Create a MAT-file which can be loaded into MATLAB.
 *
 * This program demonstrates the use of the following functions:
 *
 *  matClose
 *  matGetVariable
 *  matOpen
 *  matPutVariable
 *  matPutVariableAsGlobal
 *
 * Copyright 1984-2007 The MathWorks, Inc.
 */

//#define BUFSIZE 256
/*
  MATFile *pmat;
  mxArray *pa1, *pa2, *pa3;
  std::vector<int> myInts;
  myInts.push_back(1);
  myInts.push_back(2);
  printf("Accessing a STL vector: %d\n", myInts[1]);

  double data[9] = { 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0 };
  const char *file = "mattest.mat";
  char str[BUFSIZE];
  int status;

  printf("Creating file %s...\n\n", file);
  pmat = matOpen(file, "w");
  if (pmat == NULL) {
    printf("Error creating file %s\n", file);
    printf("(Do you have write permission in this directory?)\n");
    return(EXIT_FAILURE);
  }

  pa1 = mxCreateDoubleMatrix(3,3,mxREAL);
  if (pa1 == NULL) {
      printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
      printf("Unable to create mxArray.\n");
      return(EXIT_FAILURE);
  }

  pa2 = mxCreateDoubleMatrix(3,3,mxREAL);
  if (pa2 == NULL) {
      printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
      printf("Unable to create mxArray.\n");
      return(EXIT_FAILURE);
  }
  memcpy((void *)(mxGetPr(pa2)), (void *)data, sizeof(data));

  pa3 = mxCreateString("MATLAB: the language of technical computing");
  if (pa3 == NULL) {
      printf("%s :  Out of memory on line %d\n", __FILE__, __LINE__);
      printf("Unable to create string mxArray.\n");
      return(EXIT_FAILURE);
  }

  status = matPutVariable(pmat, "LocalDouble", pa1);
  if (status != 0) {
      printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
      return(EXIT_FAILURE);
  }

  status = matPutVariableAsGlobal(pmat, "GlobalDouble", pa2);
  if (status != 0) {
      printf("Error using matPutVariableAsGlobal\n");
      return(EXIT_FAILURE);
  }

  status = matPutVariable(pmat, "LocalString", pa3);
  if (status != 0) {
      printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
      return(EXIT_FAILURE);
  }


  // * Ooops! we need to copy data before writing the array.  (Well,
  // * ok, this was really intentional.) This demonstrates that
//   * matPutVariable will overwrite an existing array in a MAT-file.


  memcpy((void *)(mxGetPr(pa1)), (void *)data, sizeof(data));
  status = matPutVariable(pmat, "LocalDouble", pa1);
  if (status != 0) {
      printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
      return(EXIT_FAILURE);
  }

  /* clean up */
/*
  mxDestroyArray(pa1);
  mxDestroyArray(pa2);
  mxDestroyArray(pa3);

  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(EXIT_FAILURE);
  }


//   * Re-open file and verify its contents with matGetVariable


  pmat = matOpen(file, "r");
  if (pmat == NULL) {
    printf("Error reopening file %s\n", file);
    return(EXIT_FAILURE);
  }


  // * Read in each array we just wrote


  pa1 = matGetVariable(pmat, "LocalDouble");
  if (pa1 == NULL) {
    printf("Error reading existing matrix LocalDouble\n");
    return(EXIT_FAILURE);
  }
  if (mxGetNumberOfDimensions(pa1) != 2) {
    printf("Error saving matrix: result does not have two dimensions\n");
    return(EXIT_FAILURE);
  }

  pa2 = matGetVariable(pmat, "GlobalDouble");
  if (pa2 == NULL) {
    printf("Error reading existing matrix GlobalDouble\n");
    return(EXIT_FAILURE);
  }
  if (!(mxIsFromGlobalWS(pa2))) {
    printf("Error saving global matrix: result is not global\n");
    return(EXIT_FAILURE);
  }

  pa3 = matGetVariable(pmat, "LocalString");
  if (pa3 == NULL) {
    printf("Error reading existing matrix LocalString\n");
    return(EXIT_FAILURE);
  }

  status = mxGetString(pa3, str, sizeof(str));
  if(status != 0) {
      printf("Not enough space. String is truncated.");
      return(EXIT_FAILURE);
  }
  if (strcmp(str, "MATLAB: the language of technical computing")) {
    printf("Error saving string: result has incorrect contents\n");
    return(EXIT_FAILURE);
  }

  /* clean up before exit */
/*
  mxDestroyArray(pa1);
  mxDestroyArray(pa2);
  mxDestroyArray(pa3);

  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(EXIT_FAILURE);
  }
  printf("Done\n");
*/


MATFile * pmat; char ** you ; int ndir; int ndir2; mxArray * pa; mxArray * qa; mxArray *ra; mxArray *sa; mxArray *ta; mxArray *ua; mxArray *va; mxArray *xa;  mxArray *ya; mxArray *za; mxArray *zb; mxArray *in; const char * name; char ** dir;
char ** dir2;
pmat = matOpen ("init_net.mat", "r" ) ;
//dir = matGetDir ( pmat, &ndir ) ;
//double * paData;
double * qaData;
double * raData;
double * saData;
double * taData;
double * uaData;
double * xaData;
double * yaData;

double * zaData;
double * zbData;

/*
for ( int i = 0 ; i <ndir; i ++ )
{
   cout <<"i is: " << i << endl;
   name = dir[i];
   printf("%s", name);
   pa = matGetVariable(pmat, name);
   cout << endl;
  // pa = mxGetField ( pmat, name ) ;
   paData = ( double * ) mxGetData ( pa ) ;
   cout <<paData[1]<< endl;

//   cout << paData[0] << endl;

}
*/
//Get positions of nodes
pa = matGetVariable(pmat, "P");
qa = mxGetField(pa, 0, "states");
in = mxGetField(pa, 0, "fixed");
//Get spring variables.
ra = matGetVariable(pmat, "W");
sa = mxGetField(ra, 0, "k1");
ta = mxGetField(ra, 0, "k3");
ua = mxGetField(ra, 0, "d1");
xa = mxGetField(ra, 0, "d3");
ya = mxGetField(ra, 0, "l0");

za = mxGetField(ra, 0, "from");
zb = mxGetField(ra, 0, "to");

in = mxGetField(ra, 0, )

qaData = ( double * ) mxGetData ( qa ) ;

raData = ( double * ) mxGetData ( ra ) ;
saData = ( double * ) mxGetData ( sa ) ;
taData = ( double * ) mxGetData ( ta ) ;
uaData = ( double * ) mxGetData ( ua ) ;
xaData = ( double * ) mxGetData ( xa ) ;
yaData = ( double * ) mxGetData ( ya ) ;

zaData = ( double * ) mxGetData ( za ) ;
zbData = ( double * ) mxGetData ( zb ) ;
//X positions of nodes
vector<double> x_nodes;
vector<double> y_nodes;
//input nodes.
vector<bool> input_nodes;

//Get k1, k3, d1, d3, l0, from, to
vector<double> k1;
vector<double> k3;
vector<double> d1;
vector<double> d3;
vector<double> l0;
vector<double> node1;
vector<double> node2;

for(int i=0; i<78; i++)
{
  if(i<30)
  {
  x_nodes.push_back(qaData[i]);
  cout << x_nodes[i] << endl;
  }
  if(i>=30 && i<60)
  {
  y_nodes.push_back(qaData[i]);
  cout << y_nodes[i] << endl;
  }

  k1.push_back(saData[i]);
  cout << k1[i] << endl;
  k3.push_back(taData[i]);
  cout << k3[i] << endl;
  d1.push_back(uaData[i]);
  cout << d1[i] << endl;
  d3.push_back(xaData[i]);
  cout << d3[i] << endl;
  l0.push_back(yaData[i]);
  cout << l0[i] << endl;;

  node1.push_back(zaData[i]);
  cout << node1[i] << endl;
  node2.push_back(zbData[i]);
  cout << node2[i] << endl;
}
//cout <<"The first element of P and states is: " << paData[29] << endl;


double wash_out_time = 20000;
double learning_time = 200000;
double learning_time_test = 15000;

//

return 0;

//  const mxArray *mxTmp;
  //double  *tmp;

  //mxArray *pa;
  //pa = mxGetField(init_net, 0, "P.states");

//  mxTmp = mxGetField(init_net, 0,  "P.states");
  //tmp   = mxGetPr(mxTmp);

  //mxArray* net = matGetVariable(init_net, "P.states");
  //cout << *mxGetPr(mxGetCell(net, 1));

  //bool isNetStruct = mxIsStruct(net);

  // //reading data
//  MATFile* datafile =  matOpen("D:\\data.mat", "r");
  //mxArray* data = matGetVariable(datafile, "minidata");




    InitialDataValues input_data;
    //send from processor
    srand(rdtsc());

    vector<double> Volterra;
    vector<double> Input;
    vector<string> Input_Lines;
    vector<string> Volterre_Lines;  // Todo: ???

    //init path
    string init_path("");

    cout << "-- Start ---------------------------------------- " << endl;


    // Todo: I think it would be useful to define class for reading in files from csv, maybe best even encapsulated in a class
    // there are plenty of these out there
  //  ifstream file_Input ( "/Users/hh14185/Leverhulme_Trust_Proposal_spider_web/Xcode/Nonlinear-Mass-Spring-System-BR/input.csv" );
    string final_path("Data/inputsignal.csv");
    final_path.insert(0, init_path);
    ifstream file_Input (final_path); file_Input.precision(15);
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

  //  double wash_out_time = 20000;
    //double learning_time = 200000;
    //double learning_time_test = 15000;


    // setting parameters for simulation
    // This should be possible to read in from a text file
    input_data.N = 35;
  //  input_data.N = 17;
    input_data.ux=0;
    input_data.uy= 0;

    input_data.input_connectivity_percentage = 20;
    //data.w_in_initial = -1;
    input_data.min_input_weight = -1;
    input_data.max_input_weight = 1;
    input_data.min_x_position = 0;
    input_data.max_x_position = 10;
    input_data.min_y_position  = 0;
    input_data.max_y_position = 10;

    input_data.min_k3 = 1;
    input_data.max_k3  = 100;
    input_data.min_d3 = 1;
    input_data.max_d3  = 100;

    input_data.min_k1 = 1;
    input_data.max_k1  = 200;
    input_data.min_d1 = 1;
    input_data.max_d1  = 200;


    input_data.dt = 0.001;
    input_data.t0 = wash_out_time*input_data.dt;
    input_data.tmax = (wash_out_time+learning_time+learning_time_test)*input_data.dt;

    vector<double> Sine_Wave;




    vector<int> no_of_springs;
    vector<double> MSE_list;
    vector<double> times;

    ofstream MSE_list_out("MSE_list.csv", ofstream::out);
    ofstream no_of_springs_out("no_of_springs.csv", ofstream::out);
    ofstream times_out("time.csv", ofstream::out);





  //  Simulation sim(data, Volterra, Input, wash_out_time, learning_time, learning_time_test);

  for(int i=0; i<1; i++)
  {
    auto begin = std::chrono::high_resolution_clock::now();
   Simulation sim(input_data, Input, Volterra, wash_out_time, learning_time, learning_time_test);
    cout <<"The number of nodes is: " << input_data.N << endl;
    cout <<"The number of springs is: " << sim.Spring_List() << endl;



   auto end = std::chrono::high_resolution_clock::now();
   cout << "The time it took for the programme to run in total in milliseconds: ";
   std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms";

   //For the experiment matlab
   no_of_springs.push_back(sim.Spring_List());
   MSE_list.push_back(sim.return_MSE());
   times.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count());

   MSE_list_out << sim.return_MSE();
   no_of_springs_out << sim.Spring_List();
   times_out << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
   if(i<20) MSE_list_out <<',';
   if(i<20) no_of_springs_out <<',';
   if(i<20) times_out <<',';

 }
 cout << endl;
 cout << "The mean of the MSE's is: " << accumulate( MSE_list.begin(), MSE_list.end(), 0.0)/MSE_list.size() << endl;
 cout << "The mean no of springs is: " << accumulate( no_of_springs.begin(), no_of_springs.end(), 0.0)/no_of_springs.size() << endl;
 double mean = accumulate( no_of_springs.begin(), no_of_springs.end(), 0.0)/no_of_springs.size();
 double sq_sum = std::inner_product(no_of_springs.begin(), no_of_springs.end(), no_of_springs.begin(), 0.0);
 double stdev = std::sqrt(sq_sum / no_of_springs.size() - mean * mean);
 cout <<"The stdev of the springs is: "<< stdev << endl;
 cout << "The mean time taken is: " << ((accumulate( times.begin(), times.end(), 0.0)/times.size())/1000) << endl;







   return 0;
}
