#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>
#include "Simulation.h"
#include "Eigen/Dense"
#include "Eigen/QR"
#include <sstream>
#include <fstream>
#include <string>
using namespace std;
using namespace Eigen;


Simulation::Simulation(InitialDataValues &data, vector<double> &TS)
{

  //rand(); rand(); rand();
  this->N = data.N;
  //Step 2: The positions of the nodes were initialized and 20% of the nodes are connected to the input.
  this->input_connectivity = data.input_connectivity;
  this->total_input_nodes = (data.input_connectivity)*N;
  this->w_out_initial = data.w_out_initial;
  this->w_out_final = data.w_out_final;
  this->w_in_initial = data.w_in_initial;
  this->w_in_final = data.w_in_final;
  this->initial_log_uniform= data.initial_log_uniform;
  this->final_log_uniform = data.final_log_uniform;
  this->initial_uniform = data.initial_uniform;
  this->final_uniform = data.final_uniform;
  this->t0 = data.t0;
  this->tmax = data.tmax;
  this->dt = data.dt;

  this->range0x = data.range0x;
  this->range0y = data.range0y;
  this->range1x = data.range1x;
  this->range1y = data.range1y;

  this->ux = data.ux;
  this->uy = data.uy;

  this->maxtimesteps = ((data.tmax - data.t0)/data.dt);


  Target_Signal = TS;

  bool learning_phase = 0;
  //Learning phase
  Initialize_Nodes(range0x, range1x, range0y, range1y);
  Delaunay_Triangulation_and_Spring_Creation();
}

Simulation::Simulation(double radius, int rounds, int no_of_points_per_round, InitialDataValues &data, vector<double> &Lvx)
{
  this->input_connectivity = data.input_connectivity;
  this->total_input_nodes = (data.input_connectivity)*rounds*no_of_points_per_round;

  this->N = data.N;
  //Step 2: The positions of the nodes were initialized and 20% of the nodes are connected to the input.
  this->w_out_initial = data.w_out_initial;
  this->w_out_final = data.w_out_final;
  this->w_in_initial = data.w_in_initial;
  this->w_in_final = data.w_in_final;

  this->x = x;
  this->y = 1-x;

  //Learning phase
  Initialize_Nodes(radius, rounds, no_of_points_per_round, data);

  this->t0 = data.t0;
  this->tmax = data.tmax;
  this->dt = data.dt;
  this->maxtimesteps = ((data.tmax - data.t0)/data.dt);

  Target_Signal = Lvx;

  Execute_In_Time();
  Output_For_Plot();
}

void Simulation::Initialize_Nodes(double range0x, double range1x, double range0y, double range1y)
   {
    ofstream fixed("fixednode.txt");
 		ofstream Initialnodes("initial.txt");


 		double x;
 		double y;

 		double x1 =range0x;
 		double x0 =range1x;

     //for fixed nodes.
 		int j=0;
 		int k=0;

      for(int i=0; i<N; i++)
 		 {

 			 x=Uniform(range0x, range1x);
 			 if(x1<x)
 			 {
 			  x1=x;
 				j=i;
 			}
 			 y=Uniform(range0y, range1y);
 			 if(x0>x)
 			 {
 				 x0=x;
 				 k=i;
 			 }
 		   Nodes p(x, y, 0);
 			 Initialnodes <<x <<"," <<y << endl;
      //The first input_connectivity percent of the nodes are marked as input nodes, and the
        n.push_back(p);
 	   }

 		 fixed <<j << endl;
 		 fixed <<k << endl;

 		 n[j].FixedNode();
 		 n[k].FixedNode();
 		//Just one node for test;
 		 //Fixed the leftmost and rightmost nodes.
   }

void Simulation::Initialize_Nodes(double radius, int rounds, int no_of_points_per_round, InitialDataValues &data)
{
  double angle = ((2*M_PI)/no_of_points_per_round);
  double x_position;
  double y_position;

  double xprevious = 0;
  double yprevious = 0;

  double xnext =0;
  double ynext = 0;

  double k1;
  double d1;
  double k3;
  double d3;

  double x0;
  double x1;

  double y0;
  double y1;

  double l0;
  double wout;
  double win;
  double BeforeRand = 0;


  int k =0;

 for(int j=0; j<rounds; j++)
 {
   x0 = (j+1)*radius*cos((0));
   y0 = (j+1)*radius*sin((0));


   for(int i=0; i<no_of_points_per_round; i++)
   {
      //So this
      x_position = (j+1)*radius*cos((i*angle));
      y_position = (j+1)*radius*sin((i*angle));


      Nodes node(x_position, y_position, 0);




      //THIS IS NOT GOOD CODING PRACTICE. FIX IT.
      //FIX IT
      //FIX IT
      win = Uniform(data.w_in_initial, data.w_in_final);
      cout << win << endl;

       /*
      win = Uniform(0,2);
      win -= 1;
      BeforeRand = Uniform(0,1);
      */

      n.push_back(node);

      if(BeforeRand<=input_connectivity)
      {
      n[k].Input_Node(data.ux, data.uy, win);
      cout << endl;
      }

      if(i>0)
      {
      k1 = log10(Uniform(data.initial_log_uniform, data.final_log_uniform));
      d1 = log10(Uniform(data.initial_log_uniform, data.final_log_uniform));
      k3 = Uniform(data.initial_uniform, data.final_uniform);
      d3 = Uniform(data.initial_uniform, data.final_uniform);

      l0 = Eucl_Dist(n[k].X_Position(), n[k].Y_Position(), n[k-1].X_Position(), n[k-1].Y_Position());
      cout << "lo is: " << l0 << endl;
      wout = Uniform(data.w_out_initial, data.w_out_final);

      s.push_back(Springs(k1, d1, k3, d3, l0, k, k-1, wout));
      }

      //For radial pattern.
      if(j>0)
      {
      k1 = log10(Uniform(data.initial_log_uniform, data.final_log_uniform));
      d1 = log10(Uniform(data.initial_log_uniform, data.final_log_uniform));
      k3 = Uniform(data.initial_uniform, data.final_uniform);
      d3 = Uniform(data.initial_uniform, data.final_uniform);

      l0 = Eucl_Dist(n[k].X_Position(), n[k].Y_Position(), n[ k-no_of_points_per_round].X_Position(), n[ k-no_of_points_per_round].Y_Position());
         cout << "lo is: " << l0 << endl;
      wout = Uniform(data.w_out_initial, data.w_out_final);

      s.push_back(Springs(k1, d1, k3, d3, l0, k, k-no_of_points_per_round, wout));
      }

      x0 = x_position;
      y0 = y_position;
      k++;

   }
   x_position = (j+1)*radius*cos((0));
   y_position = (j+1)*radius*sin((0));

   k1 = log10(Uniform(data.initial_log_uniform, data.final_log_uniform));
   d1 = log10(Uniform(data.initial_log_uniform, data.final_log_uniform));
   k3 = Uniform(data.initial_uniform, data.final_uniform);
   d3 = Uniform(data.initial_uniform, data.final_uniform);
   l0 = Eucl_Dist(n[ k-no_of_points_per_round].X_Position(), n[ k-no_of_points_per_round].Y_Position(),n[k-1].X_Position(), n[k-1].Y_Position());
   cout << "lo is: " << l0 << endl;
   wout = Uniform(data.w_out_initial, data.w_out_final);
   cout << "wout is: " << wout << endl;

   if(no_of_points_per_round>2)
   {
   s.push_back(Springs(k1, d1, k3, d3, l0, (k-no_of_points_per_round), k-1, wout));
   }
 }



n[no_of_points_per_round*(rounds-1)].FixedNode();
n[(no_of_points_per_round/2)+no_of_points_per_round*(rounds-1)].FixedNode();

}

void Simulation::FixNode(int i)
{
  n[i].FixedNode();
}

void Simulation::Delaunay_Triangulation_and_Spring_Creation()
{
    DelaunayTriangulation DT(abs(range1x-range0x), abs(range1y-range0y));
    double win = 0;
    double BeforeRand = 0;

    //This offset is to counteract the negativity.
  //  double offset = abs(w_in_initial)+abs(w_in_final);

  //  cout <<"offset is: " << offset;
  //  cout << endl << w_in_initial;
  //  cout <<endl << w_in_final;

    for(int i=0; i<N; i++)
    {
      win = Uniform(w_in_initial, w_in_final);
  //    cout <<"win is: "<< win << endl;
  //    win -= offset;
      cout <<"win is: "<< win << endl;
      if(BeforeRand<=input_connectivity) n[i].Input_Node(ux, uy, win);
      DT.AddPoint(Point(n[i].X_Position(),n[i].Y_Position(),0));
    //  DT.AddPoint(Point(n[i].X_Position(), Point(n[i].Y_Position());
    }
    DT.print();

    Get_Triangles(DT);
    Initialize_Springs();
    Execute_In_Time();
    Output_For_Plot();

}

void Simulation::Execute_In_Time()
{
  double Fsum =0;
  double Fx_nodea =0;
  double Fy_nodea =0;
  double Fx_nodeb =0;
  double Fy_nodeb =0;
  double Fz_nodea =0;
  double Fz_nodeb =0;
  double theta = 0;
  double l = 0;

  double nodea =0;
  double nodeb = 0;

  double x0 = 0;
  double x1 = 0;
  double y0 = 0;
  double y1 = 0;

  //z coordinates
  double z0 = 0;
  double z1 = 0;


  double lx = 0;
  double ly = 0;
  double lz = 0;

  //Direction cosine triangles

  double alpha = 0;
  double beta = 0;
  double gamma =0;




  ofstream ofs("Node1.csv");
  ofstream ofs2("Node2.csv");
  ofstream ofs3("SampleForce.csv");
  ofstream bad("Badcoefficients.csv");

  double currentlength = 0;

  int y = maxtimesteps - x;





  MatrixXd LearningMatrix(Target_Signal.size(), s.size());
  MatrixXd TargetSignal(Target_Signal.size(), s.size());

  int p = 0.66666*Target_Signal.size();

  MatrixXd LearningMatrix2(p, s.size());
  MatrixXd TargetSignal2(p, s.size());




  double outputsignal = 0;
  for(int j=0; j<s.size(); j++)
  {
    LearningMatrix(0,j)=s[j].Return_Original_Length();
    LearningMatrix2(0, j) = s[j].Return_Original_Length();
    TargetSignal(0,j) = Target_Signal[0];
    TargetSignal2(0,j) = Target_Signal[0];
  }


  outputsignal = 0;
  double currenttime = 0;
  double currenttime1 = 0;
  double currenttime2 = 0;



  for(int i=1; i<Target_Signal.size(); i++)
  {
    currenttime = t0+ i*dt;
  //  Learning_Matrix_3.push_back(vector<double>());

  for(int j=0;  j<s.size(); j++)
    {
       s[j].ForceEq(Fsum);
       ofs3 << i*dt<<"," <<Fsum << endl;

       nodea = s[j].Nodea();
       nodeb = s[j].Nodeb();

       x0 = n[nodea].X_Position();
       x1 = n[nodeb].X_Position();
       y0 = n[nodea].Y_Position();
       y1 = n[nodeb].Y_Position();
       z0 = n[nodea].Z_Position();
       z1 = n[nodeb].Z_Position();


       //Change position of first node
       //theta = abs(Angle(x0, x1, y0, y1));

       //Lengths in x,y and z directions
       l = Eucl_Dist(x0, y0, z0, x1, y1, z1);



       cout <<"There shuold be oscillations in x: " << x0 << endl;
       cout <<"There shuold be oscillations in x: " << x1 << endl;

       //Unit vector for two points

       lx = x1 - x0;
       ly = y1 - y0;
       lz = z1 - z0;

      //Direction cosines
       alpha = lx/l;
       beta = ly/l;
       gamma = lz/l;



       cout<<"Fsum is: "<< Fsum << endl;
       cout <<"Fx is: " << Fsum*alpha << endl;
       cout <<"Fy is: " << Fsum*beta << endl;
       cout <<"Fz is: " << Fsum*gamma << endl;

       cout <<"alpha is: " << alpha << endl;
       cout <<"beta is: " << beta << endl;
       cout <<"gamma is: " << gamma << endl;

       Fx_nodeb = Fsum*alpha;
       Fy_nodeb = Fsum*beta;
       Fz_nodeb = Fsum*gamma;

       Fx_nodea = -Fsum*alpha;
       Fy_nodea = -Fsum*beta;
       Fz_nodea = -Fsum*gamma;


       cout <<Fx_nodea << endl;
       cout <<Fx_nodeb << endl;

       cout <<Fy_nodea << endl;
       cout <<Fy_nodeb << endl;

       cout <<Fz_nodea << endl;
       cout <<Fz_nodeb << endl;


       cout <<"Alpha is: " << alpha << endl;
       cout <<"Beta is: " << beta << endl;
       cout <<"Gamma is: " << gamma << endl;


      // Fz_nodea += Uniform(-0.01,0.01)*sin(currenttime);
    //   Fz_nodeb += Uniform(-0.01,0.01)*sin(currenttime);


       //Square web temp
       //if(j<=4) Fz_nodea += Uniform(-0.01,0.01);
       //if(j<=4) Fz_nodeb += Uniform(-0.01,0.01);

       //Square wave on first web

       //if(j==0 && currenttime == 0.001) Fz_nodea += 1;
       //if(j==0 && currenttime == 0.001) Fz_nodeb += 1;

       //if(j==0 && currenttime == 0.004) Fz_nodea -= 1;
      // if(j==0 && currenttime == 0.004) Fz_nodeb -= 1;

       //if(currenttime<currenttime2 && currenttime>currenttime1) n[nodea].Change_Z_Position(1*(currenttime)*currenttime + 1*(currenttime) + 0.5);
      // if(currenttime<currenttime2 && currenttime>currenttime1) n[nodeb].Change_Z_Position(1*(currenttime)*currenttime + 1*(currenttime) + 0.


       n[nodea].Change_Position(Fx_nodea, Fy_nodea, Fz_nodea, dt);
       n[nodeb].Change_Position(Fx_nodeb, Fy_nodeb, Fz_nodeb, dt);





       x0 = n[nodea].X_Position();
       x1 = n[nodeb].X_Position();

       ofs <<dt*i <<","<< x0;
       ofs2 <<dt*i <<"," << x1;


       y0 = n[nodea].Y_Position();
       y1 = n[nodeb].Y_Position();

       ofs <<"," << y0 << endl;
       ofs2 <<"," << y1 << endl;

       //Be very careful with the lengths here.
       l = Eucl_Dist(x0, y0, z0, x1,y1,z1);
       currentlength = l;

       cout << l << endl;


      // cout <<"Is this running?" << endl;
       LearningMatrix(i,j) = currentlength;
       if(i<p) LearningMatrix2(i,j) = currentlength;
  //     if(i>=p) Learning_Matrix_3.push_back(currentlength);

       TargetSignal(i,j) = Target_Signal[i];
       if(i<p) TargetSignal2(i,j) = Target_Signal[i];



       s[j].Change_Length_And_Velocity(dt, l);
       Fsum = 0;
      }
      outputsignal = 0;
    }

   cout<< Target_Signal.size() << endl;
   cout << p << endl;
    cout <<"The size of Learning Matrix" << LearningMatrix.rows() << endl;

    cout <<"The size of Learning Matrix 2" << LearningMatrix2.rows() << endl;
  //  cout <<"The size of Learning Matrix 3" << Learning_Matrix_3.size() << endl;


//   cout <<LearningMatrix.rows();
  // cout << endl;
  // cout <<  LearningMatrix.block(0, 0, x, LearningMatrix.cols()).rows();
  // cout << endl;

    LM = LearningMatrix;
    //This is the protocl for 2/3 learning weights, 1/3 learnt signal.
  //  LM = LearningMatrix1;

  //LearningMatrix = LearningMatrix.block(0, 0, x, LearningMatrix.cols());

  //    cout <<"TS is now " << T
    Moore_Penrose_Pseudoinverse(LearningMatrix2);
  //  cout << LearningMatrix;
  //  TargetSignal = TargetSignal.block(0,0,x,TargetSignal.cols());
    LearningMatrix2= LearningMatrix2 * TargetSignal2;

//    cout << LearningMatrix;
//    cout << TargetSignal;
    Populate_Learning_Weights(LearningMatrix2);

    cout << Learning_Weights[0] << endl;
    cout << Learning_Weights[1] << endl;
    cout << Learning_Weights[2] << endl;
    cout << Learning_Weights[3] << endl;
    cout << Learning_Weights[4] << endl;

}

void Simulation::Moore_Penrose_Pseudoinverse(MatrixXd& L)
{
  //First one third of signal.
  L = L.completeOrthogonalDecomposition().pseudoInverse();
}

vector<double>& Simulation::Return_Learning_Weights()
{
  return Learning_Weights;
}

MatrixXd& Simulation::Return_Learning_Matrix()
{
  return LM;
}

void Simulation::Populate_Learning_Weights(MatrixXd& L)
{
  for(int j=0; j<s.size(); j++)
  {
    Learning_Weights.push_back(L(j,0));
  }
}

void Simulation::Output_Signal_And_MSE()
{
//  MatrixXd LM;
//  vector<double> LW;

//  LM = Return_Learning_Matrix();
//  LW = Return_Learning_Weights();



  ofstream output("outputsignal.csv");
  ofstream output2("learningweights.csv");
  ofstream output3("targetsignal.csv");

  double outputsignal = 0;
  double currenttime = 0;

  int p = 0.66666*Target_Signal.size();
  cout << p << endl;
  cout << Target_Signal.size() << endl;

  for(int i=p; i<Target_Signal.size(); i++)
  {

  for(int j=0; j<s.size(); j++)
  {
    outputsignal += Learning_Weights[j] * LM(i,j);
    if(i==p) output2 << Learning_Weights[j] << endl;
  }

  Output_Signal.push_back(outputsignal);
  currenttime = t0 + i*dt;

  output << currenttime <<"," << outputsignal;
  output << endl;
  output3 <<currenttime <<"," << Target_Signal.at(i);
  output3 << endl;
  outputsignal = 0;
  }


cout << Output_Signal.size() << endl;
vector<double> Target_Signal_23(Target_Signal.begin() + p, Target_Signal.end());
cout << Target_Signal.size();
cout << endl;
cout <<"The mean squared error of the output signal versus the target signal is: " << MSE(Output_Signal, Target_Signal_23);
//cout <<endl;
}

void Simulation::Output_Signal_And_MSE(vector<double>& External_Weights)
{
//  MatrixXd LM;
//  vector<double> LW;

//  LM = Return_Learning_Matrix();
//  LW = Return_Learning_Weights();

  ofstream output("outputsignal.csv");
  ofstream output2("learningweights.csv");
  ofstream output3("targetsignal.csv");

  double outputsignal = 0;
  double currenttime = 0;

  int finaltime = 0.66666*Target_Signal.size();


  for(int i=0; i<maxtimesteps; i++)
  {
  for(int j=0; j<External_Weights.size(); j++)
  {
    outputsignal += External_Weights[j] * LM(i, j);
    if(i==0) output2 << External_Weights[j] << endl;
  }

  Output_Signal.push_back(outputsignal);
  currenttime = t0 + i*dt;

  output << currenttime <<"," << Output_Signal.at(i);
  output << endl;
  output3 <<currenttime <<"," << Target_Signal.at(i);
  output3 << endl;
  outputsignal = 0;
  }

  cout <<"The mean squared error of the output signal versus the target signal is: " << MSE(Output_Signal, Target_Signal);
  cout <<endl;
}

double Simulation::MSE(vector<double>& A, vector<double>& Ahat)
{

double MSEsum = 0;
double MSE = 0;
double Total;

for(int i =0; i<A.size(); i++)
{
  MSEsum += (A[i]-Ahat[i])*(A[i]-Ahat[i]);
}

Total = (double)A.size();

//cout <<"Inverse total is: " << Total << endl;
MSE = (1/Total)*MSEsum;

return MSE;
}

int Simulation::Random_Input_Nodes(int N)
{
  return N*((int) rand() / (RAND_MAX));
}

double Simulation::Uniform(double M, double N)
{
  return M + (rand() / ( RAND_MAX / (N-M) ) ) ;
}

double Simulation::Log_10_Uniform(double initial, double finalvalue)
{
  return exp(Uniform(initial, finalvalue)/(2.302585093));
}

double Simulation::Spring_And_Damping_Coefficient_1(double initial, double finalvalue)
{
  return Log_10_Uniform(initial, finalvalue);
}

double Simulation::Spring_And_Damping_Coefficient_2(double initial, double finalvalue)
{
  return Uniform(initial, finalvalue);
}

double Simulation::Eucl_Dist(double x1, double y1, double x2, double y2)
{
  return sqrt((y2-y1)*(y2-y1) + (x2-x1)*(x2-x1));
}

double Simulation::Eucl_Dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
  return sqrt((y2-y1)*(y2-y1) + (x2-x1)*(x2-x1) + (z2-z1)*(z2-z1));
}

double Simulation::Angle(double x0, double x1, double y0, double y1)
{
  return atan((y1-y0)/(x1-x0));
}

double Simulation::X_Comp(double vectorsum, double theta)
{
  return vectorsum*cos(theta);
}

double Simulation::Y_Comp(double vectorsum, double theta)
{
  return vectorsum*sin(theta);
}

void Simulation::Sort(int &a, int &b, int &c)
{
    if(a>b)
  {
        int tmp = a;
        a = b;
        b = tmp;
    }
    if(a>c)
  {
        int tmp = a;
        a=c;
        c = tmp;
    }
    if(b>c)
  {
        int tmp = b;
        b=c;
        c=tmp;
    }
}

void Simulation::Sort(int &a, int &b)
{
    if(a>b)
  {
        int tmp = a;
        a = b;
        b = tmp;
    }
}

void Simulation::RemoveDuplicates(vector<vector<double>> &x)
{
   sort (x.begin(), x.end());
   int i=1;
   while(i<x.size())
   {
    if(x[i][1] ==x[i-1][1] && x[i][0] ==x[i-1][0])
    {
      x.erase(x.begin()+(i-1));
      i-=1;
    }
      i++;
   }
}

void Simulation::Create_EdgeNodeList()
{
  for(int i=0; i<EdgeList.size(); i++)
  {
    EdgeNodeList.push_back(s[i].Nodea());
    EdgeNodeList.push_back(s[i].Nodeb());
  }
  EdgeNodeList.erase( unique( EdgeNodeList.begin(), EdgeNodeList.end() ), EdgeNodeList.end() );
}

double Simulation::Output_No_of_Edges()
{
  return EdgeList.size();
}

void Simulation::Output_Spring_And_Node_Positions()
{
  for (int i =0; i<s.size(); i++)
  {
    s[i].Output();
  }
}

void Simulation::Initialize_Springs()
{
  cout <<"The number of edges for: " <<N << " mass points is: " << EdgeList.size() << endl;
  //cout <<"The number of unduplicated edges should be " << EdgeList.size() << endl;
  //Spring and damping coefficients

  double k1 = 0;
  double d1 = 0;
  double k3 = 0;
  double d3 = 0;
  double l0 = 0;

  double x0;
  double y0;
  double x1;
  double y1;
  double wout;

  int arraysubscript1=0;
  int arraysubscript2=0;

  ofstream ofs3("k1output.csv");
  ofstream ofs4("d1output.csv");
  ofstream ofs5("k3output.csv");
  ofstream ofs6("d3output.csv");


  for(int i=0; i<EdgeList.size(); i++)
  {
      //These take the arraysubscripts and disregard the first four points

      arraysubscript1 = EdgeList[i].at(0) - 4;
      arraysubscript2 = EdgeList[i].at(1) - 4;

      x0 = n[arraysubscript1].X_Position();
      x1 = n[arraysubscript2].X_Position();
      y0 = n[arraysubscript1].Y_Position();
      y1 = n[arraysubscript2].Y_Position();

      //These spring and damping coefficients are not giving different values
      k1 = log10(Uniform(initial_log_uniform, final_log_uniform));
      d1 = log10(Uniform(initial_log_uniform, final_log_uniform));
      k3 = Uniform(initial_uniform, final_uniform);
      d3 = Uniform(initial_uniform, final_uniform);

      ofs3 <<  k1 << endl;
      ofs4 <<  d1 << endl;
      ofs5 <<  k3 << endl;
      ofs6 <<  d3  << endl;

      l0 = Eucl_Dist(x0, y0, x1, y1);
      wout = Uniform(w_out_initial, w_out_final);

      s.push_back(Springs(k1, d1, k3, d3, l0, arraysubscript1, arraysubscript2, wout));
      s[i].Output();
   }
}

void Simulation::Get_Triangles(DelaunayTriangulation &Delaunay)
{
  stringstream s1;

  vector<string> tri;
  string es;
  int k =0;

  int node1;
  int node2;
  int node3;
  char sep=',';

  double vertices = 0;
  double noofconnectingedges = 0;


  //	This takes the nodes of each triangle and parses them so only the relevant nodes a
    for(auto e: Delaunay.triangles)
    {
      //cout <<"Node1:" <<get<0>(e) <<" " << "Node2:"<<" " <<get<1>(e)<< endl;
      s1 << e;
    //	cout << e << endl;
      tri.push_back(s1.str());
      tri.at(k) = tri.at(k).substr(1, tri.at(k).size()-3);
      //cout << tri.at(k) << endl;
  //		cout <<endl;
          s1.str("");
          istringstream iss(tri.at(k));

          iss >> node1;
        //  cout <<"node1: " <<node1 <<endl;
          iss >>sep;
          iss >>node2;
        //  cout <<"node2: "<<node2 << endl;
          iss >>sep;
          iss >>node3;
      //    cout <<"node3: " << node3 << endl;
          iss >>sep;

          if(node1!=0 && node1!=1 && node1!=2 && node1!=3) vertices++;
          if(node2!=0 && node2!=1 && node2!=2 && node2!=3) vertices++;
          if(node3!=0 && node3!=1 && node3!=2 && node3!=3) vertices++;

          if(vertices ==2)
          {
      noofconnectingedges++;
    //	cout << "Connect two nodes!: " <<node1 <<" " <<node2 <<" "<< node3 <<" ";
      Sort(node1, node2, node3);
      if(node1!=0 && node1!=1 && node1!=2 && node1!=3)
      {
        NodeList.push_back(node1);
        NodeList.push_back(node2);
        EdgeList.push_back(NodeList);
        NodeList.clear();
        }
      else if(node2!=0 && node2!=1 && node2!=2 && node2!=3)
      {
        NodeList.push_back(node2);
        NodeList.push_back(node3);
        EdgeList.push_back(NodeList);
        NodeList.clear();
      }

      //Add a spring
    //	Springs = new Spring()
        }

        if(vertices ==3)
          {
      noofconnectingedges+=2;
      //cout << "Connect all nodes: " <<node1 <<" " <<node2 <<" "<< node3 <<" ";
      Sort(node1, node2, node3);
      if(node1!=0 && node1!=1 && node1!=2 && node1!=3)
      {
        NodeList.push_back(node1);
        NodeList.push_back(node2);
        NodeList.push_back(node2);
        NodeList.push_back(node3);
        EdgeList.push_back(NodeList);
        NodeList.clear();
      }
      else if(node2!=0 && node2!=1 && node2!=2 && node2!=3)
      {
        NodeList.push_back(node2);
        NodeList.push_back(node3);
        EdgeList.push_back(NodeList);
        NodeList.clear();
      }
      //Add a spring
        }
          vertices = 0;
      k++;
   }

      RemoveDuplicates(EdgeList);
 //Remove Duplicates from EdgeNodeList
}

double Simulation::Sine_Wave(double currenttime)
{
  return sin(currenttime);
}

void Simulation::Output_For_Plot()
{
  string str;
  str = "s.csv";
 // if(i%10 == 0) str.erase(str.length()-4);
  ofstream EdgesS(str);
  string str2;
  str2 = "t.csv";
  ofstream EdgesT(str2);
  vector<int>::iterator NodeNums;

  string str3;

  //In Matlab, it does not accept indices of 0 for node graphs.
  //Replace EdgeList with s for consistency.
  for(int j=0; j<s.size()-1; j++)
  {
    EdgesS <<s[j].Nodea()+1 <<",";
    EdgesT <<s[j].Nodeb()+1 <<",";
  }

  EdgesS <<s[s.size()-1].Nodea()+1;
  EdgesT <<s[s.size()-1].Nodeb()+1;


  for(int i=0; i<maxtimesteps; i++)
 {
   str = to_string(i*dt);
  // if(i%10 == 0) str.erase(str.length()-4);
   str.erase(str.length()-3);
   str.append("X.csv");
   ofstream nodesX(str);

   str2 = to_string(i*dt);
  // if(i%10 == 0) str.erase(str.length()-4);
   str2.erase(str.length()-5);
   str2.append("Y.csv");
   ofstream nodesY(str2);

   str3 = to_string(i*dt);
   str3.erase(str.length()-5);
   str3.append("Z.csv");
   ofstream nodesZ(str3);

  int j=0;
  while(j<n.size())
  {
    nodesX << n[j].X_Position();
    if(j<n.size()-1) nodesX<<",";
    nodesY << n[j].Y_Position();
    if(j<n.size()-1) nodesY<<",";
    nodesZ << n[j].Z_Position();
    if(j<n.size()-1) nodesZ<<",";
    j++;
  }

 }

}
Springs Simulation::Spring_Return(int i)
{
  return s[i];
}

Nodes Simulation::NodeReturn(int i)
{
  return n[i];
}

double Simulation::Spring_List()
{
  return EdgeList.size();
}

//The helper functions for the Dynamical Systems class

DynamicalSystems::DynamicalSystems(double t0, double tmax, double dt)
{
  this->maxtimesteps = (int)((tmax - t0)/dt);
  this->tmax = tmax;
  this->dt = dt;
  this->t0 = t0;
}

void DynamicalSystems::SineWave(vector<double> &Sine_Wave)
{
  double currenttime =0;
  for(int i =0; i<maxtimesteps; i++)
  {
     Sine_Wave.push_back(sin(t0+i*dt));
     cout <<sin(t0+i*dt) << endl;
  }
}

void DynamicalSystems::LotkaVolterra(vector<double> &LVx, vector<double> &LVy)
{
  //USe Euler's to get quick Lotka Volterra. Parameters = 1, 1 just for speed.
  double x0;
  double xnext;
  double y0;
  double ynext;

  //Initial conditions
  x0 = 1.05;
  y0 = 1.05;

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
