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


Simulation::Simulation(InitialDataValues &data, vector<double> &IS, vector<double> &TS)
{

  //rand(); rand(); rand();
  this->N = data.N;
  //Step 2: The positions of the nodes were initialized and 20% of the nodes are connected to the input.
  this->input_connectivity_percentage = data.input_connectivity_percentage;
  this->total_input_nodes = (data.input_connectivity_percentage)*N;

  this->input_weight_smallest_value = data.input_weight_smallest_value;
  this->input_weight_largest_value = data.input_weight_largest_value;

  this->log_uniform_smallest_value=  data.log_uniform_smallest_value;
  this->log_uniform_largest_value = data.log_uniform_largest_value;

  this->uniform_smallest_value= data.uniform_smallest_value;
  this->uniform_largest_value = data.uniform_largest_value;


  this->t0 = data.t0;
  this->tmax = data.tmax;
  this->dt = data.dt;

  this->smallest_x_position = data.smallest_x_position;
  this->smallest_y_position = data.smallest_y_position;

  this->largest_x_position = data.largest_x_position;
  this->largest_y_position = data.largest_y_position;

  this->ux = data.ux;
  this->uy = data.uy;

  this->maxtimesteps = ((data.tmax - data.t0)/data.dt);

  Target_Signal = TS;
  Input_Signal = IS;

//  bool learning_phase = 0;
  //Learning phase
  Initialize_Nodes(smallest_x_position, largest_x_position, smallest_y_position, largest_y_position);
  Delaunay_Triangulation_and_Spring_Creation();
}



Simulation::Simulation(double radius, int rounds, int no_of_points_per_round, InitialDataValues &data, vector<double> &Lvx, vector<double> &Lvy)
{
  this->input_connectivity_percentage = data.input_connectivity_percentage;
  this->total_input_nodes = (data.input_connectivity_percentage)*rounds*no_of_points_per_round;

  this->N = data.N;
  //Step 2: The positions of the nodes were initialized and 20% of the nodes are connected to the input.
  this->input_weight_smallest_value = data.input_weight_smallest_value;
  this->input_weight_largest_value = data.input_weight_largest_value;

  this->smallest_x_position = data.smallest_x_position;     // smallest_x_position ?? name is not very descriptive (same for the others below)
  this->largest_x_position = data.largest_x_position;
  this->smallest_y_position=data.smallest_y_position;
  this->largest_y_position = data.largest_y_position;

  //Learning phase
  Initialize_Nodes(radius, rounds, no_of_points_per_round, data);

  this->t0 = data.t0;
  this->tmax = data.tmax;
  this->dt = data.dt;
  this->maxtimesteps = ((data.tmax - data.t0)/data.dt);

  Target_Signal = Lvx;

  execute();
//  Output_For_Plot();
}

void Simulation::Initialize_Nodes(double smallest_x_position, double largest_x_position, double smallest_y_position, double largest_y_position)
   {
    // todo: still needed for debugging
    ofstream fixed("fixednode.txt");
 		ofstream Initialnodes("initial.txt");


 		double x;
 		double y;

 		double x1 =smallest_x_position;
 		double x0 =largest_x_position;


     //for fixed nodes.
 		int j=0;
 		int k=0;

      for(int i=0; i<N; i++)
 		 {

 			 x=Uniform(smallest_x_position, largest_x_position);
 			 if(x1<x)
 			 {
 			  x1=x;
 				j=i;
 			}
 			 y=Uniform(smallest_y_position, largest_y_position);
 			 if(x0>x)
 			 {
 				 x0=x;
 				 k=i;
 			 }
 		   Nodes p(x, y);
 			 Initialnodes <<x <<"," <<y << endl;
      //The first input_connectivity percent of the nodes are marked as input nodes, and the
        n.push_back(p);
 	   }

 		 fixed <<j << endl;
 		 fixed <<k << endl;

 		 n[j].set_Fixed_Node();
 		 n[k].set_Fixed_Node();
 		//Just one node for test;
 		 //Fixed the leftmost and rightmost nodes.
   }

void Simulation::Initialize_Nodes(double radius, int rounds, int no_of_points_per_round, InitialDataValues &data)
{
  double angle = ((2*M_PI)/no_of_points_per_round);
  double x_position;
  double y_position;

  double k1;
  double d1;
  double k3;
  double d3;

  double x0;
//  double x1;

  double y0;
//  double y1;

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


      Nodes node(x_position, y_position);

      //THIS IS NOT GOOD CODING PRACTICE. FIX IT.
      //FIX IT
      //FIX IT
      win = Uniform(data.input_weight_smallest_value, data.input_weight_largest_value);
      cout << win << endl;

       /*
      win = Uniform(0,2);
      win -= 1;
      BeforeRand = Uniform(0,1);
      */

      n.push_back(node);

      if(BeforeRand<=input_connectivity_percentage)
      {
      n[k].init_Input_Node(data.ux, data.uy, win);
      cout << endl;
      }

      if(i>0)
      {
      k1 = log10(Uniform(data.log_uniform_smallest_value, data.log_uniform_largest_value));
      d1 = log10(Uniform(data.log_uniform_smallest_value, data.log_uniform_largest_value));
      k3 = Uniform(data.uniform_smallest_value, data.uniform_largest_value);
      d3 = Uniform(data.uniform_smallest_value, data.uniform_largest_value);

      l0 = Eucl_Dist(x0, y0, x_position, y_position);
    //  wout = Uniform(data.w_out_initial, data.w_out_final);
      wout = 0;

      s.push_back(Springs(k1, d1, k3, d3, l0, k, k-1, wout));
      }

      //For radial pattern.
      if(j>0)
      {
      k1 = log10(Uniform(data.log_uniform_smallest_value, data.log_uniform_largest_value));
      d1 = log10(Uniform(data.log_uniform_smallest_value, data.log_uniform_largest_value));
      k3 = Uniform(data.uniform_smallest_value, data.uniform_largest_value);
      d3 = Uniform(data.uniform_smallest_value, data.uniform_largest_value);

      l0 = Eucl_Dist(x0, y0, x_position, y_position);
    //  wout = Uniform(data.w_out_initial, data.w_out_final);
       wout = 0;

      s.push_back(Springs(k1, d1, k3, d3, l0, k, k-no_of_points_per_round, wout));
      }

      x0 = x_position;
      y0 = y_position;
      k++;

   }
   x_position = (j+1)*radius*cos((0));
   y_position = (j+1)*radius*sin((0));

   k1 = log10(Uniform(data.log_uniform_smallest_value, data.log_uniform_largest_value));
   d1 = log10(Uniform(data.log_uniform_smallest_value, data.log_uniform_largest_value));
   k3 = Uniform(data.uniform_smallest_value, data.uniform_largest_value);
   d3 = Uniform(data.uniform_smallest_value, data.uniform_largest_value);
   l0 = Eucl_Dist(x0, y0, x_position, y_position);
  // wout = Uniform(data.w_out_initial, data.w_out_final);
  //Temporarily remove
  wout = 0;

   if(no_of_points_per_round>2)
   {
   s.push_back(Springs(k1, d1, k3, d3, l0, (k-no_of_points_per_round), k-1, wout));
   }
 }


 n[no_of_points_per_round*(rounds-1)].set_Fixed_Node();
 n[2+no_of_points_per_round*(rounds-1)].set_Fixed_Node();
}

void Simulation::Delaunay_Triangulation_and_Spring_Creation()
{
    DelaunayTriangulation DT(abs(largest_x_position-smallest_x_position), abs(largest_y_position-smallest_y_position));
    double win = 0;
    double BeforeRand = 0;

    //This offset is to counteract the negativity.
  //  double offset = abs(w_in_initial)+abs(w_in_final);

  //  cout <<"offset is: " << offset;
  //  cout << endl << w_in_initial;
  //  cout <<endl << w_in_final;

    for(int i=0; i<N; i++)
    {
      win = Uniform(input_weight_smallest_value, input_weight_largest_value);
      BeforeRand = Uniform(0,1);
  //    cout <<"win is: "<< win << endl;
  //    win -= offset;
      cout <<"win is: "<< win << endl;
      if(BeforeRand<=input_connectivity_percentage) n[i].init_Input_Node(ux, uy, win);
      DT.AddPoint(Point(n[i].get_x_position(),n[i].get_y_position(),0));
    //  DT.AddPoint(Point(n[i].X_Position(), Point(n[i].Y_Position());
    }
    DT.print();

    Get_Triangles(DT);
    Initialize_Springs();
    execute();
    Output_For_Plot();

}

void Simulation::execute()
{
  double Fsum =0;
  double Fx_nodea =0;
  double Fy_nodea =0;
  double Fx_nodeb =0;
  double Fy_nodeb =0;
  double theta = 0;
  double l = 0;

  double nodea = 0;
  double nodeb = 0;

  double x0 = 0;
  double x1 = 0;
  double y0 = 0;
  double y1 = 0;

  // Todo: still needed for debugging DEBUG
  ofstream ofs("Node1.csv");
  ofstream ofs2("Node2.csv");
  ofstream ofs3("SampleForce.csv");
  ofstream bad("Badcoefficients.csv");

  double currentlength = 0;

  MatrixXd LearningMatrix(maxtimesteps, s.size());
//  MatrixXd TargetSignal(maxtimesteps, s.size()); // Todo: to fill TargetSignal matrix at init/constructor
 VectorXd TargetSignal(maxtimesteps);

  double outputsignal = 0;
  for(int j=0; j<s.size(); j++)
  {

    LearningMatrix(0,j)=s[j].return_Initial_Length();
    TargetSignal(0) = Target_Signal[0]; // Todo: Target matrix is of size #outputs x #timestep!!
  }


    // double uniform1 = 0;
    // double uniform2 = 0;

    //bool nodea1;
    // bool nodeb1;

  /////////////////////////////////
  //  SIMULATION LOOP
  /////////////////////////////////
  for(int i=1; i<maxtimesteps; i++)
  {
    //Each node has to be changed
    n[nodea].zero_updatecheck();
    n[nodeb].zero_updatecheck();

  for(int j=0;  j<s.size(); j++)
    {
       s[j].get_Force(Fsum);

       nodea = s[j].Nodea();
       nodeb = s[j].Nodeb();

       x0 = n[nodea].get_x_position();
       x1 = n[nodeb].get_x_position();
       y0 = n[nodea].get_y_position();
       y1 = n[nodeb].get_y_position();

       //Change position of first node
       // Todo: Check if "abs" is really needed
       theta = abs(Angle(x0, x1, y0, y1));

        // Todo: Check if this has to be so complicated
       if(x1>x0)
       {
           Fx_nodeb = X_Comp(Fsum, theta);
           Fx_nodea = -X_Comp(Fsum, theta);
       }
       if(y1>y0)
       {
           Fy_nodeb = Y_Comp(Fsum, theta);
           Fy_nodea = -Y_Comp(Fsum, theta);
       }

       if(x0>x1)
       {
       Fx_nodeb = -X_Comp(Fsum, theta);
       Fx_nodea = X_Comp(Fsum, theta);
       }

       if(y0>y1)
       {
       Fy_nodeb = -Y_Comp(Fsum, theta);
       Fy_nodea = Y_Comp(Fsum, theta);
       }

       // Update for input signal.
       // Todo: PROBLEM - adds up input force to same nodes over and over again
       Fx_nodea += n[nodea].return_Win()*Input_Signal[i];
       Fx_nodeb += n[nodeb].return_Win()*Input_Signal[i];

       //cout <<n[nodea].return_updatecheck()<< endl;
    //   cout <<n[nodeb].return_Win() << endl;

       if((!n[nodea].return_updatecheck())
       {
        n[nodea].Update(Fx_nodea, Fy_nodea, dt);
        n[nodea].change_updatecheck();
       }

       if((!n[nodeb].return_updatecheck())
       {
        n[nodeb].Update(Fx_nodeb, Fy_nodeb, dt);
        n[nodeb].change_updatecheck();
       }

       x0 = n[nodea].get_x_position();
       x1 = n[nodeb].get_x_position();

       y0 = n[nodea].get_y_position();
       y1 = n[nodeb].get_y_position();


       //Be very careful with the lengths here.
       l = Eucl_Dist(x0, y0, x1, y1);
       currentlength = l;



      // cout <<"Is this running?" << endl;
       LearningMatrix(i,j) = currentlength;
       TargetSignal(i) = Target_Signal[i];  // Todo: update is not needed for target signal


       s[j].update_Spring_State(dt, l);

       Fsum = 0;
       Fx_nodea =0;
       Fx_nodeb =0;
       Fy_nodea =0;
       Fy_nodeb =0;

      }


      outputsignal = 0;
    }

  //  TempMat = LearningMatrix;

  // TempMat.transposeInPlace();

    //cout <<TempMat;



      // TempMat =  ((TempMat*LearningMatrix).inverse())*TempMat;
    //   TempMat = TempMat*TargetSignal;

       //Populate the learning matrix and get weights
    //   LM = LearningMatrix;
       //Moore_Penrose_Pseudoinverse(LearningMatrix);
     //  LearningMatrix= LearningMatrix * TargetSignal;
    //   Populate_Learning_Weights(TempMat);

    LM = LearningMatrix;
    Moore_Penrose_Pseudoinverse(LearningMatrix);
    LearningMatrix= LearningMatrix * TargetSignal;
    Populate_Learning_Weights(LearningMatrix);
}

void Simulation::Moore_Penrose_Pseudoinverse(MatrixXd& L)
{
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

// Todo: The name does not tell us what this fucntion is doing
// Also, writing data to hard disk during simualtion slow it down a lot.
// Btw. any graphical output (even to the terminal) slows the process down a lot
// However, you could have every 1000 points and update message to show the use the simualtion is still going
// Btw. it is good to have a functionality to switch off any of these things by the user
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

  for(int i=0; i<maxtimesteps; i++)
  {
      for(int j=0; j<LM.cols(); j++)
      {
          outputsignal += Learning_Weights[j] * LM(i, j);
          //  if(i==0) output2 << Learning_Weights[j] << endl;
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
}  // end simulation loop



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
  return exp(Uniform(initial, finalvalue)/(2.302585093 ));
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

void Simulation::Remove_Duplicates(vector<vector<double>> &x)
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

unsigned int Simulation::Output_No_of_Edges()
{
  return EdgeList.size();
}

void Simulation::Output_Spring_And_Node_Positions()
{
  for (int i =0; i<s.size(); i++)
  {
    s[i].print_output();
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

      x0 = n[arraysubscript1].get_x_position();
      x1 = n[arraysubscript2].get_x_position();
      y0 = n[arraysubscript1].get_y_position();
      y1 = n[arraysubscript2].get_y_position();

      //These spring and damping coefficients are not giving different values
      k1 = log10(Uniform(log_uniform_smallest_value, log_uniform_largest_value));
      d1 = log10(Uniform(log_uniform_smallest_value,log_uniform_largest_value));
      k3 = Uniform(uniform_smallest_value, uniform_largest_value);
      d3 = Uniform(uniform_smallest_value, uniform_largest_value);

      ofs3 <<  k1 << endl;
      ofs4 <<  d1 << endl;
      ofs5 <<  k3 << endl;
      ofs6 <<  d3  << endl;

      l0 = Eucl_Dist(x0, y0, x1, y1);
      //Initial value for the output weights. I believe this was never used.
      //wout = Uniform(input_weight_smallest_value, input_weight_largest_value);
      wout = 0;

      s.push_back(Springs(k1, d1, k3, d3, l0, arraysubscript1, arraysubscript2, wout));
      s[i].print_output();
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

      Remove_Duplicates(EdgeList);
 //Remove Duplicates from EdgeNodeList
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

  int j=0;
  while(j<n.size())
  {
    nodesX << n[j].get_x_position();
    if(j<n.size()-1) nodesX<<",";
    nodesY << n[j].get_y_position();
    if(j<n.size()-1) nodesY<<",";
    j++;
  }

 }

}


Springs Simulation::Spring_Return(int i)
{
  return s[i];
}

Nodes Simulation::Node_Return(int i)
{
  return n[i];
}

unsigned int Simulation::Spring_List()
{
  return EdgeList.size();
}

//The helper functions for the Dynamical Systems class

DataSet::DataSet(double t0, double tmax, double dt)
{
  this->maxtimesteps = (int)((tmax - t0)/dt);
  this->tmax = tmax;
  this->dt = dt;
  this->t0 = t0;
}

void DataSet::SineWave(vector<double> &Sine_Wave)
{
  //double currenttime =0;
  for(int i =0; i<maxtimesteps; i++)
  {
     Sine_Wave.push_back(sin(t0+i*dt));
     cout <<sin(t0+i*dt) << endl;
  }
}
