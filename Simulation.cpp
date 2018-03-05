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
#include "Eigen/SVD"
#include "Eigen/LU"
#include <sstream>
#include <fstream>
#include <string>
using namespace std;
using namespace Eigen;


Simulation::Simulation(InitialDataValues &data, vector<double> &IS, vector<double> &TS, int wash_out_time, int learning_time, int learning_time_test)
{

  // read out all the data from the data structure
  this->N = data.N;

  this->input_connectivity_percentage = data.input_connectivity_percentage;
  this->num_input_nodes = (data.input_connectivity_percentage)*N;

  this->input_weight_smallest_value = data.min_input_weight;
  this->input_weight_largest_value = data.max_input_weight;


  this->t0 = data.t0;
  this->tmax = data.tmax;
  this->dt = data.dt;

  this->smallest_x_position = data.min_x_position;
  this->smallest_y_position = data.min_y_position;

  this->largest_x_position = data.max_x_position;
  this->largest_y_position = data.max_y_position;

  this-> min_k1 = data.min_k1;
  this-> min_k3 = data.min_k3;
  this-> min_d1 = data.min_d1;
  this-> min_d3 = data.min_d3;

  this-> max_k1 = data.max_k1;
  this-> max_k3 = data.max_k3;
  this-> max_d1 = data.max_d1;
  this-> max_d3 = data.max_d3;

  this->ux = data.ux;
  this->uy = data.uy;

  this->wash_out_time = wash_out_time;
  this->learning_time = learning_time;
  this->learning_time_test = learning_time_test;

  this->maxtimesteps = wash_out_time + learning_time + learning_time_test;

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
  this->num_input_nodes = (data.input_connectivity_percentage)*rounds*no_of_points_per_round;

  this->N = data.N;
  //Step 2: The positions of the nodes were initialized and 20% of the nodes are connected to the input.
  this->input_weight_smallest_value = data.min_input_weight;
  this->input_weight_largest_value = data.max_input_weight;

  this->smallest_x_position = data.min_x_position;     // smallest_x_position ?? name is not very descriptive (same for the others below)
  this->largest_x_position = data.max_x_position;
  this->smallest_y_position=data.min_y_position;
  this->largest_y_position = data.max_y_position;

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
      win = Uniform(data.min_input_weight, data.max_input_weight);
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
        //This has to be done from main, this is not good here.
      k1 = log10(Uniform(data.min_k1, data.max_k1));
      d1 = log10(Uniform(data.min_d1, data.max_d1));
      k3 = log10(Uniform(data.min_k3, data.max_k3));
      d3 = log10(Uniform(data.min_d1, data.max_d1));

      l0 = Eucl_Dist(x0, y0, x_position, y_position);
    //  wout = Uniform(data.w_out_initial, data.w_out_final);
      wout = 0;

      s.push_back(Springs(k1, d1, k3, d3, l0, k, k-1, wout));
      }

      //For radial pattern.
      if(j>0)
      {
        k1 = log10(Uniform(data.min_k1, data.max_k1));
        d1 = log10(Uniform(data.min_d1, data.max_d1));
        k3 = log10(Uniform(data.min_k3, data.max_k3));
        d3 = log10(Uniform(data.min_d3, data.max_d3));

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

   k1 = log10(Uniform(data.min_k1, data.max_k1));
   d1 = log10(Uniform(data.min_d1, data.max_d1));
   k3 = log10(Uniform(data.min_k3, data.max_k3));
   d3 = log10(Uniform(data.min_d3, data.max_d3));
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
    //Why is that abs? Double check this.
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
      if(BeforeRand<input_connectivity_percentage) win = Uniform(input_weight_smallest_value, input_weight_largest_value);
      else win = 0;
  //    win -= offset;
      n[i].init_Input_Node(ux, uy, win);
      DT.AddPoint(Point(n[i].get_x_position(),n[i].get_y_position(),0));
    //  DT.AddPoint(Point(n[i].X_Position(), Point(n[i].Y_Position());
    }
    DT.print();

    Get_Triangles(DT);
    Initialize_Springs();
    execute();
//    Output_For_Plot();

}

void Simulation::input_Magnitude_of_Chaos_Force(double k, const std::string& input, const std::string& input2)
{
  this->k = k;
  str = input;
  str2 = input2;
}

void Simulation::Reset_Simulation()
{
  for(int i=0; i<s.size(); i++)
  {
    s[i].set_original_length();
    s[i].set_x1(0);
    s[i].set_x2(0);
    n[s[i].Nodea()].original_positions();
    n[s[i].Nodeb()].original_positions();

    n[s[i].Nodea()].print_position();
    //n[s[i].Nodeb()].original_positions();
  }
}


void Simulation::execute()
{
  double Fsum =0;
  double Fx_nodea =0;
  double Fy_nodea =0;
  double Fx_nodeb =0;
  double Fy_nodeb =0;
  double l = 0;

  double nodea = 0;
  double nodeb = 0;

  double x0 = 0;
  double x1 = 0;
  double y0 = 0;
  double y1 = 0;

  double k1 = 0;
  double k3 = 0;
  double d1 = 0;
  double d3 = 0;

  //set k3 and d3 to 0.

  double x1spring = 0;
  double x2spring = 0;

  double x1new = 0;

  double vector_x = 0;
  double vector_y = 0;

  double alpha = 0;
  double beta = 0;

  double theta =0;

  //Input


  // Todo: still needed for debugging DEBUG
  ofstream Node1("Node1.csv");
  ofstream Node2("Node2.csv");
  ofstream ForceVec("SampleForce.csv");

  ofstream OutputPositionsx("NodePositionsx.csv");
  ofstream OutputVelocitiesx("NodeVelocitiesx.csv");
  ofstream OutputAccelerationsx("NodeAccelerationsx.csv");

  ofstream OutputPositionsy("NodePositionsy.csv");
  ofstream OutputVelocitiesy("NodeVelocitiesy.csv");
  ofstream OutputAccelerationsy("NodeAccelerationsy.csv");



//  MatrixXd LearningMatrix(maxtimesteps, s.size());


//Matrix for entire run
  MatrixXd LearningMatrix(maxtimesteps, s.size());
 //Matrix for learning phase
//  MatrixXd LearningMatrix1(learning_time, s.size());
  //MatrixXd LearningMatrix2(maxtimesteps, s.size());
//  MatrixXd LearningMatrix3(maxtimesteps, s.size());
//  MatrixXd TargetSignal(maxtimesteps, s.size()); // Todo: to fill TargetSignal matrix at init/constructor
//Target signal for entire run
  VectorXd TargetSignal(maxtimesteps);
 //Target Signal for learning_phase
//  VectorXd TargetSignal1(learning_time);

  double outputsignal = 0;

    // double uniform1 = 0;
    // double uniform2 = 0;

    //bool nodea1;
    // bool nodeb1;

  /////////////////////////////////
  //  SIMULATION LOOP
  /////////////////////////////////
//<<<<<<< HEAD
  //for(int i=0; i<maxtimesteps; i++)
//=======


//Inputs: Target Signal, nodes n, springs s

   for(int i=0; i<maxtimesteps; i++)
    {
  //      cout << "Time step " << i << endl;
        //
        TargetSignal(i) = Target_Signal[i];
      ////  if(i>=wash_out_time && i<(learning_time+wash_out_time)) TargetSignal1(i-wash_out_time) = Target_Signal[i];
    //  for(int l=0; l<s.size(); l++) n[l].Zero_Force();



        for(int j=0;  j<s.size(); j++)
        {
            nodea = s[j].Nodea();
            nodeb = s[j].Nodeb();
            //get_node_numbers(j);

            x0 = n[nodea].get_x_position();
            x1 = n[nodeb].get_x_position();

            //get_positions();


            y0 = n[nodea].get_y_position();
            y1 = n[nodeb].get_y_position();

            //get_positions();


            vector_x = x1 - x0;
            vector_y = y1 - y0;

            //get_x_component();
            //get_y_component();

            l = sqrt(vector_x*vector_x + vector_y*vector_y);

            //get_current_length();

            //Reintroduce theta temporarily
          //  theta = atan(vector_y/vector_x);

             //Reintroduce at later stage.
            alpha = vector_x/l;
            beta = vector_y/l;

            //get_direction_cosines();

            //theta = atan(vector_y/vector_x);

          //  cout <<"alpha minus cos theta" << alpha - cos(atan(theta)) << endl;
          //  cout <<"beta minus sin theta" << beta - sin(atan(theta)) << endl;;




            LearningMatrix(i,j) = l;  // Todo: update is not needed for target signal


            //s[j].update_Spring_State(dt, l);
            k1 = s[j].get_k1();
            d1 = s[j].get_d1();
            k3 = s[j].get_k3();
            d3 = s[j].get_d3();

          //  get_spring_coefficients();
          //  get_damping_coefficients();

            //Set k3 and d3 to 0 temporarily.
          //  k1 = 0;
        //    d1 = 0;

            x1new = l - s[j].return_Initial_Length();
            x1spring = s[j].return_x1();
            x2spring = ((x1new - x1spring)/dt);

            //get_x1_new(j);
            //get_x2_new(j)



          //  if(x1new == 0) cout <<"0 here." << endl;


            s[j].set_x2(x2spring);
            s[j].set_x1(x1new);

            //End k3 and d3

            Fsum =-k3*x1new*x1new*x1new - k1*x1new - d3*x2spring*x2spring*x2spring - d1*x2spring;
        //    if(i<450 && j==0) cout <<"Fsum is: "<< Fsum << endl;


          //  s[j].get_Force(Fsum);

        //    ForceVec << Fsum << endl;

            //if(i>=wash_out_time && i<(learning_time+wash_out_time)) LearningMatrix1(i-wash_out_time, j) = l;

            Fx_nodeb = Fsum*alpha;
            Fx_nodea = -Fx_nodeb;

            Fy_nodeb = Fsum*beta;
            Fy_nodea = -Fy_nodeb;

            n[nodea].Input_Force(Fx_nodea, Fy_nodea);
            n[nodeb].Input_Force(Fx_nodeb, Fy_nodeb);

          //  cout <<"l is: " << l;
          //  cout << endl;

            //   l = Eucl_Dist(x0, y0, x1, y1);
            // theta = abs(Angle(x0, x1, y0, y1));

            Fsum = 0;
            Fx_nodea =0;
            Fx_nodeb =0;
            Fy_nodea =0;
            Fy_nodeb =0;
            x1new = 0;
            x2spring =0;

            // s[j].update_Spring_State(dt, l)
        }

        //At the end of every timestep, the net force should be zero

          for(int l=0; l<n.size(); l++)
          {
              //Output NodePositions and Velocity and Acceleration
              if(l==1) OutputPositionsx << n[l].get_x_position() <<",";
              if(l==1) OutputVelocitiesx << n[l].get_x_velocity() <<",";
              if(l==1) OutputAccelerationsx << n[l].get_x_acceleration() <<",";
              //Output NodePositions and Velocity and Acceleration
              if(l==1) OutputPositionsy << n[l].get_y_position();
              if(l==1) OutputVelocitiesy << n[l].get_y_velocity();
              if(l==1) OutputAccelerationsy << n[l].get_y_acceleration();

              //This puts in the forces due to the input force
              if(i>0) n[l].Input_Force(0.001, 0);
              //if(i>0) n[l].Input_Force(0, 0);
              //This puts in the forces due to the input force
              //n[k].Input_Force(0, return_Win();
              //Update force on node due to dt
              n[l].Update(dt);
              //Zero Force
              n[l].Zero_Force();
         }

         OutputPositionsx << endl;
         OutputVelocitiesx << endl;
         OutputAccelerationsx << endl;

         OutputPositionsy << endl;
         OutputVelocitiesy << endl;
         OutputAccelerationsy << endl;

       }


      //Jacobian singular value decomposition for Moore Penrose pseudoinverse

  //    LM = LearningMatrix;
  //    JacobiSVD<MatrixXd> svd(LearningMatrix, ComputeThinU | ComputeThinV);
  //    MatrixXd Cp = svd.matrixV() * (svd.singularValues().asDiagonal()).inverse() * svd.matrixU().transpose();
  //    MatrixXd original = svd.matrixU() * (svd.singularValues().asDiagonal()) * svd.matrixV().transpose();
  //    Cp = Cp * TargetSignal;
  //    Populate_Learning_Weights(Cp);



      LM = LearningMatrix;
    //  MatrixXd LeftInverse = ((LearningMatrix.transpose()*LearningMatrix).inverse())*LearningMatrix.transpose();
      //LeftInverse = LeftInverse * TargetSignal;
      //Populate_Learning_Weights(LeftInverse);



     //left inverse

  //    LM = LearningMatrix;
  //    MatrixXd LeftInverse = ((LearningMatrix.transpose()*LearningMatrix).inverse())*LearningMatrix.transpose();
  //    LeftInverse = LeftInverse * TargetSignal;
  //    Populate_Learning_Weights(LeftInverse);









        //SVF decomp using fast method
    //  LM = LearningMatrix;
    //  BDCSVD<MatrixXd> svd(LM, ComputeThinU | ComputeThinV);
  //    MatrixXd Cp = svd.matrixV() * (svd.singularValues().asDiagonal()).inverse() * svd.matrixU().transpose();
    //  MatrixXd original = svd.matrixU() * (svd.singularValues().asDiagonal()) * svd.matrixV().transpose();
    //  Cp = Cp * TargetSignal;
    //  Populate_Learning_Weights(Cp);









    //     LM = LearningMatrix;
//    Moore_Penrose_Pseudoinverse(LearningMatrix);
//    LearningMatrix= LearningMatrix * TargetSignal;
//    Populate_Learning_Weights(LearningMatrix);
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
  //  cout << L(j,0) << endl;
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


  ofstream outputsignal("outputsignal.csv");
  ofstream learningweights("learningweights.csv");
  ofstream targetsignal("targetsignal.csv");
  //ofstream learningmatrix("learningmatrix.csv");
  ofstream learningmatrix("learningmatrix.csv");
  ofstream inputsignalcheck("inputsignalcheck.csv");

  ofstream chaoscheck(str);

//  ofstream mse(str);

  //double outputsignal = 0;

  double wjej = 0;


  double currenttime = 0;


  //For normalised LearningMatrix
  double average = 0;
  double std = 0;
  //double outputLM = 0;

  //(learningmatrix/mean(learningmatrix) - 1 )*(1/(std(learningmatrix/mean(learningmatrix))))

//Temp stop

//  for(int i=0; i<learning_time_test; i++)
   for(int i = 0; i<maxtimesteps; i++)
  {
      for(int j=0; j<s.size(); j++)
      {
      //    outputsignal += Learning_Weights[j] * LM(i+wash_out_time+learning_time, j);
        //  outputsignal += Learning_Weights[j] * LM(i, j);
            //outputLM =  LM(i, j);
          //  wjej += LM(i, j);

            learningmatrix<<LM(i, j) <<",";

        //  if(i==0) output2 << Learning_Weights[j] << endl;
      }

    //  chaoscheck << wjej << endl;
    //  wjej = 0;

      learningmatrix <<endl;

    //  Output_Signal.push_back(outputsignal);
      currenttime = t0 + i*dt;

    //  output << currenttime <<"," << Output_Signal.at(i);
    //  output << endl;
      //output3 <<currenttime <<"," << Target_Signal.at(i);
    //  output3 << endl;
      inputsignalcheck << Input_Signal.at(i) << endl;
      targetsignal << Target_Signal.at(i) << endl;
    //  outputsignal = 0;
  }

//  cout <<"The mean squared error of the output signal versus the target signal is: " << MSE(Output_Signal, Target_Signal);
//  cout <<endl;
//  mse << MSE(Output_Signal, Target_Signal);
//  mse << endl;
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

double Simulation::Rand_In_Range_Exp_k1()
{
  cout << "min is:" << min_k1 << endl;
  cout << "max is:" << max_k1 << endl;

  double log10min = log10(min_k1);
  double log10max = log10(max_k1);

  cout << log10min << endl;
  cout << log10max << endl;

  double return_value = ((log10max-log10min)*Uniform(0,1))+log10min;
  return pow(10, return_value);
}

double Simulation::Rand_In_Range_Exp_d1()
{
  cout << "min is:" << min_d1 << endl;
  cout << "max is:" << max_d1 << endl;

  double log10min = log10(min_d1);
  double log10max = log10(max_d1);

  cout << log10min << endl;
  cout << log10max << endl;

  double return_value = ((log10max-log10min)*Uniform(0,1))+log10min;
  return pow(10, return_value);;
}

double Simulation::Rand_In_Range_Exp_k3()
{
  cout << "min is:" << min_k3 << endl;
  cout << "max is:" << max_k3 << endl;

  double log10min = log10(min_k3);
  double log10max = log10(max_k3);

  cout << log10min << endl;
  cout << log10max << endl;
  double return_value = ((log10max-log10min)*Uniform(0,1))+log10min;
  return pow(10, return_value);
}

double Simulation::Rand_In_Range_Exp_d3()
{

  double log10min = log10(min_d3);
  double log10max = log10(max_d3);
  double return_value = ((log10max-log10min)*Uniform(0,1))+log10min;
  return pow(10, return_value);
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
  //cout <<"The number of edges for: " <<N << " mass points is: " << EdgeList.size() << endl;
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

    ofstream k1output("k1output.csv");  // Todo: ofs3, etc. is really bad!
    ofstream d1output("d1output.csv");
    ofstream k3output("k3output.csv");
    ofstream d3output("d3output.csv");
    ofstream originallengthoutput("originallengthoutput.csv");


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
      k1 = Rand_In_Range_Exp_k1();
      d1 = Rand_In_Range_Exp_d1();
      k3 = Rand_In_Range_Exp_k3();
      d3 = Rand_In_Range_Exp_d3();

      k1output << k1 <<endl;
      k3output << d1 <<endl;
      d1output << k3 <<endl;
      d3output << d3 <<endl;

      l0 = Eucl_Dist(x0, y0, x1, y1);
      originallengthoutput << l0 << endl;
      //Initial value for the output weights. I believe this was never used.
      //wout = Uniform(input_weight_smallest_value, input_weight_largest_value);
      wout = 0;

      s.push_back(Springs(k1, d1, k3, d3, l0, arraysubscript1, arraysubscript2, wout));
      s[i].print_output();

      cout <<"Edgelist is: " << EdgeList.size() << endl;
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
  for(int i =0; i<maxtimesteps; i++)
  {
     Sine_Wave.push_back(sin(t0+i*dt));
     cout <<sin(t0+i*dt) << endl;
  }
}
