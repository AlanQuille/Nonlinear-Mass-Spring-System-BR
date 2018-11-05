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
#include <iomanip>
using namespace std;
using namespace Eigen;


Simulation::Simulation(InitialDataValues &data, vector<double> &IS, vector<double> &TS, int wash_out_time, int learning_time, int learning_time_test)
{

  // read out all the data from the data structure
  this->N = data.N;

  this->input_connectivity = data.input_connectivity;
  this->num_input_nodes = (data.input_connectivity)*N;

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

  //Total time
  this->maxtimesteps = wash_out_time + learning_time + learning_time_test;

  Target_Signal = TS;
  Input_Signal = IS;

  Initialize_Nodes(smallest_x_position, largest_x_position, smallest_y_position, largest_y_position);
  Delaunay_Triangulation_and_Spring_Creation();

  Initialize_Springs(data);
  //update(true, false);
 // Mean_Squared_Error = output_LearningMatrix_and_MeanSquaredError();
//  Output_For_Plot();
}




Simulation::Simulation(double radius, int rounds, int no_of_points_per_round, InitialDataValues &data, vector<double> &IS, vector<double> &TS, int wash_out_time, int learning_time, int learning_time_test, bool springs_identical, bool random_node_positions, double mean, double stdev)
{
  //this->input_connectivity_percentage = data.input_connectivity_percentage;
  this->num_input_nodes = (data.input_connectivity)*rounds*no_of_points_per_round;

  this->N = data.N;
  //Step 2: The positions of the nodes were initialized and 20% of the nodes are connected to the input.
  this->input_weight_smallest_value = data.min_input_weight;
  this->input_weight_largest_value = data.max_input_weight;

  this-> min_k1 = data.min_k1;
  this-> min_k3 = data.min_k3;
  this-> min_d1 = data.min_d1;
  this-> min_d3 = data.min_d3;
  
  this-> max_k1 = data.max_k1;
  this-> max_k3 = data.max_k3;
  this-> max_d1 = data.max_d1;
  this-> max_d3 = data.max_d3;
  
  if(springs_identical) identical = 1;
  else identical = 0;

  //Learning phase
  Initialize_Nodes(radius, rounds, no_of_points_per_round, data);

  this->t0 = data.t0;
  this->tmax = data.tmax;
  this->dt = data.dt;

  this->wash_out_time = wash_out_time;
  this->learning_time = learning_time;
  this->learning_time_test = learning_time_test;

  this->maxtimesteps = wash_out_time + learning_time + learning_time_test;
  
  Target_Signal = TS;
  Input_Signal = IS;
  
  random_positions = random_node_positions;
  mu = mean;
  sigma = stdev;

  //update(true);
  //Mean_Squared_Error = output_LearningMatrix_and_MeanSquaredError();
 // Output_For_Plot();
}

//Need this function to change input_connectivity input input_connectivity



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

     //Test to see whether the reason why you're getting those 0 springs is because of fixed nodes.

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

     double random_factor_x = 0;
     double random_factor_y = 0;

     bool odd_even_check = 1;


     int k =0;
     
     
     //If the springs/threads are all identical, only do the spring and damping coefficients once.
    k1 = Rand_In_Range_Exp(data.min_k1, data.max_k1);
    d1 = Rand_In_Range_Exp(data.min_d1, data.max_d1);

    k3 = Rand_In_Range_Exp(data.min_k3, data.max_k3);
    d3 = Rand_In_Range_Exp(data.min_d3, data.max_d3);
    
    k1_identical = k1;
    d1_identical = d1;
    
    k3_identical = k3;
    d3_identical = d3;
    
     

    for(int j=0; j<rounds; j++)
    {
      x0 = (j+1)*radius*cos((0));
      y0 = (j+1)*radius*sin((0));

      for(int i=0; i<no_of_points_per_round; i++)
      {
         //So this
         x_position = (j+1)*radius*cos((i*angle));
         y_position = (j+1)*radius*sin((i*angle));
         
         //Add in gaussian noise
         x_position += generate_Gaussian_Noise(mu, sigma);
         y_position += generate_Gaussian_Noise(mu, sigma);


         Nodes node(x_position, y_position);
         
         
         //change mass of each node
         node.change_Mass(data.mass_of_nodes);

         n.push_back(node);
         


         if(i>0)
         {

         if(!identical) k1 = Rand_In_Range_Exp(data.min_k1, data.max_k1);
         if(!identical) d1 = Rand_In_Range_Exp(data.min_d1, data.max_d1);

         cout << k1 << endl;
         cout << d1 << endl;

         if(!identical) k3 = Rand_In_Range_Exp(data.min_k3, data.max_k3);
         if(!identical) d3 = Rand_In_Range_Exp(data.min_d3, data.max_d3);


      cout << k3 << endl;
      cout << d3 << endl;


         l0 = Eucl_Dist(x0, y0, x_position, y_position);
         wout = 0;

         s.push_back(Springs(k1, d1, k3, d3, l0, k, k-1, wout));
         }

         //For radial pattern.
         if(j>0)
         {

         if(!identical) k1 = Rand_In_Range_Exp(data.min_k1, data.max_k1);
         if(!identical) d1 = Rand_In_Range_Exp(data.min_d1, data.max_d1);
         if(!identical) k3 = Rand_In_Range_Exp(data.min_k3, data.max_d3);
         if(!identical) d3 = Rand_In_Range_Exp(data.min_d3, data.max_d3);
         cout << "d3 is: " << d3 << endl;

      //   k3 = 0;
      //   d3 = 0;

      //   k1 = data.min_k1;
      //   d1 = data.min_d1;


         //k3 = 0;
         //d3 = 0;
        // k3 = data.max_k3;
    //     d3 = data.max_d3;

         //l0 = Eucl_Dist(x0, y0, x_position, y_position);
         l0 = radius;
       //  wout = Uniform(data.w_out_initial, data.w_out_final);
          wout = 0;

/*
         if(odd_even_check)
         {
         if(j%2==0 && i%2==0) s.push_back(Springs(k1, d1, k3, d3, l0, k, k-no_of_points_per_round, wout));
         if(j%2==1 && i%2==1) s.push_back(Springs(k1, d1, k3, d3, l0, k, k-no_of_points_per_round, wout));
         odd_even_check = false;
         }

         else
         {
         if(j%2==0 && i%2==1) s.push_back(Springs(k1, d1, k3, d3, l0, k, k-no_of_points_per_round, wout));
         if(j%2==1 && i%2==0) s.push_back(Springs(k1, d1, k3, d3, l0, k, k-no_of_points_per_round, wout));
      //   odd_even_check = true;
         }


*/

         s.push_back(Springs(k1, d1, k3, d3, l0, k, k-no_of_points_per_round, wout));



         }

         x0 = x_position;
         y0 = y_position;
         k++;

      }
      x_position = (j+1)*radius*cos((0));
      y_position = (j+1)*radius*sin((0));

      if(!identical) k1 = Rand_In_Range_Exp(data.min_k1, data.max_k1);
      cout <<"k1 is: " <<k1 << endl;
      if(!identical) d1 = Rand_In_Range_Exp(data.min_d1, data.max_d1);
      cout <<"d1 is: " <<d1 << endl;
      if(!identical) k3 = Rand_In_Range_Exp(data.min_k3, data.max_k3);
      cout <<"k3 is: " <<k3 << endl;
      if(identical) d3 = Rand_In_Range_Exp(data.min_d3, data.max_d3);
      cout <<"d3 is: " <<d3 << endl;

    //  k3 = 0;
    //  d3 = 0;


      l0 = Eucl_Dist(x0, y0, x_position, y_position);
      cout <<"l0 is: " <<l0 << endl;
     // wout = Uniform(data.w_out_initial, data.w_out_final);
     //Temporarily remove
     wout = 0;

      if(no_of_points_per_round>2)
      {
      s.push_back(Springs(k1, d1, k3, d3, l0, (k-no_of_points_per_round), k-1, wout));
      }
    }

    vector<int> nos;

    for(int i=0; i<n.size()-1; i++)
    {
    nos.push_back(i);
    random_shuffle(nos.begin(), nos.end());

    }


    int N = rounds * no_of_points_per_round + 1;

    int input_node_nums = 0;

    input_node_nums = data.input_connectivity*(N);

    cout << input_node_nums << endl;

    int randomnum = 0;
    double win;

    cout << "No of input nodes is: " << input_node_nums << endl;
    cout << "No of input nodes is: " << input_node_nums << endl;



    n[no_of_points_per_round*(rounds-1)].set_Fixed_Node();
    n[(no_of_points_per_round/2)+no_of_points_per_round*(rounds-1)].set_Fixed_Node();


    cout << "Outer fixed node 1 " << no_of_points_per_round*(rounds-1) << endl;
    cout << "Outer fixed node 2 " <<(no_of_points_per_round/2) + no_of_points_per_round*(rounds-1) << endl;




//Temporarily, input nodes will be fixed for this spider web.


//This fixed the outer ring of the nodes. Will get user access for this.
  //  for(int i=0; i<N; i++)
  //   {
      //Input weights for the number of input_connectivitiy nodes.
  //    win = Uniform(data.min_input_weight, data.max_input_weight);
    //  randomnum = (int)Uniform(0, N);

//      randomnum = (int)Uniform(0, no_of_points_per_round*(rounds-1));

      //Just

//      if(i<input_node_nums)
//      {
//      n[randomnum].init_Input_Node(data.ux, data.uy, win);
//      cout <<"The input node is:" << randomnum << endl;
//      }



      //Make sure fixed nodes are not input nodes

  //    if(i<input_node_nums)
  //    {
    //  If it is not a fixed node.
      
   //   while(n[randomnum].is_Fixed_Node())
   //     {
        //C++ typecasting rounds down (truncates) but this is fine going from 0 to N-1.
    //    randomnum = (int)Uniform(0, N);
    //    }

    //  n[randomnum].init_Input_Node(data.ux, data.uy, win);
      

   //   }
// for non fixed, get rid of this.

  //    if(i>=no_of_points_per_round*(rounds-1) && i<((no_of_points_per_round)+no_of_points_per_round*(rounds-1)))
  //    {
  //      n[i].set_Fixed_Node();
  //      cout <<"The: " <<i <<"th" << " fixed node is fixed." << endl;
  //    }
    

    //        if(i<no_of_points_per_round*(rounds-1) )
      //      {

    //        randomnum = (int)Uniform(0, no_of_points_per_round*(rounds-1));
    //        n[randomnum].init_Input_Node(data.ux, data.uy, win);
         //If it is not a fixed node.
            /*
            while(n[randomnum].is_Fixed_Node())
              {
              //C++ typecasting rounds down (truncates) but this is fine going from 0 to N-1.
              randomnum = (int)Uniform(0, N);
              }

            n[randomnum].init_Input_Node(data.ux, data.uy, win);
            */

      //       }


   //  }

     // add in central node here.

     cout << randomnum << endl;




     Nodes central_node(0, 0);
     n.push_back(central_node);

     win = Uniform(data.min_input_weight, data.max_input_weight);
    // n[5].init_Input_Node(data.ux, data.uy, win);
    //for another round
    // randomnum = (int)Uniform(0, n.size());
       randomnum = 10;
       
      while(randomnum!=no_of_points_per_round*(rounds-1) && randomnum!=((no_of_points_per_round/2)+no_of_points_per_round*(rounds-1)))
     {
	 randomnum = (int)Uniform(0, n.size());
     }
     randomnum = 6;
     
     //n[10].init_Input_Node(data.ux, data.uy, win);
     n[randomnum].init_Input_Node(data.ux, data.uy, win);
     //n[10].init_Input_Node(data.ux, data.uy, win);
     
    
  
     
     
    // n[randomnum].init_Input_Node(data.ux, data.uy, win);
     
     no_of_input_node = randomnum;

     cout << "number of nodes is: " << n.size() << endl;

    // win = Uniform(data.min_input_weight, data.max_input_weight);
  //   n[16].init_Input_Node(data.ux, data.uy, win);




     for(int i=0; i<no_of_points_per_round; i++)
     {
       k1 = Rand_In_Range_Exp(data.min_k1, data.max_k1);
       cout <<"k1 is: " <<k1 << endl;
       d1 = Rand_In_Range_Exp(data.min_d1, data.max_d1);
       cout <<"d1 is: " <<d1 << endl;
       k3 = Rand_In_Range_Exp(data.min_k3, data.max_k3);
       cout << "k3 is: " <<k3 << endl;
       d3 = Rand_In_Range_Exp(data.min_d3, data.max_d3);
       cout << "d3 is: " <<d3 << endl;

       l0 = Eucl_Dist(0, 0, n[i].get_x_Position(), n[i].get_y_Position());

       wout = 0;

       s.push_back(Springs(k1, d1, k3, d3, l0, n.size()-1, i, wout));

     }



   }
   
int Simulation::return_input_Node()
{
   	return no_of_input_node;
}

void Simulation::Delaunay_Triangulation_and_Spring_Creation()
{
    //Why is that abs? Double check this.
    //I think the delaunay triangulation needs to take in the absolute difference between largest and smallest.
    DelaunayTriangulation DT(abs(largest_x_position-smallest_x_position), abs(largest_y_position-smallest_y_position));

    double win = 0;
    int input_node_nums = (int)(input_connectivity*((int)N));

    cout << input_node_nums << endl;

    int randomnum;

    int num_of_input_nodes = 0;
    int num_of_fixed_nodes = 0;


    for(int i=0; i<N; i++)
    {
      //Input weights for the number of input_connectivitiy nodes.
      win = Uniform(input_weight_smallest_value, input_weight_largest_value);
      randomnum = (int)Uniform(0, N);

      //Make sure fixed nodes are not input nodes
      if(i<input_node_nums)
      {
      //If it is not a fixed node.
      while(n[randomnum].is_Fixed_Node())
      {
        //C++ typecasting rounds down (truncates) but this is fine going from 0 to N-1.
        randomnum = (int)Uniform(0, N);
      }

      n[randomnum].init_Input_Node(ux, uy, win);
      num_of_input_nodes++;

      }





      DT.AddPoint(Point(n[i].get_x_Position(),n[i].get_y_Position(),0));

    }

    cout <<"The total number of input nodes is: " << num_of_input_nodes << endl;

    DT.print();
    Get_Triangles(DT);

}

void Simulation::input_Magnitude_of_Chaos_Force(double k, const std::string& input, const std::string& input2)
{
  this->k = k;
  str = input;
  str2 = input2;
}



void Simulation::Reset_Simulation()
{
	
	
  //Springs
  for(int i=0; i<s.size(); i++)
  {
    s[i].set_Original_Length();

    s[i].set_x1(0);
    s[i].set_x2(0);

   // s[i].set_Force_0();
  //  n[s[i].Nodea()].print_position();

  }
  
  for(int j=0; j<n.size(); j++)
  {
    n[j].original_Positions();
    n[j].zero_Accel_and_Vel();
  }


}


void Simulation::update(bool bias_learning, bool impulse_response_or_input_signal)
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


  double x1spring = 0;
  double x2spring = 0;

  double x1new = 0;

  double vector_x = 0;
  double vector_y = 0;

  double alpha = 0;
  double beta = 0;

  double theta =0;
  
  this->bias_learning = bias_learning;

  //Input
  
  //Global stability check. if l is nAN for instance it is automatically unstable. Default stability = true
 // bool stability = true;
 //Now moved to simulation.h


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
  
  ofstream InputSignal("InputSignal.csv");

  //Learning matrix for entire run, learning phase and testing phase.
  MatrixXd LearningMatrix(maxtimesteps, s.size());
  MatrixXd LearningMatrix2(learning_time, s.size());
  MatrixXd LearningMatrix3(learning_time_test, s.size());

  //Detto but for Target Signal
  VectorXd TargetSignal(maxtimesteps);
  VectorXd TargetSignal2(learning_time);
  VectorXd TargetSignal3(learning_time_test); 
  

  
  vector<double> Bounded_Max_vector;
  double Bounded_Max = 0;
  
  vector<double> Sum_vector;
  double Sum = 0;
  
 // vector<double> Max_Bounded;
  
 // for(int k=0; k<s.size(); k++) Max_Bounded.at(k) = 0;
  

  double outputsignal = 0;
  

  /////////////////////////////////
  //  SIMULATION LOOP
  /////////////////////////////////
//<<<<<<< HEAD
  //for(int i=0; i<maxtimesteps; i++)
//=======

//Inputs: Input SignalTarget Signal, nodes n, springs s
//Outputs: LearningMatrix

//For some reason the number of threads/webs is gone in simulation. I will save them again.

number_of_threads_or_webs = s.size();
cout << "The number of springs is: " << s.size() << endl;

bool wcheck = 0;

double l0 =0;

   for(int i=0; i<maxtimesteps; i++)
    {
  //      cout << "Time step " << i << endl;
        //
        TargetSignal(i) = Target_Signal[i];
        if(i>=wash_out_time && i<(wash_out_time+learning_time)) TargetSignal2(i-wash_out_time) = Target_Signal[i];
        if(i>=(wash_out_time+learning_time)) TargetSignal3(i-wash_out_time-learning_time) = Target_Signal[i];

        for(int j=0;  j<s.size(); j++)
        {
            nodea = s[j].Nodea();
            nodeb = s[j].Nodeb();

            x0 = n[nodea].get_x_Position();
            x1 = n[nodeb].get_x_Position();

            y0 = n[nodea].get_y_Position();
            y1 = n[nodeb].get_y_Position();

            vector_x = x1 - x0;
            vector_y = y1 - y0;

            l = sqrt(vector_x*vector_x + vector_y*vector_y);

            alpha = vector_x/l;
            beta = vector_y/l;


            LearningMatrix(i,j) = l;  // Todo: update is not needed for target signal

          //  cout << l << endl;
            if(i>=wash_out_time && i<(wash_out_time+learning_time)) LearningMatrix2(i-wash_out_time, j) = l;
            if(i>=(wash_out_time+learning_time)) LearningMatrix3(i-wash_out_time-learning_time, j) = l;

            k1 = s[j].get_k1();
            d1 = s[j].get_d1();
            k3 = s[j].get_k3();
            d3 = s[j].get_d3();
            l0 = s[j].return_Initial_Length();

            if(i==(maxtimesteps-1))
            {
            cout << "k1 is: " <<s[j].get_k1() << endl;
            cout << "d1 is: " <<s[j].get_d1() << endl;

            cout << "k3 is: " <<s[j].get_k3() << endl;
            cout << "d3 is: " <<s[j].get_d3()<< endl;
            }


            x1new = l - s[j].return_Initial_Length();




            x1spring = s[j].return_x1();


        //    cout <<"x1spring is: " << x1spring << endl;
      //      cout <<"x1new is: " << x1new << endl;

            //cout << "For spring" <<" " <<j <<"x1new is : "<< x1new<< endl;
            x2spring = ((x1new - x1spring)/dt);
            //cout << "For spring" <<" " <<j <<"x2spring is : "<< x2spring<< endl;

            s[j].set_x2(x2spring);
            s[j].set_x1(x1new);

            Fsum =-k3*x1new*x1new*x1new - k1*x1new - d3*x2spring*x2spring*x2spring - d1*x2spring;
            
          //  cout << Fsum << endl;


            Fx_nodeb = Fsum*alpha;
            Fx_nodea = -Fx_nodeb;

            Fy_nodeb = Fsum*beta;
            Fy_nodea = -Fy_nodeb;

            n[nodea].input_Force(Fx_nodea, Fy_nodea);
            n[nodeb].input_Force(Fx_nodeb, Fy_nodeb);
            
            
           // if(nodea == 9 || nodeb ==9) cout << l << endl;  



            Fsum = 0;
            Fx_nodea =0;
            Fx_nodeb =0;
            Fy_nodea =0;
            Fy_nodeb =0;
            x1new = 0;
            x2spring =0;
            
            //isnan output
            if(isnan(l)) stability = false;
             
            
            
            //Sum_vector.at(j) = +l;
            
            
          //  if(i<= half_time) Max_Bounded.at(i) = abs(l);
        }




          for(int l=0; l<n.size(); l++)
          {

              //Input force to input nodes from input signal. If impulse_response_or_input_signal is true or 1, than input impulse response. If false, input inputsisgnal
            if(impulse_response_or_input_signal)      
			{
		        if(n[l].is_Input_Node()==true && i==0) n[l].input_Force(1,0);	
			} 
			else 
			{
				if(n[l].is_Input_Node()==true) n[l].input_Force(n[l].return_Win()*Input_Signal[i],0);	
			}
			 
			 
			 
			
			//if(n[l].is_Input_Node()==true && i==0) n[l].input_Force(1,0);
         //   if(n[l].is_Input_Node()==true) n[l].Input_Force(n[l].return_Win()*Treble_Sine_Function(2.11, 3.73, 4.33, 0.001, i, scaling_factor),0);
           // Input impulse esponse.
            //if(n[l].is_Input_Node()==true && i==0) n[l].input_Force(1,0);
            
              //Change the node position, velocity and acceleration in response.
              n[l].update(dt);
              //At the end of the loop, each node has no force acting on it.
              n[l].zero_Force();
         }

       }
       
      

      LM = LearningMatrix;
      LM2 = LearningMatrix2;
      LM3 = LearningMatrix3;
      
      TS = TargetSignal;
      TS2 = TargetSignal2;
      TS3 = TargetSignal3;
      
      //Check for stability and than perform Moore_Penrose pseudoinverse
      //stability_Check();
      if(!impulse_response_or_input_signal) Moore_Penrose_Pseudoinverse_and_Learning_Weights();
}


void Simulation::stability_Check()
{
	   int half_time = maxtimesteps/2;
	   
	   MatrixXd Lee1(maxtimesteps, s.size());
  
       MatrixXd Lee(half_time, s.size());
       
       MatrixXd Le(half_time, s.size());
       
       
       
       MatrixXd S_A = LM.colwise().mean().replicate(maxtimesteps, 1);
       
       Lee1 = LM - S_A;
       
       
       Lee = Lee1.block(0, 0, half_time, s.size());
       
       Le = Lee1.block(half_time, 0, half_time, s.size());
       
       Lee =Lee.cwiseAbs();
       Le = Le.cwiseAbs();
       
       double minCoeff = (Lee.colwise().maxCoeff() - Le.colwise().maxCoeff()).minCoeff();
       
       cout << "Minimum coefficient" <<" " <<  minCoeff << endl;
       
       
       if(minCoeff<0 || isnan(minCoeff)) stability = false;
       
       if(stability) cout <<"The structure is stable" << endl;
       else cout <<"The structure is not stable." << endl;
	
}




void Simulation::Moore_Penrose_Pseudoinverse_and_Learning_Weights()
{
	  MatrixXd LearningMatrix(maxtimesteps, s.size());
      MatrixXd LearningMatrix2(learning_time, s.size());
      MatrixXd LearningMatrix3(learning_time_test, s.size());
      
      VectorXd TargetSignal(maxtimesteps);
      VectorXd TargetSignal2(learning_time);
      VectorXd TargetSignal3(learning_time_test); 
      
      LearningMatrix = LM;
      LearningMatrix2 = LM2;
      LearningMatrix3 = LM3;
      
      TargetSignal = TS;
      TargetSignal2 = TS2;
      TargetSignal3 = TS3;

	
	 //Add in column for bias learning, linear regression weight.
	 if(bias_learning)
      {
      LearningMatrix2.conservativeResize(LearningMatrix2.rows(), LearningMatrix2.cols()+1);
      LearningMatrix2.col(LearningMatrix2.cols() - 1) = VectorXd::Ones(learning_time);
      }

      //Jacobian singular value decomposition for pseudoinverse.
      JacobiSVD<MatrixXd> svd(LearningMatrix2, ComputeThinU | ComputeThinV);
      MatrixXd Cp = svd.matrixV() * (svd.singularValues().asDiagonal()).inverse() * svd.matrixU().transpose();;
      VectorXd LearningWeightsVector = Cp *TargetSignal2;

      cout << LearningWeightsVector << endl;
      
      //Add in column for bias learning, linear regression weight.
      if(bias_learning)
      {
      LearningMatrix3.conservativeResize(LearningMatrix3.rows(), LearningMatrix3.cols()+1);
      LearningMatrix3.col(LearningMatrix2.cols()-1) = VectorXd::Ones(learning_time_test);
      }

      //Output vector.
      Output = LearningMatrix3*LearningWeightsVector;	
}

void Simulation::Moore_Penrose_Pseudoinverse(MatrixXd& L)
{
  L = L.completeOrthogonalDecomposition().pseudoInverse();
}

double Simulation::return_MSE()
{
  return Mean_Squared_Error;
}

vector<double>& Simulation::Return_Learning_Weights()
{
  return Learning_Weights;
}

MatrixXd& Simulation::Return_Learning_Matrix()
{
  return LM;
}

void Simulation::Populate_Learning_Weights(VectorXd& L)
{
  for(int j=0; j<s.size(); j++)
  {
  //  cout << L(j,0) << endl;
    Learning_Weights.push_back(L(j));
  }
}

// Also, writing data to hard disk during simualtion slow it down a lot.
// Btw. any graphical output (even to the terminal) slows the process down a lot
// However, you could have every 1000 points and update message to show the use the simualtion is still going
// Btw. it is good to have a functionality to switch off any of these things by the user
double Simulation::output_LearningMatrix_and_MeanSquaredError()
{
  ofstream output("outputsignal.csv"); output.precision(15);
  ofstream learningweights("learningweights.csv"); learningweights.precision(15);
  ofstream targetsignal("targetsignal.csv");  targetsignal.precision(15);
  ofstream targetsignal2("targetsignal2.csv");  targetsignal.precision(15);
  ofstream targetsignal3("targetsignal3.csv");  targetsignal.precision(15);

  ofstream learningmatrix("learningmatrix.csv");  learningmatrix.precision(15);
  ofstream learningmatrix2("learningmatrix2.csv");  learningmatrix.precision(15);
  ofstream learningmatrix3("learningmatrix3.csv");  learningmatrix.precision(15);
  ofstream outputsignal("outputsignal.csv");  learningmatrix.precision(15);

  ofstream inputsignalcheck("inputsignalcheck.csv");  inputsignalcheck.precision(15);

  ofstream chaoscheck(str);

  //double outputsignal = 0;
  double wjej = 0;
  double currenttime = 0;
  double currentvalue = 0;
  double average = 0;
  double std = 0;
  double Mean_squared_error = 0;

  vector<double> Output_Signal;
  vector<double> Test_Data;

  VectorXd Target_Here(learning_time_test);

//  for(int i=0; i<learning_time_test; i++)
  cout << "Is this working " << endl;
   for(int i=0; i<maxtimesteps; i++)
  // for(int i=0; i<learning_time_test; i++)
  {

      if(i>=(wash_out_time+learning_time))
      {
    //    outputsignal << Output(i-wash_out_time-learning_time);
    //    outputsignal << endl;

    //    targetsignal << Target_Signal.at(i);
    //    targetsignal << endl;
        //(i-wash_out_time-learning_time) = Target_Signal.at(i);
        Test_Data.push_back(Target_Signal.at(i));
        Output_Signal.push_back(Output(i-wash_out_time-learning_time));
      }

  }
  Mean_squared_error = MSE(Output_Signal, Test_Data);

  cout <<"The mean squared error of the output signal versus the target signal is: " << Mean_squared_error;
  cout <<endl;
  return Mean_squared_error;
}

void Simulation::output_Output_Signal(string& s)
{
  string str = s + "_" + "outputsignal.csv";
  ofstream learningweights("learningweights.csv"); learningweights.precision(15);
  ofstream targetsignal("targetsignal.csv");  targetsignal.precision(15);
  ofstream targetsignal2("targetsignal2.csv");  targetsignal.precision(15);
  ofstream targetsignal3("targetsignal3.csv");  targetsignal.precision(15);

  ofstream learningmatrix("learningmatrix.csv");  learningmatrix.precision(15);
  ofstream learningmatrix2("learningmatrix2.csv");  learningmatrix.precision(15);
  ofstream learningmatrix3("learningmatrix3.csv");  learningmatrix.precision(15);
  ofstream outputsignal(str);  learningmatrix.precision(15);

  ofstream inputsignalcheck("inputsignalcheck.csv");  inputsignalcheck.precision(15);

  ofstream chaoscheck(str);

  for(int i=0; i<maxtimesteps; i++)
 // for(int i=0; i<learning_time_test; i++)
 {

     //Use number_of_threads_or_webs isntead of s.size()
     //for(int j = 0; j<s.size(); j++)
     for(int j = 0; j<number_of_threads_or_webs; j++)
     {
       learningmatrix << LM(i, j) <<",";
       //cout << LM(i,j) << endl;
       if(i>=wash_out_time && i<(wash_out_time+learning_time)) learningmatrix2 << LM2(i-wash_out_time, j) <<",";
       if(i>=(wash_out_time+learning_time)) learningmatrix3 << LM3(i-wash_out_time-learning_time, j) << ",";
     }

     learningmatrix << endl;
     if(i>=wash_out_time && i<(wash_out_time+learning_time)) learningmatrix2 << endl;
     if(i>=wash_out_time && i<(wash_out_time+learning_time)) targetsignal2 << Target_Signal.at(i) <<",";
     if(i>=(wash_out_time+learning_time)) learningmatrix3 << endl;


     if(i>=(wash_out_time+learning_time))
     {
       outputsignal << Output(i-wash_out_time-learning_time);
       outputsignal << endl;

       targetsignal << Target_Signal.at(i);
       targetsignal << endl;
       //(i-wash_out_time-learning_time) = Target_Signal.at(i);
     }

 }

}


void Simulation::output_Learning_Matrix_CSVFile()
{
	ofstream learningmatrix("learningmatrix.csv");  learningmatrix.precision(15);
    ofstream learningmatrix2("learningmatrix2.csv");  learningmatrix.precision(15);
    ofstream learningmatrix3("learningmatrix3.csv");  learningmatrix.precision(15);
	
	
	for(int i=0; i<maxtimesteps; i++)
 // for(int i=0; i<learning_time_test; i++)
 {

     //Use number_of_threads_or_webs isntead of s.size()
    for(int j = 0; j<s.size(); j++)
    // for(int j = 0; j<number_of_threads_or_webs; j++)
     {
       learningmatrix << LM(i, j) <<",";
       //cout << LM(i,j) << endl;
       if(i>=wash_out_time && i<(wash_out_time+learning_time)) learningmatrix2 << LM2(i-wash_out_time, j) <<",";
       if(i>=(wash_out_time+learning_time)) learningmatrix3 << LM3(i-wash_out_time-learning_time, j) << ",";
     }

     learningmatrix << endl;
     if(i>=wash_out_time && i<(wash_out_time+learning_time)) learningmatrix2 << endl;
     if(i>=(wash_out_time+learning_time)) learningmatrix3 << endl;
 }

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


double Simulation::Rand_In_Range_Exp(double min, double max)
{
  double log10min = log10(min);
  double log10max = log10(max);
  double return_value = ((log10max-log10min)*Uniform(0,1))+log10min;
  return pow(10, return_value);
}

double Simulation::Treble_Sine_Function(double f1, double f2, double f3, double dt, double t, double T)
{
	return sin((2*M_PI*f1*dt*t)/T)*sin((2*M_PI*f2*dt*t)/T)*sin((2*M_PI*f3*dt*t)/T);
}


int Simulation::Random_Input_Nodes(int N)
{
  return N*((int) rand() / (RAND_MAX));
}


double Simulation::generate_Gaussian_Noise(double mu, double sigma)
{
	static const double epsilon = std::numeric_limits<double>::min();
	static const double two_pi = 2.0*3.14159265358979323846;

	thread_local double z1;
	thread_local bool generate;
	generate = !generate;

	if (!generate)
	   return z1 * sigma + mu;

	double u1, u2;
	do
	 {
	   u1 = Uniform(0, 1);
	   u2 = Uniform(0, 1);
	 }
	while ( u1 <= epsilon );

	double z0;
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
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

void Simulation::Initialize_Springs(InitialDataValues &data)
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

  double dist = 0;
  double dist2 = 0;
  double perp_dist_new = 0;
  //maximum possible distance for the initial value as this is a minimisation problem
  double perp_dist_old =largest_x_position + largest_y_position;
  double x2 = 0;
  double y2 = 0;
  int connect_node;
  int connect_node2;

  int arraysubscript1=0;
  int arraysubscript2=0;

    ofstream k1output("k1output.csv");  // Todo: ofs3, etc. is really bad!
    ofstream d1output("d1output.csv");
    ofstream k3output("k3output.csv");
    ofstream d3output("d3output.csv");
    ofstream originallengthoutput("originallengthoutput.csv");

    vector<int> node_list;
    vector<int> unconnected_nodes;


  for(int i=0; i<EdgeList.size(); i++)
  {
      //These take the arraysubscripts and disregard the first four points

      arraysubscript1 = EdgeList[i].at(0) - 4;
      arraysubscript2 = EdgeList[i].at(1) - 4;

      x0 = n[arraysubscript1].get_x_Position();
      x1 = n[arraysubscript2].get_x_Position();
      y0 = n[arraysubscript1].get_y_Position();
      y1 = n[arraysubscript2].get_y_Position();

      //These spring and damping coefficients are not giving different values
      k1 = Rand_In_Range_Exp(data.min_k1, data.max_k1);
      d1 = Rand_In_Range_Exp(data.min_d1, data.max_d1);
      k3 = Rand_In_Range_Exp(data.min_k3, data.max_d3);
      d3 = Rand_In_Range_Exp(data.min_d3, data.max_d3);

      //Just for Test
    //  d3 = 0;
  //    k3 = 0;



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
      node_list.push_back(s[i].Nodea());
      node_list.push_back(s[i].Nodeb());
    //  Sort()
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
      tri.push_back(s1.str());
      tri.at(k) = tri.at(k).substr(1, tri.at(k).size()-3);

      s1.str("");
      istringstream iss(tri.at(k));

          iss >> node1;
      //    cout <<"node1: " <<node1 <<endl;
          iss >>sep;
          iss >>node2;
        //  cout <<"node2: "<<node2 << endl;
          iss >>sep;
          iss >>node3;

          iss >>sep;


        //  cout <<" node1: " << node1 << " node2: " <<node2 << " node3: " << node3 << endl;
          Sort(node1, node2, node3);
        //  cout <<"node1: " << node1 << " node2: " <<node2 << " node3: " << node3 << endl;


          if(node1>3)
         {
            NodeList.push_back(node1);
            NodeList.push_back(node2);
            EdgeList.push_back(NodeList);
            NodeList.clear();

            NodeList.push_back(node2);
            NodeList.push_back(node3);
            EdgeList.push_back(NodeList);
            NodeList.clear();
          }

          else if(node2>3)
          {
            NodeList.push_back(node2);
            NodeList.push_back(node3);
            EdgeList.push_back(NodeList);
            NodeList.clear();
          }

          k++;
       }

      Remove_Duplicates(EdgeList);

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


   str = "X.csv";
   ofstream nodesX(str);


   str = "Y.csv";
   ofstream nodesY(str);

  int j=0;
  while(j<n.size())
  {
    nodesX << n[j].get_x_Position();
    if(j<n.size()-1) nodesX<<",";
    nodesY << n[j].get_y_Position();
    if(j<n.size()-1) nodesY<<",";
    j++;
  }

// }

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

bool Simulation::Stability_return()
{
	return stability;
}

double Simulation::return_k1()
{
	return k1_identical;
}

double Simulation::return_k3()
{
	return k3_identical;
}

double Simulation::return_d1()
{
	return d1_identical;
}

double Simulation::return_d3()
{
	return d3_identical;
}
/*
Spider_Web_Simulation::Spider_Web_Simulation(double radius, int rounds, int no_of_points_per_round)
{
	this->radius = radius;
	this->rounds = rounds;
	this->no_of_points_per_round = no_of_points_per_round;

  update(true);
  Mean_Squared_Error = output_LearningMatrix_and_MeanSquaredError();
  
  
}
/*
void Spider_Web_Simulation::Initialize_Nodes()
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

     double random_factor_x = 0;
     double random_factor_y = 0;

     bool odd_even_check = 1;


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
         
         
         //change mass of each node
         node.change_Mass(data.mass_of_nodes);

         n.push_back(node);
         


         if(i>0)
         {

         k1 = Rand_In_Range_Exp(data.min_k1, data.max_k1);
         d1 = Rand_In_Range_Exp(data.min_d1, data.max_d1);

         cout << k1 << endl;
         cout << d1 << endl;

         k3 = Rand_In_Range_Exp(data.min_k3, data.max_k3);
         d3 = Rand_In_Range_Exp(data.min_d3, data.max_d3);


      cout << k3 << endl;
      cout << d3 << endl;


         l0 = Eucl_Dist(x0, y0, x_position, y_position);
         wout = 0;

         s.push_back(Springs(k1, d1, k3, d3, l0, k, k-1, wout));
         }

         //For radial pattern.
         if(j>0)
         {

         k1 = Rand_In_Range_Exp(data.min_k1, data.max_k1);
         d1 = Rand_In_Range_Exp(data.min_d1, data.max_d1);
         k3 = Rand_In_Range_Exp(data.min_k3, data.max_d3);
         d3 = Rand_In_Range_Exp(data.min_d3, data.max_d3);
         cout << "d3 is: " << d3 << endl;

      //   k3 = 0;
      //   d3 = 0;

      //   k1 = data.min_k1;
      //   d1 = data.min_d1;


         //k3 = 0;
         //d3 = 0;
        // k3 = data.max_k3;
    //     d3 = data.max_d3;

         //l0 = Eucl_Dist(x0, y0, x_position, y_position);
         l0 = radius;
       //  wout = Uniform(data.w_out_initial, data.w_out_final);
          wout = 0;

/*
         if(odd_even_check)
         {
         if(j%2==0 && i%2==0) s.push_back(Springs(k1, d1, k3, d3, l0, k, k-no_of_points_per_round, wout));
         if(j%2==1 && i%2==1) s.push_back(Springs(k1, d1, k3, d3, l0, k, k-no_of_points_per_round, wout));
         odd_even_check = false;
         }

         else
         {
         if(j%2==0 && i%2==1) s.push_back(Springs(k1, d1, k3, d3, l0, k, k-no_of_points_per_round, wout));
         if(j%2==1 && i%2==0) s.push_back(Springs(k1, d1, k3, d3, l0, k, k-no_of_points_per_round, wout));
      //   odd_even_check = true;
         }


*/
/*
         s.push_back(Springs(k1, d1, k3, d3, l0, k, k-no_of_points_per_round, wout));



         }

         x0 = x_position;
         y0 = y_position;
         k++;

      }
      x_position = (j+1)*radius*cos((0));
      y_position = (j+1)*radius*sin((0));

      k1 = Rand_In_Range_Exp(data.min_k1, data.max_k1);
      cout <<"k1 is: " <<k1 << endl;
      d1 = Rand_In_Range_Exp(data.min_d1, data.max_d1);
      cout <<"d1 is: " <<d1 << endl;
      k3 = Rand_In_Range_Exp(data.min_k3, data.max_k3);
      cout <<"k3 is: " <<k3 << endl;
      d3 = Rand_In_Range_Exp(data.min_d3, data.max_d3);
      cout <<"d3 is: " <<d3 << endl;

    //  k3 = 0;
    //  d3 = 0;


      l0 = Eucl_Dist(x0, y0, x_position, y_position);
      cout <<"l0 is: " <<l0 << endl;
     // wout = Uniform(data.w_out_initial, data.w_out_final);
     //Temporarily remove
     wout = 0;

      if(no_of_points_per_round>2)
      {
      s.push_back(Springs(k1, d1, k3, d3, l0, (k-no_of_points_per_round), k-1, wout));
      }
    }

    vector<int> nos;

    for(int i=0; i<n.size()-1; i++)
    {
    nos.push_back(i);
    random_shuffle(nos.begin(), nos.end());

    }


    int N = rounds * no_of_points_per_round + 1;

    int input_node_nums = 0;

    input_node_nums = data.input_connectivity*(N);

    cout << input_node_nums << endl;

    int randomnum = 0;
    double win;

    cout << "No of input nodes is: " << input_node_nums << endl;
    cout << "No of input nodes is: " << input_node_nums << endl;



    n[no_of_points_per_round*(rounds-1)].set_Fixed_Node();
    n[(no_of_points_per_round/2)+no_of_points_per_round*(rounds-1)].set_Fixed_Node();


    cout << "Outer fixed node 1 " << no_of_points_per_round*(rounds-1) << endl;
    cout << "Outer fixed node 2 " <<(no_of_points_per_round/2) + no_of_points_per_round*(rounds-1) << endl;




//Temporarily, input nodes will be fixed for this spider web.


//This fixed the outer ring of the nodes. Will get user access for this.
    for(int i=0; i<N; i++)
     {
      //Input weights for the number of input_connectivitiy nodes.
  //    win = Uniform(data.min_input_weight, data.max_input_weight);
    //  randomnum = (int)Uniform(0, N);

//      randomnum = (int)Uniform(0, no_of_points_per_round*(rounds-1));

      //Just

//      if(i<input_node_nums)
//      {
//      n[randomnum].init_Input_Node(data.ux, data.uy, win);
//      cout <<"The input node is:" << randomnum << endl;
//      }



      //Make sure fixed nodes are not input nodes

//      if(i<input_node_nums)
//      {
      //If it is not a fixed node.
      /*
      while(n[randomnum].is_Fixed_Node())
        {
        //C++ typecasting rounds down (truncates) but this is fine going from 0 to N-1.
        randomnum = (int)Uniform(0, N);
        }

      n[randomnum].init_Input_Node(data.ux, data.uy, win);
      */

    //  }
// for non fixed, get rid of this.

   //   if(i>=no_of_points_per_round*(rounds-1) && i<((no_of_points_per_round)+no_of_points_per_round*(rounds-1)))
   //   {
  //      n[i].set_Fixed_Node();
  //      cout <<"The: " <<i <<"th" << " fixed node is fixed." << endl;
  //    }
    

    //        if(i<no_of_points_per_round*(rounds-1) )
      //      {

    //        randomnum = (int)Uniform(0, no_of_points_per_round*(rounds-1));
    //        n[randomnum].init_Input_Node(data.ux, data.uy, win);
         //If it is not a fixed node.
            /*
            while(n[randomnum].is_Fixed_Node())
              {
              //C++ typecasting rounds down (truncates) but this is fine going from 0 to N-1.
              randomnum = (int)Uniform(0, N);
              }

            n[randomnum].init_Input_Node(data.ux, data.uy, win);
            */

      //       }


//     }

     // add in central node here.

//     cout << randomnum << endl;

/*


     Nodes central_node(0, 0);
     n.push_back(central_node);

     win = Uniform(data.min_input_weight, data.max_input_weight);
    // n[5].init_Input_Node(data.ux, data.uy, win);
    //for another round
     n[10].init_Input_Node(data.ux, data.uy, win);

     cout << "number of nodes is: " << n.size() << endl;

    // win = Uniform(data.min_input_weight, data.max_input_weight);
  //   n[16].init_Input_Node(data.ux, data.uy, win);




     for(int i=0; i<no_of_points_per_round; i++)
     {
       k1 = Rand_In_Range_Exp(data.min_k1, data.max_k1);
       cout <<"k1 is: " <<k1 << endl;
       d1 = Rand_In_Range_Exp(data.min_d1, data.max_d1);
       cout <<"d1 is: " <<d1 << endl;
       k3 = Rand_In_Range_Exp(data.min_k3, data.max_k3);
       cout << "k3 is: " <<k3 << endl;
       d3 = Rand_In_Range_Exp(data.min_d3, data.max_d3);
       cout << "d3 is: " <<d3 << endl;

       l0 = Eucl_Dist(0, 0, n[i].get_x_Position(), n[i].get_y_Position());

       wout = 0;

       s.push_back(Springs(k1, d1, k3, d3, l0, n.size()-1, i, wout));

     }
}

*/



