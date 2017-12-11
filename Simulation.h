#include <vector>
#include "Eigen/Dense"
#include "Delaunay.hpp"
#include "Springs.cpp"
#include "Nodes.cpp"
using namespace Eigen;


struct InitialDataValues
{
	double N;
	double ux;
	double uy;
	double input_connectivity;
	double w_in_initial;
	double w_in_final;
	double w_out_initial;
	double w_out_final;
	double range0x;
	double range1x;
	double range0y;
	double range1y;
	double initial_log_uniform;
	double final_log_uniform;
	double initial_uniform;
	double final_uniform;
  double t0;
	double tmax;
	double dt;
};

class DynamicalSystems
{
private:
  double t0;
	double dt;
	double tmax;
	int maxtimesteps;
public:
	//This loads in the initial values for the signal, t0, tmax and dt
	DynamicalSystems(double t0, double tmax, double dt);

	//This is the LotkaVolterra System
	void LotkaVolterra(vector<double> &LVx, vector<double> &LVy);

	//A simple sinewave to test target signal;
	void SineWave(vector<double> &SineWave);
};

class Simulation
{

private:
  //No of mass points
  int N;
  double input_connectivity;
  int total_input_nodes;
  vector<Nodes> n;
  vector<Springs> s;
	vector<vector<double>> Loutput;
  vector<vector<double>> EdgeList;
  vector<double> NodeList;
	vector<double> EdgeNodeList;

  double w_out_initial;
  double w_out_final;
  double w_in_initial;
  double w_in_final;
  double initial_log_uniform;
  double final_log_uniform;
  double initial_uniform;
  double final_uniform;
  double t0;
  double tmax;
  double dt;

	int maxtimesteps;

  bool input_node;

  //Ranges for the delaunay triangulation
  double range0x;
  double range0y;
  double range1x;
  double range1y;

  //These are constant horizontal forces on the input nodes. If this changes so each input node receives a unique force we will have to modify the code
	double ux;
	double uy;

  //This is for the first phase of reservoir computing
	bool learning_phase_over = 0;

	//This is the final function which SHOULD be like the target signal.
	vector<double> Output_Signal;

  //This is the learning_weights vector
	vector<double> Learning_Weights;

//This is the target signal for learning
	vector<double> Target_Signal;

	//LearningMatrix for learning weight multiplication
	MatrixXd LM;


public:

  //Default constructor
  Simulation(InitialDataValues &data, vector<double> &Target_Signal);

  //This is an overloaded default constructor. This is not randomly initialized mass spring system, this is a determined one.
	Simulation(double radius, int rounds, int no_of_points_per_round, InitialDataValues &data, vector<double> &Lvx, vector<double> &Lvy);

  //This creates the nodes for the reservoir computer implementation
	void Initialize_Nodes(double radius, int rounds, int no_of_points_per_round, InitialDataValues &data);

	//This initializes the nodes and puts in appropriate values for the ranges and the weights
	void Initialize_Nodes(double range0x, double range1x, double range0y, double range1y);

	//This changes position of springs and nodes dynamically in time.
  void Execute_In_Time();

 //Overloaded Execute_in_Time, will use this for 3D
	void Execute_In_Time_2();

  //This does the delaunay triangulation for the two dimensional case and creates the springs for the reservoir computer, not the radial spider web
	void Delaunay_Triangulation_and_Spring_Creation();

  //Create EdgeNodeList, defunct function don't want to get rid of it.
	void Create_EdgeNodeList();

  //After delaunay triangulation, how many connecting edges hence springs.
	double Output_No_of_Edges();

  //This outputs the spring and node positions for the reservoir computer
	void Output_Spring_And_Node_Positions();

 //Get the triangles from the Delaunay Triangulation.
  void Get_Triangles(DelaunayTriangulation &Delaunay);

  //Create springs for reservoir computer nonlinear mass spring system
  void Initialize_Springs();

  //Does the Moore-Penrose pseudoinverse from Eigen library
	void Moore_Penrose_Pseudoinverse(MatrixXd& L);

	//Sine wave to test target signal.
	double Sine_Wave(double currenttime);

  //Output for Matlab plot
	void Output_For_Plot();

  //Mean Squared Error between vector A and estimator Ahat
	double MSE(vector<double>& A, vector<double>& Ahat);

  //Return learning matrix, MatrixXd defined in Eigen library
	MatrixXd& Return_Learning_Matrix();

  //Return Output Signal after learning phase and the mean squared error.
	void Output_Signal_And_MSE();

  //Return the weights after learning phase as a vector
	vector<double>& Return_Learning_Weights();

	//Populate the learning matrix with weights.
	void Populate_Learning_Weights(MatrixXd& L);

  //Have N random input nodes.
  int Random_Input_Nodes(int N);

  //Randomly chosen number between N and M
  double Uniform(double M, double N);

  //Randomy chosen log to the base 10 uniform number between initial and final values.
	double Log_10_Uniform(double initial, double finalvalue);

  //Expression of log 10 uniform for spring and damping coefficients
  double Spring_And_Damping_Coefficient_1(double initial, double finalvalue);
  double Spring_And_Damping_Coefficient_2(double initial, double finalvalue);

  //Euclidean distance between two points on x y plane, will overload for 3D
  double Eucl_Dist(double x1, double y1, double x2, double y2);

	//Angle for line between two points.
  double Angle(double x0, double x1, double y0, double y1);

  //X component of vector from one point to another
  double X_Comp(double vectorsum, double theta);

  //Y component of vector from one point to another.
  double Y_Comp(double vectorsum, double theta);

	//Sorting three input numbers.
  void Sort(int &a, int &b, int &c);

  //Sorting two input numbers
  void Sort(int &a, int &b);


  //Remove duplicates from two dimensional vector.
  void RemoveDuplicates(vector<vector<double>> &x);

  //Return number of springs
	Springs Spring_Return(int i);

	//Return number of nodes
	Nodes NodeReturn(int i);

  //Return number of edges from the triangle.
	double Spring_List();
};
