#include <vector>
#include "Eigen/Dense"
#include "Delaunay.hpp"
#include "Springs.cpp"
#include "Nodes.cpp"
using namespace Eigen;


struct InitialDataValues
{
    int     N;   // number of Nodes
    double ux;   // First input values in x direction  TODO: Really needed?
    double uy;   // First input values in y direction

    double input_connectivity_percentage;  // [0,1] percentage of nodes that receive input
    // lower and upper range for input weights
    double min_input_weight;
    double max_input_weight;

    // Parameters to set area where nodes can be placed in
    // Todo: maybe change to min_x... and max... also below with probabilties
    double min_x_position;     // range0x ?? name is not very descriptive (same for the others below)
    double max_x_position;
    double min_y_position;
    double max_y_position;

    // lower and upper limits for log-uniform distribution
    double min_log_uniform;
    double max_log_uniform;

    // lower and upper limits for uniform distribution
    double min_uniform;
    double max_uniform;

    double t0;      // time for first time step [s]
    double tmax;    // maximum time step [s]
    double dt;      // time step in seconds

};

// Todo: You make a class DynamicalSystem or class DynSysData, but then you have functions like Volterra in there
// I would have thought Voleterra would be a derived class from DynamicalSystem
// Also, is the name DynamicalSystem descriptive? At the end it seems to be only a container of data
// Would it be better called DataSet or something along these lines?

class DataSet
{
    private:
        double t0;
        double dt;
        double tmax;
        int maxtimesteps;

    public:
        //This loads in the initial values for the signal, t0, tmax and dt
        DataSet(double t0, double tmax, double dt);

        //A simple sinewave to test target signal;
        void SineWave(vector<double> &SineWave);
};



class Simulation
{

    private:
        int N;                                  // Number of mass points
        double input_connectivity_percentage;
        int num_input_nodes;    // Number of input nodes
        vector<Nodes> n;        // List of all nodes
        vector<Springs> s;      // List of all springs
        vector< vector<double> > EdgeList;   // Todo: Is that really part of Simulation class -
        vector<double> NodeList;             // Todo: Is that really part of Simulation class -
        vector<double> EdgeNodeList;         // Todo: Is that really part of Simulation class -

        double log_uniform_smallest_value;    // ??
		    double log_uniform_largest_value;     // ??
        double uniform_smallest_value;      // ??
	    	double uniform_largest_value;          // ??

		    double input_weight_smallest_value;
	    	double input_weight_largest_value;

        double t0;
        double tmax;
        double dt;

        //These time variables are for the washout, learning phase and test data for the weights calcualated in learning phase.
        int wash_out_time;
        int learning_time;
        int learning_time_test;

        int maxtimesteps;

        bool input_node;

        //Ranges for the delaunay triangulation
        double smallest_x_position;     // range0x ?? name is not very descriptive (same for the others below)
	    	double largest_x_position;
	    	double smallest_y_position;
	     	double largest_y_position;

        //These are constant horizontal forces on the input nodes. If this changes so each input node receives a unique force we will have to modify the code
        // Todo: Better to call them Fx and Fy
        double ux;
        double uy;

        //This is for the first phase of reservoir computing
        //bool is_learning_phase_over = false;

        //This is the final function which SHOULD be like the target signal.
        vector<double> Output_Signal;

        //This is the learning_weights vector
        vector<double> Learning_Weights;

        //This is the target signal for learning
        vector<double> Target_Signal;

        //Input signal for system to be simulated.
        vector<double> Input_Signal;

        //LearningMatrix for learning weight multiplication
        // For collecting data for learning
        MatrixXd LM;

    public:

        //Default constructor
        Simulation(InitialDataValues &data, vector<double> &Input_Signal, vector<double> &Target_Signal, int wash_out_time, int learning_time, int learning_time_test);

        //This is an overloaded default constructor. This is not randomly initialized mass spring system, this is a determined one.
        // Todo: Maybe derive a class for spiderweb simulation
        Simulation(double radius, int rounds, int no_of_points_per_round, InitialDataValues &data, vector<double> &Input_Signal, vector<double> &Target_Signal);

        //This creates the nodes for the reservoir computer implementation
        // Todo: Maybe derive a class for spiderweb simulation
        void Initialize_Nodes(double radius, int rounds, int no_of_points_per_round, InitialDataValues &data);

        //This initializes the nodes and puts in appropriate values for the ranges and the weights
        // Todo: change variable names here
        void Initialize_Nodes(double range0x, double range1x, double range0y, double range1y);

        //This changes position of springs and nodes dynamically in time.
        // Todo: Name is not ideal. Better would be to call it update() or similar
        void execute();


        //This does the delaunay triangulation for the two dimensional case and creates the springs for the reservoir computer, not the radial spider web
        void Delaunay_Triangulation_and_Spring_Creation();

        //Create EdgeNodeList, defunct function don't want to get rid of it.
        void Create_EdgeNodeList();

        //After delaunay triangulation, how many connecting edges hence springs.
        // Todo: Should the ouput be really double and not unsigned int?
        unsigned int Output_No_of_Edges();

        //This outputs the spring and node positions for the reservoir computer
        void Output_Spring_And_Node_Positions();

        //Get the triangles from the Delaunay Triangulation.
        void Get_Triangles(DelaunayTriangulation &Delaunay);

        //Create springs for reservoir computer nonlinear mass spring system
        void Initialize_Springs();

        //Does the Moore-Penrose pseudoinverse from Eigen library
        void Moore_Penrose_Pseudoinverse(MatrixXd& L);


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
        void Remove_Duplicates(vector<vector<double>> &x);

        //Return number of springs
        // Todo: Does it return number of springs (int) or an object of class spring?
        // Thi
        Springs Spring_Return(int i);

        //Return number of nodes
        //Todp: Same question here
        Nodes Node_Return(int i);

        //Return number of edges from the triangle.
        // Todo: Is there really a double neeed and not an unsigend int?
        unsigned int Spring_List();
};
