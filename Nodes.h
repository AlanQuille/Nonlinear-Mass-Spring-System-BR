#include <iostream>
using namespace std;


class Nodes
{

private:
    // mass is 1 by default
    double m = 1;
    
	//Original position of node in Cartesian coordinates
	double original_px;
	double original_py;

	//Current position of node in Cartesian coordinates
	double px;
	double py;


    //Current velocity of node in Cartesian coordinates. It is initially 0.
    double pxdot = 0;
    double pydot = 0;

    //Current acceleration of node in Cartesian coordinates. It is initially 0.
    double pxdotdot = 0;
    double pydotdot = 0;

	//Every node is by default not an input node, it must be changed to true at some stage in the simulation
    bool input_node = false;

	//Input force applied to the node in x direction
    double F_in_x = 0;
    //Input force applied to the node in x direction
    double F_in_y = 0;

    // Weight to scale input force
    // intitialised to 0
    double win = 0;

    //Node is globally fixed or not
    bool fixednode = false;

	//A bool check to determine whether a node has been updated or not.
	bool updatecheck = 0;


  public:

    //Default constructor, with Caretesian coordinates
    Nodes(double x_position, double y_position);

    //Function to make node into an input node and also determines what the horizontal and vertical input force is. It also determines the input weight
	void init_Input_Node(double ux, double uy, double win);

    //Returns true if node is input node, false otherwise.
	bool is_Input_Node();

    //Returns true if node is input node, false otherwise.
	bool is_Fixed_Node();

    //Save original position in x and y of nodes.
	void original_Positions();

    // Return input weight w_in for the node
	double return_Win();

    // Make node globally fixed
	void set_Fixed_Node();

    //Show current position 
	void print_Position();

    //Show x position of node
	double get_x_Position();

    //Show y position of node
	double get_y_Position();

    //Show x and y velocities of node
	double get_x_Velocity();
	double get_y_Velocity();

    //Show x and y accelerations of node
	double get_x_Acceleration();
	double get_y_Acceleration();

	//This is the function that incrementally changes the force on the node and hence acceleration, velocity and position in the next timestep;
    void input_Force(double Fx, double Fy);

    //make force 0.
	void zero_Force();

    //This function applies Euler's method to change the position, velocity and acceleration of the node at every timestep.
	void update(double dt);

    //halt the node, no velocity or acceleration after this is applied.
	void zero_Accel_and_Vel();

	//Change mass of springs to new_m
	void change_Mass(double new_m);
	
    //Returns mass of nodes
	double return_Mass();
};
