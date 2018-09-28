#include <iostream>
using namespace std;


class Nodes
{

private:
    // mass is 1 by default
    double m = 1;
	//double m = 0.12;

		//Original positions
		double original_px;
		double original_py;

	// Cartesian coordinates
	  double px;
	  double py;


    // Velocity
    double pxdot = 0;
    double pydot = 0;

    // Acceleration
    double pxdotdot = 0;
    double pydotdot = 0;

	// The velocity of the nodes initially it is zero.
    // int connections;  // TODO: Still needed
    bool input_node = false;

	// Input force applied to the Node
    double F_in_x = 0;
    // Force in y direction
    double F_in_y = 0;

		//Saved F_in_x
		double F_in_x_saved = 0;
		//Saved F_in_y
		double F_in_y_saved = 0;
    // Weight to scale input to force
    // intitialised to 0
    double win = 0;

    // Node is globally fixed or not
    bool fixednode = false;

		//A bool check to determine whether a node has been updated or not.
		bool updatecheck = 0;


  public:

  //Default constructor, with caretesian coordinates
  Nodes(double x_position, double y_position);

  // Function to convert Node to node that receives and an input force
	void init_Input_Node(double ux, double uy, double win);

  // Checks if node is input node
	bool is_Input_Node();

  //Returns fixed node
	bool is_Fixed_Node();

  //Save original position in x and y of nodes.
	void original_positions();

    // Return input weight w_in for the node
	double return_Win();

    // Make node globally fixed
	void set_Fixed_Node();

    //Show current position
	void print_position();

    //Show xposition
    // Todo: Same here: get_x_Position()
	double get_x_position();

    //Show yposition
    // Todo: detto
	double get_y_position();

	double get_x_velocity();
	double get_y_velocity();

	double get_x_acceleration();
	double get_y_acceleration();

	//This is the function that incrementally changes the nodes position in the next timestep;
    // maybe instead of "change" use "update"
	void Input_Force(double Fx, double Fy);

	void Zero_Force();

	void Update(double dt);

	void zero_Accel_and_Vel();

  //At every timestep, a node should be changed only once. A node should not be changed at every tiem
	void change_updatecheck();

	//Check to see if the node has been modified or not
	bool return_updatecheck();

	//Check to see if the node has been modified or not
	void zero_updatecheck();
	
	//Change mass of springs
	void change_Mass(double new_m);
	
    //Returns mass of nodes
	double return_Mass();
};
