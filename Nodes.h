#include <iostream>
using namespace std;

class Nodes
{
// Noe mass is 1 by default
private:
	double m = 1;
	// Cartesian coordinates
    // Todo: for consitency rename to px an py (see below → pydot...etc.)
	double px;
	double py;

	// The velocity of the nodes initially it is zero.
    int connections;
    bool input_node = 0;

	// Force in x direction
    // Todo: Wouldn't it make more sense to call if fx and fy
		//I want to keep this one consistent with your code.
		// Fx and Fy are the total forces on the spring due to external forces such as ux, uy, mabye others.
    double ux = 0;

    // Force in y direction
    double uy = 0;

	// Input weight intitialised to 0
    // Weight to scale input to force
    double win = 0;

    // No of connected springs to the node
	int noofconnectedsprings;

    // Node is globally fixed or not
    bool fixednode = false;

	// Acceleration
	double pxdotdot = 0;
	double pydotdot = 0;

    // Velocity
	double pxdot = 0;
	double pydot = 0;

  public:

    // Todo: Naming? What does that

    //Default constructor, which position the node is in
	Nodes(double x_position, double y_position);

    // If this function is activated than the node is input node
		//This function is activated when initialising the nodes. The input force and the input weights are input
    // Todo: How is this function acitvated???
	void Input_Node(double ux, double uy, double win);

    // Return node
    // Todo: Better would be to call it something like "is_Input_Node()"
	bool is_Input_Node();

    // Return win
    // Todo: How is this a boolean? The weight should be a doubel or real value
	int Return_Win();

    //Whether the node is fixed or not. In the reservoir computer the left and right node are fixed
    // Todo:Better "isFixedNode" or is_Fixed_Node" - be consistent Alan!
	void is_Fixed_Node();

    //Show current position
    // Todo: Should it be get_Position()
	void get_position();

    //Show xposition
    // Todo: Same here: get_x_Position()
	double get_x_position();

    //Show yposition
    // Todo: detto
	double get_y_position();

	//This is the function that incrementally changes the nodes position in the next timestep;
    // maybe instead of "change" use "update"
	void Update(double Fx, double  Fy, double  dt);
};
