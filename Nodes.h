#include <iostream>
using namespace std;

class Nodes
{
	//the node has mass M=1
private:
	double m =1;
	//Cartestian coordinates
	double x_position;
	double y_position;

	//The velocity of the nodes initially it is zero.

	int connections;
  bool input_node=0;

	//Force in x direction
  double ux =0;
	//Force in y direction
  double uy = 0;

	//Input weight
  double win = 0;

	int noofconnectedsprings;

	bool fixednode = 0;

	//Acceleration
	double pxdotdot = 0;
	double pydotdot = 0;

  //Velocity
	double pxdot = 0;
	double pydot = 0;

  public:

  int test = 0;

  //Default constructor, which position the node is in
	Nodes(double x_position, double y_position);

 //If this function is activated than the node is input node
	void Input_Node(double ux, double uy, double win);

  //Whether the node is fixed or not. In the reservoir computer the left and right node are fixed
	void FixedNode();

  //Show current position
	void Output_Position();

  //Show xposition
	double X_Position();

  //Show yposition
	double Y_Position();

	//This is the function that incrementally changes the nodes position in the next timestep;
	void Change_Position(double Fx, double  Fy, double  dt);
};
