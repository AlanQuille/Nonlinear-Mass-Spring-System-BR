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
	double pxdot = 0;
	double pydot = 0;

	int connections;

	//spring constants, not approriate for spring constnats.
//	double k;
//	double d;
    //wheter node is connected to input or not.
    bool input_node=0;
    //If it is an input node, than input force u in the x and y directions
    //For the calculations it's not relevant.
    double ux;
    double uy;
    //if it is an input node, than what is the input weight?
    double win;
    //L non crossing
    //Each node has a certain number of springs attached to it. Each spring acts as a force on the node, which changes its position, and hence the position of the spring.
	int noofconnectedsprings;

  public:

  int test = 0;

	Nodes(double x_position, double y_position)
	 {
		this->x_position=x_position;
	  this->y_position=y_position;
   }

 //If this function is activated than the node is input node
	void Input_Node(double ux, double uy, double win)
	{
		this->ux = ux;
		this->uy = uy;
		this->input_node = 1;
		this->win = win;
	}

	void Output_Position()
	{
		cout <<x_position << " " <<y_position <<" " << input_node;
	}

	double X_Position()
	{
		return x_position;
	}

	double Y_Position()
	{
		return y_position;
	}

	//This is the function that incrementally changes the nodes position in the next timestep;
	void Change_Position(double Fx, double  Fy, double  dt)
	{


	    double pxdotdot;
	    double pydotdot;

	    double pxdot;
	    double pydot;

		  pxdotdot= (Fx + win*ux)/m;
      pydotdot= (Fy/m);

        //You want to calculate the velocity so that initial velocity is 0 and than calculate the corresponding change in position
      pydot = EulerMethod(pydotdot, pydotdot, dt);
      this->y_position = EulerMethod(y_position, pydot, dt);

      pxdot = EulerMethod(pxdotdot, pxdotdot, dt);
      this->x_position = EulerMethod(x_position, pxdot, dt);
	}

	double EulerMethod(double x0, double f, double dt)
	{
		double xnew = x0 + dt*f;
		return xnew;
	}



    //whether node is connected to output or not;
    //bool output_node;
    //That doesn't make a difference per se whether the node is an output node.


	//the node has position (x,y) which is determined by a uniform distribution. Keep it simple for the moment and have it determined by a uniform distribution between 0 and 1.
	/*
	void initializeposition()
	{
	    uniform_real_distribution<double> distribution(0.0,1.0);
	    cout <<
	}
	*/
};
