#include <iostream>
#include "Nodes.h"
using namespace std;

Nodes::Nodes(double x_position, double y_position, double z_position)
 {
  this->x_position=x_position;
  this->y_position=y_position;
  this->z_position=z_position;
 }

//If this function is activated than the node is input node
void Nodes::Input_Node(double ux, double uy, double win)
{
  this->ux = ux;
  this->uy = uy;
  this->input_node = 1;
  this->win = win;
}

void Nodes::FixedNode()
{
  fixednode = 1;
}

void Nodes::Output_Position()
{
  cout <<"X_position" << x_position << " " <<"Y_Position" << y_position <<" " << "Z_position" <<" " <<z_position << input_node;
}

double Nodes::X_Position()
{
  return x_position;
}

double Nodes::Y_Position()
{
  return y_position;
}

double Nodes::Z_Position()
{
  return z_position;
}

//This is the function that incrementally changes the nodes position in the next timestep;
//Z Position will be directed by a different function
void Nodes::Change_Position(double Fx, double  Fy, double Fz, double  dt)
{
    if(fixednode ==0)
    {
    pxdotdot= (Fx + (win*ux))/m;
    pydotdot= (Fy/m);
    pzdotdot = (Fz/m);
      //You want to calculate the velocity so that initial velocity is 0 and than calculate the corresponding change in position
    pzdot += dt*pzdotdot;
    z_position += dt*pzdot;

    pydot += dt*pydotdot;
    y_position += dt*pydot;

    pxdot += dt*pxdotdot;
    x_position += dt*pxdot;
    }
    else
    {
      pxdotdot =0;
      pxdot = 0;

      pydotdot = 0;
      pydot = 0;
    }
}
void Nodes::Change_X_Position(double input)
{
  x_position += input;
}

void Nodes::Change_Y_Position(double input)
{
  y_position += input;
}

void Nodes::Change_Z_Position(double input)
{
  z_position += input;
}
