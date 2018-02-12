#include <iostream>
#include "Nodes.h"
using namespace std;

Nodes::Nodes(double x_position, double y_position)
 {
  this->px = x_position;
  this->py = y_position;
 }

//If this function is activated than the node is input node
void Nodes::init_Input_Node(double ux, double uy, double win)
{
  this->F_in_x = ux;
  this->F_in_y = uy;
  this->input_node = 1;
  this->win = win;
}

bool Nodes::is_Input_Node()
{
  return input_node;
}

double Nodes::return_Win()
{
  return win;
}


void Nodes::set_Fixed_Node()
{
  this->fixednode = true;
}

void Nodes::print_position()
{
  cout <<"X_position" << px << " " <<"Y_Position" << py <<" " << input_node;
}

double Nodes::get_x_position()
{
  return px;
}

double Nodes::get_y_position()
{
  return py;
}

//This is the function that incrementally changes the nodes position in the next timestep;
void Nodes::Update(double Fx, double  Fy, double  dt)
{
    if(fixednode == false)
    {
    pxdotdot = (Fx + (win*F_in_x))/m;
    pydotdot = (Fy/m);
      //You want to calculate the velocity so that initial velocity is 0 and than calculate the corresponding change in position
    pydot += dt*pydotdot;
    py += dt*pydot;

    pxdot += dt*pxdotdot;
    px += dt*pxdot;
    }
    else
    {
      pxdotdot = 0;
      pxdot = 0;
    }
}
