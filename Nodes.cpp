#include <iostream>
#include "Nodes.h"
using namespace std;

Nodes::Nodes(double x_position, double y_position)
 {
  this->px = x_position;
  this->py = y_position;
 }

//If this function is activated than the node is input node
void Nodes::Input_Node(double ux, double uy, double win)
{
  this->ux = ux;
  this->uy = uy;
  this->input_node = 1;
  this->win = win;
}

bool Nodes::is_Input_Node()
{
  return input_node;
}

int Nodes::Return_Win()
{
  return win;
}


void Nodes::is_Fixed_Node()
{
  this->fixednode = 1;
}

void Nodes::get_position()
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
    if(fixednode ==0)
    {
    pxdotdot= (Fx + (win*ux))/m;
    pydotdot= (Fy/m);
      //You want to calculate the velocity so that initial velocity is 0 and than calculate the corresponding change in position
    pydot += dt*pydotdot;
    py += dt*pydot;

    pxdot += dt*pxdotdot;
    px += dt*pxdot;
    }
    else
    {
      pxdotdot =0;
      pxdot = 0;
    }
}
