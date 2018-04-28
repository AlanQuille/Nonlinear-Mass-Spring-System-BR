#include <iostream>
#include "Nodes.h"
using namespace std;

Nodes::Nodes(double x_position, double y_position)
 {
  this->px = x_position;
  this->py = y_position;

  this->original_px = x_position;
  this->original_py = y_position;
 }

//If this function is activated than the node is input node
void Nodes::init_Input_Node(double ux, double uy, double win)
{
  this->F_in_x = ux;
  this->F_in_y = uy;
  this->input_node = 1;
  this->win = win;
  cout <<"win is: " << win;
}

bool Nodes::is_Input_Node()
{
  return input_node;
}

bool Nodes::is_Fixed_Node()
{
  return fixednode;
}

double Nodes::return_Win()
{
  return win;
}



void Nodes::original_positions()
{
  px = original_px;
  py = original_py;
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

double Nodes::get_x_velocity()
{
  return pxdot;
}

double Nodes::get_y_velocity()
{
  return pydot;
}

double Nodes::get_x_acceleration()
{
  return pxdotdot;
}

double Nodes::get_y_acceleration()
{
  return pydotdot;
}


//This is the function that incrementally changes the nodes position in the next timestep;

void Nodes::Input_Force(double Fx, double Fy)
{
  if(fixednode == false)
  {
    F_in_x += Fx;
    F_in_y += Fy;
  }
}

void Nodes::Zero_Force()
{
  F_in_x = 0;
  F_in_y = 0;
}

void Nodes::Update(double dt)
{

  //Standard Euler's
    pxdotdot = (F_in_x)/m;
//    cout <<"acceleration x is: " << pxdotdot << endl;
    pydotdot = (F_in_y)/m;
  //  cout <<"acceleration y is: " << pydotdot << endl;
      //You want to calculate the velocity so that initial velocity is 0 and than calculate the corresponding change in position
    pydot += dt*pydotdot;
  //  cout <<"velocity y is: " << pydot << endl;
    py += dt*pydot;
  //  cout <<"position y is: " << py << endl;

    pxdot += dt*pxdotdot;
//    cout <<"velocity x is: " << pxdot << endl;
    px += dt*pxdot;
//    cout <<"position x is: " <<px << endl;



/*
  double k1x, k1v, k2x, k2v, k3x, k3v, k4x, k4v;

  k1x = pxdot;
//  k1v = feval(func, t(i)    , x(i)         , v(i)        ,F );
  k1v = F_in_x_saved;

  k2x = pxdot+k1v*dt/2;
  k2v = (F_in_x_saved+F_in_x)/2;

  k3x = pxdot+k2v*dt/2;
  k3v =(F_in_x_saved+F_in_x)/2;

  k4x = pxdot+k3v*dt;
  k4v = F_in_x;

  px += (k1x + 2*k2x + 2*k3x + k4x)*dt/6;
  pxdot += (k1v + 2*k2v + 2*k3v + k4v)*dt/6;

  F_in_x_saved = F_in_x;
  */


}

void Nodes::change_updatecheck()
{
    updatecheck = 1;
}

bool Nodes::return_updatecheck()
{
    return updatecheck;
}

void Nodes::zero_updatecheck()
{
   updatecheck = 0;
}
