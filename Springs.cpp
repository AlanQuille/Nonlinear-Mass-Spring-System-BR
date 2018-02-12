#include <iostream>
#include "Springs.h"
using namespace std;

Springs::Springs(double k1, double k3, double d1, double d3, double l0, int nodea, int nodeb, double wout)
{
    this->l0=l0;
    this->x1 = 0;
    this->k1=k1;
    this->k3=k3;
    this->d1=d1;
    this->d3=d3;
    this->nodea = nodea;
    this->nodeb = nodeb;
    this->wout = wout;

    cout <<nodea << endl;
    cout <<nodeb << endl;

    //The initial x1 = l0 - l0 =0, and initial spring velocity = 0
};

void Springs::get_Force(double &Fsum)
{
    // Todo: Clean up, simplify
    p=k3*x1*x1*x1 + k1*x1;
    q=d3*x2*x2*x2 + d1*x2;
    F_total = -p-q;
    Fsum = F_total;
};

void Springs::update_Spring_State(double &dt, double &l)
{
    double x1new = l-l0;
    x2 = ((x1new - x1)/dt);
    this->x1 = x1new;  // Todo: Check if this-> is needed
    //RungeKutta2ndOrder(dt, x2);
};



double Springs::return_Initial_Length()
{
    return l0;
};


void Springs::print_output()
{
    cout <<"Node 1 is: " << nodea << endl;
    cout <<"Node 2 is: " << nodeb << endl;
    cout <<"Length of spring at rest is: " << l0 << endl;
}

int Springs::Nodea()
{
    return nodea;
}

int Springs::Nodeb()
{
    return nodeb;
}

double Springs::get_Output_Weight()
{
return wout;
}

double Springs::get_k1()
{
return k1;
}

double Springs::get_k3()
{
return k3;
}

double Springs::get_d1()
{
return d1;
}

double Springs::get_d3()
{
return d3;
}
