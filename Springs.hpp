#include <iostream>
using namespace std;


class Springs
{

  //natural length of each spring, i.e., distance between mass nodes as time 0. This might not be necessary to include.
  private:

	double l0;
	double l;
	//l is the current length.
	double x1;
	//l-l0 is the x1 at time(t)
	double x2;
    //spring constants, not approriate for spring constnats.
	double k1;
	double d1;
	double k3;
	double d3;

	//which nodes are each spring connected to
	int nodea;
	int nodeb;

	//the weight on each spring.
	double wout;

	public:

	Springs(double k1, double d1, double k3, double d3, double l0, double nodea, double nodeb, double wout)
	{
		this->l0=l0;
		this->k1=k1;
	  this->k3=k3;
		this->d1=d1;
		this->d3=d3;
		this->nodea = nodea;
		this->nodeb = nodeb;
		this->wout = wout;

		//The initial x1 = l0 - l0 =0, and initial spring velocity = 0
		this->x1 = 0;
		this->x2 = 0;
	};

	void ForceEq(double &Fsum)
	{
		double p;
		double q;
		p=k3*x1*x1*x1 + k1*x1;
        q=d3*x2*x2*x2 + d1*x2;
        Fsum = -p-q;
	};

	void Change_Length_And_Velocity(double &l, double &dt)
	{
		x1 = l - l0;
		double x2force = -(k3*x1*x1*x1 + k1*x1) - (d3*x2*x2*x2 + d1*x2);
		cout << x2force << endl;

		x2 += dt*x2force;
	}


	void Output()
	{
		cout <<"Node 1 is: " << nodea << endl;
		cout <<"Node 2 is: " << nodeb << endl;
		cout <<"Length of spring at rest is: " << l0 << endl;
	}

	double Nodea()
	{
		return nodea;
	}

	double Nodeb()
	{
		return nodeb;
	}

//	void InputNextX1()
//	{
//		x1
//	}
};
