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
    //This is the velocity of the spring relative to l0.
    double x2 = 0;
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

    //Total force on spring
    double x2force =0;
    double p;
    double q;

    public:

    //Default constructor to load in spring and damping coefficeints
    Springs(double k1, double d1, double k3, double d3, double l0, int nodea, int nodeb, double wout);

    //The equation to change the force due to the spring
    void ForceEq(double &Fsum);

    //Change the length and velocity of the spring
    void Change_Length_And_Velocity(double &dt, double &l);

    //Output l0, initial length of spring without extension.
    double Return_Original_Length();

    //Output current length of spring
    void Output();

    //Output first node number
    double Nodea();

    //Output second node number
    double Nodeb();

    //Output the output weight
    double Output_Weight();

    //Output damping and spring coefficients.
    double Outputk1();
    double Outputk3();
    double Outputd1();
    double Outputd3();

};
