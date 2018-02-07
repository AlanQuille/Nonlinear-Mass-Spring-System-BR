#include <iostream>
using namespace std;


class Springs
{

  //natural length of each spring, i.e., distance between mass nodes as time 0. This might not be necessary to include.
    private:

        double l0;  // resting lenght
        double l;   // current length.
        double x1;  // difference between current and resting lenght, i.e., l-l0 i
        double x2 = 0; // Velocity of the spring relative to l0 (= dx1/dt)
    
        // spring constants
        // force from the spring stiffness based on x1 is
        // p = x1^3*k3 + x1*k1
        double k1;
        double k3;
        // damping constants
        // force from the damping based on x2 (change of length over time)
        // q = x2^3*d3 + x3*d1
        double d1;
        double d3;

        // Two nodes that are connected by this spring
        int nodea;
        int nodeb;

        // output weight assigned to this spring
        double wout;

        // Total force on spring
        // Todo: Needs to be explained better
        double x2force = 0;
        double p;
        double q;

    
    public:

        // Default constructor to load in spring and damping coefficeints
        Springs(double k1, double k3, double d1, double d3, double l0, int nodea, int nodeb, double wout);

        //The equation to change the force due to the spring
        // Todo: Better name is needed, what does it do acutally? get_force? update_force, etc..
        void ForceEq(double &Fsum);

        //Change the length and velocity of the spring
        // Todo: Maybe update_spring_state would be better
        void Change_Length_And_Velocity(double &dt, double &l);

        //Output l0, initial length of spring without extension.
        // TODO: Why is it called "original" length
        double Return_Original_Length();

        // Output current length of spring
        void Output();

        //Output first node number
        // Todo: is a double needed? Would be an integer enough?
        double Nodea();

        //Output second node number
        // Todo: is a double needed? Would be an integer enough?
        double Nodeb();

        //Output the output weight
        double Output_Weight();

        //Output damping and spring coefficients.
        // Todo: get_k1, get_k3 etc..
        double Outputk1();
        double Outputk3();
        double Outputd1();
        double Outputd3();

};
