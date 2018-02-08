#include <iostream>
using namespace std;

// TODO: Maybe it would make sense to make a base class called connection and dervie spring from that
class Springs
{

    //Todo: natural length of each spring, i.e., distance between mass nodes as time 0. This might not be necessary to include.
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
        // Todo: We should make this classe more general - the stiffness and damping functiosn should be overloaded
        // So we can implement variations of that
        Springs(double k1, double k3, double d1, double d3, double l0, int nodea, int nodeb, double wout);

        //The equation to change the force due to the spring
        // Todo: Better name is needed, what does it do acutally? get_force? update_force, etc..
<<<<<<< HEAD
        void Update_Force(double &Fsum);
=======
        // This could be overloaded
        void ForceEq(double &Fsum);
>>>>>>> 0c3e96e802ad7dca06ea4b8fca8d37d1c48af347

        //Change the length and velocity of the spring
        // Todo: Maybe update_spring_state would be better
        void Update_Spring_State(double &dt, double &l);

        //Output l0, initial length of spring without extension.
        // TODO: Why is it called "original" length
        //This is
        double Return_Initial_Length();

        // Output current length of spring
        void Output();

        //Output first node number
        // Todo: is a double needed? Would be an integer enough?
        int Nodea();

        //Output second node number
        // Todo: is a double needed? Would be an integer enough?
        int Nodeb();

        //Output the output weight
        double Output_Weight();

        //Output damping and spring coefficients.
        // Todo: get_k1, get_k3 etc..
        double get_k1();
        double get_k3();
        double get_d1();
        double get_d3();

};
