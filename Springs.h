#include <iostream>
using namespace std;

// TODO: Maybe it would make sense to make a base class called connection and dervie spring from that
class Springs
{

    //Todo: natural length of each spring, i.e., distance between mass nodes as time 0. This might not be necessary to include.
    private:

        double l0;  // resting length
        double l;   // current length.
        double x1 = 0;  // difference between current and resting lenght, i.e., x1 = l-l0
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
        // Todo: Clean up code
        double F_total = 0;
        double p; // Todo: still needed?
        double q; // Todo: still needed?


    public:

        // Default constructor to load in spring and damping coefficeints
        // Todo: We should make this classe more general - the stiffness and damping functiosn should be overloaded
        // So we can implement variations of that
        Springs(double k1, double k3, double d1, double d3, double l0, int nodea, int nodeb, double wout);

        //The equation to change the force due to the spring
        // Todo: Better name is needed, what does it do acutally? get_force? update_force, etc..

        void get_Force(double &Fsum);


        // This could be overloaded
        // void ForceEq(double &Fsum);


        // Change the length and velocity of the spring
        // Todo: Maybe update_spring_state would be better
      //;  void update_Spring_State(double &dt, double &l);


        // get resting length l0
        double return_Initial_Length();

        void set_original_length();

        //Make total force 0 for reset.
        void set_Force_0();

        // Output current length of spring
        void print_output();

        //Output first node number
        // Todo: is a double needed? Would be an integer enough?
        int Nodea();

        //Output second node number
        // Todo: is a double needed? Would be an integer enough?
        int Nodeb();

        //Output the output weight
        double get_Output_Weight();

        //Output damping and spring coefficients.
        // Todo: get_k1, get_k3 etc.
        double get_k1();
        double get_k3();
        double get_d1();
        double get_d3();

        //Return x1 and x2
        double return_x1();
        double return_x2();

        //Set x1 and x2
        void set_x1(double x1);
        void set_x2(double x2);

};
