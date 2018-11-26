#include <iostream>
using namespace std;


class Springs
{

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

    public:

        // Default constructor to load in spring and damping coefficeints
        Springs(double k1, double k3, double d1, double d3, double l0, int nodea, int nodeb, double wout);

        //The equation to change the force due to the spring
        void get_Force(double &Fsum);

        // get resting length l0
        double return_Initial_Length();

        //return all springs to the original length l0.
        void set_Original_Length();

        //Make total force 0 for reset.
        void set_Force_0();

        // Output current length of spring
        void print_output();

        //Output first node number
        int Nodea();

        //Output second node number
        int Nodeb();

        //Output the output weight
        double get_Output_Weight();

        //Output damping and spring coefficients.
        double get_k1();
        double get_k3();
        double get_d1();
        double get_d3();

        //These are for changing the k1, k3, d1 and d3.
        void set_k1(double k1);
        void set_k3(double k3);
        void set_d1(double d1);
        void set_d3(double d3);

        //Return x1 and x2
        double return_x1();
        double return_x2();

        //Set x1 and x2
        void set_x1(double x1);
        void set_x2(double x2);

};
