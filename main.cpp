#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>
#include "Delaunay.hpp"
#include "Springs.hpp"
#include "Nodes.hpp"
#include <sstream>
#include <set>
#include <fstream>
using namespace std;


int Random_Input_Nodes(int N)
{
	return N*((int) rand() / (RAND_MAX));
}

double Uniform(double M, double N)
{

	return M + (rand() / ( RAND_MAX / (N-M) ) ) ;
}


double Log_Uniform(double initial, double final)
{
	return exp(Uniform(initial, final));

}

double Spring_And_Damping_Coefficient_1(double initial, double final)
{
	return Log_Uniform(initial, final);
}

double Spring_And_Damping_Coefficient_2(double initial, double final)
{
	return Uniform(initial, final);
}

double Eucl_Dist(double x1, double y1, double x2, double y2)
{
	return sqrt((y2-y1)*(y2-y1) + (x2-x1)*(x2-x1));
}

double Angle(double x1, double x2, double y1, double y2)
{
	return atan((y2-y1)/(x2-x1));
}

double X_Comp(double vectorsum, double theta)
{
	return vectorsum*cos(theta);
}

double Y_Comp(double sum, double theta)
{
	return sum*sin(theta);
}

void Sort(int &a, int &b, int &c){
    if(a>b)
	{
        int tmp = a;
        a = b;
        b = tmp;
    }
    if(a>c)
	{
        int tmp = a;
        a=c;
        c = tmp;
    }
    if(b>c)
	{
        int tmp = b;
        b=c;
        c=tmp;
    }
}

void Sort(int &a, int &b)
{
    if(a>b)
	{
        int tmp = a;
        a = b;
        b = tmp;
    }
}

void RemoveDuplicates(vector<vector<double>> &x)
{
	// sort (x.begin(), x.end());
	 int i=1;
	 while(i<x.size())
	 {
	 	if(x[i][1] ==x[i-1][1] && x[i][0] ==x[i-1][0])
	 	{
		  x.erase(x.begin()+(i-1));
		  i-=1;
    	}
      i++;
	 }
}


int main(int argc, char** argv)
{
//	initializeposition();
   //Step 1: The nodes were initialized according to an initial pattern and number of nodes N
   srand (time(NULL));
   //Total number of mass points
   int N = 3;
   // HOrizontal force
   double ux = 1;
   //Vertical force
   double uy = 0;
   Nodes* n = new Nodes[N];
   cout << n[0].test;

   cout <<Random_Input_Nodes(0.2) <<endl;

   //Step 2: The positions of the nodes were initialized and 20% of the nodes are connected to the input.
   double input_connectivity = 0.2;
   double total_input_nodes = input_connectivity * ((double)N);
 //  addinputnodes(&n, N, input_connectivity);

   //First x% of the nodes are given input status tahn shuffled.
   for(int i=0; i<N; i++)
   {
   n[i].Declare_Position_Variables(Uniform(0,1), Uniform(0,1));
   if(i < (int)total_input_nodes) n[i].Input_Nodes(1, ux, uy, Uniform(-1, 1));
   }

   random_shuffle(&n[0], &n[N-1]);
   //Do Delaunay Triangulation
   DelaunayTriangulation DT(1, 1);
   for(int i=0; i<N; i++)
   {
   n[i].Output_Position();
   cout <<endl;
   DT.AddPoint(Point(n[i].X_Position(), n[i].Y_Position()));
   }

    //Delaunay Triangulation
    /*
    DelaunayTriangulation DT(1, 1);
    DT.AddPoint(Point(0.1, 0.41));
    DT.AddPoint(Point(0.15, 0.41));
    DT.AddPoint(Point(0.125, 0.38));
    */
    cout << "Reporting" << endl;

		DT.print();

		cout << endl;


	char i;
	int j;
	stringstream s1;


	int iter=0;
	int total = 0;
	int nodea = 0;
	int nodeb = 0;
	//vector<Springs> s;
	vector<string> tri;
	string es;
	int k =0;

	int node1;
	int node2;
	int node3;
	char sep=',';

	int noofconnectingedges = 0;
	int vertices = 0;

	bool oneortwo;
	bool twoorthree;
	bool oneorthree;

	vector<vector<double>> EdgeList;
	vector<double> NodeList;

//	This takes the nodes of each triangle and parses them so only the relevant nodes a
	for(auto e: DT.triangles)
	{
		//cout <<"Node1:" <<get<0>(e) <<" " << "Node2:"<<" " <<get<1>(e)<< endl;
		s1 << e;
	//	cout << e << endl;
		tri.push_back(s1.str());
		tri.at(k) = tri.at(k).substr(1, tri.at(k).size()-3);
		cout << tri.at(k) << endl;
//		cout <<endl;
        s1.str("");
        istringstream iss(tri.at(k));

        iss >> node1;
      //  cout <<"node1: " <<node1 <<endl;
        iss >>sep;
        iss >>node2;
      //  cout <<"node2: "<<node2 << endl;
        iss >>sep;
        iss >>node3;
    //    cout <<"node3: " << node3 << endl;
        iss >>sep;

        if(node1!=0 && node1!=1 && node1!=2 && node1!=3) vertices++;
        if(node2!=0 && node2!=1 && node2!=2 && node2!=3) vertices++;
        if(node3!=0 && node3!=1 && node3!=2 && node3!=3) vertices++;

        if(vertices ==2)
        {
		noofconnectingedges++;
	//	cout << "Connect two nodes!: " <<node1 <<" " <<node2 <<" "<< node3 <<" ";
		Sort(node1, node2, node3);
		if(node1!=0 && node1!=1 && node1!=2 && node1!=3)
		{
			NodeList.push_back(node1);
			NodeList.push_back(node2);
			EdgeList.push_back(NodeList);
			NodeList.clear();
	    }
		else if(node2!=0 && node2!=1 && node2!=2 && node2!=3)
		{
			NodeList.push_back(node2);
			NodeList.push_back(node3);
			EdgeList.push_back(NodeList);
			NodeList.clear();
		}

		//Add a spring
	//	Springs = new Spring()
	    }

	    if(vertices ==3)
        {
		noofconnectingedges+=2;
		//cout << "Connect all nodes: " <<node1 <<" " <<node2 <<" "<< node3 <<" ";
		Sort(node1, node2, node3);
		if(node1!=0 && node1!=1 && node1!=2 && node1!=3)
		{
			NodeList.push_back(node1);
			NodeList.push_back(node2);
			NodeList.push_back(node2);
			NodeList.push_back(node3);
			EdgeList.push_back(NodeList);
			NodeList.clear();
		}
		else if(node2!=0 && node2!=1 && node2!=2 && node2!=3)
		{
			NodeList.push_back(node2);
			NodeList.push_back(node3);
			EdgeList.push_back(NodeList);
			NodeList.clear();
		}
		//Add a spring
	    }
        vertices = 0;
		k++;

	}

	//I have each node number. I need to take the ones that connect up to each node. Neglect points 0, 1, 2 and 3.

	sort (EdgeList.begin(), EdgeList.end());
	RemoveDuplicates(EdgeList);

	cout <<"The number of edges for: " <<N << "Mass points is: " << EdgeList.size() << endl;
	//cout <<"The number of unduplicated edges should be " << EdgeList.size() << endl;
  //Spring and damping coefficients
	double k1 = 0;
	double d1 = 0;
	double k3 = 0;
	double d3 = 0;
	double l0 = 0;

	double x0;
	double y0;
	double x1;
	double y1;
	double wout;

	int arraysubscript1=0;
	int arraysubscript2=0;;
//	double l0 = ;

	//for(int i=0; i<)
	vector<Springs> s;

	//For uniform and log normal distribution
	double initialvalue = 0;
	double finalvalue = 1;
	double mean = 0;
	double stdev = 1;

	for(int i=0; i<EdgeList.size(); i++)
	{
		  //These take the arraysubscripts and disregard the first four points
		  arraysubscript1 = EdgeList[i].at(0) - 4;
		  arraysubscript2 = EdgeList[i].at(1) - 4;

		  x0 = n[arraysubscript1].X_Position();
		  x1 = n[arraysubscript2].X_Position();
		  y0 = n[arraysubscript1].Y_Position();
		  y1 = n[arraysubscript2].Y_Position();

		  //These damping and spring coefficients are not working properly, i will have to fix this somehow.
		  k1 = Spring_And_Damping_Coefficient_1(initialvalue, finalvalue);
      d1 = Spring_And_Damping_Coefficient_1(initialvalue, finalvalue);
	    k3 = Spring_And_Damping_Coefficient_2(initialvalue, finalvalue);
	    d3 = Spring_And_Damping_Coefficient_2(initialvalue, finalvalue);

		  l0 = Eucl_Dist(x0, y0, x1, y1);
		  wout = Uniform(-1, 1);

	    s.push_back(Springs(k1, d1, k3, d3, l0, arraysubscript1, arraysubscript2, wout));
	    s[i].Output();
   }

   double tmax = 1;
   double t0 = 0;
   double Fsum =0;
   double Fx =0;
   double Fy =0;
   double theta = 0;
   double dt = 0.001;
   double l = 0;

   ofstream ofs ("test.csv", ofstream::out);

   int maxtimesteps = (int)(tmax/dt);

   nodea =0;
   nodeb = 0;

   for(int i=0; i<maxtimesteps; i++)
   {

		for(int j=0; j<EdgeList.size(); j++)
   	{

   	   s[j].ForceEq(Fsum);

   	   nodea = s[j].Nodea();
   	   nodeb = s[j].Nodeb();

   	   x0 = n[nodea].X_Position();
   	   x1 = n[nodeb].X_Position();
   	   y0 = n[nodea].Y_Position();
   	   y1 = n[nodeb].Y_Position();

   	   theta = Angle(x0, x1, y0, y1);
   	   Fx = X_Comp(Fsum, theta);
   	   Fy = Y_Comp(Fsum, theta);

   	   n[nodea].Change_Position(Fx, Fy, dt);
   	   n[nodeb].Change_Position(Fx, Fy, dt);

   	   x0 = n[nodea].X_Position();
   	   x1 = n[nodeb].X_Position();

   	   y0 = n[nodea].Y_Position();
   	   y1 = n[nodeb].Y_Position();

   	   l = Eucl_Dist(x0, x1, y0, y1);

       s[j].Change_Length_And_Velocity(l, dt);

    }

       ofs << i*dt;
   	   ofs <<",";
   	   ofs << n[nodea].X_Position();
   	   ofs <<endl;
  }

   ofs.close();
	 return 0;
}
