/*-------------------------------------------------------------------------------//
Main Program for looking at convergence behavior for Laplace Problem.
Finds maximum change in solution for w values of 1, 1.3, 1.5 for 20x20 array


Jerin Roberts 2016
compiled using g++/gcc version 5.4.0 on Ubuntu 16.04.02 and are available for clone 
via the link provided: url{https://github.com/j16out/
//-------------------------------------------------------------------------------*/

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <math.h> 
#include "TApplication.h"
#include "vroot/root.hpp"
#include "numerical/numerical.hpp"

#define BIG 1000000
#define E06 0.000001
#define E07 0.0000001
#define E08 0.00000001
#define E09 0.000000001



using namespace std;



int main(int argc, char **argv)
{

float diff = 1;
carray gsarraySOR1;//my main arrays
carray gsarraySOR2;
carray gsarraySOR3;
carray analytic1;
carray analytic2;
carray analytic3;

//set size
set_array_size(gsarraySOR1, 20, 20, 1.0);//array, xsize, ysize, dimension
set_array_size(gsarraySOR2, 20, 20, 1.0);
set_array_size(gsarraySOR3, 20, 20, 1.0);

set_array_size(analytic1, 20, 20, 1.0);
set_array_size(analytic2, 20, 20, 1.0);
set_array_size(analytic3, 20, 20, 1.0);

//set analytic solution
set_analytic(analytic1);
set_analytic(analytic2);
set_analytic(analytic3);


//set ghost cells as boundary conditions
set_zero(gsarraySOR1);
set_ghostcells(gsarraySOR1);

set_zero(gsarraySOR2);
set_ghostcells(gsarraySOR2);

set_zero(gsarraySOR3);
set_ghostcells(gsarraySOR3);



//---------------------GS SOR w=1 loop----------------------//

solve_arraySOR(gsarraySOR1, E09, 1.0);
float n = get_l2norm(gsarraySOR1, analytic1);

//---------------------GS SOR w=1.3 loop----------------------//


solve_arraySOR(gsarraySOR2, E09, 1.3);
n = get_l2norm(gsarraySOR2, analytic2);
//---------------------GS SOR w=1.5 loop----------------------//


solve_arraySOR(gsarraySOR3, E09, 1.5);
n = get_l2norm(gsarraySOR3, analytic3);

//----------------------Draw Data---------------------//

if(1)//start root app
{
	TApplication theApp("App", &argc, argv);
	draw_graph_diff3(gsarraySOR1, gsarraySOR2, gsarraySOR3);//Draw data using root
	theApp.Run();
}



//end
}
