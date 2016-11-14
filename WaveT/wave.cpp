/*-------------------------------------------------------------------------------//
Main Program for finding pressure for imcompressible flows using Poisson equations.
Finds solution at P(1/2,1/2) for Land descretization error for w values of 1 for 
20x20,40x40 and 60x60 array

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

using namespace std;


#define E07 0.0000001
#define E08 0.00000001
#define E09 0.000000001
#define E10 0.0000000001
#define E11 0.00000000001




int main(int argc, char **argv)
{


carray wave1;//my main array
carray analytic;


//set size


//set array size or default used 162x162
set_array_size(wave1, 1200, 1, 1.0);//array, xsize, ysize, dimension
set_array_size(analytic, 1200, 1, 1.0);

//set analytic solution
set_analytic(analytic);

//print_array(analytic);



//set ghost cells as boundary conditions


set_zero(wave1);
set_intial_cond(wave1);//set ghost cells/boundaries
print_array(wave1);//print array in terminal




//---------------------GS SOR w=1.3 loop 1----------------------//

solve_arrayRK2(wave1, 1.0, 0.0001);
print_array(wave1);
get_l2norm(wave1, analytic);
//cout << "Solution: " << get_solution(poisson1) << "\n";


//----------------------Draw Data---------------------//



if(1)//start root application
{
	TApplication theApp("App", &argc, argv);
	draw_graph_wave1(analytic, wave1);//draw 3d graph
	theApp.Run();
}



//end
}



