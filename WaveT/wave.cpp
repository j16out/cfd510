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
carray analytic1;
carray wave2;//my main array
carray analytic2;
carray wave3;//my main array
carray analytic3;


//set size


//set array size or default used 162x162
set_array_size(wave1, 20, 1, 1.0);//array, xsize, ysize, dimension
set_array_size(wave2, 40, 1, 1.0);
set_array_size(wave3, 80, 1, 1.0);
set_array_size(analytic1, 20, 1, 1.0);
set_array_size(analytic2, 40, 1, 1.0);
set_array_size(analytic3, 80, 1, 1.0);



//print_array(analytic);



//set intial conditions


set_zero(wave1);
set_intial_cond(wave1);
//print_array(wave1);//print array in terminal
set_zero(wave2);
set_intial_cond(wave2);

set_zero(wave3);
set_intial_cond(wave3);




//---------------------solve array1----------------------//
solve_arrayRK2(wave1, 1.0, 0.4);
set_analytic(analytic1, wave1);
get_l2norm(wave1, analytic1);
//cout << "Solution: " << get_solution(poisson1) << "\n";



//---------------------solve array2----------------------//
solve_arrayRK2(wave2, 1.0, 0.4);
set_analytic(analytic2, wave2);
get_l2norm(wave2, analytic2);


//---------------------solve array2----------------------//
solve_arrayRK2(wave3, 1.0, 0.4);
set_analytic(analytic3, wave3);
get_l2norm(wave3, analytic3);



//----------------------Draw Data---------------------//

if(1)//start root application
{
	TApplication theApp("App", &argc, argv);
	draw_graph_wave1(analytic1, wave1);//draw 3d graph
	theApp.Run();
}



//end
}



