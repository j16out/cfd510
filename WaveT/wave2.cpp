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
set_array_size(wave1, 500, 1, 20.0, 0);//array, xsize, ysize, dimension
set_zero(wave1);
set_intial_cond(wave1);
solve_arrayRK2(wave1, 8.0, 0.1);//def stable

//---------------------plot unstable cfl----------------------//
float cfl = 0.4;
float factor = 0.0025;

while(cfl <= 0.7)
{
carray wave2;
set_array_size(wave2, 500, 1, 20.0, 0);
set_zero(wave2);
set_intial_cond(wave2);
solve_arrayRK2(wave2, 8.0, cfl);
//set_analytic(analytic1, wave1);
float l2 = get_l2norm(wave1, wave2);

if(l2>BIG)
break;
else
{
wave1.l2norm.push_back(l2);
wave1.diff.push_back(cfl);
}
if(l2>0.052)
factor = 0.001;

cfl = cfl + factor;
}

//-----------------------
carray wave2;
set_array_size(wave2, 500, 1, 20.0, 0);
set_zero(wave2);
set_intial_cond(wave2);
solve_arrayRK2(wave2, 8.0, 0.5000);

carray wave3;//my main array
set_array_size(wave3, 500, 1, 20.0, 0);//array, xsize, ysize, dimension
set_zero(wave3);
set_intial_cond(wave3);
solve_arrayRK2(wave3, 8.0, 0.5015);//def stable




//----------------------Draw Data---------------------//

if(1)//start root application
{
	TApplication theApp("App", &argc, argv);//no more than two subs
    draw_graph_wave1_p2(wave1, wave2, wave3);
	//draw_graph_q1(wave1, wave2, wave3, analytic1, analytic2, analytic3);
    //draw_graph_q1a(wave1, wave2, wave3, analytic1, analytic2, analytic3);
	theApp.Run();
}



//end
}



