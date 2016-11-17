/*-------------------------------------------------------------------------------//
Main Program for finding solution for wave equation. Employs a RK2 time advance 
with 2nd order upwind flux scheme.


Jerin Roberts 2016
compiled using g++/gcc version 5.4.0 on Ubuntu 16.04.02 and are available for clone 
via the link provided: url{https://github.com/j16out/
//-------------------------------------------------------------------------------*/

//***NOTE need to change intial condition in numerical.cpp***//



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




int main(int argc, char **argv)
{


carray wave1;//my main array
set_array_size(wave1, 500, 1, 20.0, 0);//array, xsize, ysize, dimension
set_zero(wave1);
set_intial_cond(wave1);
solve_arrayRK2(wave1, 8.0, 0.1);//def stable solution

//---------------------plot l2 vs cfl----------------------//
float cfl = 0.1;
float factor = 0.0025;

while(cfl <= 1.0)
{
carray wave2;//solve array for different cfl
set_array_size(wave2, 500, 1, 20.0, 0);
set_zero(wave2);
set_intial_cond(wave2);
solve_arrayRK2(wave2, 8.0, cfl);
float l2 = get_l2norm(wave1, wave2);//compare to stable solution

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

//-----------------------look at particular cfl------------------//
carray wave2;
set_array_size(wave2, 500, 1, 20.0, 0);
set_zero(wave2);
set_intial_cond(wave2);
solve_arrayRK2(wave2, 8.0, 0.5000);

carray wave3;//my main array
set_array_size(wave3, 500, 1, 20.0, 0);//array, xsize, ysize, dimension
set_zero(wave3);
set_intial_cond(wave3);
solve_arrayRK2(wave3, 8.0, 0.5015);




//----------------------Draw Data---------------------//

if(1)//start root application
{
	TApplication theApp("App", &argc, argv);//no more than two subs
    draw_graph_wave1_p2(wave1, wave2, wave3);
	draw_graph_q1(wave1, wave2, wave3, analytic1, analytic2, analytic3);
	theApp.Run();
}



//end
}



