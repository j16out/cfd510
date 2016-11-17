/*-------------------------------------------------------------------------------//
Main Program for finding solution for wave equation. Employs a RK2 time advance 
with 2nd order upwind flux scheme.


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




int main(int argc, char **argv)
{

//used as storage for this example
carray wave11;//my main array1
carray wave21;//my main array2
carray wave31;//my main array3




int smesh = 20;

while(smesh <= 2000)
{
//arrays being solved
carray wave1;
carray analytic1;
carray wave2;
carray analytic2;
carray wave3;
carray analytic3;


//set array size 
set_array_size(wave1, smesh, 1, 1.0, 0);//array, xsize, ysize, dimension, scheme
set_array_size(wave2, smesh, 1, 1.0, 1);
set_array_size(wave3, smesh, 1, 1.0, 2);

//set analytic size
set_array_size(analytic1, smesh, 1, 1.0, 0);
set_array_size(analytic2, smesh, 1, 1.0, 0);
set_array_size(analytic3, smesh, 1, 1.0, 0);

//set intial conditions
set_zero(wave1);
set_intial_cond(wave1);

set_zero(wave2);
set_intial_cond(wave2);

set_zero(wave3);
set_intial_cond(wave3);


float l2 = 0;
float l1 = 0;
float lin = 0;


//---------------------solve array1----------------------//
solve_arrayRK2(wave1, 1.0, 0.4);//array,time,cfl
set_analytic(analytic1, wave1);
l2 = get_l2norm(wave1, analytic1);
l1 = get_l1norm(wave1, analytic1);
lin = get_linf_norm(wave1, analytic1);
wave11.l1norm.push_back(l1);
wave11.linfnorm.push_back(lin);
//cout << "Solution: " << get_solution(poisson1) << "\n";



//---------------------solve array2----------------------//
solve_arrayRK2(wave2, 1.0, 0.4);
set_analytic(analytic2, wave2);
l2 = get_l2norm(wave2, analytic2);
l1 = get_l1norm(wave2, analytic2);
lin = get_linf_norm(wave2, analytic2);
wave21.l1norm.push_back(l1);
wave21.linfnorm.push_back(lin);


//---------------------solve array2----------------------//
solve_arrayRK2(wave3, 1.0, 0.4);
set_analytic(analytic3, wave3);
l2 = get_l2norm(wave3, analytic3);
l1 = get_l1norm(wave3, analytic3);
lin = get_linf_norm(wave3, analytic3);
wave31.l1norm.push_back(l1);
wave31.linfnorm.push_back(lin);

smesh = smesh + 50;

}

//----------------------Draw Data---------------------//

if(1)//start root application
{
	TApplication theApp("App", &argc, argv);//no more than two subs
	draw_graph_wave1_p3(wave11, wave21, wave31); 
	//draw_graph_q1(wave1, wave2, wave3, analytic1, analytic2, analytic3);   
	theApp.Run();
}

 //draw_graph_wave1(wave1, wave2, wave3);


//end
}



