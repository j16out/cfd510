/*-------------------------------------------------------------------------------//
Main Program for finding solutions for energy equation. Employs a implicit euler
time advance with 2nd order centered flux scheme.


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
cdata mydata;

carray flow1;//my main array
carray flow2;
carray flow3;


//set array size or default used 162x162
set_array_size(flow1, 200, 80, 40.0, 1.0, 0);//array, xsize, ysize, dimension
set_array_size(flow2, 100, 40, 40.0, 1.0, 0);
set_array_size(flow3, 50, 20, 40.0, 1.0, 0);


set_zero(flow1);
set_zero(flow2);
set_zero(flow3);


//---------------------solve IE----------------------//
solve_array_IE(flow1, 14.5, 0.1);
solve_array_IE(flow2, 14.5, 0.1);
solve_array_IE(flow3, 14.5, 0.1);

//print_array(flow1);
double dx = 0.0;

find_max(flow1, dx);
get_discrete_Error(flow1, flow2, flow3);
//----------------------Draw Data---------------------//
if(1)//start root application
{
	TApplication theApp("App", &argc, argv);//no more than two subs  
	//draw_3Dgraph(flow1, flow2);
	find_maxvalues(flow1, flow2, flow3);
	theApp.Run();
}



//end
}




