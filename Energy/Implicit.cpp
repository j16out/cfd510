/*-------------------------------------------------------------------------------//
Main Program for finding solutions for wave equation. Employs a RK2 time advance 
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
cdata mydata;

carray flow1;//my main array
carray flow2;
carray analytic;

//set array size or default used 162x162
set_array_size(flow1, 25, 10, 5.0, 1.0, 0);//array, xsize, ysize, dimension
set_array_size(flow2, 25, 10, 5.0, 1.0, 0);
set_array_size(analytic, 20, 20, 5.0, 1.0, 0);

set_zero(flow1);
set_zero(flow2);
set_zero(analytic);
//print_array(flow2);
//---------------------solve EE----------------------//

solve_array_EE(flow1, 5.101, 0.01);
set_analytic(analytic, flow1);


//---------------------solve IE----------------------//
solve_array_IE(flow2, 5.101, 0.1);
//print_array(flow2);


//----------------------Draw Data---------------------//
if(1)//start root application
{
	TApplication theApp("App", &argc, argv);//no more than two subs  
	draw_3Dgraph(flow1, flow2);
	theApp.Run();
}



//end
}



