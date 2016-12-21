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
cdata mydata2;

int N = 10;

carray flow1;//my main array
carray flow2;
carray flow3;//my main array
carray flow4;


set_array_size(flow1, 40, 120, 1.0, 3.0, 0);//array, xsize, ysize, dimension
//set_array_size(flow2, 3*N, 3*N, 1.0, 1.0, 0);
//set_array_size(flow3, 9*N, 9*N, 1.0, 1.0, 0);//array, xsize, ysize, dimension
//set_array_size(flow4, N, N, 1.0, 1.0, 0);

set_zero(flow1);
//set_zero(flow2);
//set_zero(flow3);
//set_zero(flow4);



//---------------------solve IE----------------------//

//solve:( array, maxtime, cfl/timestep, lid velocity, data array)
solve_array_IE(flow1, 100.0, 0.1, 1.0, mydata); 

//solve_array_IE(flow2, 100.0, 0.1, 1.0, mydata); 

//solve_array_IE(flow3, 100.0, 0.1, 1.0, mydata); 

//solve_array_IE(flow4, 100.0, 0.1, 1.0, mydata); 




//----------------------Draw Data---------------------//
if(1)//start root application
{
	TApplication theApp("App", &argc, argv);//no more than two subs  
	//draw_stab_l2(mydata);
	draw_3Dgraph_s(flow1, flow2);
	draw_u(flow1);
	//draw_order_l2(mydata);
	theApp.Run();
}



//end
}
