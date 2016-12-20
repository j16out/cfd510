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

int N = 20;

carray flow1;//my main array
carray analytic;


set_array_size(flow1, N, N, 1.0, 1.0, 0);//array, xsize, ysize, dimension
set_array_size(analytic, N, N, 1.0, 1.0, 0);

set_zero(flow1);
set_zero(analytic);
set_init_cond(analytic);


//---------------------solve IE----------------------//

solve_array_IE(flow1, 10.0, 0.05, mydata);




//----------------------Draw Data---------------------//
if(1)//start root application
{
	TApplication theApp("App", &argc, argv);//no more than two subs  
	draw_stab_l2(mydata);
	draw_3Dgraph_s(flow1, analytic);
	//draw_order_l2(mydata);
	theApp.Run();
}



//end
}
