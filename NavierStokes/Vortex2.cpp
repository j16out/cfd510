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

int N = 5;

for(int i = 0; N < 60; ++i){
carray flow1;//my main array
carray flow2;



set_array_size(flow1, N, N, 1.0, 1.0, 0);//array, xsize, ysize, dimension
set_array_size(flow2, N*2, N*2, 1.0, 1.0, 0);


set_zero(flow1);
set_zero(flow2);



//---------------------solve IE----------------------//

//solve:( array, maxtime, cfl/timestep, lid velocity, data array)
solve_array_IE(flow1, 100.0, 0.1, 1.0, mydata); 

solve_array_IE(flow2, 100.0, 0.1, 1.0, mydata); 

get_l2norm_2array(flow1, flow2, mydata2, N);

N = N+10;
}


//----------------------Draw Data---------------------//
if(1)//start root application
{
	TApplication theApp("App", &argc, argv);//no more than two subs 
	draw_order_l2(mydata2);
	//draw_3Dgraph_s(flow1, flow2); 
	theApp.Run();
}
    //get_vortex(flow1);
	//draw_u(flow1);
	
	//draw_stab_l2(mydata);
	//draw_3Dgraph_s(flow1, flow2);

//end
}
