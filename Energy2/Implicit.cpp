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


//--------Gets l2 and order for either IE or EE---------//

int main(int argc, char **argv)
{
cdata mydata;
double l2;
int n = 25;
int p = 10;
int nt = 200;


//loop for different array sizes
while(n < nt)
{
carray flow1;//my main array
carray flow2;

carray analytic;

//set array size or default used 162x162
set_array_size(flow1, n, p, 5.0, 1.0, 0);//array, xsize, ysize, dimension
set_array_size(flow2, n, p, 5.0, 1.0, 0);
set_array_size(analytic, n, p, 5.0, 1.0, 0);

//set zero
set_zero(flow1);
set_zero(flow2);
set_zero(analytic);


//---------------------solve IE----------------------//
solve_array_IE(flow1, 4.5, 0.1);
solve_array_IE(flow2, 4.5, 0.1);
set_analytic(analytic, flow1);

//get l2 norm and save
l2 = get_l2norm(flow1, analytic);
mydata.l2norm.push_back(l2);

n +=5;
p +=2;
}


//----------------------Draw Data---------------------//
if(1)//start root application
{
	TApplication theApp("App", &argc, argv);//no more than two subs  
	//draw_3Dgraph(flow1, analytic);
	draw_order_l2(mydata); 
	theApp.Run();
}



//end
}




