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

carray flux;//my main array
carray analytic;

//set array size or default used 162x162
set_array_size(flux, 20, 20, 1.0, 1.0, 0);//array, xsize, ysize, dimension
set_array_size(analytic, 20, 20, 1.0, 1.0, 0);

set_zero(flux);
set_zero(analytic);

//---------------------solve array1----------------------//

solve_array(flux, 1.0, 1.0);
set_analytic(analytic, flux);









int n = 10;
for(int i = 0; i < 8; ++i)
{
carray flux1;//my main array
carray analytic1;


//set array size or default used 162x162
set_array_size(flux1, n, n, 1.0, 1.0, 0);//array, xsize, ysize, dimension
set_array_size(analytic1, n, n, 1.0, 1.0, 0);


set_zero(flux1);
set_zero(analytic1);

//print_array(wave1);//print array in terminal


//---------------------solve array1----------------------//

solve_array(flux1, 1.0, 1.0);
set_analytic(analytic1, flux1);

//----------------------get error------------------------//

double lf1 = get_l2norm(flux1, analytic1);
mydata.l2norm.push_back(lf1);
n = n+10;
}


//----------------------Draw Data---------------------//
if(1)//start root application
{
	TApplication theApp("App", &argc, argv);//no more than two subs 
	draw_order_l2(mydata); 
	draw_3Dgraph(flux, analytic);
	theApp.Run();
}



//end
}



