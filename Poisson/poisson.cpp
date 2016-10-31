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


carray poisson1;//my main array
carray poisson2;
carray poisson3;

//set array size or default used 162x162
set_array_size(poisson1, 20, 20, 1.0);//array, xsize, ysize, dimension
set_array_size(poisson2, 40, 40, 1.0);
set_array_size(poisson3, 80, 80, 1.0);


//set ghost cells as boundary conditions


set_zero(poisson1);
set_ghostcells(poisson1);//set ghost cells/boundaries
//print_array(poisson1);//print array in terminal


set_zero(poisson2);
set_ghostcells(poisson2);
//print_array(poisson2);

set_zero(poisson3);
set_ghostcells(poisson3);
//print_array(poisson3);



//---------------------GS SOR w=1.3 loop 1----------------------//

solve_arraySOR(poisson1, E11, 1.3);
cout << "Solution: " << get_solution(poisson1) << "\n";


//---------------------GS SOR w=1.3 loop 2----------------------//

solve_arraySOR(poisson2, E11, 1.3);
cout << "Solution: " << get_solution(poisson2) << "\n";

//---------------------GS SOR w=1 loop 3----------------------//

solve_arraySOR(poisson3, E11, 1.3);
cout << "Solution: " << get_solution(poisson3) << "\n";

//---------------------calc error based on ASME---------------//

get_discrete_Error(poisson3, poisson2, poisson1, 1.0);


//----------------------Draw Data---------------------//



if(1)//start root application
{
	TApplication theApp("App", &argc, argv);
	draw_3DgraphP(poisson3);//draw 3d graph
	theApp.Run();
}



//end
}



