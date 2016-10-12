/*-------------------------------------------------------------------------------//
Main Program for looking at Point Guass-Seidel Scheme for Laplace Problem.
Finds maximum change in solution and L2 norm for w values of 1, 1.5 for10x10 array


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

#define BIG 1000000
#define E06 0.000001
#define E07 0.0000001
#define E08 0.00000001
#define E09 0.000000001


using namespace std;



int main(int argc, char **argv)
{

float diff = 1;
carray gsarray;//my main array
carray gsarraySOR;


//set size
set_array_size(gsarray, 10, 10, 1.0);
set_array_size(gsarraySOR, 10, 10, 1.0);


//set ghost cells as boundary conditions
set_zero(gsarray);
set_ghostcells(gsarray);

set_zero(gsarraySOR);
set_ghostcells(gsarraySOR);

print_array(gsarraySOR);



cout << "done!\n";


//---------------------GS loop----------------------//

solve_arraySOR(gsarray, E09, 1.0);

//---------------------GS SOR loop----------------------//

solve_arraySOR(gsarraySOR, E09, 1.5);

//----------------------Draw Data---------------------//

if(1)
{
	TApplication theApp("App", &argc, argv);
	draw_graph(gsarray, gsarraySOR);
	draw_3Dgraph(gsarray, gsarraySOR);
	theApp.Run();
}



//end
}
