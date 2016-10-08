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
#include "laplace/laplace.hpp"


//g++ rocketSIM.cpp calc.cpp -Wall -o2 -o test1 `root-config --cflags --glibs` -std=c++0x -pthread

using namespace std;



int main(int argc, char **argv)
{

float diff = 1;
carray poisson1;//my main array





//set ghost cells as boundary conditions


set_zero(poisson1);
set_ghostcells(poisson1);



//---------------------GS SOR w=1 loop----------------------//
diff = 1;

while(diff > 0.0000001)
{
diff = gs_iter_SOR(poisson1, 1.3);
//cout << "difference " << diff << "\n"; 
if(diff > 100000000)
break;
	if(poisson1.iterations > 100000){
	break;
	cout << "solution failed to converge\n";
	}
}

cout << "Iterations: " << poisson1.iterations << "\n";





//----------------------Draw Data---------------------//

if(1)
{
	TApplication theApp("App", &argc, argv);
	draw_3Dgraph(poisson1);
	theApp.Run();
}



//end
}
