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
carray gsarraySOR1;//my main array
carray gsarraySOR13;
carray gsarraySOR15;




//set ghost cells as boundary conditions
set_zero(gsarraySOR1);
set_ghostcells(gsarraySOR1);

set_zero(gsarraySOR13);
set_ghostcells(gsarraySOR13);

set_zero(gsarraySOR15);
set_ghostcells(gsarraySOR15);



//---------------------GS SOR w=1 loop----------------------//
diff = 1;

while(diff > 0.00000001)
{
diff = gs_iter_SOR(gsarraySOR1, 1.0);
//cout << "difference " << diff << "\n"; 
if(diff > 100000000)
break;
	if(gsarraySOR1.iterations > 100000){
	break;
	cout << "solution failed to converge\n";
	}
}

cout << "Iterations: " << gsarraySOR1.iterations << "\n";



//---------------------GS SOR w=1.3 loop----------------------//
diff = 1;

while(diff > 0.00000001)
{
diff = gs_iter_SOR(gsarraySOR13, 1.3);
//cout << "difference " << diff << "\n"; 
if(diff > 100000000)
break;
	if(gsarraySOR13.iterations > 100000){
	break;
	cout << "solution failed to converge\n";
	}
}

cout << "Iterations: " << gsarraySOR13.iterations << "\n";



//---------------------GS SOR w=1.5 loop----------------------//
diff = 1;

while(diff > 0.00000001)
{
diff = gs_iter_SOR(gsarraySOR15, 1.5);
//cout << "difference " << diff << "\n"; 
if(diff > 100000000)
break;
	if(gsarraySOR15.iterations > 100000){
	break;
	cout << "solution failed to converge\n";
	}

}

cout << "Iterations: " << gsarraySOR15.iterations << "\n";




//----------------------Draw Data---------------------//

if(1)
{
	TApplication theApp("App", &argc, argv);
	draw_graph_l2norm3(gsarraySOR1, gsarraySOR13, gsarraySOR15);
	theApp.Run();
}



//end
}
