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
#define E07 0.0000001
#define E08 0.00000001
#define E09 0.000000001


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

while(diff > E08)
{
diff = gs_iter_SOR(gsarraySOR1, 1.0);
//cout << "difference " << diff << "\n"; 
if(diff > BIG)
break;
	if(gsarraySOR1.iterations > 10000){
	cout << "solution failed to converge\n";
	break;
	}
}

cout << "Iterations: " << gsarraySOR1.iterations << "\n";



//---------------------GS SOR w=1.3 loop----------------------//
diff = 1;

while(diff > E08)
{
diff = gs_iter_SOR(gsarraySOR13, 1.3);
//cout << "difference " << diff << "\n"; 
if(diff > BIG)
break;
	if(gsarraySOR13.iterations > 10000){
	cout << "solution failed to converge\n";
	break;
	}
}

cout << "Iterations: " << gsarraySOR13.iterations << "\n";



//---------------------GS SOR w=1.5 loop----------------------//
diff = 1;

while(diff > E08)
{
diff = gs_iter_SOR(gsarraySOR15, 1.5);
//cout << "difference " << diff << "\n"; 
if(diff > BIG)
break;
	if(gsarraySOR15.iterations > 10000){
	cout << "solution failed to converge\n";
	break;
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
