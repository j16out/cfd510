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
carray gsarray;//my main array
carray gsarraySOR;





//set ghost cells as boundary conditions
set_zero(gsarray);
set_ghostcells(gsarray);

set_zero(gsarraySOR);
set_ghostcells(gsarraySOR);
cout << "done!\n";


//---------------------GS loop----------------------//

while(diff > 0.0000001)
{
diff = gs_iter(gsarray);
//cout << "difference " << diff << "\n"; 
if(diff > 100000000)
break;

}

cout << "Iterations: " << gsarray.iterations << "\n";

//---------------------GS SOR loop----------------------//
float diff2 = 1;

while(diff2 > 0.0000001)
{
diff2 = gs_iter_SOR(gsarraySOR, 1.5);
//cout << "difference " << diff << "\n"; 
if(diff2 > 100000000)
break;

}

cout << "Iterations: " << gsarraySOR.iterations << "\n";



//----------------------Draw Data---------------------//

if(1)
{
	TApplication theApp("App", &argc, argv);
	draw_graph(gsarray, gsarraySOR);
	theApp.Run();
}



//end
}
