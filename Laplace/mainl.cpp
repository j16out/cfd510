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
carray gsarray;//my main array
carray gsarraySOR;





//set ghost cells as boundary conditions
set_zero(gsarray);
set_ghostcells(gsarray);

set_zero(gsarraySOR);
set_ghostcells(gsarraySOR);

print_mcell(gsarraySOR);



cout << "done!\n";


//---------------------GS loop----------------------//

while(diff > E07)
{
diff = gs_iter_SOR(gsarray, 1.0);
//cout << "difference " << diff << "\n"; 
if(diff > BIG)
break;

}

print_mcell(gsarray);
cout << "Iterations: " << gsarray.iterations << "\n";

//---------------------GS SOR loop----------------------//
float diff2 = 1;

while(diff2 > E07)
{
diff2 = gs_iter_SOR(gsarraySOR, 1.5);
//cout << "difference " << diff << "\n"; 
if(diff2 > BIG)
break;

}
print_mcell(gsarraySOR);
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
