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

#define BIG 1000000
#define E07 0.0000001
#define E08 0.00000001
#define E09 0.000000001


int main(int argc, char **argv)
{

float diff = 1;
carray poisson1;//my main array





//set ghost cells as boundary conditions


set_zero(poisson1);
set_ghostcells(poisson1);



//---------------------GS SOR w=1 loop----------------------//
diff = 1;
int update = 0;
int update2 = 100;

while(diff > E07)
{
diff = gs_iter_SOR(poisson1, 1.4);



if(diff > BIG)
break;
	
if(poisson1.iterations > 100000){
break;
cout << "solution failed to converge\n";
}

if(update >= update2)
{cout << "Update: step " << update << " divergence " << diff << " \n"; 
 update2 = update2 + 100;
}	

++update;	
}

cout << "Iterations: " << poisson1.iterations << "\n";





//----------------------Draw Data---------------------//

if(1)
{
	TApplication theApp("App", &argc, argv);
	draw_3DgraphP(poisson1);
	theApp.Run();
}



//end
}
