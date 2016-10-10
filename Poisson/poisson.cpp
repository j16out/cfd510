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
carray poisson2;
carray poisson3;

//set array size or default used 162x162
set_array_size(poisson1, 20, 20, 1.0);
set_array_size(poisson2, 40, 40, 1.0);
set_array_size(poisson3, 60, 60, 1.0);


//set ghost cells as boundary conditions


set_zero(poisson1);
set_ghostcells(poisson1);
print_array(poisson1);

set_zero(poisson2);
set_ghostcells(poisson2);
print_array(poisson2);

set_zero(poisson3);
set_ghostcells(poisson3);
print_array(poisson3);


//---------------------GS SOR w=1 loop 1----------------------//
diff = 1;
int update = 0;
int update2 = 100;

while(diff > E07)
{
diff = gs_iter_SOR(poisson1, 1.3);



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


//---------------------GS SOR w=1 loop 2----------------------//
diff = 1;
update = 0;
update2 = 100;

while(diff > E07)
{
diff = gs_iter_SOR(poisson2, 1.3);



if(diff > BIG)
break;
	
if(poisson2.iterations > 100000){
break;
cout << "solution failed to converge\n";
}

if(update >= update2)
{cout << "Update: step " << update << " divergence " << diff << " \n"; 
 update2 = update2 + 100;
}	

++update;	
}

cout << "Iterations: " << poisson2.iterations << "\n";

//---------------------GS SOR w=1 loop 3----------------------//
/*
diff = 1;
update = 0;
update2 = 100;

while(diff > E07)
{
diff = gs_iter_SOR(poisson3, 1.3);



if(diff > BIG)
break;
	
if(poisson3.iterations > 100000){
break;
cout << "solution failed to converge\n";
}

if(update >= update2)
{cout << "Update: step " << update << " divergence " << diff << " \n"; 
 update2 = update2 + 100;
}	

++update;	
}

cout << "Iterations: " << poisson3.iterations << "\n";*/


//----------------------Draw Data---------------------//

if(1)
{
	TApplication theApp("App", &argc, argv);
	draw_3DgraphP(poisson2);

	theApp.Run();
}



//end
}
