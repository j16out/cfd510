#include "numerical.hpp"





//**************************************************************************//
//---------------------------Setting Array----------------------------------//
//**************************************************************************//


//----------set array size (working area excluding ghost)---------------//

void set_array_size(carray & myarray, int x, int y, float DIM)
{
	if(x < 160 && y < 160)
	{
	myarray.sizex = x+2;
	myarray.sizey = y+2;
	myarray.DIM1 = DIM/(x);
	}
	else
	cout << "Array size to big, setting to default 160" << "\n";

}

//--------------------------Print array in terminal----------------------------//

void print_array(carray & myarray)
{
cout << "\n";

	for(int j = 0; j < myarray.sizey; ++j)
	{
	cout << "\n|";	
		for(int i = 0; i < myarray.sizex; ++i)
		{
		if(myarray.mcellSOL[i][j] >= 0)
		cout << setprecision(3) << fixed << myarray.mcellSOL[i][j] <<"|";
		if(myarray.mcellSOL[i][j] < 0)
		cout << setprecision(2) << fixed << myarray.mcellSOL[i][j] <<"|";
		}
	
	}
cout << "\n";
}


//--------------------------zero array----------------------------//

void set_zero(carray & myarray)
{
	for(int i = 0; i < myarray.sizex; ++i)
	{
		for(int j = 0; j < myarray.sizey; ++j)
		{
		myarray.mcellSOL[i][j] = 0;//set everything to zero

		}
	}
}


//--------------------------set ghost cells for Wave----------------------------//


void set_ghostcells(carray & myarray)
{
float DIM1 = myarray.DIM1;

//set boundary conditions in ghost cells
for(int i = 1; i < myarray.sizex-1; ++i)
	{float d = (i-0.5)*DIM1;
	
	//neumann boundaries	
	myarray.mcellSOL[i][0] = myarray.mcellSOL[i][1];//top ghost cells 
	myarray.mcellSOL[0][i] = myarray.mcellSOL[1][i];//left ghost cells
		
	//dirichlet boundaries	
	myarray.mcellSOL[i][myarray.sizex-1] = 2.0*(5.0-((1.0/2.0)*pow((1.0+pow(d, 2)),3))) - myarray.mcellSOL[i][myarray.sizex-2];
	myarray.mcellSOL[myarray.sizey-1][i] = 2.0*(5.0-((1.0/2.0)*pow((1.0+pow(d, 2)),3))) - myarray.mcellSOL[myarray.sizex-2][i];
	
        }
}

void set_intial_cond(carray & myarray)
{
float DIM1 = myarray.DIM1;
float dx =0.0;
float f;

for(int j = 1; j < myarray.sizey-1; ++j)
{

    for(int i = 1; i < myarray.sizex-1; ++i)
    {
    dx = (i-0.5)*DIM1;
    f = -sin(2.0*PI*dx);
    myarray.mcellSOL[i][j] = f;
    printf("f: %f  dx: %f\n", f, dx);
    }

}
}


//**************************************************************************//
//---------------------------Array Solving----------------------------------//
//**************************************************************************//


//--------------------------Set FI values for array mcellFI Interior----------------------------//

int get_FIarray(carray & myarray, int stage)
{
float DIM1 = myarray.DIM1;//get dimensions of array (not grid size)
float l2sum = 0.0;
float sx = myarray.sizex-2;
float sy = myarray.sizey-2;

float dx = 0.0;//change in x
float dy = 0.0;//change in y

carray oldarray=myarray;

for(int j = 1; j < myarray.sizey-1; ++j)
{

    for(int i = 2; i < myarray.sizex-2; ++i)
    {
    dx = (i-0.5)*DIM1;
    dy = (j-0.5)*DIM1;

    //----get surrounding cells and compute new cell-------//
    get_surcells(myarray, i, j, stage);
    float newcell = calc_FI(myarray); 

    //-----update current cell----//
    if(stage == 1)
    myarray.mcellFI[i][j] = newcell;

    if(stage == 2)
    myarray.mcellFI2[i][j] = newcell;
    }

}
	

return 1;
}


//--------------------Calculate new cell value from neighbors ---------------------//


float calc_FI(carray & myarray)
{
float DIM1 = myarray.DIM1;
float chx = DIM1;
//float chy = DIM1;
float temp = 1.0;
float newcell = (2.0)*(3.0*myarray.Ti_j+4.0*myarray.Tim1_j+myarray.Tim2_j)/(2.0*chx);

return newcell;
}




//--------------------------Get current cell values----------------------------//

void get_surcells(carray & myarray, int i, int j, int stage)
{
float fcon = false;
float sizex = myarray.sizex;
float sizey = myarray.sizey;

if(stage == 1)
{
myarray.Tim1_j = myarray.mcellSOL[i-1][j];
myarray.Tim2_j = myarray.mcellSOL[i-2][j];
myarray.Ti_jm1 = myarray.mcellSOL[i][j-1];
myarray.Ti_jm2 = myarray.mcellSOL[i][j-2];
myarray.Ti_j = myarray.mcellSOL[i][j]; 
} 

if(stage == 2)
{
myarray.Tim1_j = myarray.mcellSOL2[i-1][j];
myarray.Tim2_j = myarray.mcellSOL2[i-2][j];
myarray.Ti_jm1 = myarray.mcellSOL2[i][j-1];
myarray.Ti_jm2 = myarray.mcellSOL2[i][j-2];
myarray.Ti_j = myarray.mcellSOL2[i][j]; 
} 
        

}




//--------------------------Solve array using RK2 and 2ndUW----------------------------//

void solve_arrayRK2(carray & myarray)
{
int tomp;
float tmax = myarray.tmax;
float ctime = myarray.ctime;
while(ctime < tmax)
{

tomp = get_FIarray(myarray, 1);
tomp = get_RK2(myarray, 1);

tomp = get_FIarray(myarray, 2);
tomp = get_RK2(myarray, 2);


ctime = myarray.ctime;
}

}

//--------------------------Solve RK2 interation----------------------------//


int get_RK2(carray & myarray, int stage)
{

if(stage == 1)
{
for(int j = 1; j < myarray.sizey-1; ++j)
{

    for(int i = 2; i < myarray.sizex-2; ++i)
    {
      myarray.mcellSOL2[i][j] = myarray.mcellSOL[i][j]-myarray.tstep*(myarray.mcellFI[i][j]);
    }

}
}

if(stage == 2)
{
for(int j = 1; j < myarray.sizey-1; ++j)
{

    for(int i = 2; i < myarray.sizex-2; ++i)
    {
     myarray.mcellSOL2[i][j] = myarray.mcellSOL[i][j]-myarray.tstep*((myarray.mcellFI2[i][j]+myarray.mcellFI[i][j])/2.0);
    }

}
}


}



//**************************************************************************//
//---------------------------Error Checking---------------------------------//
//**************************************************************************//



//--------------------------Get descrete error----------------------------//

void get_discrete_Error(carray ray1, carray ray2, carray ray3, float DIM)
{
//Calculating error as described in paper "procedure for estimation and reporting of uncertainty due to discretization in CFD applications"//

printf("\nCalculating Error...\n");

float h1 = DIM/ray1.sizex;
float h2 = DIM/ray2.sizex;
float h3 = DIM/ray3.sizex;


float sol1 = get_solution(ray1);
float sol2 = get_solution(ray2);
float sol3 = get_solution(ray3);



printf("h1: %f \nh2: %f \nh3: %f, \nsol1: %f \nsol2: %f \nsol3: %f\n",h1, h2, h3, sol1, sol2, sol3);

float r21 = h2/h1;
float r32 = h3/h2;

printf("\nr32: %f \nr21: %f\n",r32, r21);

float e32 = sol3-sol2;
float e21 = sol2-sol1;

float s = (e32/e21);
if(s >= 0)
s = 1;
else
s = -1;

float p_n = 0;
float p = (1/log(r21))*(abs(log(abs(e32/e21))+0));

printf("intial guess: %f \n", p);

float diff = 1;

	while(diff > 0.0000001)
	{

	float p_n = (1/log(r21))*(abs(log(abs(e32/e21))+log((pow(r21,p)-s)/(pow(r32,p)-s)) ));
	diff = abs(p_n -p);
	//printf("p_n: %f p: %f diff: %f\n",p_n, p, diff);

	p = p_n;
	}
 
//
float sol_ext21 = (pow(r21, p)*sol1-sol2)/(pow(r21,p)-1.0);
float sol_ext32 = (pow(r32, p)*sol2-sol3)/(pow(r32,p)-1.0);

printf("order: %f \nphi_ext21: %f \nphi_ext32 %f\n",p, sol_ext21, sol_ext32);

float ea21 = abs((sol1-sol2)/sol1);

float e_ext21 = abs((sol_ext21-sol1)/sol_ext21);

float GCI_21 = (1.25*ea21)/(pow(r21,p)-1.0);


printf("ea21: %f  \ne_ext21: %f  \nGC121 %f \n", ea21, e_ext21, GCI_21);

}



//-------------------------Get L2 nrom for unknown analytical----------------------//

float get_l2norm(carray & myarray, carray myarray2)
{
float l2sum =0;
float sx = myarray.sizex-2;
float sy = myarray.sizey-2;

for(int j = 1; j < myarray.sizey-1; ++j)
{	
	for(int i = 1; i < myarray.sizex-1; ++i)
	{

	float P = myarray.mcellSOL[i][j];
	float T = myarray2.mcellSOL[i][j];
	l2sum =  l2sum + pow((P-T),2);

	}

}

float l2 = sqrt(l2sum/(sx*sy));
cout << "L2 norm: " << l2 << "\n";
return l2;
}



//----------------------------Set a Analytical Solution------------------------------//
void set_analytic(carray & myarray)
{
float DIM1 = myarray.DIM1;
for(int j = 1; j < myarray.sizey-1; ++j)
{

	for(int i = 1; i < myarray.sizex-1; ++i)
	{
	float dx = (i-0.5)*DIM1;
	float dy = (j-0.5)*DIM1;
	float T = sin(2*PI*(2*(1)-dx));
	myarray.mcellSOL[i][j] = T;
	}

}
}


//-----------------------Get average solution at point (1/2)(1/2)--------------------//

float get_solution(carray & myarray)
{
float DIM1 = myarray.DIM1;
int sx = (myarray.sizex)/2.0;
int sy = (myarray.sizey)/2.0;
float sol = (myarray.mcellSOL[sx-1][sy]+myarray.mcellSOL[sx][sy]+myarray.mcellSOL[sx][sy-1]+myarray.mcellSOL[sx-1][sy-1])/4.0;

printf("cell 1: %f cell 2: %f cell 3: %f cell 4: %f\n",myarray.mcellSOL[sx-1][sy],myarray.mcellSOL[sx][sy],myarray.mcellSOL[sx][sy-1],myarray.mcellSOL[sx-1][sy-1]);
//for Poisson problem only, finds value based on average of four surrounding cells

return sol;
}





