#include "numerical.hpp"





//**************************************************************************//
//---------------------------Setting Array----------------------------------//
//**************************************************************************//


//----------set array size (working area excluding ghost)---------------//

void set_array_size(carray & myarray, int x, int y, double DIMx, double DIMy, int scheme)
{
	if(x <= maxx && y <= maxy)
	{
	myarray.sizex = x+2;
	myarray.sizey = y+2;
	myarray.DIMx = DIMx/(x);
	myarray.DIMy = DIMy/(y);
	myarray.scheme = scheme;
	}
	else
	cout << "Array size to big, setting to default 100" << "\n";

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
		if(myarray.T1[i][j] >= 0)
		cout << setprecision(3) << fixed << myarray.T1[i][j] <<"|";
		if(myarray.T1[i][j] < 0)
		cout << setprecision(2) << fixed << myarray.T1[i][j] <<"|";
		}
	
	}
cout << "\n";
}


//--------------------------zero array----------------------------//

void set_zero(carray & myarray)
{
	for(int j = 0; j < myarray.sizey; ++j)
	{
		for(int i = 0; i < myarray.sizex; ++i)
		{
		myarray.T1[i][j] = 0;//set everything to zero
        myarray.f1[i][j] = 0;
        myarray.v1[i][j] = 0;
        myarray.u1[i][j] = 0;

		}
	}
}


//--------------------------set ghost cells for Wave----------------------------//


void set_ghostcells(carray & myarray)
{
//set boundary conditions in ghost cells


}

//--------------------------set intial condition---------------------------------//

void set_intial_cond(carray & myarray)
{
double DIMx = myarray.DIMx;
double DIMy = myarray.DIMy;
double dx = 0.0;
double dy = 0.0;

for(int j = 0; j < myarray.sizey; ++j)
{

    for(int i = 0; i < myarray.sizex; ++i)
    {
    dx = (i-0.5)*DIMx;
    dy = (j-0.5)*DIMy;
    

    myarray.T1[i][j] = T0*cos(PI*dx)*sin(PI*dy);
    myarray.u1[i][j] = U0*dy*sin(PI*dx);
    myarray.v1[i][j] = V0*dx*cos(PI*dy);
    //printf("f: %f  dx: %f\n", f, dx);
    }

}
}


//**************************************************************************//
//----------------------------- Array Solving-------------------------------//
//**************************************************************************//

void solve_array(carray & myarray, double tmax, double cfl)
{


double ctime = myarray.ctime;

set_intial_cond(myarray);



compute_Flux(myarray);



printf("Solved numeric at %f time\n",ctime);
}


//**************************************************************************//
//---------------------------Compute Flux-----------------------------------//
//**************************************************************************//


void compute_Flux(carray & myarray)
{
surr mysurr;

for(int j = 1; j < myarray.sizey-1; ++j)
{
    for(int i = 1; i < myarray.sizex-1; ++i)
    {
    //----get surrounding cells and compute new cell----//
    get_nsurcells(myarray, i, j, mysurr);
    //-----update current cell----//
    double newcell = calc_newcell(myarray, mysurr);    
    myarray.f1[i][j] = newcell;   
    }

}
}

//--------------------------get surrounding cells-----------------------//

void get_nsurcells(carray & myarray, int i, int j, surr & mysurr)
{

mysurr.Tim1_j = myarray.T1[i-1][j];
mysurr.Tip1_j = myarray.T1[i+1][j];
mysurr.Ti_j = myarray.T1[i][j]; 
mysurr.Ti_jm1 = myarray.T1[i][j-1];
mysurr.Ti_jp1 = myarray.T1[i][j+1];

mysurr.uim1_j = myarray.u1[i-1][j];
mysurr.uip1_j = myarray.u1[i+1][j];
 

mysurr.vi_jm1 = myarray.v1[i][j-1];
mysurr.vi_jp1 = myarray.v1[i][j+1];
}


//-----------------------calculate flux for new cell--------------------//

double calc_newcell(carray & myarray, surr & s1)
{
double chx = myarray.DIMx;
double chy = myarray.DIMy;
//float chy = DIM1;
double source = 0.0;
double a = (s1.uip1_j * s1.Tip1_j  -  s1.uim1_j * s1.Tim1_j)/(2.0);
double b = (s1.Tip1_j  - 2.0*s1.Ti_j + s1.Tim1_j)/(chx*(RE*PR));
double c = (s1.vi_jp1 * s1.Ti_jp1  -  s1.vi_jm1 * s1.Ti_jm1)/(2.0);
double d = (s1.Ti_jp1  - 2.0*s1.Ti_j + s1.Ti_jm1)/(chy*(RE*PR));

double newcell = (-1.0/chx)*(a-b)   +   (-1.0/chy)*(c-d)   +   source;
return -newcell;
}



//**************************************************************************//
//---------------------------Error Checking---------------------------------//
//**************************************************************************//





//----------------------------Set a Analytical Solution------------------------------//
void set_analytic(carray & myarray, carray & numarray)
{
double DIMx = myarray.DIMx;
double DIMy = myarray.DIMy;
double ctime = numarray.ctime;
for(int j = 1; j < myarray.sizey-1; ++j)
{

	for(int i = 1; i < myarray.sizex-1; ++i)
	{
	double dx = (i-0.5)*DIMx;
	double dy = (j-0.5)*DIMy;
	
	double a = U0*T0*PI*cos(2.0*PI*dx)*dy*sin(PI*dy);
	double b = V0*T0*PI*dx*cos(PI*dx)*cos(2.0*PI*dy);
	double c = (2.0*T0*pow(PI,2)*cos(PI*dx)*sin(PI*dy))/(RE*PR);
	myarray.f1[i][j] = a + b + c;
	
	}

}

printf("setting analytic at %f time\n",ctime);
}






