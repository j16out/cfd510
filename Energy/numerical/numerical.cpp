#include "numerical.hpp"





//**************************************************************************//
//---------------------------Setting Array----------------------------------//
//**************************************************************************//


//----------set array size (working area excluding ghost)---------------//

void set_array_size(carray & myarray, int x, int y, double DIMx, double DIMy, int scheme)
{
	if(x <= 100 && y <= 100)
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
double DIMx = myarray.DIMx;

//set boundary conditions in ghost cells
if(myarray.scheme == 0)//2nd order upwind
{
myarray.T1[0][1] = -2.0*(sin(4.0*PI*myarray.ctime)) + 3.0*myarray.T1[1][1];
myarray.T1[1][1] = 2.0*(sin(4.0*PI*myarray.ctime)) - myarray.T1[2][1];
}

if(myarray.scheme == 1)//1st order upwind
{
myarray.T1[1][1] = 2.0*(sin(4.0*PI*myarray.ctime)) - myarray.T1[2][1];
}

if(myarray.scheme == 2)//2nd order centered
{
myarray.T1[1][1] = 2.0*(sin(4.0*PI*myarray.ctime)) - myarray.T1[2][1];
}


	
}

//--------------------------set intial condition---------------------------------//

void set_intial_cond(carray & myarray)
{
double DIMx = myarray.DIMx;
double DIMy = myarray.DIMy;
double dx = 0.0;
double dy = 0.0;

for(int j = 0; j < myarray.sizey-0; ++j)
{

    for(int i = 0; i < myarray.sizex-0; ++i)
    {
    dx = (i-0.5)*DIMx;
    dy = (i-0.5)*DIMy;
    

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
double a = (s1.uip1_j * s1.Tip1_j  -  s1.uim1_j * s1.Tim1_j)/(2.0*chx);
double b = (s1.Tip1_j  - 2.0*s1.Ti_j + s1.Tim1_j)/(chx*chx*(RE*PR));
double c = (s1.vi_jp1 * s1.Ti_jp1  -  s1.vi_jm1 * s1.Ti_jm1)/(2.0*chy);
double d = (s1.Ti_jp1  - 2.0*s1.Ti_j + s1.Ti_jm1)/(chy*chy*(RE*PR));

double newcell = (-a+b)   +   (-c+d)   +   source;
return newcell*-1.0;
}



//**************************************************************************//
//---------------------------Error Checking---------------------------------//
//**************************************************************************//





//-------------------------Get L1 norm for unknown analytical----------------------//

double get_l1norm(carray & myarray, carray myarray2)
{
double l1sum =0;
double sx = myarray.sizex-2;
double sy = myarray.sizey-2;

for(int j = 1; j < myarray.sizey-1; ++j)
{	
	for(int i = 2; i < myarray.sizex; ++i)
	{

	double P = myarray.T1[i][j];
	double T = myarray2.T1[i][j];
	l1sum =  l1sum + abs(P-T);

	}

}

double l1 = l1sum/(sx);
cout << setprecision(8) << fixed << "L1 norm: " << l1 << "\n";
return l1;
}


//-------------------------Get L infinty norm for unknown analytical----------------------//

double get_linf_norm(carray & myarray, carray myarray2)
{
double error =0;
double maxerror = -1;
double sx = myarray.sizex-2;
double sy = myarray.sizey-2;

for(int j = 1; j < myarray.sizey-1; ++j)
{	
	for(int i = 2; i < myarray.sizex; ++i)
	{

	double P = myarray.T1[i][j];
	double T = myarray2.T1[i][j];
	error =  abs(P-T);
    if(error > maxerror)
    maxerror = error;

	}

}


cout << setprecision(8) << fixed << "L infinity norm: " << maxerror << "\n";
return maxerror;
}



//-------------------------Get L2 nrom for unknown analytical----------------------//

double get_l2norm(carray & myarray, carray myarray2)
{
double l2sum =0;
double sx = myarray.sizex-2;
double sy = myarray.sizey-2;

for(int j = 1; j < myarray.sizey-1; ++j)
{	
	for(int i = 2; i < myarray.sizex; ++i)
	{

	double P = myarray.T1[i][j];
	double T = myarray2.T1[i][j];
	l2sum =  l2sum + pow((P-T),2);

	}

}

double l2 = sqrt(l2sum/(sx));
cout << setprecision(8) << fixed << "L2 norm: " << l2 << "\n";
return l2;
}



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






