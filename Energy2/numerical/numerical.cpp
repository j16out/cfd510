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


//--------------------------set ghost cells for Energy----------------------------//


void set_ghostcells(carray & myarray)
{
double DIMx = myarray.DIMx;
double DIMy = myarray.DIMy;
double dx = 0.0;
double dy = 0.0;
    //set ghost cells top/bottom
	for(int i = 0; i < myarray.sizex; ++i)
	{
	myarray.T1[i][0] = 3.0*(0.0) - (5.0/2.0)*myarray.T1[i][1] + (1.0/2.0)*myarray.T1[i][2];
	myarray.T1[i][myarray.sizey-1] =3.0*(1.0)- (5.0/2.0)*myarray.T1[i][myarray.sizey-2] + (1.0/2.0)*myarray.T1[i][myarray.sizey-3];
	}	
	
    //set ghost cells inflow/outflow	
	for(int j = 0; j < myarray.sizey; ++j)
	{ 
    dy = (j-0.5)*DIMy;
    
	myarray.T1[0][j] = 2.0*( dy+((3.0/4.0)*PR*EC*(pow(U0,2))*(1.0-pow((1.0-2.0*dy),4))) ) - myarray.T1[1][j];
	
	
	myarray.T1[myarray.sizex-1][j] = myarray.T1[myarray.sizex-2][j];

	}
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
    dy = (j-0.5)*DIMy;
    //set intial conditions for array    
    myarray.T1[i][j] = dy;
    myarray.u1[i][j] = 6.0*U0*dy*(1.0-dy);
    myarray.v1[i][j] = 0;
    //printf("f: %f  dx: %f\n", f, dx);
    }
}
}

//**************************************************************************//
//-----------------------------IE Array Solving-----------------------------//
//**************************************************************************//

void solve_array_IE(carray & myarray, double tmax, double cfl)
{

//double tstep = (cfl*(myarray.DIMx))/2.0;
double tstep = cfl;
double ctime = 0.0;

set_intial_cond(myarray);
set_ghostcells(myarray);

printf("\n\nRunning size: %d time step: %f\n",myarray.sizex,tstep);


int n = 0;
int nt = 10;
double mdiff = 0;

while(ctime < tmax-tstep)
{

//implicit time advance
ctime = ctime+tstep;
compute_Flux(myarray);
solve_LinSys(myarray, tstep, mdiff);
set_ghostcells(myarray);
//print_array(myarray);

if(n >= nt)//status
{
printf("Run: %d time: %f diff %f\n",n,ctime, mdiff);
nt = 10+n;
} 
if(mdiff<0.000001)
break; 

++n;
}

myarray.ctime = ctime;
printf("Solved numeric at %f time with tstep %f\n",ctime, tstep);

}

//-------------------------------LHS approx factor----------------------//

void solve_LinSys(carray & myarray, double tstep, double & mdiff)
{
mdiff = -1.0;
//--solve linear system num 1---//
crow myrow;
//print_array(myarray);
//print_arrayu(myarray);

for(int j = 0; j < myarray.sizey; ++j)
{
//load row data via appox fact then tomas to solve
load_row(myarray, myrow, j, tstep);
solve_thomas(myrow, myarray.sizex);

    for(int i = 0; i < myarray.sizex; ++i)
    {
    myarray.f1[i][j] = myrow.RHS[i];
    }
}

//print_array(myarray);
//print_arrayu(myarray);

//--linear system num 2---//
ccol mycol;
for(int i = 0; i < myarray.sizex; ++i)
{
//load col data via appox fact then tomas to solve
load_col(myarray, mycol, i, tstep);
solve_thomas(mycol, myarray.sizey);

    for(int j = 0; j < myarray.sizey; ++j)
    {
    double temp = mycol.RHS[j];
    if(abs(temp) > mdiff)
    mdiff = abs(temp);
    myarray.f1[i][j] = temp;
    myarray.T1[i][j] = temp + myarray.T1[i][j];
    }
}
//print_arrayu(myarray);
}



//--------------------------Load row for Thomson---------------------------//

void load_row(carray & myarray, crow & myrow, int j, double tstep)
{
double chx = myarray.DIMx;
for(int i = 0; i < myarray.sizex; ++i)
{

double alpha = tstep/(RE*PR*pow(chx,2.0));
double beta = (myarray.u1[i][j]*tstep)/(2.0*chx);
if(i == 0)
{
myrow.LHS[i][0] = 0.0;//apply dirichlet
myrow.LHS[i][1] = 1.0;
myrow.LHS[i][2] = 1.0;
}
else if(i == myarray.sizex-1)
{
myrow.LHS[i][0] =  -1.0;//apply newmann
myrow.LHS[i][1] =  1.0;
myrow.LHS[i][2] =  0.0;
}
else
{
myrow.LHS[i][0] = (-alpha - beta);
myrow.LHS[i][1] = (1.0 + 2.0*alpha);
myrow.LHS[i][2] = (-alpha + beta);
}
//load the flux

myrow.RHS[i] = tstep * myarray.f1[i][j];
}

}



//--------------------------Load column for Thomson---------------------------//

void load_col(carray & myarray, ccol & mycol, int i, double tstep)
{
double chy = myarray.DIMy;
for(int j = 0; j < myarray.sizey; ++j)
{
double alpha = tstep/(RE*PR*pow(chy,2));
double beta = (myarray.v1[i][j]*tstep)/(2.0*chy);


if(j == 0)//apply dirichlet
{
mycol.LHS[j][0] =  0.0;
mycol.LHS[j][1] =  1.0;
mycol.LHS[j][2] =  1.0;
}
else if(j == myarray.sizey-1)//apply dirichlet
{
mycol.LHS[j][0] =  1.0;
mycol.LHS[j][1] =  1.0;
mycol.LHS[j][2] =  0.0;
}
else
{
mycol.LHS[j][0] = (-alpha - beta);
mycol.LHS[j][1] = (1.0 + 2.0*alpha);
mycol.LHS[j][2] = (-alpha + beta);
}
//load the solution from first linsys solve
mycol.RHS[j] = myarray.f1[i][j]; 
}
mycol.RHS[0] = 0.0; 
mycol.RHS[myarray.sizey] = 0.0; 

}


//------------------------------Thomas solving----------------------------------//


void solve_thomas(crow & r, int iSize)
{
  int i;
  /* This next line actually has no effect, but it -does- make clear that
     the values in those locations have no impact. */
  r.LHS[0][0] = r.LHS[iSize-1][2] = 0;
  /* Forward elimination */
  for (i = 0; i < iSize-1; i++) {
    r.LHS[i][2] /= r.LHS[i][1];
    r.RHS[i] /= r.LHS[i][1];
    r.LHS[i+1][1] -= r.LHS[i][2]*r.LHS[i+1][0];
    r.RHS[i+1] -= r.LHS[i+1][0]*r.RHS[i];
  }
  /* Last line of elimination */
  r.RHS[iSize-1] /= r.LHS[iSize-1][1];

  /* Back-substitution */
  for (i = iSize-2; i >= 0; i--) {
    r.RHS[i] -= r.RHS[i+1]*r.LHS[i][2];
  }
  
}

void solve_thomas(ccol & r, int iSize)
{
  int i;
  /* This next line actually has no effect, but it -does- make clear that
     the values in those locations have no impact. */
  r.LHS[0][0] = r.LHS[iSize-1][2] = 0;
  /* Forward elimination */
  for (i = 0; i < iSize-1; i++) {
    r.LHS[i][2] /= r.LHS[i][1];
    r.RHS[i] /= r.LHS[i][1];
    r.LHS[i+1][1] -= r.LHS[i][2]*r.LHS[i+1][0];
    r.RHS[i+1] -= r.LHS[i+1][0]*r.RHS[i];
  }
  /* Last line of elimination */
  r.RHS[iSize-1] /= r.LHS[iSize-1][1];

  /* Back-substitution */
  for (i = iSize-2; i >= 0; i--) {
    r.RHS[i] -= r.RHS[i+1]*r.LHS[i][2];
  }
  
}

//**************************************************************************//
//-----------------------------EE Array Solving-----------------------------//
//**************************************************************************//

void solve_array_EE(carray & myarray, double tmax, double cfl)
{
double tstep = (cfl*(myarray.DIMx))/2.0;
double ctime = 0.0;
double mdiff = 0.0;
set_intial_cond(myarray);
set_ghostcells(myarray);

printf("\n\nRunning size: %d time step: %f\n",myarray.sizex,tstep);

int n = 0;
int nt = 10;

while(ctime < tmax-tstep)
{
  
ctime = ctime+tstep;
compute_Flux(myarray);
time_advance_EE(myarray, tstep, mdiff);

set_ghostcells(myarray);


if(n >= nt)//status
{
printf("Run: %d time: %f diff %f\n",n,ctime, mdiff);
nt = 10+n;
}
if(mdiff<0.000001)
break; 

++n;
}

myarray.ctime = ctime;
printf("Solved numeric at %f time with tstep %f\n",ctime, tstep);
}

//------------------------------EE time advance----------------------------//

void time_advance_EE(carray & myarray, double tstep, double & mdiff)
{
mdiff = -1.0;
for(int j = 1; j < myarray.sizey-1; ++j)
{
    for(int i = 1; i < myarray.sizex-1; ++i)
    {  
    double temp = tstep*(myarray.f1[i][j]); 
    if(abs(temp)> mdiff)
    mdiff = abs(temp);
       
    myarray.T1[i][j] = myarray.T1[i][j] + temp;
    }
}
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
mysurr.ui_jm1 = myarray.u1[i][j-1];
mysurr.ui_jp1 = myarray.u1[i][j+1];
 

mysurr.vi_jm1 = myarray.v1[i][j-1];
mysurr.vi_jp1 = myarray.v1[i][j+1];
mysurr.vim1_j = myarray.v1[i-1][j];
mysurr.vip1_j = myarray.v1[i+1][j];
}


//-----------------------calculate flux for new cell--------------------//

double calc_newcell(carray & myarray, surr & s1)
{
double chx = myarray.DIMx;
double chy = myarray.DIMy;
//float chy = DIM1;

double a = (s1.uip1_j * s1.Tip1_j  -  s1.uim1_j * s1.Tim1_j)/(2.0);
double b = (s1.Tip1_j  - 2.0*s1.Ti_j + s1.Tim1_j)/(chx*(RE*PR));
double c = (s1.vi_jp1 * s1.Ti_jp1  -  s1.vi_jm1 * s1.Ti_jm1)/(2.0);
double d = (s1.Ti_jp1  - 2.0*s1.Ti_j + s1.Ti_jm1)/(chy*(RE*PR));

double e = (2.0*pow((s1.uip1_j-s1.uim1_j)/(2.0*chx),2)+2.0*pow((s1.vi_jp1-s1.vi_jm1)/(2.0*chy),2));
double f = pow(((s1.vip1_j-s1.vim1_j)/(2.0*chx))+((s1.ui_jp1-s1.ui_jm1)/(2.0*chy)), 2);
double source = (EC/RE)*(e+f);

double newcell = 1.0*((-1.0/chx)*(a-b)   +   (-1.0/chy)*(c-d))   +   source;
return newcell;
}



//**************************************************************************//
//---------------------------Error Checking---------------------------------//
//**************************************************************************//
//--------------------------Get descrete error----------------------------//

void get_discrete_Error(carray ray1, carray ray2, carray ray3)
{
//Calculating error as described in paper "procedure for estimation and reporting of uncertainty due to discretization in CFD applications"//

printf("\nCalculating Error...\n");

float h1 = ray1.DIMx;
float h2 = ray2.DIMx;
float h3 = ray3.DIMx;

/*
float sol1 = ray1.T1[40][20];
float sol2 = ray2.T1[20][10];
float sol3 = ray3.T1[10][5];
*/
float sol1 = ray1.T1[44][20];
float sol2 = ray2.T1[22][10];
float sol3 = ray3.T1[11][5];



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
//-------------------------------------l2norm--------------------------//

double get_l2norm(carray & myarray, carray myarray2)
{
double l2sum =0;

double sx = myarray.sizex-2;
double sy = myarray.sizey-2;

for(int j = 1; j < myarray.sizey-1; ++j)
{	
	for(int i = 1; i < myarray.sizex-1; ++i)
	{
	double P = myarray.T1[i][j];
	double T = myarray2.T1[i][j];
	l2sum =  l2sum + pow((P-T),2);

	}


}

double l2 = sqrt(l2sum/(sx*sy));
cout << setprecision(20) << fixed << "L2 norm: " << l2 << "\n";
return l2;
}

//-------------------------------------l2norm between diff--------------------------//

double get_l1normD(carray & myarray, carray myarray2)
{
double l2sum =0;

double sx = myarray.sizex-2;
double sy = myarray.sizey-2;

for(int j = 1; j < myarray.sizey-1; ++j)
{	
	for(int i = 1; i < myarray.sizex-1; ++i)
	{
	double P = myarray.f1[i][j];
	double T = myarray2.f1[i*2][j*2];
	l2sum =  l2sum + abs(P-T);

	}


}

double l2 = l2sum/(sx*sy);
cout << setprecision(9) << fixed << "L1 norm: " << l2 << "\n";
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
	
	
	myarray.T1[i][j] =  ( dy+((3.0/4.0)*PR*EC*(pow(U0,2))*(1.0-pow((1.0-2.0*dy),4))) );
	
	}

}

printf("setting analytic at %f time\n",ctime);
}

//**************************************************************************//
//---------------------------Print Functions--------------------------------//
//**************************************************************************//


//--------------------------Print array in terminal----------------------------//

void print_array(carray & myarray)
{
cout << "Solution:\n       |";
for(int i = 0; i < myarray.sizex; ++i)
		{
		if(i < 10)
		cout << "   i: " << i <<"|";
		if(i > 9)
		cout << "   i:" << i <<"|";
        }
cout << "\n";
	for(int j = 0; j < myarray.sizey; ++j)
	{
	if(j > 9)
	cout << "\nj:" << j << "|  |";
	if(j < 10)	
	cout << "\nj: " << j << "|  |";
		for(int i = 0; i < myarray.sizex; ++i)
		{
		if(myarray.T1[i][j] >= 0)
		cout << setprecision(5) << fixed << myarray.T1[i][j] <<"|";
		if(myarray.T1[i][j] < 0)
		cout << setprecision(4) << fixed << myarray.T1[i][j] <<"|";
		}
	
	}
cout << "\n\n";
}

//--------------------------Print u array in terminal----------------------------//

void print_arrayu(carray & myarray)
{

cout << "Flux:\n       |";
for(int i = 0; i < myarray.sizex; ++i)
		{
		if(i < 10)
		cout << "   i: " << i <<"|";
		if(i > 9)
		cout << "   i:" << i <<"|";
        }
cout << "\n";
	for(int j = 0; j < myarray.sizey; ++j)
	{
	if(j > 9)
	cout << "\nj:" << j << "|  |";
	if(j < 10)	
	cout << "\nj: " << j << "|  |";
		for(int i = 0; i < myarray.sizex; ++i)
		{
		if(myarray.f1[i][j] >= 0)
		cout << setprecision(5) << fixed << myarray.f1[i][j] <<"|";
		if(myarray.f1[i][j] < 0)
		cout << setprecision(4) << fixed << myarray.f1[i][j] <<"|";
		}
	
	}
cout << "\n\n";
}


//-------------------------------print row------------------------------------//

void print_row(crow & myrow, carray & myarray)
{
cout << "\nprint rows\n";
for(int i = 0; i < myarray.sizex; ++i)
{
cout << "|" << myrow.LHS[i][0] << "|" << myrow.LHS[i][1] << "|" << myrow.LHS[i][2] << "|   |" << myrow.RHS[i] << "|\n";

}

}



void print_col(ccol & mycol, carray & myarray)
{
cout << "\nprint rows\n";
for(int i = 0; i < myarray.sizey; ++i)
{
cout << "|" << mycol.LHS[i][0] << "|" << mycol.LHS[i][1] << "|" << mycol.LHS[i][2] << "|   |" << mycol.RHS[i] << "|\n";

}

}

/*
void set_ghostcells(carray & myarray)
{
double DIMx = myarray.DIMx;
double DIMy = myarray.DIMy;
double dx = 0.0;
double dy = 0.0;
    //set ghost cells top/bottom
	for(int i = 0; i < myarray.sizex; ++i)
	{
	myarray.T1[i][0] = 2.0*(0.0) - myarray.T1[i][1];
	myarray.T1[i][myarray.sizey-1] =2.0*(1.0)-myarray.T1[i][myarray.sizey-2];
	}	
	
    //set ghost cells inflow/outflow	
	for(int j = 0; j < myarray.sizey; ++j)
	{ 
    dy = (j-0.5)*DIMy;
    
	myarray.T1[0][j] = 2.0*( dy+((3.0/4.0)*PR*EC*(pow(U0,2))*(1.0-pow((1.0-2.0*dy),4))) ) - myarray.T1[1][j];
	
	
	myarray.T1[myarray.sizex-1][j] = myarray.T1[myarray.sizex-2][j];

	}
}
*/


