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
	myarray.T1[i][0] = 2.0*(0.0) - myarray.T1[i][1];
	myarray.T1[i][myarray.sizey] =2.0*(1.0)-myarray.T1[i][myarray.sizey-1];
	}	
	
    //set ghost cells inflow/outflow	
	for(int j = 1; j < myarray.sizey-1; ++j)
	{ 
    dy = (j-0.5)*DIMy;
    
	myarray.T1[0][j] = 2.0*( dy+((3.0/4.0)*PR*EC*(pow(U0,2))*(1.0-pow((1.0-2.0*dy),4))) ) - myarray.T1[1][j];
	
	
	myarray.T1[myarray.sizex][j] = myarray.T1[myarray.sizex-1][j];//set everything to zero

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

double tstep = (cfl*(myarray.DIMx))/2.0;
double ctime = 0.0;

set_intial_cond(myarray);
set_ghostcells(myarray);


printf("\n\nRunning size: %d time step: %f\n",myarray.sizex,tstep);

int n = 0;
int nt = 1000;

while(ctime < tmax-tstep)
{

if(n >= nt)//status
{
printf("Run: %d time: %f\n",n,ctime);
nt = 1000+n;
}   

//implicit time advance
ctime = ctime+tstep;
compute_Flux(myarray);
solve_LinSys(myarray, tstep);
set_ghostcells(myarray);

print_array(myarray);
print_arrayu(myarray);

++n;
}

myarray.ctime = ctime;
printf("Solved numeric at %f time\n",ctime);
}

//-------------------------------LHS approx factor----------------------//

void solve_LinSys(carray & myarray, double tstep)
{

//--linear system num 1---//
crow myrow;

print_array(myarray);
print_arrayu(myarray);
for(int j = 0; j < myarray.sizey; ++j)
{

load_row(myarray, myrow, j, tstep);
solve_thomas(myrow, myarray.sizex);

    for(int i = 0; i < myarray.sizex; ++i)
    {
    myarray.f1[i][j] = myrow.RHS[i];
    }
}
print_array(myarray);
print_arrayu(myarray);

//--linear system num 2---//
ccol mycol;
for(int i = 0; i < myarray.sizex; ++i)
{
load_col(myarray, mycol, i, tstep);
solve_thomas(mycol, myarray.sizey);

    for(int j = 0; j < myarray.sizey; ++j)
    {
    myarray.f1[i][j] = mycol.RHS[i];
    myarray.T1[i][j] = mycol.RHS[j] + myarray.T1[i][j];
    }
}



}



//--------------------------Load row for Thomson---------------------------//

void load_row(carray & myarray, crow & myrow, int j, double tstep)
{
double chx = myarray.DIMx;
for(int i = 0; i < myarray.sizex; ++i)
{
double alpha = tstep/(RE*PR*pow(chx,2));
double beta = (myarray.u1[i][j]*tstep)/(2.0*chx);

myrow.LHS[i][0] = (-alpha - beta);
myrow.LHS[i][1] = (1.0 + 2.0*alpha);
myrow.LHS[i][2] = (-alpha + beta);
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

mycol.LHS[j][0] = (-alpha - beta);
mycol.LHS[j][1] = (1.0 + 2.0*alpha);
mycol.LHS[j][2] = (-alpha + beta);
//load the solution from first linsys solve
mycol.RHS[j] = myarray.f1[i][j]; 
}

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

set_intial_cond(myarray);
set_ghostcells(myarray);

printf("\n\nRunning size: %d time step: %f\n",myarray.sizex,tstep);

int n = 0;
int nt = 1000;

while(ctime < tmax-tstep)
{

if(n >= nt)//status
{
printf("Run: %d time: %f\n",n,ctime);
nt = 1000+n;
}   

ctime = ctime+tstep;
compute_Flux(myarray);
time_advance_EE(myarray, tstep);
set_ghostcells(myarray);

++n;
}

myarray.ctime = ctime;
printf("Solved numeric at %f time\n",ctime);
}

//------------------------------EE time advance----------------------------//

void time_advance_EE(carray & myarray, double tstep)
{
for(int j = 1; j < myarray.sizey-1; ++j)
{
    for(int i = 1; i < myarray.sizex-1; ++i)
    {
      myarray.T1[i][j] = myarray.T1[i][j]+tstep*(myarray.f1[i][j]);
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


double get_l2norm(carray & myarray, carray myarray2)
{
double l2sum =0;

double sx = myarray.sizex-2;
double sy = myarray.sizey-2;

for(int j = 1; j < myarray.sizey-1; ++j)
{	
	for(int i = 1; i < myarray.sizex-1; ++i)
	{
	double P = myarray.f1[i][j];
	double T = myarray2.f1[i][j];
	l2sum =  l2sum + pow((P-T),2);

	}


}

double l2 = sqrt(l2sum/(sx*sy));
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
	
	double d = pow((U0*sin(PI*dx)+V0*cos(PI*dy)), 2);
	double source = (EC/RE)*((2.0*pow((U0*PI*cos(PI*dx)*dy),2)+(2.0*pow((V0*PI*sin(PI*dy)*dx),2)))+d);
	myarray.f1[i][j] =  a + b + c + source;
	
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
		cout << " i: " << i <<"|";
		if(i > 9)
		cout << " i:" << i <<"|";
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
		cout << setprecision(3) << fixed << myarray.T1[i][j] <<"|";
		if(myarray.T1[i][j] < 0)
		cout << setprecision(2) << fixed << myarray.T1[i][j] <<"|";
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
		cout << " i: " << i <<"|";
		if(i > 9)
		cout << " i:" << i <<"|";
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
		cout << setprecision(3) << fixed << myarray.f1[i][j] <<"|";
		if(myarray.f1[i][j] < 0)
		cout << setprecision(2) << fixed << myarray.f1[i][j] <<"|";
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


