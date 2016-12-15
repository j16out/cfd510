#include "numerical.hpp"





//==========================================================================//
//---------------------------Setting Array----------------------------------//
//==========================================================================//


//-----------set array size (working area excluding ghost)-----------------//

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




//-----------------------------zero array----------------------------------//

void set_zero(carray & myarray)
{
	for(int j = 0; j < myarray.sizey; ++j)
	{
		for(int i = 0; i < myarray.sizex; ++i)
		{
		//set solution array zero
		myarray.s1[i][j].P = 0;
        myarray.s1[i][j].u = 0;
        myarray.s1[i][j].v = 0;
        //set flux array zero
        myarray.f1[i][j].P = 0;
        myarray.f1[i][j].u = 0;
        myarray.f1[i][j].v = 0;

		}
	}
}



//-------------------------set intial condition----------------------------//

void set_init_cond(carray & myarray)
{
double DIMx = myarray.DIMx;
double DIMy = myarray.DIMy;


for(int j = 0; j < myarray.sizey; ++j)
{

    for(int i = 0; i < myarray.sizex; ++i)
    {
    double x = (i-0.5)*DIMx;
	double y = (j-0.5)*DIMy;
    
    myarray.s1[i][j].P = P0*cos(PI*x)*cos(PI*y);
    myarray.s1[i][j].u = U0*sin(PI*x)*sin(2.0*PI*y);
    myarray.s1[i][j].v = V0*sin(2.0*PI*x)*sin(PI*y);;
    //printf("f: %f  dx: %f\n", f, dx);
    }

}
}

//==========================================================================//
//-----------------------------IE Array Solving-----------------------------//
//==========================================================================//

void solve_array_IE(carray & myarray, double tmax, double cfl)
{

//double tstep = (cfl*(myarray.DIMx))/2.0;
double tstep = cfl;
double ctime = 0.0;

set_init_cond(myarray);
//set_ghostcells(myarray);

update_flux(myarray);

}

//---------------------------LHS approx factor----------------------------//













//==========================================================================//
//---------------------------Compute Flux-----------------------------------//
//==========================================================================//


void update_flux(carray & myarray)
{
surr mysurr;
vec ftemp;

for(int j = 1; j < myarray.sizey-1; ++j)
{
    for(int i = 1; i < myarray.sizex-1; ++i)
    {
    //----get surrounding cells and compute new cell----//
    get_nsurcells(myarray, i, j, mysurr);    
    calc_flux(myarray, mysurr, ftemp); 
    
    //-----update current cell----//  
    myarray.f1[i][j] = ftemp;   
    }

}
}

//--------------------------get surrounding cells-----------------------//

void get_nsurcells(carray & myarray, int i, int j, surr & mysurr)
{

mysurr.Pim1_j = myarray.s1[i-1][j].P;
mysurr.Pip1_j = myarray.s1[i+1][j].P;
mysurr.Pi_j = myarray.s1[i][j].P; 
mysurr.Pi_jm1 = myarray.s1[i][j-1].P;
mysurr.Pi_jp1 = myarray.s1[i][j+1].P;

mysurr.uim1_j = myarray.s1[i-1][j].u;
mysurr.uip1_j = myarray.s1[i+1][j].u;
mysurr.ui_j = myarray.s1[i][j].u; 
mysurr.ui_jm1 = myarray.s1[i][j-1].u;
mysurr.ui_jp1 = myarray.s1[i][j+1].u;
 

mysurr.vi_jm1 = myarray.s1[i][j-1].v;
mysurr.vi_jp1 = myarray.s1[i][j+1].v;
mysurr.vi_j = myarray.s1[i][j].v; 
mysurr.vim1_j = myarray.s1[i-1][j].v;
mysurr.vip1_j = myarray.s1[i+1][j].v;
}


//-----------------------calculate flux for new cell--------------------//

void calc_flux(carray & myarray, surr & s1, vec & ftemp)
{
double chx = myarray.DIMx;
double chy = myarray.DIMy;
ftemp.P = 0;
ftemp.u = 0;
ftemp.v = 0;

//x direction

double fp = (((s1.uip1_j + s1.ui_j)/(2.0*BETA)) - ((s1.ui_j + s1.uim1_j)/(2.0*BETA))) /chx;

double fu = (   (pow((s1.uip1_j+s1.ui_j)/2.0, 2) + ((s1.Pip1_j+s1.Pi_j)/2.0) - ((s1.uip1_j-s1.ui_j)/(chx*RE)))
            - (pow((s1.ui_j+s1.uim1_j)/2.0, 2) + ((s1.Pi_j+s1.Pim1_j)/2.0) - ((s1.ui_j-s1.uim1_j)/(chx*RE)))  )/chx; 
                        
double fv =   (   (((s1.uip1_j+s1.ui_j)/2.0)*((s1.vip1_j+s1.vi_j)/2.0) - ((s1.vip1_j-s1.vi_j)/(chx*RE)))
            - (((s1.ui_j+s1.uim1_j)/2.0)*((s1.vi_j+s1.vim1_j)/2.0) - ((s1.vi_j-s1.vim1_j)/(chx*RE)))  )/chx;

double gp = (((s1.vi_jp1 + s1.vi_j)/(2.0*BETA)) - ((s1.vi_j + s1.vi_jm1)/(2.0*BETA))) /chy;

                        
double gu =   (   (((s1.ui_jp1+s1.ui_j)/2.0)*((s1.vi_jp1+s1.vi_j)/2.0) - ((s1.ui_jp1-s1.ui_j)/(chy*RE)))
            - (((s1.ui_j+s1.ui_jm1)/2.0)*((s1.vi_j+s1.vi_jm1)/2.0) - ((s1.ui_j-s1.ui_jm1)/(chy*RE)))  )/chy;

double gv = (   (pow((s1.vi_jp1+s1.vi_j)/2.0, 2) + ((s1.Pi_jp1+s1.Pi_j)/2.0) - ((s1.vi_jp1-s1.vi_j)/(chy*RE)))
            - (pow((s1.vi_j+s1.vi_jm1)/2.0, 2) + ((s1.Pi_j+s1.Pi_jm1)/2.0) - ((s1.vi_j-s1.vi_jm1)/(chy*RE)))  )/chy ;
            
            
ftemp.P = -fp-gp;
ftemp.u = -fu-gu;
ftemp.v = -fv-gv;
}



//==========================================================================//
//---------------------------Error Checking---------------------------------//
//==========================================================================//



//--------------------------Get descrete error----------------------------//
/*
void get_discrete_Error(carray ray1, carray ray2, carray ray3)
{
//Calculating error as described in paper "procedure for estimation and reporting of uncertainty due to discretization in CFD applications"//

printf("\nCalculating Error...\n");

double h1 = ray1.DIMx;
double h2 = ray2.DIMx;
double h3 = ray3.DIMx;


double sol1 = ray1.P1[40][20];
double sol2 = ray2.P1[20][10];
double sol3 = ray3.P1[10][5];

double sol1 = 1.114209;
double sol2 = 1.112606;
double sol3 = 1.111794;




printf("h1: %f \nh2: %f \nh3: %f, \nsol1: %f \nsol2: %f \nsol3: %f\n",h1, h2, h3, sol1, sol2, sol3);

double r21 = h2/h1;
double r32 = h3/h2;

printf("\nr32: %f \nr21: %f\n",r32, r21);

double e32 = sol3-sol2;
double e21 = sol2-sol1;

double s = (e32/e21);
if(s >= 0)
s = 1;
else
s = -1;

double p_n = 0;
double p = (1/log(r21))*(abs(log(abs(e32/e21))+0));

printf("intial guess: %f \n", p);

double diff = 1;

	while(diff > 0.0000001)
	{

	double p_n = (1/log(r21))*(abs(log(abs(e32/e21))+log((pow(r21,p)-s)/(pow(r32,p)-s)) ));
	diff = abs(p_n -p);
	//printf("p_n: %f p: %f diff: %f\n",p_n, p, diff);

	p = p_n;
	}
 
//
double sol_ext21 = (pow(r21, p)*sol1-sol2)/(pow(r21,p)-1.0);
double sol_ext32 = (pow(r32, p)*sol2-sol3)/(pow(r32,p)-1.0);

printf("order: %f \nphi_ext21: %f \nphi_ext32 %f\n",p, sol_ext21, sol_ext32);

double ea21 = abs((sol1-sol2)/sol1);

double e_ext21 = abs((sol_ext21-sol1)/sol_ext21);

double GCI_21 = (1.25*ea21)/(pow(r21,p)-1.0);


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
	double P = myarray.f1[i][j];
	double T = myarray2.f1[i][j];
	l2sum =  l2sum + pow((P-T),2);

	}


}

double l2 = sqrt(l2sum/(sx*sy));
cout << setprecision(8) << fixed << "L2 norm: " << l2 << "\n";
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
*/
//----------------------------Set a Analytical Solution------------------------------//
void set_analytic(carray & myarray)
{
double DIMx = myarray.DIMx;
double DIMy = myarray.DIMy;



for(int j = 1; j < myarray.sizey-1; ++j)
{

	for(int i = 1; i < myarray.sizex-1; ++i)
	{
	double x = (i-0.5)*DIMx;
	double y = (j-0.5)*DIMy;
	
	double cx = cos(PI*x);
    double sx = sin(PI*x);
    double cy = cos(PI*y);
    double sy = sin(PI*y);
	double c2x = cos(2.0*PI*x);
    double s2x = sin(2.0*PI*x);
    double c2y = cos(2.0*PI*y);
    double s2y = sin(2.0*PI*y);
	
	
	myarray.f1[i][j].P =  (-PI/BETA)*(U0*cx*s2y   +   V0*s2x*cy);
	myarray.f1[i][j].u =  (P0*PI*sx*cy)  
	                       -  (pow(U0,2)*PI*s2x*pow(s2y,2)) 
	                       -  (U0*V0*PI*sx*s2x*(cy*s2y+2.0*c2y*sy))
	                       -U0*((5.0*pow(PI,2)*sx*s2y)/(RE));
	                       
    myarray.f1[i][j].v =  (P0*PI*cx*sy)  
	                       -  (pow(V0,2)*PI*pow(s2x,2)*s2y) 
	                       -  (U0*V0*PI*sy*s2y*(cx*s2x+2.0*c2x*sx))
	                       -V0*((5.0*pow(PI,2)*s2x*sy)/(RE));
	}

}

printf("setting analytic\n");
}

//==========================================================================//
//---------------------------Print Functions--------------------------------//
//==========================================================================//




//--------------------------Print array in terminal----------------------------//

void print_array_sP(carray & myarray)
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
		if(myarray.s1[i][j].P >= 0)
		cout << setprecision(5) << fixed << myarray.s1[i][j].P <<"|";
		if(myarray.s1[i][j].P < 0)
		cout << setprecision(4) << fixed << myarray.s1[i][j].P <<"|";
		}
	
	}
cout << "\n\n";
}

//--------------------------Print u array in terminal----------------------------//

void print_array_fP(carray & myarray)
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
		if(myarray.f1[i][j].P >= 0)
		cout << setprecision(5) << fixed << myarray.f1[i][j].P <<"|";
		if(myarray.f1[i][j].P < 0)
		cout << setprecision(4) << fixed << myarray.f1[i][j].P <<"|";
		}
	
	}
cout << "\n\n";
}

/*
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
*/
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
	myarray.P1[i][0] = 2.0*(0.0) - myarray.P1[i][1];
	myarray.P1[i][myarray.sizey-1] =2.0*(1.0)-myarray.P1[i][myarray.sizey-2];
	}	
	
    //set ghost cells inflow/outflow	
	for(int j = 0; j < myarray.sizey; ++j)
	{ 
    dy = (j-0.5)*DIMy;
    
	myarray.P1[0][j] = 2.0*( dy+((3.0/4.0)*PR*EC*(pow(U0,2))*(1.0-pow((1.0-2.0*dy),4))) ) - myarray.P1[1][j];
	
	
	myarray.P1[myarray.sizex-1][j] = myarray.P1[myarray.sizex-2][j];

	}
}
*/


