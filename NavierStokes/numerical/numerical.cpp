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
	cout << "Array size to big, array not set***" << "\n";

}




//-----------------------------zero array----------------------------------//

void set_zero(carray & myarray)
{
	for(int j = 0; j < myarray.sizey; ++j)
	{
		for(int i = 0; i < myarray.sizex; ++i)
		{
		//set solution array zero
		myarray.s1[i][j].P = 0.0;
        myarray.s1[i][j].u = 0.0;
        myarray.s1[i][j].v = 0.0;
        //set flux array zero
        myarray.f1[i][j].P = 0.0;
        myarray.f1[i][j].u = 0.0;
        myarray.f1[i][j].v = 0.0;
        

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

//--------------------------set ghost cells for Energy----------------------------//


void set_ghostcells(carray & myarray)
{
double DIMx = myarray.DIMx;
double DIMy = myarray.DIMy;
double dx = 0.0;
double dy = 0.0;
double Uwall = 0.0;
    //set ghost cells top/bottom
	for(int i = 0; i < myarray.sizex; ++i)
	{
	//top
	myarray.s1[i][0].P = myarray.s1[i][1].P ;
	myarray.s1[i][0].u = (2.0*UW)-myarray.s1[i][1].u ;
	myarray.s1[i][0].v = -myarray.s1[i][1].v ;
	//bottom
	myarray.s1[i][myarray.sizey-1].P = myarray.s1[i][myarray.sizey-2].P ;
	myarray.s1[i][myarray.sizey-1].u = -myarray.s1[i][myarray.sizey-2].u ;
	myarray.s1[i][myarray.sizey-1].v = -myarray.s1[i][myarray.sizey-2].v ;
	
	}	
	
    //set ghost cells inflow/outflow	
	for(int j = 0; j < myarray.sizey; ++j)
	{ 
	//left
    myarray.s1[0][j].P = myarray.s1[1][j].P;
    myarray.s1[0][j].u = -myarray.s1[1][j].u;
    myarray.s1[0][j].v = -myarray.s1[1][j].v;
    //right
    myarray.s1[myarray.sizex-1][j].P = myarray.s1[myarray.sizex-2][j].P ;
    myarray.s1[myarray.sizex-1][j].u = -myarray.s1[myarray.sizex-2][j].u ;
    myarray.s1[myarray.sizex-1][j].v = -myarray.s1[myarray.sizex-2][j].v ;

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

//set intial conditions/boundaries
set_init_cond(myarray);
set_ghostcells(myarray);

//print_array_sP(myarray);


printf("\n\nRunning size: %d time step: %f\n",myarray.sizex,tstep);


int n = 0;
int nt = 10;
double mdiff = -1.0;

while(ctime < tmax-tstep)
{
ctime = ctime+tstep;

//update flux
update_flux(myarray, tstep);
//print_array_fP(myarray, 1.0);

//implicit time advance
solve_LinSys(myarray, tstep, mdiff);

//update boundaries
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




//---------------------------LHS approx factor----------------------------//

//==========================================================================//
//---------------------------LHS Jacobian-----------------------------------//
//==========================================================================//




void solve_LinSys(carray & myarray, double tstep, double & mdiff)
{
mdiff = -1.0;
//--linear system num 1---//
crow myrow;
//print_arrayu(myarray);

for(int j = 1; j < myarray.sizey-1; ++j)
{
load_row(myarray, myrow, j, tstep);
solve_block_thomas(myarray, myrow, myarray.sizex, j);

}

//print_array_fP(myarray, 1.0);

ccol mycol;
for(int i = 1; i < myarray.sizex-1; ++i)
{
load_col(myarray, mycol, i, tstep);
solve_block_thomas(myarray, mycol, myarray.sizey, i, mdiff);
}
//print_array_sP(myarray);
//print_array_fP(myarray, 1.0);

update_sol(myarray);
//print_array_sP(myarray);
}

void update_sol(carray & myarray)
{
 for(int j = 0; j < myarray.sizey; ++j)
{
    for(int i = 0; i < myarray.sizex; ++i)
    {
    myarray.s1[i][j].P = myarray.f1[i][j].P + myarray.s1[i][j].P;
    myarray.s1[i][j].u = myarray.f1[i][j].u + myarray.s1[i][j].u;
    myarray.s1[i][j].v = myarray.f1[i][j].v + myarray.s1[i][j].v;
    }
}

}
//--------------------------Load row for Thomson---------------------------//

void load_row(carray & myarray, crow & myrow, int j, double tstep)
{
LHScX c1;

for(int i = 0; i < myarray.sizex; ++i)
{
calc_LHS_constX(myarray, c1, i, j, tstep);

if(i == 0)
{
set_wall(c1, 1);
myrow.LHS[i] = c1;
myrow.RHS[i].P = 0.0;
myrow.RHS[i].u = 0.0;
myrow.RHS[i].v = 0.0;
}
else if(i == myarray.sizex-1)
{
set_wall(c1, 0);
myrow.LHS[i] = c1;
myrow.RHS[i].P = 0.0;
myrow.RHS[i].u = 0.0;
myrow.RHS[i].v = 0.0;
}
else
{
myrow.LHS[i] = c1;
myrow.RHS[i].P = myarray.f1[i][j].P;
myrow.RHS[i].u = myarray.f1[i][j].u;
myrow.RHS[i].v = myarray.f1[i][j].v;
}


}

}


void load_col(carray & myarray, ccol & mycol, int i, double tstep)
{
LHScY c1;

for(int j = 0; j < myarray.sizey; ++j)
{
calc_LHS_constY(myarray, c1, i, j, tstep);

if(j == 0)
{
set_wall(c1, 1);
mycol.LHS[j] = c1;
mycol.RHS[j].P = 0.0;
mycol.RHS[j].u = 0.0;
mycol.RHS[j].v = 0.0;
}
else if(j == myarray.sizey-1)
{
set_wall(c1, 0);
mycol.LHS[j] = c1;
mycol.RHS[j].P = 0.0;
mycol.RHS[j].u = 0.0;
mycol.RHS[j].v = 0.0;
}
else
{
mycol.LHS[j] = c1;
mycol.RHS[j].P = myarray.f1[i][j].P;
mycol.RHS[j].u = myarray.f1[i][j].u;
mycol.RHS[j].v = myarray.f1[i][j].v;
}


}

}




//-------------------------------get x constants----------------------------------//

void calc_LHS_constX(carray & a1, LHScX & c1, int i, int j, double tstep)
{
surr s1;
get_nsurcells(a1, i, j, s1); 

double chx = a1.DIMx;
double chy = a1.DIMy;

tstep = -tstep;

c1.Ax.rP.P = 0.0;
c1.Ax.rP.u = tstep*(1.0/chx)*(1.0/2.0);
c1.Ax.rP.v = 0.0;
c1.Ax.ru.P = tstep*(1.0/chx)*(1.0/(2.0*BETA));
c1.Ax.ru.u = tstep*(1.0/chx)*( ((s1.ui_j+s1.uim1_j)/2.0) + (1.0/(RE*chx)) );
c1.Ax.ru.v = tstep*(1.0/chx)*( ((s1.vi_j+s1.vim1_j)/4.0) );
c1.Ax.rv.P = 0.0;
c1.Ax.rv.u = 0.0;
c1.Ax.rv.v = tstep*(1.0/chx)*( ((s1.ui_j+s1.uim1_j)/4.0) + (1.0/(RE*chx)) );


c1.Bx.rP.P = 1.0; 
c1.Bx.rP.u = 0.0;
c1.Bx.rP.v = 0.0;
c1.Bx.ru.P = 0.0;
c1.Bx.ru.u = 1.0 + (tstep*(-1.0/chx)*( ((s1.ui_j+s1.uip1_j)/2.0) - ((s1.ui_j+s1.uim1_j)/2.0) + (2.0/(RE*chx)) ));
c1.Bx.ru.v = tstep*(-1.0/chx)*( ((s1.vi_j+s1.vip1_j)/4.0) - ((s1.vi_j+s1.vim1_j)/4.0));
c1.Bx.rv.P = 0.0;
c1.Bx.rv.u = 0.0; 
c1.Bx.rv.v = 1.0 + (tstep*(-1.0/chx)*( ((s1.ui_j+s1.uip1_j)/4.0) - ((s1.ui_j+s1.uim1_j)/4.0)+ (2.0/(RE*chx)) ));


c1.Cx.rP.P = 0.0;
c1.Cx.rP.u = tstep*(-1.0/chx)*(1.0/2.0);
c1.Cx.rP.v = 0.0;
c1.Cx.ru.P = tstep*(-1.0/chx)*(1.0/(2.0*BETA));
c1.Cx.ru.u = tstep*(-1.0/chx)*( ((s1.ui_j+s1.uip1_j)/2.0) - (1.0/(RE*chx)) );
c1.Cx.ru.v = tstep*(-1.0/chx)*( ((s1.vi_j+s1.vip1_j)/4.0) );
c1.Cx.rv.P = 0.0;
c1.Cx.rv.u = 0.0;
c1.Cx.rv.v = tstep*(-1.0/chx)*( ((s1.ui_j+s1.uip1_j)/4.0) - (1.0/(RE*chx)) );


}

//-------------------------------get y constants----------------------------------//

void calc_LHS_constY(carray & a1, LHScY & c2, int i, int j, double tstep)
{
surr s1;
get_nsurcells(a1, i, j, s1); 

double chx = a1.DIMx;
double chy = a1.DIMy;

tstep = -tstep;

c2.Ay.rP.P = 0.0;
c2.Ay.rP.u = 0.0;
c2.Ay.rP.v = tstep*(1.0/chy)*(1.0/2.0);
c2.Ay.ru.P = 0.0;
c2.Ay.ru.u = tstep*(1.0/chy)*( ((s1.vi_j+s1.vi_jm1)/4.0) + (1.0/(RE*chy)) );
c2.Ay.ru.v = 0.0;
c2.Ay.rv.P = tstep*(1.0/chy)*(1.0/(2.0*BETA));
c2.Ay.rv.u = tstep*(1.0/chy)*( ((s1.ui_j+s1.ui_jm1)/4.0) );
c2.Ay.rv.v = tstep*(1.0/chy)*( ((s1.vi_j+s1.vi_jm1)/2.0) + (1.0/(RE*chy)) );

c2.By.rP.P = 1.0;
c2.By.rP.u = 0.0;
c2.By.rP.v = 0.0;
c2.By.ru.P = 0.0;
c2.By.ru.u = 1.0 + tstep*(-1.0/chy)*( ((s1.vi_j+s1.vi_jp1)/4.0) - ((s1.vi_j+s1.vi_jm1)/4.0) + (2.0/(RE*chy)) );
c2.By.ru.v = 0.0;
c2.By.rv.P = 0.0;
c2.By.rv.u = tstep*(-1.0/chy)*( ((s1.ui_j+s1.ui_jp1)/4.0) - ((s1.ui_j+s1.ui_jm1)/4.0) );
c2.By.rv.v = 1.0 + tstep*(-1.0/chy)*( ((s1.vi_j+s1.vi_jp1)/2.0) - ((s1.vi_j+s1.vi_jm1)/2.0) + (2.0/(RE*chy)) );


c2.Cy.rP.P = 0.0;
c2.Cy.rP.u = 0.0;
c2.Cy.rP.v = tstep*(-1.0/chy)*(1.0/2.0);
c2.Cy.ru.P = 0.0;
c2.Cy.ru.u = tstep*(-1.0/chy)*( ((s1.vi_j+s1.vi_jp1)/4.0) - (1.0/(RE*chy)) );
c2.Cy.ru.v = 0.0;
c2.Cy.rv.P = tstep*(-1.0/chy)*(1.0/(2.0*BETA));
c2.Cy.rv.u = tstep*(-1.0/chy)*( ((s1.ui_j+s1.ui_jp1)/4.0) );
c2.Cy.rv.v = tstep*(-1.0/chy)*( ((s1.vi_j+s1.vi_jp1)/2.0) - (1.0/(RE*chy)) );

}





//==========================================================================//
//-----------------------RHS Compute Flux-----------------------------------//
//==========================================================================//


void update_flux(carray & myarray, double tstep)
{
surr mysurr;
vec ftemp;
double artP = 1.0;
for(int j = 1; j < myarray.sizey-1; ++j)
{
    for(int i = 1; i < myarray.sizex-1; ++i)
    {
    //----get surrounding cells and compute new cell----//
    get_nsurcells(myarray, i, j, mysurr);    
    calc_flux(myarray, mysurr, ftemp, tstep); 
    //-----update current cell----//  
    myarray.f1[i][j].P = tstep * ftemp.P;
    myarray.f1[i][j].u = tstep * ftemp.u;
    myarray.f1[i][j].v = tstep * ftemp.v;   
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

void calc_flux(carray & myarray, surr & s1, vec & ftemp, double tstep)
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
  
double artP = ARTVIS*((s1.Pip1_j-2.0*s1.Pi_j+s1.Pim1_j)/(2.0*pow(chx,2)) + (s1.Pi_jp1-2.0*s1.Pi_j+s1.Pi_jm1)/(2.0*pow(chy,2)))*chx*chy;            
            
ftemp.P = (-fp-gp)+(artP/tstep);
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
*/
void get_l2norm(carray & myarray, carray myarray2, cdata & mydata)
{
double l2sumP =0;
double l2sumu =0;
double l2sumv =0;


double sx = myarray.sizex-2;
double sy = myarray.sizey-2;

for(int j = 1; j < myarray.sizey-1; ++j)
{	
	for(int i = 1; i < myarray.sizex-1; ++i)
	{
	double P = myarray.f1[i][j].P;
	double T = myarray2.f1[i][j].P;
	l2sumP =  l2sumP + pow((P-T),2);
	
	P = myarray.f1[i][j].u;
	T = myarray2.f1[i][j].u;
	l2sumu =  l2sumu + pow((P-T),2);
	
	P = myarray.f1[i][j].v;
	T = myarray2.f1[i][j].v;
	l2sumv =  l2sumv + pow((P-T),2);

	}

}

double l2P = sqrt(l2sumP/(sx*sy));
double l2u = sqrt(l2sumu/(sx*sy));
double l2v = sqrt(l2sumv/(sx*sy));

cout << setprecision(8) << fixed << "L2 norm (P|u|v): " << l2P << " | " << l2u<< " | "  <<  l2v << "\n";
mydata.l2normP.push_back(l2P);
mydata.l2normu.push_back(l2u);
mydata.l2normv.push_back(l2v);
}


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
cout << "------------------------------------------------------------------------------------------------\n";

cout << "Solution Pressure:\n       |";
for(int i = 0; i < myarray.sizex; ++i)
		{
		if(i < 10)
		cout << "   i:  " << i <<"|";
		if(i > 9)
		cout << "   i: " << i <<"|";
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
		cout << setprecision(6) << fixed << myarray.s1[i][j].P <<"|";
		if(myarray.s1[i][j].P < 0)
		cout << setprecision(5) << fixed << myarray.s1[i][j].P <<"|";
		}
	
	}
	
	
	cout << "\nSolution u:\n       |";
for(int i = 0; i < myarray.sizex; ++i)
		{
		if(i < 10)
		cout << "   i:  " << i <<"|";
		if(i > 9)
		cout << "   i: " << i <<"|";
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
		cout << setprecision(6) << fixed << myarray.s1[i][j].u <<"|";
		if(myarray.s1[i][j].P < 0)
		cout << setprecision(5) << fixed << myarray.s1[i][j].u <<"|";
		}
	
	}
	
	
	cout << "\nSolution v:\n       |";
for(int i = 0; i < myarray.sizex; ++i)
		{
		if(i < 10)
		cout << "   i:  " << i <<"|";
		if(i > 9)
		cout << "   i: " << i <<"|";
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
		cout << setprecision(6) << fixed << myarray.s1[i][j].v <<"|";
		if(myarray.s1[i][j].P < 0)
		cout << setprecision(5) << fixed << myarray.s1[i][j].v <<"|";
		}
	
	}
cout << "\n\n";
}

//--------------------------Print u array in terminal----------------------------//

void print_array_fP(carray & myarray, double tstep)
{
cout << "------------------------------------------------------------------------------------------------\n";

cout << "Flux Pressure:\n       |";
for(int i = 0; i < myarray.sizex; ++i)
		{
		if(i < 10)
		cout << "   i:  " << i <<"|";
		if(i > 9)
		cout << "   i: " << i <<"|";
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
		cout << setprecision(6) << fixed << myarray.f1[i][j].P <<"|";
		if(myarray.f1[i][j].P < 0)
		cout << setprecision(5) << fixed << myarray.f1[i][j].P <<"|";
		}
	
	}
	
cout << "\nFlux u:\n       |";
for(int i = 0; i < myarray.sizex; ++i)
		{
		if(i < 10)
		cout << "   i:  " << i <<"|";
		if(i > 9)
		cout << "   i: " << i <<"|";
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
		cout << setprecision(6) << fixed << myarray.f1[i][j].u <<"|";
		if(myarray.f1[i][j].P < 0)
		cout << setprecision(5) << fixed << myarray.f1[i][j].u <<"|";
		}
	
	}
	
cout << "\nFlux v:\n       |";
for(int i = 0; i < myarray.sizex; ++i)
		{
		if(i < 10)
		cout << "   i:  " << i <<"|";
		if(i > 9)
		cout << "   i: " << i <<"|";
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
		cout << setprecision(6) << fixed << myarray.f1[i][j].v <<"|";
		if(myarray.f1[i][j].P < 0)
		cout << setprecision(5) << fixed << myarray.f1[i][j].v <<"|";
		}
	
	}		
cout << "\n\n";
}


/* Code for solving block-tridiagonal matrix problems, using the Thomas
   algorithm.  The only subroutine in here that you'll -need- to call is
   SolveThomas, although things like Add3x3 or AddVec might be useful,
   too. */

/* LHS array is sized as [*][3][3][3].  The last two indices identify
   the element within a block; the third index is the row and the fourth
   is the column of the Jacobian matrix.  The second index tells which
   block it is: 0 is below the main diagonal, 1 is on the main diagonal,
   2 is above the main diagonal.  The first index tells which block row
   you're looking at (the i or j index from the discretization). */

/* RHS array is [*][3].  The second index tells which element of the
   solution vector, and the first is the block row. */

/* Before linking this with your own code, you'll want to remove the
   main program included here as a test. */
#include "numerical.hpp"

static void SpewMatrix(double Source[3][3])
{
  printf("%10.6f %10.6f %10.6f\n", Source[0][0], Source[1][0], Source[2][0]);
  printf("%10.6f %10.6f %10.6f\n", Source[0][1], Source[1][1], Source[2][1]);
  printf("%10.6f %10.6f %10.6f\n", Source[0][2], Source[1][2], Source[2][2]);
}

static void SpewVector(double Source[3])
{
  printf("%10.6f %10.6f %10.6f\n", Source[0], Source[1], Source[2]);
}

static inline void CopyVec(const double Source[3],
			   double Target[3])
{
  Target[0] = Source[0];
  Target[1] = Source[1];
  Target[2] = Source[2];
}

static inline void Copy3x3(double Source[3][3],
			   double Target[3][3])
{
  Target[0][0] = Source[0][0];
  Target[0][1] = Source[0][1];
  Target[0][2] = Source[0][2];

  Target[1][0] = Source[1][0];
  Target[1][1] = Source[1][1];
  Target[1][2] = Source[1][2];

  Target[2][0] = Source[2][0];
  Target[2][1] = Source[2][1];
  Target[2][2] = Source[2][2];
}

static inline void Mult3x3(double A[3][3],
			   double B[3][3],
			   double C[3][3])
{
  C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
  C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1]; 
  C[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2]; 

  C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0]; 
  C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1]; 
  C[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2]; 

  C[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0]; 
  C[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1]; 
  C[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2]; 
}

static inline void MultVec(double A[3][3],
			   const double Vec[3],
			   double Result[3])
{
  Result[0] = A[0][0]*Vec[0] + A[0][1]*Vec[1] + A[0][2]*Vec[2]; 
  Result[1] = A[1][0]*Vec[0] + A[1][1]*Vec[1] + A[1][2]*Vec[2]; 
  Result[2] = A[2][0]*Vec[0] + A[2][1]*Vec[1] + A[2][2]*Vec[2]; 
}

static inline void Add3x3(double A[3][3],
			  double B[3][3],
			  const double Factor,
			  double C[3][3])
{
  C[0][0] = A[0][0] + Factor * B[0][0];
  C[0][1] = A[0][1] + Factor * B[0][1];
  C[0][2] = A[0][2] + Factor * B[0][2];

  C[1][0] = A[1][0] + Factor * B[1][0];
  C[1][1] = A[1][1] + Factor * B[1][1];
  C[1][2] = A[1][2] + Factor * B[1][2];

  C[2][0] = A[2][0] + Factor * B[2][0];
  C[2][1] = A[2][1] + Factor * B[2][1];
  C[2][2] = A[2][2] + Factor * B[2][2];
}

static inline void AddVec(const double A[3],
			  const double B[3],
			  const double Factor,
			  double C[3])
{
  C[0] = A[0] + Factor * B[0];
  C[1] = A[1] + Factor * B[1];
  C[2] = A[2] + Factor * B[2];
}

static inline void Invert3x3(double Block[3][3],
			     double Inverse[3][3])
{
  double DetInv = 1. / (+ Block[0][0]*Block[1][1]*Block[2][2]
			+ Block[0][1]*Block[1][2]*Block[2][0]
			+ Block[0][2]*Block[1][0]*Block[2][1]
			- Block[0][2]*Block[1][1]*Block[2][0]
			- Block[0][1]*Block[1][0]*Block[2][2]
			- Block[0][0]*Block[1][2]*Block[2][1]);

  /* Expand by minors to compute the inverse */
  Inverse[0][0] = + DetInv * (Block[1][1]*Block[2][2] -
			      Block[2][1]*Block[1][2]); 
  Inverse[1][0] = - DetInv * (Block[1][0]*Block[2][2] -
			      Block[2][0]*Block[1][2]); 
  Inverse[2][0] = + DetInv * (Block[1][0]*Block[2][1] -
			      Block[2][0]*Block[1][1]); 
  Inverse[0][1] = - DetInv * (Block[0][1]*Block[2][2] -
			      Block[2][1]*Block[0][2]); 
  Inverse[1][1] = + DetInv * (Block[0][0]*Block[2][2] -
			      Block[2][0]*Block[0][2]); 
  Inverse[2][1] = - DetInv * (Block[0][0]*Block[2][1] -
			      Block[2][0]*Block[0][1]); 
  Inverse[0][2] = + DetInv * (Block[0][1]*Block[1][2] -
			      Block[1][1]*Block[0][2]); 
  Inverse[1][2] = - DetInv * (Block[0][0]*Block[1][2] -
			      Block[1][0]*Block[0][2]); 
  Inverse[2][2] = + DetInv * (Block[0][0]*Block[1][1] -
			      Block[1][0]*Block[0][1]); 
}

void SolveBlockTri(double LHS[MAXSIZE][3][3][3],
		   double RHS[MAXSIZE][3],
		   int iNRows)
{
  int j;
  double Inv[3][3];

  for (j = 0; j < iNRows-1; j++) {
    /* Compute the inverse of the main block diagonal. */
    Invert3x3(LHS[j][1], Inv);
    /* Scale the right-most block diagonal by the inverse. */
    {
      double Temp[3][3];
      Mult3x3(Inv, LHS[j][2], Temp);
      Copy3x3(Temp, LHS[j][2]);
    }

    /* Scale the right-hand side by the inverse. */
    {
      double Temp[3];
      MultVec(Inv, RHS[j], Temp);
      CopyVec(Temp, RHS[j]);
    }      


    /* Left-multiply the jth row by the sub-diagonal on the j+1st row
       and subtract from the j+1st row.  This involves the
       super-diagonal term and the RHS of the jth row. */
    {
      /* First the LHS manipulation */
#define A LHS[j+1][0]
#define B LHS[j+1][1]
#define C LHS[j  ][2]
      double Temp[3][3], Temp2[3][3];
      double TVec[3], TVec2[3];
      Mult3x3(A, C, Temp);
      Add3x3(B, Temp, -1., Temp2);
      Copy3x3(Temp2, B);

      /* Now the RHS manipulation */
      MultVec(A, RHS[j], TVec);
      AddVec(RHS[j+1], TVec, -1., TVec2);
      CopyVec(TVec2, RHS[j+1]);
#undef A
#undef B
#undef C
    }
  } /* Done with forward elimination loop */
  /* Compute the inverse of the last main block diagonal. */
  j = iNRows-1;
  Invert3x3(LHS[j][1], Inv);

  /* Scale the right-hand side by the inverse. */
  {
    double Temp[3];
    MultVec(Inv, RHS[j], Temp);
    CopyVec(Temp, RHS[j]);
  }      
  
  /* Now do the back-substitution. */
  for (j = iNRows-2; j >= 0; j--) {
    /* Matrix-vector multiply and subtract. */
#define C LHS[j][2]
    RHS[j][0] -= (C[0][0]*RHS[j+1][0] +
		  C[0][1]*RHS[j+1][1] +
		  C[0][2]*RHS[j+1][2]);
    RHS[j][1] -= (C[1][0]*RHS[j+1][0] +
		  C[1][1]*RHS[j+1][1] +
		  C[1][2]*RHS[j+1][2]);
    RHS[j][2] -= (C[2][0]*RHS[j+1][0] +
		  C[2][1]*RHS[j+1][1] +
		  C[2][2]*RHS[j+1][2]);
#undef C
  }
}


/* LHS array is sized as [*][3][3][3].  The last two indices identify
   the element within a block; the third index is the row and the fourth
   is the column of the Jacobian matrix.  The second index tells which
   block it is: 0 is below the main diagonal, 1 is on the main diagonal,
   2 is above the main diagonal.  The first index tells which block row
   you're looking at (the i or j index from the discretization). */

void solve_block_thomas(carray & myarray, crow & r1, int NRows, int j)
{
  double LHS[100][3][3][3], RHS[100][3];
  int i;
          //|row|col
  for (i = 0; i < NRows; i++) {
  //A
    LHS[i][0][0][0] = r1.LHS[i].Ax.rP.P;
    LHS[i][0][0][1] = r1.LHS[i].Ax.ru.P;
    LHS[i][0][0][2] = r1.LHS[i].Ax.rv.P;   
    LHS[i][0][1][0] = r1.LHS[i].Ax.rP.u;
    LHS[i][0][1][1] = r1.LHS[i].Ax.ru.u;
    LHS[i][0][1][2] = r1.LHS[i].Ax.rv.u;    
    LHS[i][0][2][0] = r1.LHS[i].Ax.rP.v;
    LHS[i][0][2][1] = r1.LHS[i].Ax.ru.v;
    LHS[i][0][2][2] = r1.LHS[i].Ax.rv.v;
 //B
    LHS[i][1][0][0] = r1.LHS[i].Bx.rP.P;
    LHS[i][1][0][1] = r1.LHS[i].Bx.ru.P;
    LHS[i][1][0][2] = r1.LHS[i].Bx.rv.P;   
    LHS[i][1][1][0] = r1.LHS[i].Bx.rP.u;
    LHS[i][1][1][1] = r1.LHS[i].Bx.ru.u;
    LHS[i][1][1][2] = r1.LHS[i].Bx.rv.u;   
    LHS[i][1][2][0] = r1.LHS[i].Bx.rP.v;
    LHS[i][1][2][1] = r1.LHS[i].Bx.ru.v;
    LHS[i][1][2][2] = r1.LHS[i].Bx.rv.v;
 //C
    LHS[i][2][0][0] = r1.LHS[i].Cx.rP.P;
    LHS[i][2][0][1] = r1.LHS[i].Cx.ru.P;
    LHS[i][2][0][2] = r1.LHS[i].Cx.rv.P;   
    LHS[i][2][1][0] = r1.LHS[i].Cx.rP.u;
    LHS[i][2][1][1] = r1.LHS[i].Cx.ru.u;
    LHS[i][2][1][2] = r1.LHS[i].Cx.rv.u;    
    LHS[i][2][2][0] = r1.LHS[i].Cx.rP.v;
    LHS[i][2][2][1] = r1.LHS[i].Cx.ru.v;
    LHS[i][2][2][2] = r1.LHS[i].Cx.rv.v;

    RHS[i][0] = r1.RHS[i].P;
    RHS[i][1] = r1.RHS[i].u;
    RHS[i][2] = r1.RHS[i].v;
  }


  SolveBlockTri(LHS, RHS, NRows);

  /*for (i = 0; i < NRows; i++) {
    printf("%d %15.10f %15.10f %15.10f\n",
	   i, RHS[i][0], RHS[i][1], RHS[i][2]);
  }*/
  
  for (i = 1; i < NRows-1; i++) {
    myarray.f1[i][j].P = RHS[i][0];
    myarray.f1[i][j].u = RHS[i][1];
    myarray.f1[i][j].v = RHS[i][2];
  }
}


void solve_block_thomas(carray & myarray, ccol & r1, int NRows, int i, double & mdiff)
{
  double LHS[100][3][3][3], RHS[100][3];
  double tomp = 0;
  int j;
          //|row|col
  for (j = 0; j < NRows; j++) {
  //A
    LHS[j][0][0][0] = r1.LHS[j].Ay.rP.P;
    LHS[j][0][0][1] = r1.LHS[j].Ay.ru.P;
    LHS[j][0][0][2] = r1.LHS[j].Ay.rv.P;   
    LHS[j][0][1][0] = r1.LHS[j].Ay.rP.u;
    LHS[j][0][1][1] = r1.LHS[j].Ay.ru.u;
    LHS[j][0][1][2] = r1.LHS[j].Ay.rv.u;    
    LHS[j][0][2][0] = r1.LHS[j].Ay.rP.v;
    LHS[j][0][2][1] = r1.LHS[j].Ay.ru.v;
    LHS[j][0][2][2] = r1.LHS[j].Ay.rv.v;
 //B
    LHS[j][1][0][0] = r1.LHS[j].By.rP.P;
    LHS[j][1][0][1] = r1.LHS[j].By.ru.P;
    LHS[j][1][0][2] = r1.LHS[j].By.rv.P;   
    LHS[j][1][1][0] = r1.LHS[j].By.rP.u;
    LHS[j][1][1][1] = r1.LHS[j].By.ru.u;
    LHS[j][1][1][2] = r1.LHS[j].By.rv.u;   
    LHS[j][1][2][0] = r1.LHS[j].By.rP.v;
    LHS[j][1][2][1] = r1.LHS[j].By.ru.v;
    LHS[j][1][2][2] = r1.LHS[j].By.rv.v;
 //C
    LHS[j][2][0][0] = r1.LHS[j].Cy.rP.P;
    LHS[j][2][0][1] = r1.LHS[j].Cy.ru.P;
    LHS[j][2][0][2] = r1.LHS[j].Cy.rv.P;   
    LHS[j][2][1][0] = r1.LHS[j].Cy.rP.u;
    LHS[j][2][1][1] = r1.LHS[j].Cy.ru.u;
    LHS[j][2][1][2] = r1.LHS[j].Cy.rv.u;    
    LHS[j][2][2][0] = r1.LHS[j].Cy.rP.v;
    LHS[j][2][2][1] = r1.LHS[j].Cy.ru.v;
    LHS[j][2][2][2] = r1.LHS[j].Cy.rv.v;

    RHS[j][0] = r1.RHS[j].P;
    RHS[j][1] = r1.RHS[j].u;
    RHS[j][2] = r1.RHS[j].v;
  }


  SolveBlockTri(LHS, RHS, NRows);

  /*for (i = 0; i < NRows; i++) {
    printf("%d %15.10f %15.10f %15.10f\n",
	   i, RHS[i][0], RHS[i][1], RHS[i][2]);
  }*/
  
    for (j = 1; j < NRows-1; j++) {
    tomp = (abs(RHS[j][0]) + abs(RHS[j][1]) + abs(RHS[j][2]))/3.0;
    if(tomp > mdiff)
    mdiff = tomp;
  }
  
  for (j = 1; j < NRows-1; j++) {
    myarray.f1[i][j].P = RHS[j][0];
    myarray.f1[i][j].u = RHS[j][1];
    myarray.f1[i][j].v = RHS[j][2];
  }
}




//----------------------implicit wall boundaries-------------------------//


void set_wall(LHScX & temp, int par)
{
if(par == 1){
temp.Ax.rP.P = 0.;
temp.Ax.rP.u = 0.;
temp.Ax.rP.v = 0.;
temp.Ax.ru.P = 0.;
temp.Ax.ru.u = 0.;
temp.Ax.ru.v = 0.;
temp.Ax.rv.P = 0.;
temp.Ax.rv.u = 0.;
temp.Ax.rv.v = 0.;

temp.Bx.rP.P = 1.;
temp.Bx.rP.u = 0.;
temp.Bx.rP.v = 0.;
temp.Bx.ru.P = 0.;
temp.Bx.ru.u = 1.;
temp.Bx.ru.v = 0.;
temp.Bx.rv.P = 0.;
temp.Bx.rv.u = 0.;
temp.Bx.rv.v = 1.;

temp.Cx.rP.P = -1.;
temp.Cx.rP.u = 0.;
temp.Cx.rP.v = 0.;
temp.Cx.ru.P = 0.;
temp.Cx.ru.u = 1.;
temp.Cx.ru.v = 0.;
temp.Cx.rv.P = 0.;
temp.Cx.rv.u = 0.;
temp.Cx.rv.v = 1.0;
}

if(par == 0){
temp.Ax.rP.P = -1.;
temp.Ax.rP.u = 0.;
temp.Ax.rP.v = 0.;
temp.Ax.ru.P = 0.;
temp.Ax.ru.u = 1.;
temp.Ax.ru.v = 0.;
temp.Ax.rv.P = 0.;
temp.Ax.rv.u = 0.;
temp.Ax.rv.v = 1.;

temp.Bx.rP.P = 1.;
temp.Bx.rP.u = 0.;
temp.Bx.rP.v = 0.;
temp.Bx.ru.P = 0.;
temp.Bx.ru.u = 1.;
temp.Bx.ru.v = 0.;
temp.Bx.rv.P = 0.;
temp.Bx.rv.u = 0.;
temp.Bx.rv.v = 1.;

temp.Cx.rP.P = 0.;
temp.Cx.rP.u = 0.;
temp.Cx.rP.v = 0.;
temp.Cx.ru.P = 0.;
temp.Cx.ru.u = 0.;
temp.Cx.ru.v = 0.;
temp.Cx.rv.P = 0.;
temp.Cx.rv.u = 0.;
temp.Cx.rv.v = 0.;
}




}


void set_wall(LHScY & temp, int par)
{
if(par == 1){
temp.Ay.rP.P = 0.;
temp.Ay.rP.u = 0.;
temp.Ay.rP.v = 0.;
temp.Ay.ru.P = 0.;
temp.Ay.ru.u = 0.;
temp.Ay.ru.v = 0.;
temp.Ay.rv.P = 0.;
temp.Ay.rv.u = 0.;
temp.Ay.rv.v = 0.;

temp.By.rP.P = 1.;
temp.By.rP.u = 0.;
temp.By.rP.v = 0.;
temp.By.ru.P = 0.;
temp.By.ru.u = 1.;
temp.By.ru.v = 0.;
temp.By.rv.P = 0.;
temp.By.rv.u = 0.;
temp.By.rv.v = 1.;

temp.Cy.rP.P = -1.;
temp.Cy.rP.u = 0.;
temp.Cy.rP.v = 0.;
temp.Cy.ru.P = 0.;
temp.Cy.ru.u = 1.;
temp.Cy.ru.v = 0.;
temp.Cy.rv.P = 0.;
temp.Cy.rv.u = 0.;
temp.Cy.rv.v = 1.0;
}

if(par == 0){
temp.Ay.rP.P = -1.;
temp.Ay.rP.u = 0.;
temp.Ay.rP.v = 0.;
temp.Ay.ru.P = 0.;
temp.Ay.ru.u = 1.;
temp.Ay.ru.v = 0.;
temp.Ay.rv.P = 0.;
temp.Ay.rv.u = 0.;
temp.Ay.rv.v = 1.;

temp.By.rP.P = 1.;
temp.By.rP.u = 0.;
temp.By.rP.v = 0.;
temp.By.ru.P = 0.;
temp.By.ru.u = 1.;
temp.By.ru.v = 0.;
temp.By.rv.P = 0.;
temp.By.rv.u = 0.;
temp.By.rv.v = 1.;

temp.Cy.rP.P = 0.;
temp.Cy.rP.u = 0.;
temp.Cy.rP.v = 0.;
temp.Cy.ru.P = 0.;
temp.Cy.ru.u = 0.;
temp.Cy.ru.v = 0.;
temp.Cy.rv.P = 0.;
temp.Cy.rv.u = 0.;
temp.Cy.rv.v = 0.;
}




}


