#ifndef numerical_INCLUDED
#define numerical_INCLUDED


#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <math.h> 
#include <iomanip>

using namespace std;

#define BIG 10000
#define maxx 120
#define maxy 120
#define PI 3.141592654

#define U0 3
#define T0 1
#define V0 1

#define RE 50
#define PR 0.7
#define EC 0.1



struct carray{
//arrays

double f1 [maxx][maxy];//second stage solution mesh
double T1 [maxx][maxy];//first stage and solution mesh
double v1 [maxx][maxy];//first stage and solution mesh
double u1 [maxx][maxy];
//array attributes
int sizex = maxx;
int sizey = maxy;
double DIMx = 0.0;
double DIMy = 0.0;
//scheme
int scheme = 0;
double ctime = 0;
};



struct cdata{
//data storage specific to array
vector<double> l2norm;
vector<double> l1norm;
vector<double> linfnorm;
};

struct crow{
double LHS [maxx][3];
double RHS [maxx];
};

struct ccol{
double LHS [maxy][3];
double RHS [maxy];
};

struct surr{
//surrounding cells
double Tim1_j = 0;
double Tip1_j = 0;
double Ti_j = 0; 
double Ti_jm1 = 0;
double Ti_jp1 = 0;

double uim1_j = 0;
double uip1_j = 0;
double ui_j = 0;
double ui_jm1 = 0;
double ui_jp1 = 0; 


double vim1_j = 0;
double vip1_j = 0;
double vi_j = 0; 
double vi_jm1 = 0;
double vi_jp1 = 0;
};


//--------------------Init Arrays-----------------------------------------//

//set array size

void set_zero(carray & myarray);//zero entire array

void print_array(carray & myarray);//print array in terminal

void print_arrayu(carray & myarray);

void print_row(crow & myrow, carray & myarray);

void print_col(ccol & mycol, carray & myarray);

void set_array_size(carray & myarray, int x, int y, double DIMx, double DIMy, int scheme);

//-------------------Boundary and Intial Conditions------------------------//

void set_ghostcells(carray & myarray);//set ghost cells

void set_intial_cond(carray & myarray);


//--------------------Solve Implicit Euler------------------------------------//

void solve_array_IE(carray & myarray, double tmax, double cfl);

void solve_LinSys(carray & myarray, double tstep);

void load_row(carray & myarray, crow & myrow, int j, double tstep);

void load_col(carray & myarray, ccol & mycol, int i, double tstep);

void solve_thomas(crow & r, int iSize);

void solve_thomas(ccol & r, int iSize);

//--------------------Solve Explicit Euler------------------------------------//

void solve_array_EE(carray & myarray, double tmax, double cfl);

void time_advance_EE(carray & myarray, double tstep);

//--------------------Flux calculation------------------------------------//

void compute_Flux(carray & myarray);

void get_nsurcells(carray & myarray, int i, int j, surr & mysurr);

double calc_newcell(carray & myarray, surr & s1);



//-----------------------Error calc related functions---------------------------//

double get_l2norm(carray & myarray, carray myarray2);//get estimated vale for l2 norm between arrays

void set_analytic(carray & myarray, carray & numarray);//set analytic solution to a mesh



#endif


















//
