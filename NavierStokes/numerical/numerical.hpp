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

#define P0 1.0
#define U0 1.0
#define V0 1.0
#define BETA 1.0

#define RE 10
#define PR 0.7
#define EC 0.001


struct vec{
double P = 0.0;
double v = 0.0;
double u = 0.0;
};

struct carray{

//flux and solution vector arrays
struct vec f1 [maxx][maxy];
struct vec s1 [maxx][maxy];
struct vec lhs [maxx][maxy];


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
vector<double> l2normP;
vector<double> l2normu;
vector<double> l2normv;
vector<double> l1norm;
vector<double> linfnorm;
vector<double> time1;
vector<double> time2;
};

struct crow{
//data for loaded row
double LHS [maxx][3];
double RHS [maxx];
};

struct ccol{
//data for loaded col
double LHS [maxy][3];
double RHS [maxy];
};

struct surr{
//surrounding cells
double Pim1_j = 0;
double Pip1_j = 0;
double Pi_j = 0; 
double Pi_jm1 = 0;
double Pi_jp1 = 0;

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

//------------------------LHS--------------------//
struct vec3x3{
vec rP;
vec ru;
vec rv;
};


struct LHSconst{
vec3x3 Ax;
vec3x3 Ay;
vec3x3 Bx;
vec3x3 By;
vec3x3 Cx;
vec3x3 Cy;

};
//--------------------Init Arrays-----------------------------------------//

void set_zero(carray & myarray);//zero entire array

void set_array_size(carray & myarray, int x, int y, double DIMx, double DIMy, int scheme);//set array size


//-------------------Boundary and Intial Conditions------------------------//


void set_init_cond(carray & myarray);


//--------------------Solve Implicit Euler------------------------------------//

void solve_array_IE(carray & myarray, double tmax, double cfl);


//------------------LHS calculation--------------------------------------//

void calc_LHS_const(carray & myarray, LHSconst & c1, int i, int j);

void solve_array_LHS(carray & myarray);


//--------------------Flux calculation------------------------------------//

void update_flux(carray & myarray);

void get_nsurcells(carray & myarray, int i, int j, surr & mysurr);

void calc_flux(carray & myarray, surr & s1, vec & ftemp);



//-----------------------Error calc related functions---------------------------//
/*

double get_l1normD(carray & myarray, carray myarray2);

void get_discrete_Error(carray ray1, carray ray2, carray ray3);
*/
void get_l2norm(carray & myarray, carray myarray2, cdata & mydata);

void set_analytic(carray & myarray);//set analytic solution to a mesh


//-----------------------Print related functions-------------------------------//

void print_array_sP(carray & myarray);//print array in terminal

void print_array_fP(carray & myarray);

/*
void print_row(crow & myrow, carray & myarray);

void print_col(ccol & mycol, carray & myarray);
*/

#endif


















//
