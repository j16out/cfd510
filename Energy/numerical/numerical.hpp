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
#define maxx 8000
#define maxy 3
#define PI 3.141592654
#define U0 0
#define T0 1
#define V0 0
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
double DIMx = 0;
double DIMy = 0;
//scheme
int scheme = 0;
double tstep = 0;
double ctime = 0;
};



struct cdata{
//data storage specific to array
vector<double> l2norm;
vector<double> l1norm;
vector<double> linfnorm;
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

double vi_j = 0; 
double vi_jm1 = 0;
double vi_jp1 = 0;
};


//--------------------Init Arrays-----------------------------------------//

//set array size

void set_zero(carray & myarray);//zero entire array

void print_array(carray & myarray);//print array in terminal

void set_array_size(carray & myarray, int x, int y, double DIMx, double DIMy, int scheme);

//-------------------Boundary and Intial Conditions------------------------//

void set_ghostcells(carray & myarray);//set ghost cells

void set_intial_cond(carray & myarray);


//--------------------Flux calculation------------------------------------//
void solve_array(carray & myarray, double tmax, double cfl);

void compute_Flux(carray & myarray);

void get_nsurcells(carray & myarray, int i, int j, surr & mysurr);

double calc_newcell(carray & myarray, surr & s1);



//-----------------------Error calc related functions---------------------------//

double get_l2norm(carray & myarray, carray myarray2);//get estimated vale for l2 norm between arrays
double get_linf_norm(carray & myarray, carray myarray2);
double get_l1norm(carray & myarray, carray myarray2);

void set_analytic(carray & myarray, carray & numarray);//set analytic solution to a mesh



#endif


















//
