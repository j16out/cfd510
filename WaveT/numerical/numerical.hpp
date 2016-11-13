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

#define BIG 1000000
#define maxx 160
#define maxy 160
#define PI 3.141592654



struct carray{
//arrays
float mcellSOL [maxx][maxy];
float mcellSOL2 [maxx][maxy];
float mcellFI [maxx][maxy];
float mcellFI2 [maxx][maxy];

//array attributes
int sizex = maxx;
int sizey = maxy;
int iterations = 0;
float DIM1 = 0;

//current time 
float tstep = 0;
float ctime = 0;
float tmax = 0; 

//data specific to array
vector<float> l2norm;
vector<float> diff;

//temporary cells to store
float Tim1_j=0.0;
float Tim2_j=0.0;
float Ti_jm1=0.0;
float Ti_jm2=0.0;
float Ti_j=0.0;



};


void set_array_size(carray & myarray, int x, int y, float DIM);//set array size

void set_ghostcells(carray & myarray);//set ghost cells

void set_intial_cond(carray & myarray);


void set_zero(carray & myarray);//zero entire array


void print_array(carray & myarray);//print array in terminal


void solve_arrayRK2(carray & myarray, float E0);//solve the array



//fI related functions

int get_FIarray(carray & myarray, int stage);//get all FI for array for specific stage

float calc_FI(carray & myarray);//calculate new cell value based on 2nd order scheme

void get_surcells(carray & myarray, int i, int j, int stage);//obtain values of surrounding cells stage defines were result stored

int get_RK2(carray & myarray, int stage);


//Error calc related functions

void get_discrete_Error(carray ray1, carray ray2, carray ray3, float DIM);//get error using 3 arrays based on ASME solution accuarcy handout

float get_l2norm(carray & myarray, carray myarray2);//get estimated vale for l2 norm between arrays

void set_analytic(carray & myarray);

float get_solution(carray & myarray);


#endif


















//
