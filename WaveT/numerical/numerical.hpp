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



struct carray{
//arrays
float mcellSOL [maxx][maxy];//first stage and solution mesh
float mcellSOL2 [maxx][maxy];//second stage solution mesh
float mcellFI [maxx][maxy];//first stage flux
float mcellFI2 [maxx][maxy];//second stage flux

//array attributes
int sizex = maxx;
int sizey = maxy;
float DIM1 = 0;

//current time 
float tstep = 0;
float ctime = 0;

//data specific to array
vector<float> l2norm;
vector<float> diff;

//temporary cells to store
float Tim1_j=0.0;
float Tim2_j=0.0;
float Ti_j=0.0;

};


void set_array_size(carray & myarray, int x, int y, float DIM);//set array size

void set_ghostcells(carray & myarray);//set ghost cells

void set_intial_cond(carray & myarray);


void set_zero(carray & myarray);//zero entire array


void print_array(carray & myarray);//print array in terminal


void solve_arrayRK2(carray & myarray, float tmax, float cfl);//solve the array



//fI related functions

void get_FIarray(carray & myarray, int stage);//get all FI for array for specific stage

void get_FIarray_1stcell(carray & myarray, int stage);

float calc_2nd_UW(carray & myarray);//calculate new cell value based on 2nd order scheme

void get_surcells(carray & myarray, int i, int j, int stage);//obtain values of surrounding cells stage defines were result stored

void get_RK2(carray & myarray, int stage);

void mv_SOL2_to_SOL1(carray & myarray);


//Error calc related functions

void get_discrete_Error(carray ray1, carray ray2, carray ray3, float DIM);//get error using 3 arrays based on ASME solution accuarcy handout

float get_l2norm(carray & myarray, carray myarray2);//get estimated vale for l2 norm between arrays

void set_analytic(carray & myarray, carray & numarray);

float get_solution(carray & myarray);


#endif


















//
