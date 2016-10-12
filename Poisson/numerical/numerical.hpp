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
#define maxx 162
#define maxy 162
#define PI 3.141592654



struct carray{
float mcell [maxx][maxy];
vector<float> l2norm;
vector<float> diff;
int sizex = maxx;
int sizey = maxy;
int iterations = 0;
float DIM1 = 0;

};


void set_array_size(carray & myarray, int x, int y, float DIM);

void set_ghostcells(carray & myarray);

void set_zero(carray & myarray);

void print_array(carray & myarray);

float gs_iter_SOR(carray & myarray, float omega);

float calc_source(carray & myarray, int i, int j);

float calc_newcell(carray & myarray, float source, float Tip1_j, float Tim1_j, float Ti_jp1 ,float Ti_jm1);

void get_surcells(carray & myarray, float & Tip1_j, float & Tim1_j, float & Ti_jp1 ,float & Ti_jm1, int i, int j);

void solve_arraySOR(carray & myarray, float E0, float w);

void get_discrete_Error(carray ray1, carray ray2, carray ray3, float DIM);

float get_solution(carray & myarray);



























#endif
