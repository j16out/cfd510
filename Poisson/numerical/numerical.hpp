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


#define maxx 42
#define maxy 42
#define DIM 1.00000000
#define PI 3.141592654



struct carray{
float mcell [maxx][maxy];
vector<float> l2norm;
vector<float> diff;
int sizex = maxx;
int sizey = maxy;
int iterations = 0;

};



void set_ghostcells(carray & myarray);

void set_zero(carray & myarray);

void print_mcell(carray & myarray);

float gs_iter_SOR(carray & myarray, float omega);

float calc_source(int i, int j);

float calc_newcell(float source, float Tip1_j, float Tim1_j, float Ti_jp1 ,float Ti_jm1);

void get_surcells(carray & myarray, float & Tip1_j, float & Tim1_j, float & Ti_jp1 ,float & Ti_jm1, int i, int j);
































#endif
