#ifndef laplace_INCLUDED
#define laplace_INCLUDED


#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <math.h> 

using namespace std;


#define maxx 22
#define maxy 22
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




































#endif
