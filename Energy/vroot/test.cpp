
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


int main()
{
//Calculating error as described in paper "procedure for estimation and reporting of uncertainty due to discretization in CFD applications"//

printf("\nCalculating Error...\n");

double h1 =1/sqrt(18000.0);//= ray1.DIMx;
double h2 =1/sqrt(8000.0);//= ray2.DIMx;
double h3 =1/sqrt(4500.0);//= ray3.DIMx;


float sol1 = 6.063;
float sol2 = 5.972;
float sol3 = 5.863;


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

return 1;
}
