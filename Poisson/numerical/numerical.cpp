#include "numerical.hpp"





//----------set array size (working area excluding ghost)---------------//

void set_array_size(carray & myarray, int x, int y, float DIM)
{
	if(x < 160 && y < 160)
	{
	myarray.sizex = x+2;
	myarray.sizey = y+2;
	myarray.DIM1 = DIM/(x);
	}
	else
	cout << "Array size to big, setting to default 160" << "\n";

}


//--------------------------zero array----------------------------//

void set_zero(carray & myarray)
{
	for(int i = 0; i < myarray.sizex; ++i)
	{
		for(int j = 0; j < myarray.sizey; ++j)
		{
		myarray.mcell[i][j] = 0;//set everything to zero

		}
	}
}


//--------------------------set ghost cells for Poisson----------------------------//


void set_ghostcells(carray & myarray)
{
float DIM1 = myarray.DIM1;

//set boundary conditions in ghost cells
for(int i = 1; i < myarray.sizex-1; ++i)
	{float d = (i-0.5)*DIM1;
	
	//neumann boundaries	
	myarray.mcell[i][0] = myarray.mcell[i][1];//top ghost cells 
	myarray.mcell[0][i] = myarray.mcell[1][i];//left ghost cells
		
	//dirichlet boundaries	
	myarray.mcell[i][myarray.sizex-1] = 2.0*(5.0-((1.0/2.0)*pow((1.0+pow(d, 2)),3))) - myarray.mcell[i][myarray.sizex-2];
	myarray.mcell[myarray.sizey-1][i] = 2.0*(5.0-((1.0/2.0)*pow((1.0+pow(d, 2)),3))) - myarray.mcell[myarray.sizex-2][i];
	
        }
}



//--------------------------Guass-Seidel SOR----------------------------//

float gs_iter_SOR(carray & myarray, float omega)
{
float DIM1 = myarray.DIM1;//get dimensions of array (not grid size)
float Tip1_j, Tim1_j, Ti_jp1, Ti_jm1;//define values for surrounding cells
float l2sum = 0.0;
float sx = myarray.sizex-2;
float sy = myarray.sizey-2;

float dx = 0.0;//change in x
float dy = 0.0;//change in y

carray oldarray=myarray;


//-----iterate through all x for steps y -----//
set_ghostcells(myarray);

for(int j = 1; j < myarray.sizey-1; ++j)
{


	for(int i = 1; i < myarray.sizex-1; ++i)
	{
	
	dx = (i-0.5)*DIM1;
	dy = (j-0.5)*DIM1;
	
	//----get surrounding cells and compute new cell-------//
	get_surcells(myarray, Tip1_j, Tim1_j, Ti_jp1 , Ti_jm1, i, j);		
	float source = calc_source(myarray, i, j);
	float newcell = calc_newcell(myarray, source, Tip1_j, Tim1_j, Ti_jp1 , Ti_jm1); 
	
	//----apply over-relaxation-----//
	float delta = newcell - myarray.mcell[i][j];
	float newcellSOR = myarray.mcell[i][j] + omega*(delta);
	
	//-----update current cell----//
	myarray.mcell[i][j] = newcellSOR;
	

	}

}
	
set_ghostcells(myarray);	
//-----iterate through all y for steps x -----//		
for(int i = 1; i < myarray.sizey-1; ++i)
{


	for(int j = 1; j < myarray.sizex-1; ++j)
	{
	
	dx = (i-0.5)*DIM1;
	dy = (j-0.5)*DIM1;
	
	//----get surrounding cells and compute new cell-------//
	get_surcells(myarray, Tip1_j, Tim1_j, Ti_jp1 , Ti_jm1, i, j);		
	float source = calc_source(myarray, i, j);
	float newcell = calc_newcell(myarray, source, Tip1_j, Tim1_j, Ti_jp1 , Ti_jm1); 
	
	//----apply over-relaxation-----//
	float delta = newcell - myarray.mcell[i][j];
	float newcellSOR = myarray.mcell[i][j] + omega*(delta);

	
	//-----update current cell----//
	myarray.mcell[i][j] = newcellSOR;
	
	}

}


float maxdiff = -1.0;

for(int i = 2; i < myarray.sizey-2; ++i)
{	
	for(int j = 2; j < myarray.sizex-2; ++j)
	{
	float diff = abs(oldarray.mcell[i][j] - myarray.mcell[i][j]);
		if(diff > maxdiff)
		{
		maxdiff = diff;
		//cout << "coord " << maxdiff << " " << i << "  " << j << "\n";
		}
	}	

}

//for obtaining l2norm convergence

/*carray analytic;
set_array_size(analytic, 10, 10, 1.0);
set_analytic(analytic);
//float norm = get_l2norm(myarray, analytic);

myarray.l2norm.push_back(norm);	*/
myarray.diff.push_back(maxdiff);
++myarray.iterations; 

return maxdiff;
}
//-------------------------Get L2 nrom for unknown analytical----------------------//

float get_l2norm(carray & myarray, carray myarray2)
{
float l2sum =0;
float sx = myarray.sizex-2;
float sy = myarray.sizey-2;

for(int j = 1; j < myarray.sizey-1; ++j)
{	
	for(int i = 1; i < myarray.sizex-1; ++i)
	{

	float P = myarray.mcell[i][j];
	float T = myarray2.mcell[i][j];
	l2sum =  l2sum + pow((P-T),2);

	}

}

float l2 = sqrt(l2sum/(sx*sy));
cout << "L2 norm: " << l2 << "\n";
return l2;
}

//--------------------------Solve array using GS-iterations----------------------------//

void solve_arraySOR(carray & myarray, float E0, float w)
{
printf("\n\nSolving Grid size: %d Relaxation: %f\n", myarray.sizex, w);
bool relax_on = true; 
float diff = 1;// current difference
float ldiff = BIG;// previous difference
int div = 0;
int update = 0;
int update2 = 100;

while(diff >= E0)
{
diff = gs_iter_SOR(myarray, w);


	if(diff > BIG)//avoid infinite loops if diverges
	break;

	if(update >= update2)//report difference every 100 steps
	{cout << "Update: step " << update << " Solution Change: " << setprecision(9) << fixed << diff << " \n"; 
	//print_array(myarray);
	 update2 = update2 + 100;
	}
	
	if(ldiff == diff)//checks for repeated values indication of instability for high w
	++div;
	else
	div = 0;

	if(div > 3 && w > 1.1)//reduces over-relaxation for high w when unstable
	{
	w = 1.0;
	cout << "Relaxation Reduced to "<<w<<" @ " << myarray.iterations << " \n";
	}
		
ldiff = diff;
++update;	
}

cout << "Iterations: " << myarray.iterations <<  "\n";
}


//--------------------------Print array in terminal----------------------------//

void print_array(carray & myarray)
{
cout << "\n";

	for(int j = 0; j < myarray.sizey; ++j)
	{
	cout << "\n|";	
		for(int i = 0; i < myarray.sizex; ++i)
		{
		if(myarray.mcell[i][j] >= 0)
		cout << setprecision(3) << fixed << myarray.mcell[i][j] <<"|";
		if(myarray.mcell[i][j] < 0)
		cout << setprecision(2) << fixed << myarray.mcell[i][j] <<"|";
		}
	
	}
cout << "\n";
}

//--------------------Calculate new cell value from neighbors ---------------------//


float calc_newcell(carray & myarray, float source, float Tip1_j, float Tim1_j, float Ti_jp1 ,float Ti_jm1)
{
float DIM1 = myarray.DIM1;
float chx = DIM1;
float chy = DIM1;
float temp = (pow(chx,2)*pow(chy,2)) /  (2*(pow(chx,2)+pow(chy,2)));
float newcell = ((  ((Tip1_j+Tim1_j)/pow(chx,2))  +  ((Ti_jp1+Ti_jm1)/pow(chy,2)) - source ) * temp) ;

return newcell;
}
//---------------------Get source term for poisson problem----------------------//
float calc_source(carray & myarray, int i, int j)
{
float DIM1 = myarray.DIM1;
float dx = (i-0.5)*DIM1;
float dy = (j-0.5)*DIM1;
float source = -1.0*(pow(3*pow(dx,2)-3*pow(dy,2),2)+72.0*(pow(dx,2)*pow(dy,2))+pow(3*pow(dy,2)-3.0*pow(dx,2),2));
//source = 0 for Laplace problem

return source;
}

//-----------------------Get average solution at point (1/2)(1/2)--------------------//

float get_solution(carray & myarray)
{
float DIM1 = myarray.DIM1;
int sx = (myarray.sizex)/2.0;
int sy = (myarray.sizey)/2.0;
float sol = (myarray.mcell[sx-1][sy]+myarray.mcell[sx][sy]+myarray.mcell[sx][sy-1]+myarray.mcell[sx-1][sy-1])/4.0;

printf("cell 1: %f cell 2: %f cell 3: %f cell 4: %f\n",myarray.mcell[sx-1][sy],myarray.mcell[sx][sy],myarray.mcell[sx][sy-1],myarray.mcell[sx-1][sy-1]);
//for Poisson problem only, finds value based on average of four surrounding cells

return sol;
}

//--------------------------Get current cell values----------------------------//

void get_surcells(carray & myarray, float & Tip1_j, float & Tim1_j, float & Ti_jp1 ,float & Ti_jm1, int i, int j)
{
float fcon = false;
float sizex = myarray.sizex;
float sizey = myarray.sizey;

		Tip1_j = myarray.mcell[i+1][j];
		Tim1_j = myarray.mcell[i-1][j];
		Ti_jp1 = myarray.mcell[i][j+1];
		Ti_jm1 = myarray.mcell[i][j-1];

}

//--------------------------Get descrete error----------------------------//

void get_discrete_Error(carray ray1, carray ray2, carray ray3, float DIM)
{
//Calculating error as described in paper "procedure for estimation and reporting of uncertainty due to discretization in CFD applications"//

printf("\nCalculating Error...\n");

float h1 = DIM/ray1.sizex;
float h2 = DIM/ray2.sizex;
float h3 = DIM/ray3.sizex;


float sol1 = get_solution(ray1);
float sol2 = get_solution(ray2);
float sol3 = get_solution(ray3);



printf("h1: %f \nh2: %f \nh3: %f, \nsol1: %f \nsol2: %f \nsol3: %f\n",h1, h2, h3, sol1, sol2, sol3);

float r21 = h2/h1;
float r32 = h3/h2;

printf("\nr32: %f \nr21: %f\n",r32, r21);

float e32 = sol3-sol2;
float e21 = sol2-sol1;

float s = (e32/e21);
if(s >= 0)
s = 1;
else
s = -1;

float p_n = 0;
float p = (1/log(r21))*(abs(log(abs(e32/e21))+0));

printf("intial guess: %f \n", p);

float diff = 1;

	while(diff > 0.0000001)
	{

	float p_n = (1/log(r21))*(abs(log(abs(e32/e21))+log((pow(r21,p)-s)/(pow(r32,p)-s)) ));
	diff = abs(p_n -p);
	//printf("p_n: %f p: %f diff: %f\n",p_n, p, diff);

	p = p_n;
	}
 
//
float sol_ext21 = (pow(r21, p)*sol1-sol2)/(pow(r21,p)-1.0);
float sol_ext32 = (pow(r32, p)*sol2-sol3)/(pow(r32,p)-1.0);

printf("order: %f \nphi_ext21: %f \nphi_ext32 %f\n",p, sol_ext21, sol_ext32);

float ea21 = abs((sol1-sol2)/sol1);

float e_ext21 = abs((sol_ext21-sol1)/sol_ext21);

float GCI_21 = (1.25*ea21)/(pow(r21,p)-1.0);


printf("ea21: %f  \ne_ext21: %f  \nGC121 %f \n", ea21, e_ext21, GCI_21);

}













