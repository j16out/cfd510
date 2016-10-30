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
		myarray.mcell[i][j] = 0;

		}
	}
}


//--------------------------set ghost cells for Laplace----------------------------//

void set_ghostcells(carray & myarray)//boundary conditions
{
float DIM1 = myarray.DIM1;
 
for(int i = 1; i < myarray.sizex-1; ++i)
	{
	float d = i*DIM1;
	myarray.mcell[i][0] = 0;//top
	myarray.mcell[0][i] = myarray.mcell[1][i];//left
	myarray.mcell[i][myarray.sizex-1] = cos(PI*d);//bottom
	myarray.mcell[myarray.sizey-1][i] = myarray.mcell[myarray.sizey-2][i];//right
        }
}



//--------------------------Guass-Seidel SOR----------------------------//

float gs_iter_SOR(carray & myarray, float omega)
{
float DIM1 = myarray.DIM1;
float maxdiff = -1;
float Tip1_j, Tim1_j, Ti_jp1, Ti_jm1;
float l2sum = 0;
float sx = myarray.sizex-2;
float sy = myarray.sizey-2;

float dx = 0;
float dy = 0;



//-----iterate through all y for steps x -----//

       for(int j = 1; j < myarray.sizey-1; ++j)
	{
	
	
		for(int i = 1; i < myarray.sizex-1; ++i)
		{
		dx = DIM1*i;
		dy = DIM1*j;
		
		//----get surrounding cells and compute new cell-------//
		get_surcells(myarray, Tip1_j, Tim1_j, Ti_jp1 , Ti_jm1, i, j);		
		float source = calc_source(myarray, i, j);
		float newcell = calc_newcell(myarray, source, Tip1_j, Tim1_j, Ti_jp1 , Ti_jm1); 
		
		//----apply over-relaxation-----//
		float delta = newcell - myarray.mcell[i][j];
		float newcellSOR = myarray.mcell[i][j] + omega*(delta);
		
		//-----find difference between new and old cell-----//
		float diff = abs(newcell - myarray.mcell[i][j]);
		if(diff > maxdiff)
		maxdiff = diff;
		
		
		//-----update current cell----//
		myarray.mcell[i][j] = newcellSOR;
		set_ghostcells(myarray);

		}
	
	}
	
//-----iterate through all x for steps y -----//	
	
       for(int i = 1; i < myarray.sizey-1; ++i)
	{
	
	
		for(int j = 1; j < myarray.sizex-1; ++j)
		{
		dx = DIM1*i;
		dy = DIM1*j;
		
		//----get surrounding cells and compute new cell-------//
		get_surcells(myarray, Tip1_j, Tim1_j, Ti_jp1 , Ti_jm1, i, j);		
		float source = calc_source(myarray, i, j);
		float newcell = calc_newcell(myarray, source, Tip1_j, Tim1_j, Ti_jp1 , Ti_jm1); 
		
		//----apply over-relaxation-----//
		float delta = newcell - myarray.mcell[i][j];
		float newcellSOR = myarray.mcell[i][j] + omega*(delta);
		
		//-----find difference between new and old cell-----//
		float diff = abs(newcell - myarray.mcell[i][j]);
		if(diff > maxdiff)
		maxdiff = diff;
		
		
		//-----update current cell----//
		myarray.mcell[i][j] = newcellSOR;
		set_ghostcells(myarray);
		
		//-------error compare-------//
		
		float T = (cos(PI*dx)*sinh(PI*dy))/sinh(PI);
		l2sum =  l2sum + pow((newcellSOR-T),2);
		
		 
		}
	
	}
//---get l2norm value------//	
float l2 = sqrt(l2sum/(sx*sy));

//----store data-------//
myarray.l2norm.push_back(l2);
myarray.diff.push_back(maxdiff);
++myarray.iterations; 

return maxdiff;
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

	while(diff > E0)
	{
	diff = gs_iter_SOR(myarray, w);


		if(diff > BIG)//avoid infinite loops if diverges
		break;
	
		if(update >= update2)//report difference every 100 steps
		{cout << "Update: step " << update << " Solution Change: " << setprecision(9) << fixed << diff << " \n"; 
		 update2 = update2 + 100;
		}
		
		if(ldiff == diff)//checks for repeated values indication of instability for high w
		++div;
		else
		div = 0;
	
		if(div > 20 && relax_on && w > 1.3)//reduces over-relaxation for high w when unstable
		{
		w = 1.3;
		relax_on = false;
		cout << "Relaxation Reduced to 1.0 @ " << myarray.iterations << " \n";
		}
			
	ldiff = diff;
	++update;	
	}

cout << "Iterations: " << myarray.iterations << "   L2 norm: " << setprecision(6) << fixed << myarray.l2norm.at(myarray.iterations-1) <<  "\n";
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

//--------------------Calculate new cell value from neibhors ---------------------//


float calc_newcell(carray & myarray, float source, float Tip1_j, float Tim1_j, float Ti_jp1 ,float Ti_jm1)
{
float DIM1 = myarray.DIM1;
float chx = DIM1;
float chy = DIM1;
float temp = (pow(chx,2)*pow(chy,2)) /  (2*(pow(chx,2)+pow(chy,2)));
float newcell = ((  ((Tip1_j+Tim1_j)/pow(chx,2))  +  ((Ti_jp1+Ti_jm1)/pow(chy,2))  ) * temp) - source;

return newcell;
}
//---------------------Get source term for Laplace problem----------------------//
float calc_source(carray & myarray, int i, int j)
{
float DIM1 = myarray.DIM1;
float dx = DIM1*i;
float dy = DIM1*j;
float source = 0;

return source;
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









//-----------------------------BC checking-------------------------------------//

#define MAX(a, b) ((a) > (b) ? (a) : (b))

double BCfunc(const double coord) 
{
  return cos(M_PI*coord);
}

int approxEqual(const double a, const double b)// checks if equal within a relative tolerance, truncating only works if certain your values fall within a range, this avoids that
{
  double ave = (a + b) * 0.5;
  double tol = MAX(ave*1.e-10, 1.e-10);
  return (fabs(a-b) < tol ? 1 : 0);
}

int BCtest(const double coord, const double interiorValue,
	    const double ghostValue)
{
  double wallValue = BCfunc(coord);
  double aveValue = (interiorValue + ghostValue) * 0.5;
  return approxEqual(wallValue, aveValue);
}

#define IMAX 20
#define JMAX 20

int testAllBCs(double soln[IMAX+2][JMAX+2]) 
{
  double dx = 1./IMAX;
  int i;
  int OK = 1;
  for (i = 1; i <= IMAX && OK; i++) {
    double x = (i-0.5) * dx;
    OK = BCtest(x, soln[i][JMAX], soln[i][JMAX+1]);
  }

  if (i == IMAX+1 && OK == 1) {
    printf("BC test passed\n");
    return 1;
  }
  else {
    printf("BC test failed\n");
    return 0;
  }
}

void testmyBC(carray & myarray)
{

  double soln[IMAX+2][JMAX+2];
  double dx = 1./IMAX;
  double dy = 1./JMAX;
  int i, j;
  
  for (i = 1; i <= IMAX; i++) {
    for (j = 1; j <= JMAX; j++) {
      soln[i][j] = 0;
    }
  }

  for (i = 1; i <= IMAX; i++) {
    double x = (i-0.5) * dx;
    soln[i][JMAX+1] = 2*BCfunc(x) - soln[i][JMAX];
  }

  testAllBCs(soln);

}





