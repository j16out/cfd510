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


//--------------------------set ghost cells----------------------------//

void set_ghostcells(carray & myarray)//boundary conditions
{
float DIM1 = myarray.DIM1;


for(int i = 1; i < myarray.sizex-1; ++i)
	{float d = i*DIM1;	
	myarray.mcell[i][0] = myarray.mcell[i][1];// y = 0 ie top
	myarray.mcell[0][i] = myarray.mcell[1][i];// x = 0 ie left
	myarray.mcell[i][myarray.sizex-1] = 5+((1/2)*pow((1+pow(d, 2)),3));//y = 1 ie bottom
	myarray.mcell[myarray.sizey-1][i] = 5+((1/2)*pow((1+pow(d, 2)),3));//x = 1 right
        }
}



//--------------------------Guass-Seidel SOR----------------------------//

float gs_iter_SOR(carray & myarray, float omega)
{
float DIM1 = myarray.DIM1;
float maxdiff = -1;
float Tip1_j, Tim1_j, Ti_jp1, Ti_jm1;
float l2sum = 0;
float sx = maxx-2;
float sy = maxy-2;

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
//---------------------Get source term for poisson problem----------------------//
float calc_source(carray & myarray, int i, int j)
{
float DIM1 = myarray.DIM1;
float dx = DIM1*i;
float dy = DIM1*j;
float source = pow(3*pow(dx,2)-3*pow(dy,2),2)+72*(pow(dx,2)*pow(dy,2))+pow(3*pow(dy,2)-3*pow(dx,2),2);

return source;
}

//--------------------------Get average solution at point------------------------//

float get_solution(carray & myarray)
{
int sx = (myarray.sizex-2)/2;
int sy = (myarray.sizey-2)/2;
float sol = (myarray.mcell[sx-1][sy+1]+myarray.mcell[sx+1][sy+1]+myarray.mcell[sx+1][sy-1]+myarray.mcell[sx+1][sy-1])/4;



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

