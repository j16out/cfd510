#include "laplace.hpp"


float DIM1 = DIM/(maxx-2);

bool get_cells(carray & myarray, float & Tip1_j, float & Tim1_j, float & Ti_jp1 ,float & Ti_jm1, int i, int j);


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
 


for(int i = 0; i < myarray.sizex-1; ++i)
	{
	float d = i*DIM1;
	myarray.mcell[i][0] = 0;//top
	myarray.mcell[0][i] = myarray.mcell[1][i];//left
	myarray.mcell[i][myarray.sizex-1] = cos(PI*d);//bottom
	myarray.mcell[myarray.sizey-1][i] = myarray.mcell[myarray.sizey-2][i];//right
        }
}

//--------------------------Guass-Seidel----------------------------//

float gs_iter(carray & myarray)
{
float maxdiff = -1;
float Tip1_j, Tim1_j, Ti_jp1, Ti_jm1;
float l2sum = 0;
float sx = maxx-2;
float sy = maxy-2;

//------over then down------//

       for(int j = 1; j < myarray.sizey-1; ++j)
	{
	
	
		for(int i = 1; i < myarray.sizex-1; ++i)
		{
		
		bool fcon = get_cells(myarray, Tip1_j, Tim1_j, Ti_jp1 , Ti_jm1, i, j);
		
		if(fcon == false)
		{
		cout << "****ERROR****  couldn't find match";
		break;
		} 
		
		float newcell = ((Tip1_j+Tim1_j+Ti_jp1+Ti_jm1)/(4));
		
		float diff = abs(newcell - myarray.mcell[i][j]);
		if(diff > maxdiff)
		maxdiff = diff;
		
		myarray.mcell[i][j] = newcell;
		set_ghostcells(myarray);

		}
	
	}
	
//------down then over------//	
	
       for(int i = 1; i < myarray.sizey-1; ++i)
	{
	
	
		for(int j = 1; j < myarray.sizex-1; ++j)
		{
		
		bool fcon = get_cells(myarray, Tip1_j, Tim1_j, Ti_jp1 , Ti_jm1, i, j);
		
		if(fcon == false)
		{
		cout << "****ERROR****  couldn't find match";
		break;
		} 

		float newcell = ((Tip1_j+Tim1_j+Ti_jp1+Ti_jm1)/(4));
		
		float diff = abs(newcell - myarray.mcell[i][j]);
		if(diff > maxdiff)
		maxdiff = diff;
		
		myarray.mcell[i][j] = newcell;		
		set_ghostcells(myarray);
		
		//-----error compare----//
		float dx = DIM1*i;
		float dy = DIM1*j;
		float T = (cos(PI*dx)*sinh(PI*dy))/sinh(PI);

		l2sum =  l2sum + pow((newcell-T),2);
		 
		}
	
	}
float l2 = sqrt(l2sum/(sx*sy));
myarray.l2norm.push_back(l2);
myarray.diff.push_back(maxdiff);
++myarray.iterations; 

return maxdiff;
}


//--------------------------Guass-Seidel SOR----------------------------//

float gs_iter_SOR(carray & myarray, float omega)
{
float maxdiff = -1;
float Tip1_j, Tim1_j, Ti_jp1, Ti_jm1;
float l2sum = 0;
float sx = maxx-2;
float sy = maxy-2;

//-----over then down-----//

       for(int j = 1; j < myarray.sizey-1; ++j)
	{
	
	
		for(int i = 1; i < myarray.sizex-1; ++i)
		{
		
		bool fcon = get_cells(myarray, Tip1_j, Tim1_j, Ti_jp1 , Ti_jm1, i, j);
		
		if(fcon == false)
		{
		cout << "****ERROR****  couldn't find match";
		break;
		} 
		
		float newcell = ((Tip1_j+Tim1_j+Ti_jp1+Ti_jm1)/(4));
		
		float diff = abs(newcell - myarray.mcell[i][j]);
		float delta = newcell - myarray.mcell[i][j];
		
		if(diff > maxdiff)
		maxdiff = diff;
		
		float newcellSOR = myarray.mcell[i][j] + omega*(delta);
		
		myarray.mcell[i][j] = newcellSOR;
		set_ghostcells(myarray);

		}
	
	}
	
//------down then over-----//	
	
       for(int i = 1; i < myarray.sizey-1; ++i)
	{
	
	
		for(int j = 1; j < myarray.sizex-1; ++j)
		{
		
		bool fcon = get_cells(myarray, Tip1_j, Tim1_j, Ti_jp1 , Ti_jm1, i, j);
		
		if(fcon == false)
		{
		cout << "****ERROR****  couldn't find match";
		break;
		} 

		float newcell = ((Tip1_j+Tim1_j+Ti_jp1+Ti_jm1)/(4));
		
		float diff = abs(newcell - myarray.mcell[i][j]);
		float delta = newcell - myarray.mcell[i][j];
		
		if(diff > maxdiff)
		maxdiff = diff;
		
		float newcellSOR = myarray.mcell[i][j] + omega*(delta);
		
		myarray.mcell[i][j] = newcellSOR;
		set_ghostcells(myarray);
		
		//-------error compare-------//
		float dx = DIM1*i;
		float dy = DIM1*j;
		float T = (cos(PI*dx)*sinh(PI*dy))/sinh(PI);

		l2sum =  l2sum + pow((newcellSOR-T),2);
		 
		}
	
	}
float l2 = sqrt(l2sum/(sx*sy));
myarray.l2norm.push_back(l2);
myarray.diff.push_back(maxdiff);
++myarray.iterations; 

return maxdiff;
}


//--------------------------Print array in terminal----------------------------//

void print_mcell(carray & myarray)
{
for(int j = 0; j < myarray.sizey; ++j)
	{
	cout << "\n";
	
		for(int i = 0; i < myarray.sizex; ++i)
		{
		cout << myarray.mcell[i][j] <<",";
		

		}
	
	}

}


//--------------------------Get current cell values----------------------------//

bool get_cells(carray & myarray, float & Tip1_j, float & Tim1_j, float & Ti_jp1 ,float & Ti_jm1, int i, int j)
{
float fcon = false;
float sizex = myarray.sizex;
float sizey = myarray.sizey;

		Tip1_j = myarray.mcell[i+1][j];
		Tim1_j = myarray.mcell[i-1][j];
		Ti_jp1 = myarray.mcell[i][j+1];
		Ti_jm1 = myarray.mcell[i][j-1];
		fcon = true;


return fcon;
}

