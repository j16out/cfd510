#include "root.hpp"

using namespace std;



//-----------------------------------part 3---------------------------------------//

void draw_graph_wave1_p3(carray & myarray1, carray myarray2, carray myarray3)
{
TCanvas *c16 = new TCanvas("c16","The FillRandom example",200,50,900,700);   
string titlefile;
const char* c;  


c16->cd();
titlefile = "Solutions; Mesh Size; T Error";
c = titlefile.c_str();
TGraph *gr1 = new TGraph();	
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(24);
	gr1->SetLineColor(1);
	gr1->SetTitle(c);  
TGraph *gr2 = new TGraph();	
	gr2->SetMarkerColor(3);
	gr2->SetMarkerStyle(25);
	gr2->SetLineColor(1);
	gr2->SetTitle(c);
TGraph *gr3 = new TGraph();	
	gr3->SetMarkerColor(2);
	gr3->SetMarkerStyle(26);
	gr3->SetLineColor(1);
	gr3->SetTitle(c); 
	
	TGraph *gr11 = new TGraph();	
	gr11->SetMarkerColor(4);
	gr11->SetMarkerStyle(20);
	gr11->SetLineColor(1);
	gr11->SetTitle(c);  
TGraph *gr21 = new TGraph();	
	gr21->SetMarkerColor(3);
	gr21->SetMarkerStyle(21);
	gr21->SetLineColor(1);
	gr21->SetTitle(c);
TGraph *gr31 = new TGraph();	
	gr31->SetMarkerColor(2);
	gr31->SetMarkerStyle(22);
	gr31->SetLineColor(1);
	gr31->SetTitle(c); 
	
	TGraph *gr = new TGraph();
	gr->SetTitle(c); 
	
int p = 50;
int size = myarray1.l1norm.size();
int n = 20;
printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
float temp = myarray1.l1norm.at(i);
gr1->SetPoint(i,n,temp);
n = n+p;
//printf("temp: %f i: %d\n", temp, n);
}



size = myarray2.l1norm.size();
n = 20;
//printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
float temp = myarray2.l1norm.at(i);
gr2->SetPoint(i,n,temp);
n = n+p;
//printf("temp: %f i: %d\n", temp, n);
}



size = myarray3.l1norm.size();
n = 20;
//printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
float temp = myarray3.l1norm.at(i);
gr3->SetPoint(i,n,temp);
n = n+p;
//printf("temp: %f i: %d\n", temp, n);
}



size = myarray1.linfnorm.size();
n = 20;
//printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
float temp = myarray1.linfnorm.at(i);
gr11->SetPoint(i,n,temp);
n = n+p;
//printf("temp: %f i: %d\n", temp, n);
}



size = myarray2.linfnorm.size();
n = 20;
//printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
float temp = myarray2.linfnorm.at(i);
gr21->SetPoint(i,n,temp);
n = n+p;
//printf("temp: %f i: %d\n", temp, n);
}


size = myarray3.linfnorm.size();
n = 20;
//printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
float temp = myarray3.linfnorm.at(i);
gr31->SetPoint(i,n,temp);
n = n+p;
//printf("temp: %f i: %d\n", temp, n);
}

 c16->SetLogy();
gr->SetPoint(0,0,10.0);

gr->SetPoint(1,3000,0.00001);
 
	

gr->Draw("AP");
gr21->Draw("samePl");
gr1->Draw("samePl");
gr3->Draw("samePl");

gr11->Draw("samePl");
gr2->Draw("samePl");
gr31->Draw("samePl");

	 

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);

leg1->AddEntry(gr1,"L1 2nd Upwind","APl");
leg1->AddEntry(gr2,"L1 1st Upwind","APl");
leg1->AddEntry(gr3,"L1 2nd Centered","APl");
leg1->AddEntry(gr11,"L Infinity 2nd Upwind","APl");
leg1->AddEntry(gr21,"L Infinity 1st Upwind","APl");
leg1->AddEntry(gr31,"L Infinity 2nd Centered","APl");


//leg->AddEntry(fitb,"this one","l");
leg1->Draw();
}





//-----------------------------------------prt 2----------------------------------------------------//


void draw_graph_wave1_p2(carray & myarray1, carray myarray2, carray myarray3)
{
TCanvas *c14 = new TCanvas("c14","The FillRandom example",200,50,900,700);   
string titlefile;
const char* c;  


c14->cd();
titlefile = "Solutions; x; T";
c = titlefile.c_str();
TGraph *gr1 = new TGraph();	
	gr1->SetMarkerColor(1);
	//gr1->SetMarkerStyle(24);
	gr1->SetLineColor(4);
	gr1->SetTitle(c);  
TGraph *gr2 = new TGraph();	
	gr2->SetMarkerColor(1);
	//gr2->SetMarkerStyle(25);
	gr2->SetLineColor(1);
	gr2->SetTitle(c);
TGraph *gr3 = new TGraph();	
	gr3->SetMarkerColor(1);
	//gr3->SetMarkerStyle(26);
	gr3->SetLineColor(2);
	gr3->SetTitle(c); 
	
	


for (int i = 2; i < myarray1.sizex; i++) 
{
float DIM2 = myarray1.DIM1; 
float dx = DIM2*(i-1.5);
float T = myarray1.mcellSOL[i][1];
//printf("T: %f dx: %f\n", T, dx);

      gr1->SetPoint(i-2,dx,T);
}  

for (int i = 2; i < myarray2.sizex; i++) 
{
float DIM2 = myarray2.DIM1; 
float dx = DIM2*(i-1.5);
float T = myarray2.mcellSOL[i][1];
//printf("T: %f dx: %f\n", T, dx);

      gr2->SetPoint(i-2,dx,T);
}  

for (int i = 2; i < myarray3.sizex; i++) 
{
float DIM2 = myarray3.DIM1; 
float dx = DIM2*(i-1.5);
float T = myarray3.mcellSOL[i][1];
//printf("T: %f dx: %f\n", T, dx);

      gr3->SetPoint(i-2,dx,T);
}  


 
 //c1->SetLogy();	

gr3->Draw("APc");
gr2->Draw("samePc");
gr1->Draw("samePc");

	 

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);

leg1->AddEntry(gr1,"CFL 0.4000 Very Stable ","Pl");
leg1->AddEntry(gr2,"CFL 0.5000 Borderline Stable","Pl");
leg1->AddEntry(gr3,"CFL 0.5015 Slightly Unstable ","Pl");


//leg->AddEntry(fitb,"this one","l");
leg1->Draw();
//---------------sub section for order----------//
TCanvas *c5 = new TCanvas("c5","The FillRandom example",200,50,900,700); 
c5->cd();
titlefile = "Order Evaluation; CFL Number ; L2 Error";
c = titlefile.c_str();


TGraph *gr5 = new TGraph();	
	gr5->SetMarkerColor(4);
	gr5->SetMarkerStyle(23);
	gr5->SetLineColor(1);
	gr5->SetTitle(c);  

c5->SetLogy();	
int size = myarray1.l2norm.size();
int n = 10;
printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
float temp = myarray1.l2norm.at(i);
float cfl = myarray1.diff.at(i);
gr5->SetPoint(i,cfl,temp);

printf("temp: %f i: %f\n", temp, cfl);
}


gr5->Draw("APl");


TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);
leg2->AddEntry(gr5,"L2 Error","AP");
leg2->Draw();
}


//--------------------------wave1 graphs----------------------------//

void draw_graph_wave1(carray & myarray1, carray myarray2, carray myarray3)
{
TCanvas *c17 = new TCanvas("c17","The FillRandom example",200,50,900,700);   
string titlefile;
const char* c;  


c17->cd();
titlefile = "Solutions; x; T";
c = titlefile.c_str();
TGraph *gr1 = new TGraph();	
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(24);
	gr1->SetLineColor(1);
	gr1->SetTitle(c);  
TGraph *gr2 = new TGraph();	
	gr2->SetMarkerColor(3);
	gr2->SetMarkerStyle(25);
	gr2->SetLineColor(1);
	gr2->SetTitle(c);
TGraph *gr3 = new TGraph();	
	gr3->SetMarkerColor(2);
	gr3->SetMarkerStyle(26);
	gr3->SetLineColor(1);
	gr3->SetTitle(c); 
TGraph *gr4 = new TGraph();	
//	gr4->SetMarkerColor(2);
//	gr4->SetMarkerStyle(7);
//	gr4->SetLineColor(1);
	gr4->SetTitle(c); 		
	


for (int i = 1; i < 500; i++) 
{ float dx = (1.0/500.0)*(i-0.5);
  float T = sin(2*PI*(2*(0.5)-dx));

      gr4->SetPoint(i,dx,T);
}

for (int i = 2; i < myarray1.sizex; i++) 
{
float DIM2 = myarray1.DIM1; 
float dx = DIM2*(i-1.5);
float T = myarray1.mcellSOL[i][1];
//printf("T: %f dx: %f\n", T, dx);

      gr1->SetPoint(i-2,dx,T);
}  

for (int i = 2; i < myarray2.sizex; i++) 
{
float DIM2 = myarray2.DIM1; 
float dx = DIM2*(i-1.5);
float T = myarray2.mcellSOL[i][1];
//printf("T: %f dx: %f\n", T, dx);

      gr2->SetPoint(i-2,dx,T);
}  

for (int i = 2; i < myarray3.sizex; i++) 
{
float DIM2 = myarray3.DIM1; 
float dx = DIM2*(i-1.5);
float T = myarray3.mcellSOL[i][1];
//printf("T: %f dx: %f\n", T, dx);

      gr3->SetPoint(i-2,dx,T);
}  


 
 //c1->SetLogy();	

gr4->Draw("Ac");
gr2->Draw("sameP");
gr3->Draw("sameP");
gr1->Draw("sameP");
	 

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);
/*
leg1->AddEntry(gr1,"20x1 mesh","AP");
leg1->AddEntry(gr2,"40x1 mesh","AP");
leg1->AddEntry(gr3,"80x1 mesh","AP");
*/
leg1->AddEntry(gr1,"2nd Upwind","APl");
leg1->AddEntry(gr2,"1st Upwind","APl");
leg1->AddEntry(gr3,"2nd Centered","APl");
leg1->AddEntry(gr4,"Analytic Solution","l");


//leg->AddEntry(fitb,"this one","l");
leg1->Draw();

//---------------------------------------------------------------------------------------------------

TCanvas *c5 = new TCanvas("c5","The FillRandom example",200,50,900,700); 
c5->cd();
titlefile = "Order Evaluation; Log(dx); Log(L2)";
c = titlefile.c_str();


TGraph *gr5 = new TGraph();	
	gr5->SetMarkerColor(4);
	gr5->SetMarkerStyle(23);
	gr5->SetLineColor(1);
	gr5->SetTitle(c);  


int size = myarray1.l2norm.size();
int n = 10;
printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
float temp = myarray1.l2norm.at(i);

n = n*2;
gr5->SetPoint(i,log(1.0/n),log(temp));

printf("temp: %f change x: %d\n", temp, n);
}





TF1 *tfit2 = new TF1("tfit2", "pol1");
tfit2->SetLineColor(2);
tfit2->SetLineWidth(1);
gr5->Fit(tfit2, "", "", -10, 10);
gStyle->SetOptFit();
gr5->Draw("AP");
TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(gr5,"Data","AP");
leg2->AddEntry(tfit2,"Linear Fit","l");
leg2->Draw();

}



void draw_graph_q1(carray & myarray1, carray & myarray2, carray & myarray3, carray analytic1, carray analytic2, carray analytic3)
{TCanvas *c11 = new TCanvas("c11","The FillRandom example",200,50,900,700);   
string titlefile;
const char* c;  
float DIM2 = myarray2.DIM1;

c11->cd();
titlefile = "Convergence Behavior for a 20 x 20 mesh; X; Error       ";
c = titlefile.c_str();
TGraph *gr1 = new TGraph();	
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(20);
	gr1->SetLineColor(1);
	gr1->SetTitle(c);  
TGraph *gr2 = new TGraph();	
	gr2->SetMarkerColor(3);
	gr2->SetMarkerStyle(21);
	gr2->SetLineColor(1);
	gr2->SetTitle(c);
TGraph *gr3 = new TGraph();	
	gr3->SetMarkerColor(2);
	gr3->SetMarkerStyle(22);
	gr3->SetLineColor(1);
	gr3->SetTitle(c); 	
	
	
  

for (int i = 2; i < myarray1.sizex; i++) 
{
float DIM2 = myarray1.DIM1; 
float dx = DIM2*(i-1.5);
float T = myarray1.mcellSOL[i][1]-analytic1.mcellSOL[i][1];
//printf("T: %f dx: %f\n", T, dx);

      gr1->SetPoint(i-2,dx,T);
}  

for (int i = 2; i < myarray2.sizex; i++) 
{
float DIM2 = myarray2.DIM1; 
float dx = DIM2*(i-1.5);
float T = myarray2.mcellSOL[i][1]-analytic2.mcellSOL[i][1];
//printf("T: %f dx: %f\n", T, dx);

      gr2->SetPoint(i-2,dx,T);
}  

for (int i = 2; i < myarray3.sizex; i++) 
{
float DIM2 = myarray3.DIM1; 
float dx = DIM2*(i-1.5);
float T = myarray3.mcellSOL[i][1]-analytic3.mcellSOL[i][1];
//printf("T: %f dx: %f\n", T, dx);

      gr3->SetPoint(i-2,dx,T);
}  


 
 //c1->SetLogy();	

gr2->Draw("APl");
gr1->Draw("samePl");
gr3->Draw("samePl");
	 

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);

leg1->AddEntry(gr1,"2nd Upwind","APl");
leg1->AddEntry(gr2,"1st Upwind","APl");
leg1->AddEntry(gr3,"2nd Centered","APl");


//leg->AddEntry(fitb,"this one","l");
leg1->Draw();

TCanvas *c14 = new TCanvas("c14","The FillRandom example",200,50,900,700); 	
c14->cd();	

titlefile = "Convergence Behavior for a 20 x 20 mesh; Iterations (N); T(K)-T(K+1)";
c = titlefile.c_str();
TGraph *gr11 = new TGraph();	
	gr11->SetMarkerColor(4);
	gr11->SetMarkerStyle(20);
	gr11->SetLineColor(1);
	gr11->SetTitle(c);  
TGraph *gr21 = new TGraph();	
	gr21->SetMarkerColor(3);
	gr21->SetMarkerStyle(21);
	gr21->SetLineColor(1);
	gr21->SetTitle(c);
TGraph *gr31 = new TGraph();	
	gr31->SetMarkerColor(2);
	gr31->SetMarkerStyle(22);
	gr31->SetLineColor(1);
	gr31->SetTitle(c); 
  

for (int i = 2; i < myarray1.sizex; i++) 
{
float DIM2 = myarray1.DIM1; 
float dx = DIM2*(i-1.5);
float T = abs(myarray1.mcellSOL[i][1]-analytic1.mcellSOL[i][1]);
//printf("T: %f dx: %f\n", T, dx);

      gr11->SetPoint(i-2,dx,T);
}  

for (int i = 2; i < myarray2.sizex; i++) 
{
float DIM2 = myarray2.DIM1; 
float dx = DIM2*(i-1.5);
float T = abs(myarray2.mcellSOL[i][1]-analytic2.mcellSOL[i][1]);
//printf("T: %f dx: %f\n", T, dx);

      gr21->SetPoint(i-2,dx,T);
}  

for (int i = 2; i < myarray3.sizex; i++) 
{
float DIM2 = myarray3.DIM1; 
float dx = DIM2*(i-1.5);
float T = abs(myarray3.mcellSOL[i][1]-analytic3.mcellSOL[i][1]);
//printf("T: %f dx: %f\n", T, dx);

      gr31->SetPoint(i-2,dx,T);
}  


 
 //c1->SetLogy();	

gr11->Draw("APl");
gr21->Draw("samePl");
gr31->Draw("samePl");
	 

TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(gr11,"2nd Upwind","AP");
leg2->AddEntry(gr21,"1st Upwind","AP");
leg2->AddEntry(gr31,"2nd Centered","AP");


//leg->AddEntry(fitb,"this one","l");
leg2->Draw();






}



//--------------------------3D graphs----------------------------//

void draw_3DgraphP(carray myarray)
{
float DIM2 = myarray.DIM1;

TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);



string titlefile = "Poisson Numerical; x; y; z";
const char* c; 
c = titlefile.c_str();	

TGraph2D *gr = new TGraph2D();


gr->SetTitle(c);

titlefile = "Poisson Analytical; x; y; z";
c = titlefile.c_str();
	


int N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIM2*(i-0.5);
	  float dy = DIM2*(j-0.5);
	  float T = myarray.mcellSOL[i][j];

	      gr->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}


N = 0;

c1->cd();

   gStyle->SetPalette(1);
   //gr->Draw("colz");
   gr->Draw("surf1z"); 

   

   
}

//--------------------------3D graphs----------------------------//

void draw_3Dgraph(carray myarray, carray myarray2)
{
float DIM2 = myarray.DIM1;

TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);
TCanvas *c2 = new TCanvas("c2","The FillRandom example",200,50,900,700);
TCanvas *c4 = new TCanvas("c4","The FillRandom example",200,50,900,700);

c1->cd();
string titlefile = "Laplace Numerical; x; y; z";
const char* c; 
c = titlefile.c_str();	

TGraph2D *gr = new TGraph2D();
TGraph2D *gr1 = new TGraph2D();
TGraph2D *gr5 = new TGraph2D();
gr->SetTitle(c);

titlefile = "Laplace Analytical; x; y; z";
c = titlefile.c_str();
gr1->SetTitle(c);	

int N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIM2*(i-1);
	  float dy = DIM2*(j-1);
	      float T = (cos(PI*dx)*sinh(PI*dy))/sinh(PI);

	      gr1->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}

N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 0; j < myarray.sizey; j++) 
	{ float dx = DIM2*(i-1);
	  float dy = DIM2*(j-1);
	  float T = myarray.mcellSOL[i][j];

	      gr->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}


N = 0;
for (int i = 1; i < myarray2.sizex-1; i++) 
{
	
	for (int j = 0; j < myarray2.sizey; j++) 
	{ float dx = DIM2*(i-1);
	  float dy = DIM2*(j-1);
	  float T = myarray2.mcellSOL[i][j];

	      gr5->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}
c1->cd();

   gStyle->SetPalette(1);
   gr->Draw("colz");
c2->cd();   
   gStyle->SetPalette(1);
   gr1->Draw("colz");
   
c4->cd();   
   gStyle->SetPalette(1);
   gr5->Draw("colz");   
   
}


//--------------------------diff graphs----------------------------//

void draw_graph(carray myarray, carray myarray2)
{TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);  
 TCanvas *c2 = new TCanvas("c2","The FillRandom example",200,50,900,700); 
string titlefile;
const char* c;  
float DIM2 = myarray.DIM1;

c1->cd();
titlefile = "Convergence Behavior for a 10 x 10 mesh; Iterations (N); T(K)-T(K+1)";
c = titlefile.c_str();
TGraph *gr1 = new TGraph();	
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(7);
	gr1->SetTitle(c);  
TGraph *gr2 = new TGraph();	
	gr2->SetMarkerColor(3);
	gr2->SetMarkerStyle(7);
	gr2->SetTitle(c); 	
	
for (int i = 0; i < myarray.diff.size(); ++i)
	{	
	float tomp = myarray.diff.at(i);
	gr1->SetPoint(i,i,tomp);		
	}
for (int i = 0; i < myarray2.diff.size(); ++i)
	{	
	float temp = myarray2.diff.at(i);
	gr2->SetPoint(i,i,temp);		
	}

 c1->SetLogy();	
gr1->Draw("ALP");
gr2->Draw("sameLP");
	 

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);

leg1->AddEntry(gr1,"Gauss-Seidel","AP");
leg1->AddEntry(gr2,"Gauss-Seidel SOR 1.5","AP");

//leg->AddEntry(fitb,"this one","l");
leg1->Draw();

c2->cd();
titlefile = "L2 norm Converged Solution for 10 x 10 mesh; Iterations (N); L2";
c = titlefile.c_str();
TGraph *gr3 = new TGraph();	
	gr3->SetMarkerColor(4);
	gr3->SetMarkerStyle(7);
	gr3->SetTitle(c);  
TGraph *gr4 = new TGraph();	
	gr4->SetMarkerColor(3);
	gr4->SetMarkerStyle(7);
	gr4->SetTitle(c); 	
	
for (int i = 0; i < myarray.l2norm.size(); ++i)
	{	
	float tomp = myarray.l2norm.at(i);
	gr3->SetPoint(i,i,tomp);		
	}
for (int i = 0; i < myarray2.l2norm.size(); ++i)
	{	
	float temp = myarray2.l2norm.at(i);
	gr4->SetPoint(i,i,temp);		
	}


gr3->Draw("ALP");
gr4->Draw("sameLP");
	 

TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(gr3,"Gauss-Seidel","AP");
leg2->AddEntry(gr4,"Gauss-Seidel SOR 1.5","AP");

//leg->AddEntry(fitb,"this one","l");
leg2->Draw();


}


//--------------------------diff graphs----------------------------//

void draw_graph_diff3(carray myarray, carray myarray2, carray myarray3)
{TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);   
string titlefile;
const char* c;  
float DIM2 = myarray.DIM1;

c1->cd();
titlefile = "Convergence Behavior for a 20 x 20 mesh; Iterations (N); T(K)-T(K+1)";
c = titlefile.c_str();
TGraph *gr1 = new TGraph();	
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(7);
	gr1->SetLineColor(1);
	gr1->SetTitle(c);  
TGraph *gr2 = new TGraph();	
	gr2->SetMarkerColor(3);
	gr2->SetMarkerStyle(7);
	gr2->SetLineColor(1);
	gr2->SetTitle(c); 
TGraph *gr3 = new TGraph();	
	gr3->SetMarkerColor(2);
	gr3->SetMarkerStyle(7);
	gr3->SetLineColor(1);
	gr3->SetTitle(c); 		
	
for (int i = 0; i < myarray.diff.size(); ++i)
	{	
	float tomp = myarray.diff.at(i);
	gr1->SetPoint(i,i,tomp);		
	}
for (int i = 0; i < myarray2.diff.size(); ++i)
	{	
	float temp = myarray2.diff.at(i);
	gr2->SetPoint(i,i,temp);		
	}
for (int i = 0; i < myarray3.diff.size(); ++i)
	{	
	float temp = myarray3.diff.at(i);
	gr3->SetPoint(i,i,temp);		
	}
 
 c1->SetLogy();	

gr1->Draw("AP");
gr2->Draw("sameP");
gr3->Draw("sameP");
	 

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);

leg1->AddEntry(gr1,"SOR w=1","AP");
leg1->AddEntry(gr2,"SOR w=1.3","AP");
leg1->AddEntry(gr3,"SOR w=1.5","AP");

//leg->AddEntry(fitb,"this one","l");
leg1->Draw();



}

//--------------------------l2 graphs----------------------------//

void draw_graph_l2norm3(carray myarray, carray myarray2, carray myarray3)
{TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);   
string titlefile;
const char* c;  
float DIM2 = myarray.DIM1;

c1->cd();
titlefile = "Accuracy Behavior for a 40 x 40 mesh; Iterations (N); L2 norm                 ";
c = titlefile.c_str();
TGraph *gr1 = new TGraph();	
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(7);
	gr1->SetLineColor(1);
	gr1->SetTitle(c);  
TGraph *gr2 = new TGraph();	
	gr2->SetMarkerColor(3);
	gr2->SetMarkerStyle(7);
	gr2->SetLineColor(1);
	gr2->SetTitle(c); 
TGraph *gr3 = new TGraph();	
	gr3->SetMarkerColor(2);
	gr3->SetMarkerStyle(7);
	gr3->SetLineColor(1);
	gr3->SetTitle(c); 		
	
for (int i = 0; i < myarray.l2norm.size(); ++i)
	{	
	float tomp = myarray.l2norm.at(i);
	gr1->SetPoint(i,i,tomp);		
	}
for (int i = 0; i < myarray2.l2norm.size(); ++i)
	{	
	float temp = myarray2.l2norm.at(i);
	gr2->SetPoint(i,i,temp);		
	}
for (int i = 0; i < myarray3.l2norm.size(); ++i)
	{	
	float temp = myarray3.l2norm.at(i);
	gr3->SetPoint(i,i,temp);		
	}
 
 //c1->SetLogy();	

gr1->Draw("ACP");
gr2->Draw("sameCP");
gr3->Draw("sameCP");
	 

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);

leg1->AddEntry(gr1,"SOR w=1","AP");
leg1->AddEntry(gr2,"SOR w=1.3","AP");
leg1->AddEntry(gr3,"SOR w=1.5","AP");

//leg->AddEntry(fitb,"this one","l");
leg1->Draw();



}






