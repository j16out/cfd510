#include "root.hpp"

using namespace std;


//--------------------------3D graphs----------------------------//

void draw_3Dgraph(carray myarray, carray myarray2)
{

float DIM2 = myarray.DIM1;
TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700);
TCanvas *c2 = new TCanvas("c2","The FillRandom example",200,50,900,700);
TCanvas *c4 = new TCanvas("c4","The FillRandom example",200,50,900,700);

c1->cd();
string titlefile = "Numerical Solution 10x10 mesh; x; y; z";
const char* c; 
c = titlefile.c_str();	

TGraph2D *gr = new TGraph2D();
TGraph2D *gr1 = new TGraph2D();
TGraph2D *gr5 = new TGraph2D();
gr->SetTitle(c);

titlefile = "Solution Error 10x10 mesh; x; y; z";
c = titlefile.c_str();
gr1->SetTitle(c);	

int N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIM2*i;
	  float dy = DIM2*j;
	      float T = (cos(PI*dx)*sinh(PI*dy))/sinh(PI);
		float nT = myarray.mcell[i][j];		
	      gr1->SetPoint(N,dx,dy,abs(nT-T));
	      ++N;
	 }     
}

N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 0; j < myarray.sizey; j++) 
	{ float dx = DIM2*i;
	  float dy = DIM2*j;
	  float T = myarray.mcell[i][j];

	      gr->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}


N = 0;
for (int i = 1; i < myarray2.sizex-1; i++) 
{
	
	for (int j = 0; j < myarray2.sizey; j++) 
	{ float dx = DIM2*i;
	  float dy = DIM2*j;
	  float T = myarray2.mcell[i][j];

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
   gr5->Draw("surf1");   
   
}


//--------------------------error graphs----------------------------//

void draw_graph(carray myarray, carray myarray2)
{TCanvas *c6 = new TCanvas("c6","The FillRandom example",200,50,900,700);  
 TCanvas *c7 = new TCanvas("c7","The FillRandom example",200,50,900,700); 
string titlefile;
const char* c;  
float DIM2 = myarray.DIM1;

c6->cd();
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

 c6->SetLogy();	
gr1->Draw("ALP");
gr2->Draw("sameLP");
	 

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);

leg1->AddEntry(gr1,"Gauss-Seidel","AP");
leg1->AddEntry(gr2,"Gauss-Seidel SOR 1.5","AP");

//leg->AddEntry(fitb,"this one","l");
leg1->Draw();

c7->cd();
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


//--------------------------error graphs----------------------------//

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

gr1->Draw("AcP");
gr2->Draw("samecP");
gr3->Draw("samecP");
	 

TLegend *leg1 = new TLegend(0.75,0.9,0.9,0.8);

leg1->AddEntry(gr1,"SOR w=1","AP");
leg1->AddEntry(gr2,"SOR w=1.3","AP");
leg1->AddEntry(gr3,"SOR w=1.5","AP");

//leg->AddEntry(fitb,"this one","l");
leg1->Draw();



}


void draw_graph_l2norm3(carray myarray, carray myarray2, carray myarray3)
{TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,50,900,700); 
//c1->SetLeftMargin(0.15);
  
string titlefile;
const char* c;  
float DIM2 = myarray.DIM1;

c1->cd();
titlefile = "Accuracy Behavior for a 40 x 40 mesh; Iterations (N); L2 norm (10^-2)";
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






