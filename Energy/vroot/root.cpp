#include "root.hpp"

using namespace std;






void find_maxvalues(carray & myarray, carray & myarray2, carray & myarray3)
{

TCanvas *c4 = new TCanvas("c4","The FillRandom example",200,50,900,700);


string titlefile = "Energy Explicit flow1; Channel Bottom Length; Temperature Gradient";
const char* c; 
c = titlefile.c_str();	
TGraph *gr6 = new TGraph();	
	gr6->SetMarkerColor(4);
	gr6->SetMarkerStyle(24);
	gr6->SetLineColor(1);
	gr6->SetTitle(c); 
	TGraph *gr7 = new TGraph();	
	gr7->SetMarkerColor(3);
	gr7->SetMarkerStyle(25);
	gr7->SetLineColor(1);
	gr7->SetTitle(c); 
	TGraph *gr8 = new TGraph();	
	gr8->SetMarkerColor(2);
	gr8->SetMarkerStyle(26);
	gr8->SetLineColor(1);
	gr8->SetTitle(c); 

double temp;
 double mdif = -1000.0;
 double q = 0;
    double DIMx = myarray.DIMx;
    double DIMy = myarray.DIMy;  
	for(int j = 1; j < 2; ++j)
	{
		for(int i = 0; i < myarray.sizex; ++i)
		{
		TGraph *gr5 = new TGraph();	
		 
		 
		 	gr5->SetPoint(0,DIMy*(4*j-0.5),myarray.T1[i][4*j]);
		 gr5->SetPoint(1,DIMy*(4*(j+1)-0.5),myarray.T1[i][4*(j+1)]);
		 gr5->SetPoint(2,DIMy*(4*(j+2)-0.5),myarray.T1[i][4*(j+2)]);
		 //gr5->SetPoint(3,DIMy*(4*(j+3)-0.5),myarray.T1[i][4*(j+3)]);

		
		TF1 *tfit2 = new TF1("tfit2", "pol1");
        gr5->Fit(tfit2, "BRQE", "", -10, 10);
        float x = tfit2->GetParameter(1);
         if(x > mdif)
        {
        mdif = x;
        q = DIMx*(i-0.5);
        }
        //printf("x: %f\n", x);
        gr6->SetPoint(i,DIMx*(i-0.5),x);   

		}
		
	}
	
	printf("max slope: %f at %f\n", mdif, q);
	
	DIMx = myarray2.DIMx;
    DIMy = myarray2.DIMy; 
    temp;
    mdif = -1000.0;
	
		for(int j = 1; j < 2; ++j)
	{
		for(int i = 0; i < myarray2.sizex; ++i)
		{
		TGraph *gr5 = new TGraph();	
		 gr5->SetPoint(0,DIMy*(2*j-0.5),myarray2.T1[i][2*j]);
		 gr5->SetPoint(1,DIMy*(2*(j+1)-0.5),myarray2.T1[i][2*(j+1)]);
		 gr5->SetPoint(2,DIMy*(2*(j+2)-0.5),myarray2.T1[i][2*(j+2)]);
		 //gr5->SetPoint(3,DIMy*(2*(j+3)-0.5),myarray2.T1[i][2*(j+3)]);

		
		TF1 *tfit2 = new TF1("tfit2", "pol1");
        gr5->Fit(tfit2, "BRQE", "", -10, 10);
        float x = tfit2->GetParameter(1);
         if(x > mdif)
        {
        mdif = x;
        q = DIMx*(i-0.5);
        }
        //printf("x: %f\n", x);
        gr7->SetPoint(i,DIMx*(i-0.5),x);   

		}
		
	}
	
	printf("max slope: %f at %f\n", mdif, q);
	
	
	DIMx = myarray3.DIMx;
    DIMy = myarray3.DIMy;
    temp;
    mdif = -1000.0; 
	
		for(int j = 1; j < 2; ++j)
	{
		for(int i = 0; i < myarray3.sizex; ++i)
		{
		TGraph *gr5 = new TGraph();	
         gr5->SetPoint(0,DIMy*(j-0.5),myarray3.T1[i][j]);
		 gr5->SetPoint(1,DIMy*(j-0.5+1),myarray3.T1[i][j+1]);
		 gr5->SetPoint(2,DIMy*(j-0.5+2),myarray3.T1[i][j+2]);
		// gr5->SetPoint(3,DIMy*(j-0.5+3),myarray3.T1[i][j+3]);

		
		TF1 *tfit2 = new TF1("tfit2", "pol1");
        gr5->Fit(tfit2, "BRQE", "", -10, 10);
        float x = tfit2->GetParameter(1);
        if(x > mdif)
        {
        mdif = x;
        q = DIMx*(i-0.5);
        }
        
        gr8->SetPoint(i,DIMx*(i-0.5),x);   

		}
		
	}
	printf("max slope: %f at %f\n", mdif, q);
	

c4->SetGridx();
c4->SetGridy();
c4->SetTickx(1);
c4->SetTicky(1);


gr6->Draw("AP");
gr7->Draw("sameP");
gr8->Draw("sameP");
TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(gr6,"40x1 channel","AP");
leg2->AddEntry(gr7,"30x1 channel","AP");
leg2->AddEntry(gr8,"10x1 channel","AP");
leg2->Draw();
}
//-----------------------------------part 3---------------------------------------//

void draw_3Dgraph(carray myarray, carray myarray2)
{
float DIMx = myarray.DIMx;
float DIMy = myarray.DIMy;


TCanvas *c2 = new TCanvas("c2","The FillRandom example",200,50,900,700);
TCanvas *c4 = new TCanvas("c4","The FillRandom example",200,50,900,700);


string titlefile = "Energy Explicit flow1; x; y; z";
const char* c; 
c = titlefile.c_str();	


TGraph2D *gr1 = new TGraph2D();
TGraph2D *gr5 = new TGraph2D();
gr1->SetTitle(c);

titlefile = "Energy Implicit flow2; x; y; z";
c = titlefile.c_str();
gr5->SetTitle(c);	


int N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < 5; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray.T1[i][j];

	      gr1->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}


N = 0;
for (int i = 1; i < myarray2.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray2.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray2.T1[i][j];

	      gr5->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}
/*
c2->cd();   
   gStyle->SetPalette(1);
   gr1->Draw("colz");
   
c4->cd();   
   gStyle->SetPalette(1);
   gr5->Draw("colz");   */
   
   c2->cd();   
   gStyle->SetPalette(1);
   gr1->Draw("colz");
   
c4->cd();   
   gStyle->SetPalette(1);
   gr5->Draw("surf1z");   
   
}


//---------------------------------------------------------------------------------------

void draw_order_l2(cdata & mydata)
{
string titlefile;
const char* c;  

TCanvas *c5 = new TCanvas("c5","The FillRandom example",200,50,900,700); 
c5->cd();
titlefile = "Order Evaluation; Log(dx); Log(L2)";
c = titlefile.c_str();


TGraph *gr5 = new TGraph();	
	gr5->SetMarkerColor(4);
	gr5->SetMarkerStyle(23);
	gr5->SetLineColor(1);
	gr5->SetTitle(c);  


int size = mydata.l2norm.size();
int n = 5;
printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
double temp = mydata.l2norm.at(i);


gr5->SetPoint(i,log(1.0/n),log(temp));

printf("l2: %f change x: %d\n", temp, n);
n = n+5;
}


c5->SetGridx();
c5->SetGridy();
c5->SetTickx(1);
c5->SetTicky(1);


TF1 *tfit2 = new TF1("tfit2", "pol1");
tfit2->SetLineColor(2);
tfit2->SetLineWidth(1);
gr5->Fit(tfit2, "", "", -10, 10);
gStyle->SetOptFit();
gr5->Draw("AP");
TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(gr5,"log(l2)","AP");
leg2->AddEntry(tfit2,"Linear Fit","l");
leg2->Draw();


}



void draw_order_eff(cdata & mydata)
{
string titlefile;
const char* c;  

TCanvas *c5 = new TCanvas("c5","The FillRandom example",200,50,900,700); 
c5->cd();
titlefile = "Order Evaluation; Array Size; Time (us)";
c = titlefile.c_str();


TGraph *gr5 = new TGraph();	
	gr5->SetMarkerColor(4);
	gr5->SetMarkerStyle(23);
	gr5->SetLineColor(1);
	gr5->SetTitle(c); 
	
TGraph *gr6 = new TGraph();	
	gr6->SetMarkerColor(2);
	gr6->SetMarkerStyle(21);
	gr6->SetLineColor(1);
	gr6->SetTitle(c);  	 


int size = mydata.time1.size();
int n = 5;
printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
double temp = mydata.time1.at(i);
gr5->SetPoint(i,n,temp);

printf("time: %f change x: %d\n", temp, n);
n = n+5;
}



size = mydata.time2.size();
n = 5;
printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
double temp = mydata.time2.at(i);
gr6->SetPoint(i,n,temp);

printf("time: %f change x: %d\n", temp, n);
n = n+5;
}


c5->SetGridx();
c5->SetGridy();
c5->SetTickx(1);
c5->SetTicky(1);


gr6->Draw("AP");
gr5->Draw("sameP");
TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(gr5,"Implicit Euler","AP");
leg2->AddEntry(gr6,"Explicit Euler","AP");
leg2->Draw();


}

