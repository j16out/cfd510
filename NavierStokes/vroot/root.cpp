#include "root.hpp"

using namespace std;






//-----------------------------------part 3---------------------------------------//

void draw_3Dgraph_s(carray myarray, carray myarray2)
{
float DIMx = myarray.DIMx;
float DIMy = myarray.DIMy;


TCanvas *c2 = new TCanvas("c2","The FillRandom example",200,50,500,300);
TCanvas *c4 = new TCanvas("c4","The FillRandom example",800,50,500,300);

TCanvas *c22 = new TCanvas("c22","The FillRandom example",200,350,500,300);
TCanvas *c44 = new TCanvas("c44","The FillRandom example",800,350,500,300);

TCanvas *c222 = new TCanvas("c222","The FillRandom example",200,650,500,300);
TCanvas *c444 = new TCanvas("c444","The FillRandom example",800,650,500,300);


string titlefile;
const char* c; 


//----------------------Pressure Flux---------------------//

TGraph2D *gr1 = new TGraph2D();
TGraph2D *gr5 = new TGraph2D();

titlefile = "P flow1; x; y; z";
c = titlefile.c_str();	
gr1->SetTitle(c);

titlefile = "P analytic; x; y; z";
c = titlefile.c_str();
gr5->SetTitle(c);		


int N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray.s1[i][j].P;

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
	  float T = myarray2.s1[i][j].P;

	      gr5->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}


//-------------u-------------//

TGraph2D *gr11 = new TGraph2D();
TGraph2D *gr51 = new TGraph2D();

titlefile = "u flow1; x; y; z";
c = titlefile.c_str();	
gr11->SetTitle(c);

titlefile = "u analytic; x; y; z";
c = titlefile.c_str();
gr51->SetTitle(c);	


N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray.s1[i][j].u;

	      gr11->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}


N = 0;
for (int i = 1; i < myarray2.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray2.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray2.s1[i][j].u;

	      gr51->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}

//--------------------v-------------------//

TGraph2D *gr111 = new TGraph2D();
TGraph2D *gr511 = new TGraph2D();

titlefile = "v flow1; x; y; z";
c = titlefile.c_str();	
gr111->SetTitle(c);

titlefile = "v analytic; x; y; z";
c = titlefile.c_str();
gr511->SetTitle(c);	


 N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray.s1[i][j].v;

	      gr111->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}


N = 0;
for (int i = 1; i < myarray2.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray2.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray2.s1[i][j].v;

	      gr511->SetPoint(N,dx,dy,T);
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
   gr1->Draw("surf1z");
   
c4->cd();   
   gStyle->SetPalette(1);
   gr5->Draw("surf1z");  
   
   c22->cd();   
   gStyle->SetPalette(1);
   gr11->Draw("surf1z");
   
c44->cd();   
   gStyle->SetPalette(1);
   gr51->Draw("surf1z"); 
   
   
   c222->cd();   
   gStyle->SetPalette(1);
   gr111->Draw("surf1z");
   
c444->cd();   
   gStyle->SetPalette(1);
   gr511->Draw("surf1z");  
   
}


//-------------------------------flux--------------------------------------//

void draw_3Dgraph_f(carray myarray, carray myarray2)
{
float DIMx = myarray.DIMx;
float DIMy = myarray.DIMy;


TCanvas *c2 = new TCanvas("c2","The FillRandom example",200,50,500,300);
TCanvas *c4 = new TCanvas("c4","The FillRandom example",800,50,500,300);

TCanvas *c22 = new TCanvas("c22","The FillRandom example",200,350,500,300);
TCanvas *c44 = new TCanvas("c44","The FillRandom example",800,350,500,300);

TCanvas *c222 = new TCanvas("c222","The FillRandom example",200,650,500,300);
TCanvas *c444 = new TCanvas("c444","The FillRandom example",800,650,500,300);


string titlefile;
const char* c; 


//----------------------Pressure Flux---------------------//

TGraph2D *gr1 = new TGraph2D();
TGraph2D *gr5 = new TGraph2D();

titlefile = "P flow1; x; y; z";
c = titlefile.c_str();	
gr1->SetTitle(c);

titlefile = "P analytic; x; y; z";
c = titlefile.c_str();
gr5->SetTitle(c);		


int N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray.f1[i][j].P;

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
	  float T = myarray2.f1[i][j].P;

	      gr5->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}


//-------------u-------------//

TGraph2D *gr11 = new TGraph2D();
TGraph2D *gr51 = new TGraph2D();

titlefile = "u flow1; x; y; z";
c = titlefile.c_str();	
gr11->SetTitle(c);

titlefile = "u analytic; x; y; z";
c = titlefile.c_str();
gr51->SetTitle(c);	


N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray.f1[i][j].u;

	      gr11->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}


N = 0;
for (int i = 1; i < myarray2.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray2.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray2.f1[i][j].u;

	      gr51->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}

//--------------------v-------------------//

TGraph2D *gr111 = new TGraph2D();
TGraph2D *gr511 = new TGraph2D();

titlefile = "v flow1; x; y; z";
c = titlefile.c_str();	
gr111->SetTitle(c);

titlefile = "v analytic; x; y; z";
c = titlefile.c_str();
gr511->SetTitle(c);	


 N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray.f1[i][j].v;

	      gr111->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}


N = 0;
for (int i = 1; i < myarray2.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray2.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray2.f1[i][j].v;

	      gr511->SetPoint(N,dx,dy,T);
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
   gr1->Draw("surf1z");
   
c4->cd();   
   gStyle->SetPalette(1);
   gr5->Draw("surf1z");  
   
   c22->cd();   
   gStyle->SetPalette(1);
   gr11->Draw("surf1z");
   
c44->cd();   
   gStyle->SetPalette(1);
   gr51->Draw("surf1z"); 
   
   
   c222->cd();   
   gStyle->SetPalette(1);
   gr111->Draw("surf1z");
   
c444->cd();   
   gStyle->SetPalette(1);
   gr511->Draw("surf1z");  
   
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

TGraph *gr = new TGraph();
gr->SetTitle(c);	

TGraph *gr5 = new TGraph();	
	gr5->SetMarkerColor(1);
	gr5->SetMarkerStyle(24);
	gr5->SetLineColor(1);
	gr5->SetTitle(c); 
	
	TGraph *gr51 = new TGraph();	
	gr51->SetMarkerColor(1);
	gr51->SetMarkerStyle(25);
	gr51->SetLineColor(1);
	gr51->SetTitle(c);  
	
	TGraph *gr52 = new TGraph();	
	gr52->SetMarkerColor(1);
	gr52->SetMarkerStyle(26);
	gr52->SetLineColor(1);
	gr52->SetTitle(c);   


int size = mydata.l2normP.size();
int n = 5;
printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
double temp = mydata.l2normP.at(i);
gr5->SetPoint(i,log(1.0/n),log(temp));

temp = mydata.l2normu.at(i);
gr51->SetPoint(i,log(1.0/n),log(temp));

temp = mydata.l2normv.at(i);
gr52->SetPoint(i,log(1.0/n),log(temp));

printf("l2P: %f change x: %d\n", temp, n);
n = n+5;
}


gr->SetPoint(0,-4,0);
gr->SetPoint(1,0,-10);

c5->SetGridx();
c5->SetGridy();
c5->SetTickx(1);
c5->SetTicky(1);

gr->Draw("AP");

TF1 *tfit2 = new TF1("tfit2", "pol1");
tfit2->SetLineColor(2);
tfit2->SetLineWidth(1);
gr5->Fit(tfit2, "", "", -10, 10);
gStyle->SetOptFit();
gr5->Draw("sameP");

TF1 *tfit21 = new TF1("tfit21", "pol1");
tfit21->SetLineColor(3);
tfit21->SetLineWidth(2);
gr51->Fit(tfit21, "", "", -10, 10);
gStyle->SetOptFit();
gr51->Draw("sameP");

TF1 *tfit22 = new TF1("tfit22", "pol1");
tfit22->SetLineColor(4);
tfit22->SetLineWidth(1);
gr52->Fit(tfit22, "", "", -10, 10);
gStyle->SetOptFit();
gr52->Draw("sameP");


TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(gr5,"P log(l2)","AP");
leg2->AddEntry(tfit2,"P Linear Fit","l");

leg2->AddEntry(gr51,"u log(l2)","AP");
leg2->AddEntry(tfit21,"u Linear Fit","l");

leg2->AddEntry(gr52,"v log(l2)","AP");
leg2->AddEntry(tfit22,"v Linear Fit","l");
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

