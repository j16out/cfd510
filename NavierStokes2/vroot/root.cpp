#include "root.hpp"

using namespace std;





//first one is finest array
float get_discrete_Error(carray ray1, carray ray2, carray ray3, double sol1, double sol2, double sol3)
{
//Calculating error as described in paper "procedure for estimation and reporting of uncertainty due to discretization in CFD applications"//

printf("\nCalculating Error...\n");

double h1 = 1.0/sqrt(ray1.sizex-2);
double h2 = 1.0/sqrt(ray2.sizex-2);
double h3 = 1.0/sqrt(ray3.sizex-2);


//sol1 = ray1.P1[40][20];
// sol2 = ray2.P1[20][10];
// sol3 = ray3.P1[10][5];


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

return GCI_21;
}
//---------------------------------------------------------------------------------------------3.3

void draw_um(carray myarray, carray myarray2, carray myarray3, carray myarray4)
{


string titlefile;
const char* c;  

TCanvas *c5 = new TCanvas("c5","The FillRandom example",200,50,900,700); 
//TCanvas *c7 = new TCanvas("c7","The FillRandom example",200,50,900,700);
c5->cd();
titlefile = "Order Evaluation; y vertical; Velocity u";
c = titlefile.c_str();

TGraph *gr = new TGraph();
gr->SetTitle(c);	

TGraph *gr5 = new TGraph();	
	gr5->SetMarkerColor(4);
	gr5->SetMarkerStyle(21);
	gr5->SetLineColor(1);
	gr5->SetTitle(c); 
	
	TGraph *gr51 = new TGraph();	
	gr51->SetMarkerColor(2);
	gr51->SetMarkerStyle(22);
	gr51->SetLineColor(1);
	gr51->SetTitle(c);  
	
	TGraph *gr515 = new TGraph();	
	gr515->SetMarkerColor(3);
	gr515->SetMarkerStyle(23);
	gr515->SetLineColor(1);
	gr515->SetTitle(c); 


titlefile = "Order Evaluation; y vertical; GCI Velocity u";
c = titlefile.c_str();
	
	TGraph *gr5155 = new TGraph();	
	gr5155->SetMarkerColor(4);
	gr5155->SetMarkerStyle(24);
	gr5155->SetLineColor(1);
	gr5155->SetTitle(c);
	
//c7->cd();

vector<float> gc1;
vector<float> gc2;
vector<float> gc3;


for(int j = 1; j < myarray.sizey-1; ++j)
{
float DIMx = myarray.DIMx;
float DIMy = myarray.DIMy;
int i = (myarray.sizex-2)/2;
float dx = DIMx*(i-0.5);
float dy = DIMy*(j-0.5);
 
double t1 = myarray.s1[i-1][j].u;
double t2 = myarray.s1[i][j].u;
double t3 = myarray.s1[i+1][j].u;
double t4 = myarray.s1[i+2][j].u;

TGraph *gr6 = new TGraph();	
         gr6->SetPoint(0,DIMx*(i-0.5-1),t1);
		 gr6->SetPoint(1,DIMx*(i-0.5),t2);
		 gr6->SetPoint(2,DIMx*(i-0.5+1),t3);
		 gr6->SetPoint(3,DIMx*(i-0.5+2),t4);

//printf("t1 %f t2 %f t3 %f t4 %f\n", t1, t2, t3, t4);
		
		TF1 *tfit2 = new TF1("tfit2", "pol2");
        gr6->Fit(tfit2, "BRQE", "", 0.4, 0.6);
        float p0 = tfit2->GetParameter(0);
        float p1 = tfit2->GetParameter(1);
        float p2 = tfit2->GetParameter(2);
        
        float u = p2*pow(0.5,2)+p1*(0.5)+p0;
//"BRQE"        
        gr5->SetPoint(j-1,DIMy*(j-0.5),u);   
gc1.push_back(u);
//gr6->Draw("A*c");
}


for(int j = 1; j < myarray2.sizey-1; ++j)
{
float DIMx = myarray2.DIMx;
float DIMy = myarray2.DIMy;
int i = (myarray2.sizex-2)/2;
float dx = DIMx*(i-0.5);
float dy = DIMy*(j-0.5);
 
double t1 = myarray2.s1[i-1][j].u;
double t2 = myarray2.s1[i][j].u;
double t3 = myarray2.s1[i+1][j].u;
double t4 = myarray2.s1[i+2][j].u;

TGraph *gr6 = new TGraph();	
         gr6->SetPoint(0,DIMx*(i-0.5-1),t1);
		 gr6->SetPoint(1,DIMx*(i-0.5),t2);
		 gr6->SetPoint(2,DIMx*(i-0.5+1),t3);
		 gr6->SetPoint(3,DIMx*(i-0.5+2),t4);

//printf("t1 %f t2 %f t3 %f t4 %f\n", t1, t2, t3, t4);
		
		TF1 *tfit2 = new TF1("tfit2", "pol2");
        gr6->Fit(tfit2, "BRQE", "", 0.4, 0.6);
        float p0 = tfit2->GetParameter(0);
        float p1 = tfit2->GetParameter(1);
        float p2 = tfit2->GetParameter(2);
        
        float u = p2*pow(0.5,2)+p1*(0.5)+p0;
//"BRQE"        
        gr51->SetPoint(j-1,DIMy*(j-0.5),u);   
gc2.push_back(u);
//gr6->Draw("A*c");
}


for(int j = 1; j < myarray3.sizey-1; ++j)
{
float DIMx = myarray3.DIMx;
float DIMy = myarray3.DIMy;
int i = (myarray3.sizex-2)/2;
float dx = DIMx*(i-0.5);
float dy = DIMy*(j-0.5);
 
double t1 = myarray3.s1[i-1][j].u;
double t2 = myarray3.s1[i][j].u;
double t3 = myarray3.s1[i+1][j].u;
double t4 = myarray3.s1[i+2][j].u;

TGraph *gr6 = new TGraph();	
         gr6->SetPoint(0,DIMx*(i-0.5-1),t1);
		 gr6->SetPoint(1,DIMx*(i-0.5),t2);
		 gr6->SetPoint(2,DIMx*(i-0.5+1),t3);
		 gr6->SetPoint(3,DIMx*(i-0.5+2),t4);

//printf("t1 %f t2 %f t3 %f t4 %f\n", t1, t2, t3, t4);
		
		TF1 *tfit2 = new TF1("tfit2", "pol2");
        gr6->Fit(tfit2, "BRQE", "", 0.4, 0.6);
        float p0 = tfit2->GetParameter(0);
        float p1 = tfit2->GetParameter(1);
        float p2 = tfit2->GetParameter(2);
        
        float u = p2*pow(0.5,2)+p1*(0.5)+p0;
//"BRQE"        
        gr515->SetPoint(j-1,DIMy*(j-0.5),u);   
gc3.push_back(u);
//gr6->Draw("A*c");
}


for(int j = 1; j < myarray4.sizey-1; ++j)
{
float DIMx = myarray4.DIMx;
float DIMy = myarray4.DIMy;
int i = (myarray4.sizex-2)/2;
float dx = DIMx*(i-0.5);
float dy = DIMy*(j-0.5);
 
double t1 = myarray4.s1[i-1][j].u;
double t2 = myarray4.s1[i][j].u;
double t3 = myarray4.s1[i+1][j].u;
double t4 = myarray4.s1[i+2][j].u;

TGraph *gr6 = new TGraph();	
         gr6->SetPoint(0,DIMx*(i-0.5-1),t1);
		 gr6->SetPoint(1,DIMx*(i-0.5),t2);
		 gr6->SetPoint(2,DIMx*(i-0.5+1),t3);
		 gr6->SetPoint(3,DIMx*(i-0.5+2),t4);

//printf("t1 %f t2 %f t3 %f t4 %f\n", t1, t2, t3, t4);
		
		TF1 *tfit2 = new TF1("tfit2", "pol2");
        gr6->Fit(tfit2, "BRQE", "", 0.4, 0.6);
        float p0 = tfit2->GetParameter(0);
        float p1 = tfit2->GetParameter(1);
        float p2 = tfit2->GetParameter(2);
        
        float u = p2*pow(0.5,2)+p1*(0.5)+p0;
//"BRQE"        
        //gr5155->SetPoint(j-1,DIMy*(j-0.5),u);   

//gr6->Draw("A*c");
}

for(int i = 0; i < myarray.sizey-2; ++i)
{
float DIMy = myarray.DIMy;
float sol1 = gc1.at(i);
float sol2 = gc2.at(i*3);
float sol3 = gc3.at(i*9);

float GCI = get_discrete_Error(myarray3, myarray2, myarray, sol3, sol2,  sol1);
gr5155->SetPoint(i,DIMy*(i-0.5+1), abs(GCI));
}

c5->cd();

gr->SetPoint(0,-4,0);
gr->SetPoint(1,0,-10);

c5->SetGridx();
c5->SetGridy();
c5->SetTickx(1);
c5->SetTicky(1);
//c5->SetLogy();
//gr->Draw("AP");



gr5->Draw("APc");
gr51->Draw("samePc");
gr515->Draw("samePc");

TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(gr5,"Velocity u 10x10","AP");
leg2->AddEntry(gr51,"Velocity u 30x30","AP");
leg2->AddEntry(gr515,"Velocity u 90x90","AP");

//leg2->AddEntry(gr51,"velocity u log(l2)","AP");
leg2->Draw();

TCanvas *c55 = new TCanvas("c55","The FillRandom example",200,50,900,700); 
c55->cd();
gr5155->Draw("aPl");


//gr51->Draw("samePc");





}

//------------------------------------------------------------------------3.22

void draw_3Dgraph_san(carray myarray, carray myarray2)
{
float DIMx = myarray.DIMx;
float DIMy = myarray.DIMy;




TCanvas *c22 = new TCanvas("c22","The FillRandom example",200,350,500,300);
TCanvas *c44 = new TCanvas("c44","The FillRandom example",800,350,500,300);
TCanvas *c222 = new TCanvas("c222","The FillRandom example",200,650,500,300);




string titlefile;
const char* c; 



//-------------u-------------//

TGraph2D *gr11 = new TGraph2D();
TGraph2D *gr51 = new TGraph2D();
TGraph2D *gr111 = new TGraph2D();

titlefile = "u flow1; x; y; u";
c = titlefile.c_str();	
gr11->SetTitle(c);

titlefile = "u analytic; x; y; u";
c = titlefile.c_str();
gr51->SetTitle(c);

titlefile = "v flow1; x; y; z";
c = titlefile.c_str();	
gr111->SetTitle(c);	


int N = 0;
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
	{ float dx = (DIMx*(i-0.5));
	  float dy = DIMy*(j-0.5);
	  float T = myarray2.s1[abs(myarray2.sizex-i-1)][j].u;

	      gr51->SetPoint(N,dx,dy,T);
	      ++N;
	 }     
}

N = 0;
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T = myarray.s1[i][j].u;
	  float T2 = myarray2.s1[abs(myarray2.sizex-i-1)][j].u;

	      gr111->SetPoint(N,dx,dy,T+T2);
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
   
 
   
   c22->cd();   
   gStyle->SetPalette(1);
   gr11->Draw("cont4z");
   
c44->cd();   
   gStyle->SetPalette(1);
   gr51->Draw("cont4z"); 
   
 


   c222->cd();   
   gStyle->SetPalette(1);
   gr111->Draw("cont4z");  

   
}


//-------------------------------------------------------------------------------------------------------3.21

void draw_u(carray myarray)
{


string titlefile;
const char* c;  

TCanvas *c5 = new TCanvas("c5","The FillRandom example",200,50,900,700); 
//TCanvas *c7 = new TCanvas("c7","The FillRandom example",200,50,900,700);
c5->cd();
titlefile = "Order Evaluation; y vertical; Velocity u";
c = titlefile.c_str();

TGraph *gr = new TGraph();
gr->SetTitle(c);	

TGraph *gr5 = new TGraph();	
	gr5->SetMarkerColor(4);
	gr5->SetMarkerStyle(21);
	gr5->SetLineColor(1);
	gr5->SetTitle(c); 
	
	TGraph *gr51 = new TGraph();	
	gr51->SetMarkerColor(2);
	gr51->SetMarkerStyle(22);
	gr51->SetLineColor(1);
	gr51->SetTitle(c);  
	
//c7->cd();

for(int j = 1; j < myarray.sizey-1; ++j)
{
float DIMx = myarray.DIMx;
float DIMy = myarray.DIMy;
int i = (myarray.sizex-2)/2;
float dx = DIMx*(i-0.5);
float dy = DIMy*(j-0.5);
 
double t1 = myarray.s1[i-1][j].u;
double t2 = myarray.s1[i][j].u;
double t3 = myarray.s1[i+1][j].u;
double t4 = myarray.s1[i+2][j].u;

TGraph *gr6 = new TGraph();	
         gr6->SetPoint(0,DIMx*(i-0.5-1),t1);
		 gr6->SetPoint(1,DIMx*(i-0.5),t2);
		 gr6->SetPoint(2,DIMx*(i-0.5+1),t3);
		 gr6->SetPoint(3,DIMx*(i-0.5+2),t4);

//printf("t1 %f t2 %f t3 %f t4 %f\n", t1, t2, t3, t4);
		
		TF1 *tfit2 = new TF1("tfit2", "pol2");
        gr6->Fit(tfit2, "BRQE", "", 0.4, 0.6);
        float p0 = tfit2->GetParameter(0);
        float p1 = tfit2->GetParameter(1);
        float p2 = tfit2->GetParameter(2);
        
        float u = p2*pow(0.5,2)+p1*(0.5)+p0;
//"BRQE"        
        gr5->SetPoint(j-1,DIMy*(j-0.5),u);   

//gr6->Draw("A*c");
}

c5->cd();

gr->SetPoint(0,-4,0);
gr->SetPoint(1,0,-10);

c5->SetGridx();
c5->SetGridy();
c5->SetTickx(1);
c5->SetTicky(1);
//c5->SetLogy();
//gr->Draw("AP");



gr5->Draw("APc");


//gr51->Draw("samePc");


TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(gr5,"Velocity u at x = 0.5","AP");
//leg2->AddEntry(gr51,"velocity u log(l2)","AP");


leg2->Draw();


}

void draw_stab_l2(cdata & mydata)
{
string titlefile;
const char* c;  

TCanvas *c5 = new TCanvas("c5","The FillRandom example",200,50,900,700); 
c5->cd();
titlefile = "Order Evaluation; Iteration; Log(L2)";
c = titlefile.c_str();

TGraph *gr = new TGraph();
gr->SetTitle(c);	

TGraph *gr5 = new TGraph();	
	gr5->SetMarkerColor(4);
	gr5->SetMarkerStyle(24);
	gr5->SetLineColor(1);
	gr5->SetTitle(c); 
	
	TGraph *gr51 = new TGraph();	
	gr51->SetMarkerColor(2);
	gr51->SetMarkerStyle(25);
	gr51->SetLineColor(1);
	gr51->SetTitle(c);  
	
	TGraph *gr52 = new TGraph();	
	gr52->SetMarkerColor(3);
	gr52->SetMarkerStyle(26);
	gr52->SetLineColor(1);
	gr52->SetTitle(c);   


int size = mydata.l2normP.size();

printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
double temp = mydata.l2normP.at(i);
gr5->SetPoint(i,i,log(temp));

temp = mydata.l2normu.at(i);
gr51->SetPoint(i,i,log(temp));

temp = mydata.l2normv.at(i);
gr52->SetPoint(i,i,log(temp));


}


gr->SetPoint(0,-4,0);
gr->SetPoint(1,0,-10);

c5->SetGridx();
c5->SetGridy();
c5->SetTickx(1);
c5->SetTicky(1);
//c5->SetLogy();
//gr->Draw("AP");



gr5->Draw("APc");


gr51->Draw("samePc");


gr52->Draw("samePc");


TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(gr5,"Pressure log(l2)","AP");
leg2->AddEntry(gr51,"velocity u log(l2)","AP");
leg2->AddEntry(gr52,"velocity v log(l2)","AP");

leg2->Draw();


}

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

titlefile = "P flow1; x; y; P";
c = titlefile.c_str();	
gr1->SetTitle(c);

titlefile = "P analytic; x; y; P";
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

titlefile = "u flow1; x; y; u";
c = titlefile.c_str();	
gr11->SetTitle(c);

titlefile = "u analytic; x; y; u";
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

titlefile = "v flow1; x; y; v";
c = titlefile.c_str();	
gr111->SetTitle(c);


titlefile = "v analytic; x; y; v";
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
   gr1->Draw("cont4z");
   
c4->cd();   
   gStyle->SetPalette(1);
   gr1->Draw("surf1z");  
   
   c22->cd();   
   gStyle->SetPalette(1);
   gr11->Draw("cont4z");
   
c44->cd();   
   gStyle->SetPalette(1);
   gr11->Draw("surf1z"); 
   
   
   c222->cd();   
   gStyle->SetPalette(1);
   gr111->Draw("cont4z");
   
c444->cd();   
   gStyle->SetPalette(1);
   gr111->Draw("surf1z");  
   
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

