#include "root.hpp"

using namespace std;










//---------------------------------------------------------------------------------------------4/3


void draw_um(carray myarray, carray myarray2, carray myarray3, carray myarray4)
{

 

TCanvas *c59 = new TCanvas("c59","The FillRandom example",200,50,900,700); 
TCanvas *c535 = new TCanvas("c535","The FillRandom example",200,50,900,700); 
//TCanvas *c7 = new TCanvas("c7","The FillRandom example",200,50,900,700);
c59->cd();


string titlefile;
const char* c; 

titlefile = "Order Evaluation; y vertical; Velocity u";
c = titlefile.c_str();

TGraph *gr2 = new TGraph();
gr2->SetTitle(c);	

TGraph *gr52 = new TGraph();	
	gr52->SetMarkerColor(4);
	gr52->SetMarkerStyle(21);
	gr52->SetLineColor(1);
	gr52->SetTitle(c); 
	
	TGraph *gr512 = new TGraph();	
	gr512->SetMarkerColor(2);
	gr512->SetMarkerStyle(22);
	gr512->SetLineColor(1);
	gr512->SetTitle(c);  
	
	TGraph *gr5152 = new TGraph();	
	gr5152->SetMarkerColor(3);
	gr5152->SetMarkerStyle(23);
	gr5152->SetLineColor(1);
	gr5152->SetTitle(c); 


titlefile = "Order Evaluation; y vertical; GCI Velocity u";
c = titlefile.c_str();
	
	TGraph *gr51552 = new TGraph();	
	gr51552->SetMarkerColor(4);
	gr51552->SetMarkerStyle(24);
	gr51552->SetLineColor(1);
	gr51552->SetTitle(c);
	
//c7->cd();

vector<float> gc1;
vector<float> gc2;
vector<float> gc3;


for(int j = 1; j < myarray.sizey-1; ++j)
{
float DIMx = myarray.DIMx;
float DIMy = myarray.DIMy;
int i = (myarray.sizex-2)/2.0;
float dx = DIMx*(i-0.5);
float dy = DIMy*(j-0.5);
 
double t1 = myarray.s1[i-1][j].u;
double t2 = myarray.s1[i][j].u;
double t3 = myarray.s1[i+1][j].u;
double t4 = myarray.s1[i+2][j].u;

TGraph *gr62 = new TGraph();	
         gr62->SetPoint(0,DIMx*(i-0.5-1),t1);
		 gr62->SetPoint(1,DIMx*(i-0.5),t2);
		 gr62->SetPoint(2,DIMx*(i-0.5+1),t3);
		 gr62->SetPoint(3,DIMx*(i-0.5+2),t4);

//printf("t1 %f t2 %f t3 %f t4 %f\n", t1, t2, t3, t4);
		
		TF1 *tfit2 = new TF1("tfit2", "pol2");
        gr62->Fit(tfit2, "BRQE", "", 0.4, 0.6);
        float p0 = tfit2->GetParameter(0);
        float p1 = tfit2->GetParameter(1);
        float p2 = tfit2->GetParameter(2);
        
        float u = p2*pow(0.5,2)+p1*(0.5)+p0;
//"BRQE"        
        gr52->SetPoint(j-1,DIMy*(j-0.5),u);   
gc1.push_back(u);

}


for(int j = 1; j < myarray2.sizey-1; ++j)
{
float DIMx = myarray2.DIMx;
float DIMy = myarray2.DIMy;
int i = (myarray2.sizex-2)/2.0;
float dx = DIMx*(i-0.5);
float dy = DIMy*(j-0.5);
 
double t1 = myarray2.s1[i-1][j].u;
double t2 = myarray2.s1[i][j].u;
double t3 = myarray2.s1[i+1][j].u;
double t4 = myarray2.s1[i+2][j].u;

TGraph *gr61 = new TGraph();	
         gr61->SetPoint(0,DIMx*(i-0.5-1),t1);
		 gr61->SetPoint(1,DIMx*(i-0.5),t2);
		 gr61->SetPoint(2,DIMx*(i-0.5+1),t3);
		 gr61->SetPoint(3,DIMx*(i-0.5+2),t4);

//printf("t1 %f t2 %f t3 %f t4 %f\n", t1, t2, t3, t4);
		
		TF1 *tfit2 = new TF1("tfit2", "pol2");
        gr61->Fit(tfit2, "BRQE", "", 0.4, 0.6);
        float p0 = tfit2->GetParameter(0);
        float p1 = tfit2->GetParameter(1);
        float p2 = tfit2->GetParameter(2);
        
        float u = p2*pow(0.5,2)+p1*(0.5)+p0;
//"BRQE"        
        gr512->SetPoint(j-1,DIMy*(j-0.5),u);   
gc2.push_back(u);
//gr6->Draw("A*c");
}


for(int j = 1; j < myarray3.sizey-1; ++j)
{
float DIMx = myarray3.DIMx;
float DIMy = myarray3.DIMy;
int i = (myarray3.sizex-2)/2.0;
float dx = DIMx*(i-0.5);
float dy = DIMy*(j-0.5);
 
double t1 = myarray3.s1[i-1][j].u;
double t2 = myarray3.s1[i][j].u;
double t3 = myarray3.s1[i+1][j].u;
double t4 = myarray3.s1[i+2][j].u;

TGraph *gr63 = new TGraph();	
         gr63->SetPoint(0,DIMx*(i-0.5-1),t1);
		 gr63->SetPoint(1,DIMx*(i-0.5),t2);
		 gr63->SetPoint(2,DIMx*(i-0.5+1),t3);
		 gr63->SetPoint(3,DIMx*(i-0.5+2),t4);

//printf("t1 %f t2 %f t3 %f t4 %f\n", t1, t2, t3, t4);
		
		TF1 *tfit2 = new TF1("tfit2", "pol2");
        gr63->Fit(tfit2, "BRQE", "", 0.4, 0.6);
        float p0 = tfit2->GetParameter(0);
        float p1 = tfit2->GetParameter(1);
        float p2 = tfit2->GetParameter(2);
        
        float u = p2*pow(0.5,2)+p1*(0.5)+p0;
//"BRQE"        
        gr5152->SetPoint(j-1,DIMy*(j-0.5),u);   
gc3.push_back(u);
//gr6->Draw("A*c");
}


for(int j = 1; j < myarray4.sizey-1; ++j)
{
float DIMx = myarray4.DIMx;
float DIMy = myarray4.DIMy;
int i = (myarray4.sizex-2)/2.0;
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

c59->cd();

gr2->SetPoint(0,-4,0);
gr2->SetPoint(1,0,-10);

c59->SetGridx();
c59->SetGridy();
c59->SetTickx(1);
c59->SetTicky(1);
//c5->SetLogy();
//gr->Draw("AP");



gr52->Draw("APc");
gr512->Draw("samePc");
gr5152->Draw("samePc");

TLegend *leg2 = new TLegend(0.75,0.9,0.9,0.8);

leg2->AddEntry(gr52,"Velocity u 10x30","AP");
leg2->AddEntry(gr512,"Velocity u 20x60","AP");
leg2->AddEntry(gr5152,"Velocity u 40x120","AP");

//leg2->AddEntry(gr51,"velocity u log(l2)","AP");
leg2->Draw();


c535->cd();
gr51552->Draw("aPl");

}
//--------------------------------------------------------------------------------------------
 
void get_vortex(carray myarray)
{
TCanvas *ca = new TCanvas("ca","The FillRandom example",200,50,900,700); 
string titlefile;
const char* c; 


titlefile = "Order Evaluation; x; y";
c = titlefile.c_str();
	TGraph *gr = new TGraph();
	gr->SetTitle(c);
	
	TGraph *gr5a = new TGraph();	
	gr5a->SetMarkerColor(4);
	gr5a->SetMarkerStyle(24);
	gr5a->SetLineColor(1);
	gr5a->SetTitle(c);
	
	TGraph *gr5b = new TGraph();	
	gr5b->SetMarkerColor(2);
	gr5b->SetMarkerStyle(24);
	gr5b->SetLineColor(1);
	gr5b->SetTitle(c);
	
float q = 1.0;
//adjust j to shorten
for(int j = q*109; j < q*114; ++j)
{
float DIMx = myarray.DIMx;
float DIMy = myarray.DIMy;


float dy = DIMy*(j-0.5);
 
 int n = 0;
TGraph *gr6 = new TGraph();	 
for(int i = q*16; i < q*24; ++i)
{
 float dx = DIMx*(i-0.5);
 float T1 = myarray.s1[i][j].v;
 float T2 = myarray.s1[i][j].u;
 float T = sqrt( pow(T1,2) + pow(T2,2));
gr6->SetPoint(n,dx,T);
n++;
}

//printf("t1 %f t2 %f t3 %f t4 %f\n", t1, t2, t3, t4);
TF1 *tfit2 = new TF1("tfit2", "pol2");
tfit2->SetParameter(2, 0.5);
tfit2->SetParLimits(2, 0, 50);
gr6->Fit(tfit2, "", "", 0.0, 3.0);
float p0 = tfit2->GetParameter(0);
float p1 = tfit2->GetParameter(1);
float p2 = tfit2->GetParameter(2);

float apex = -p1/(2.0*p2);	
	
//"BRQE"
        
gr5a->SetPoint(j-q*109,apex, DIMy*(j-0.5));   


}

//adjust i to shorten


for(int i = q*17; i < q*24; ++i)
{
float DIMx = myarray.DIMx;
float DIMy = myarray.DIMy;


float dx = DIMx*(i-0.5);
 
 int n = 0;
TGraph *gr6 = new TGraph();	 
for(int j = q*106; j < q*114; ++j)
{
 float dy = DIMy*(j-0.5);
 float T1 = myarray.s1[i][j].v;
 float T2 = myarray.s1[i][j].u;
 float T = sqrt( pow(T1,2) + pow(T2,2));
gr6->SetPoint(n,dy,T);
n++;
}

//printf("t1 %f t2 %f t3 %f t4 %f\n", t1, t2, t3, t4);
TF1 *tfit2 = new TF1("tfit2", "pol2");
tfit2->SetParameter(2, 0.5);
tfit2->SetParLimits(2, 0.0001, 2.0);

gr6->Fit(tfit2, "", "", 2.5, 3.0);
float p0 = tfit2->GetParameter(0);
float p1 = tfit2->GetParameter(1);
float p2 = tfit2->GetParameter(2);

float apex = -p1/(2.0*p2);	
	
//"BRQE"        
gr5b->SetPoint(i-q*17,DIMx*(i-0.5), apex);   


}

gr->SetPoint(0,0.4,2.8);
gr->SetPoint(1,0.6,2.6);

TF1 *tfit1 = new TF1("tfit1", "pol1");
gr5a->Fit(tfit1, "", "", 0.4, 0.6);
float cc = tfit1->GetParameter(0);
float a = tfit1->GetParameter(1);


TF1 *tfit11 = new TF1("tfit11", "pol1");
gr5b->Fit(tfit11, "", "", 0.0, 3.0);
float d = tfit11->GetParameter(0);
float b = tfit11->GetParameter(1);
gr->Draw("ap");
gr5a->Draw("samePc");
gr5b->Draw("samePc");

float x = (d-cc)/(a-b);
float y = a*((d-cc)/(a-b))+cc;


printf("intercept x = %f y = %f\n",x, y); 







}

//first one is finest array




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

vector<double> u1;
vector<double> d1;


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
   u1.push_back(u); 
   d1.push_back(DIMy*(j-0.5));     
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

//--------------------------------------------------vortice strength
while(1)
{
float rr1 = 0.0;
float rr2 = 0.8;
int r1 = (myarray.sizex-2)*rr1;
int r2 = (myarray.sizex-2)*rr2;

float n = 0;
float p = 0;

float sum = 0;
float sum2 = 0;
for(int i = r1; i < r2; ++i)
{

float temp = u1.at(i)*myarray.DIMy;
//printf("i %d temp %f\n", i , temp );
if(temp < 0)
{
sum2 = temp + sum2;
++n;
}
else
{
sum = temp + sum;
++p;
}
}
float avep = sum/p;

float aven = sum/n;
printf("vortice #1 strength p %f strength n %f averages p %f, n %f\n",sum, sum2, avep, aven); 
break;
}

while(1)
{
float rr1 = 0.8;
float rr2 = 1.8;
int r1 = (myarray.sizex-2)*rr1;
int r2 = (myarray.sizex-2)*rr2;

float n = 0;
float p = 0;

float sum = 0;
float sum2 = 0;
for(int i = r1; i < r2; ++i)
{

float temp = u1.at(i)*myarray.DIMy;
//printf("i %d temp %f\n", i , temp );
if(temp < 0)
{
sum2 = temp + sum2;
++n;
}
else
{
sum = temp + sum;
++p;
}
}
float avep = sum/p;

float aven = sum/n;
printf("vortice #2 strength p %f strength n %f averages p %f, n %f\n",sum, sum2, avep, aven); 
break;
}

while(1)
{
float rr1 = 1.8;
float rr2 = 3.0;
int r1 = (myarray.sizex-2)*rr1;
int r2 = (myarray.sizex-2)*rr2;

float n = 0;
float p = 0;

float sum = 0;
float sum2 = 0;
for(int i = r1; i < r2; ++i)
{

float temp = u1.at(i)*myarray.DIMy;
//printf("i %d temp %f\n", i , temp );
if(temp < 0)
{
sum2 = temp + sum2;
++n;
}
else
{
sum = temp + sum;
++p;
}
}
float avep = sum/p;

float aven = sum/n;
printf("vortice #3 strength p %f strength n %f averages p %f, n %f\n",sum, sum2, avep, aven); 
break;
}


}







//--------------------------------------------------------------------------




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

	      //gr51->SetPoint(N,dx,dy,T);
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
for (int i = 1; i < myarray.sizex-1; i++) 
{
	
	for (int j = 1; j < myarray.sizey-1; j++) 
	{ float dx = DIMx*(i-0.5);
	  float dy = DIMy*(j-0.5);
	  float T1 = myarray.s1[i][j].v;
      float T2 = myarray.s1[i][j].u;
      float T = sqrt( pow(T1,2) + pow(T2,2));
	      gr511->SetPoint(N,dx,dy,T);
	      gr51->SetPoint(N,dx,dy,T);
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
   gr51->Draw("surf1z");  
   
   c22->cd();   
   gStyle->SetPalette(1);
   gr11->Draw("cont4z");
   
c44->cd();   
   gStyle->SetPalette(1);
   gr51->Draw("cont1z"); 
   
   
   c222->cd();   
   gStyle->SetPalette(1);
   gr111->Draw("cont4z");
   
c444->cd();   
   gStyle->SetPalette(1);
   gr511->Draw("cont4z");  
   
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
float n = 5;
printf("array size %d\n", size);

for(int i = 0; i < size; ++i)
{
double temp = mydata.l2normP.at(i);
gr5->SetPoint(i,log(1.0/(n)),log(temp));

temp = mydata.l2normu.at(i);
gr51->SetPoint(i,log(1.0/(n)),log(temp));

temp = mydata.l2normv.at(i);
gr52->SetPoint(i,log(1.0/(n)),log(temp));

printf("l2P: %f change x: %f\n", temp, n);
n = n+1;
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
gr5->Fit(tfit2, "", "", -100, 100);
gStyle->SetOptFit();
gr5->Draw("sameP");

TF1 *tfit21 = new TF1("tfit21", "pol1");
tfit21->SetLineColor(3);
tfit21->SetLineWidth(2);
gr51->Fit(tfit21, "", "", -100, 100);
gStyle->SetOptFit();
gr51->Draw("sameP");

TF1 *tfit22 = new TF1("tfit22", "pol1");
tfit22->SetLineColor(4);
tfit22->SetLineWidth(1);
gr52->Fit(tfit22, "", "", -100, 100);
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

