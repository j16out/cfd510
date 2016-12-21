#ifndef ROOT_INCLUDED
#define ROOT_INCLUDED

#include "TF1.h"
#include <TAxis.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include "TF1.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TApplication.h"
#include <TLatex.h>
#include <TImage.h>
#include <TRandom3.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TPad.h>
#include <TFrame.h>
#include "TH3.h"
#include "TNtuple.h"

#include <TRandom3.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <math.h> 
#include <thread>
#include <stdio.h>
#include <iterator>
#include <sys/stat.h>
#include <unistd.h>
#include "/home/jerin/cfd510/NavierStokes/numerical/numerical.hpp"

using namespace std;




void draw_3Dgraph_s(carray myarray, carray myarray2);
void draw_3Dgraph_f(carray myarray, carray myarray2);

void draw_order_l2(cdata & mydata);

void draw_order_eff(cdata & mydata);

void find_maxvalues(carray & myarray, carray & myarray2, carray & myarray3);

void draw_stab_l2(cdata & mydata);

void draw_u(carray myarray);

void draw_3Dgraph_san(carray myarray, carray myarray2);

void draw_um(carray myarray, carray myarray2, carray myarray3, carray myarray4);

float get_discrete_Error(carray ray1, carray ray2, carray ray3, double sol1, double sol2, double sol3);

void get_vortex(carray myarray);















#endif
