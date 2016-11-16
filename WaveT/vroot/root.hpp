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
#include "/home/jerin/cfd510/WaveT/numerical/numerical.hpp"

using namespace std;



void draw_graph(carray myarray, carray myarray2);

void draw_3DgraphP(carray myarray);

void draw_3Dgraph(carray myarray, carray myarray2);

void draw_graph_diff3(carray myarray, carray myarray2, carray myarray3);

void draw_graph_l2norm3(carray myarray, carray myarray2, carray myarray3);

//------------------------------------new=---------------------------//

void draw_graph_wave1_p2(carray & myarray1, carray myarray2, carray myarray3);

void draw_graph_wave1_p3(carray & myarray1, carray myarray2, carray myarray3);

void draw_graph_wave1(carray & myarray1, carray myarray2, carray myarray3);

void draw_graph_q1(carray & myarray1, carray & myarray2, carray & myarray3, carray analytic1, carray analytic2, carray analytic3);












#endif
