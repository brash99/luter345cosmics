#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TROOT.h"




//Global Variables
bool remakepedfile= false;
const int bin = 100;
const Int_t Nadc = 16;
const Int_t Ntdc = 16;
const Int_t pedrun = 172;
const Int_t bl=Ntdc-16+0;//Start at channel 0
const Int_t br=Ntdc-16+1;
const Int_t tl=Ntdc-16+2;
const Int_t tr=Ntdc-16+3;
const Double_t adjadcto=1400.0;//value to ADJust ADC TO
Double_t t_fullscale = 140.0E-09; // full scale TDC range in seconds
Double_t t_convert=t_fullscale/4096.0;

const Int_t nthetabins = 51;
const Double_t thetalow = -100.5;
const Double_t thetahigh = 100.5;


Double_t vn = 2.997E08/nscint;
Double_t resolution = 0.0232*nscint*nscint-0.1061*nscint+0.1617;
resolution = resolution*0.934;
Double_t granularity = t_convert*vn/2.0;
Double_t xpos_range = 0.30;
Double_t dscint = 0.105; // distance between scintillators in metres
const int xposbin = 2.0*xpos_range/granularity;

