#include <iostream>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TLinearFitter.h>
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

using namespace std;
using RNode = ROOT::RDF::RNode;

std::vector<RNode> v;

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
const Int_t tdc_min = 500;
const Int_t tdc_max = 3500;
const Double_t adjadcto=1400.0;//value to ADJust ADC TO
Double_t t_fullscale = 140.0E-09; // full scale TDC range in seconds
Double_t t_convert=t_fullscale/4096.0;

const Int_t nthetabins = 51;
const Double_t thetalow = -100.5;
const Double_t thetahigh = 100.5;

const Double_t nscint=1.50;

Double_t vn = 2.997E08/nscint;
Double_t resolution = 0.934*(0.0232*nscint*nscint-0.1061*nscint+0.1617);
Double_t granularity = t_convert*vn/2.0;
Double_t xpos_range = 0.30;
Double_t dscint = 0.105; // distance between scintillators in metres
const int xposbin = 2.0*xpos_range/granularity;

Char_t tdcnames[][Ntdc]={"Bottom Left","Bottom Right","Top Left","Top Right","4","5","6","7","8","9","10","11","12","13","14","15"};
Char_t adcnames[][Nadc]={"Bottom Left","Bottom Right","Top Left","Top Right","4","5","6","7","8","9","10","11","12","13","14","15"};

//Set correction values
Double_t tdccorrect[Ntdc];
Double_t ped[Nadc];
Double_t gain[Nadc];

TRandom r;
Double_t rnd;

bool getTrigger(int* tdc, int* adc) {

	bool tdc_tl = false;
	bool tdc_tr = false;
	bool tdc_bl = false;
	bool tdc_br = false;
	bool trigger = false;

	if (tdc[tl]>tdc_min&&tdc[tl]<tdc_max) {
		tdc_tl = true;
	}
	if (tdc[tr]>tdc_min&&tdc[tr]<tdc_max) {
		tdc_tr = true;
	}
	if (tdc[bl]>tdc_min&&tdc[bl]<tdc_max) {
		tdc_bl = true;
	}
	if (tdc[br]>tdc_min&&tdc[br]<tdc_max) {
		tdc_br = true;
	}

	if (tdc_tl && tdc_tr && tdc_bl && tdc_br) {
		trigger = true;
	}

	return trigger;
}

TCanvas* calculateADCFactors() {
	
	ROOT::RDF::RResultPtr<TH1D> hadcraw[Ntdc];
	TCanvas *c3 = new TCanvas("c3", "c3", 100,100,600,600);
	c3->Divide(2,2, 0.01, 0.01, 0);

	//=====================================GET PED for ADC==========================
        if(remakepedfile){cout << "Executing luterpedestals.C" << endl; gROOT->ProcessLine(Form(".x luterpedestals.C(%d)",pedrun));cout<<"pedestal file made"<<endl;}

        //cout << "Opening pedestal values file ... " << endl;
        FILE *adcpeds = fopen(Form("./pedestalfiles/pedestalrun%d.dat",pedrun),"r");

        //cout << "Filling ADC pedestal array ..." << endl;
        for( int i = 0; i < Nadc ; i++) {//Start ADC filling loop
                  fscanf(adcpeds,"%lf\n",&ped[i]);
	}
	
	for (int i = bl; i <= tr ; i++) {
		auto varname = "adc[" + std::to_string(i) + "]";
		hadcraw[i] = v[0].Define(Form("hadc%02d",i),varname).Histo1D({Form("hadc%02d",i), adcnames[i], bin, 0, 3000.0},Form("hadc%02d",i));
		gain[i] = adjadcto/(hadcraw[i]->GetMean()-ped[i]);
		cout << i << " " << hadcraw[i]->GetMean() << " " << ped[i] << " " << gain[i] << endl;
		c3->cd(i+1);
		hadcraw[i]->DrawClone();
	}
	
	return c3;

}

TCanvas* calculateTOffsets(){

	ROOT::RDF::RResultPtr<TH1D> htdcraw[Ntdc];
	TCanvas *c1 = new TCanvas("c1", "c1", 75,75,600,600);
	c1->Divide(2,2, 0.01, 0.01, 0);

	for (int i = bl; i <= tr ; i++) {
		auto varname = "tdc[" + std::to_string(i) + "]";
		htdcraw[i] = v[0].Define(Form("htdc%02d",i),varname).Histo1D({Form("htdc%02d",i), tdcnames[i], bin, 1700, 2100.0},Form("htdc%02d",i));
		cout << i << " " << htdcraw[i]->GetMean() << " " << htdcraw[i]->GetStdDev() << endl;
		c1->cd(i+1);
		int highbin=5;
		for (int ii=0;ii<=bin;ii++) {
			if ((htdcraw[i]->GetBinContent(ii))>(htdcraw[i]->GetBinContent(highbin))&&ii>=5){
                                highbin=ii;
                        }
                }
                int median = htdcraw[i]->GetBinCenter(highbin);
                htdcraw[i]->Fit("gaus","Q","",(median-200),(median+200));
                tdccorrect[i] = 2000 - (htdcraw[i]->GetFunction("gaus")->GetParameter(1));
		htdcraw[i]->DrawClone();
	}

	return c1;

}

float getTheta(bool trigger, float xtop, float xbottom) {

	double rtod=180.0/3.14159265;
	if (trigger) {
		rnd = r.Gaus(0.0,1.5);
		return rtod*atan((xbottom-xtop)/dscint)+rnd;
	}

	return -1000;

}

float getTheta2(bool trigger, float xtop, float xbottom) {

	double rtod=180.0/3.14159265;
        Double_t rnd_top_pos = r.Gaus(0.0,resolution);
        Double_t rnd_bottom_pos = r.Gaus(0.0,resolution);

        TF1 *fphi = new TF1("fsin", "x", 0, 2*TMath::Pi());
        TF1 *fcos = new TF1("fcos", "cos(x)*cos(x)",-TMath::Pi()/2.0,TMath::Pi()/2.0);
        Double_t phisim, cossim, rndxt, rndyt, rndrt, rndrr, rndxb, rndyb;
        while (true) {
            phisim = fphi->GetRandom();
            cossim = fcos->GetRandom();
            rndxt = r.Uniform(-0.10,0.10);
            rndyt = r.Uniform(-0.15,0.15);
            rndrt = sqrt(rndxt*rndxt+rndyt*rndyt)+rnd_top_pos;
            rndrr = dscint*tan(cossim);
            rndxb = rndxt+rndrr*cos(phisim);
            rndyb = rndyt+rndrr*sin(phisim);
            if (rndxb > -0.10 && rndxb < 0.10 && rndyb > -0.15 && rndyb < 0.15) break;
        }

        Double_t rndrb = sqrt(rndxb*rndxb+rndyb*rndyb)+rnd_bottom_pos;
	Double_t theta2 = 0.8281-rtod*atan((rndrt-rndrb)/dscint);	

	if (trigger) {
		return theta2;
	}

	return -1000;

}

float getAdcTL(bool trigger, int* adc) {

	//std::vector<float> v;

	if(trigger) {
		//v.push_back(gain[tl]*(adc[tl]-ped[tl]));
		return gain[tl]*(adc[tl]-ped[tl]);
        }

	return -1;
}

float getAdcTR(bool trigger, int* adc) {

	//std::vector<float> v;

	if(trigger) {
		return gain[tr]*(adc[tr]-ped[tr]);
        }

	return -1;
}

float getAdcBL(bool trigger, int* adc) {

	//std::vector<float> v;

	if(trigger) {
		return gain[bl]*(adc[bl]-ped[bl]);
        }

	return -1;
}

float getAdcBR(bool trigger, int* adc) {

	//std::vector<float> v;

	if(trigger) {
		return gain[br]*(adc[br]-ped[br]);
        }

	return -1;
}

float getTdcTL(bool trigger, int* tdc) {

	//std::vector<float> v;

	if(trigger) {
		//v.push_back(tdc[tl]+tdccorrect[tl]);
		return tdc[tl]+tdccorrect[tl];
        }

	return -1;
}

float getTdcTR(bool trigger, int* tdc) {

	//std::vector<float> v;

	if(trigger) {
		//v.push_back(tdc[tr]+tdccorrect[tr]);
		return tdc[tr]+tdccorrect[tr];
        }

	return -1;
}

float getTdcBL(bool trigger, int* tdc) {

	//std::vector<float> v;

	if(trigger) {
		//v.push_back(tdc[bl]+tdccorrect[bl]);
		return tdc[bl]+tdccorrect[bl];
        }

	return -1;
}

float getTdcBR(bool trigger, int* tdc) {

	//std::vector<float> v;

	if(trigger) {
		//v.push_back(tdc[br]+tdccorrect[br]);
		return tdc[br]+tdccorrect[br];
        }

	return -1;
}

float getXTop(bool trigger, double ttl, double ttr, double tbl, double tbr) {

	//std::vector<float> v;

	if(trigger) {
		return (ttl-ttr)/2.0*t_convert*vn;
        }

	return -1000;
}

float getXBottom(bool trigger, double ttl, double ttr, double tbl, double tbr) {

	//std::vector<float> v;

	if(trigger) {
		return (tbl-tbr)/2.0*t_convert*vn;
        }

	return -1000;
}

float getERatio(bool trigger, double etop, double ebottom) {

	if (trigger) {
		return etop/ebottom;
	}

	return -1000;
}

float getETop(bool trigger, double adctl, double adctr) {

	if (trigger) {
		return (adctr+adctl)/2.0;
	}

	return -1000;
}

float getEBottom(bool trigger, double adcbl, double adcbr) {

	if (trigger) {
		return (adcbr+adcbl)/2.0;
	}

	return -1000;
}

TCanvas* plotTDCvsXpos(){

	TCanvas *c6 = new TCanvas("c6", "c6", 112,112,600,600);
	c6->Divide(2,4, 0.01, 0.01, 0);

	auto hTDCTLvsXT = v[1].Histo2D({"h1","TDC Top Left vs XTop",bin,-0.3,0.3,bin,1950,2050},"xtop","tdctl");
	auto hTDCTRvsXT = v[1].Histo2D({"h2","TDC Top Right vs XTop",bin,-0.3,0.3,bin,1950,2050},"xtop","tdctr");
	auto hTDCTLvsXB = v[1].Histo2D({"h3","TDC Top Left vs XBottom",bin,-0.3,0.3,bin,1950,2050},"xbottom","tdctl");
	auto hTDCTRvsXB = v[1].Histo2D({"h4","TDC Top Right vs XBottom",bin,-0.3,0.3,bin,1950,2050},"xbottom","tdctr");
	auto hTDCBLvsXT = v[1].Histo2D({"h1","TDC Bottom Left vs XTop",bin,-0.3,0.3,bin,1950,2050},"xtop","tdcbl");
	auto hTDCBRvsXT = v[1].Histo2D({"h2","TDC Bottom Right vs XTop",bin,-0.3,0.3,bin,1950,2050},"xtop","tdcbr");
	auto hTDCBLvsXB = v[1].Histo2D({"h3","TDC Bottom Left vs XBottom",bin,-0.3,0.3,bin,1950,2050},"xbottom","tdcbl");
	auto hTDCBRvsXB = v[1].Histo2D({"h4","TDC Bottom Right vs XBottom",bin,-0.3,0.3,bin,1950,2050},"xbottom","tdcbr");

	c6->cd(1);
 	hTDCTLvsXT->Draw("COLZ");
 	c6->cd(2);
 	hTDCTRvsXT->Draw("COLZ");
  	c6->cd(3);
  	hTDCTLvsXB->Draw("COLZ");
  	c6->cd(4);
  	hTDCTRvsXB->Draw("COLZ");
	c6->cd(5);
 	hTDCBLvsXT->Draw("COLZ");
 	c6->cd(6);
 	hTDCBRvsXT->Draw("COLZ");
  	c6->cd(7);
  	hTDCBLvsXB->Draw("COLZ");
  	c6->cd(8);
  	hTDCBRvsXB->Draw("COLZ");

	c6->DrawClone();

	return c6;

}
TCanvas* plotADCvsXpos(){

	TCanvas *c5 = new TCanvas("c5", "c5", 125,125,600,600);
	c5->Divide(2,2, 0.01, 0.01, 0);

	auto hADCvsXTL = v[1].Histo2D({"h1","ADC Top Left vs XTop",bin,-0.3,0.3,bin,0,2500},"xtop","adctl");
	auto hADCvsXTR = v[1].Histo2D({"h2","ADC Top Right vs XTop",bin,-0.3,0.3,bin,0,2500},"xtop","adctr");
	auto hADCvsXBL = v[1].Histo2D({"h3","ADC Bottom Left vs XBottom",bin,-0.3,0.3,bin,0,2500},"xbottom","adcbl");
	auto hADCvsXBR = v[1].Histo2D({"h4","ADC Bottom Right vs XBottom",bin,-0.3,0.3,bin,0,2500},"xbottom","adcbr");

	c5->cd(1);
 	hADCvsXTL->Draw("COLZ");
 	c5->cd(2);
 	hADCvsXTR->Draw("COLZ");
  	c5->cd(3);
  	hADCvsXBL->Draw("COLZ");
  	c5->cd(4);
  	hADCvsXBR->Draw("COLZ");

	c5->DrawClone();

	return c5;

}

TCanvas* plotTDCAdjusted(){

	TCanvas *c2 = new TCanvas("c2", "c2", 150,150,600,600);
	c2->Divide(2,2, 0.01, 0.01, 0);

	auto hTDCTL = v[1].Histo1D({"h1","TDC Top Left",bin,1900,2100},"tdctl");
	auto hTDCTR = v[1].Histo1D({"h2","TDC Top Right",bin,1900,2100},"tdctr");
	auto hTDCBL = v[1].Histo1D({"h3","TDC Bottom Left",bin,1900,2100},"tdcbl");
	auto hTDCBR = v[1].Histo1D({"h4","TDC Bottom Right",bin,1900,2100},"tdcbr");

	c2->cd(1);
 	hTDCTL->Draw();
 	c2->cd(2);
 	hTDCTR->Draw();
  	c2->cd(3);
  	hTDCBL->Draw();
  	c2->cd(4);
  	hTDCBR->Draw();

	c2->DrawClone();

	return c2;

}

TCanvas* plotAngles(){

	TCanvas *c7 = new TCanvas("c7", "c7", 175,175,600,600);
	c7->Divide(2,2, 0.01, 0.01, 0);

	auto hXTop = v[1].Histo1D({"h1","Top Position",bin,-0.3,0.3},"xtop");
	auto hXBottom = v[1].Histo1D({"h2","Bottom Position",bin,-0.3,0.3},"xbottom");
	auto hTheta = v[1].Histo1D({"h3","Angle 1",bin,-100,100},"theta");
	auto hTheta2 = v[1].Histo1D({"h4","Angle 2",bin,-100,100},"theta2");

	c7->cd(1);
 	hXTop->Draw();
 	c7->cd(2);
 	hXBottom->Draw();
  	c7->cd(3);
  	hTheta->Draw();
  	c7->cd(4);
  	hTheta2->Draw();

	c7->DrawClone();

	return c7;

}

TCanvas* plotADCSums(){

	TCanvas *c8 = new TCanvas("c8", "c8", 200,200,600,600);
	c8->Divide(1,2, 0.01, 0.01, 0);

	auto hADCT = v[1].Histo1D({"h1","EDep Top",bin,0,6000},"etop");
	auto hADCB = v[1].Histo1D({"h4","EDep Bottom",bin,0,6000},"ebottom");

	c8->cd(1);
 	hADCT->Draw();
 	c8->cd(2);
 	hADCB->Draw();

	c8->DrawClone();

	return c8;

}

TCanvas* plotADCRatio(){

	TCanvas *c9 = new TCanvas("c9", "c9", 225,225,600,600);
	c9->Divide(1,2, 0.01, 0.01, 0);

	auto hADCT = v[1].Histo1D({"h1","EDep (Top: Blue, Bottom: Red)",bin,0,6000},"etop");
	auto hADCB = v[1].Histo1D({"h4","EDep Bottom",bin,0,6000},"ebottom");
	auto hADCR = v[1].Histo1D({"h4","EDep Ratio",bin,0,3.0},"eratio");

	c9->cd(1);
 	hADCT->Draw();
	hADCB->SetLineColor(kRed);
	hADCB->Draw("SAME");
 	c9->cd(2);
 	hADCR->Draw();

	c9->DrawClone();

	return c9;

}

TCanvas* plotADCTheta(){

	TCanvas *c10 = new TCanvas("c10", "c10", 250,250,600,600);
	c10->Divide(2,2, 0.01, 0.01, 0);

	auto hADCTopvsTheta = v[1].Histo2D({"h1","ADC Top vs Theta",bin,-60,60,bin,0,4500},"theta","etop");
	auto hADCBottomvsTheta = v[1].Histo2D({"h2","ADC Bottom vs Theta",bin,-60,60,bin,0,4500},"theta","ebottom");

	c10->cd(1);
 	hADCTopvsTheta->Draw("COLZ");
 	c10->cd(2);
 	hADCBottomvsTheta->Draw("COLZ");

  	c10->cd(3);
        TF1 *myLeBronFit = new TF1("myLeBronFit","[0]*(1.0+[1]*cos(x*3.14159/180.0))",-60.0,60.0);
        TProfile *prof = hADCTopvsTheta->ProfileX();
        prof->Fit("myLeBronFit","QR");
        prof->Draw();

  	c10->cd(4);
        TProfile *prof2 = hADCBottomvsTheta->ProfileX();
        prof2->Fit("myLeBronFit","QR");
        prof2->Draw();

	c10->DrawClone();

	return c10;

}

TCanvas* plotADCAdjusted(){

	TCanvas *c4 = new TCanvas("c4", "c4", 300,300,600,600);
	c4->Divide(2,2, 0.01, 0.01, 0);

	auto hADCTL = v[1].Histo1D({"h1","ADC Top Left",bin,0,3000},"adctl");
	auto hADCTR = v[1].Histo1D({"h2","ADC Top Right",bin,0,3000},"adctr");
	auto hADCBL = v[1].Histo1D({"h3","ADC Bottom Left",bin,0,3000},"adcbl");
	auto hADCBR = v[1].Histo1D({"h4","ADC Bottom Right",bin,0,3000},"adcbr");

	c4->cd(1);
 	hADCTL->Draw();
 	c4->cd(2);
 	hADCTR->Draw();
  	c4->cd(3);
  	hADCBL->Draw();
  	c4->cd(4);
  	hADCBR->Draw();

	c4->DrawClone();

	return c4;

}

void rdfluterplots(int run_number = 42) {

	cout << run_number << endl;
        auto fileName = "rootfiles/test"+std::to_string(run_number)+".root";
        auto treeName = "tdata";
	cout << fileName << " " << treeName << endl;

        //TFile* f = new TFile((TString)fileName,"READ");
        //TTree* t = (TTree*)f->Get(treeName);

	ROOT::EnableImplicitMT();
        ROOT::RDataFrame d(treeName,fileName);
	cout << "Opened RDataFrame" << endl;

	auto fdf = d.Define("trigger","getTrigger(&tdc[0],&adc[0])")
	;
	cout << "Defined addtional variables" << endl;

	auto entries = d.Count();
	cout << *entries << " entries in Tree with no filter" << endl;
	
	auto triggers = fdf.Filter("trigger==true").Count();
	cout << *triggers << " entries passed Main trigger" << endl;

	auto fdft = fdf.Filter("trigger==true");
	v.push_back(fdft); 
	
	calculateTOffsets();
	calculateADCFactors();

	auto fdftt = fdft.Define("tdctl","getTdcTL(trigger,&tdc[0])")
			 .Define("tdctr","getTdcTR(trigger,&tdc[0])")
			 .Define("tdcbr","getTdcBR(trigger,&tdc[0])")
			 .Define("tdcbl","getTdcBL(trigger,&tdc[0])")
			 .Define("adctl","getAdcTL(trigger,&adc[0])")
			 .Define("adctr","getAdcTR(trigger,&adc[0])")
			 .Define("adcbl","getAdcBL(trigger,&adc[0])")
			 .Define("adcbr","getAdcBR(trigger,&adc[0])")
			 .Define("xtop","getXTop(trigger,tdctl,tdctr,tdcbl,tdcbr)")
			 .Define("xbottom","getXBottom(trigger,tdctl,tdctr,tdcbl,tdcbr)")
			 .Define("theta","getTheta(trigger,xtop,xbottom)")
			 .Define("theta2","getTheta2(trigger,xtop,xbottom)")
			 .Define("etop","getETop(trigger,adctl,adctr)")
			 .Define("ebottom","getEBottom(trigger,adcbl,adcbr)")
			 .Define("eratio","getERatio(trigger,etop,ebottom)")
	;
	v.push_back(fdftt);
}

