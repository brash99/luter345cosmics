//===========Program to plot TDC and ADC data from root trees by Davison, Hesse, Losasda, Kerver Brash==========
#include <iostream>
#include <fstream>

//Global Variables
const int tbin = 1000; 
const int abin = 500; 
const Int_t Nadc = 16;
const Int_t Ntdc = 16;
const Int_t Ntdc2 = 4;
//const Int_t pedrun = 89;
//const Int_t bl=Ntdc-16+0;//Start at channel 0
//const Int_t br=Ntdc-16+1;
//const Int_t tl=Ntdc-16+2;
//const Int_t tr=Ntdc-16+3;



//=============================================Method starts here for plotting===================================

void luterpedestals(Int_t nrun) {
	FILE *f = fopen(Form("./pedestalfiles/pedestalrun%d.dat",nrun),"w");	
	TCanvas *adc_canvas2 = new TCanvas("adc_canvas2","ADC pedestals",50,50,600,600);
        adc_canvas2->Divide(2,2);
	Char_t tdcnames[][Ntdc]={"Bottom Left","Bottom Right","Top Left","Top Right","4","5","6","7","8","9","10","11","12","13","14","15"};	
	Char_t adcnames[][Nadc]={"Bottom Left","Bottom Right","Top Left","Top Right","4","5","6","7","8","9","10","11","12","13","14","15"};

//Set correction values
    double adcpedestal[Nadc];
	
//Read files and tree branches 
	gStyle->SetOptStat(1);
	TFile *froot =  new TFile(Form("./rootfiles/test%d.root",nrun));
	TTree *troot = (TTree*)froot->Get("tdata");
	Int_t tdc[Ntdc]; 
	Int_t adc[Nadc]; 
	troot->SetBranchAddress("tdc",&tdc);
	troot->SetBranchAddress("adc",&adc);

//Create Histograms
	TH1F *htdcraw[Ntdc], *hadcraw[Nadc];
	
	for ( int i = 0; i < Ntdc ; i++) {
		htdcraw[i] = new TH1F(Form("htdcraw%02d", i),Form("%s   raw tdc",tdcnames[i]),tbin,1500,2500);
	}
    
	for ( int i = 0; i < Nadc ; i++) {
        hadcraw[i] = new TH1F(Form("hadcraw%02d", i),Form("%s   raw adc",adcnames[i]),abin,0,500);
    }
    
	
//================================ Start looping through the number of events ==========================
    int nevents_in_file = (int)troot->GetEntries();
	//cout << " Nentries = "<< nevents_in_file << endl;
//===================================FILL RAW HISTOS==============================

    for (int ie = 0; ie < nevents_in_file; ie++) {
		troot->GetEntry(ie);
		for( int i = 0; i < Nadc ; i++) {//Start ADC filling loop
            		if(adc[i] > 0) {hadcraw[i]->Fill(adc[i]);}
		}
		for( int i = 0; i < Ntdc ; i++) {//Start TDC filling loop
            		if(tdc[i] > 0) {htdcraw[i]->Fill(tdc[i]);}
		}
	}
	for( int i = 0; i < 4 ; i++) {//Start ADC filling loop
		adc_canvas2->cd(i+1);
		hadcraw[i]->Draw();
		cout << "Mean = " << hadcraw[i]->GetMean() << "  RMS = " << hadcraw[i]->GetRMS() << endl;
		adcpedestal[i] = hadcraw[i]->GetMean()+3.0*hadcraw[i]->GetRMS();
//		hadcraw[i]->Fit("gaus","Q",0,500);
//		adcpedestal[i] = (hadcraw[i]->GetFunction("gaus")->GetParameter(1))+(2*hadcraw[i]->GetFunction("gaus")->GetParameter(2));
//		printf("%f\n",adcpedestal[i]);
		fprintf(f,"%f\n",adcpedestal[i]);	

	}
	fclose(f);
//	adc_canvas2->Close();
	
}
