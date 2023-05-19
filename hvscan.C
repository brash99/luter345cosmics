//===========Program to plot TDC and ADC data from root trees by Davison, Hesse, Losasda, Kerver Brash==========
#include <iostream>
#include <fstream>

//Global Variables
bool remakepedfile= false;
const int bin = 100; 
const Int_t Nadc = 16;
const Int_t Ntdc = 16;
const Int_t pedrun = 89;
const Double_t adjadcto=2000.0;//value to ADJust ADC TO
const Int_t numscan = 14;



//=============================================Method starts here for plotting===================================
//91 93 94 96 98 99 100 101 102
void hvscan() {
const Int_t runlist[numscan] = {91,  109, 106, 108, 104, 110 ,103, 111, 100, 117, 118, 114, 115, 116 };
const Int_t *hv[numscan]={      1850,1800,1750,1700,1650,1625,1600,1575,1550,1525,1500,1475,1450,1425};
TCanvas *adc_canvas1[numscan];
Double_t *adc_means[numscan][Nadc];
FILE *f = fopen("./scan_means.dat","w");	

for(int irun=0;irun<numscan;irun++){
int nrun = runlist[irun];


	char *adcnames[Nadc]={"Bottom Left","Bottom Right","Top Left","Top Right","4","5","6","7","8","9","10","11","12","13","14","15"};


//Set correction values

    Double_t ped[Nadc];
    Double_t gain[Nadc];

	
//Read files and tree branches 
	gStyle->SetOptStat(1);
	TFile *froot =  new TFile(Form("./rootfiles/test%d.root",nrun));
	TTree *troot = (TTree*)froot->Get("tdata");
	Int_t tdc[Ntdc]; 
	Int_t adc[Nadc]; 
	troot->SetBranchAddress("tdc",&tdc);
	troot->SetBranchAddress("adc",&adc);
	const int nevents_in_file = (int)troot->GetEntries();

//Create Histograms
	TH1F *hadcraw[Nadc], *hadccut[Nadc], *hadcadjusted[Nadc];
    
	for ( int i = 0; i < Nadc ; i++) {
        	hadcraw[i] = new TH1F(Form("hadcraw%02d", i),Form("%s   raw adc",adcnames[i]),bin,0,5000);
        	hadccut[i] = new TH1F(Form("hadccut%02d", i),Form("%s  cut adc",adcnames[i]),bin,0,5000);
        	//hadcadjusted[i] = new TH1F(Form("hadcadjusted%02d", i),Form("%s   adjusted adc",adcnames[i]),bin,0,8500);
        }
    
	
//=====================================GET PED for ADC==========================
	if(remakepedfile){gROOT->ProcessLine(Form(".x luterpedestals.C(%d)",pedrun));cout<<"pedestal file made"<<endl;}
	FILE *adcpeds = fopen(Form("./pedestalfiles/pedestalrun%d.dat",pedrun),"r");
	for( int i = 0; i < Nadc ; i++) {//Start ADC filling loop
		fscanf(adcpeds,"%lf\n",&ped[i]);
//		printf("%lf\n",ped[i]);
	}

//=====================================GET ADJ FACTORS for ADC&TDC==========================
    for (int ie = 0; ie < nevents_in_file; ie++) {
		troot->GetEntry(ie);
		if(ie%130000==0.0&&ie!=0){cout<<"Progress: "<<((double)ie)/((double)nevents_in_file)*100<<"%"<<endl;}
		
		for( int i = 0; i < Nadc ; i++) {//Start ADC-ADJUSTED filling loop
			hadcraw[i]->Fill(adc[i]);
			hadccut[i]->Fill(adc[i]-ped[i]);
		}//End ADC for loop
	}



//==============LOOP THROUGH EVENTS IN FILE AND DO MATHS=============================
/*	for (int ie = 0; ie < nevents_in_file; ie++) {
		troot->GetEntry(ie);  
	//=========================================ADC===============================
        	if(ie%130000==0.0&&ie!=0){cout<<"Progress: "<<((double)ie)/((double)nevents_in_file)*100<<"%"<<endl;}
		for( int i = 0; i < Nadc ; i++) {//Start ADC-ADJUSTED filling loop
			gain[i]=adjadcto/(hadccut[i]->GetMean());	 
			hadcadjusted[i]->Fill(gain[i]*(adc[i]-ped[i]));
		}//End ADC for loop	    	
	} //End of loop over events

//*/   
//=======================================Canvases=============================

//Draw ADC raw data
        adc_canvas1[irun] = new TCanvas(Form("adc_canvas1_%d",irun),Form("%d",hv[irun]),125+irun*25,125+irun*25,600,600);
	adc_canvas1[irun]->Divide(2,2);
        for ( int i = 0; i <= 3 ; i++) {
                adc_canvas1[irun]->cd(i+1);
		adc_means[irun][i] = hadcraw[i]->GetMean();
		if(irun==0){
			hadcraw[i]->Draw();
		}
		else{
			hadcraw[i]->SetMarkerColor(1+irun);
			hadcraw[i]->Draw("same");

		}

        }
	

 
}//ends irun loop

for(int irun=0;irun<numscan;irun++){
	printf("%i %f %f %f %f\n",hv[irun],adc_means[irun][0],adc_means[irun][1],adc_means[irun][2],adc_means[irun][3]);
	fprintf(f,"%i %f %f %f %f\n",hv[irun],adc_means[irun][0],adc_means[irun][1],adc_means[irun][2],adc_means[irun][3]);    
} 

fclose(f);  
}
