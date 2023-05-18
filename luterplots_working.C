//===========Program to plot TDC and ADC data from root trees by Davison, Hesse, Losasda, Kerver Brash==========
#include <iostream>
#include <fstream>

//Global Variables
const int bin = 50; 
const Int_t Nadc = 16;
const Int_t Ntdc = 16;
const Int_t ped = 80;
const Int_t bl=Ntdc-16+0;//Start at channel 0
const Int_t br=Ntdc-16+1;
const Int_t tl=Ntdc-16+2;
const Int_t tr=Ntdc-16+3;

//=============================================Method starts here for plotting===================================

void luterplots(Int_t nrun) {
	
//Read files and tree branches 
	gStyle->SetOptStat(1);
	TFile *froot =  new TFile(Form("./rootfiles/test%d.root",nrun));
	TTree *troot = (TTree*)froot->Get("tdata");
	Int_t tdc[Ntdc]; 
	Int_t adc[Nadc]; 
	troot->SetBranchAddress("tdc",&tdc);
	troot->SetBranchAddress("adc",&adc);

	

//Create Histograms
	TH1F *htdcraw[Ntdc], *hadcraw[Nadc], *hadccut[Nadc];
	htpos = new TH1F("htpos","Top Position",1000,-0.50,0.50); //Histogram for top scintillator
	hbpos = new TH1F("hbpos","Bottom Position",1000,-0.50,0.50); //Histogram for bottom scintillator 
	htheta = new TH1F("htheta","Angle (Degrees)",bin,-100.5,100.5); //Histogram for incidence angel 

	for ( int i = 0; i < Ntdc ; i++) {
		htdcraw[i] = new TH1F(Form("htdcraw%02d", i),Form("%02d   raw tdc",i),1000,1500,2500);
	}
    
	for ( int i = 0; i < Nadc ; i++) {
        hadcraw[i] = new TH1F(Form("hadcraw%02d", i),Form("%02d   raw adc",i),1000,0,2500);
        hadccut[i] = new TH1F(Form("hadccut%02d", i),Form("%02d   cut adc",i),1000,0,2500);
        }

	
    
//================================ Start looping through the number of events ==========================
    
    int nevents_in_file = (int)troot->GetEntries();
    


	//cout << " Nentries = "<< nevents_in_file << endl;
    
    
	for (int ie = 0; ie < nevents_in_file; ie++) {
        
		troot->GetEntry(ie);
        
		for ( int i = 0; i < Ntdc ; i++) {//Start TDC filling loop
            
            htdcraw[i]->Fill(tdc[i]);
            
		}//End TDC for loop
        
        
    for ( int i = 0; i < Nadc ; i++) {//Start ADC filling loop
            
            hadcraw[i]->Fill(adc[i]);
            
            if(adc[i] > ped){
                
                 hadccut[i]->Fill(adc[i]);
                 hadccut[i]->GetMean();
                            }// end if for ped cut

                }//End ADC for loop

//
		Double_t tbl = tdc[bl]+164.0;
		Double_t tbr = tdc[br]+43.0;
		Double_t ttl = tdc[tl]+111.0;
		Double_t ttr = tdc[tr]-13.0;
		Double_t rtod = 180.0/3.14159265;
		Double_t nscint = 1.50; // index of refraction of scintillator
		Double_t vn = 2.997E08/nscint;
		Double_t dscint = 0.105; // distance between scintillators in metres
		Double_t t_fullscale = 140.0E-09; // full scale TDC range in seconds
		Double_t t_convert=t_fullscale/4096.0;
		Double_t tdiff = (ttl-ttr)/2.0;
		Double_t bdiff = (tbl-tbr)/2.0;
		Double_t xtop = tdiff*t_convert*vn;
		Double_t xbottom = bdiff*t_convert*vn;
		Double_t theta = rtod*atan((xbottom-xtop)/dscint);

		htpos->Fill(xtop);
		hbpos->Fill(xbottom);
	        htheta->Fill(theta); 

		
	} //end of loop over events 

	/* TFile *meanfile =  new TFile(("./analyzer/adcmeans.txt","CREATE"));

        Double_t means[4];
	
	for(int j=0;j<4;j++)
	{
            means[j]=hadccut[j]->GetMean();  
		meanfile<<"jhdfskajhfdj"<<endl;
        }*/

	

 	TCanvas *traw = new TCanvas("traw","Raw TDC spectra",100,100,600,600);
   	traw->Divide(2,2);
   	for ( int i = bl; i <= tr ; i++) {
      		traw->cd(i-bl+1);
      		htdcraw[i]->Draw();
	}

        TCanvas *araw = new TCanvas("araw","Raw ADC spectra",125,125,600,600);
        araw->Divide(2,2);
        for ( int i = 0; i <= 3 ; i++) {
                araw->cd(i+1);
                hadcraw[i]->Draw();
        }
/* 	TCanvas *traw2 = new TCanvas("traw2","X Positions",150,150,600,900);
	traw2->Divide(1,3);
	gStyle->SetOptFit(1);
	traw2->cd(1);
	htpos->Draw();
	traw2->cd(2);
	hbpos->Draw();
	traw2->cd(3);
        TF1 *myfit = new TF1("myfit","[0]*pow(cos(x/180.0*3.14159),2)",-70,70);
        myfit->SetParameter(0,140);
	htheta->Fit("myfit","R");
*/
}
