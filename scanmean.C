#include <string>
//#include <ostream>
#include <sstream>
#define _USE_MATH_DEFINES


//**************USER OPTIONS*******************************\\
/***/ const int m=2;    //ORDER of fit			
/***/// Option_t* opt = "pol2"; //replace # with value of m
/***/ const bool errors = false; //true if using error	
/***/ FILE *f = fopen("scan_means.dat","r");			
//********************************************************\\
//data arrays
Double_t hv[100000],pmt1[100000],pmt2[100000],pmt3[100000],pmt4[100000];
int n=0;
//methods
void fill();
void graph();

//MAIN
void scanmean(){
	fill();
	graph();
}


void graph(){
	TCanvas *c1 = new TCanvas("canvas","HW2",40,40,600,600);
	//set options
	gStyle->SetOptStat(1); //turns on(1)/off(0) stat box
	gStyle->SetOptFit(1);  //turns on(1)/off(0) fit box
	c1->ToggleEventStatus(); //turns on event status bar
	const TGraph *g = new TGraph(n,hv,pmt1);
	h = new TGraph(n,hv,pmt2);
	j = new TGraph(n,hv,pmt3);
	k = new TGraph(n,hv,pmt4);
	g->SetLineColor(2);		h->SetLineColor(4);		j->SetLineColor(6);		k->SetLineColor(8);
	g->SetMarkerColor(2);	h->SetMarkerColor(4);	j->SetMarkerColor(6);	k->SetMarkerColor(8);
	g->SetMarkerStyle(21);	h->SetMarkerStyle(21);	j->SetMarkerStyle(21);	k->SetMarkerStyle(21);
	
	TMultiGraph *gg = new TMultiGraph();
	gg->Add(g,"lp");gg->Add(h,"lp");gg->Add(j,"lp");gg->Add(k,"lp");	
	g->SetTitle("hv scan");
	gg->Draw("a");
	TLegend *l =  new TLegend(.2,1,.4,.8);
	l->AddEntry(g, "pmt1","lp");
	l->AddEntry(h, "pmt2","lp");
	l->AddEntry(j, "pmt3","lp");
	l->AddEntry(k, "pmt4","lp");
	l->Draw();
	c1->Update();
}
//read in file of unknown length
void fill(){
	int i=0;
	while(!feof(f)){
		fscanf(f,"%lf %lf %lf %lf %lf\n",&hv[i],&pmt1[i],&pmt2[i],&pmt3[i],&pmt4[i]);
//		printf("%i %lf %lf %lf %lf %lf\n",i,&hv[i],&pmt1[i],&pmt2[i],&pmt3[i],&pmt4[i]);
		i++;n++;
	}
	cout<<"file read."<<endl;
}
