#include<iostream>
void adcmeans()
{
	TCanvas *c1 = new TCanvas("ADC Means","ADC Means", 200, 10, 700, 500);
        c1->SetFillColor(42);
        c1->SetGrid();
     
	const int num=9;

	int voltage[num]={1200,1250,1300,1350,1400,1450,1500,1550,1600}//,1650,1700,1750,1800,1850,1900,1950,2000}

	double adc1[num]={242.5,251.2,249,259.7,282.2,308.4,349.5,397.6,456.3}
	double adc2[num]={364.1,380.8,365.2,403.7,460.4,543.6,651.9,779.3,927.7}
	double adc3[num]={293.8,310.7,309.9,348.5,416.8,505.2,614.3,746,903.6}
	double adc4[num]={291.3,317.7,322.1,364.4,422.9,493.1,596.5,713.7,850.5}

	gr1 = new TGraph(num,voltage,adc1)
        TGraph *gr2 = new TGraph(num,voltage,adc2)
        TGraph *gr3 = new TGraph(num,voltage,adc3)
        TGraph *gr4 = new TGraph(num,voltage,adc3)

	gr1->SetMarkerColor(kBlue);
	gr1->Draw("AP");

	gr2->SetMarkerColor(kGreen);
        gr2->Draw("P");

        gr3->SetMarkerColor(kOrange);
        gr3->Draw("P");

        gr4->SetMarkerColor(kRed);
        gr4->Draw("P");





}
