#include "TH1F.h"
#include "TPad.h"
#include "TROOT.h"

void test_gROOT() {
	TH1F* h = (TH1F*)gROOT->FindObject("h_gaus");

	if (h) {
		h->Reset();
	} else {
		h = new TH1F("h_gaus","Une gaussienne", 100, 0.0, 100.0);
	}

	h->Draw();

	for (Double_t x=0; x<100; x++) {
		Double_t f = 20.0*exp(-pow(x-50.,2)/2./pow(10,2));
		h->Fill(x,f);
		gPad->Modified();
		gPad->Update();
	}
}
