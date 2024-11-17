using namespace std;
using RNode = ROOT::RDF::RNode;
std::vector<RNode> v;

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
//const Double_t adjadcto=1400.0;//value to ADJust ADC TO
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

float getTheta(bool trigger, float xtop, float xbottom) {

    double rtod=180.0/3.14159265;
    if (trigger) {
        rnd = r.Gaus(0.0,1.5);
        return rtod*atan((xbottom-xtop)/dscint)+rnd;
    }

    return -1000;

}

float getXTop(bool trigger, double ttl, double ttr, double tbl, double tbr) {

        //std::vector<float> v;

        if(trigger) {
                return (ttl-ttr)/2.0*t_convert*vn;
        }

        return -1000;
}

float getXMeanTop(bool trigger, double ttl, double ttr, double tbl, double tbr) {

        //std::vector<float> v;

        if(trigger) {
                return (ttl+ttr-4000.0)/2.0*t_convert*vn;
        }

        return -1000;
}

float getXMean(bool trigger, double xmeantop, double xmeanbottom) {

        //std::vector<float> v;

        if(trigger) {
                return (xmeantop-xmeanbottom);
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

float getXMeanBottom(bool trigger, double ttl, double ttr, double tbl, double tbr) {

        //std::vector<float> v;

        if(trigger) {
                return (tbl+tbr-4000.0)/2.0*t_convert*vn;
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
