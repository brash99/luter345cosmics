////////////////////////////////////////////////////////
//
//    EEL VME test stand simple decoding macro
//    to be used with Analyzer v1.5 or more recent
//    For 1 V792 and one V775
//
//    Dec 23rd 2008 
//    
//    A. Camsonne
//    camsonne@jlab.org
//
/////////////////////////////////////////////////////////


//adjust according to the analysis computer used
#define DATAFILENAME "/home/coda/CODA/data/test_%d.dat"
#include<iostream>

void decodeLuter(int run)
{
  THaCodaFile * coda = new THaCodaFile(Form(DATAFILENAME,run));

  // variable declaration
  const int Nchan=16;
  int * buffer; // coda buffer where the event will be stored
  int adc[Nchan]; // data array for the adc
  int tdc[Nchan]; // data array for the tddc
  int offset=1; // offset where the first data word appears after the header
  int chan,val;// decoded channel and corresponding value variables
  int index=0; // index used to go through the buffer
  int evcount=0;// number of events processed
  int size; // number of channels to be read
  if (coda != NULL) {                               
    TFile * out = new TFile (Form("./rootfiles/test%d.root",run),"RECREATE"); // creation of the output file                                         
    TTree * tdata = new TTree("tdata", Form("Run %d",run)); // creation of the tree
    tdata->Branch("adc",&adc,"adc[16]/I"); // creation of the branch to hold the adc data
    tdata->Branch("tdc",&tdc,"tdc[16]/I"); // creation of the branch to hold the tdc data
    THaCodaDecoder * ev = new THaCodaDecoder;
    while (coda->codaRead()==0 )
    {
	// initialisation of variables used for each events
	index=0;
	evcount++;
	for ( int i = 0 ;i <Nchan ;i++)
	  {
	    adc[i]=0;
	    tdc[i]=0;
	  }
	buffer = coda->getEvBuffer(); // load the data of th event buffer
	ev->LoadEvent(buffer);
	int nevt=ev->GetEvLength();

	while (index<nevt && ((buffer[index]&0xff000000)!=0xfa000000))//TONY- may change
	  {
	    //cout<<hex<<buffer[index]<<endl;
 	    index++; // look for header of the adc
	  }//TONY: for our setup index==1 at this point
	size = (0x00001f00&buffer[index])>>8;//TONY must be 1f to get full size
//	size = (0xf000&buffer[index+1])>>8  ;
	//if(evcount<20){cout<<" Index = " << index << "  Size :"<<size<<endl;}
	 // decode the channels read 
	 for ( int i = 0 ; i < size ; i++)
	   {
	     chan=(((buffer[i+index+offset])&0x1e0000)>>17);
	     val=((buffer[i+index+offset])&0x0fff);
	     // val=((buffer[i+index+offset])&0x1fff);//TONY this pulls overflowbit
	     //if(evcount<20){cout<< i << " " << offset << " " << evcount<<" "<<chan<<" "<<val<<endl;}
	     adc[chan]=val;
	   }
	index++;//TONY: increment index past fa001000 for qdc 	
	 while (index<nevt && ((buffer[index]&0xff000000)!=0xfa000000)) //this moves index to word count for tdc
          {
            index++;//look for header of the tdc
	    //if (evcount <=10){cout<<dec<<index<<"l2 "<<hex<<buffer[index]<<" "<<buffer[index+1]<<endl;}
          }
	    //if (evcount <=10){cout<<dec<<index<<" "<<hex<<buffer[index]<<" "<<buffer[index+1]<<endl;}
	size = (0x00001f00&buffer[index])>>8;
//   	size = buffer[index+1] -2 ;
	//if (evcount <= 50 && size>0){cout<<"TDC size: "<< size<<endl;}
 // decode the tdc channels read 
   	for ( int i = 0 ; i < size ; i++){
	     		chan=(((buffer[i+index+offset])&0x1e0000)>>17);
	     		val=((buffer[i+index+offset])&0xfff);
	     		//if (evcount <= 50){cout<<"TDC chan: "<< chan<<"     TDC value: "<<val<<endl;}
	     		tdc[chan]=val;
     	}		

	
   	tdata->Fill();// fill the tree with the new eventB
   	if ( (evcount%100)==0)
     	{
       		cout<<evcount<<endl; // display the number of events processed every 100 events
     	}
}
    tdata->Write(); // end of the coda file ... write the tree
    out->Close(); // and close the output file

  }
 
}
