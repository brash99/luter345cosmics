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
//  adjusted for tdc only by ARL 
/////////////////////////////////////////////////////////


//adjust according to the analysis computer used
#define DATAFILENAME "/home/coda/LuterCODA/data/test_%d.dat"

void decodeLuter2(int run)
{
  THaCodaFile * coda = new THaCodaFile(Form(DATAFILENAME,run));

  // variable declaration
  const int Nchan=16;
  int * buffer; // coda buffer where the event will be stored
  int tdc[Nchan]; // data array for the tddc
  int offset=1; // offset where the first data word appears after the header
  int chan,val;// decoded channel and corresponding value variables
  int index=0; // index used to go through the buffer
  int evcount=0;// number of events processed
  int size; // number of channels to be read
  if (coda != NULL) {                               
    TFile * out = new TFile (Form("./rootfiles/test%d.root",run),"RECREATE"); // creation of the output file                                         
    TTree * tdata = new TTree("tdata", Form("Run %d",run)); // creation of the tree
    tdata->Branch("tdc",&tdc,"tdc[16]/I"); // creation of the branch to hold the tdc data
    THaCodaDecoder * ev = new THaCodaDecoder;
    while (coda->codaRead()==0)
    {
// initialisation of variables used for each events
	index=0; 
	evcount++;
	for ( int i = 0 ;i < Nchan ;i++)
	  {
	    tdc[i]=0;
	  }
	buffer = coda->getEvBuffer(); // load the data of the event buffer
	ev->LoadEvent(buffer); //get event from buffer
	int nevt=ev->GetEvLength();
	while (index<nevt && ((buffer[index]&0xff000000)!=0xfa000000)) //brash changed this from d00dd00d
          {
            index++;//look for header of the tdc
	    if (evcount <=3){cout<<dec<<index<<" l2 "<<hex<<buffer[index]<<" "<<buffer[index+1]<<endl;}
          }
	    if (evcount <=3){cout<<dec<<index<<" out of index loop "<<hex<<buffer[index]<<" "<<buffer[index+1]<<endl;}
	size = (0x0f00&buffer[index])>>8  ;
	//size = buffer[index+1] -2 ;
	//if (evcount <= 50 && size>0){cout<<"TDC size: "<< size<<endl;}
// decode the tdc channels read 
   	for ( int i = 0 ; i < size ; i++){
			/* V775 = 0x1f0000 >> 16 */
			/* V775N = 0x1e0000 >> 17 */
	     		chan=(((buffer[i+index+offset])&0x1e0000)>>17);
	     		val=((buffer[i+index+offset])&0xfff);
	     		if (evcount <= 3){cout<<"TDC chan: "<< chan<<"     TDC value: "<<val<<endl;}
	     		tdc[chan]=val;cout<<"TONY:  Eventcount:"<<evcount<<" channel:"<<chan<<" value:"<<val<<endl;
     	}		

	
   	tdata->Fill();// fill the tree with the new eventB
   	if ( (evcount%100)==0)
     	{
       		cout<<evcount<<endl; // display the number of events processed every 100 events
     	}
    }//closes  while(coda->codaRead()==0)

    tdata->Write(); // end of the coda file ... write the tree
    out->Close(); // and close the output file

  }
 
}
