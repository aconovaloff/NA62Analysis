#ifndef SACANALYSIS_h
#define SACANALYSIS_h

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "UserMethods.hh"
#include "TRecoSACEvent.hh"

class SACAnalysis 
{
  public :
 	SACAnalysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoSACEvent*, Double_t, Double_t, Double_t, bool); 
     	Int_t Veto();
        void BookHistos_Set();
	void SaveAllPlots();  
	UserMethods* fUserMethods;
  	~SACAnalysis();
  private:
    	TRecoSACEvent *event;
    	Double_t reftime;
    	Double_t tLow;
    	Double_t tHi;  
	bool MCflag; 

};

#endif
