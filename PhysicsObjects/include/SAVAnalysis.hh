#ifndef SAVANALYSIS_h
#define SAVANALYSIS_h

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "UserMethods.hh"
#include "TRecoSAVEvent.hh"

class SAVAnalysis 
{
  public :
 	SAVAnalysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoSAVEvent*, Double_t, Double_t, Double_t, bool); 
     	Int_t Veto();
        void BookHistos_Set();
	void SaveAllPlots();  
	UserMethods* fUserMethods;
  	~SAVAnalysis();
  private:
    	TRecoSAVEvent *event;
    	Double_t reftime;
    	Double_t tLow;
    	Double_t tHi;  
	bool MCflag;
};

#endif
