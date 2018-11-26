#ifndef MUV3Analysis_h
#define MUV3Analysis_h

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "UserMethods.hh"
//#include "UserMethods.hh"
#include "TRecoMUV3Event.hh"
//#include "TRecoSpectrometerCandidate.hh"
//#include "CHODCandidate.hh"

//class AnalysisTools;

class MUV3Analysis 
{
  public :
 	MUV3Analysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoMUV3Event*, Double_t, Double_t, Double_t, bool); 
     	Int_t Veto();
	Int_t VetoL(); 
	Int_t Coincidental();
	Int_t MuonID(TVector3);
	Double_t GetMinDistance();
	void BookHistos_MuonID();
        void BookHistos_Set();
	void SaveAllPlots();  
	UserMethods* fUserMethods;
  	~MUV3Analysis();
  private:
//  AnalysisTools *fTools;
//  UserMethods *fUserMethods;
  //  MUV3Candidate *fMUV3Candidate;

  private:
    	TRecoMUV3Event *event;
    	Double_t reftime;
    	Double_t tLow;
    	Double_t tHi;
	Double_t fMUV3_minD;  
	bool MCflag;
//  Int_t MUV3Veto(TRecoMUV3Event*,Double_t,Double_t,Double_t);
//  Int_t MakeExpectedCandidate(Double_t,Double_t,Double_t,Double_t*);
//  Int_t Quadrant(Int_t);
};

#endif
