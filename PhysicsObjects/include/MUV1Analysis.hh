#ifndef MUV1Analysis_h
#define MUV1Analysis_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TRecoMUV1Event.hh"
#include "UserMethods.hh"
class MUV1Analysis
{
  public :
        MUV1Analysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoMUV1Event*, Double_t, Double_t, Double_t, bool);
	void BookHistos_MatchedTrack();
	void BookHistos_MuonID(); 
        void BookHistos_Set();
	void SaveAllPlots(); 
        Int_t MuonID();
	Int_t MatchedTrack(TVector3);
	Double_t TotalEnergyHits(); 
        ~MUV1Analysis();
	UserMethods* fUserMethods;
	Double_t MatchedClusterE;
	Double_t GetMatchedClusterE(){return MatchedClusterE;};
	Int_t MUV1ID;
        Int_t GetMUV1ID(){return MUV1ID;};
	Double_t Distance; 
	Double_t GetDistance(){return Distance;}; 

  private:
  private:
        TRecoMUV1Event *event;
        Double_t reftime;
        Double_t tLow;
        Double_t tHi;
	bool MCflag;
};

#endif
