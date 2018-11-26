#ifndef MUV2Analysis_h
#define MUV2Analysis_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TRecoMUV2Event.hh"
#include "UserMethods.hh"
class MUV2Analysis
{
  public :
        MUV2Analysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoMUV2Event*, Double_t, Double_t, Double_t, bool);
        Int_t MuonID();
	Int_t MatchedTrack(TVector3);	
	void SaveAllPlots();
	void BookHistos_MuonID();
	void BookHistos_Set();
	void BookHistos_MatchedTrack(); 
	Double_t TotalEnergy();
	Double_t TotalEnergyHits();  
	Double_t MatchedClusterE;
	Double_t GetMatchedClusterE(){return MatchedClusterE;}; 
	UserMethods* fUserMethods;
	Int_t MUV2ID;
        Int_t GetMUV2ID(){return MUV2ID;};
	Double_t Distance;
	Double_t GetDistance(){return Distance;}; 
        ~MUV2Analysis();
  private:
        TRecoMUV2Event *event;
        Double_t reftime;
        Double_t tLow;
        Double_t tHi;
	bool MCflag; 
};	

#endif
