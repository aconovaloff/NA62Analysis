#ifndef RUNCALCULATIONS_HH
#define RUNCALCULATIONS_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include "MCSimple.hh"
#include "DetectorAcceptance.hh"
#include <TCanvas.h>
#include "OneTrack.hh"
#include "CalorimeterAnalysis.hh"
#include "PhotonAnalysis.hh"
#include "TrackAnalysis.hh"
#include "StrawAnalysis.hh"

class TH1I;
class TH2F;
class TGraph;
class TTree;


class RunCalculations : public NA62Analysis::Analyzer
{
	public:
		RunCalculations(NA62Analysis::Core::BaseAnalysis *ba);
		~RunCalculations();
		void InitHist();
		void InitOutput();
		void DefineMCSimple();
		void Process(int iEvent);
		void StartOfBurstUser();
		void EndOfBurstUser();
		void StartOfRunUser();
		void EndOfRunUser();
		void PostProcess();
		void DrawPlot();
		Int_t NEvents;
		Int_t TypeCounter; 
		Int_t burst_number; 
	protected:
	TrackAnalysis* Tracka;
	StrawAnalysis* Strawa; 
	TString acc_mmsq[4];
	TString mc_pi[4];
	TString mc_p[4];  
	TString mc_ps[4];
	TString mc_chambers[4];
	TString mc_chambersp[4];
};
#endif
