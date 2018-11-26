#ifndef LKRANALYSIS_H
#define LKRANALYSIS_H
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TRecoLKrEvent.hh"
#include "UserMethods.hh"

class LKrAnalysis
{
  public :
        LKrAnalysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoLKrEvent*, Double_t, Double_t, Double_t, bool);
        Int_t MuonID(TVector3);
	Int_t SinglePionCluster(TVector3, TVector3, bool);
	Int_t MatchedPion(TVector3, TVector3); 
	Int_t MatchedMuon(TVector3); 
        Int_t SinglePhoton(TVector3, TVector3, TVector3, TLorentzVector, TLorentzVector);
	Int_t TrackCellMatching(TVector3*);
        void BookHistos_SinglePhoton();
	void BookHistos_SinglePionCluster();
        void BookHistos_TwoNeutralPions();
	void BookHistos_Set();
	void BookHistos_MatchedCluster();
	void BookHistos_MuonID();  
	Int_t MatchedCluster(TVector3);   
	void SaveAllPlots(); 
        Double_t GetMinDistance();
	Int_t OneNeutralPion(TLorentzVector, Int_t);
	 Int_t OneNeutralPion_Straws(TVector3, TVector3, TLorentzVector, Int_t);
	Int_t TwoNeutralPions(TLorentzVector, TVector3); 
	void BookHistos_OneNeutralPion(); 
	void BookHistos_OneNeutralPion_Straws(); 
	UserMethods* fUserMethods;
	TLorentzVector Pi04Mom;
	TLorentzVector GetPi0Momentum() {return Pi04Mom;}; 
	TLorentzVector PiPlus4Mom_LKr;
	TLorentzVector GetPiPlus4Mom_LKr(){return PiPlus4Mom_LKr;};
	TVector3 PiPlus3Mom_LKr;
	TVector3 GetPiPlus3Mom(){return PiPlus3Mom_LKr;};
	TVector3 DecayVert;
	TVector3 GetDecayVert(){return DecayVert;} 
	bool nPion_flag;
	bool GetnPionFlag(){return nPion_flag;}; 
	Double_t PiPlus_mmsq;
	Double_t GetPiPlus_mmsq(){return PiPlus_mmsq;}; 
        Double_t clusterE;
        Double_t GetclusterE(){return clusterE;};
	TVector3 DecayVertexAve;
	TVector3 GetDecayVertexAve(){return DecayVertexAve;};
        ~LKrAnalysis();
	void Test(TRecoLKrEvent*&);
        Double_t mindist;
        Double_t Getmindist(){return mindist;};
	Int_t MatchedPion_unmatched;
	Int_t GetMatchedPion_unmatched(){return MatchedPion_unmatched;}; 
	Int_t SinglePionCluster_unmatched;
        Int_t GetSinglePionCluster_unmatched(){return SinglePionCluster_unmatched;};
	Int_t SinglePionCluster_minID;
        Int_t GetSinglePionCluster_minID(){return SinglePionCluster_minID;};
	Int_t MatchedPion_minID;
        Int_t GetMatchedPion_minID(){return MatchedPion_minID;};
	Int_t OneNeutralPion_Straws_minID;
        Int_t GetOneNeutralPion_Straws_minID(){return OneNeutralPion_Straws_minID;}; 
        Int_t MatchedMuon_unmatched;
        Int_t GetMatchedMuon_unmatched(){return MatchedMuon_unmatched;};
	TRecoLKrCandidate *fLKrCellCandidate;
	TRecoLKrCandidate* GetLKrCellCandidate(){return fLKrCellCandidate;}; 
	Double_t RMS;
        Double_t GetRMS(){return RMS;};
	Double_t MatchedCluster_RMS;
        Double_t GetMatchedCluster_RMS(){return MatchedCluster_RMS;};
	Double_t OneNeutralPion_Straws_RMS;
        Double_t GetOneNeutralPion_Straws_RMS(){return OneNeutralPion_Straws_RMS;};
	Int_t MatchedCluster_minID;
        Int_t GetMatchedCluster_minID(){return MatchedCluster_minID;};
	Int_t MatchedMuon_minID;
        Int_t GetMatchedMuon_minID(){return MatchedMuon_minID;};
  private:
        TRecoLKrEvent *event;
        Double_t reftime;
        Double_t tLow;
        Double_t tHi;
        Double_t fLKr_minD;
	bool MCflag; 

};

#endif
