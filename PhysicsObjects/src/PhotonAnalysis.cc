#include "PhotonAnalysis.hh"
#include "CHODAnalysis.hh"
#include "CedarAnalysis.hh"
#include "GigaTrackerAnalysis.hh"
#include "CHANTIAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoSpectrometerEvent.hh"
#include "Event.hh"
#include "Parameters.hh"
#include "LKrAnalysis_G.hh"
#include "LAVMatching.hh"
#include "SAVMatching.hh"
#include "BlueTubeTrackerAnalysis.hh"
#include "MUV1Analysis.hh"
#include "MUV2Analysis.hh"
#include "MUV3Analysis.hh"
#include "SAVAnalysis.hh"
#include "LAVAnalysis.hh"
#include "IRCAnalysis.hh"
#include "SACAnalysis.hh"
#include "LAVMatchingMC.hh"
#include "SAVMatchingMC.hh"

using namespace TMath;

PhotonAnalysis::PhotonAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
LKra_G = new LKrAnalysis_G(1, ba);
SAVa = new SAVAnalysis(ba);
LAVa = new LAVAnalysis(ba);
IRCa = new IRCAnalysis(ba);
SACa = new SACAnalysis(ba);  
}

void PhotonAnalysis::Set(TRecoLKrEvent* fLKrEvent, TRecoLAVEvent* fLAVEvent, TRecoIRCEvent* fIRCEvent, TRecoSACEvent* fSACEvent, TRecoSAVEvent* fSAVEvent, TRecoLKrCandidate* fLKrCandidate, Double_t fReft, Double_t fmmsq, Int_t fUnmatchedClusters, TVector3 fLKr_pos, Bool_t fMCflag){
LAVEvent = fLAVEvent;
IRCEvent = fIRCEvent;
SACEvent = fSACEvent;
SAVEvent = fSAVEvent;
LKrEvent = fLKrEvent;
LKrCand = fLKrCandidate;
LKr_pos = fLKr_pos;
Reft = fReft;
mmsq = fmmsq;
UnmatchedClusters = fUnmatchedClusters; 
MCflag = fMCflag;
// fLAVMatching = new LAVMatching();
//fSAVMatching = new SAVMatching(); 
}

Int_t PhotonAnalysis::Rejection(){
Parameters* P = new Parameters(3809);
fUserMethods->FillHisto("mmsq_photon_initial", mmsq);
if (UnmatchedClusters > 0)return 0; 
fUserMethods->FillHisto("mmsq_LKr", mmsq);
////////////////Giuseppe's Find Extra Cluster Stuff////////////////////////
TVector2 MatchedPionPos;
MatchedPionPos.Set(LKrCand->GetClusterX(), LKrCand->GetClusterY());
if(UnmatchedClusters ==0){
Int_t nlkrnew = LKra_G->FindNewClusters(Reft, LKrEvent,  &MatchedPionPos , &LKr_pos); 
Double_t timeref = Reft;
fUserMethods->FillHisto("new_clusters", nlkrnew); 
    Int_t nNewPhotonsInTime = 0;
    Int_t nNewPhotons = 0;
    Double_t emax = 0;
    for (Int_t jnew=0; jnew<nlkrnew; jnew++) {
      LKrCandidate *newphoton = LKra_G->GetNewClusterCandidate(jnew);
      if (newphoton->GetEnergy()<=1) continue;
      if (newphoton->GetEnergy()>emax) emax = newphoton->GetEnergy();
//      if (fYear==2015) {
        if (fabs((timeref-newphoton->GetTime())-3.0)<8) nNewPhotonsInTime++;
//      }
//      if (fYear==2016) {
//        if (fabs((timeref-newphoton->GetTime()))<8) nNewPhotonsInTime++;
//      }
      nNewPhotons++; //change to new photons in time. 
    }
if(nNewPhotons>0)return 0;
}
///////////////////////////////////////////////////////////////////////
fUserMethods->FillHisto("mmsq_LKr_G", mmsq);
//FillHisto("p_zvert_mid",VertZ, PionMomentum.Mag());
//FillHisto("mmsq_p_LKr",PionMomentum.Mag(),mmsq);
/*
LAVMatching* fLAVMatching = new LAVMatching();
fLAVMatching->SetReferenceTime(Reft);
Bool_t matched = fLAVMatching->LAVHasTimeMatching(LAVEvent);
if(matched)return 0;
fUserMethods->FillHisto("mmsq_LAV", mmsq);
SAVMatching *fSAVMatching = new SAVMatching();
fSAVMatching->SetReferenceTime(Reft);
Bool_t matched_SAV = fSAVMatching->SAVHasTimeMatching(IRCEvent, SACEvent);
if(matched_SAV)return 0;
fUserMethods->FillHisto("mmsq_SAV", mmsq);
SAVa->Set(SAVEvent, Reft, P->GettSAVLow(), P->GettSAVHi(), MCflag);
if (SAVa->Veto()==0)return 0;
fUserMethods->FillHisto("mmsq_SAVa", mmsq);
*/
if (OtherVetoes(LAVEvent, IRCEvent, SACEvent, SAVEvent, Reft, MCflag)==0)return 0; 
return 1;
}
void PhotonAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
SAVa->SaveAllPlots();
SACa->SaveAllPlots();
IRCa->SaveAllPlots();
LAVa->SaveAllPlots(); 
}
void PhotonAnalysis::BookHistos(){
fUserMethods->BookHisto(new TH1F("mmsq_photon_initial", "Missing Mass squared (No LKr photons)", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("mmsq_LKr", "Missing Mass squared (No LKr photons)", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("new_clusters", "found LKr clusters", 11, 0,11));
  fUserMethods->BookHisto(new TH1F("mmsq_LKr_G", "Missing Mass squared (No LKr photons)", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("mmsq_LAV", "Missing Mass squared (LAV Veto)", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("mmsq_SAV", "Missing Mass squared (SAC Veto)", 250, -.35, .15));
 fUserMethods->BookHisto(new TH1F("mmsq_SAVa", "Missing Mass squared (SAV cream veto)", 250, -.35, .15));
//SAVa->BookHistos_Set(); 
LAVa->BookHistos_Set(); 
SAVa->BookHistos_Set();
SACa->BookHistos_Set(); 
IRCa->BookHistos_Set();  
}

Int_t PhotonAnalysis::OtherVetoes(TRecoLAVEvent* rLAVEvent, TRecoIRCEvent* rIRCEvent, TRecoSACEvent* rSACEvent, TRecoSAVEvent* rSAVEvent, Double_t rReft, Bool_t rMCflag){  
Parameters* P = new Parameters(3809);
//for MC if SAV/LAV matching are not working in MC
LAVa->Set(rLAVEvent, rReft, P->GettLAVLow(), P->GettLAVHi(), rMCflag);
IRCa->Set(rIRCEvent, rReft, P->GettIRCLow(), P->GettIRCHi(), rMCflag);
SACa->Set(rSACEvent, rReft, P->GettSACLow(), P->GettSACHi(), rMCflag);

LAVMatching* fLAVMatching = new LAVMatching();
fLAVMatching->SetReferenceTime(rReft);
LAVMatchingMC* fLAVMatchingMC = new LAVMatchingMC();
fLAVMatchingMC->SetReferenceTime(rReft);
Bool_t matchedMC = fLAVMatchingMC->LAVHasTimeMatching(rLAVEvent); 
Bool_t matched = fLAVMatching->LAVHasTimeMatching(rLAVEvent);
//cout<<"mc: "<<matchedMC<<endl;
//cout<<"normal: "<<matched<<endl;//doesn't work for MC, digis? 
if(!rMCflag){if(matched)return 0;}
//else {if (LAVa->Veto()==0)return 0;}
else {if(matchedMC)return 0;}
fUserMethods->FillHisto("mmsq_LAV", mmsq);
SAVMatching *fSAVMatching = new SAVMatching();
fSAVMatching->SetReferenceTime(rReft);
Bool_t matched_SAV = fSAVMatching->SAVHasTimeMatching(rIRCEvent, rSACEvent);
SAVMatchingMC *fSAVMatchingMC = new SAVMatchingMC();
fSAVMatchingMC->SetReferenceTime(rReft);
Bool_t matched_SAVMC = fSAVMatchingMC->SAVHasTimeMatching(rIRCEvent, rSACEvent);

if(!rMCflag){if(matched_SAV)return 0;}
//else {if (SACa->Veto()==0||IRCa->Veto()==0)return 0;}
else {if(matched_SAVMC)return 0;}
fUserMethods->FillHisto("mmsq_SAV", mmsq);
SAVa->Set(rSAVEvent, rReft, P->GettSAVLow(), P->GettSAVHi(), rMCflag);
if (SAVa->Veto()==0)return 0;
fUserMethods->FillHisto("mmsq_SAVa", mmsq);
return 1;
}
