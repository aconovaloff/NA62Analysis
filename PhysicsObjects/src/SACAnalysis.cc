#include "SACAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoSACEvent.hh"
#include "Event.hh"

SACAnalysis::SACAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}

void SACAnalysis::Set(TRecoSACEvent* fevent, double freftime, double ftLow, double ftHi, bool fMCflag){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
MCflag = fMCflag;
/*
TRecoSACCandidate* SACCand;
for(int i=0; i<event->GetNCandidates(); i++){
  SACCand = (TRecoSACCandidate*)event->GetCandidate(i);
  Double_t Tdiff = reftime - SACCand->GetTime();
  fUserMethods->FillHisto("SAC_Set_reft", Tdiff);
}
*/
 TRecoSACHit* SACHit;
 TClonesArray& Hits = (*(event->GetHits())); 
for(int i=0; i<event->GetNHits(); i++){
  SACHit = (TRecoSACHit*)Hits[i];
  Double_t Tdiff = reftime - SACHit->GetTime();
  fUserMethods->FillHisto("SAC_Set_reft", Tdiff);

}
}
Int_t SACAnalysis::Veto(){
int SAC_N = 0;
/*
TRecoSACCandidate *SACCand;
for(int i=0; i<event->GetNCandidates(); i++) {
 SACCand = (TRecoSACCandidate*)event->GetCandidate(i);
 double Tdiff = reftime - SACCand->GetTime();
 if(Tdiff<tLow||Tdiff>tHi)continue;
 SAC_N++;
}
*/
 TRecoSACHit* SACHit;
 TClonesArray& Hits = (*(event->GetHits()));
for(int i=0; i<event->GetNHits(); i++) {
 SACHit = (TRecoSACHit*)Hits[i];
 double Tdiff = reftime - SACHit->GetTime();
 if(!MCflag){if(Tdiff<tLow||Tdiff>tHi)continue;}
 SAC_N++;
}

if (SAC_N>0)return 0;
if (SAC_N==0)return 1;

}

void SACAnalysis::BookHistos_Set(){
fUserMethods->BookHisto(new TH1F("SAC_Set_reft", "SAC Reference Time", 400, -50, 50)); 
}

void SACAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots(); 
}
