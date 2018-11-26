#include "SAVAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoSAVEvent.hh"
#include "Event.hh"

SAVAnalysis::SAVAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}

void SAVAnalysis::Set(TRecoSAVEvent* fevent, double freftime, double ftLow, double ftHi, bool fMCflag ){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
MCflag = fMCflag;

TRecoSAVCandidate* SAVCand;
for(int i=0; i<event->GetNCandidates(); i++){
  SAVCand = (TRecoSAVCandidate*)event->GetCandidate(i);
  Double_t Tdiff = reftime - SAVCand->GetTime();
  fUserMethods->FillHisto("SAV_Set_reft", Tdiff);
}
/*
 TRecoSAVHit* SAVHit;
 TClonesArray& Hits = (*(event->GetHits()));
for(int i=0; i<event->GetNHits(); i++){
  SAVHit = (TRecoSAVHit*)Hits[i];
  Double_t Tdiff = reftime - SAVHit->GetTime();
  fUserMethods->FillHisto("SAV_Set_reft", Tdiff);

}
*/
}
Int_t SAVAnalysis::Veto(){
int SAV_N = 0;

TRecoSAVCandidate *SAVCand;
for(int i=0; i<event->GetNCandidates(); i++) {
 SAVCand = (TRecoSAVCandidate*)event->GetCandidate(i);
 double Tdiff = reftime - SAVCand->GetTime();
 if(Tdiff<tLow||Tdiff>tHi)continue;
 SAV_N++;
}
/*
 TRecoSAVHit* SAVHit;
 TClonesArray& Hits = (*(event->GetHits()));
for(int i=0; i<event->GetNHits(); i++) {
 SAVHit = (TRecoSAVHit*)Hits[i];
 double Tdiff = reftime - SAVHit->GetTime();
if(!MCflag){ if(Tdiff<tLow||Tdiff>tHi)continue;}
 SAV_N++;
}
*/
if (SAV_N>0)return 0;
if (SAV_N==0)return 1;

}

void SAVAnalysis::BookHistos_Set(){
fUserMethods->BookHisto(new TH1F("SAV_Set_reft", "SAV Reference Time", 400, -50, 50)); 
}

void SAVAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots(); 
}
