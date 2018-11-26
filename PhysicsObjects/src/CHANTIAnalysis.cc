#include "CHANTIAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoCHANTIEvent.hh"
#include "Event.hh"


CHANTIAnalysis::CHANTIAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}

void CHANTIAnalysis::Set(TRecoCHANTIEvent* fevent, double freftime, double ftLow, double ftHi, bool fMCflag ){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
MCflag = fMCflag;
TRecoCHANTICandidate* CHANTICand;
for(int i=0; i<event->GetNCandidates(); i++){
  CHANTICand = (TRecoCHANTICandidate*)event->GetCandidate(i);
  Double_t Tdiff = reftime - CHANTICand->GetTime();
  fUserMethods->FillHisto("CHANTI_Set_reft", Tdiff);
}

}

Int_t CHANTIAnalysis::Veto(){
int CHANTI_N = 0;
TRecoCHANTICandidate *CHANTICand;
for(int i=0; i<event->GetNCandidates(); i++) {
 CHANTICand = (TRecoCHANTICandidate*)event->GetCandidate(i);
 double Tdiff = reftime - CHANTICand->GetTime();
if(!MCflag){ if(Tdiff<tLow||Tdiff>tHi)continue;}
 CHANTI_N++;
}
if (CHANTI_N>0)return 0;
if (CHANTI_N==0)return 1;

}

void CHANTIAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
}

void CHANTIAnalysis::BookHistos_Set(){
fUserMethods->BookHisto(new TH1F("CHANTI_Set_reft","reference time", 800, -50, 50));
}
