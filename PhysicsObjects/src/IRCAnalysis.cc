#include "IRCAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoIRCEvent.hh"
#include "Event.hh"

IRCAnalysis::IRCAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}

void IRCAnalysis::Set(TRecoIRCEvent* fevent, double freftime, double ftLow, double ftHi, bool fMCflag ){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
MCflag = fMCflag;
/*
TRecoIRCCandidate* IRCCand;
for(int i=0; i<event->GetNCandidates(); i++){
  IRCCand = (TRecoIRCCandidate*)event->GetCandidate(i);
  Double_t Tdiff = reftime - IRCCand->GetTime();
  fUserMethods->FillHisto("IRC_Set_reft", Tdiff);
}
*/
 TRecoIRCHit* IRCHit;
 TClonesArray& Hits = (*(event->GetHits()));
for(int i=0; i<event->GetNHits(); i++){
  IRCHit = (TRecoIRCHit*)Hits[i];
  Double_t Tdiff = reftime - IRCHit->GetTime();
  fUserMethods->FillHisto("IRC_Set_reft", Tdiff);

}
}
Int_t IRCAnalysis::Veto(){
int IRC_N = 0;
/*
TRecoIRCCandidate *IRCCand;
for(int i=0; i<event->GetNCandidates(); i++) {
 IRCCand = (TRecoIRCCandidate*)event->GetCandidate(i);
 double Tdiff = reftime - IRCCand->GetTime();
 if(Tdiff<tLow||Tdiff>tHi)continue;
 IRC_N++;
}
*/
 TRecoIRCHit* IRCHit;
 TClonesArray& Hits = (*(event->GetHits()));
for(int i=0; i<event->GetNHits(); i++) {
 IRCHit = (TRecoIRCHit*)Hits[i];
 double Tdiff = reftime - IRCHit->GetTime();
if(!MCflag){ if(Tdiff<tLow||Tdiff>tHi)continue;}
 IRC_N++;
}

if (IRC_N>0)return 0;
if (IRC_N==0)return 1;

}

void IRCAnalysis::BookHistos_Set(){
fUserMethods->BookHisto(new TH1F("IRC_Set_reft", "IRC Reference Time", 400, -50, 50)); 
}

void IRCAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots(); 
}
