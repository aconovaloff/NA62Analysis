#include "CedarAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoCedarEvent.hh"
#include "Event.hh"
using namespace TMath;
CedarAnalysis::CedarAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}

void CedarAnalysis::Set(TRecoCedarEvent* fevent, double freftime, double ftLow, double ftHi, bool fMCflag ){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
MCflag = fMCflag; 
nSectors = -1;
TRecoCedarCandidate* CedarCand;
for(int i=0; i<event->GetNCandidates(); i++){
  CedarCand = (TRecoCedarCandidate*)event->GetCandidate(i);
  Double_t Tdiff = reftime - CedarCand->GetTime();
  fUserMethods->FillHisto("Cedar_Set_reft", Tdiff);
}
}

Int_t CedarAnalysis::MatchedKaon(){
 TRecoCedarCandidate *CedarCand;
 int Cedar = 0; 
 int CedarMaxID;
 nSectors = 0;
 for(int i=0; i<event->GetNCandidates(); i++)        {
  CedarCand = (TRecoCedarCandidate*)event->GetCandidate(i);
  Double_t Tdiff = reftime - CedarCand->GetTime();
//fUserMethods->FillHisto("Cedar_MatchedKaon_reft", Tdiff); 
if(!MCflag) {if(Tdiff<tLow||Tdiff>tHi)continue;} 
   if (CedarCand->GetNSectors()> nSectors)      {
    nSectors = CedarCand->GetNSectors();
      Cedar++;
      CedarMaxID = i;
      fTime = CedarCand->GetTime();
}
}
//CedarCand = (TRecoCedarCandidate*)event->GetCandidate(CedarMaxID); 
 if (Cedar == 0)return 0;
 if (Cedar > 0)return 1;
/*
 if (Cedar != 1)return 0;
 else {
cout<<Cedar<<endl;
return 1; 
	}
*/
}

Int_t CedarAnalysis::GetnSectors(){
return nSectors; 
}

void CedarAnalysis::BookHistos_MatchedKaon(){
}
void CedarAnalysis::BookHistos_Set(){
fUserMethods->BookHisto(new TH1F("Cedar_Set_reft", "Reference Time", 800, -25, 25));
}

void CedarAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
}

Int_t CedarAnalysis::MatchedKaon_Straw(Double_t reft){
 TRecoCedarCandidate *CedarCand;
Double_t mindiff = 999999;
Int_t minID = -1;
 for(int i=0; i<event->GetNCandidates(); i++)        {
  CedarCand = (TRecoCedarCandidate*)event->GetCandidate(i);
  Double_t Tdiff = reft - CedarCand->GetTime();
if(!MCflag) {if(Tdiff<-10||Tdiff>10)continue;}////////get these values
   if (CedarCand->GetNSectors()<5)continue;
   if (abs(Tdiff)<mindiff){
      mindiff = Tdiff;
      minID = i;        
  }
 }
 return minID;
}

