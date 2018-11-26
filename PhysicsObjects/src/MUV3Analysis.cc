#include "MUV3Analysis.hh"
#include "AnalysisTools.hh"
#include "TRecoMUV3Event.hh"
#include "Event.hh"
//#include "Parameters.hh"

// fFlag == 0: call from Track.cc
// fFlag == 1; call from TwoPhotonAnalysis.cc

MUV3Analysis::MUV3Analysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba); 
}

void MUV3Analysis::Set(TRecoMUV3Event* fevent, double freftime, double ftLow, double ftHi, bool fMCflag ){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
MCflag = fMCflag;
fMUV3_minD = -1;
TRecoMUV3Candidate* MUV3Cand;
for(int i=0; i<event->GetNCandidates(); i++){
  MUV3Cand = (TRecoMUV3Candidate*)event->GetCandidate(i);
  Double_t Tdiff = reftime - MUV3Cand->GetTime();
  fUserMethods->FillHisto("MUV3_Set_reft", Tdiff);
}
}


Int_t MUV3Analysis::Veto(){				
int MUV3_N = 0;
TRecoMUV3Candidate *MUV3Cand;
for(int i=0; i<event->GetNCandidates(); i++) {
 MUV3Cand = (TRecoMUV3Candidate*)event->GetCandidate(i);
 double Tdiff = reftime - MUV3Cand->GetTime();
 if(!MCflag){if(Tdiff<tLow||Tdiff>tHi)continue;}
 MUV3_N++;
}
if (MUV3_N>0)return 0;
if (MUV3_N==0)return 1;

}
Int_t MUV3Analysis::VetoL(){				
int MUV3_N = 0;
TRecoMUV3Candidate *MUV3Cand;
for(int i=0; i<event->GetNCandidates(); i++) {
 MUV3Cand = (TRecoMUV3Candidate*)event->GetCandidate(i);
 double Tdiff = reftime - MUV3Cand->GetTime();
 if(!MCflag){if(Tdiff<-10||Tdiff>10)continue;}
 MUV3_N++;
}
if (MUV3_N>0)return 0;
if (MUV3_N==0)return 1;

}

Int_t MUV3Analysis::MuonID(TVector3 Pos){
 TRecoMUV3Candidate *MUV3Cand;
// Double_t MUV3StartPos = 244341;
// TVector3 Pion_Prop_MUV3 = Pos;
TVector3 Pion_Prop_MUV3 = Pos;
 int MUV3 = 0;
 int MUV3_minD;
 TVector2 Pion_Prop_MUV3_2D(Pion_Prop_MUV3.X(), Pion_Prop_MUV3.Y());
 double distanceMUV3 = 99999;
 for(int i=0; i<event->GetNCandidates(); i++)        {
  MUV3Cand = (TRecoMUV3Candidate*)event->GetCandidate(i);
 double Tdiff_MUV3 = reftime-MUV3Cand->GetTime();
 if(!MCflag){if(Tdiff_MUV3 < tLow || Tdiff_MUV3 > tHi)continue;}
  TVector2 HitPos(MUV3Cand->GetPosition().X() , MUV3Cand->GetPosition().Y());
  TVector2 distMUV3 = Pion_Prop_MUV3_2D - HitPos;
   if (distMUV3.Mod() < distanceMUV3)      {
    distanceMUV3 = distMUV3.Mod();
	fMUV3_minD = distanceMUV3;
	MUV3_minD = i;
     if (distanceMUV3<200)   {
      MUV3++;
}
}
}
  fUserMethods->FillHisto("MUV3_MuonID_MinD", distanceMUV3);
  if(MUV3==0)return 0;
  if(MUV3>0)return 1;
}

Int_t MUV3Analysis::Coincidental(){
int MUV3_N = 0;
TRecoMUV3Candidate *MUV3Cand;
for(int i=0; i<event->GetNCandidates(); i++) {
 MUV3Cand = (TRecoMUV3Candidate*)event->GetCandidate(i);
 double Tdiff = reftime - MUV3Cand->GetTime();
 if(!MCflag){if(Tdiff<tLow||Tdiff>tHi)continue;}
 MUV3_N++;
}
if (MUV3_N==0||MUV3_N==1)return 0;
if (MUV3_N>=2)return 1;

}



void MUV3Analysis::BookHistos_MuonID(){
fUserMethods->BookHisto(new TH1I("MUV3_MuonID_MinD", " Distance of minimum distance track to closest matched hit", 200, 0, 200));
}

void MUV3Analysis::BookHistos_Set(){
fUserMethods->BookHisto(new TH1I("MUV3_Set_reft", " reference time", 400, -50, 50));
}
void MUV3Analysis::SaveAllPlots(){
fUserMethods->SaveAllPlots(); 
}

Double_t MUV3Analysis::GetMinDistance(){
return fMUV3_minD;
}
