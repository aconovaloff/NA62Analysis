#include "OneTrack.hh"
#include "CHODAnalysis.hh"
#include "CedarAnalysis.hh"
#include "GigaTrackerAnalysis.hh"
#include "CHANTIAnalysis.hh"
#include "LKrAnalysis.hh"
#include "RICHAnalysis_G.hh"
#include "RICHCandidate.hh"
#include "AnalysisTools.hh"
#include "TRecoSpectrometerEvent.hh"
#include "Event.hh"
#include "Parameters.hh"
#include "BlueTubeTrackerAnalysis.hh"
using namespace TMath;

OneTrack::OneTrack(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
Strawa = new StrawAnalysis(ba);
CHODa = new CHODAnalysis(ba);
Cedara = new CedarAnalysis(ba);
Giga = new GigaTrackerAnalysis(ba);
CHANTIa = new CHANTIAnalysis(ba);
LKra = new LKrAnalysis(ba); 
RICHa_G = new RICHAnalysis_G(ba); 
Bluea = new BlueTubeTrackerAnalysis(ba); 
}

void OneTrack::Set(TRecoSpectrometerEvent* fSpectrometerEvent, TRecoCHODEvent* fCHODEvent, TRecoCedarEvent* fCedarEvent, TRecoGigaTrackerEvent* fGigaTrackerEvent, TRecoCHANTIEvent* fCHANTIEvent, TRecoLKrEvent* fLKrEvent, TRecoRICHEvent* fRICHEvent, L0TPData* fL0, Bool_t fMCflag, Int_t fRunNumber){
SpectrometerEvent = fSpectrometerEvent;
CHODEvent = fCHODEvent;
CedarEvent = fCedarEvent;
GigaTrackerEvent = fGigaTrackerEvent;
CHANTIEvent = fCHANTIEvent;
LKrEvent = fLKrEvent; 
RICHEvent = fRICHEvent;
L0 = fL0;
MCflag = fMCflag;
RunNumber = fRunNumber; 
}

Int_t OneTrack::OneTrackSelection(){
Bool_t beam_bg=false;
Parameters* P= new Parameters(3809);
Reft=999999;
//L0TPData *L0 = GetL0Data();
Strawa->Set(SpectrometerEvent, 1,P->GettStrawLow(),P->GettStrawHi(), MCflag); 
if (Strawa->SingleTrack()==0)return 0;
PionMomentum = Strawa->GetPionMomentum();
fUserMethods->FillHisto("pionmomentum", PionMomentum.Mag());
////////////////////Beam Background/////////////////////////
if(beam_bg){if (PionMomentum.Mag()<74||PionMomentum.Mag()>76)return 0;}
//////////////////////////////////////////////////////////
TLorentzVector PionFourMomentum = Strawa->GetPionFourMomentum();
PionPosition = Strawa->GetPionPosition();
TRecoSpectrometerCandidate* StrawCand;
StrawCand = (TRecoSpectrometerCandidate*)SpectrometerEvent->GetCandidate(0);
CHODa->Set(CHODEvent, StrawCand->GetTime(), P->GettStrawLow(),P->GettStrawHi(), MCflag);
 
TVector3 CHOD_pos= prop(PionPosition, PionMomentum*1000, P->GetCHODStartPos());
//TRecoSpectrometerCandidate* StrawCand;
//StrawCand = (TRecoSpectrometerCandidate*)SpectrometerEvent->GetCandidate(0);
if (CHODa->MatchedTrack(CHOD_pos)==0)return 0;
//if (CHODa->MinDisc(CHOD_pos, StrawCand->GetTime())==0)return;
Reft = CHODa->GetTime(); 
Double_t unit = 24.951059536/256;
fUserMethods->FillHisto("trigger_time_CHOD", Reft - L0->GetReferenceFineTime()*unit);
if(!MCflag){if (abs(Reft - L0->GetReferenceFineTime()*unit)>25)return 0;} 
//TRecoCHODCandidate* CHODCand;
//CHODCand = (TRecoCHODCandidate*)CHODEvent->GetCandidate(CHODa->GetMinDiscID());
//Double_t Reft = CHODCand->GetTime();
Strawa->MatchedTime(Reft); 
if(Strawa->MatchedTime(Reft)==0)return 0;

Cedara->Set(CedarEvent, Reft,P->GettCedarLow(),P->GettCedarHi(), MCflag);
fUserMethods->FillHisto("nCand", CedarEvent->GetNCandidates());
////////////////////////Beam Background//////////////////////
if(!beam_bg){if (Cedara->MatchedKaon()==0)return 0;}
//Cedara->MatchedKaon(); 
///////////////////////////////////////////////////////////fUserMethods->
fUserMethods->FillHisto("nSectors", Cedara->GetnSectors());
fUserMethods->FillHisto("p_Cedar", PionMomentum.Mag());
////////////////////////Beam Background////////////////////
if(!beam_bg){if (Cedara->GetnSectors()<5)return 0;}
//////////////////////////////////////////////////////////
fUserMethods->FillHisto("nCand_post", CedarEvent->GetNCandidates());
fUserMethods->FillHisto("p_Cedar_post", PionMomentum.Mag());
fUserMethods->FillHisto("nSectors_post", Cedara->GetnSectors());
Giga->Set(GigaTrackerEvent, Reft, P->GettGTKLow(), P->GettGTKHi(), MCflag);
TLorentzVector KaonFourMomentum = Giga->GetKaonFourMomentum(); 
TVector3 KaonPosition = Giga->GetKaonPosition(); 
TVector3 KaonMomentum = Giga->GetKaonMomentum(); 
double theta = ACos(PionMomentum*KaonMomentum/(PionMomentum.Mag()*KaonMomentum.Mag()));
Double_t mK = P->GetKaonMass();
Double_t mPi = P->GetPionMass(); 
Double_t mmsq = mK*mK*(1 - PionMomentum.Mag()/KaonMomentum.Mag()) + mPi*mPi*(1 - KaonMomentum.Mag()/PionMomentum.Mag()) - KaonMomentum.Mag()*PionMomentum.Mag()*theta*theta;
fUserMethods->FillHisto("mmsq_Straw", mmsq);
fUserMethods->FillHisto("mmsq_p",PionMomentum.Mag(),mmsq);
Int_t gtkmatched = Giga->MatchedTrack(PionMomentum, PionPosition, RunNumber);
if (gtkmatched==0)return 0; //took out for MC
KaonFourMomentum=Giga->GetKaonFourMomentum();
KaonPosition = Giga->GetKaonPosition();
KaonMomentum = Giga->GetKaonMomentum();
theta = ACos(PionMomentum*KaonMomentum/(PionMomentum.Mag()*KaonMomentum.Mag()));
mmsq = mK*mK*(1 - PionMomentum.Mag()/KaonMomentum.Mag()) + mPi*mPi*(1 - KaonMomentum.Mag()/PionMomentum.Mag()) - KaonMomentum.Mag()*PionMomentum.Mag()*theta*theta;
fUserMethods->FillHisto("mmsq_gtk", mmsq);
fUserMethods->FillHisto("mmsq_p_mid",PionMomentum.Mag(),mmsq);
Double_t close = cda(PionMomentum, KaonMomentum, PionPosition, KaonPosition);
fUserMethods->FillHisto("cda", close);
TVector3 Vert = Vertex(PionMomentum, KaonMomentum, PionPosition, KaonPosition);
Double_t VertZ = Vert.Z();
Double_t closest = cda(PionMomentum, KaonMomentum, PionPosition, KaonPosition); 
fUserMethods->FillHisto("cda_vertZ", Vert.Z(), closest);
fUserMethods->FillHisto("p_zvert",VertZ, PionMomentum.Mag());

BlueTubeTracker *gTracker = new BlueTubeTracker(); 
Bluea->Set(gTracker, KaonPosition, Vert, PionPosition, PionMomentum); 
TVector3 VertexDiff = Bluea->GetVertexCorrection() - Vert;
TLorentzVector NewMomentum = Bluea->GetMomentumCorrection();
TVector3 MomentumDiff = NewMomentum.Vect() - PionMomentum;
fUserMethods->FillHisto("Vertex_corr", VertexDiff.Mag());
fUserMethods->FillHisto("p_corr", MomentumDiff.Mag());
fUserMethods->FillHisto("vertex_z", Vert.Z());
if (VertZ < 110000 || VertZ > 165000) return 0;    
///////////////blue tube cuts///////////////////////////  

//Double_t mmsq_pre = mmsq; 
//PionMomentum = pionMomentumCorr.Vect(); 
//theta = ACos(PionMomentum*KaonMomentum/(PionMomentum.Mag()*KaonMomentum.Mag()));
//mmsq = mK*mK*(1 - PionMomentum.Mag()/KaonMomentum.Mag()) + mPi*mPi*(1 - KaonMomentum.Mag()/PionMomentum.Mag()) - KaonMomentum.Mag()*PionMomentum.Mag()*theta*theta;
//FillHisto("mmsq_btube", mmsq);
//FillHisto("mmsq_corr", mmsq_pre - mmsq); 
CHANTIa->Set(CHANTIEvent, Reft, P->GettCHANTILow(), P->GettCHANTIHi(), MCflag);
if(CHANTIa->Veto()==0) return 0; 
/*
if (RICHEvent->GetNRingCandidates()>1)return 0;
Int_t Multi = RICHa_G->TrackMultiRingMatching(Reft, RICHEvent, StrawCand, PionMomentum, Vert);
Int_t Single = RICHa_G->TrackSingleRingMatching(Reft, RICHEvent, StrawCand, PionMomentum, Vert);
RICHCandidate *MultiRICHCand;
RICHCandidate *SingleRICHCand;
Double_t SingleRICHMass;
Double_t MultiRICHMass;
Bool_t isSingle=false;
Bool_t isMulti=false;
if (Multi>-1){
MultiRICHCand = RICHa_G->GetRICHMultiCandidate();
isMulti = RICHa_G->RICHMultiCandidate(*MultiRICHCand, StrawCand, PionMomentum, Vert);
if(isMulti)fUserMethods->FillHisto("multimass",RICHa_G->Getmultirichmass());
fUserMethods->FillHisto("multi_disc",MultiRICHCand->GetDiscriminant());
fUserMethods->FillHisto("DistCheck", RICHa_G->GetDistCheck());
}

if (Single>-1){
 SingleRICHCand = RICHa_G->GetRICHSingleCandidate();
 isSingle = RICHa_G->RICHSingleCandidate(*SingleRICHCand, StrawCand, PionMomentum, Vert);
if(isSingle)fUserMethods->FillHisto("singlemass", RICHa_G->Getsinglerichmass());
fUserMethods->FillHisto("single_disc", SingleRICHCand->GetDiscriminant());

}
RICHmass=999999;
if(isSingle==true)RICHmass=RICHa_G->Getsinglerichmass();
else if (isMulti==true)RICHmass=RICHa_G->Getmultirichmass();
else return 0;
*/
if (MatchedRICHCandidate(RICHEvent, StrawCand, PionMomentum, Vert, Reft)==0)return 0; 
fUserMethods->FillHisto("RICHmass", RICHmass);
LKr_pos = prop(PionPosition, PionMomentum*1000, P->GetLKrStartPos());
Bool_t PiPiZero_background = false; 
LKra->Set(LKrEvent, Reft, P->GettLKrLow(), P->GettLKrHi(), MCflag);
Int_t number = LKra->SinglePionCluster(PionMomentum, LKr_pos, PiPiZero_background);
if (number==2||number==3)return 0;
//if (number == 3)return 0; //replace above line w/ this one for ke3 mc
if (number==0){if(LKra->TrackCellMatching(&LKr_pos)==0)return 0;};
Int_t index = LKra->GetSinglePionCluster_minID();
if(number==1)LKrCand = (TRecoLKrCandidate*)LKrEvent->GetCandidate(index);
//if(number==1||number==2)LKrCand = (TRecoLKrCandidate*)LKrEvent->GetCandidate(index); //same as four lines up
else LKrCand = LKra->GetLKrCellCandidate();
if (LKrCand->GetClusterEnergy()/PionMomentum.Mag()>0.8)return 0; //take out for ke3 MC
finalmmsq = mmsq; 
//fUserMethods->FillHisto("mmsq_gtkeff", mmsq);
//if (gtkmatched==1)fUserMethods->FillHisto("mmsq_gtkeff_post", mmsq);
fUserMethods->FillHisto("OneTrack_RICHmass_final", RICHmass);
fUserMethods->FillHisto("OneTrack_pi3_v_RICHrad", PionMomentum.Mag(),RingRadius);
fUserMethods->FillHisto("OneTrack_mmsq_final", mmsq);
fUserMethods->FillHisto("OneTrack_mmsq_p_final",PionMomentum.Mag(),mmsq);

Int_t bin = 10;
if (PionMomentum.Mag()<15)bin=0;
else if (PionMomentum.Mag()>=15&&PionMomentum.Mag()<20)bin=1;
else if (PionMomentum.Mag()>=20&&PionMomentum.Mag()<25)bin=2;
else if (PionMomentum.Mag()>=25&&PionMomentum.Mag()<30)bin=3;
else if (PionMomentum.Mag()>=30&&PionMomentum.Mag()<35)bin=4;
else bin=5;
for(int i=0; i<6; i++){
if(i==bin)fUserMethods->FillHisto(acc[i].Data(), PionMomentum.Mag(), mmsq);
}

return 1; 
}
void OneTrack::SaveAllPlots(){
fUserMethods->SaveAllPlots();
Strawa->SaveAllPlots();
CHODa->SaveAllPlots();
Cedara->SaveAllPlots();
Giga->SaveAllPlots();
CHANTIa->SaveAllPlots(); 
LKra->SaveAllPlots(); 
RICHa_G->SaveAllPlots();  
}
void OneTrack::BookHistos(){
//Strawa->BookHistos_Set();
fUserMethods->BookHisto(new TH1F("pionmomentum", "pion momentum", 200, 0, 100));
 fUserMethods->BookHisto(new TH1F("trigger_time_CHOD", "CHOD time minus trigger fine time", 200, -100,100));
fUserMethods->BookHisto(new TH1F("nCand", "Cand", 9, 0, 9));
fUserMethods->BookHisto(new TH1F("nCand_post", "Cand", 9, 0, 9));
fUserMethods->BookHisto(new TH1F("nSectors", "nSectors", 9, 0, 9));
fUserMethods->BookHisto(new TH1F("nSectors_post", "nSectors", 9, 0, 9));
fUserMethods->BookHisto(new TH1F("p_Cedar", "pi momenutm", 100, 0, 100));
fUserMethods->BookHisto(new TH1F("p_Cedar_post", "pi momenutm", 100, 0, 100));
 fUserMethods->BookHisto(new TH1F("mmsq_Straw", "Missing Mass squared (matched CHOD cand and CEDAR Kaon) ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("mmsq_gtk", "Missing Mass Square (Matched GTK track)", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("mmsq_gtkeff", "Missing Mass Square ", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("mmsq_gtkeff_post", "Missing Mass Square", 250, -.35, .15));
TH2I* mmsq_p = new TH2I("mmsq_p", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
mmsq_p->SetOption("COLZ");
mmsq_p->Draw();
fUserMethods->BookHisto(mmsq_p);
TH2I* mmsq_p_mid = new TH2I("mmsq_p_mid", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
mmsq_p_mid->SetOption("COLZ");
mmsq_p_mid->Draw();
fUserMethods->BookHisto(mmsq_p_mid);
TH2I* cda_vertZ = new TH2I("cda_vertZ", "z vertex (x-axis), cda (y-axis)", 180, 0, 180000, 20, 0, 20);
cda_vertZ->SetOption("COLZ");
cda_vertZ->Draw();
fUserMethods->BookHisto(cda_vertZ);
        fUserMethods->BookHisto(new TH1F("cda", "CDA", 200, 0, 100));
TH2I* p_zvert = new TH2I("p_zvert", "z coord of vertex (x-axis), momentum (y-axis)", 180, 0, 180000, 100, 0, 100);
p_zvert->SetOption("COLZ");
p_zvert->Draw();
fUserMethods->BookHisto(p_zvert);
fUserMethods->BookHisto(new TH1F("Vertex_corr", "blue tube mmsq correction", 100, 0, 100));
fUserMethods->BookHisto(new TH1F("p_corr", "blue tube momentum correction", 100, 0, 10));
 fUserMethods->BookHisto(new TH1F("vertex_z", "decay vertex z", 180, 0, 180000));
 fUserMethods->BookHisto(new TH1I("singlemass", "single rich mass", 140, 0,1.4));
        fUserMethods->BookHisto(new TH1I("multimass", "single rich mass", 140, 0,1.4));
 fUserMethods->BookHisto(new TH1I("multi_disc", "discriminant of multi rich candidates", 220, 0,110));
        fUserMethods->BookHisto(new TH1I("single_disc", "discriminant of single rich candidates", 220, 0,110));
        fUserMethods->BookHisto(new TH1I("DistCheck", "distance check of multi rich cand", 800, 0,200));
        fUserMethods->BookHisto(new TH1I("RICHmass", "rich mass of best candidate", 140, 0,1.4));
fUserMethods->BookHisto(new TH1I("OneTrack_RICHmass_final", "RICH mass after single track selection cuts", 140, 0,1.4));
TH2I* OneTrack_pi3_v_RICHrad = new TH2I("OneTrack_pi3_v_RICHrad", "Pion Momentum vs RICH Ring Radius", 100, 0, 100, 100, 0, 300);
OneTrack_pi3_v_RICHrad->SetOption("COLZ");
OneTrack_pi3_v_RICHrad->Draw();
fUserMethods->BookHisto(OneTrack_pi3_v_RICHrad);
fUserMethods->BookHisto(new TH1F("OneTrack_mmsq_final", "Missing Mass squared after all single track cuts ", 250, -.35, .15));
TH2I* OneTrack_mmsq_p_final = new TH2I("OneTrack_mmsq_p_final", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
OneTrack_mmsq_p_final->SetOption("COLZ");
OneTrack_mmsq_p_final->Draw();
fUserMethods->BookHisto(OneTrack_mmsq_p_final);
TH2I *acc_h[6]; 
for(int i=0; i<6; i++){
  acc[i].Form("OneTrack_acc_%d", i);
  acc_h[i] = new TH2I(acc[i].Data(), "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
  acc_h[i]->SetOption("COLZ");
  acc_h[i]->Draw();
  fUserMethods->BookHisto(acc_h[i]);
}

	CHODa->BookHistos_Set();
        CHODa->BookHistos_MatchedTrack();
        Strawa->BookHistos_MatchedTime();
	Strawa->BookHistos_SingleTrack(); 
        Giga->BookHistos_Set();
        Giga->BookHistos_MatchedTrack();	
	Cedara->BookHistos_Set();
	CHANTIa->BookHistos_Set();
	LKra->BookHistos_Set();
	LKra->BookHistos_SinglePionCluster();  
}

Int_t OneTrack::MatchedRICHCandidate(TRecoRICHEvent* rRICHEvent, TRecoSpectrometerCandidate* rStrawCand, TVector3 rPionMomentum, TVector3 rVert, Double_t rReft){
if (RICHEvent->GetNRingCandidates()>1)return 0;
Int_t Multi = RICHa_G->TrackMultiRingMatching(rReft, rRICHEvent, rStrawCand, rPionMomentum, rVert);
Int_t Single = RICHa_G->TrackSingleRingMatching(rReft, rRICHEvent, rStrawCand, rPionMomentum, rVert);
RICHCandidate *MultiRICHCand;
RICHCandidate *SingleRICHCand;
Double_t SingleRICHMass;
Double_t MultiRICHMass;
Bool_t isSingle=false;
Bool_t isMulti=false;
if (Multi>-1){
MultiRICHCand = RICHa_G->GetRICHMultiCandidate();
isMulti = RICHa_G->RICHMultiCandidate(*MultiRICHCand, rStrawCand, rPionMomentum, rVert);
if(isMulti)fUserMethods->FillHisto("multimass",RICHa_G->Getmultirichmass());
fUserMethods->FillHisto("multi_disc",MultiRICHCand->GetDiscriminant());
fUserMethods->FillHisto("DistCheck", RICHa_G->GetDistCheck());
}

if (Single>-1){
 SingleRICHCand = RICHa_G->GetRICHSingleCandidate();
 isSingle = RICHa_G->RICHSingleCandidate(*SingleRICHCand, rStrawCand, rPionMomentum, rVert);
if(isSingle)fUserMethods->FillHisto("singlemass", RICHa_G->Getsinglerichmass());
fUserMethods->FillHisto("single_disc", SingleRICHCand->GetDiscriminant());

}
RICHmass=999999;
if(isSingle==true){
 RICHmass=RICHa_G->Getsinglerichmass();
 RingRadius=SingleRICHCand->GetRadius();
} 
else if (isMulti==true){
 RICHmass=RICHa_G->Getmultirichmass();
 RingRadius=MultiRICHCand->GetRadius();
}
else return 0;
if (RICHmass>1.4)return 0;
fUserMethods->FillHisto("RICHmass", RICHmass);
return 1; 
}

Int_t OneTrack::BackgroundSelection(){
Parameters* P= new Parameters(3809);
//L0TPData *L0 = GetL0Data();
if(CedarEvent->GetNCandidates()!=0)return 0;
Strawa->Set(SpectrometerEvent, 1,P->GettStrawLow(),P->GettStrawHi(), MCflag); 
if (Strawa->SingleTrack()==0)return 0;
PionMomentum = Strawa->GetPionMomentum();
fUserMethods->FillHisto("pionmomentum", PionMomentum.Mag());
if(PionMomentum.Mag()<74||PionMomentum.Mag()>76)return 0;
TLorentzVector PionFourMomentum = Strawa->GetPionFourMomentum();
PionPosition = Strawa->GetPionPosition();
TRecoSpectrometerCandidate* StrawCand;
StrawCand = (TRecoSpectrometerCandidate*)SpectrometerEvent->GetCandidate(0);
CHODa->Set(CHODEvent, StrawCand->GetTime(), P->GettStrawLow(),P->GettStrawHi(), MCflag);
 
TVector3 CHOD_pos= prop(PionPosition, PionMomentum*1000, P->GetCHODStartPos());
//TRecoSpectrometerCandidate* StrawCand;
//StrawCand = (TRecoSpectrometerCandidate*)SpectrometerEvent->GetCandidate(0);
if (CHODa->MatchedTrack(CHOD_pos)==0)return 0;
//if (CHODa->MinDisc(CHOD_pos, StrawCand->GetTime())==0)return;
Reft = CHODa->GetTime();
Double_t unit = 24.951059536/256;
fUserMethods->FillHisto("trigger_time_CHOD", Reft - L0->GetReferenceFineTime()*unit);
if (abs(Reft - L0->GetReferenceFineTime()*unit)>25)return 0; 
//TRecoCHODCandidate* CHODCand;
//CHODCand = (TRecoCHODCandidate*)CHODEvent->GetCandidate(CHODa->GetMinDiscID());
//Double_t Reft = CHODCand->GetTime();
Strawa->MatchedTime(Reft); 
if(Strawa->MatchedTime(Reft)==0)return 0;

Cedara->Set(CedarEvent, Reft,P->GettCedarLow(),P->GettCedarHi(), MCflag);
fUserMethods->FillHisto("nCand", CedarEvent->GetNCandidates());
fUserMethods->FillHisto("nSectors", Cedara->GetnSectors());
fUserMethods->FillHisto("p_Cedar", PionMomentum.Mag());
fUserMethods->FillHisto("nCand_post", CedarEvent->GetNCandidates());
fUserMethods->FillHisto("p_Cedar_post", PionMomentum.Mag());
fUserMethods->FillHisto("nSectors_post", Cedara->GetnSectors());
Giga->Set(GigaTrackerEvent, Reft, P->GettGTKLow(), P->GettGTKHi(), MCflag);
TLorentzVector KaonFourMomentum = Giga->GetKaonFourMomentum(); 
TVector3 KaonPosition = Giga->GetKaonPosition(); 
TVector3 KaonMomentum = Giga->GetKaonMomentum(); 
double theta = ACos(PionMomentum*KaonMomentum/(PionMomentum.Mag()*KaonMomentum.Mag()));
Double_t mK = P->GetKaonMass();
Double_t mPi = P->GetPionMass(); 
Double_t mmsq = mK*mK*(1 - PionMomentum.Mag()/KaonMomentum.Mag()) + mPi*mPi*(1 - KaonMomentum.Mag()/PionMomentum.Mag()) - KaonMomentum.Mag()*PionMomentum.Mag()*theta*theta;
fUserMethods->FillHisto("mmsq_Straw", mmsq);
fUserMethods->FillHisto("mmsq_p",PionMomentum.Mag(),mmsq);
//if (Giga->MatchedTrack(PionMomentum, PionPosition, RunNumber)==0)return 0;
KaonFourMomentum=Giga->GetKaonFourMomentum();
KaonPosition = Giga->GetKaonPosition();
KaonMomentum = Giga->GetKaonMomentum();
theta = ACos(PionMomentum*KaonMomentum/(PionMomentum.Mag()*KaonMomentum.Mag()));
mmsq = mK*mK*(1 - PionMomentum.Mag()/KaonMomentum.Mag()) + mPi*mPi*(1 - KaonMomentum.Mag()/PionMomentum.Mag()) - KaonMomentum.Mag()*PionMomentum.Mag()*theta*theta;
fUserMethods->FillHisto("mmsq_gtk", mmsq);
fUserMethods->FillHisto("mmsq_p_mid",PionMomentum.Mag(),mmsq);
Double_t close = cda(PionMomentum, KaonMomentum, PionPosition, KaonPosition);
fUserMethods->FillHisto("cda", close);
TVector3 Vert = Vertex(PionMomentum, KaonMomentum, PionPosition, KaonPosition);
Double_t VertZ = Vert.Z();
Double_t closest = cda(PionMomentum, KaonMomentum, PionPosition, KaonPosition); 
fUserMethods->FillHisto("cda_vertZ", Vert.Z(), closest);
fUserMethods->FillHisto("p_zvert",VertZ, PionMomentum.Mag());

BlueTubeTracker *gTracker = new BlueTubeTracker(); 
Bluea->Set(gTracker, KaonPosition, Vert, PionPosition, PionMomentum); 
TVector3 VertexDiff = Bluea->GetVertexCorrection() - Vert;
TLorentzVector NewMomentum = Bluea->GetMomentumCorrection();
TVector3 MomentumDiff = NewMomentum.Vect() - PionMomentum;
fUserMethods->FillHisto("Vertex_corr", VertexDiff.Mag());
fUserMethods->FillHisto("p_corr", MomentumDiff.Mag());
fUserMethods->FillHisto("vertex_z", Vert.Z());
if (VertZ < 110000 || VertZ > 165000) return 0;    
///////////////blue tube cuts///////////////////////////  

//Double_t mmsq_pre = mmsq; 
//PionMomentum = pionMomentumCorr.Vect(); 
//theta = ACos(PionMomentum*KaonMomentum/(PionMomentum.Mag()*KaonMomentum.Mag()));
//mmsq = mK*mK*(1 - PionMomentum.Mag()/KaonMomentum.Mag()) + mPi*mPi*(1 - KaonMomentum.Mag()/PionMomentum.Mag()) - KaonMomentum.Mag()*PionMomentum.Mag()*theta*theta;
//FillHisto("mmsq_btube", mmsq);
//FillHisto("mmsq_corr", mmsq_pre - mmsq); 
CHANTIa->Set(CHANTIEvent, Reft, P->GettCHANTILow(), P->GettCHANTIHi(), MCflag);
if(CHANTIa->Veto()==0) return 0; 
/*
if (RICHEvent->GetNRingCandidates()>1)return 0;
Int_t Multi = RICHa_G->TrackMultiRingMatching(Reft, RICHEvent, StrawCand, PionMomentum, Vert);
Int_t Single = RICHa_G->TrackSingleRingMatching(Reft, RICHEvent, StrawCand, PionMomentum, Vert);
RICHCandidate *MultiRICHCand;
RICHCandidate *SingleRICHCand;
Double_t SingleRICHMass;
Double_t MultiRICHMass;
Bool_t isSingle=false;
Bool_t isMulti=false;
if (Multi>-1){
MultiRICHCand = RICHa_G->GetRICHMultiCandidate();
isMulti = RICHa_G->RICHMultiCandidate(*MultiRICHCand, StrawCand, PionMomentum, Vert);
if(isMulti)fUserMethods->FillHisto("multimass",RICHa_G->Getmultirichmass());
fUserMethods->FillHisto("multi_disc",MultiRICHCand->GetDiscriminant());
fUserMethods->FillHisto("DistCheck", RICHa_G->GetDistCheck());
}

if (Single>-1){
 SingleRICHCand = RICHa_G->GetRICHSingleCandidate();
 isSingle = RICHa_G->RICHSingleCandidate(*SingleRICHCand, StrawCand, PionMomentum, Vert);
if(isSingle)fUserMethods->FillHisto("singlemass", RICHa_G->Getsinglerichmass());
fUserMethods->FillHisto("single_disc", SingleRICHCand->GetDiscriminant());

}
RICHmass=999999;
if(isSingle==true)RICHmass=RICHa_G->Getsinglerichmass();
else if (isMulti==true)RICHmass=RICHa_G->Getmultirichmass();
else return 0;
*/
if (MatchedRICHCandidate(RICHEvent, StrawCand, PionMomentum, Vert, Reft)==0)return 0; 
fUserMethods->FillHisto("RICHmass", RICHmass);
LKr_pos = prop(PionPosition, PionMomentum*1000, P->GetLKrStartPos());
Bool_t PiPiZero_background = false; 
LKra->Set(LKrEvent, Reft, P->GettLKrLow(), P->GettLKrHi(), MCflag);
Int_t number = LKra->SinglePionCluster(PionMomentum, LKr_pos, PiPiZero_background);
if (number==2||number==3)return 0;
if (number==0){if(LKra->TrackCellMatching(&LKr_pos)==0)return 0;};
Int_t index = LKra->GetSinglePionCluster_minID();
if(number==1)LKrCand = (TRecoLKrCandidate*)LKrEvent->GetCandidate(index);
else LKrCand = LKra->GetLKrCellCandidate();
if (LKrCand->GetClusterEnergy()/PionMomentum.Mag()>0.8)return 0;
finalmmsq = mmsq; 

return 1; 
}
