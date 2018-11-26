#include "TrackAnalysis.hh"
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
#include "MUV3Analysis.hh"
#include "Parameters.hh"
#include "BlueTubeTrackerAnalysis.hh"
#include "CalorimeterAnalysis.hh"
#include "CalorimeterPlotGenerator.hh"
#include "OneTrack.hh"
#include "PhotonAnalysis.hh"
using namespace TMath;

TrackAnalysis::TrackAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
Calo = new CalorimeterAnalysis(ba);
One = new OneTrack(ba);
Photon = new PhotonAnalysis(ba);
LKra = new LKrAnalysis(ba);  
MUV3a = new MUV3Analysis(ba);
CaloPlot = new CalorimeterPlotGenerator(ba);  
}

void TrackAnalysis::Set(TRecoCedarEvent* fCedarEvent, TRecoGigaTrackerEvent* fGigaTrackerEvent, TRecoCHANTIEvent* fCHANTIEvent, TRecoSpectrometerEvent* fSpectrometerEvent, TRecoCHODEvent* fCHODEvent, TRecoLKrEvent* fLKrEvent, TRecoRICHEvent* fRICHEvent, TRecoMUV1Event* fMUV1Event, TRecoMUV2Event* fMUV2Event, TRecoMUV3Event* fMUV3Event, TRecoLAVEvent* fLAVEvent,  TRecoIRCEvent* fIRCEvent, TRecoSACEvent* fSACEvent, TRecoSAVEvent* fSAVEvent, L0TPData* fL0, Bool_t fMCflag, Int_t fRunNumber){
SpectrometerEvent = fSpectrometerEvent; 
MUV1Event = fMUV1Event;
MUV2Event = fMUV2Event;
MUV3Event = fMUV3Event;
LAVEvent = fLAVEvent; 
IRCEvent = fIRCEvent; 
SACEvent = fSACEvent;
SAVEvent = fSAVEvent;
CHODEvent = fCHODEvent;
CedarEvent = fCedarEvent;
GigaTrackerEvent = fGigaTrackerEvent;
CHANTIEvent = fCHANTIEvent;
LKrEvent = fLKrEvent; 
RICHEvent = fRICHEvent;
L0 = fL0;
MCflag = fMCflag;
RunNumber = fRunNumber; 
One->Set(SpectrometerEvent, CHODEvent, CedarEvent, GigaTrackerEvent, CHANTIEvent, LKrEvent, RICHEvent,  L0, MCflag, RunNumber);
}

Int_t TrackAnalysis::PhotonRejectionFactor(){
Double_t mmsq = One->Getmmsq();
TVector3 PionMomentum = One->GetPionMomentum();
fUserMethods->FillHisto("PhotonRejectionFactor_mmsq_onetrack", mmsq); 
fUserMethods->FillHisto("PhotonRejectionFactor_mmsq_p_onetrack",PionMomentum.Mag(),mmsq);
if (Calo->Rejection()==0)return 0;
fUserMethods->FillHisto("PhotonRejectionFactor_mmsq_calo", mmsq);
if (One->GetRICHmass()<0.13)return 0;
fUserMethods->FillHisto("PhotonRejectionFactor_mmsq_RICH", mmsq);
fUserMethods->FillHisto("PhotonRejectionFactor_mmsq_p_muon",PionMomentum.Mag(),mmsq);
if ((PionMomentum.Mag()>15&&PionMomentum.Mag()<35)&&(mmsq>0.01&&mmsq<0.026))fUserMethods->FillHisto("PhotonRejectionFactor_mmsq_den", mmsq);
if(Photon->Rejection()==0)return 0;
fUserMethods->FillHisto("PhotonRejectionFactor_mmsq_photon", mmsq);
fUserMethods->FillHisto("PhotonRejectionFactor_mmsq_p_photon",PionMomentum.Mag(),mmsq);
if ((PionMomentum.Mag()>15&&PionMomentum.Mag()<35)&&(mmsq>0.01&&mmsq<0.026))fUserMethods->FillHisto("PhotonRejectionFactor_mmsq_num", mmsq);

fUserMethods->FillHisto("PhotonRejectionFactor_mmsq_end", mmsq); 
if (PionMomentum.Mag()>15&&PionMomentum.Mag()<35)fUserMethods->FillHisto("PhotonRejectionFactor_mmsq_end_signal", mmsq);
if ((PionMomentum.Mag()>15&&PionMomentum.Mag()<35)&&(mmsq>0&&mmsq<0.01||mmsq>0.026&&mmsq<0.068))return 1; 
}

void TrackAnalysis::CalorimeterRejectionFactor(){
if (L0->GetTriggerFlags()!=32)return; //took out, not sure why i put it here. 
Double_t mmsq = One->Getmmsq();
 TVector3 PionMomentum = One->GetPionMomentum();
fUserMethods->FillHisto("CalorimeterRejectionFactor_mmsq_onetrack", mmsq);
fUserMethods->FillHisto("CalorimeterRejectionFactor_mmsq_p_onetrack",PionMomentum.Mag(),mmsq);
if(Photon->Rejection()==0)return;
fUserMethods->FillHisto("CalorimeterRejectionFactor_mmsq_photon", mmsq);
fUserMethods->FillHisto("CalorimeterRejectionFactor_mmsq_p_photon",PionMomentum.Mag(),mmsq);
if ((PionMomentum.Mag()>15&&PionMomentum.Mag()<35)&&(mmsq>-0.02&&mmsq<0))fUserMethods->FillHisto("CalorimeterRejectionFactor_mmsq_den", mmsq);
if (Calo->Rejection()==0)return;
fUserMethods->FillHisto("CalorimeterRejectionFactor_mmsq_calo", mmsq);
fUserMethods->FillHisto("CalorimeterRejectionFactor_mmsq_p_calo",PionMomentum.Mag(),mmsq); 
if ((PionMomentum.Mag()>15&&PionMomentum.Mag()<35)&&(mmsq>-0.02&&mmsq<0))fUserMethods->FillHisto("CalorimeterRejectionFactor_mmsq_num", mmsq);
if (One->GetRICHmass()<0.13)return;
fUserMethods->FillHisto("CalorimeterRejectionFactor_mmsq_RICH", mmsq);
fUserMethods->FillHisto("CalorimeterRejectionFactor_mmsq_p_RICH",PionMomentum.Mag(),mmsq);
fUserMethods->FillHisto("MuonRejectionFactor_mmsq_end", mmsq);
if (PionMomentum.Mag()>15&&PionMomentum.Mag()<35)fUserMethods->FillHisto("MuonRejectionFactor_mmsq_end_signal", mmsq);
}

void TrackAnalysis::RICHRejectionFactor(){
Double_t mmsq = One->Getmmsq();
TVector3 PionMomentum = One->GetPionMomentum();
fUserMethods->FillHisto("RICHRejectionFactor_mmsq_onetrack", mmsq);
fUserMethods->FillHisto("RICHRejectionFactor_mmsq_p_onetrack",PionMomentum.Mag(),mmsq);
if(Photon->Rejection()==0)return;
fUserMethods->FillHisto("RICHRejectionFactor_mmsq_photon", mmsq);
fUserMethods->FillHisto("RICHRejectionFactor_mmsq_p_photon",PionMomentum.Mag(),mmsq);
if ((PionMomentum.Mag()>15&&PionMomentum.Mag()<35)&&(mmsq>-0.02&&mmsq<0))fUserMethods->FillHisto("RICHRejectionFactor_mmsq_den", mmsq);
if (One->GetRICHmass()<0.13)return;
fUserMethods->FillHisto("RICHRejectionFactor_mmsq_RICH", mmsq);
fUserMethods->FillHisto("RICHRejectionFactor_mmsq_p_RICH",PionMomentum.Mag(),mmsq);
if ((PionMomentum.Mag()>15&&PionMomentum.Mag()<35)&&(mmsq>-0.02&&mmsq<0))fUserMethods->FillHisto("RICHRejectionFactor_mmsq_num", mmsq);
if (Calo->Rejection()==0)return;
fUserMethods->FillHisto("RICHRejectionFactor_mmsq_calo", mmsq);
fUserMethods->FillHisto("RICHRejectionFactor_mmsq_p_muon",PionMomentum.Mag(),mmsq);
}

void TrackAnalysis::PionKineFactor(){
Double_t mmsq = One->Getmmsq();
 TVector3 PionMomentum = One->GetPionMomentum();
fUserMethods->FillHisto("PionKineFactor_mmsq_onetrack", mmsq);
if (Calo->Rejection()==0)return;
fUserMethods->FillHisto("PionKineFactor_mmsq_calo", mmsq);
if (One->GetRICHmass()<0.13)return;
fUserMethods->FillHisto("PionKineFactor_mmsq_RICH", mmsq);
if(Photon->OtherVetoes(LAVEvent, IRCEvent, SACEvent, SAVEvent, One->GetReft(), MCflag)==0)return;
fUserMethods->FillHisto("PionKineFactor_mmsq_photon", mmsq);
if (LKrEvent->GetNCandidates()>3)return;
TRecoLKrCandidate* LKrCand = One->GetLKrCand();
if(LKrCand->GetClusterEnergy()<2)return;
if (PionMomentum.Mag()>15&&PionMomentum.Mag()<35)fUserMethods->FillHisto("PionKineFactor_mmsq_den", mmsq);
if ((PionMomentum.Mag()>15&&PionMomentum.Mag()<35)&&((mmsq>0&&mmsq<0.01)||(mmsq>0.026&&mmsq<0.068)))fUserMethods->FillHisto("PionKineFactor_mmsq_num", mmsq);
}

void TrackAnalysis::RICHAnalysisPion(){
if (L0->GetTriggerFlags()!=32)return; //this is for !MUV3 trigger efficiency
Double_t mmsq = One->Getmmsq();
TVector3 PionMomentum = One->GetPionMomentum();
TVector3 PionPosition = One->GetPionPosition();
Parameters* P = new Parameters(3809); 
LKra->Set(LKrEvent, One->GetReft(), P->GettLKrLow(), P->GettLKrHi(), MCflag);
if(LKra->OneNeutralPion_Straws(One->GetPionMomentum(),One->GetPionPosition(), One->Giga->GetKaonFourMomentum(), 0)==0)return;
Int_t index = LKra->GetOneNeutralPion_Straws_minID();
if (index ==-1)return;
Int_t index2 = One->LKra->GetSinglePionCluster_minID();
if (index!=index2)return; 
//LKrCand = (TRecoLKrCandidate*)LKrEvent->GetCandidate(index);
Double_t mmsq_LKr = LKra->GetPiPlus_mmsq();
if(Photon->OtherVetoes(LAVEvent, IRCEvent, SACEvent, SAVEvent, One->GetReft(), MCflag)==0)return; 
fUserMethods->FillHisto("RICHAnalysisPion_mmsq_LKr", mmsq_LKr);
if(mmsq<0.01 || mmsq>0.026)return;
//if(mmsq<0.01 || mmsq>0.026)return; //took out after g's plots
fUserMethods->FillHisto("RICHAnalysisPion_mmsq_LKr_all", mmsq_LKr);
if(mmsq_LKr<0 || mmsq_LKr>.04)return;
//if(mmsq<0.01 || mmsq>0.026)return; //took out after g's plots
if (MUV3a->Veto()==1)fUserMethods->FillHisto("RICHAnalysisPion_mmsq_pion_selection", mmsq);
if (MUV3a->VetoL()==1)fUserMethods->FillHisto("RICHAnalysisPion_mmsq_pion_selection2", mmsq);
fUserMethods->FillHisto("PiPi0_mmissVSP_kine_3gtk_2",PionMomentum.Mag(),mmsq);
if (PionMomentum.Mag()>15&&PionMomentum.Mag()<35)fUserMethods->FillHisto("RICHAnalysisPion_mmsq_den", mmsq);
//if (Calo->Rejection()==0)return;
fUserMethods->FillHisto("PiPi0_mmissVSP_kine_3gtk_3",PionMomentum.Mag(),mmsq);
if (PionMomentum.Mag()>15&&PionMomentum.Mag()<35)fUserMethods->FillHisto("RICHAnalysisPion_mmsq_num_calo", mmsq);
fUserMethods->FillHisto("RICHAnalysisPion_RICHmass", One->GetRICHmass()); 
if(One->GetPionMomentum().Mag()>=15&&One->GetPionMomentum().Mag()<35)fUserMethods->FillHisto("RICHAnalysisPion_RICHmass_signal", One->GetRICHmass());
Int_t bin = 10; 
if (One->GetPionMomentum().Mag()<15)bin=0;
else if (One->GetPionMomentum().Mag()>=15&&One->GetPionMomentum().Mag()<20)bin=1;
else if (One->GetPionMomentum().Mag()>=20&&One->GetPionMomentum().Mag()<25)bin=2;
else if (One->GetPionMomentum().Mag()>=25&&One->GetPionMomentum().Mag()<30)bin=3;
else if (One->GetPionMomentum().Mag()>=30&&One->GetPionMomentum().Mag()<35)bin=4;
else bin=5; 
for(int i=0; i<6; i++){
if(i==bin)fUserMethods->FillHisto(pp_bin[i].Data(), One->GetRICHmass());
}
if (One->GetRICHmass()<0.13)return; 
///to use, remove Calo->Rejection function above, make sure not using other calo plot generator in main program, this is for pion selection,
//for muon see CalorimeterPlots function all the way at the bottom(could have stuck in RICH analysis muon function, don't know why i didn't.
CaloPlot->Set(MUV1Event, MUV2Event, MUV3Event, One->GetLKrCand(), One->GetReft(), One->Getmmsq(), One->GetPionPosition(), One->GetPionMomentum(), MCflag);
CaloPlot->Rejection();

fUserMethods->FillHisto("PiPi0_mmissVSP_kine_3gtk",PionMomentum.Mag(),mmsq);
fUserMethods->FillHisto("PiPi0_mmissVSP_kine_nogtk",PionMomentum.Mag(),mmsq);
for(int i=0; i<6; i++){
if(i==bin)fUserMethods->FillHisto(pp_bin_post[i].Data(), One->GetRICHmass());
}
if(One->GetPionMomentum().Mag()>=15&&One->GetPionMomentum().Mag()<35)fUserMethods->FillHisto("RICHAnalysisPion_RICHmass_signal_post", One->GetRICHmass());
if (PionMomentum.Mag()>15&&PionMomentum.Mag()<35)fUserMethods->FillHisto("RICHAnalysisPion_mmsq_num", mmsq);
}

void TrackAnalysis::RICHAnalysisMuon(){
Parameters* P = new Parameters(3809);  
Double_t mmsq = One->Getmmsq();
TVector3 PionMomentum = One->GetPionMomentum(); 
TVector3 PionPosition = One->GetPionPosition(); 
fUserMethods->FillHisto("RICHAnalysisMuon_mmsq_onetrack", mmsq);
if (mmsq>0||mmsq<-0.02)return;// take out for MC
TLorentzVector MuMom(PionMomentum.Px(), PionMomentum.Py(), PionMomentum.Pz(),sqrt(PionMomentum*PionMomentum+(105.66/1000)*(105.66/1000)));
Photon->Set(LKrEvent, LAVEvent, IRCEvent, SACEvent, SAVEvent, One->GetLKrCand(), One->GetReft(), One->Getmmsq(), One->LKra->GetSinglePionCluster_unmatched(), One->GetLKrPosition(), MCflag);
if(Photon->Rejection()==0)return;
TLorentzVector NeutrinoMom = One->Giga->GetKaonFourMomentum() - MuMom;
Double_t mmsq_n = NeutrinoMom.M2();
fUserMethods->FillHisto("RICHAnalysisMuon_mmsq_n", mmsq_n);
//if(mmsq_n<-0.004||mmsq_n>0.004)return;
fUserMethods->FillHisto("RICHAnalysisMuon_mmsq_n_post", mmsq_n);
TVector3 MUV3_pos = prop(PionPosition, PionMomentum*1000, P->GetMUV3StartPos()); 
if (MUV3_pos.Y()>1220||MUV3_pos.Y()<-1220)return; //from what i can tell, nominal parameter, say +/-1320 for outer, +/- 76 for inner (x and y)
if (MUV3_pos.Y()<100&&MUV3_pos.Y()>-100)return;
if (MUV3_pos.X()>1220||MUV3_pos.X()<-1220)return;
if (MUV3_pos.X()<100&&MUV3_pos.X()>-100)return;
MUV3a->Set(MUV3Event, One->GetReft(), P->GettMUV3Low(), P->GettMUV3Hi(), MCflag);
if (MUV3a->Veto()==1)return;
 fUserMethods->FillHisto("RICHAnalysisMuon_mmsq_n_all", mmsq_n);
if(mmsq_n<-0.004||mmsq_n>0.004)return;
fUserMethods->FillHisto("RICHAnalysisMuon_mmsq_MuonSelection", mmsq);
fUserMethods->FillHisto("RICHAnalysisMuon_mmsq_p_selection",PionMomentum.Mag(),mmsq);
//////////////end muon selection///////////////////
fUserMethods->FillHisto("RICHAnalysisMuon_RICHmass", One->GetRICHmass()); 
if(One->GetPionMomentum().Mag()>=15&&One->GetPionMomentum().Mag()<35)fUserMethods->FillHisto("RICHAnalysisMuon_RICHmass_signal", One->GetRICHmass());
Int_t bin = 10; 
if (One->GetPionMomentum().Mag()<15)bin=0;
else if (One->GetPionMomentum().Mag()>=15&&One->GetPionMomentum().Mag()<20)bin=1;
else if (One->GetPionMomentum().Mag()>=20&&One->GetPionMomentum().Mag()<25)bin=2;
else if (One->GetPionMomentum().Mag()>=25&&One->GetPionMomentum().Mag()<30)bin=3;
else if (One->GetPionMomentum().Mag()>=30&&One->GetPionMomentum().Mag()<35)bin=4;
else bin=5; 
for(int i=0; i<6; i++){
if(i==bin)fUserMethods->FillHisto(mp_bin[i].Data(), One->GetRICHmass());
}
if (One->GetRICHmass()>0.13)return; 
for(int i=0; i<6; i++){
if(i==bin)fUserMethods->FillHisto(mp_bin_post[i].Data(), One->GetRICHmass());
}
if(One->GetPionMomentum().Mag()>=15&&One->GetPionMomentum().Mag()<35)fUserMethods->FillHisto("RICHAnalysisMuon_RICHmass_signal_post", One->GetRICHmass());

//fUserMethods->FillHisto("RICHAnalysisMuon_mmsq_p_selection",PionMomentum.Mag(),mmsq);
}

void TrackAnalysis::MuonAccidentalPhotonFactor(){ 
Parameters* P = new Parameters(3809);  
Double_t mmsq = One->Getmmsq();
TVector3 PionMomentum = One->GetPionMomentum(); 
TVector3 PionPosition = One->GetPionPosition(); 
fUserMethods->FillHisto("MuonAccidentalPhotonFactor_mmsq_onetrack", mmsq);
//if (mmsq>0||mmsq<-0.02)return; 
TLorentzVector MuMom(PionMomentum.Px(), PionMomentum.Py(), PionMomentum.Pz(),sqrt(PionMomentum*PionMomentum+(105.66/1000)*(105.66/1000)));
Photon->Set(LKrEvent, LAVEvent, IRCEvent, SACEvent, SAVEvent, One->GetLKrCand(), One->GetReft(), One->Getmmsq(), One->LKra->GetSinglePionCluster_unmatched(), One->GetLKrPosition(), MCflag);
TLorentzVector NeutrinoMom = One->Giga->GetKaonFourMomentum() - MuMom;
Double_t mmsq_n = NeutrinoMom.M2();
fUserMethods->FillHisto("MuonAccidentalPhotonFactor_mmsq_n", mmsq_n);
if(mmsq_n<-0.005||mmsq_n>0.005)return;//changed from 0.004 to match kaon integral calculation
fUserMethods->FillHisto("MuonAccidentalPhotonFactor_mmsq_n_post", mmsq_n);
TVector3 MUV3_pos = prop(PionPosition, PionMomentum*1000, P->GetMUV3StartPos()); 
if (MUV3_pos.Y()>1220||MUV3_pos.Y()<-1220)return; //from what i can tell, nominal parameter, say +/-1320 for outer, +/- 76 for inner (x and y)
if (MUV3_pos.Y()<100&&MUV3_pos.Y()>-100)return;
if (MUV3_pos.X()>1220||MUV3_pos.X()<-1220)return;
if (MUV3_pos.X()<100&&MUV3_pos.X()>-100)return;
if (Calo->MuonSelection()==0)return;
fUserMethods->FillHisto("MuonAccidentalPhotonFactor_mmsq_MuonSelection", mmsq);
//////////////end muon selection///////////////////
fUserMethods->FillHisto("MuonAccidentalPhotonFactor_RICHmass", One->GetRICHmass()); 
if (One->GetRICHmass()>0.13)return; 
fUserMethods->FillHisto("MuonAccidentalPhotonFactor_mmsq_den", mmsq);
if(Photon->Rejection()==0)return;
fUserMethods->FillHisto("MuonAccidentalPhotonFactor_mmsq_num", mmsq);
}

void TrackAnalysis::MuonKineFactor(){
TVector3 PionMomentum = One->GetPionMomentum();
Double_t mmsq = One->Getmmsq();
fUserMethods->FillHisto("MuonKineFactor_mmsq_onetrack", mmsq);
if(Photon->Rejection()==0)return;
fUserMethods->FillHisto("MuonKineFactor_mmsq_photon", mmsq);
if (One->GetRICHmass()>0.13)return;
fUserMethods->FillHisto("MuonKineFactor_mmsq_RICH", mmsq);
if (Calo->MuonSelection()==0)return;
fUserMethods->FillHisto("MuonKineFactor_mmsq_calo", mmsq);
if (PionMomentum.Mag()>15&&PionMomentum.Mag()<35)fUserMethods->FillHisto("MuonKineFactor_mmsq_den", mmsq);
if (PionMomentum.Mag()>15&&PionMomentum.Mag()<35&&mmsq>0&&mmsq<0.01)fUserMethods->FillHisto("MuonKineFactor_mmsq_num", mmsq);
}

void TrackAnalysis::BeamAccidentalPhotonFactor(){
Double_t mmsq = One->Getmmsq();
fUserMethods->FillHisto("BeamAccidentalPhotonFactor_mmsq_onetrack", mmsq);
if (Calo->Rejection()==0)return;
fUserMethods->FillHisto("BeamAccidentalPhotonFactor_mmsq_calo", mmsq);
if (One->GetRICHmass()<0.13)return;
fUserMethods->FillHisto("BeamAccidentalPhotonFactor_mmsq_RICH", mmsq);
/*if (PionMomentum.Mag()>15&&PionMomentum.Mag()<35)*/fUserMethods->FillHisto("BeamAccidentalPhotonFactor_mmsq_den", mmsq);
if(Photon->Rejection()==0)return;
fUserMethods->FillHisto("BeamAccidentalPhotonFactor_mmsq_photon", mmsq);
/*if (PionMomentum.Mag()>15&&PionMomentum.Mag()<35)*/fUserMethods->FillHisto("BeamAccidentalPhotonFactor_mmsq_num", mmsq);
}

void TrackAnalysis::ep(){
Double_t mmsq = One->Getmmsq();
TVector3 PionMomentum = One->GetPionMomentum();
TVector3 PionPosition = One->GetPionPosition();
Parameters* P = new Parameters(3809);
//(!MCflag){ 
LKra->Set(LKrEvent, One->GetReft(), P->GettLKrLow(), P->GettLKrHi(), MCflag);
if(Photon->OtherVetoes(LAVEvent, IRCEvent, SACEvent, SAVEvent, One->GetReft(), MCflag)==0)return;
//if(!MCflag){
if(LKra->OneNeutralPion_Straws(One->GetPionMomentum(),One->GetPionPosition(), One->Giga->GetKaonFourMomentum(), 0)==0)return;
Int_t index = LKra->GetOneNeutralPion_Straws_minID();
if (index ==-1)return;
Int_t index2 = One->LKra->GetSinglePionCluster_minID();
if (index!=index2)return; 
//LKrCand = (TRecoLKrCandidate*)LKrEvent->GetCandidate(index);
Double_t mmsq_LKr = LKra->GetPiPlus_mmsq();
//if(Photon->OtherVetoes(LAVEvent, IRCEvent, SACEvent, SAVEvent, One->GetReft(), MCflag)==0)return;//now put in for everything 
//fUserMethods->FillHisto("RICHAnalysisPion_mmsq_LKr", mmsq_LKr);
if(mmsq<0.01 || mmsq>0.026)return;
//if(mmsq<0.01 || mmsq>0.026)return; //took out after g's plots
//fUserMethods->FillHisto("RICHAnalysisPion_mmsq_LKr_all", mmsq_LKr);
fUserMethods->FillHisto("ep_mmsq_LKr_pre", mmsq_LKr);
if(mmsq_LKr<0 || mmsq_LKr>.04)return;
//if(mmsq<0.01 || mmsq>0.026)return; //took out after g's plots
//fUserMethods->FillHisto("RICHAnalysisPion_mmsq_pion_selection", mmsq);
//if (PionMomentum.Mag()>15&&PionMomentum.Mag()<35)fUserMethods->FillHisto("RICHAnalysisPion_mmsq_den", mmsq);
//if (Calo->Rejection()==0)return;
//if (PionMomentum.Mag()>15&&PionMomentum.Mag()<35)fUserMethods->FillHisto("RICHAnalysisPion_mmsq_num_calo", mmsq);
//fUserMethods->FillHisto("RICHAnalysisPion_RICHmass", One->GetRICHmass()); 
//if(One->GetPionMomentum().Mag()>=15&&One->GetPionMomentum().Mag()<35)fUserMethods->FillHisto("RICHAnalysisPion_RICHmass_signal", One->GetRICHmass());
fUserMethods->FillHisto("ep_mmsq_LKr", mmsq_LKr);
//}//put in for pi+pi0 mc, take out for pinn, put in for pip0 purity study

Int_t bin = 10; 
if (One->GetPionMomentum().Mag()<15)bin=0;
else if (One->GetPionMomentum().Mag()>=15&&One->GetPionMomentum().Mag()<20)bin=1;
else if (One->GetPionMomentum().Mag()>=20&&One->GetPionMomentum().Mag()<25)bin=2;
else if (One->GetPionMomentum().Mag()>=25&&One->GetPionMomentum().Mag()<30)bin=3;
else if (One->GetPionMomentum().Mag()>=30&&One->GetPionMomentum().Mag()<35)bin=4;
else bin=5; 
for(int i=0; i<6; i++){
if(i==bin)fUserMethods->FillHisto(epRICH[i].Data(), mmsq);
if(i==bin)fUserMethods->FillHisto(epCal[i].Data(), mmsq); 
}
if (One->GetRICHmass()>0.13){
 for(int i=0; i<6; i++){
  if(i==bin)fUserMethods->FillHisto(epRICH_post[i].Data(), mmsq);
 }
}
//if (One->GetRICHmass()<0.13)return; 
CaloPlot->Set(MUV1Event, MUV2Event, MUV3Event, One->GetLKrCand(), One->GetReft(), One->Getmmsq(), One->GetPionPosition(), One->GetPionMomentum(), MCflag);  
CaloPlot->Rejection();

if (Calo->Rejection()==1){
 for(int i=0; i<6; i++){
  if(i==bin)fUserMethods->FillHisto(epCal_post[i].Data(), mmsq);
 }
}
}

void TrackAnalysis::MUV3Efficiency(){
Parameters* P = new Parameters(3809);  
Double_t mmsq = One->Getmmsq();
TVector3 PionMomentum = One->GetPionMomentum(); 
TVector3 PionPosition = One->GetPionPosition(); 
//if (mmsq>0||mmsq<-0.02)return; // leave out for now
TLorentzVector MuMom(PionMomentum.Px(), PionMomentum.Py(), PionMomentum.Pz(),sqrt(PionMomentum*PionMomentum+(105.66/1000)*(105.66/1000)));
Photon->Set(LKrEvent, LAVEvent, IRCEvent, SACEvent, SAVEvent, One->GetLKrCand(), One->GetReft(), One->Getmmsq(), One->LKra->GetSinglePionCluster_unmatched(), One->GetLKrPosition(), MCflag);
if(Photon->Rejection()==0)return;
TLorentzVector NeutrinoMom = One->Giga->GetKaonFourMomentum() - MuMom;
Double_t mmsq_n = NeutrinoMom.M2();
if(mmsq_n<-0.004||mmsq_n>0.004)return;
fUserMethods->FillHisto("RICHAnalysisMuon_mmsq_n_post", mmsq_n);
TVector3 MUV3_pos = prop(PionPosition, PionMomentum*1000, P->GetMUV3StartPos()); 
if (MUV3_pos.Y()>1220||MUV3_pos.Y()<-1220)return; //from what i can tell, nominal parameter, say +/-1320 for outer, +/- 76 for inner (x and y)
if (MUV3_pos.Y()<100&&MUV3_pos.Y()>-100)return;
if (MUV3_pos.X()>1220||MUV3_pos.X()<-1220)return;
if (MUV3_pos.X()<100&&MUV3_pos.X()>-100)return;
MUV3a->Set(MUV3Event, One->GetReft(), P->GettMUV3Low(), P->GettMUV3Hi(), MCflag);
//////////////end muon selection/////////////////// 
//if(!MCflag){
if (One->GetRICHmass()>0.13)return; 
if (Calo->MUV3Efficiency()==0)return;
//}
fUserMethods->FillHisto("MUV3Efficiency_den", mmsq_n); 
if (MUV3a->Veto()==1) return; 
fUserMethods->FillHisto("MUV3Efficiency_num", mmsq_n);
}

void TrackAnalysis::Kmu2Normalization(){
if(!MCflag){if (L0->GetTriggerFlags()!=32)return;} 
Parameters* P = new Parameters(3809);  
Double_t mmsq = One->Getmmsq();
TVector3 PionMomentumBeforeMagnet = One->GetPionMomentumBeforeMagnet();
TVector3 PionMomentum = One->GetPionMomentum(); 
TVector3 PionPosition = One->GetPionPosition(); 
//if (mmsq>0||mmsq<-0.02)return; // leave out for now
TLorentzVector MuMom(PionMomentumBeforeMagnet.Px(), PionMomentumBeforeMagnet.Py(), PionMomentumBeforeMagnet.Pz(),sqrt(PionMomentumBeforeMagnet*PionMomentumBeforeMagnet+(105.66/1000)*(105.66/1000)));
Photon->Set(LKrEvent, LAVEvent, IRCEvent, SACEvent, SAVEvent, One->GetLKrCand(), One->GetReft(), One->Getmmsq(), One->LKra->GetSinglePionCluster_unmatched(), One->GetLKrPosition(), MCflag);
if(MCflag)fUserMethods->FillHisto("Kmu2Normalization_photonveto_counter", mmsq);
if(!MCflag){if(Photon->Rejection()==0)return;}//NOT SURE about this! maybe the same as usual, MC and data
if(MCflag&&(Photon->Rejection()==0))fUserMethods->FillHisto("Kmu2Normalization_photonveto_counter_post", mmsq); 
TLorentzVector NeutrinoMom = One->Giga->GetKaonFourMomentum() - MuMom;
Double_t mmsq_n = NeutrinoMom.M2();
//if(mmsq_n<-0.004||mmsq_n>0.004)return;
TVector3 MUV3_pos = prop(PionPosition, PionMomentum*1000, P->GetMUV3StartPos()); 
if (MUV3_pos.Y()>1220||MUV3_pos.Y()<-1220)return; //from what i can tell, nominal parameter, say +/-1320 for outer, +/- 76 for inner (x and y)
if (MUV3_pos.Y()<100&&MUV3_pos.Y()>-100)return;
if (MUV3_pos.X()>1220||MUV3_pos.X()<-1220)return;
if (MUV3_pos.X()<100&&MUV3_pos.X()>-100)return;
MUV3a->Set(MUV3Event, One->GetReft(), P->GettMUV3Low(), P->GettMUV3Hi(), MCflag);
if (MUV3a->Veto()==1)return; 
fUserMethods->FillHisto("Kmu2Normalization_number_before", mmsq_n);
if(mmsq_n<-0.005||mmsq_n>0.005)return;
//if(mmsq_n<-0.004||mmsq_n>0.004)return;//changed to compare to .005

fUserMethods->FillHisto("Kmu2Normalization_number", mmsq_n);
fUserMethods->FillHisto("Kmu2Normalization_count", PionMomentum.Mag(), mmsq);
}

void TrackAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
LKra->SaveAllPlots();
CaloPlot->SaveAllPlots();
}

void TrackAnalysis::BookHistos(){
	LKra->BookHistos_Set();
	MUV3a->BookHistos_Set(); 
	LKra->BookHistos_OneNeutralPion_Straws();  
	CaloPlot->BookHistos(); 
        fUserMethods->BookHisto(new TH1F("PhotonRejectionFactor_mmsq_den", "Missing Mass squared before rej", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("PhotonRejectionFactor_mmsq_num", "Missing Mass squared after rej", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("PhotonRejectionFactor_mmsq_onetrack", "Missing Mass squared, one track selection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("PhotonRejectionFactor_mmsq_calo", "Missing Mass squared, calo rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("PhotonRejectionFactor_mmsq_RICH", "Missing Mass squared, RICH rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("PhotonRejectionFactor_mmsq_photon", "Missing Mass squared, photon rejection", 250, -.35, .15));
	fUserMethods->BookHisto(new TH1F("PhotonRejectionFactor_mmsq_end", "Missing Mass squared, photon rejection", 250, -.35, .15));
	fUserMethods->BookHisto(new TH1F("PhotonRejectionFactor_mmsq_end_signal", "Missing Mass squared, photon rejection", 250, -.35, .15));
TH2I* PhotonRejectionFactor_mmsq_p_onetrack = new TH2I("PhotonRejectionFactor_mmsq_p_onetrack", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
PhotonRejectionFactor_mmsq_p_onetrack->SetOption("COLZ");
PhotonRejectionFactor_mmsq_p_onetrack->Draw();
fUserMethods->BookHisto(PhotonRejectionFactor_mmsq_p_onetrack);
TH2I* PhotonRejectionFactor_mmsq_p_photon = new TH2I("PhotonRejectionFactor_mmsq_p_photon", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
PhotonRejectionFactor_mmsq_p_photon->SetOption("COLZ");
PhotonRejectionFactor_mmsq_p_photon->Draw();
fUserMethods->BookHisto(PhotonRejectionFactor_mmsq_p_photon);
TH2I* PhotonRejectionFactor_mmsq_p_muon = new TH2I("PhotonRejectionFactor_mmsq_p_muon", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
PhotonRejectionFactor_mmsq_p_muon->SetOption("COLZ");
PhotonRejectionFactor_mmsq_p_muon->Draw();
fUserMethods->BookHisto(PhotonRejectionFactor_mmsq_p_muon);

TH2I* PiPi0_mmissVSP_kine_3gtk = new TH2I("PiPi0_mmissVSP_kine_3gtk", "momentum (x-axis), mmiss (y-axis)", 200, 0, 100, 375, -0.1, 0.14375);
PiPi0_mmissVSP_kine_3gtk->SetOption("COLZ");
PiPi0_mmissVSP_kine_3gtk->Draw();
fUserMethods->BookHisto(PiPi0_mmissVSP_kine_3gtk);

TH2I* PiPi0_mmissVSP_kine_nogtk = new TH2I("PiPi0_mmissVSP_kine_nogtk", "momentum (x-axis), mmiss (y-axis)", 200, 0, 100, 375, -0.1, 0.14375);
PiPi0_mmissVSP_kine_nogtk->SetOption("COLZ");
PiPi0_mmissVSP_kine_nogtk->Draw();
fUserMethods->BookHisto(PiPi0_mmissVSP_kine_nogtk);

TH2I* PiPi0_mmissVSP_kine_3gtk_2 = new TH2I("PiPi0_mmissVSP_kine_3gtk_2", "momentum (x-axis), mmiss (y-axis)", 200, 0, 100, 375, -0.1, 0.14375);
PiPi0_mmissVSP_kine_3gtk_2->SetOption("COLZ");
PiPi0_mmissVSP_kine_3gtk_2->Draw();
fUserMethods->BookHisto(PiPi0_mmissVSP_kine_3gtk_2);

TH2I* PiPi0_mmissVSP_kine_3gtk_3 = new TH2I("PiPi0_mmissVSP_kine_3gtk_3", "momentum (x-axis), mmiss (y-axis)", 200, 0, 100, 375, -0.1, 0.14375);
PiPi0_mmissVSP_kine_3gtk_3->SetOption("COLZ");
PiPi0_mmissVSP_kine_3gtk_3->Draw();
fUserMethods->BookHisto(PiPi0_mmissVSP_kine_3gtk_3);


	fUserMethods->BookHisto(new TH1F("CalorimeterRejectionFactor_mmsq_den", "Missing Mass squared before rej", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("CalorimeterRejectionFactor_mmsq_num", "Missing Mass squared after rej", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("CalorimeterRejectionFactor_mmsq_onetrack", "Missing Mass squared, one track selection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("CalorimeterRejectionFactor_mmsq_calo", "Missing Mass squared, calo rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("CalorimeterRejectionFactor_mmsq_RICH", "Missing Mass squared, RICH rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("CalorimeterRejectionFactor_mmsq_photon", "Missing Mass squared, Calorimeter rejection", 250, -.35, .15));
TH2I* CalorimeterRejectionFactor_mmsq_p_onetrack = new TH2I("CalorimeterRejectionFactor_mmsq_p_onetrack", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
CalorimeterRejectionFactor_mmsq_p_onetrack->SetOption("COLZ");
CalorimeterRejectionFactor_mmsq_p_onetrack->Draw();
fUserMethods->BookHisto(CalorimeterRejectionFactor_mmsq_p_onetrack);
TH2I* CalorimeterRejectionFactor_mmsq_p_photon = new TH2I("CalorimeterRejectionFactor_mmsq_p_photon", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
CalorimeterRejectionFactor_mmsq_p_photon->SetOption("COLZ");
CalorimeterRejectionFactor_mmsq_p_photon->Draw();
fUserMethods->BookHisto(CalorimeterRejectionFactor_mmsq_p_photon);
TH2I* CalorimeterRejectionFactor_mmsq_p_RICH = new TH2I("CalorimeterRejectionFactor_mmsq_p_RICH", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
CalorimeterRejectionFactor_mmsq_p_RICH->SetOption("COLZ");
CalorimeterRejectionFactor_mmsq_p_RICH->Draw();
fUserMethods->BookHisto(CalorimeterRejectionFactor_mmsq_p_RICH);
TH2I* CalorimeterRejectionFactor_mmsq_p_calo = new TH2I("CalorimeterRejectionFactor_mmsq_p_calo", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
CalorimeterRejectionFactor_mmsq_p_calo->SetOption("COLZ");
CalorimeterRejectionFactor_mmsq_p_calo->Draw();
fUserMethods->BookHisto(CalorimeterRejectionFactor_mmsq_p_calo);
	fUserMethods->BookHisto(new TH1F("MuonRejectionFactor_mmsq_end", "Missing Mass squared, photon rejection", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("MuonRejectionFactor_mmsq_end_signal", "Missing Mass squared, photon rejection", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("RICHRejectionFactor_mmsq_den", "Missing Mass squared before rej", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHRejectionFactor_mmsq_num", "Missing Mass squared after rej", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHRejectionFactor_mmsq_onetrack", "Missing Mass squared, one track selection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHRejectionFactor_mmsq_calo", "Missing Mass squared, calo rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHRejectionFactor_mmsq_photon", "Missing Mass squared, RICH rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHRejectionFactor_mmsq_RICH", "Missing Mass squared, RICH rejection", 250, -.35, .15));
TH2I* RICHRejectionFactor_mmsq_p_onetrack = new TH2I("RICHRejectionFactor_mmsq_p_onetrack", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
RICHRejectionFactor_mmsq_p_onetrack->SetOption("COLZ");
RICHRejectionFactor_mmsq_p_onetrack->Draw();
fUserMethods->BookHisto(RICHRejectionFactor_mmsq_p_onetrack);

TH2I* RICHRejectionFactor_mmsq_p_muon = new TH2I("RICHRejectionFactor_mmsq_p_muon", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
RICHRejectionFactor_mmsq_p_muon->SetOption("COLZ");
RICHRejectionFactor_mmsq_p_muon->Draw();
fUserMethods->BookHisto(RICHRejectionFactor_mmsq_p_muon);
TH2I* RICHRejectionFactor_mmsq_p_photon = new TH2I("RICHRejectionFactor_mmsq_p_photon", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
RICHRejectionFactor_mmsq_p_photon->SetOption("COLZ");
RICHRejectionFactor_mmsq_p_photon->Draw();
fUserMethods->BookHisto(RICHRejectionFactor_mmsq_p_photon);
TH2I* RICHRejectionFactor_mmsq_p_RICH = new TH2I("RICHRejectionFactor_mmsq_p_RICH", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
RICHRejectionFactor_mmsq_p_RICH->SetOption("COLZ");
RICHRejectionFactor_mmsq_p_RICH->Draw();
fUserMethods->BookHisto(RICHRejectionFactor_mmsq_p_RICH);
	fUserMethods->BookHisto(new TH1F("PionKineFactor_mmsq_den", "Missing Mass squared before rej", 1250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("PionKineFactor_mmsq_num", "Missing Mass squared after rej", 1250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("PionKineFactor_mmsq_onetrack", "Missing Mass squared, one track selection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("PionKineFactor_mmsq_calo", "Missing Mass squared, calo rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("PionKineFactor_mmsq_RICH", "Missing Mass squared, RICH rejection ", 250, -.35, .15));
 fUserMethods->BookHisto(new TH1F("PionKineFactor_mmsq_photon", "Missing Mass squared, RICH rejection", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("MuonKineFactor_mmsq_den", "Missing Mass squared before rej", 1250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("MuonKineFactor_mmsq_num", "Missing Mass squared after rej", 1250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("MuonKineFactor_mmsq_onetrack", "Missing Mass squared, one track selection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("MuonKineFactor_mmsq_calo", "Missing Mass squared, calo rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("MuonKineFactor_mmsq_photon", "Missing Mass squared, RICH rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("MuonKineFactor_mmsq_RICH", "Missing Mass squared, RICH rejection", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("BeamAccidentalPhotonFactor_mmsq_den", "Missing Mass squared before rej", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("BeamAccidentalPhotonFactor_mmsq_num", "Missing Mass squared after rej", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("BeamAccidentalPhotonFactor_mmsq_onetrack", "Missing Mass squared, one track selection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("BeamAccidentalPhotonFactor_mmsq_calo", "Missing Mass squared, calo rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("BeamAccidentalPhotonFactor_mmsq_RICH", "Missing Mass squared, RICH rejection ", 250, -.35, .15));
      fUserMethods->BookHisto(new TH1F("BeamAccidentalPhotonFactor_mmsq_photon", "Missing Mass squared, RICH rejection", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("MuonAccidentalPhotonFactor_mmsq_den", "Missing Mass squared before rej", 500, -.35, .15));//changed to accomadte .005
        fUserMethods->BookHisto(new TH1F("MuonAccidentalPhotonFactor_mmsq_num", "Missing Mass squared after rej", 500, -.35, .15));
        fUserMethods->BookHisto(new TH1F("MuonAccidentalPhotonFactor_mmsq_onetrack", "Missing Mass squared, one track selection ", 500, -.35, .15));
        fUserMethods->BookHisto(new TH1F("MuonAccidentalPhotonFactor_mmsq_calo", "Missing Mass squared, calo rejection ", 500, -.35, .15));
        fUserMethods->BookHisto(new TH1F("MuonAccidentalPhotonFactor_mmsq_RICH", "Missing Mass squared, RICH rejection ", 500, -.35, .15));
      fUserMethods->BookHisto(new TH1F("MuonAccidentalPhotonFactor_mmsq_photon", "Missing Mass squared, RICH rejection", 500, -.35, .15));
fUserMethods->BookHisto(new TH1F("MuonAccidentalPhotonFactor_mmsq_n", "Missing Mass squared", 1000, -.35, .15));
fUserMethods->BookHisto(new TH1F("MuonAccidentalPhotonFactor_mmsq_n_post", "Missing Mass squared", 500, -.35, .15));
fUserMethods->BookHisto(new TH1F("MUV3Efficiency_den", "MUV3Efficiency Denominator", 1000, -.35, .15));
fUserMethods->BookHisto(new TH1F("MUV3Efficiency_num", "MUV3Efficiency Numerator", 1000, -.35, .15));
fUserMethods->BookHisto(new TH1F("MuonAccidentalPhotonFactor_mmsq_MuonSelection", "Missing Mass squared", 500, -.35, .15));
 fUserMethods->BookHisto(new TH1I("MuonAccidentalPhotonFactor_RICHmass", "RICH mass",  280, 0,1.4));
        fUserMethods->BookHisto(new TH1F("RICHAnalysisPion_mmsq_den", "Missing Mass squared before rej", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHAnalysisPion_mmsq_num", "Missing Mass squared after rej", 250, -.35, .15));
	fUserMethods->BookHisto(new TH1F("RICHAnalysisPion_mmsq_num_calo", "Missing Mass squared after rej", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHAnalysisPion_mmsq_onetrack", "Missing Mass squared, one track selection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHAnalysisPion_mmsq_calo", "Missing Mass squared, calo rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHAnalysisPion_mmsq_RICH", "Missing Mass squared, RICH rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHAnalysisPion_mmsq_photon", "Missing Mass squared, RICH rejection", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("RICHAnalysisPion_mmsq_LKr", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("ep_mmsq_LKr_pre", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("ep_mmsq_LKr", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("RICHAnalysisPion_mmsq_LKr_all", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("RICHAnalysisPion_mmsq_pion_selection", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("RICHAnalysisPion_mmsq_pion_selection2", "Missing Mass squared", 250, -.35, .15));
 fUserMethods->BookHisto(new TH1I("RICHAnalysisPion_RICHmass", "RICH mass",  280, 0,1.4));
 fUserMethods->BookHisto(new TH1I("RICHAnalysisPion_RICHmass_signal", "RICH mass",  280, 0,1.4));
fUserMethods->BookHisto(new TH1I("RICHAnalysisPion_RICHmass_signal_post", "RICH mass",  280, 0,1.4));
for(int i=0; i<6; i++){
  pp_bin[i].Form("RICHAnalysisPion_RICHmass_%d", i);
  pp_bin_post[i].Form("RICHAnalysisPion_RICHmass_%d_post", i);
  fUserMethods->BookHisto(new TH1I(pp_bin[i].Data(), "RICH mass (muons)", 450, 0, 0.45));
  fUserMethods->BookHisto(new TH1I(pp_bin_post[i].Data(), "RICH mass (muons) post cut", 450, 0, 0.45));
}
        fUserMethods->BookHisto(new TH1F("RICHAnalysisMuon_mmsq_den", "Missing Mass squared before rej", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHAnalysisMuon_mmsq_num", "Missing Mass squared after rej", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHAnalysisMuon_mmsq_onetrack", "Missing Mass squared, one track selection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHAnalysisMuon_mmsq_calo", "Missing Mass squared, calo rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHAnalysisMuon_mmsq_RICH", "Missing Mass squared, RICH rejection ", 250, -.35, .15));
        fUserMethods->BookHisto(new TH1F("RICHAnalysisMuon_mmsq_photon", "Missing Mass squared, RICH rejection", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("RICHAnalysisMuon_mmsq_n", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("RICHAnalysisMuon_mmsq_n_all", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("RICHAnalysisMuon_mmsq_n_post", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("RICHAnalysisMuon_mmsq_MuonSelection", "Missing Mass squared", 250, -.35, .15));
 fUserMethods->BookHisto(new TH1I("RICHAnalysisMuon_RICHmass", "RICH mass",  280, 0,1.4));
 fUserMethods->BookHisto(new TH1I("RICHAnalysisMuon_RICHmass_signal", "RICH mass",  280, 0,1.4));
 fUserMethods->BookHisto(new TH1I("RICHAnalysisMuon_RICHmass_signal_post", "RICH mass",  280, 0,1.4));
for(int i=0; i<6; i++){
  mp_bin[i].Form("RICHAnalysisMuon_RICHmass_%d", i);
  mp_bin_post[i].Form("RICHAnalysisMuon_RICHmass_%d_post", i);
 fUserMethods->BookHisto(new TH1I(mp_bin[i].Data(), "RICH mass (muons)", 450, 0, 0.45));  
fUserMethods->BookHisto(new TH1I(mp_bin_post[i].Data(), "RICH mass (muons) post cut", 450, 0, 0.45));
}
TH2I* RICHAnalysisMuon_mmsq_p_selection = new TH2I("RICHAnalysisMuon_mmsq_p_selection", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
RICHAnalysisMuon_mmsq_p_selection->SetOption("COLZ");
RICHAnalysisMuon_mmsq_p_selection->Draw();
fUserMethods->BookHisto(RICHAnalysisMuon_mmsq_p_selection);
fUserMethods->BookHisto(new TH1F("CalorimeterPlotGenerator_mmsq_n", "Missing Mass squared", 1250, -.35, .15));
for(int i=0; i<6; i++){
  epCal[i].Form("epCal_%d", i);
  epCal_post[i].Form("epCal_post_%d", i);
  fUserMethods->BookHisto(new TH1I(epCal[i].Data(), "Pions before Calorimeters (mmsq)", 250, -0.35, 0.15));
  fUserMethods->BookHisto(new TH1I(epCal_post[i].Data(), "Pions after Calorimeters (mmsq)", 250, -0.35, 0.15));
}
for(int i=0; i<6; i++){
  epRICH[i].Form("epRICH_%d", i);
  epRICH_post[i].Form("epRICH_post_%d", i);
  fUserMethods->BookHisto(new TH1I(epRICH[i].Data(), "Pions before RICH (mmsq)", 250, -0.35, 0.15));
  fUserMethods->BookHisto(new TH1I(epRICH_post[i].Data(), "Pions after RICH (mmsq)", 250, -0.35, 0.15));
}
TH2I* Kmu2Normalization_count = new TH2I("Kmu2Normalization_count", "momentum (x-axis), mmsq (y-axis)", 100, 0, 100, 125, -0.1, 0.15);
Kmu2Normalization_count->SetOption("COLZ");
Kmu2Normalization_count->Draw();
fUserMethods->BookHisto(Kmu2Normalization_count);
fUserMethods->BookHisto(new TH1F("Kmu2Normalization_photonveto_counter", "Missing Mass squared", 500, -.35, .15));//changed to accomdate .005
fUserMethods->BookHisto(new TH1F("Kmu2Normalization_photonveto_counter_post", "Missing Mass squared before rej", 500, -.35, .15));
fUserMethods->BookHisto(new TH1F("Kmu2Normalization_number_before", "Number before", 500, -.35, .15));
fUserMethods->BookHisto(new TH1F("Kmu2Normalization_number", "Number", 500, -.35, .15));
}


void TrackAnalysis::CalorimeterPlots(){
//if you want a muon sample for plots put cuts in, otherwise take out for muon and pions in same plots
CaloPlot->Set(MUV1Event, MUV2Event, MUV3Event, One->GetLKrCand(), One->GetReft(), One->Getmmsq(), One->GetPionPosition(), One->GetPionMomentum(), MCflag);
Parameters* P = new Parameters(3809);  
//CaloPlot->Rejection(); //move down for data
Double_t mmsq = One->Getmmsq();
if (mmsq>0||mmsq<-0.02)return;
TVector3 PionMomentum = One->GetPionMomentum(); 
TVector3 PionPosition = One->GetPionPosition(); 
//if (mmsq>0||mmsq<-0.02)return;// take out for MC
TLorentzVector MuMom(PionMomentum.Px(), PionMomentum.Py(), PionMomentum.Pz(),sqrt(PionMomentum*PionMomentum+(105.66/1000)*(105.66/1000)));
Photon->Set(LKrEvent, LAVEvent, IRCEvent, SACEvent, SAVEvent, One->GetLKrCand(), One->GetReft(), One->Getmmsq(), One->LKra->GetSinglePionCluster_unmatched(), One->GetLKrPosition(), MCflag);
if(Photon->Rejection()==0)return;
TLorentzVector NeutrinoMom = One->Giga->GetKaonFourMomentum() - MuMom;
Double_t mmsq_n = NeutrinoMom.M2();
fUserMethods->FillHisto("CalorimeterPlotGenerator_mmsq_n", mmsq_n);
if(mmsq_n<-0.004||mmsq_n>0.004)return;
//CaloPlot->Rejection();//move up for pi nu nu mc
Calo->Phil(); 
}

