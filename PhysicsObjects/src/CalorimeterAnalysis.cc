#include "CalorimeterAnalysis.hh"
#include "CHODAnalysis.hh"
#include "CedarAnalysis.hh"
#include "GigaTrackerAnalysis.hh"
#include "CHANTIAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoSpectrometerEvent.hh"
#include "Event.hh"
#include "Parameters.hh"
#include "BlueTubeTrackerAnalysis.hh"
#include "MUV1Analysis.hh"
#include "MUV2Analysis.hh"
#include "MUV3Analysis.hh"

using namespace TMath;

CalorimeterAnalysis::CalorimeterAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
MUV1a = new MUV1Analysis(ba);
MUV2a = new MUV2Analysis(ba);
MUV3a = new MUV3Analysis(ba);
}

void CalorimeterAnalysis::Set(TRecoMUV1Event* fMUV1Event, TRecoMUV2Event* fMUV2Event, TRecoMUV3Event* fMUV3Event, TRecoLKrCandidate* fLKrCandidate, Double_t fReft, Double_t fmmsq, TVector3 fPionPosition, TVector3 fPionMomentum, Bool_t fMCflag){
MUV1Event = fMUV1Event;
MUV2Event = fMUV2Event;
MUV3Event = fMUV3Event;
LKrCand = fLKrCandidate;
Reft = fReft;
mmsq = fmmsq;
PionPosition = fPionPosition;
PionMomentum = fPionMomentum; 
MCflag = fMCflag;
}

Int_t CalorimeterAnalysis::Rejection(){
Parameters* P = new Parameters(3809);
MUV3a->Set(MUV3Event, Reft, P->GettMUV3Low(), P->GettMUV3Hi(), MCflag);
MUV1a->Set(MUV1Event, Reft, P->GettMUV1Low(), P->GettMUV1Hi(), MCflag);
MUV2a->Set(MUV2Event, Reft, P->GettMUV2Low(), P->GettMUV2Hi(), MCflag);
TVector3 MUV3_pos= prop(PionPosition, PionMomentum*1000, P->GetMUV3StartPos());
TVector3 MUV1_pos= prop(PionPosition, PionMomentum*1000, P->GetMUV1StartPos());
TVector3 MUV2_pos= prop(PionPosition, PionMomentum*1000, P->GetMUV2StartPos());

Double_t MUV1E;
Double_t MUV2E; 
Double_t e = LKrCand->GetClusterEnergy(); 
Double_t ratio = (LKrCand->GetClusterSeedEnergy()/1000)/LKrCand->GetClusterEnergy();
if (ratio<0)ratio=0; 
Int_t Matched_MUV1= MUV1a->MatchedTrack(MUV1_pos);
if(Matched_MUV1)MUV1E=MUV1a->GetMatchedClusterE();
else  MUV1E = MUV1a->TotalEnergyHits();
Int_t Matched_MUV2 = MUV2a->MatchedTrack(MUV2_pos);
if(Matched_MUV2)MUV2E=MUV2a->GetMatchedClusterE();
else  MUV2E = MUV2a->TotalEnergyHits();
Bool_t pixel=false;
if (MUV1E==0&&MUV2E==0){
//FillHisto("MUV_matched", Matched_MUV2, Matched_MUV1);
//FillHisto("noMUV_LKrE", LKrCand->GetClusterEnergy()*1000); 
//FillHisto("noMUV_EoverP", LKrCand->GetClusterEnergy()/PionMomentum.Mag()); 
//FillHisto("noMUV_mmsq", mmsq);
//FillHisto("noMUV_CellEnergy", LKrCand->GetClusterEnergy(), LKrCand->GetNCells());
if(LKrCand->GetNCells()>(LKrCand->GetClusterEnergy()*1.66+10))pixel=true;
if(e>0){ 
//FillHisto("noMUV_seed", ratio, LKrCand->GetNCells()/LKrCand->GetClusterEnergy());
//FillHisto("noMUV_CellEnergy_number", LKrCand->GetNCells()/LKrCand->GetClusterEnergy());
}
}
if(e>0){
//FillHisto("seed", ratio, LKrCand->GetNCells()/LKrCand->GetClusterEnergy());
//FillHisto("CellEnergy_number", LKrCand->GetNCells()/LKrCand->GetClusterEnergy());
}
//FillHisto("CellEnergy", LKrCand->GetClusterEnergy(), LKrCand->GetNCells());
//Double_t MUV1E = 0;
//if (MUV1a->MatchedTrack(MUV1_pos) ==1) MUV1E = MUV1a->GetMatchedClusterE(); 
//Double_t MUV1E = MUV1a->GetMatchedClusterE();
fUserMethods->FillHisto("CalorimeterAnalsyis_Veto_mmsq_start", mmsq);
Double_t TotalE = 0;
TotalE = MUV1E + LKrCand->GetClusterEnergy()*1000 + MUV2E;
fUserMethods->FillHisto("LKrE", LKrCand->GetClusterEnergy()*1000);
fUserMethods->FillHisto("MUV1E", MUV1E);
fUserMethods->FillHisto("MUV2E", MUV2E);
fUserMethods->FillHisto("ECalo", TotalE);
if (TotalE>0)fUserMethods->FillHisto("MUV_energyratio_muon", MUV2E/TotalE, MUV1E/TotalE);
if (LKrCand->GetNCells()/LKrCand->GetClusterEnergy() < 1.66)return 0;
fUserMethods->FillHisto("CalorimeterAnalsyis_Veto_mmsq_cellenergyratio", mmsq);
//if(LKrCand->GetNCells()<(LKrCand->GetClusterEnergy()*1.66+10))return;
if (LKrCand->GetClusterEnergy()/PionMomentum.Mag()<0.07)return 0;
fUserMethods->FillHisto("CalorimeterAnalsyis_Veto_mmsq_LKrEPratio", mmsq);
Double_t RMS = sqrt(LKrCand->GetClusterRMSY()*LKrCand->GetClusterRMSY()+LKrCand->GetClusterRMSX()*LKrCand->GetClusterRMSX());
fUserMethods->FillHisto("CalorimeterAnalsyis_Veto_RMS", RMS);
//if (RMS < 8.25)return 0;
fUserMethods->FillHisto("CalorimeterAnalsyis_Veto_mmsq_RMS", mmsq);
if (LKrCand->GetClusterEnergy()<2)return 0;
fUserMethods->FillHisto("CalorimeterAnalsyis_Veto_mmsq_LKre", mmsq);

/*
if (index != -1)FillHisto("cluster_type", LKrCand->GetClusterEnergy(), LKrCand->GetNCells());
if (index!=-1){if (LKrCand->GetNCells()/LKrCand->GetClusterEnergy() < 1.66)return;}
FillHisto("mmsq_LKrRatio", mmsq); 
if (index!=-1){if (LKra->GetRMS() < 8.25)return;}
if (index!=-1){if (LKra->GetclusterE()/PionMomentum.Mag()<0.07)return;} 
*/
fUserMethods->FillHisto("mmsq_LKrRMS", mmsq);
if (TotalE>0)fUserMethods->FillHisto("MUV_energyratio", MUV2E/TotalE, MUV1E/TotalE);
//FillHisto("p1_Eratio",PionMomentum.Mag(),  MUV1E/TotalE);
//FillHisto("p2_Eratio",PionMomentum.Mag(),  MUV2E/TotalE);
//if (TotalE<4000)return;
if (TotalE<4000)return 0;
//FillHisto("ECalo_post", TotalE);
fUserMethods->FillHisto("mmsq_ECalo", mmsq);
Double_t EP = TotalE/(PionMomentum.Mag()*1000);
fUserMethods->FillHisto("EP", EP);
if (EP<0.1)return 0;
//FillHisto("EP_post", EP);
//FillHisto("mmsq_p_calo",PionMomentum.Mag(),mmsq);
fUserMethods->FillHisto("mmsq_EP", mmsq);
//FillHisto("mass_RICH", RICHa->Getmass_RICH());


if(Matched_MUV1&&Matched_MUV2&&MUV1a->GetDistance()<100&&MUV2a->GetDistance()<100)return 0;
fUserMethods->FillHisto("CalorimeterAnalsyis_Veto_mmsq_MUVmatched", mmsq);
/*
TRecoMUV1Candidate* MUV1Cand;
TRecoMUV2Candidate* MUV2Cand;
if(Matched_MUV1){
MUV1Cand = (TRecoMUV1Candidate*)MUV1Event->GetCandidate(MUV1a->GetMUV1ID());
if(MUV1Cand->GetShowerWidth()<1)return;
}
if(Matched_MUV2){
MUV2Cand = (TRecoMUV2Candidate*)MUV2Event->GetCandidate(MUV2a->GetMUV2ID());
if (MUV2Cand->GetShowerWidth()<1)return;
}
*/
Double_t E2 =  MUV2E/TotalE;
Double_t E1 = MUV1E/TotalE;
if((E1>0&&E1<(0.4*E2/0.1))||E1>(1-1*E2/0.3)/*||(E1==0&&E2>0.4)*/)return 0;
fUserMethods->FillHisto("CalorimeterAnalsyis_Veto_mmsq_MUVratio", mmsq);
//if(MUV1a->TotalEnergyHits()>0&&MUV2a->TotalEnergyHits()>0)return;

if(MUV3a->Veto()==0)return 0;
fUserMethods->FillHisto("mmsq_MUV3", mmsq);
Double_t Total1 = MUV1a->TotalEnergyHits()/1000;
Double_t Total2 = MUV2a->TotalEnergyHits()/1000;
Int_t n;
Double_t L = LKrCand->GetClusterEnergy();
if (L>Total1&&Total1>Total2)n=0;
else if (L>Total2&&Total2>Total1)n=1;
else if (Total1>L&&L>Total2)n=2;
else if (Total1>Total2&&Total2>L)n=3;
else if (Total2>Total1&&Total1>L)n=4;
else if (Total2>L&&L>Total1)n=5;
//if (n>3)return;

if ((PionMomentum.Mag()>15 && PionMomentum.Mag()<35)&&mmsq<0&&mmsq>-0.02){
//FillHisto("mmsq_num_nselection", mmsq);
//FillHisto("MUV_energyratio_diff", MUV2E/TotalE, MUV1E/TotalE);
}
fUserMethods->FillHisto("mmsq_Eratio", mmsq);
//FillHisto("mmsq_p_Eratio",PionMomentum.Mag(),mmsq);
//if (TotalE>0)FillHisto("MUV_energyratio_post", MUV2E/TotalE, MUV1E/TotalE);
//if ((PionMomentum.Mag()>15 && PionMomentum.Mag()<35)&&mmsq<0&&mmsq>-0.02)FillHisto("mmsq_num", mmsq); 
if(MUV1E==0&&MUV2E==0&&pixel==false)return 0; 
//if ((PionMomentum.Mag()>15 && PionMomentum.Mag()<35)&&mmsq<0&&mmsq>-0.02)FillHisto("mmsq_num_noMUV", mmsq);
if (TotalE>0)fUserMethods->FillHisto("MUV_energyratio_pion", MUV2E/TotalE, MUV1E/TotalE);
return 1;
}
void CalorimeterAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
MUV1a->SaveAllPlots();
MUV2a->SaveAllPlots();
MUV3a->SaveAllPlots();
}
void CalorimeterAnalysis::BookHistos(){
fUserMethods->BookHisto(new TH1I("CalorimeterAnalsyis_Phil", "passed", 5, 0, 5));
fUserMethods->BookHisto(new TH1I("CalorimeterAnalsyis_Veto_RMS", "LKr Cluster RMS", 120, 0, 12));
fUserMethods->BookHisto(new TH1I("LKrE", "LKr Cluster Energy", 100, 0, 50000));
fUserMethods->BookHisto(new TH1I("MUV1E", "MUV1 energy", 100, 0, 50000));
        fUserMethods->BookHisto(new TH1I("MUV2E", "MUV2 E", 100, 0, 50000));
 fUserMethods->BookHisto(new TH1I("ECalo", "Total Calorimeter Energy", 100, 0, 50000));
fUserMethods->BookHisto(new TH1F("mmsq_LKrRMS", "Did RMS do anything???", 250, -.35, .15));
 TH2I* MUV_energyratio = new TH2I("MUV_energyratio", "MUV2/ECalo(x) MUV1/ECalo (y)", 110, 0, 1.1, 110, 0, 1.1);
        MUV_energyratio->SetOption("COLZ");
        MUV_energyratio->Draw();
        fUserMethods->BookHisto(MUV_energyratio);
 TH2I* MUV_energyratio_muon = new TH2I("MUV_energyratio_muon", "MUV2/ECalo(x) MUV1/ECalo (y)", 110, 0, 1.1, 110, 0, 1.1);
        MUV_energyratio_muon->SetOption("COLZ");
        MUV_energyratio_muon->Draw();
        fUserMethods->BookHisto(MUV_energyratio_muon);
 TH2I* MUV_energyratio_pion = new TH2I("MUV_energyratio_pion", "MUV2/ECalo(x) MUV1/ECalo (y)", 110, 0, 1.1, 110, 0, 1.1);
        MUV_energyratio_pion->SetOption("COLZ");
        MUV_energyratio_pion->Draw();
        fUserMethods->BookHisto(MUV_energyratio_pion);
fUserMethods->BookHisto(new TH1F("CalorimeterAnalsyis_Veto_mmsq_start", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("CalorimeterAnalsyis_Veto_mmsq_cellenergyratio", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("CalorimeterAnalsyis_Veto_mmsq_LKrEPratio", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("CalorimeterAnalsyis_Veto_mmsq_RMS", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("CalorimeterAnalsyis_Veto_mmsq_LKre", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("mmsq_ECalo", "Missing Mass squared (Calo Energy)", 250, -.35, .15));
fUserMethods->BookHisto(new TH1I("EP", "Calorimeter Energy to Momentum ratio", 120, 0,1.2));
fUserMethods->BookHisto(new TH1F("mmsq_EP", "Missing Mass squared (E/P ratio)", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("CalorimeterAnalsyis_Veto_mmsq_MUVmatched", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("CalorimeterAnalsyis_Veto_mmsq_MUVratio", "Missing Mass squared", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("mmsq_MUV3", "Missing Mass squared (MUV3 Veto)", 250, -.35, .15));
fUserMethods->BookHisto(new TH1F("mmsq_Eratio", "Missing Mass squared", 250, -.35, .15));
MUV1a->BookHistos_Set();
MUV2a->BookHistos_Set();
MUV3a->BookHistos_Set();
MUV1a->BookHistos_MatchedTrack();
MUV2a->BookHistos_MatchedTrack();  
}


Int_t CalorimeterAnalysis::MuonSelection(){
Parameters* P = new Parameters(3809);
MUV3a->Set(MUV3Event, Reft, P->GettMUV3Low(), P->GettMUV3Hi(), MCflag);
MUV1a->Set(MUV1Event, Reft, P->GettMUV1Low(), P->GettMUV1Hi(), MCflag);
MUV2a->Set(MUV2Event, Reft, P->GettMUV2Low(), P->GettMUV2Hi(), MCflag);
TVector3 MUV3_pos= prop(PionPosition, PionMomentum*1000, P->GetMUV3StartPos());
TVector3 MUV1_pos= prop(PionPosition, PionMomentum*1000, P->GetMUV1StartPos());
TVector3 MUV2_pos= prop(PionPosition, PionMomentum*1000, P->GetMUV2StartPos());

Double_t MUV1E;
Double_t MUV2E; 
Double_t e = LKrCand->GetClusterEnergy(); 
Double_t ratio = (LKrCand->GetClusterSeedEnergy()/1000)/LKrCand->GetClusterEnergy();
if (ratio<0)ratio=0; 
Int_t Matched_MUV1= MUV1a->MatchedTrack(MUV1_pos);
if(Matched_MUV1)MUV1E=MUV1a->GetMatchedClusterE();
else  MUV1E = MUV1a->TotalEnergyHits();
Int_t Matched_MUV2 = MUV2a->MatchedTrack(MUV2_pos);
if(Matched_MUV2)MUV2E=MUV2a->GetMatchedClusterE();
else  MUV2E = MUV2a->TotalEnergyHits();

Double_t TotalE = 0;
TotalE = MUV1E + LKrCand->GetClusterEnergy()*1000 + MUV2E;
//fUserMethods->FillHisto("LKrE", LKrCand->GetClusterEnergy()*1000);
//fUserMethods->FillHisto("MUV1E", MUV1E);
//fUserMethods->FillHisto("MUV2E", MUV2E);
//fUserMethods->FillHisto("ECalo", TotalE);
if (LKrCand->GetNCells()/LKrCand->GetClusterEnergy() > 1.66)return 0;
//if(LKrCand->GetNCells()<(LKrCand->GetClusterEnergy()*1.66+10))return;
if (LKrCand->GetClusterEnergy()/PionMomentum.Mag()>0.07)return 0;
Double_t RMS = sqrt(LKrCand->GetClusterRMSY()*LKrCand->GetClusterRMSY()+LKrCand->GetClusterRMSX()*LKrCand->GetClusterRMSX());
if (RMS > 8.25)return 0;
if (LKrCand->GetClusterEnergy()>2)return 0;
//fUserMethods->FillHisto("mmsq_LKrRMS", mmsq);
if (TotalE>4000)return 0;
//fUserMethods->FillHisto("mmsq_ECalo", mmsq);
Double_t EP = TotalE/(PionMomentum.Mag()*1000);
//fUserMethods->FillHisto("EP", EP);
if (EP>0.1)return 0;
//fUserMethods->FillHisto("mmsq_EP", mmsq);
//if(Matched_MUV1&&Matched_MUV2&&MUV1a->GetDistance()<100&&MUV2a->GetDistance()<100)return 0;
//Double_t E2 =  MUV2E/TotalE;
//Double_t E1 = MUV1E/TotalE;
//if((E1>0&&E1<(0.4*E2/0.1))||E1>(1-1*E2/0.3)/*||(E1==0&&E2>0.4)*/)return 0;
//if(MUV1a->TotalEnergyHits()>0&&MUV2a->TotalEnergyHits()>0)return;
if(MUV3a->Veto()==1)return 0;
return 1;
}

Int_t CalorimeterAnalysis::MUV3Efficiency(){
Parameters* P = new Parameters(3809);
MUV3a->Set(MUV3Event, Reft, P->GettMUV3Low(), P->GettMUV3Hi(), MCflag);
MUV1a->Set(MUV1Event, Reft, P->GettMUV1Low(), P->GettMUV1Hi(), MCflag);
MUV2a->Set(MUV2Event, Reft, P->GettMUV2Low(), P->GettMUV2Hi(), MCflag);
TVector3 MUV3_pos= prop(PionPosition, PionMomentum*1000, P->GetMUV3StartPos());
TVector3 MUV1_pos= prop(PionPosition, PionMomentum*1000, P->GetMUV1StartPos());
TVector3 MUV2_pos= prop(PionPosition, PionMomentum*1000, P->GetMUV2StartPos());
Double_t MUV1E;
Double_t MUV2E; 
Double_t e = LKrCand->GetClusterEnergy(); 
Double_t ratio = (LKrCand->GetClusterSeedEnergy()/1000)/LKrCand->GetClusterEnergy();
if (ratio<0)ratio=0; 
Int_t Matched_MUV1= MUV1a->MatchedTrack(MUV1_pos);
if(Matched_MUV1)MUV1E=MUV1a->GetMatchedClusterE();
else  MUV1E = MUV1a->TotalEnergyHits();
Int_t Matched_MUV2 = MUV2a->MatchedTrack(MUV2_pos);
if(Matched_MUV2)MUV2E=MUV2a->GetMatchedClusterE();
else  MUV2E = MUV2a->TotalEnergyHits();
Double_t TotalE = 0;
TotalE = MUV1E + LKrCand->GetClusterEnergy()*1000 + MUV2E;
//if (LKrCand->GetNCells()/LKrCand->GetClusterEnergy() > 1.66)return 0;
//if (LKrCand->GetClusterEnergy()/PionMomentum.Mag()>0.07)return 0;
Double_t RMS = sqrt(LKrCand->GetClusterRMSY()*LKrCand->GetClusterRMSY()+LKrCand->GetClusterRMSX()*LKrCand->GetClusterRMSX());
//if (RMS > 8.25)return 0;
if (LKrCand->GetClusterEnergy()>2)return 0;
if (TotalE>4000)return 0;
Double_t EP = TotalE/(PionMomentum.Mag()*1000);
//if (EP>0.1)return 0;
return 1;
}



void CalorimeterAnalysis::Phil(){
Parameters* P = new Parameters(3809);
MUV3a->Set(MUV3Event, Reft, P->GettMUV3Low(), P->GettMUV3Hi(), MCflag);
MUV1a->Set(MUV1Event, Reft, P->GettMUV1Low(), P->GettMUV1Hi(), MCflag);
MUV2a->Set(MUV2Event, Reft, P->GettMUV2Low(), P->GettMUV2Hi(), MCflag);
TVector3 MUV3_pos= prop(PionPosition, PionMomentum*1000, P->GetMUV3StartPos());
TVector3 MUV1_pos= prop(PionPosition, PionMomentum*1000, P->GetMUV1StartPos());
TVector3 MUV2_pos= prop(PionPosition, PionMomentum*1000, P->GetMUV2StartPos());
Double_t MUV1E;
Double_t MUV2E; 
Double_t e = LKrCand->GetClusterEnergy(); 
Double_t ratio = (LKrCand->GetClusterSeedEnergy()/1000)/LKrCand->GetClusterEnergy();
if (ratio<0)ratio=0; 
Int_t Matched_MUV1= MUV1a->MatchedTrack(MUV1_pos);
if(Matched_MUV1)MUV1E=MUV1a->GetMatchedClusterE();
else  MUV1E = MUV1a->TotalEnergyHits();
Int_t Matched_MUV2 = MUV2a->MatchedTrack(MUV2_pos);
if(Matched_MUV2)MUV2E=MUV2a->GetMatchedClusterE();
else  MUV2E = MUV2a->TotalEnergyHits();
Double_t TotalE = 0;
TotalE = MUV1E + LKrCand->GetClusterEnergy()*1000 + MUV2E;
fUserMethods->FillHisto("LKrE", LKrCand->GetClusterEnergy()*1000);
fUserMethods->FillHisto("MUV1E", MUV1E);
fUserMethods->FillHisto("MUV2E", MUV2E);
fUserMethods->FillHisto("ECalo", TotalE);
if (LKrCand->GetClusterEnergy()<2)return;
//if (TotalE>0)fUserMethods->FillHisto("MUV_energyratio", MUV2E/TotalE, MUV1E/TotalE);
fUserMethods->FillHisto("CalorimeterAnalsyis_Phil", 0);
if (TotalE>0)fUserMethods->FillHisto("MUV_energyratio", MUV2E/TotalE, MUV1E/TotalE);
Double_t E2 =  MUV2E/TotalE;
Double_t E1 = MUV1E/TotalE;
if((E1>0&&E1<(0.4*E2/0.1)))return;
fUserMethods->FillHisto("CalorimeterAnalsyis_Phil", 1);
//if((E1>0&&E1<(0.4*E2/0.1))/*||E1>(1-1*E2/0.3))return 0;
//fUserMethods->FillHisto("CalorimeterAnalsyis_Veto_mmsq_MUVratio", mmsq);
//if ((PionMomentum.Mag()>15 && PionMomentum.Mag()<35)&&mmsq<0&&mmsq>-0.02)
return;
}

