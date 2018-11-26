#include "LAVAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoLAVEvent.hh"
#include "Event.hh"

LAVAnalysis::LAVAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}

void LAVAnalysis::Set(TRecoLAVEvent* fevent, double freftime, double ftLow, double ftHi, bool fMCflag){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
MCflag = fMCflag;
TRecoLAVCandidate* LAVCand;
for(int i=0; i<event->GetNCandidates(); i++){
  LAVCand = (TRecoLAVCandidate*)event->GetCandidate(i);
  Double_t Tdiff = reftime - LAVCand->GetTime();
  fUserMethods->FillHisto("LAV_Set_reft", Tdiff);
}
}

Int_t LAVAnalysis::Veto(){
TRecoLAVCandidate* LAVCand;
				double Emin =150;
				Int_t nLAVRecoCandi = event->GetNCandidates();
                                TClonesArray& clusArray = (* (event->GetCandidates()));
                                int LAV_N = 0;
                                int MIP = 0;
                                for (Int_t iCl = 0; iCl < nLAVRecoCandi; iCl++) {
//                                TRecoLAVCandidate* LAVCand = (TRecoLAVCandidate*) clusArray[iCl];
				LAVCand = (TRecoLAVCandidate*)event->GetCandidate(iCl);
                                double Tdiff = reftime - LAVCand->GetTime();
				LAVCand->GetClusterType(); 
                                if(!MCflag){if(Tdiff<tLow||Tdiff>tHi)continue;}
                                if(LAVCand->GetEnergy()>Emin)LAV_N++;
                                if(LAVCand->GetEnergy()<Emin) MIP++;
                                                                                }

if (LAV_N>0)return 0;
if (LAV_N==0)return 1;

}

void LAVAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
}

void LAVAnalysis::BookHistos_Set(){
fUserMethods->BookHisto(new TH1F("LAV_Set_reft","reference time", 800, -50, 50));
}
/*
Int_t LAVAnalysis::SensitiveRegion(TVector3 PhotonMomentum, TVector3 DecayVertex){
Double_t LAV_start[12];
Double_t LAV_end[12];
Double_t LAV_in[12];
Double_t LAV_out[12];
LAV_start[0] = 
LAV_end[0] = 
LAV_in[0] = 
LAV_out[0] =
 LAV_start[1] = 
LAV_end[1] = 
LAV_in[1] = 
LAV_out[1] = 
LAV_start[2] = 
LAV_end[2] = 
LAV_in[2] = 
LAV_out[2] = 
LAV_start[3] = 
LAV_end[3] = 
LAV_in[3] = 
LAV_out[3] =
 LAV_start[4] = 
LAV_end[4] = 
LAV_in[4] = 
LAV_out[4] = 
LAV_start[5] = 
LAV_end[5] = 
LAV_in[5] = 
LAV_out[5] =
 LAV_start[0] = 
LAV_end[0] = 
LAV_in[0] = 
LAV_out[0] =
 LAV_start[0] = 
LAV_end[0] = 
LAV_in[0] = 
LAV_out[0] =
 LAV_start[0] = 
LAV_end[0] = 
LAV_in[0] = 
LAV_out[0] =
 LAV_start[0] = 
LAV_end[0] = 
LAV_in[0] = 
LAV_out[0] =
 LAV_start[0] = 
LAV_end[0] = 
LAV_in[0] = 
LAV_out[0] = 
//LAV_start[1]=
//Double_t LAV1_start = 100;
//Double_t LAV1_end = 200;
//Double_t LAV1_in = 5;
//Double_t LAV1_out = 10; 
Int_t station = 0:
for (int i = 0; i=11; i++){
TVector3 LAV_ent = Propagate(0, &PhotonMomentum, &DecayVertex, LAV_start[i]);
TVector3 LAV_ex = Propagate(0, &PhotonMomentum, &LAV_ent, LAV_end[i]);
Double_t TwoD_front = sqrt(LAV_ent.X()*LAV_ent.X() + LAV_ent.Y()*LAV_ent.Y());
Double_t TwoD_end = sqrt(LAV_ex.X()*LAV_ex.X() + LAV_ex.Y()*LAV_ex.Y();
if (TwoD_front>LAV_in[i]&&TwoD_front<LAV_out[i]){
station = i+1;
break;
  }
if (TwoD_end>LAV_in[i]){
station = i+1;
break;
  }
 }
return station; 
//TVector3 K_pos= Propagate(1,&KaonMomentum,&KaonPosition, Vert.Z());
}
*/
