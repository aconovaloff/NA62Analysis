#include "StrawAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoSpectrometerEvent.hh"
#include "Event.hh"

StrawAnalysis::StrawAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}

void StrawAnalysis::Set(TRecoSpectrometerEvent* fevent, double freftime, double ftLow, double ftHi, bool fMCflag){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
MCflag = fMCflag;
Pi3.SetXYZ(0,0,0);
PiPos.SetXYZ(0,0,0);
Pi4.SetPxPyPzE(0,0,0,0);
}

Int_t StrawAnalysis::SingleTrack(){
TRecoSpectrometerCandidate *StrawCand;
        Double_t Straw1 = 183508.0;
        Double_t Straw2 = 194063.5;
        Double_t Straw3 = 204459.0;
        Double_t Straw4 = 218885.0;
Int_t track = 0;
for(int j=0; j<event->GetNCandidates(); j++)  {
if (event->GetNCandidates() != 1)continue;
StrawCand = (TRecoSpectrometerCandidate*)event->GetCandidate(j);
PiPos = ComputePosition(StrawCand, Straw1);
Pi4 = ComputeMomentum(StrawCand,0.13957018);
Pi3 = Pi4.Vect();
Pi4Prop = ComputeMomentumAfterMagnet(StrawCand,0.13957018);
Pi3Prop = Pi4Prop.Vect();
PiPosProp = ComputePosition(StrawCand, Straw4);

Double_t quality_momentum = fabs(StrawCand->GetMomentumBeforeFit())-StrawCand->GetMomentum();
fUserMethods->FillHisto("StrawAnalysis_SingleTrack_QualityMomentum", quality_momentum);
if (StrawCand->GetChi2()>20) continue;
if (fabs(quality_momentum)>20000) continue;
if (StrawCand->GetNChambers()<4) continue;
TVector3 At_Straw1 = ComputePosition(StrawCand, Straw1);
TVector3 At_Straw2 = ComputePosition(StrawCand, Straw2);
TVector3 At_Straw3 = ComputePosition(StrawCand, Straw3);
TVector3 At_Straw4 = ComputePosition(StrawCand, Straw4);

TVector2 At_Straw1_2D (At_Straw1.X(), At_Straw1.Y());
TVector2 At_Straw2_2D (At_Straw2.X(), At_Straw2.Y());
TVector2 At_Straw3_2D (At_Straw3.X(), At_Straw3.Y());
TVector2 At_Straw4_2D (At_Straw4.X(), At_Straw4.Y());
if (At_Straw1_2D.Mod() < 75 || At_Straw1_2D.Mod() > 1000) continue;
if (At_Straw2_2D.Mod() < 75 || At_Straw2_2D.Mod() > 1000) continue;
if (At_Straw3_2D.Mod() < 75 || At_Straw3_2D.Mod() > 1000) continue;
if (At_Straw4_2D.Mod() < 75 || At_Straw4_2D.Mod() > 1000) continue;
track++;
ftime = StrawCand->GetTime(); 
}
if (track==0) return 0;
if (track>0) return 1;

}

TVector3 StrawAnalysis::GetPionMomentum(){
return Pi3;
}

TVector3 StrawAnalysis::GetPionPosition(){
return PiPos; 
}

TVector3 StrawAnalysis::GetPionMomentumProp(){
return Pi3Prop;
}

TVector3 StrawAnalysis::GetPionPositionProp(){
return PiPosProp;
}

TLorentzVector StrawAnalysis::GetPionFourMomentum(){
return Pi4;
}

TLorentzVector StrawAnalysis::ComputeMomentum(TRecoSpectrometerCandidate *cand, Double_t mass) {
    TLorentzVector pmom;
    Double_t thetaX = cand->GetSlopeXBeforeMagnet();
    Double_t thetaY = cand->GetSlopeYBeforeMagnet();
    Double_t pmag = cand->GetMomentum();
    Double_t pmomz = pmag/sqrt(1.+thetaX*thetaX+thetaY*thetaY);
    Double_t pmomx = pmomz*thetaX;
    Double_t pmomy = pmomz*thetaY;
    pmom.SetXYZM(pmomx/1000,pmomy/1000,pmomz/1000,mass);
    return pmom;
}

TLorentzVector StrawAnalysis::ComputeMomentumAfterMagnet(TRecoSpectrometerCandidate *cand, Double_t mass) {
    TLorentzVector pmom;
    Double_t thetaX = cand->GetSlopeXAfterMagnet();
    Double_t thetaY = cand->GetSlopeYAfterMagnet();
    Double_t pmag = cand->GetMomentum();
    Double_t pmomz = pmag/sqrt(1.+thetaX*thetaX+thetaY*thetaY);
    Double_t pmomx = pmomz*thetaX;
    Double_t pmomy = pmomz*thetaY;
    pmom.SetXYZM(pmomx/1000,pmomy/1000,pmomz/1000,mass);
    return pmom;
}


TVector3 StrawAnalysis::ComputePosition(TRecoSpectrometerCandidate *cand, Double_t zpos) {
  Double_t posx;
  Double_t posy;
if (zpos<196345) {
 posx = cand->GetPositionBeforeMagnet().X()+cand->GetSlopeXBeforeMagnet()*(zpos-cand->GetPositionBeforeMagnet().Z());
    posy = cand->GetPositionBeforeMagnet().Y()+cand->GetSlopeYBeforeMagnet()*(zpos-cand->GetPositionBeforeMagnet().Z());
  } else {
    posx = cand->GetPositionAfterMagnet().X()+cand->GetSlopeXAfterMagnet()*(zpos-cand->GetPositionAfterMagnet().Z());
    posy = cand->GetPositionAfterMagnet().Y()+cand->GetSlopeYAfterMagnet()*(zpos-cand->GetPositionAfterMagnet().Z());
  }
  return TVector3(posx,posy,zpos);
}

void StrawAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
}

Int_t StrawAnalysis::MatchedTime(Double_t reft){
Double_t ReferenceTime = reft;
TRecoSpectrometerCandidate *StrawCand;
Int_t Straw=0;
for(int i=0; i<event->GetNCandidates(); i++){
  StrawCand = (TRecoSpectrometerCandidate*)event->GetCandidate(i);
  Double_t Tdiff = ReferenceTime - StrawCand->GetTime();
  fUserMethods->FillHisto("Straw_MatchedTime_reft", Tdiff);
  Double_t Tdiff2 = ReferenceTime - ftime; //ftime will be time of "best track" when using multiple tracks (should I ever get there) 
  if(!MCflag){if(Tdiff2 < tLow || Tdiff2 > tHi)continue;}
  Straw++;
}
if (Straw==0)return 0;
if (Straw>0)return 1; 
}

void StrawAnalysis::BookHistos_MatchedTime(){
fUserMethods->BookHisto(new TH1F("Straw_MatchedTime_reft", "Reference Time", 100, -50, 50));
}

void StrawAnalysis::BookHistos_SingleTrack(){
fUserMethods->BookHisto(new TH1I("StrawAnalysis_SingleTrack_QualityMomentum", "Quality Momentum",1200,-60000, 60000));
}
