#include "BlueTubeTrackerAnalysis.hh"
#include "AnalysisTools.hh"
#include "BlueTubeTracker.hh"
#include "Event.hh"

BlueTubeTrackerAnalysis::BlueTubeTrackerAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}

void BlueTubeTrackerAnalysis::Set(BlueTubeTracker* ftracker, TVector3 fK_pos, TVector3 fDecayVertex, TVector3 fPionPosition, TVector3 fPionMomentum){
Tracker = ftracker; 
K_pos = fK_pos;
DecayVertex = fDecayVertex;
PionPosition = fPionPosition;
PionMomentum = fPionMomentum; 
 Double_t zvertex = DecayVertex.Z();
 Double_t zmax = 183311.;
  Double_t zmin = 101800.;
  Double_t gtom = 1000.;
  TVector3 initPos = PionPosition;
  if (initPos.Z()>zmax) {
    initPos.SetX(PionPosition.X()+(PionMomentum.X()/PionMomentum.Z())*(zmax-PionPosition.Z()));
    initPos.SetY(PionPosition.Y()+(PionMomentum.Y()/PionMomentum.Z())*(zmax-PionPosition.Z()));
    initPos.SetZ(zmax);
  }
  Tracker->SetCharge(1);
  Tracker->SetInitialPosition(initPos); // Unit: mm
  Tracker->SetInitialMomentum(PionMomentum.X()*gtom,PionMomentum.Y()*gtom,PionMomentum.Z()*gtom); // Unit: MeV
  Double_t zEndTracker = zvertex;
  if (zvertex<zmin) zEndTracker = zmin;
  if (zvertex>=zmax) zEndTracker = zmax;
  Tracker->SetZFinal(zEndTracker); // mm
  Tracker->TrackParticle(); // 
   vertexCorr = Tracker->GetFinalPosition(); // New vertex
   TVector3 correctedMomentum = Tracker->GetFinalMomentum();
pionMomentumCorr.SetXYZM(correctedMomentum.X()/gtom,correctedMomentum.Y()/gtom,correctedMomentum.Z()/gtom,0.13957018); // New pion momentum
  if (zvertex<zmin) vertexCorr = vertexCorr+correctedMomentum*(1./correctedMomentum.Z())*(zvertex-vertexCorr.Z());
//TVector2 cda_new_vect(K_pos.X()-vertexCorr.X(), K_pos.Y()-vertexCorr.Y()); //pretty sure this is wrong
}
//void LAVAnalysis::SaveAllPlots(){
//fUserMethods->SaveAllPlots();
//}

//void LAVAnalysis::BookHistos_Set(){
//fUserMethods->BookHisto(new TH1F("LAV_Set_reft","reference time", 200, -50, 50));
//}

