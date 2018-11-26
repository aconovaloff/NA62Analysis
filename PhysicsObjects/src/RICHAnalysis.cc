#include "RICHAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoRICHEvent.hh"
#include "Event.hh"
#include "UserMethods.hh"
using namespace TMath;

RICHAnalysis::RICHAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}
 
void RICHAnalysis::Set(TRecoRICHEvent* fevent, double freftime, double ftLow, double ftHi, bool fMCflag){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
MCflag = fMCflag;
mass_RICH = -999;
pick_RICH_mmsq= -999;
TRecoRICHCandidate* RICHCand;
for(int i=0; i<event->GetNCandidates(); i++){
  RICHCand = (TRecoRICHCandidate*)event->GetCandidate(i);
  Double_t Tdiff = reftime - RICHCand->GetTime();
  fUserMethods->FillHisto("RICH_Set_reft", Tdiff);
}
}

Int_t RICHAnalysis::PionID(TVector3 Pos, TVector3 Trackp){
bool pion = false;
TRecoRICHCandidate *RICHCand;
if (event->GetNRingCandidates() == 1){
		TVector3 Pi3 = Trackp; 
                TVector3 Pion_Prop_RICH_test = Pos;
                Double_t xcenter_2 = 1633-(Pion_Prop_RICH_test.X()+17000*tan(atan((1633-Pion_Prop_RICH_test.X())/17000)-(Pi3.Px()-0.270)/Pi3.Pz()));
                Double_t ycenter_2 = -(Pion_Prop_RICH_test.Y()+17000*tan(atan((-Pion_Prop_RICH_test.Y())/17000)-(Pi3.Py())/Pi3.Pz()));
                for(int i=0; i<event->GetNRingCandidates(); i++)    {
                RICHCand = (TRecoRICHCandidate*)event->GetRingCandidate(i);
                TVector2 dist_2(RICHCand->GetRingCenter().X() - xcenter_2, RICHCand->GetRingCenter().Y() - ycenter_2);
                if(dist_2.Mod() < 10 || dist_2.Mod() > 90) continue;
                if(dist_2.Mod() > 40 && dist_2.Mod()<75) continue;
                if (RICHCand->GetRingRadius()>180||(RICHCand->GetRingRadius()>((90/25)*Pi3.Mag()+90)))continue;
		fUserMethods->FillHisto("RICH_PionID_pi3_rad", Pi3.Mag(), RICHCand->GetRingRadius());
                pion = true;

}
}
if (!pion) return 0;
if (pion) return 1; 
}

Int_t RICHAnalysis::MuonID(TVector3 Pos, TVector3 Trackp){
bool pion = false;
TRecoRICHCandidate *RICHCand;
if (event->GetNRingCandidates() == 1){
		TVector3 Pi3 = Trackp; 
                TVector3 Pion_Prop_RICH_test = Pos;
                Double_t xcenter_2 = 1633-(Pion_Prop_RICH_test.X()+17000*tan(atan((1633-Pion_Prop_RICH_test.X())/17000)-(Pi3.Px()-0.270)/Pi3.Pz()));
                Double_t ycenter_2 = -(Pion_Prop_RICH_test.Y()+17000*tan(atan((-Pion_Prop_RICH_test.Y())/17000)-(Pi3.Py())/Pi3.Pz()));
                for(int i=0; i<event->GetNRingCandidates(); i++)    {
                RICHCand = (TRecoRICHCandidate*)event->GetRingCandidate(i);
                TVector2 dist_2(RICHCand->GetRingCenter().X() - xcenter_2, RICHCand->GetRingCenter().Y() - ycenter_2);
                if(dist_2.Mod() < 10 || dist_2.Mod() > 90) continue;
                if(dist_2.Mod() > 40 && dist_2.Mod()<75) continue;
                if (RICHCand->GetRingRadius()>180||(RICHCand->GetRingRadius()>((90/25)*Pi3.Mag()+90)))continue;
                pion = true;

}
}
if (pion) return 0;
if (!pion) return 1; 
}

Int_t RICHAnalysis::MatchedTrack(TVector3 Pos, TLorentzVector Pion4, TLorentzVector Kaon4){
        	double trim = 0.0012;
		TLorentzVector Pi4 = Pion4; 		
		TVector3 Pi3 = Pi4.Vect(); 
		TLorentzVector K4 = Kaon4;
		TRecoRICHCandidate *RICHCand;
//                double mass_RICH;
                double mmsq_RICH;
		mass_RICH = 9999; 
		double mK = 0.493677; //in GeV
        	double mPi = 0.13957018; //in GeV
        	double f = 17020; //focal length of rich
        	double rMax = 190.2; //maximum ring radius in RICH in mm
                TVector3 Pion_Prop_RICH = Pos;
//		TVector3 Pi3 = Trackp;
                int s = 0;
                Double_t xcenter = 1633-(Pion_Prop_RICH.X()+17020*tan(atan((1633-Pion_Prop_RICH.X())/17020)-(Pi3.Px()-0.270)/Pi3.Pz()));
                Double_t ycenter = -(Pion_Prop_RICH.Y()+17020*tan(atan((-Pion_Prop_RICH.Y())/17020)-(Pi3.Py())/Pi3.Pz()));
                double min_dist_RICH = 999999;
                int RICH_minD;
                for(int i=0; i<event->GetNRingCandidates(); i++)    {
                RICHCand = (TRecoRICHCandidate*)event->GetRingCandidate(i);
		if (event->GetNRingCandidates() != 1)continue;
		double Tdiff_RICH = reftime-RICHCand->GetTime();
                if(!MCflag){if(Tdiff_RICH < tLow || Tdiff_RICH > tHi) continue;}
                TVector2 dist(RICHCand->GetRingCenter().X() - xcenter, RICHCand->GetRingCenter().Y() - ycenter);
                fUserMethods->FillHisto("RICH_MatchedTrack_D", dist.Mod());
//                if(dist.Mod() < 10 || dist.Mod() > 110) continue;
		if(dist.Mod() > 85) continue;
//                if(dist.Mod() > 40  && dist.Mod()<90) continue; very confused, check different runs and number of bursts
		if (dist.Mod()>40&&dist.Mod()<66)continue;
                fUserMethods->FillHisto("RICH_MatchedTrack_D_post", dist.Mod());
		if(RICHCand->GetRingRadius()>rMax)continue;
                double Pi3_RICH = (mPi*f)/sqrt(rMax*rMax - RICHCand->GetRingRadius()*RICHCand->GetRingRadius());
//                mass_RICH = (Pi3.Mag()/f)*sqrt(rMax*rMax - RICHCand->GetRingRadius()*RICHCand->GetRingRadius());
  Double_t changle = atan(RICHCand->GetRingRadius()/f);
  mass_RICH = 1.000062*1.000062*cos(changle)*cos(changle)-1>0 ? Pi3.Mag()*sqrt(1.000062*1.000062*cos(changle)*cos(changle)-1) : 9999.;

                s++;
                double m_pion =  Pi4.Mag2();
                double E_RICH = sqrt(Pi3_RICH*Pi3_RICH + m_pion*m_pion);
                double rad = RICHCand->GetRingRadius();
                TLorentzVector Pi4_RICH(Pi3.Unit().X()*Pi3_RICH, Pi3.Unit().Y()*Pi3_RICH, Pi3.Unit().Z()*Pi3_RICH, E_RICH);
                TLorentzVector mmsqRICHV = K4 -Pi4_RICH;
                mmsq_RICH = mmsqRICHV.Mag2();
		fUserMethods->FillHisto("RICH_MatchedTrack_pi3_rad", Pi3.Mag(), RICHCand->GetRingRadius());
//                FillHisto("RICH_MatchedTrack_mmsq_RICHmmsq", MassMissSq, mmsq_RICH);
//                        if (dist.Mod() < min_dist_RICH) {
  //                              min_dist_RICH = dist.Mod();
//                                pick_RICH_mmsq = mmsq_RICH;
//                                RICH_minD = i;
  //                                                      }

                                                                        }
//                if(s == 0)return 0; 
if(s==0/*||mmsq_RICH <-0.03*/) return 0; //don't think i should use this anymore
//		if(s>0)return 1; 

//                                RICHCand = (TRecoRICHCandidate*)event->GetRingCandidate(RICH_minD);
   //                             double Tdiff_RICH = reftime-RICHCand->GetTime();
 //                               if(Tdiff_RICH < tLow || Tdiff_RICH > tHi) return 0;
//if (pick_RICH_mmsq <-0.03) return 0;
		else return 1;

}

void RICHAnalysis::BookHistos_PionID(){
        TH2I* RICH_PionID_pi3_rad = new TH2I("RICH_PionID_pi3_rad", "Pion Momentum vs RICH Ring Radius (with 15<p<35", 100, 0, 100, 100, 0, 300);
        RICH_PionID_pi3_rad->SetOption("COLZ");
        RICH_PionID_pi3_rad->Draw();
        fUserMethods->BookHisto(RICH_PionID_pi3_rad);
}

void RICHAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
}
void RICHAnalysis::BookHistos_MatchedTrack(){
        fUserMethods->BookHisto(new TH1I("RICH_MatchedTrack_D", "Difference b/w ring center and expected position", 1000, 0, 250));
        fUserMethods->BookHisto(new TH1I("RICH_MatchedTrack_D_post", "Difference b/w ring center and expected position", 1000, 0, 250));
        TH2I* RICH_MatchedTrack_pi3_rad = new TH2I("RICH_MatchedTrack_pi3_rad", "Pion Momentum vs RICH Ring Radius (Not Matched)", 100, 0, 100, 100, 0, 300);
        RICH_MatchedTrack_pi3_rad->SetOption("COLZ");
        RICH_MatchedTrack_pi3_rad->Draw();
        fUserMethods->BookHisto(RICH_MatchedTrack_pi3_rad);
}
void RICHAnalysis::BookHistos_Set(){
        fUserMethods->BookHisto(new TH1F("RICH_Set_reft", "Reference Time", 200, -25, 25));
}


