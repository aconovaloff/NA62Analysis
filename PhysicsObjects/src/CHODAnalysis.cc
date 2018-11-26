#include "CHODAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoCHODEvent.hh"
#include "Event.hh"

CHODAnalysis::CHODAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}

void CHODAnalysis::Set(TRecoCHODEvent* fevent, double freftime, double ftLow, double ftHi, bool fMCflag ){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
minDCandCHOD = -1; 
TRecoCHODCandidate *CHODCand;
MCflag = fMCflag; 
 for(int i=0; i<event->GetNCandidates(); i++){
  CHODCand = (TRecoCHODCandidate*)event->GetCandidate(i);
  Double_t Tdiff = CHODCand->GetTime() -reftime; 
  fUserMethods->FillHisto("CHOD_Set_reft", Tdiff);
} 
}

Int_t CHODAnalysis::MatchedTrack(TVector3 Pos){
		TRecoCHODCandidate *CHODCand;
                TVector3 Pion_Prop_CHOD = Pos;
                int p = 0;
                TVector2 Pion_Prop_CHOD2D(Pion_Prop_CHOD.X(), Pion_Prop_CHOD.Y());//check units
                double distanceCHOD = 99999;
                for(int i=0; i<event->GetNCandidates(); i++)        {
               CHODCand = (TRecoCHODCandidate*)event->GetCandidate(i);
                TVector2 HitPos(CHODCand->GetHitPosition().X() , CHODCand->GetHitPosition().Y());
                TVector2 distCHOD = Pion_Prop_CHOD2D - HitPos;
                if (distCHOD.Mod() < distanceCHOD)      {
                        distanceCHOD = distCHOD.Mod();
                        minDCandCHOD = i;
                        if (distanceCHOD<80) p++;
                                                        }
                                                                        }
		if(distanceCHOD<99999)fUserMethods->FillHisto("CHOD_MatchedTrack_MinD", distanceCHOD); 
                        fUserMethods->FillHisto("CHOD_MatchedTrack_MinD_mirror", distanceCHOD);
                        fUserMethods->FillHisto("CHOD_MatchedTrack_MinD_mirror",-distanceCHOD);
                if(p==0)return 0;
		if(p>0)return 1;
}

Int_t CHODAnalysis::MatchedTrack_t(TVector3 Pos){
                TRecoCHODCandidate *CHODCand;
                TVector3 Pion_Prop_CHOD = Pos;
                int p = 0;
                TVector2 Pion_Prop_CHOD2D(Pion_Prop_CHOD.X(), Pion_Prop_CHOD.Y());//check units
                double distanceCHOD = 99999;
                for(int i=0; i<event->GetNCandidates(); i++)        {		
                CHODCand = (TRecoCHODCandidate*)event->GetCandidate(i);
		Double_t Tdiff = reftime - CHODCand->GetTime();
		if (!MCflag){if(Tdiff<tLow||Tdiff>tHi)continue;} 
                TVector2 HitPos(CHODCand->GetHitPosition().X() , CHODCand->GetHitPosition().Y());
                TVector2 distCHOD = Pion_Prop_CHOD2D - HitPos;
                if (distCHOD.Mod() < distanceCHOD)      {
                        distanceCHOD = distCHOD.Mod();
                        minDCandCHOD = i;
                        if (distanceCHOD<80) p++;
                                                        }
                                                                        }
                if(distanceCHOD<99999){
			fUserMethods->FillHisto("CHOD_MatchedTrack_MinD", distanceCHOD); 
//			fUserMethods->FillHisto("CHOD_MatchedTrack_MinD_mirror", distanceCHOD);
//			fUserMethods->FillHisto("CHOD_MatchedTrack_MinD_mirror",-distanceCHOD);
                if(p==0)return 0;
                if(p>0)return 1;
}
}
Double_t CHODAnalysis::GetTime(){
TRecoCHODCandidate *CHODCand;
CHODCand = (TRecoCHODCandidate*)event->GetCandidate(minDCandCHOD);
Double_t time = CHODCand->GetTime();
return time; 
}

void CHODAnalysis::BookHistos_MatchedTrack(){
fUserMethods->BookHisto(new TH1F("CHOD_MatchedTrack_MinD", "Minimum Distance", 500, 0, 500)); 
fUserMethods->BookHisto(new TH1F("CHOD_MatchedTrack_MinD_mirror", "Minimum Distance", 1000, -500, 500));
}
void CHODAnalysis::BookHistos_Set(){
fUserMethods->BookHisto(new TH1F("CHOD_Set_reft", "CHOD Reference Time", 100, -50, 50));
}

void CHODAnalysis::BookHistos_MinDisc(){
fUserMethods->BookHisto(new TH1F("CHOD_MinDisc_disc", "CHOD Discriminate", 100, 0, 10));
}

void CHODAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
}

Int_t CHODAnalysis::MinDisc(TVector3 Pos, Double_t TrackTime){
                TRecoCHODCandidate *CHODCand;
                TVector3 Pion_Prop_CHOD = Pos;
                Int_t d = 0;
                TVector2 Pion_Prop_CHOD2D(Pion_Prop_CHOD.X(), Pion_Prop_CHOD.Y());//check units
//                double distanceCHOD = 99999;
                Double_t mindisc = 999999; ///////////check 1.1
                Double_t sigmaT = 5.1;
                Double_t sigmaTHdiff= 5.1;
                Double_t sigmasep = 80;
		MinDiscID = -1; 
                for(int i=0; i<event->GetNCandidates(); i++)        {
                CHODCand = (TRecoCHODCandidate*)event->GetCandidate(i);
                Double_t Tdiff = CHODCand->GetTime() - TrackTime;
                 if (!MCflag){if(Tdiff<tLow||Tdiff>tHi)continue;}
                TVector2 HitPos(CHODCand->GetHitPosition().X() , CHODCand->GetHitPosition().Y());
                TVector2 distCHOD = Pion_Prop_CHOD2D - HitPos;
                if (distCHOD.Mod()>80)continue;
                Double_t sep = distCHOD.Mod();
//              Double_t disc = (Tdiff*Tdiff)/(3*sigmaT*sigmaT)+(THdiff*THdiff)/(3*sigmaTHdiff*sigmaTHdiff)+(sep*sep)/(3*sigmasep*sigmasep);
                Double_t disc = (Tdiff*Tdiff)/(2*sigmaT*sigmaT)+(sep*sep)/(sigmasep*sigmasep);
		fUserMethods->FillHisto("CHOD_MinDisc_disc", disc);
                if (disc < mindisc)      {
                        d++;
                        mindisc = disc;
                        MinDiscID = i;
                                         }
                                                                        }
//                if(distanceCHOD<99999)fUserMethods->FillHisto("CHOD_MatchedTrack_MinD", distanceCHOD);
                if(mindisc>1.1)return 0;
                if(mindisc<1.1)return 1;
}

