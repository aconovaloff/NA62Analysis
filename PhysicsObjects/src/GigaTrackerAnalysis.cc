#include "GigaTrackerAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoGigaTrackerEvent.hh"
#include "Event.hh"
#include "UserMethods.hh"
using namespace TMath; 

GigaTrackerAnalysis::GigaTrackerAnalysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}

void GigaTrackerAnalysis::Set(TRecoGigaTrackerEvent* fevent, double freftime, double ftLow, double ftHi, bool fMCflag ){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
MCflag = fMCflag; 
TRecoGigaTrackerCandidate* GTKCand;
for(int i=0; i<event->GetNCandidates(); i++)         {
GTKCand = (TRecoGigaTrackerCandidate*)event->GetCandidate(i);
Double_t Tdiff = reftime - GTKCand->GetTime();
fUserMethods->FillHisto("GTK_Set_reft", Tdiff);
}
K3.SetXYZ(Sin(0.0012)*75, 0, Cos(0.0012)*75);
KPos.SetXYZ(0,0,101800);
K4.SetPxPyPzE(Sin(0.0012)*75, 0, Cos(0.0012)*75, sqrt(75*75+.493677*.493677));
/*
         double trim = 0.0012;
                TVector3 KPos(0, 0, 101800);
                TVector3 K3(Sin(trim)*75, 0, Cos(trim)*75);
                TLorentzVector K4(Sin(trim)*75, 0, Cos(trim)*75, sqrt(75*75+.493677*.493677));

*/
}

Int_t GigaTrackerAnalysis::MatchedTrack(TVector3 PionMom, TVector3 PionPos, Int_t RunNumber){
TVector3 PionMomentum = PionMom;
TVector3 PionPosition = PionPos;
TRecoGigaTrackerCandidate *GTKCand;
Double_t close = 999999;
Int_t minID_t;
Int_t minID_cda;
Int_t min_t = 0;
Int_t min_cda = 0;
Double_t Tdiff_min = 99999;
Double_t min_d = 99999;
Double_t min_d2 = 99999;
Double_t min_d123 = 99999;
Int_t minID_d=-1;
Int_t minID_d2;
Int_t minID_d123;
Int_t dis = 0;
Int_t dis2 = 0;
Int_t dis123 = 0;
Int_t goodCand = 0;
Double_t sigmaCDA; //two GTK stations
Double_t sigmaCDA123;
Double_t sigmaT;
///how many burst needed ????? do one burst, then after full run, fine tune, doesn't look like it changes it much
if (RunNumber==3809){
sigmaCDA=2.601; ////////not sure what run, prob 3809
sigmaCDA123=1.572; //////not sure what run, prob 3809
sigmaT=0.3365;
}
else if (RunNumber==3821){
sigmaCDA=2.31; //run 3821
sigmaCDA123=1.63; //run 3821
//Double_t sigmaT=0.3565; //not sure what run, prob 3809
sigmaT=0.3038;//run 3821
}
else if (RunNumber==3810){
sigmaCDA=2.579;
sigmaCDA123=1.556;
sigmaT=0.3366;
}
else if (RunNumber==3811){
sigmaCDA=2.612;
sigmaCDA123=1.567;
sigmaT=0.3404;
}
else if (RunNumber==3813){
sigmaCDA=2.673;
sigmaCDA123=1.695;
sigmaT=0.3426;
}

else if (RunNumber==3818){
sigmaCDA=2.280;
sigmaCDA123=1.677;
sigmaT=0.3478;
}

else {
sigmaCDA=2.7; ////////not sure what run, prob 3809
sigmaCDA123=1.69; //////not sure what run, prob 3809
sigmaT=0.3565;
}
for(int i=0; i<event->GetNCandidates(); i++)         {
GTKCand = (TRecoGigaTrackerCandidate*)event->GetCandidate(i);
Double_t Tdiff = reftime - GTKCand->GetTime();
TVector3 Kmom(GTKCand->GetMomentum().X()/1000, GTKCand->GetMomentum().Y()/1000, GTKCand->GetMomentum().Z()/1000);
if (Kmom.Mag()<72||Kmom.Mag()>77.5)continue;
Double_t xdiff = GTKCand->GetPosition(0).X()- GTKCand->GetPosition(2).X();
Double_t ydiff = GTKCand->GetPosition(0).Y()- GTKCand->GetPosition(2).Y();
//Double_t zdiff= GTKCand->GetPosition(0).Z()- GTKCand->GetPosition(2).Z();
Double_t zdiff=22800; 
Double_t zeta = ACos(zdiff/sqrt(zdiff*zdiff+xdiff*xdiff));
Double_t eta = ACos(zdiff/sqrt(zdiff*zdiff+ydiff*ydiff));
fUserMethods->FillHisto("zeta", zeta);
fUserMethods->FillHisto("eta", eta);
if (zeta>0.000200)continue;
if(eta>0.000200)continue; 

//TVector3 Kpos(GTKCand->GetPosition(0).X(), GTKCand->GetPosition(0).Y(), GTKCand->GetPosition(0).Z());
//TVector3 Kpos(GTKCand->GetPosition(0).X(), GTKCand->GetPosition(0).Y(),79600);
//TVector3 Kpos(GTKCand->GetPosition(1).X(), GTKCand->GetPosition(1).Y(),92800);
TVector3 Kpos(GTKCand->GetPosition(2).X(), GTKCand->GetPosition(2).Y(),102400);
 if(fabs(Tdiff) < Tdiff_min)     {
                        Tdiff_min = Tdiff;
                        minID_t = i;
                        min_t++;
                                 }
Double_t cldiap = cda(PionMomentum, Kmom, PionPosition,Kpos);
fUserMethods->FillHisto("cldiap", cldiap);
fUserMethods->FillHisto("cldiap_mirror", cldiap);
fUserMethods->FillHisto("cldiap_mirror", -cldiap);
fUserMethods->FillHisto("GTK_MatchedTrack_disc", Tdiff, cldiap);
if (cldiap<close){
close = cldiap;
minID_cda = i;
min_cda++;
}
/////could get sigmas from "best" cdas later, not sure if this will do anything (e.g., don't do it for timing)
if(GTKCand->GetType()==12||GTKCand->GetType()==13||GTKCand->GetType()==23){
fUserMethods->FillHisto("cldiap2", cldiap);
fUserMethods->FillHisto("cldiap2", -cldiap);
}
if(GTKCand->GetType()==123){
fUserMethods->FillHisto("cldiap123", cldiap);
fUserMethods->FillHisto("cldiap123", -cldiap);
}

if(!MCflag){
if (Tdiff<0){ if(cldiap>(20*Tdiff+20))continue;}
if (Tdiff>0){if(cldiap>(20-20*Tdiff))continue;}
}
else {if(cldiap>20)continue;} 
fUserMethods->FillHisto("GTK_MatchedTrack_goodtime", Tdiff);
fUserMethods->FillHisto("GTK_MatchedTrack_goodcda", cldiap);
fUserMethods->FillHisto("GTK_MatchedTrack_goodcda_mirror", cldiap);
fUserMethods->FillHisto("GTK_MatchedTrack_goodcda_mirror", -cldiap);
fUserMethods->FillHisto("GTK_MatchedTrack_disc_post", Tdiff, cldiap);
goodCand++;

if(GTKCand->GetType()==12||GTKCand->GetType()==13||GTKCand->GetType()==23){
double d2;
if(!MCflag) d2= sqrt(Tdiff*Tdiff/(sigmaT*sigmaT) + (cldiap/sigmaCDA)*(cldiap/sigmaCDA));
else d2 = sqrt(cldiap*cldiap); 
if (d2<min_d2){
min_d2 = d2;
minID_d2 = i;
dis2++;
}
}

if(GTKCand->GetType()==123){
double d123;
if(!MCflag) d123 = sqrt(Tdiff*Tdiff/(sigmaT*sigmaT) + (cldiap/sigmaCDA123)*(cldiap/sigmaCDA123));
else d123 = sqrt(cldiap*cldiap); 

if (d123<min_d123){
min_d123 = d123;
minID_d123 = i;
dis123++;
}
}
/*
double d = sqrt(Tdiff*Tdiff/(0.3565*0.3565) + (cldiap/2.64)*(cldiap/2.64));
if (d<min_d){
min_d = d;
minID_d = i;
dis++;
}
*/
}
fUserMethods->FillHisto("close", close);
fUserMethods->FillHisto("close_mirror", close);
fUserMethods->FillHisto("close_mirror", -close);
fUserMethods->FillHisto("GTK_MatchedTrack_goodCand", goodCand);
fUserMethods->FillHisto("GTK_MatchedTrack_disc_best", Tdiff_min, close);
if (dis2==0&&dis123==0) return 0; 
if (dis123>0)minID_d = minID_d123;
else minID_d = minID_d2;
//else return 0;  //for g's plot
//if (dis==0) return 0;
GTKCand = (TRecoGigaTrackerCandidate*)event->GetCandidate(minID_d);
//cout<<(GTKCand->GetPosition(1).Theta()-GTKCand->GetPosition(2).Theta());
K3.SetXYZ(GTKCand->GetMomentum().X()/1000, GTKCand->GetMomentum().Y()/1000, GTKCand->GetMomentum().Z()/1000);
KPos.SetXYZ(GTKCand->GetPosition(2).X(), GTKCand->GetPosition(2).Y(), /*GTKCand->GetPosition(2).Z()*/102400);
K4.SetPxPyPzE(GTKCand->GetMomentum().X()/1000, GTKCand->GetMomentum().Y()/1000, GTKCand->GetMomentum().Z()/1000,sqrt(K3*K3+.493677*.493677));
Double_t cda_final = cda(PionMomentum, K3, PionPosition,KPos);
fUserMethods->FillHisto("GTK_MatchedTrack_cda_final",cda_final);

return 1; 
}

void GigaTrackerAnalysis::BookHistos_MatchedTrack(){
        TH2I* GTK_MatchedTrack_disc = new TH2I("GTK_MatchedTrack_disc", "cda (y-axis), time (x-axis)", 50, -5, 5, 50, 0, 100);
        GTK_MatchedTrack_disc->SetOption("COLZ");
        GTK_MatchedTrack_disc->Draw();
        fUserMethods->BookHisto(GTK_MatchedTrack_disc);
        TH2I* GTK_MatchedTrack_disc_post = new TH2I("GTK_MatchedTrack_disc_post", "cda (y-axis), time (x-axis)", 50, -5, 5, 50, 0, 100);
        GTK_MatchedTrack_disc_post->SetOption("COLZ");
        GTK_MatchedTrack_disc_post->Draw();
        fUserMethods->BookHisto(GTK_MatchedTrack_disc_post);
        TH2I* GTK_MatchedTrack_disc_best = new TH2I("GTK_MatchedTrack_disc_best", "cda vs t for gtk", 50, -5, 5, 50, 0, 100);
        GTK_MatchedTrack_disc_best->SetOption("COLZ");
        GTK_MatchedTrack_disc_best->Draw();
        fUserMethods->BookHisto(GTK_MatchedTrack_disc_best);
        fUserMethods->BookHisto(new TH1F("GTK_MatchedTrack_goodCand", "good GTK candidates", 11, 0, 11));
	fUserMethods->BookHisto(new TH1F("GTK_MatchedTrack_goodtime", "GTK reference time", 16, -2, 2));
	fUserMethods->BookHisto(new TH1F("GTK_MatchedTrack_goodcda", "cda", 120, 0, 30));
        fUserMethods->BookHisto(new TH1F("GTK_MatchedTrack_goodcda_mirror", "cda", 240, -30, 30));
	fUserMethods->BookHisto(new TH1F("zeta", "zeta", 100, 0,0.000500));
	fUserMethods->BookHisto(new TH1F("eta", "eta", 100, 0, 0.000500));
        fUserMethods->BookHisto(new TH1F("cldiap", "cda", 100, 0, 50));
        fUserMethods->BookHisto(new TH1F("cldiap_mirror", "cda", 200,-50, 50));
        fUserMethods->BookHisto(new TH1F("close", "cda", 100, 0, 50));
        fUserMethods->BookHisto(new TH1F("close_mirror", "cda", 200,-50, 50));
        fUserMethods->BookHisto(new TH1F("cldiap2", "cda", 200,-50, 50));
        fUserMethods->BookHisto(new TH1F("cldiap123", "cda", 200,-50, 50));
	fUserMethods->BookHisto(new TH1F("GTK_MatchedTrack_cda_final", "cda", 200, 0, 100));
}
void GigaTrackerAnalysis::BookHistos_Set(){
        fUserMethods->BookHisto(new TH1F("GTK_Set_reft", "GTK reference time", 200, -50, 50));

}
void GigaTrackerAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
}

TVector3 GigaTrackerAnalysis::GetKaonMomentum(){
return K3;
}

TVector3 GigaTrackerAnalysis::GetKaonPosition(){
return KPos;
}

TLorentzVector GigaTrackerAnalysis::GetKaonFourMomentum(){
return K4;
}

