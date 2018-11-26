#include "PreAnalyzers.hh"
#include "AnalysisTools.hh"
#include "TRecoLKrEvent.hh"
#include "Event.hh"

PreAnalyzers::PreAnalyzers(){
}

void PreAnalyzers::LKrCorrections(TRecoLKrEvent*& fevent){ //use for runs below 1600 (i think)
LKrEvent = fevent;

  for (Int_t iClus=0; iClus<LKrEvent->GetNCandidates(); iClus++) {

    TRecoLKrCandidate* Lcand = (TRecoLKrCandidate*)LKrEvent->GetCandidate(iClus);

Lcand->SetClusterEnergy(Lcand->GetClusterEnergy()/1000);
                                Double_t ue = Lcand->GetClusterEnergy();
                                Double_t cells = Lcand->GetNCells();
                                Double_t ce = LKrZeroSup(cells, ue);
Lcand->SetClusterEnergy(ce);

  }

}
void PreAnalyzers::LKrMToG(TRecoLKrEvent*& fevent){
LKrEvent = fevent;
for (Int_t iClus=0; iClus<LKrEvent->GetNCandidates(); iClus++) {
	TRecoLKrCandidate* Lcand = (TRecoLKrCandidate*)LKrEvent->GetCandidate(iClus);
	Lcand->SetClusterEnergy(Lcand->GetClusterEnergy()/1000);
 }

}
void PreAnalyzers::StrawCorrections(TRecoSpectrometerEvent*& fevent){
SpectrometerEvent = fevent; 
  Double_t alpha  = 2.15e-8; // [MeV^{-1}]
  Double_t beta   = 1.10e-3; // dimensionless

  for (Int_t iTrack=0; iTrack<SpectrometerEvent->GetNCandidates(); iTrack++) {
    TRecoSpectrometerCandidate* Scand =
      (TRecoSpectrometerCandidate*)SpectrometerEvent->GetCandidate(iTrack);

    Double_t Ptrack0 = Scand->GetMomentum();
    Int_t    Q       = Scand->GetCharge();
    Double_t Ptrack  = Ptrack0 * (1+beta) * (1+alpha*Q*Ptrack0);
    Scand->SetMomentum(Ptrack); 
  }
}
