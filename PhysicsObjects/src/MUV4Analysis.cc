#include "MUV4Analysis.hh"
#include "AnalysisTools.hh"
#include "TRecoMUV2Event.hh"
#include "Event.hh"

MUV4Analysis::MUV4Analysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}
void MUV4Analysis::BookHistos(){
             fUserMethods->BookHisto(new TH1F("a", "Missing Mass squared after MUV3 veto", 5, 0, 5));

}
void MUV4Analysis::FillHistos(){
fUserMethods->FillHisto("a", 1); 
}

void MUV4Analysis::SaveHistos(){
fUserMethods->SaveAllPlots();
}

