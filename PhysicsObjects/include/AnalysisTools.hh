#ifndef AnalysisTools_HH
#define AnalysisTools_HH  

#include <stdlib.h>
#include <vector>

#include <iostream>
#include <TChain.h>

//#include "MCSimple.hh"
//#include "functions.hh"


#include "Event.hh"
//#include "Persistency.hh"

//#include "Analyzer.hh"

#include <TCanvas.h>
using namespace std;



//class AnalysisTools
//{

//public:
// AnalysisTools();

//public:

  TVector3 Propagate(Int_t, TVector3*, TVector3*, Double_t);
  TVector3 Vertex(TVector3, TVector3, TVector3, TVector3);
  TVector3 prop(TVector3, TVector3, double);
  TVector3 prop_neg(TVector3, TVector3, double);
  Double_t cda(TVector3, TVector3, TVector3, TVector3);
  Double_t LKrZeroSup(Double_t, Double_t);
TVector3 backprop(TVector3, TVector3, double);
//TVector3 BlueTubeProp(TVector3, TVector3, Double_t);   

//private:



//};
#endif


                                                                      
