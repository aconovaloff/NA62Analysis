#ifndef Parameters_h
#define Parameters_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "UserMethods.hh"
class Parameters
{
  public :
        Parameters(Int_t);
        ~Parameters();
/*
 but these are private
    Double_t GetEnergy()       {return fEnergy;      };
    Int_t    GetNCells()       {return fNCells;      };
    Int_t    GetNMerged()      {return fNMerged;     };
    Double_t GetEMerged()      {return fEMerged;     };
*/
        Double_t GettCedarHi() {return tCedarHi;};
        Double_t GettCedarLow() {return tCedarLow;};
        double GettLAVLow() {return tLAVLow;};
        double GettLAVHi() {return tLAVHi;};
        double GettMUV3Low() {return tMUV3Low;};
        double GettMUV3Hi() {return tMUV3Hi;};
        Double_t GettCHANTILow() {return tCHANTILow;};
        Double_t GettCHANTIHi() {return tCHANTIHi;};
        Double_t GettMUV1Low() {return tMUV1Low;};
        Double_t GettMUV1Hi() {return tMUV1Hi;};
        Double_t GettMUV2Low() {return tMUV2Low;};
        Double_t GettMUV2Hi() {return tMUV2Hi;};
        Double_t GettRICHLow() {return tRICHLow;};
        Double_t GettRICHHi() {return tRICHHi;};
        Double_t GettLKrLow() {return tLKrLow;};
        Double_t GettLKrHi() {return tLKrHi;};
        Double_t GettStrawLow() {return tStrawLow;};
        Double_t GettStrawHi() {return tStrawHi;};
        Double_t GettSACLow() {return tSACLow;};
        Double_t GettSACHi() {return tSACHi;};
	Double_t GettSAVLow() {return tSAVLow;};
        Double_t GettSAVHi() {return tSAVHi;};
        Double_t GettIRCLow() {return tIRCLow;};
        Double_t GettIRCHi() {return tIRCHi;};
        Double_t GettGTKLow() {return tGTKLow;};
        Double_t GettGTKHi() {return tGTKHi;};
        double GetLKrStartPos() {return LKrStartPos;};
        Double_t GetMUV3StartPos() {return MUV3StartPos;};
        double GetRICHStartPos() {return RICHStartPos;};
        Double_t GetMUV1StartPos() {return MUV1StartPos;};
	Double_t GetMUV2StartPos() {return MUV2StartPos;};
        double GetCHODStartPos() {return CHODStartPos;};
        double GetKaonMass() {return mK;};
        double GetPionMass() {return mPi;};

        Double_t tCedarHi;
        Double_t tCedarLow;
        Double_t tCedarmid;
        double tLAVLow;
        double tLAVHi;
        double tMUV3Low;
        double tMUV3Hi;
        Double_t tCHANTILow;
        Double_t tCHANTIHi;
        Double_t tMUV1Low;
        Double_t tMUV1Hi;
        Double_t tMUV2Low;
        Double_t tMUV2Hi;
        Double_t tRICHLow;
        Double_t tRICHHi;
       Double_t tLKrLow;
        Double_t tLKrHi;
       Double_t tStrawLow;
        Double_t tStrawHi;
       Double_t tSACLow;
        Double_t tSACHi;
      Double_t tSAVLow;
        Double_t tSAVHi;
       Double_t tIRCLow;
        Double_t tIRCHi;
	Double_t tGTKLow;
	Double_t tGTKHi;
        double LKrStartPos;
        Double_t MUV3StartPos;
        double RICHStartPos;
        Double_t MUV1StartPos;
	Double_t MUV2StartPos;
        double CHODStartPos;
        double mK;
        double mPi;

  private:

};

#endif

