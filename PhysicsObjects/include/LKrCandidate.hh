#ifndef LKrCandidate_h
#define LKrCandidate_h

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include "DetectorAcceptance.hh"
#include <TCanvas.h>

class LKrCandidate 
{
  public :
    LKrCandidate();
    virtual ~LKrCandidate();
    void Clear();

  public:
    void SetIsLKrCandidate(Bool_t val) {fIsLKrCandidate=val;};
    void SetDiscriminant(Double_t val) {fDiscriminant=val;};
    void SetDeltaTime(Double_t val)    {fDeltaTime=val;   };
    void SetX(Double_t val)            {fX=val;           };
    void SetY(Double_t val)            {fY=val;           };
    void SetEnergy(Double_t val)       {fEnergy=val;      };
    void SetSeedEnergy(Double_t val)   {fSeedEnergy=val;  };
    void SetNCells(Int_t val)          {fNCells=val;      };
    void SetTime(Double_t val)         {fTime=val;        };

    Bool_t   GetIsLKrCandidate() {return fIsLKrCandidate;};
    Double_t GetDiscriminant() {return fDiscriminant;};
    Double_t GetDeltaTime()    {return fDeltaTime;   };
    Double_t GetX()            {return fX;           };
    Double_t GetY()            {return fY;           };
    Double_t GetEnergy()       {return fEnergy;      };
    Double_t GetSeedEnergy()   {return fSeedEnergy;  };
    Int_t    GetNCells()       {return fNCells;      };
    Double_t GetTime()         {return fTime;        };

  private:
    Bool_t fIsLKrCandidate;
    Double_t fDiscriminant;
    Double_t fDeltaTime;
    Double_t fX;
    Double_t fY;
    Double_t fEnergy;
    Double_t fSeedEnergy;
    Int_t fNCells;
    Double_t fTime;
};

#endif
