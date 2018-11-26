#ifndef PhotonVetoCandidate_h
#define PhotonVetoCandidate_h

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
//#include "DetectorAcceptance.hh"
#include <TCanvas.h>

class PhotonVetoCandidate 
{
  public :
    PhotonVetoCandidate();
//  virtual ~PhotonVetoCandidate();
 ~PhotonVetoCandidate();
    void Clear();

  public:
    void SetTime(Double_t val)    {fTime=val;};
    void SetID(Int_t val)         {fID=val;};

    Double_t GetTime()       {return fTime;};
    Double_t GetID()         {return fID;};

  private:
    Double_t fTime;
    Int_t fID;
};

#endif
