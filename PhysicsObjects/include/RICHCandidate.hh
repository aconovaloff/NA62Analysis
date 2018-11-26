#ifndef RICHCandidate_h
#define RICHCandidate_h

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include "DetectorAcceptance.hh"
#include <TCanvas.h>

class RICHCandidate 
{
  public :
    RICHCandidate();
    virtual ~RICHCandidate();
    void Clear();

  public:
    void SetIsRICHCandidate(Bool_t val) {fIsRICHCandidate=val;};
    void SetDiscriminant(Double_t val) {fDiscriminant=val;};
    void SetDX(Double_t val)           {fDX=val;};
    void SetDY(Double_t val)           {fDY=val;};
    void SetXcenter(Double_t val)      {fXcenter=val;};
    void SetYcenter(Double_t val)      {fYcenter=val;};
    void SetRadius(Double_t val)       {fRadius=val;};
    void SetDeltaTime(Double_t val)    {fDeltaTime=val;};
    void SetChi2(Double_t val)         {fChi2=val;};
    void SetNHits(Int_t val)           {fNHits=val;};
    void SetMirrorID(Int_t val)        {fMirrorID=val;};
  void SetTimeCandInd(Int_t val)        {fTimeCandInd=val;};
  void SetHitsMaxAngle(Int_t val)        {fHitsMaxAngle=val;};




    Bool_t   GetIsRICHCandidate() {return fIsRICHCandidate;};
    Double_t GetDiscriminant() {return fDiscriminant;};
    Double_t GetDX()           {return fDX;};
    Double_t GetDY()           {return fDY;};
    Double_t GetXcenter()      {return fXcenter;};
    Double_t GetYcenter()      {return fYcenter;};
    Double_t GetRadius()       {return fRadius;};
    Double_t GetDeltaTime()    {return fDeltaTime;};
    Double_t GetChi2()         {return fChi2;};
    Int_t GetNHits()           {return fNHits;};
    Int_t GetMirrorID()        {return fMirrorID;};
  Int_t GetTimeCandInd()       {return fTimeCandInd;}; 
  Double_t GetHitsMaxAngle(Int_t val)        {return fHitsMaxAngle;};


  private:
    Bool_t fIsRICHCandidate;
    Double_t fDiscriminant;
    Double_t fDX;
    Double_t fDY;
    Double_t fXcenter;
    Double_t fYcenter;
    Double_t fRadius;
    Double_t fDeltaTime;
    Double_t fChi2;
    Int_t fNHits;
    Int_t fMirrorID;
  Int_t fTimeCandInd;
  Double_t fHitsMaxAngle;
};

#endif
