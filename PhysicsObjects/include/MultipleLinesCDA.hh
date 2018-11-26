// ---------------------------------------------------------------
// History:
//
// Created by Plamen Petrov (plpetrov@cern.ch) 2015-04-20
//
// ---------------------------------------------------------------

#ifndef MULTIPLELINESCDA_H
#define MULTIPLELINESCDA_H

#include "TVector3.h"
#include "TClonesArray.h"


/// \class MultipleLinesCDA
/// \Brief
/// Computation of CDAs for multiple pairs of lines and a common vertex 
/// \EndBrief
/// \Detailed
/// The method is based on TwoLinesCDA. It calculates the CDA between each pair 
/// of tracks and their coordinates at the CDA. The method returns also a vertex estimate 
/// based on the mean of the centers of all CDA vectors.
/// \EndDetailed

class MultipleLinesCDA {

public:

  MultipleLinesCDA();
  ~MultipleLinesCDA();

  void ComputeVertex();
  void Print();
  void Reset();

  void AddLinePoint1Dir (TVector3 p1, TVector3 dir);
  void RemoveLine(Int_t n);

  void SetLinePoint1 (Int_t n, Double_t x, Double_t y, Double_t z);
  void SetLineDir    (Int_t n, Double_t x, Double_t y, Double_t z);

  TVector3 GetLinePoint1 (Int_t n);
  TVector3 GetLineDir(Int_t n);

  TVector3 GetVertex()             { return fVertex; }
  void     SetVertex(TVector3 val) { fVertex = val;  }
  Double_t GetCDA(Int_t n1, Int_t n2);
  Int_t    GetNumberOfLines()      { return fNLines; }
 
private:

  TClonesArray fLineP1;
  TClonesArray fDir;
  TVector3 fVertex;
  Int_t fNLines;

};

#endif
