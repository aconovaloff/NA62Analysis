#include "AnalysisTools_G.hh"
#include "Parameters_G.hh"
#include "LKrAnalysis_G.hh"

LKrAnalysis_G::LKrAnalysis_G(Int_t flag, NA62Analysis::Core::BaseAnalysis *ba) { 

  fFlag = flag;
  fTools = AnalysisTools_G::GetInstance();
  fUserMethods = new UserMethods(ba);

  if (flag==0) { // Track - cluster matching
//    Parameters_G *par = Parameters_G::GetInstance(); 
//    par->LoadParameters_G("/afs/cern.ch/user/r/ruggierg/workspace/public/run2015/database/lkrt0.dat");
//    fLKrCellT0 = (Double_t **)par->GetLKrCellT0(); 
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_distvst_cluster","",200,-100,100,300,0,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_distvstchod_cluster","",400,-100,100,300,0,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_xy_cluster","",300,-300,300,300,-300,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_mindistvst_cluster","",200,-100,100,300,0,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_mindistvstchod_cluster","",400,-100,100,300,0,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_minxy_cluster","",300,-300,300,300,-300,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_energyvsdist_cluster","",300,0,300,400,0,100));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_mindistvst_cell","",200,-100,100,300,0,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_mindistvstchod_cell","",400,-100,100,300,0,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_minxy_cell","",300,-300,300,300,-300,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_energyvsdist_cell","",300,0,300,400,0,100));
    fLKrClusterCandidate = new LKrCandidate();
    fLKrCellCandidate = new LKrCandidate();
  }

  if (flag==1) { // Extra clusters finding
//    Parameters_G *par = Parameters_G::GetInstance(); 
//    par->LoadParameters_G("/afs/cern.ch/user/r/ruggierg/workspace/public/run2015/database/lkrt0.dat");
//    fLKrCellT0 = (Double_t **)par->GetLKrCellT0(); 
//    par->LoadParameters_G(t0finename.Data());
//    par->StoreLKrParameters_G();
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_not_matched","",200,0,2000,1200,-150,150));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_not_matched_all","",200,0,2000,1200,-150,150));
    fUserMethods->BookHisto(new TH1F("LKrAnalysis_G_not_matched_highe","",500,-50,50));
    fUserMethods->BookHisto(new TH1F("LKrAnalysis_G_not_matched_lowe","",500,-50,50));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_not_matched_photon","",200,0,2000,1200,-150,150));
    for (Int_t kk=0; kk<10; kk++) fPhotonCandidate[kk] = new LKrCandidate();
    for (Int_t kk=0; kk<5; kk++) fNewClusterCandidate[kk] = new LKrCandidate();
  }

  if (flag==2) { // Two photon analysis 
//    Parameters_G *par = Parameters_G::GetInstance(); 
//    par->LoadParameters_G("/afs/cern.ch/user/r/ruggierg/workspace/public/run2015/database/lkrt0.dat");
//    fLKrCellT0 = (Double_t **)par->GetLKrCellT0(); 
//    par->LoadParameters_G(t0finename.Data());
//    par->StoreLKrParameters_G();
    fTwoEMClustersID = new Int_t[45];
    fUserMethods->BookHisto(new TH2F("hbo","",40,0,40,400,-200,200));
    fUserMethods->BookHisto(new TH1F("Clusters_n","",10,0,10));
    fUserMethods->BookHisto(new TH1F("Clusters_all_time","",400,0,500));
    fUserMethods->BookHisto(new TH2F("NCellsvsEClusters_raw","",300,0,100,150,0,150));
    fUserMethods->BookHisto(new TH2F("NCellEnergyratiovsSeedEnergyRatio_raw","",110,0,1.1,200,0,20));
    fUserMethods->BookHisto(new TH1F("MipClusterEnergy","",200,0,100));
    fUserMethods->BookHisto(new TH1F("MipClusterEnergy_zoom","",100,0,2));
    fUserMethods->BookHisto(new TH2F("MipClusterXY","",128,-1263.2,1263.2,128,-1263.2,1263.2));
    fUserMethods->BookHisto(new TH2F("NCellEnergyratiovsSeedEnergyRatio_nomip","",110,0,1.1,400,-20,20));
    fUserMethods->BookHisto(new TH1F("EMClusterEnergy","",200,0,100));
    fUserMethods->BookHisto(new TH2F("EMClusterXY","",128,-1263.2,1263.2,128,-1263.2,1263.2));
    fUserMethods->BookHisto(new TH2F("EMClusterXY_high_pnn","",128,-1263.2,1263.2,128,-1263.2,1263.2));
    fUserMethods->BookHisto(new TH2F("EMClusterXY_high_ctr","",128,-1263.2,1263.2,128,-1263.2,1263.2));
    fUserMethods->BookHisto(new TH2F("EMClusterXY_low_pnn","",128,-1263.2,1263.2,128,-1263.2,1263.2));
    fUserMethods->BookHisto(new TH2F("EMClusterXY_low_ctr","",128,-1263.2,1263.2,128,-1263.2,1263.2));
    fUserMethods->BookHisto(new TH1F("OtherClusterEnergy","",200,0,100));
    fUserMethods->BookHisto(new TH2F("OtherClusterXY","",128,-1263.2,1263.2,128,-1263.2,1263.2));
    fUserMethods->BookHisto(new TH1F("Clusters_nMips","",10,0,10));
    fUserMethods->BookHisto(new TH1F("Clusters_nEM","",10,0,10));
    fUserMethods->BookHisto(new TH1F("Clusters_nOther","",10,0,10));
    fUserMethods->BookHisto(new TH1F("TwoEMClusterTime","",200,-10,10));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_2gammadistvst_cluster","",200,-100,100,300,0,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_2gammadistvstchod_cluster","",400,-100,100,300,0,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_2gammaxy_cluster","",300,-300,300,300,-300,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_2gammamindistvst_cluster","",200,-100,100,300,0,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_2gammamindistvstchod_cluster","",400,-100,100,300,0,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_2gammaminxy_cluster","",300,-300,300,300,-300,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_2gammaenergyvsdist_cluster","",300,0,300,400,0,100));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_2gammamindistvst_cell","",200,-100,100,300,0,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_2gammamindistvstchod_cell","",400,-100,100,300,0,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_2gammaminxy_cell","",300,-300,300,300,-300,300));
    fUserMethods->BookHisto(new TH2F("LKrAnalysis_G_2gammaenergyvsdist_cell","",300,0,300,400,0,100));
    fLKrClusterCandidate = new LKrCandidate();
    fLKrCellCandidate = new LKrCandidate();
  }

}

LKrAnalysis_G::~LKrAnalysis_G()
{ }

void LKrAnalysis_G::StartBurst(Int_t year, Int_t runid) {
  Parameters_G *par = Parameters_G::GetInstance(); 
  TString t0finename;
  t0finename.Form("/afs/cern.ch/user/r/ruggierg/workspace/public/na62git/database/run%d/lkrt0.dat",runid); 
  par->LoadParameters_G(t0finename.Data()); 
  par->StoreLKrParameters_G();
  fLKrCellT0 = (Double_t **)par->GetLKrCellT0(); 
  cout << t0finename << " loaded" << endl; 
  fYear = year;
}

void LKrAnalysis_G::Clear(Int_t flag) {

  if (flag==0) {
    fLKrClusterCandidate->Clear();
    fLKrCellCandidate->Clear(); 
  }

  if (flag==1) {
    for (Int_t kk=0; kk<10; kk++) fPhotonCandidate[kk]->Clear(); 
    for (Int_t kk=0; kk<5; kk++) fNewClusterCandidate[kk]->Clear(); 
  }

  if (flag==2) {
    fNEMClusters = 0;
    fNMIPClusters = 0;
    fNOtherClusters = 0;
    fill(fEMClusterID,fEMClusterID+10,-1);
    fill(fMIPClusterID,fMIPClusterID+10,-1);
    fill(fOtherClusterID,fOtherClusterID+10,-1);
    fNEMDoublets = 0;
    for (Int_t j=0; j<45; j++) fTwoEMClustersID[j]=-1;
    fLKrClusterCandidate->Clear();
    fLKrCellCandidate->Clear(); 
  }

}

///////////////////////////////////////////////////////
// Cluster analysis for pip0 selection using photons //
///////////////////////////////////////////////////////
Int_t LKrAnalysis_G::PhotonAnalysis(Bool_t mcflag, Int_t runid, TRecoLKrEvent* event, L0TPData* l0tp) {
  fL0Data = l0tp;
  fUserMethods->FillHisto("Clusters_n",event->GetNCandidates());
  ClusterAnalysis(mcflag,runid,event);
  fUserMethods->FillHisto("Clusters_nMips",fNMIPClusters);
  fUserMethods->FillHisto("Clusters_nEM",fNEMClusters);
  fUserMethods->FillHisto("Clusters_nOther",fNOtherClusters);
  Int_t fIsTwoPhotonEvent = 0;
  if (fNEMClusters>=2) fIsTwoPhotonEvent = TwoPhotonCandidate(event);
  return fIsTwoPhotonEvent;
} 
void LKrAnalysis_G::ClusterAnalysis(Bool_t mcflag, Int_t runid, TRecoLKrEvent* event) {
  for (Int_t jclus=0; jclus<event->GetNCandidates(); jclus++) {
    TRecoLKrCandidate *thisCluster = (TRecoLKrCandidate *)event->GetCandidate(jclus);
    if (!mcflag) fTools->ClusterCorrections(runid,thisCluster);
    fUserMethods->FillHisto("Clusters_all_time",thisCluster->GetClusterTime());
    fUserMethods->FillHisto("NCellsvsEClusters_raw",thisCluster->GetClusterEnergy()/1000.,thisCluster->GetNCells());
    fUserMethods->FillHisto("NCellEnergyratiovsSeedEnergyRatio_raw",thisCluster->GetClusterSeedEnergy()/thisCluster->GetClusterEnergy(),thisCluster->GetNCells()/(thisCluster->GetClusterEnergy()/1000.));
    Bool_t isMIPCluster = SelectMIPCluster(jclus,thisCluster);
    Bool_t isEMCluster = SelectEMCluster(jclus,thisCluster);
    Bool_t isOtherCluster = SelectOtherCluster(jclus,thisCluster,isMIPCluster,isEMCluster);
  }
}
Bool_t LKrAnalysis_G::TwoPhotonCandidate(TRecoLKrEvent* event) {
  for (Int_t iemclus=0; iemclus<fNEMClusters-1; iemclus++) {
    TRecoLKrCandidate *iCluster = (TRecoLKrCandidate *)event->GetCandidate(fEMClusterID[iemclus]);
    if (CloseCluster(fEMClusterID[iemclus],iCluster->GetClusterX(),iCluster->GetClusterY(),event)) continue;
    for (Int_t jemclus=iemclus+1; jemclus<fNEMClusters; jemclus++) {
      TRecoLKrCandidate *jCluster = (TRecoLKrCandidate *)event->GetCandidate(fEMClusterID[jemclus]);
      if (CloseCluster(fEMClusterID[jemclus],jCluster->GetClusterX(),jCluster->GetClusterY(),event)) continue;
      Double_t dtime = iCluster->GetClusterTime()-jCluster->GetClusterTime();
      fUserMethods->FillHisto("TwoEMClusterTime",dtime);  
      if (fabs(dtime)>2.5) continue;
      fTwoEMClustersID[2*fNEMDoublets+0] = fEMClusterID[iemclus];
      fTwoEMClustersID[2*fNEMDoublets+1] = fEMClusterID[jemclus];
      fNEMDoublets++;
      if (fNEMDoublets==45) return 1;
    }   
  }   
  return (fNEMDoublets)?1:0;
}
Bool_t LKrAnalysis_G::SelectMIPCluster(Int_t jclus, TRecoLKrCandidate *thisCluster) {
  if (thisCluster->GetNCells()>=6) return 0;
  if (thisCluster->GetClusterDDeadCell()<20.) return 0;
  Double_t seedRatio = thisCluster->GetClusterSeedEnergy()/thisCluster->GetClusterEnergy();
  Double_t cellRatio = thisCluster->GetNCells()/(thisCluster->GetClusterEnergy()/1000.);
  if (cellRatio<0.2) return 0;
  fUserMethods->FillHisto("MipClusterEnergy",thisCluster->GetClusterEnergy()/1000.);
  fUserMethods->FillHisto("MipClusterEnergy_zoom",thisCluster->GetClusterEnergy()/1000.);
  fUserMethods->FillHisto("MipClusterXY",thisCluster->GetClusterX(),thisCluster->GetClusterY());
  if (fNMIPClusters<10) SetMIPClusterID(jclus);
  AddMIPClusters();
  return 1;
}
Bool_t LKrAnalysis_G::SelectEMCluster(Int_t jclus, TRecoLKrCandidate *thisCluster) {
  if (thisCluster->GetClusterEnergy()<3000) return 0;
  if (thisCluster->GetClusterDDeadCell()<20) return 0;
  if (!fTools->LKrAcceptance(thisCluster->GetClusterX(),thisCluster->GetClusterY(),150.,1130.)); 
  Double_t seedRatio = thisCluster->GetClusterSeedEnergy()/thisCluster->GetClusterEnergy();
  Double_t cellRatio = (thisCluster->GetNCells()-2.21)/(0.93*thisCluster->GetClusterEnergy()/1000.);
  fUserMethods->FillHisto("NCellEnergyratiovsSeedEnergyRatio_nomip",seedRatio,cellRatio); 
  Bool_t isNotEM = true;
  if (seedRatio>0.23&&seedRatio<0.46&&cellRatio<2&&cellRatio>=1) isNotEM = false;
  if (isNotEM) return 0; 
  fUserMethods->FillHisto("EMClusterEnergy",thisCluster->GetClusterEnergy()/1000.);

  // Debug LKr trigger
  fUserMethods->FillHisto("EMClusterXY",thisCluster->GetClusterX(),thisCluster->GetClusterY());
  if (thisCluster->GetClusterEnergy()/1000.>25) {
    if ((fL0Data->GetDataType()>>0)&1 && (fL0Data->GetTriggerFlags()>>1)&1) { // pnn trigger
      fUserMethods->FillHisto("EMClusterXY_high_pnn",thisCluster->GetClusterX(),thisCluster->GetClusterY());
////      Double_t xclus = thisCluster->GetClusterX();
////      Double_t yclus = thisCluster->GetClusterY();
////      if (xclus<-160.&&xclus>-240.&&yclus<0&&yclus>-160) {
////        Int_t idcell = (Int_t)((xclus+240.)/19.74)+(Int_t)(4*((yclus+160)/19.74));
////        if (fL0Data->GetPrimitive(0,6).GetPrimitiveID()!=0) {
////          Double_t primtime = fL0Data->GetPrimitive(0,6).GetFineTime();
////          Double_t dtime = thisCluster->GetClusterTime()-primtime*25./256.;
////          fUserMethods->FillHisto("hbo",idcell,dtime);
////        }
////      }
    }
    if ((fL0Data->GetDataType()>>4)&1) { // control trigger
      fUserMethods->FillHisto("EMClusterXY_high_ctr",thisCluster->GetClusterX(),thisCluster->GetClusterY());
    }
  }
  if (thisCluster->GetClusterEnergy()/1000.<16) {
    if ((fL0Data->GetDataType()>>0)&1 && (fL0Data->GetTriggerFlags()>>1)&1) { // pnn trigger
      fUserMethods->FillHisto("EMClusterXY_low_pnn",thisCluster->GetClusterX(),thisCluster->GetClusterY());
    }
    if ((fL0Data->GetDataType()>>4)&1) { // control trigger
      fUserMethods->FillHisto("EMClusterXY_low_ctr",thisCluster->GetClusterX(),thisCluster->GetClusterY());
    }
  }

  if (fNEMClusters<10) SetEMClusterID(jclus);
  AddEMClusters();
  return 1;
}
Bool_t LKrAnalysis_G::SelectOtherCluster(Int_t jclus, TRecoLKrCandidate *thisCluster, Bool_t isMIPCluster, Bool_t isEMCluster) {
  if (isMIPCluster) return 0;
  if (isEMCluster) return 0;
  fUserMethods->FillHisto("OtherClusterEnergy",thisCluster->GetClusterEnergy()/1000.);
  fUserMethods->FillHisto("OtherClusterXY",thisCluster->GetClusterX(),thisCluster->GetClusterY());
  if (fNOtherClusters<10) SetOtherClusterID(jclus);
  AddOtherClusters();
  return 1;
}
Bool_t LKrAnalysis_G::CloseCluster(Int_t clid, Double_t xlkr, Double_t ylkr, TRecoLKrEvent* event) {
  for (Int_t jclus=0; jclus<event->GetNCandidates(); jclus++) {
    if (jclus==clid) continue;
    TRecoLKrCandidate *clus = (TRecoLKrCandidate *)event->GetCandidate(jclus);
    Double_t thisx = clus->GetClusterX();
    Double_t thisy = clus->GetClusterY();
    Double_t dist = (xlkr-thisx)*(xlkr-thisx)+(ylkr-thisy)*(ylkr-thisy);
    if (dist<150.*150.) return 1;
  }
  return 0;
}

//////////////////////////////////////////////////
// Track - cluster matching for pinunu analysis //
//////////////////////////////////////////////////
Int_t LKrAnalysis_G::TrackClusterMatching(Bool_t mcflag, Int_t runid, Double_t chodTime, Double_t trackTime, TRecoLKrEvent* event, TVector3 *posatlkr) {
  fLKrEvent = event;
  Double_t xpart = posatlkr->X();
  Double_t ypart = posatlkr->Y();
  Int_t minid = -1;

  // Look for clusters matching the track
  Double_t mindist = 99999999.;
  Double_t mintime = -99999.;
  Double_t mintimechod = 9999999.;
  Double_t mindx = -99999.;
  Double_t mindy = -99999.;
  for (Int_t jclus=0; jclus<fLKrEvent->GetNCandidates(); jclus++) {
    TRecoLKrCandidate *thisCluster = (TRecoLKrCandidate *)fLKrEvent->GetCandidate(jclus);
    if (!mcflag) fTools->ClusterCorrections(runid,thisCluster);
    TVector2 posclust(thisCluster->GetClusterX(),thisCluster->GetClusterY());
    TVector2 postrack(xpart,ypart);
    Double_t dist = (posclust-postrack).Mod();
    Double_t dtime = thisCluster->GetClusterTime()-trackTime-1.2;
    Double_t dtimechod = thisCluster->GetClusterTime()-chodTime-1.2;
    Double_t dx = (posclust-postrack).X();
    Double_t dy = (posclust-postrack).Y();
    if (fFlag==0) {
      fUserMethods->FillHisto("LKrAnalysis_G_distvst_cluster",dtime,dist);
      fUserMethods->FillHisto("LKrAnalysis_G_distvstchod_cluster",dtimechod,dist);
      fUserMethods->FillHisto("LKrAnalysis_G_xy_cluster",dx,dy);
    }
    if (fFlag==2) {
      fUserMethods->FillHisto("LKrAnalysis_G_2gammadistvst_cluster",dtime,dist);
      fUserMethods->FillHisto("LKrAnalysis_G_2gammadistvstchod_cluster",dtimechod,dist);
      fUserMethods->FillHisto("LKrAnalysis_G_2gammaxy_cluster",dx,dy);
    }
//    if (dist<mindist && fabs(dtime)<25. && dtimechod<15 && dtimechod>-10) {
    if (dist<mindist) {
      mindist = dist;
      minid = jclus;
      mintime = dtime;
      mintimechod = dtimechod;
      mindx = dx;
      mindy = dy;
    }
  } 
  if (minid==-1) return minid;

  TRecoLKrCandidate *cluster = (TRecoLKrCandidate *)fLKrEvent->GetCandidate(minid);
  if (fFlag==0) {
    fUserMethods->FillHisto("LKrAnalysis_G_mindistvst_cluster",mintime,mindist);
    fUserMethods->FillHisto("LKrAnalysis_G_mindistvstchod_cluster",mintimechod,mindist);
    fUserMethods->FillHisto("LKrAnalysis_G_minxy_cluster",mindx,mindy);
    fUserMethods->FillHisto("LKrAnalysis_G_energyvsdist_cluster",mindist,cluster->GetClusterEnergy()/1000.);
  }
  if (fFlag==2) {
    fUserMethods->FillHisto("LKrAnalysis_G_2gammamindistvst_cluster",mintime,mindist);
    fUserMethods->FillHisto("LKrAnalysis_G_2gammamindistvstchod_cluster",mintimechod,mindist);
    fUserMethods->FillHisto("LKrAnalysis_G_2gammaminxy_cluster",mindx,mindy);
    fUserMethods->FillHisto("LKrAnalysis_G_2gammaenergyvsdist_cluster",mindist,cluster->GetClusterEnergy()/1000.);
  }
  fLKrClusterCandidate->SetDiscriminant(mindist);
  fLKrClusterCandidate->SetDeltaTime(mintime);
  fLKrClusterCandidate->SetEnergy(cluster->GetClusterEnergy()/1000.);
  fLKrClusterCandidate->SetSeedEnergy(cluster->GetClusterSeedEnergy()/1000.);
  fLKrClusterCandidate->SetX(cluster->GetClusterX());
  fLKrClusterCandidate->SetY(cluster->GetClusterY());
  fLKrClusterCandidate->SetNCells(cluster->GetNCells());
  return minid; 
}
Int_t LKrAnalysis_G::TrackCellMatching(Bool_t mcflag, Int_t runid, Double_t chodTime, Double_t trackTime, TRecoLKrEvent* event, TVector3 *posatlkr) {
  fLKrEvent = event;
  Double_t xpart = posatlkr->X();
  Double_t ypart = posatlkr->Y();
  Int_t minid = -1;

  // Look for cells associated to the track, if no cluster candidate matches the track
  Double_t eclus = 0.;
  Double_t tclus = 0.;
  Double_t xclus = 0.;
  Double_t yclus = 0.;
  Int_t ncells = 0;
  Double_t emax = -999999.;
  TClonesArray& Hits = (*(fLKrEvent->GetHits()));
  for (Int_t jHit=0; jHit<fLKrEvent->GetNHits(); jHit++) {
    TRecoLKrHit *hit = (TRecoLKrHit*)Hits[jHit];
    Int_t ixcell = (Int_t)hit->GetXCellID();
    Int_t iycell = (Int_t)hit->GetYCellID();
    Double_t xcell = hit->GetPosition().X();
    Double_t ycell = hit->GetPosition().Y();
    Double_t ecell = hit->GetEnergy()/1000.;
    if (ecell<=-999999.) continue;
    Double_t t0cell = (ixcell>=0 && iycell>=0) ? fLKrCellT0[ixcell][iycell] : 0;
    Double_t tcell = hit->GetTime()+3.2-t0cell;
    Double_t dist = sqrt((xcell-xpart)*(xcell-xpart)+(ycell-ypart)*(ycell-ypart));
    Double_t dtime = tcell-trackTime;
    if (dist<150. && fabs(dtime)<40) { // Select a region around the extrapolated track position at LKr
      eclus += ecell;
      if (emax<ecell) {
        emax = ecell;
        tclus = tcell+1.3;
      }
      xclus += xcell*ecell;
      yclus += ycell*ecell;
      ncells++;
      minid = 1;
    }
  }
  if (!ncells || emax<0.04) { // no cells found
    tclus = -99999.;
    eclus = 0;
    xclus = -99999.;
    yclus = -99999.;
    minid = -1;
    return minid;
  }
  
  xclus /= eclus;
  yclus /= eclus;
  Double_t ue = eclus;
  if (ncells>9) {
    if (ue<22) eclus = ue/(0.7666+0.0573489*log(ue));
    if (ue>=22 && ue<65) eclus = ue/(0.828962+0.0369797*log(ue));
    if (ue>=65) eclus = ue/(0.828962+0.0369797*log(65));
  } 
  eclus *= 1.03;
  eclus = eclus<15 ? (eclus+0.015)/(15+0.015)*15*0.9999 : eclus; 
  Double_t mindist = sqrt((xclus-xpart)*(xclus-xpart)+(yclus-ypart)*(yclus-ypart));
  fLKrCellCandidate->SetDiscriminant(mindist);
  fLKrCellCandidate->SetDeltaTime(tclus-trackTime);
  fLKrCellCandidate->SetEnergy(eclus);
  fLKrCellCandidate->SetSeedEnergy(emax);
  fLKrCellCandidate->SetX(xclus);
  fLKrCellCandidate->SetY(yclus);
  fLKrCellCandidate->SetNCells(ncells);
  if (fFlag==0) {
    fUserMethods->FillHisto("LKrAnalysis_G_mindistvst_cell",tclus-trackTime,mindist);
    fUserMethods->FillHisto("LKrAnalysis_G_mindistvstchod_cell",tclus-chodTime,mindist);
    fUserMethods->FillHisto("LKrAnalysis_G_minxy_cell",xclus-xpart,yclus-ypart);
    fUserMethods->FillHisto("LKrAnalysis_G_energyvsdist_cell",mindist,eclus);
  }
  if (fFlag==2) {
    fUserMethods->FillHisto("LKrAnalysis_G_2gammamindistvst_cell",tclus-trackTime,mindist);
    fUserMethods->FillHisto("LKrAnalysis_G_2gammamindistvstchod_cell",tclus-chodTime,mindist);
    fUserMethods->FillHisto("LKrAnalysis_G_2gammaminxy_cell",xclus-xpart,yclus-ypart);
    fUserMethods->FillHisto("LKrAnalysis_G_2gammaenergyvsdist_cell",mindist,eclus);
  }
  return minid;
}

///////////////////////////////////////
// Photon search for pinunu analysis //
///////////////////////////////////////
Int_t LKrAnalysis_G::ExtraClusters(Double_t reftime, TRecoLKrEvent* event, TVector2 *posmatchclust, TVector3 *posatlkr, Int_t eventnumber) {
  Int_t nPhoton = 0;
  TClonesArray& Hits = *(event->GetHits());
  for (Int_t jclus=0; jclus<event->GetNCandidates(); jclus++) {
    TRecoLKrCandidate *thisCluster = (TRecoLKrCandidate *)event->GetCandidate(jclus);
    Int_t ixcell = -1;
    Int_t iycell = -1;
    Int_t seedid = thisCluster->GetIdSeed();
    if (seedid>=0) { // LKr T0 fine 
      TRecoLKrHit *hit = (TRecoLKrHit*)Hits[seedid];
      ixcell = hit->GetXCellID();
      iycell = hit->GetXCellID();
    }
    Double_t t0cell = (ixcell>=0 && iycell>=0) ? fLKrCellT0[ixcell][iycell] : 0;
    TVector2 posclust(thisCluster->GetClusterX(),thisCluster->GetClusterY()); 
    TVector2 postrack(posatlkr->X(),posatlkr->Y());
    Double_t dist = (posclust-postrack).Mod();
    Double_t dtime = thisCluster->GetClusterTime()-reftime-t0cell;
    fUserMethods->FillHisto("LKrAnalysis_G_not_matched_all",dist,dtime);
////    if (eventnumber==48381) cout << "       " << jclus << " " << reftime << " " << dist << " " << dtime << " " << (posclust-*posmatchclust).Mod() << endl;
    if ((posclust-*posmatchclust).Mod()<5.) continue;
    fUserMethods->FillHisto("LKrAnalysis_G_not_matched",dist,dtime);
    Bool_t isPhoton = dist>150. ? SelectTiming(dtime,thisCluster->GetClusterEnergy()/1000.) : false;
    if (!isPhoton) continue;
    fUserMethods->FillHisto("LKrAnalysis_G_not_matched_photon",dist,dtime);
    if (nPhoton<10) {
      fPhotonCandidate[nPhoton]->SetTime(dtime+reftime);
      fPhotonCandidate[nPhoton]->SetX(posclust.X());
      fPhotonCandidate[nPhoton]->SetY(posclust.Y());
      fPhotonCandidate[nPhoton]->SetEnergy(thisCluster->GetClusterEnergy()/1000.);
      fPhotonCandidate[nPhoton]->SetNCells(thisCluster->GetNCells());
    }
    nPhoton++;
  }   
  if (!nPhoton) return 0;
  
  return nPhoton;
}
Bool_t LKrAnalysis_G::SelectTiming(Double_t dtime, Double_t energy) {
  Double_t dtime12 = dtime-12.5;
  Double_t dtime25 = dtime+25;
  if (fYear==2016) dtime = dtime-1.22; // Offset measured during 2016 run. To be removed when centrally fixed 
  if (energy>2) {
    fUserMethods->FillHisto("LKrAnalysis_G_not_matched_highe",dtime);
    if (dtime>-10 && dtime<15) return true; // Broad cut because of the photon rate: adjust ? (see analysis meeting 22/07/2016)
    if (fabs(dtime25)<3) return true;
    else return false;
  } else {
    fUserMethods->FillHisto("LKrAnalysis_G_not_matched_lowe",dtime);
    if (dtime>-5 && dtime<8) return true; // Strict cut because of the MIP rate: adjust ? (see analysis meeting 22/07/2016)
    else return false;
  }
  return false;
}

//////////////////////////////////////
//////////////////////////////////////
////                              ////
////  ALTERNATIVE RECONSTRUCTION  ////
////                              ////
//////////////////////////////////////
//////////////////////////////////////
Int_t LKrAnalysis_G::FindNewClusters(Double_t timeref, TRecoLKrEvent *event, TVector2 *posmatchclust, TVector3 *posatlkr)
{
  // Init
  Int_t newCluster = 0;
  ReadCell(timeref,event,posmatchclust,posatlkr);  
  OrderRawColumn();

  // Variables definition
  Int_t iread[14000],id[14000];
  for(Int_t j=0;j<Ncell;j++) 
  {
    iread[j]=0;
    id[j]=-1;
  }

  // Cascade loop over the cells
  Int_t m=0;
  for(Int_t i=0;i<Ncell;i++)
  {
    if (iread[i]) continue;
    Int_t l=0;
    Int_t k=i;
    Double_t energy = 0.;
    Double_t xpos = 0.;
    Double_t ypos = 0.;
    Double_t time = 0.;
    Double_t energyprev = 0.;
    for(Int_t j=0;j<Ncell;j++)
    {
      if(Nearby(i,j)) 
      {
        if (!iread[j])
        {
          iread[j]=1;
          id[l]=j;
          i=id[l];
          energy+=Ecells[i];
          if (Ecells[i]>energyprev) {
            time = Tcells[i];
            energyprev = Ecells[i];
          }
          xpos+=Xcells[i]*Ecells[i];
          ypos+=Ycells[i]*Ecells[i];
          j=0;
          l++;
        }
      }
      if(j==Ncell-1&&m<l)
      {
        j=0;
        i = id[m];
        m++;
      } 
    }
    
    // Store cluster
    if (l>1 && newCluster<5)
    {
      if (energy>0) 
      {
        xpos=xpos/energy;
        ypos=ypos/energy;
      } else {
        energy=0;
        time=999999.;
        xpos=999999.;
        ypos=999999.;
      }

      TVector2 postrack(posatlkr->X(),posatlkr->Y());
      TVector2 posclust(xpos,ypos); 
      if ((posclust-postrack).Mod()>=150. && (posclust-*posmatchclust).Mod()>100. && energy>0) {
        fNewClusterCandidate[newCluster]->SetEnergy(energy);
        fNewClusterCandidate[newCluster]->SetTime(time);
        fNewClusterCandidate[newCluster]->SetX(xpos);
        fNewClusterCandidate[newCluster]->SetY(ypos);
        fNewClusterCandidate[newCluster]->SetNCells(l);
        newCluster++;
      }
    }

    // Reset
    i=k;
    m=0;
  }

  return newCluster;
}

void LKrAnalysis_G::ReadCell(Double_t timeref, TRecoLKrEvent *event, TVector2 *posmatchclust, TVector3 *posatlkr)
{
  fill(xcellcoor,xcellcoor+14000,-1);
  fill(ycellcoor,ycellcoor+14000,-1);
  fill(Tcells,Tcells+14000,-99999);
  fill(Ecells,Ecells+14000,0);
  fill(Xcells,Xcells+14000,-99999);
  fill(Ycells,Ycells+14000,-99999);
  TClonesArray& Hits = *(event->GetHits());
  Int_t m = 0;
//cout<<"NHits"<<event->GetNHits()<<endl; 
  for (Int_t j=0; j<event->GetNHits(); j++) { 
    TRecoLKrHit *hit = (TRecoLKrHit*)Hits[j];
    if (hit->GetEnergy()/1000.<=0.04) continue;//took out /1000
//	cout<<"enough energy"<<endl;
    Double_t xcell = hit->GetPosition().X();
    Double_t ycell = hit->GetPosition().Y();
    TVector2 poscell(xcell,ycell); 
//    if ((poscell-*posmatchclust).Mod()<100.) continue; // far from the cluster associated with the track extrapolation
    TVector2 postrack(posatlkr->X(),posatlkr->Y());
//    if ((poscell-postrack).Mod()<150.) continue; // far from the track extrapolation
    Int_t ixcell = (Int_t)hit->GetXCellID();
    Int_t iycell = (Int_t)hit->GetYCellID();
    Double_t t0cell = (ixcell>=0 && iycell>=0) ? fLKrCellT0[ixcell][iycell] : 0;
    Double_t tcell = hit->GetTime()+3.2-t0cell;
	
    if (fabs(tcell-timeref)<40) continue; // in time with the reference, changed from > to <
// cout<<"in time"<<endl;
    if (fabs(hit->GetEnergy()/1000.)>100) continue;//took out /1000
//cout<<"too much energy?"<<endl;
    xcellcoor[m] = ixcell;
    ycellcoor[m] = iycell;
    Tcells[m] = hit->GetTime()+3.2;
    Ecells[m] = hit->GetEnergy()/1000.;//took out /1000
    
Xcells[m] = xcell;
    Ycells[m] = ycell;
    m++; 
  }
  Ncell = m;
}

void LKrAnalysis_G::OrderRawColumn()
{
  Int_t i,ic,irmin[130],icmin[130],l,k,j,step,skip,kk,skipr,written,foundmin,Nn;
  Double_t rmin=9999,cmin=9999;
  Int_t xorder[14000],yorder[14000];
  Double_t torder[14000],eorder[14000],xposorder[14000],yposorder[14000];
  Int_t xord[14000],yord[14000];
  Double_t tord[14000],eord[14000],xposord[14000],yposord[14000];

  Nn=0;
  step=0;
  for (kk=0;kk<Ncell;kk++)
  {
    foundmin=0;
    written=0;
    rmin=9999;
    for (i=0;i<Ncell;i++)
    {
      skipr=0;
      for (j=0;j<step;j++){if (ycellcoor[i]==ycellcoor[irmin[j]]) skipr=1;}
      if (skipr) continue;
      if (ycellcoor[i]<rmin)
      {
        foundmin=1;
        rmin=ycellcoor[i];
        irmin[step]=i; 
      }
    }
    if (foundmin==0) continue;
    l=0;
    for (i=0;i<Ncell;i++)
    {
      if (ycellcoor[i]==ycellcoor[irmin[step]])
      {
        xorder[l]=xcellcoor[i];
        yorder[l]=ycellcoor[i];
        torder[l]=Tcells[i];
        eorder[l]=Ecells[i];
        xposorder[l]=Xcells[i];
        yposorder[l]=Ycells[i];
        l++;
      }
    }
    for (k=0;k<l;k++)
    {
      cmin=9999;
      for (ic=0;ic<l;ic++)
      {
        skip=0;
        for (j=0;j<k;j++){if (ic==icmin[j]) skip=1;}
        if (skip) continue;
        if (xorder[ic]<cmin)
        {
          cmin=xorder[ic];
          icmin[k]=ic;
        }
      }
      xord[Nn+k]=xorder[icmin[k]];
      yord[Nn+k]=yorder[icmin[k]];
      tord[Nn+k]=torder[icmin[k]];
      eord[Nn+k]=eorder[icmin[k]];
      xposord[Nn+k]=xposorder[icmin[k]];
      yposord[Nn+k]=yposorder[icmin[k]];
      written=1;
    }
    Nn+=k;
    if (written) step++;
  }

  for (i=0;i<Ncell;i++)
  {
    xcellcoor[i]=xord[i];
    ycellcoor[i]=yord[i];
    Ecells[i]=eord[i];
    Tcells[i]=tord[i];
    Xcells[i]=xposord[i];
    Ycells[i]=yposord[i];
  }

}

Int_t LKrAnalysis_G::Nearby(Int_t i, Int_t j)
{
  for (Int_t kk=0; kk<10; kk++)
  {
    if (xcellcoor[j]==xcellcoor[i]+kk)
    {
      for (Int_t k=0; k<10; k++)
      {
        if ((ycellcoor[j]==ycellcoor[i]+k) ||
            (ycellcoor[j]==ycellcoor[i]-k)) return 1;
      }
    }
    if (xcellcoor[j]==xcellcoor[i]-kk)
    {
      for (Int_t k=0; k<10; k++)
      {
        if ((ycellcoor[j]==ycellcoor[i]+k) ||
            (ycellcoor[j]==ycellcoor[i]-k)) return 1;
      }
    }
  }
  return 0;
}
 
