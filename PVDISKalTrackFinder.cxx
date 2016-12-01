//c++
#include <cmath>

//SoLIDTracking
#include "PVDISKalTrackFinder.h"
#include "SoLKalTrackSystem.h"
#include "SoLKalTrackSite.h"
#include "SoLKalTrackState.h"
#define MAXHITGEM 1500
#define MAXSEED 10000
PVDISKalTrackFinder::PVDISKalTrackFinder(bool isMC)
:SoLKalTrackFinder(), fIsMC(isMC)
{
  fChi2PerNDFCut = 30;
  fTargetPlaneZ  =  0.1;
  fTargetCenter  =  0.1;
  fTargetLength  =  0.4;
  fGEMTracker.clear();
  fWindowHits.clear();
  fWindowHits.reserve(MAXWINDOWHIT);
}
//_____________________________________________________________________________
PVDISKalTrackFinder::~PVDISKalTrackFinder()
{
  Clear();
  delete fCoarseTracks;
}
//_____________________________________________________________________________
void PVDISKalTrackFinder::Clear( Option_t* opt )
{

  if (fCoarseTracks->GetEntries() != 0)
  fCoarseTracks->Delete();

  fCoarseTracks->Clear(opt);

  fCaloHits = nullptr;
  fNSeeds = 0;
  fNGoodTrack = 0;

  fSeedEfficiency = false;
  fMcTrackEfficiency = false;

  map< SeedType, vector<DoubletSeed> >::iterator itt;
  for (itt = fSeedPool.begin(); itt != fSeedPool.end(); itt++) { (itt->second).clear(); }

  map< Int_t, vector<SoLIDGEMHit*> >::iterator it;
  for (it = fGoodHits.begin(); it != fGoodHits.end(); it++) { (it->second).clear(); }
  fGoodHits.clear();
}
//____________________________________________________________________________
void PVDISKalTrackFinder::SetGEMDetector(vector<SoLIDGEMTracker*> thetrackers)
{
  fGEMTracker = thetrackers;
  fNTrackers  = (Int_t)thetrackers.size();
  assert(thetrackers.size() > 2 && thetrackers[2]->GetNChamber() == 1);
  fRefPhi = fGEMTracker[2]->GetChamber(0)->GetPhiInLab();
  fRefSin = sin(-1.* fRefPhi);
  fRefCos = cos(-1.* fRefPhi);
}
//____________________________________________________________________________
void PVDISKalTrackFinder::ProcessHits(TClonesArray* theTracks)
{
  assert(fGEMTracker.size() != 0);
  assert(fCaloHits == nullptr);

  fRefPhi = fGEMTracker[2]->GetChamber(0)->GetPhiInLab();
  fCaloHits = fECal->GetCaloHits();

  //finding doublet seed from last three GEM planes
  FindDoubletSeed(3, 4);
  FindDoubletSeed(2, 4);
  FindDoubletSeed(2, 3);

  //merge doublet seed to from triplets
  MergeSeed();

  //Follow the direction of seed and look for potential hits
  TrackFollow();

  //find the interaction vertex
  FindandAddVertex();


  FinalSelection(theTracks);
  fEventNum++;
}
//______________________________________________________________________________
void PVDISKalTrackFinder::FindDoubletSeed(Int_t planej, Int_t planek)
{
  //not using the front trackers to make seed
  assert(planek > planej && planek >= 1 && planej >=2);

  SeedType seedType = kMidBack;
  if (planej == 3 && planek == 4){
    seedType = kMidBack;
  }else if(planej == 2 && planek == 3){
    seedType = kFrontMid;
  }else if (planej == 2 && planek == 4){
    seedType = kFrontBack;
  }

  for (int k=0; k<fGEMTracker[planek]->GetNChamber(); k++){
    TSeqCollection* planekHitArray = fGEMTracker[planek]->GetChamber(k)->GetHits();

    int totalHitk = planekHitArray->GetLast()+1;
    if (totalHitk > MAXHITGEM) return;

    for (int nhitk = 0; nhitk < totalHitk; nhitk++){
      SoLIDGEMHit *hitk = (SoLIDGEMHit*)planekHitArray->At(nhitk);

      int ECIndexk = 0;
      if (planek >= 3 && !ECCoarseCheck(hitk, ECIndexk)) continue;
      assert(ECIndexk >= 0);

      for (int j=0; j<fGEMTracker[planej]->GetNChamber(); j++){
        TSeqCollection* planejHitArray = fGEMTracker[planej]->GetChamber(j)->GetHits();

        int totalHitj = planejHitArray->GetLast()+1;
        if (totalHitj > MAXHITGEM) return;

        for (int nhitj = 0; nhitj < totalHitj; nhitj++){
          SoLIDGEMHit *hitj = (SoLIDGEMHit*)planejHitArray->At(nhitj);

          //TODO: What if there are two very close EC hits, the two GEM hits may match to different EC hits
          if (planej >= 3 && !ECCoarseCheck(hitj, ECIndexk)) continue;
          assert(ECIndexk >= 0);

          //after coarse check with EC, we use straight line to connect to the GEM hits and see if it lead to the EC hit
          Double_t xk  = hitk->GetX();
          Double_t yk  = hitk->GetY();
          Double_t zk  = hitk->GetZ();
          Double_t xj  = hitj->GetX();
          Double_t yj  = hitj->GetY();
          Double_t zj  = hitj->GetZ();
          Double_t xec = fCaloHits->at(ECIndexk).fXPos;
          Double_t yec = fCaloHits->at(ECIndexk).fYPos;

          Double_t initMom = fCaloHits->at(ECIndexk).fEdp;
          Double_t initTheta = atan( (sqrt(xk*xk + yk*yk) - sqrt(xj*xj + yj*yj))/(zk-zj) );
          Double_t initPhi = atan2(yk - yj, xk - xj);
          Double_t charge = -1;
          ECType type = kFAEC;

          //senity check for the local theta angle
          if (initTheta > 0.7 || initTheta < 0.3) continue;

          Rotate(xk, yk);
          Rotate(xj, yj);
          Rotate(xec, yec);
          Double_t predictX = StraightLinePredict(xk, zk, xj, zj, fECal->GetECZ(kFAEC));
          Double_t predictY = StraightLinePredict(yk, zk, yj, zj, fECal->GetECZ(kFAEC));

          if (sqrt(pow(predictX - xec, 2) + pow(predictY - yec, 2)) > 0.05) continue;

          //so the hit pairs has passed all the cuts, now we can save it into a container and waiting for merge
          map< SeedType, vector<DoubletSeed> >::iterator it = fSeedPool.find(seedType);

          if (it != fSeedPool.end()){
            (it->second).push_back(DoubletSeed(seedType, hitj, hitk, initMom, initTheta, initPhi, charge, type));
          }
          else{
            cout<<"should never happen, fSeedPool should be init in constructor"<<endl;
            vector<DoubletSeed> thisVector;
            thisVector.push_back(DoubletSeed(seedType, hitj, hitk, initMom, initTheta, initPhi, charge, type));
            fSeedPool.insert(std::pair< SeedType, vector<DoubletSeed> >(seedType, thisVector));
          }

#ifdef MCDATA
          if (dynamic_cast<SoLIDMCGEMHit*>(hitj)->IsSignalHit() == 1 && dynamic_cast<SoLIDMCGEMHit*>(hitk)->IsSignalHit() == 1)
          fSeedEfficiency = true;
#endif
        }
      }
    }
  }

}
//______________________________________________________________________________
void PVDISKalTrackFinder::MergeSeed()
{
  unsigned int totalSeed = fSeedPool[kMidBack].size() + fSeedPool[kFrontMid].size()
                         + fSeedPool[kFrontBack].size();
  if (totalSeed > MAXSEED) return;

  for (unsigned int i=0; i<fSeedPool[kMidBack].size(); i++){
    for (unsigned int j=0; j<fSeedPool[kFrontMid].size(); j++){
      //see if there is a common point in between the two doublet seeds
      if (fSeedPool[kMidBack].at(i).hita == fSeedPool[kFrontMid].at(j).hitb){
        //now check the third type of seed
        for (unsigned int k=0; k<fSeedPool[kFrontBack].size(); k++){
          if (fSeedPool[kFrontBack].at(k).hita == fSeedPool[kFrontMid].at(j).hita &&
            fSeedPool[kFrontBack].at(k).hitb == fSeedPool[kMidBack].at(i).hitb){

            Double_t xa = (fSeedPool[kFrontMid].at(j).hita)->GetX();
            Double_t ya = (fSeedPool[kFrontMid].at(j).hita)->GetY();
            Double_t za = (fSeedPool[kFrontMid].at(j).hita)->GetZ();

            Double_t xb = (fSeedPool[kFrontMid].at(j).hitb)->GetX();
            Double_t yb = (fSeedPool[kFrontMid].at(j).hitb)->GetY();
            Double_t zb = (fSeedPool[kFrontMid].at(j).hitb)->GetZ();

            Double_t xc = (fSeedPool[kMidBack].at(i).hitb)->GetX();
            Double_t yc = (fSeedPool[kMidBack].at(i).hitb)->GetY();
            Double_t zc = (fSeedPool[kMidBack].at(i).hitb)->GetZ();

            Rotate(xa, ya);
            Rotate(xb, yb);
            Rotate(xc, yc);

            Double_t predictX = StraightLinePredict(xa, za, xb, zb, zc);
            Double_t predictY = StraightLinePredict(ya, za, yb, zb, zc);

            if (fabs(predictX - xc) > 0.01) continue;
            if (fabs(predictY - yc) > 0.006) continue;

            fSeedPool[kMidBack].at(i).Deactive();
            fSeedPool[kFrontMid].at(j).Deactive();
            fSeedPool[kFrontBack].at(k).Deactive();

            SoLKalTrackSite & initSite =  SiteInitWithSeed(&(fSeedPool[kFrontBack].at(k)));
            SoLKalTrackSystem *thisSystem = new ((*fCoarseTracks)[fNSeeds++]) SoLKalTrackSystem();
            thisSystem->SetMass(kElectronMass);
            thisSystem->SetCharge(fSeedPool[kFrontBack].at(k).charge);
            thisSystem->SetElectron(kTRUE);
            thisSystem->SetAngleFlag(fSeedPool[kFrontBack].at(k).flag);
            thisSystem->SetSeedType(kTriplet);
            thisSystem->SetOwner();
            thisSystem->Add(&initSite);

            //remember finding tracks always go backward
            SoLKalTrackSite& backSite = *new SoLKalTrackSite(fSeedPool[kMidBack].at(i).hitb, kMdim, kSdim, kMdim*fChi2PerNDFCut);
            if (!(thisSystem->AddAndFilter(backSite))) thisSystem->SetTrackStatus(false);

            SoLKalTrackSite& midSite = *new SoLKalTrackSite(fSeedPool[kMidBack].at(i).hita, kMdim, kSdim, kMdim*fChi2PerNDFCut);
            if (!(thisSystem->AddAndFilter(midSite))) thisSystem->SetTrackStatus(false);

            SoLKalTrackSite& frontSite = *new SoLKalTrackSite(fSeedPool[kFrontBack].at(k).hita, kMdim, kSdim, kMdim*fChi2PerNDFCut);
            if (!(thisSystem->AddAndFilter(frontSite))) thisSystem->SetTrackStatus(false);
          }
        }
      }
    }
  }

  //end of triplet seed matching and begin the remaining doublet seed init
  map< SeedType, vector<DoubletSeed> >::iterator it;
  for (it = fSeedPool.begin(); it != fSeedPool.end(); it++){
    vector<DoubletSeed> & thisVector = (it->second);

    for (unsigned int i=0; i<thisVector.size(); i++){
      if (!thisVector.at(i).isActive) continue;
      SoLKalTrackSite & initSite =  SiteInitWithSeed(&(thisVector.at(i)));
      SoLKalTrackSystem *thisSystem = new ((*fCoarseTracks)[fNSeeds++]) SoLKalTrackSystem();
      thisSystem->SetMass(kElectronMass);
      thisSystem->SetCharge(thisVector.at(i).charge);
      thisSystem->SetElectron(kTRUE);
      thisSystem->SetAngleFlag(thisVector.at(i).flag);
      thisSystem->SetSeedType(thisVector.at(i).type);
      thisSystem->SetOwner();
      thisSystem->Add(&initSite);
      //We assume that the doublet seed has already missed one hit (otherwise it is suppose to be part
      //of a triplet seed and thus be set as inactived already)
      thisSystem->AddMissingHits();

      SoLKalTrackSite& backSite = *new SoLKalTrackSite(thisVector.at(i).hitb, kMdim, kSdim, kMdim*fChi2PerNDFCut);
      if (!(thisSystem->AddAndFilter(backSite))) thisSystem->SetTrackStatus(false);

      SoLKalTrackSite& midSite = *new SoLKalTrackSite(thisVector.at(i).hita, kMdim, kSdim,  kMdim*fChi2PerNDFCut);
      if (!(thisSystem->AddAndFilter(midSite))) thisSystem->SetTrackStatus(false);
    }
  }
}
//______________________________________________________________________________
void PVDISKalTrackFinder::TrackFollow()
{
    //this function is responsible for propagating the seed track toward the next tracker, find suitable hits
  //the process stop until the track reach the first tracker upstream (track searching always go backward)

  for (Int_t i=0; i<fCoarseTracks->GetLast()+1; i++){
    SoLKalTrackSystem* thisSystem = (SoLKalTrackSystem*)(fCoarseTracks->At(i));
    thisSystem->CheckTrackStatus();
    if (!thisSystem->GetTrackStatus()) continue; //skip the bad tracks

    thisSystem->SetCurInstancePtr(thisSystem);
    Int_t currentTracker = ((thisSystem->GetCurSite()).GetHit())->GetTrackerID();

    //seed from type kMidBack will skip the front seed plane. We assume for this type of seed, the hit on the
    //front seed plane is missing, (otherwise the seed should be absorbed into the triplet seed)
    if (thisSystem->GetSeedType() == kMidBack) currentTracker--;

    while (currentTracker > 0){
      currentTracker--;

      thisSystem->CheckTrackStatus();
      if (!thisSystem->GetTrackStatus()) break; //skip the bad tracks

      SoLKalTrackState currentState = (thisSystem->GetCurSite()).GetCurState();
      currentState.InitPredictSV();
      SoLKalTrackState *predictState = currentState.PredictSVatNextZ(fGEMTracker[currentTracker]->GetZ());

      bool flag = (thisSystem->GetNHits() >= 2);

      int size = GetHitsInWindow(currentTracker, (*predictState)(kIdxX0, 0), (predictState->GetCovMat())(kIdxX0, kIdxX0),
                                (*predictState)(kIdxY0, 0), (predictState->GetCovMat())(kIdxY0, kIdxY0), flag);


      if (size <= 0){

        thisSystem->AddMissingHits();

      }
      else if (size == 1){
        SoLKalTrackSite &newSite = *new SoLKalTrackSite(fWindowHits.at(0), kMdim, kSdim, kMdim*fChi2PerNDFCut);
        newSite.Add(predictState);
        if (newSite.Filter()){
          thisSystem->Add(&newSite);
          thisSystem->IncreaseChi2(newSite.GetDeltaChi2());
          currentState.ClearAttemptSV();
        }
        else{
          thisSystem->AddMissingHits();
        }

      }
      else{
        //find the cloest one for now, should use concurrent tracking in the future
        SoLKalTrackSite &newSite = *new SoLKalTrackSite(FindCloestHitInWindow((*predictState)(kIdxX0, 0),
                                    (*predictState)(kIdxY0, 0)), kMdim, kSdim,  kMdim*fChi2PerNDFCut);
        newSite.Add(predictState);
        if (newSite.Filter()){
          thisSystem->Add(&newSite);
          thisSystem->IncreaseChi2(newSite.GetDeltaChi2());
          currentState.ClearAttemptSV();
        }
        else{
          thisSystem->AddMissingHits();
        }

      }
    }

    //now that we have all the hits selected, we can look at the chi2 per ndf and charge asymmetry to
    //get rid of some potential bad tracks, before doing other things
    if (thisSystem->GetChi2perNDF() > fChi2PerNDFCut) {
      thisSystem->SetTrackStatus(kFALSE);
      continue;
    }
    if (!CheckChargeAsy(thisSystem)){
      thisSystem->SetTrackStatus(kFALSE);
      continue;
    }
#ifdef MCDATA
    //------------check MC track efficiency--------------//
    if (!thisSystem->GetTrackStatus()) continue;
    bool allMC = true;
    for (Int_t j=1; j!=thisSystem->GetLast()+1;j++){
      SoLIDGEMHit* thisHit = (SoLIDGEMHit*)((SoLKalTrackSite*)thisSystem->At(j))->GetHit();
      if (dynamic_cast<SoLIDMCGEMHit*>(thisHit)->IsSignalHit() != 1) allMC = false;
    }
    if (allMC) fMcTrackEfficiency = true;
    //---------------------------------------------------//
#endif
  }
}
//______________________________________________________________________________
void PVDISKalTrackFinder::FindandAddVertex()
{
   for (Int_t i=0; i<fCoarseTracks->GetLast()+1; i++){
      SoLKalTrackSystem* thisSystem = (SoLKalTrackSystem*)(fCoarseTracks->At(i));
      thisSystem->SetCurInstancePtr(thisSystem);

      thisSystem->CheckTrackStatus();
      if (thisSystem->GetTrackStatus() == kFALSE) continue; //skip bad tracks

      SoLKalTrackState currentState = (thisSystem->GetCurSite()).GetCurState();
      currentState.InitPredictSV();

      SoLKalTrackState *predictState = currentState.PredictSVatNextZ(fTargetCenter);
      Double_t vertexz = FindVertexZ(predictState);

      if (fabs(vertexz - fTargetCenter) > 0.25){
        thisSystem->SetTrackStatus(kFALSE);
        continue;
      }

      //propagate the state vector to the interaction vertex that just found
      //not sure if this is the best way to add vertex
      currentState.InitPredictSV();
      predictState = currentState.PredictSVatNextZ(vertexz);

      //make a site at the interaction vertex to add to the fitting
      SoLKalTrackSite &vertexSite = *new SoLKalTrackSite(kMdim, kSdim,  kGiga);
      vertexSite.SetMeasurement(fBPMX, fBPMY);
      vertexSite.SetHitResolution(3e-4, 3e-4);
      vertexSite.Add(predictState);
      if (vertexSite.Filter()){
        //calculate vertex variables and set info to the track system
        Double_t temp_tx =  vertexSite.GetCurState()(kIdxTX, 0);
        Double_t temp_ty =  vertexSite.GetCurState()(kIdxTY, 0);
        Double_t temp_qp =  vertexSite.GetCurState()(kIdxQP, 0);

        TVector3 vertex_vdir;
        vertex_vdir.SetZ( 1./(TMath::Sqrt(temp_tx*temp_tx + temp_ty*temp_ty + 1. )) );
        vertex_vdir.SetX(temp_tx * vertex_vdir.Z());
        vertex_vdir.SetY(temp_ty * vertex_vdir.Z());
        vertex_vdir = vertex_vdir.Unit();

        thisSystem->SetMomentum(thisSystem->GetCharge()/temp_qp);
        thisSystem->SetTheta(acos(1./sqrt(1. + pow( (vertex_vdir.X()/vertex_vdir.Z()), 2)
                             + pow((vertex_vdir.Y()/vertex_vdir.Z()) , 2))));

        thisSystem->SetVertexZ(vertexz);
        thisSystem->SetPhi(atan2(vertex_vdir.Y(), vertex_vdir.X()));
      }
      currentState.ClearAttemptSV();
      delete &vertexSite;
   }
}
//______________________________________________________________________________
void PVDISKalTrackFinder::FinalSelection(TClonesArray *theTracks)
{
  fCoarseTracks->Sort();

  /*Int_t countTrack = 0;

  for (Int_t i=0; i<fCoarseTracks->GetLast()+1; i++){
    SoLKalTrackSystem *thisSystem = (SoLKalTrackSystem*)(fCoarseTracks->At(i));
    if (thisSystem->GetTrackStatus()) countTrack++;
  }
  cout<<countTrack<<endl;*/

  for (Int_t i=0; i<fCoarseTracks->GetLast()+1; i++){

    SoLKalTrackSystem *thisSystem = (SoLKalTrackSystem*)(fCoarseTracks->At(i));
    thisSystem->SetCurInstancePtr(thisSystem);

    if (!thisSystem->GetTrackStatus()) continue;
    Int_t flag = 0;

    //start from 1 because the 0th is the dummy site that we used to initialize Kalman Filter
    //TODO remember not to add the last one since later it will be the BPM, not GEM hit
    for (Int_t j=1;j!=thisSystem->GetLast()+1;j++){
      SoLIDGEMHit *thishit = (SoLIDGEMHit*)(static_cast<SoLKalTrackSite*>(thisSystem->At(j))->GetHit());
      Int_t layer = thishit->GetTrackerID();

      map< Int_t, vector<SoLIDGEMHit*> >::iterator it = fGoodHits.find(layer);

      if (it != fGoodHits.end()){
        for (UInt_t n = 0; n<(it->second).size(); n++){
          if ((thishit->GetX() == ((it->second).at(n))->GetX()) &&
              (thishit->GetY() == ((it->second).at(n))->GetY())) { flag = 1; }
        }
      }

    }
    if (flag == 0){
      SoLIDTrack* newtrack = 0;
      if (fIsMC){
#ifdef MCDATA
        newtrack = new ((*theTracks)[fNGoodTrack++]) SoLIDMCTrack();
#endif
      }
      else{
        newtrack = new ((*theTracks)[fNGoodTrack++]) SoLIDTrack();
      }
    CopyTrack(newtrack, thisSystem);
    }

  }

}
//______________________________________________________________________________
inline Bool_t PVDISKalTrackFinder::ECCoarseCheck(SoLIDGEMHit *theHit, Int_t& index)
{
  for (UInt_t ec_count=0; ec_count<fCaloHits->size(); ec_count++){
	  assert(fCaloHits->at(ec_count).fECID != kLAEC); //should never happen for PVDIS
	  Double_t ecHitPhi = TMath::ATan2(fCaloHits->at(ec_count).fYPos, fCaloHits->at(ec_count).fXPos);
	  Double_t ecHitR = TMath::Sqrt( TMath::Power(fCaloHits->at(ec_count).fXPos, 2) +
	    			  	TMath::Power(fCaloHits->at(ec_count).fYPos, 2) );
	  Double_t tmpDeltaPhi = CalDeltaPhi(ecHitPhi, theHit->GetPhi());
	  Double_t tmpDeltaR   = CalDeltaR(ecHitR, theHit->GetR());
	  if (theHit->GetTrackerID()==3 && (tmpDeltaPhi < 0.035 && tmpDeltaPhi > -0.025) && (tmpDeltaR < 0.18 && tmpDeltaR > 0.02)){
	   index = ec_count;
       return kTRUE;
     }else if (theHit->GetTrackerID()==4 && (tmpDeltaPhi<0.03 && tmpDeltaPhi>-0.03) && (tmpDeltaR<0.11 && tmpDeltaR > -0.01) ){
       index = ec_count;
       return kTRUE;
     }
  }
  return kFALSE;
}
//_______________________________________________________________________________
inline SoLKalTrackSite & PVDISKalTrackFinder::SiteInitWithSeed(DoubletSeed* thisSeed)
{
  TVector3 initDir(cos(thisSeed->initPhi), sin(thisSeed->initPhi),
                   1./tan(thisSeed->initTheta));
  initDir = initDir.Unit();

  //-----------prepare seeds for Kalman Filter track finding------------//
  SoLKalMatrix svd(kSdim,1);
  svd.Zero();
  svd(kIdxX0,0) = (thisSeed->hitb)->GetX();
  svd(kIdxY0,0) = (thisSeed->hitb)->GetY();
  svd(kIdxTX,0) = initDir.X()/initDir.Z();
  svd(kIdxTY,0) = initDir.Y()/initDir.Z();
  svd(kIdxQP,0) = thisSeed->charge/thisSeed->initMom;

  SoLKalMatrix C(kSdim,kSdim);
  C.Zero();

  Double_t phi = (thisSeed->hitb)->GetPhi();
  Double_t dr = 4.e-4;
  Double_t drphi = 4.5e-5;

  Double_t dx = sqrt( pow( cos(phi)*dr, 2) + pow( sin(phi)*drphi, 2) );
  Double_t dy = sqrt( pow( sin(phi)*dr, 2) + pow( cos(phi)*drphi, 2) );

  C(kIdxX0, kIdxX0) = 10*pow(dx, 2);
  C(kIdxY0, kIdxY0) = 10*pow(dy, 2);
  C(kIdxTX, kIdxTX) = 0.001;
  C(kIdxTY, kIdxTY) = 0.001;
  C(kIdxQP, kIdxQP) = 0.0025;

  SoLKalTrackSite& initSite = *new SoLKalTrackSite(thisSeed->hitb, kMdim, kSdim, kMdim*fChi2PerNDFCut);

  initSite.Add(new SoLKalTrackState(svd, C, initSite, SoLKalTrackSite::kPredicted));
  initSite.Add(new SoLKalTrackState(svd, C, initSite, SoLKalTrackSite::kFiltered));
  initSite.SetHitResolution(kGiga, kGiga); //give it a very large resolution (100m) since it is a virtual site

  return initSite;
}
//______________________________________________________________________________________
inline SoLIDGEMHit* PVDISKalTrackFinder::FindCloestHitInWindow(double &x, double &y){
  double minD = kGiga;
  SoLIDGEMHit *minHit = nullptr;
  for (unsigned int i=0; i<fWindowHits.size(); i++){
    double r = sqrt(pow(x - fWindowHits.at(i)->GetX(), 2) + pow(y - fWindowHits.at(i)->GetY(), 2));
     if (r < minD) {
       minHit = fWindowHits.at(i);
       minD = r;
     }
   }
   return minHit;
}
//______________________________________________________________________________________
inline int PVDISKalTrackFinder::GetHitsInWindow(int plane, double x, double wx, double y, double wy, bool flag)
{
  assert(plane >= 0);
  fWindowHits.clear();
  vector<Int_t> ChamberList;

  double thisR = sqrt(x*x + y*y);

  for (int i=0; i<fGEMTracker[plane]->GetNChamber(); i++){
    TSeqCollection* HitArray = fGEMTracker[plane]->GetChamber(i)->GetHits();
    for (int nhit = 0; nhit < HitArray->GetLast()+1; nhit++){
      SoLIDGEMHit *hit = (SoLIDGEMHit*)HitArray->At(nhit);

      if (hit->IsUsed()) continue;
      if (hit->GetR() < thisR - 0.03) continue;
      if (hit->GetR() > thisR + 0.03) break;

      bool condition;
      if (!flag) condition = sqrt( pow(hit->GetX() - x, 2) + pow(hit->GetY() - y, 2) ) < 0.05 ;
      else condition = ( fabs(hit->GetX() - x) < 10.*sqrt(wx) && fabs(hit->GetY() - y) < 10.*sqrt(wy) );

      if (condition){
        fWindowHits.push_back(hit);
        if (fWindowHits.size() > MAXWINDOWHIT) return -1; //too many hits to be considered
      }
    }
  }
  return fWindowHits.size();
}
//______________________________________________________________________________________
inline double PVDISKalTrackFinder::CalDeltaPhi(const double & phi1, const double & phi2)
{
  double deltaPhi = phi1 - phi2;
  return TVector2::Phi_mpi_pi(deltaPhi);
}
//_____________________________________________________________________________________
inline double PVDISKalTrackFinder::CalDeltaR(const double & r1, const double & r2)
{
  return r1 - r2;
}
//______________________________________________________________________________________
inline void PVDISKalTrackFinder::Rotate(double& x, double& y)
{
  Double_t tempx = x;
  x = fRefCos*x     - fRefSin*y;
  y = fRefSin*tempx + fRefCos*y;
}
//______________________________________________________________________________________
inline Double_t PVDISKalTrackFinder::StraightLinePredict
(const Double_t& x1, const Double_t& z1, const Double_t& x2, const Double_t& z2, const Double_t& targetZ)
{
  return (x1-x2)/(z1-z2)*(targetZ - z1) +x1;
}
//______________________________________________________________________________________
inline Bool_t PVDISKalTrackFinder::CheckChargeAsy(SoLKalTrackSystem* theSystem)
{
  Int_t countCharge = 0;

  for (Int_t i=1; i<theSystem->GetLast()+1; i++){
    SoLIDGEMHit* theHit = (SoLIDMCGEMHit*)(static_cast<SoLKalTrackSite*>(theSystem->At(i))->GetHit());
    if (fabs( (theHit->GetQU() - theHit->GetQV())/(theHit->GetQU() + theHit->GetQV()) ) < 0.6) countCharge++;
  }
  if (countCharge >= 3){
    return kTRUE;
  }
  else{
    return kFALSE;
  }
}
//_______________________________________________________________________________________
inline Double_t PVDISKalTrackFinder::FindVertexZ(SoLKalTrackState* thisState)
{
  assert(fabs(thisState->GetZ0() - fTargetCenter) < 0.01);
  Double_t vertexz = fTargetCenter+ (1./(pow((*thisState)(kIdxTX, 0),2) + pow((*thisState)(kIdxTY, 0),2)))*
  ((*thisState)(kIdxTX,0)*(fBPMX-(*thisState)(kIdxX0,0)) + (*thisState)(kIdxTY,0)*(fBPMY-(*thisState)(kIdxY0, 0)));

  return vertexz;
}
//_______________________________________________________________________________________
inline void PVDISKalTrackFinder::CopyTrack(SoLIDTrack* soltrack, SoLKalTrackSystem* kaltrack)
{
  soltrack->SetStatus(kaltrack->GetTrackStatus());
  soltrack->SetCoarseFitStatus(kTRUE);
  soltrack->SetAngleFlag(kaltrack->GetAngleFlag());
  soltrack->SetCharge(kaltrack->GetCharge());
  soltrack->SetNDF(kaltrack->GetNDF());
  soltrack->SetCoarseChi2(kaltrack->GetChi2());
  soltrack->SetMomentum(kaltrack->GetMomentum());
  soltrack->SetVertexZ(kaltrack->GetVertexZ());
  soltrack->SetPhi(kaltrack->GetPhi());
  soltrack->SetTheta(kaltrack->GetTheta());
  soltrack->SetMomMax(kaltrack->fDeltaECX);
  soltrack->SetMomMin(kaltrack->fDeltaECY);
  soltrack->SetThetaMin(kaltrack->fDeltaECE);

  //start from 1 because the 0th is the dummy site that we used to initialize Kalman Filter
  //TODO remember not to add the last one since later it will be the BPM, not GEM hit

  for (Int_t j=1; j!=kaltrack->GetLast()+1;j++){
    SoLIDGEMHit* thishit = 0;
    Int_t layer = 0;

    thishit = (SoLIDGEMHit*)(static_cast<SoLKalTrackSite*>(kaltrack->At(j))->GetPredInfoHit());
    //thishit->SetUsed();
    layer = thishit->GetTrackerID();

    map< Int_t, vector<SoLIDGEMHit*> >::iterator it = fGoodHits.find(layer);

    if (it != fGoodHits.end()){
      (it->second).push_back(thishit);
    }
    else{
      vector<SoLIDGEMHit*> thisVector;
      thisVector.push_back(thishit);
      fGoodHits.insert(std::pair<Int_t, vector<SoLIDGEMHit*> >(layer, thisVector));
    }

    assert(thishit != 0);
    soltrack->AddHit(thishit);
  }
}






