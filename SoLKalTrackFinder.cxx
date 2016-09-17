//c++
#include <cmath>
//ROOT
#include "TRandom.h"
//SoLIDTracking
#include "SoLKalTrackFinder.h"
#include "SoLKalTrackSystem.h"
#include "SoLKalTrackSite.h"
#include "SoLKalTrackState.h"

//these should definitely need to go to the database
#define MAXNTRACKS_FAEC 1000
#define MAXNTRACKS_LAEC 1000
#define MAXNSEEDS 2000
#define MAXWINDOWHIT 1000


SoLKalTrackFinder::SoLKalTrackFinder(bool isMC, Int_t ntrackers)
:fIsMC(isMC), fECal(nullptr), fNSeeds(0), fSeedEfficiency(false), fMcTrackEfficiency(false),
 fNTrackers(ntrackers), fEventNum(0), fBPMX(0.), fBPMY(0.), fNGoodTrack(0), fChi2PerNDFCut(30.)
{
  fGEMTracker.clear();
  fCoarseTracks = new TClonesArray("SoLKalTrackSystem", MAXNTRACKS_FAEC, kTRUE);
  
  fFieldStepper = SoLKalFieldStepper::GetInstance();
  fWindowHits.clear();
  fWindowHits.reserve(MAXWINDOWHIT);
  fTargetPlaneZ = -3.2;
  fTargetCenter = -3.5;
  fTargetLength =  0.4;
  
  vector<DoubletSeed> midBackSeed;
  midBackSeed.clear();
  vector<DoubletSeed> frontMidSeed;
  frontMidSeed.clear();
  vector<DoubletSeed> frontBackSeed;
  frontBackSeed.clear();
  
  fSeedPool[kMidBack] = midBackSeed;
  fSeedPool[kFrontMid] = frontMidSeed;
  fSeedPool[kFrontBack] = frontBackSeed;
}
//__________________________________________________________________________
SoLKalTrackFinder::~SoLKalTrackFinder()
{
  Clear();
  fGoodHits.clear();
  delete fCoarseTracks;
}
//___________________________________________________________________________
void SoLKalTrackFinder::Clear( Option_t* opt )
{
  
  if (fCoarseTracks->GetEntries() != 0) 
   fCoarseTracks->Delete();
  
  fCoarseTracks->Clear(opt);
  
  fCaloHits.clear();
  fNSeeds = 0;
  fNGoodTrack = 0;
  fSeedEfficiency = false;
  fMcTrackEfficiency = false;
  
  map< Int_t, vector<SoLIDGEMHit*> >::iterator it;
  for (it = fGoodHits.begin(); it != fGoodHits.end(); it++) { (it->second).clear(); }
  fGoodHits.clear();
  
  map< SeedType, vector<DoubletSeed> >::iterator itt;
  for (itt = fSeedPool.begin(); itt != fSeedPool.end(); itt++) { (itt->second).clear(); }
}
//___________________________________________________________________________
void SoLKalTrackFinder::SetGEMDetector(vector<SoLIDGEMTracker*> thetrackers)
{
  fGEMTracker = thetrackers;
}
//___________________________________________________________________________
void SoLKalTrackFinder::ProcessHits(TClonesArray* theTracks)
{
  if (fGEMTracker.size() == 0) return;
  
  fNSeeds = 0;
  
  //forward angle seed finding
  FindDoubletSeed(4, 5, kFAEC);
  FindDoubletSeed(3, 4, kFAEC);
  FindDoubletSeed(3, 5, kFAEC);
  MergeSeed();
  
  map< SeedType, vector<DoubletSeed> >::iterator itt;
  for (itt = fSeedPool.begin(); itt != fSeedPool.end(); itt++) { (itt->second).clear(); }
  
  //large angle seed finding
  //FindDoubletSeed(2, 3, kLAEC);
  //FindDoubletSeed(1, 2, kLAEC);
  //FindDoubletSeed(1, 3, kLAEC);
  //MergeSeed();
  
  
  TrackFollow();
  FindandAddVertex();
  ECalFinalMatch();
  FinalSelection(theTracks);
  
  fEventNum++;
}

//___________________________________________________________________________________________________________________
void SoLKalTrackFinder::FindDoubletSeed(Int_t planej, Int_t planek, ECType type)
{
  double philimit[2] = {0};
  double rlimit[2][2] = {{0, 0}, {0, 0}};
  double deltar[2] = {0};
  double dphi, dr;
  SeedType seedType = kMidBack;
  
  double charge = 0;
  if (type == kFAEC){
    if (planek == 5 && planej == 4){
      seedType = kMidBack;
      rlimit[1][0] = 0.43; rlimit[1][1] = 0.95;
      rlimit[0][0] = 0.54; rlimit[0][1] = 1.19; 
      deltar[0]=0.084;deltar[1]=0.241; 
      philimit[0]=0.022; philimit[1]=0.207;
    }
    else if(planek == 4 && planej == 3){
      seedType = kFrontMid;
      rlimit[1][0] = 0.33; rlimit[1][1] = 0.78; 
      rlimit[0][0] = 0.43; rlimit[0][1] = 0.95;
      deltar[0]=0.081; deltar[1]=0.201;
      philimit[0]=0.017; philimit[1]=0.173;
    }
    else if (planek == 5 && planej == 3){
      seedType = kFrontBack;
      rlimit[1][0] = 0.33; rlimit[1][1] = 0.78;
      rlimit[0][0] = 0.54; rlimit[0][1] = 1.19;
      deltar[0]=0.157; deltar[1]=0.439;
      philimit[0]=0.042; philimit[1]=0.378;
    } 
  }
  else{ 
    if (planek == 3 && planej == 2){
      seedType = kMidBack;
      rlimit[1][0] = 0.60; rlimit[1][1] = 1.12;
      rlimit[0][0] = 0.74; rlimit[0][1] = 1.35;
      deltar[0]=0.121; deltar[1]=0.242;
      philimit[0]=0.012; philimit[1]=0.127;
    }
    else if (planek == 2 && planej == 1){
      seedType = kFrontMid;
      rlimit[1][0] = 0.51; rlimit[1][1] = 0.98;
      rlimit[0][0] = 0.60; rlimit[0][1] = 1.12;
      deltar[0]=0.076; deltar[1]=0.151;
      philimit[0]=0.006;philimit[1]=0.07;
    }
    else if (planek == 3 && planej == 1){
      seedType = kFrontBack;
      rlimit[1][0] = 0.51; rlimit[1][1] = 0.98;
      rlimit[0][0] = 0.74; rlimit[0][1] = 1.35;
      deltar[0]=0.2; deltar[1]=0.39;
      philimit[0]=0.018;philimit[1]=0.2;
    }
    
  }
  
  
  for (int k=0; k<fGEMTracker[planek]->GetNChamber(); k++){
    TSeqCollection* planekHitArray = fGEMTracker[planek]->GetChamber(k)->GetHits();
  
    for (int nhitk = 0; nhitk < planekHitArray->GetLast()+1; nhitk++){
      SoLIDGEMHit *hitk = (SoLIDGEMHit*)planekHitArray->At(nhitk);
      if (hitk->GetR() < rlimit[0][0]) continue;
      if (hitk->GetR() > rlimit[0][1]) break; // check if the hit is within r range
      if (!TriggerCheck(hitk, type)) continue;
      
      vector<Int_t> jChamberList;
      GetHitChamberList(jChamberList, k, 3);
      
      for (int j=0; j<(int)jChamberList.size(); j++){
        TSeqCollection* planejHitArray = fGEMTracker[planej]->GetChamber(jChamberList.at(j))->GetHits();
      
        for (int nhitj = 0; nhitj < planejHitArray->GetLast()+1; nhitj++){
          SoLIDGEMHit *hitj = (SoLIDGEMHit*)planejHitArray->At(nhitj);
          if (hitj->IsUsed()) continue;
          if (hitj->GetR()<rlimit[1][0]) continue; 
          if (hitj->GetR()>rlimit[1][1]) break;
          
          dr = CalDeltaR(hitk->GetR(), hitj->GetR());
          if (dr > deltar[1]) continue;
          if (dr < deltar[0]) break;
          
          charge = 0;
	        dphi = CalDeltaPhi(hitj->GetPhi(), hitk->GetPhi());
	        if(((dphi >philimit[0]&& dphi <philimit[1])||(dphi < -1*philimit[0]&& dphi > -1*philimit[1]))){
	         
	          if (dphi >philimit[0]&& dphi <philimit[1]){
	              charge = 1;
	          }else{
	              charge = -1;
	          }
          }
          else continue;
          
          assert(charge != 0); //should never happen
          
          //using correction function to calculate initial momentum and angles of the particle at plane k
          double initTheta = 0;
          double initMom   = 0;
          double initPhi   = 0;
          if (!CalInitParForPair(hitj, hitk, charge, initMom, initTheta, initPhi, type)) continue;
          
          if (type == kFAEC && (initTheta > 0.3 || initTheta < 0.1)) continue;
          if (type == kLAEC && (initTheta > 0.5 || initTheta < 0.24)) continue;
          if (initMom > 8. || initMom < 0.8) continue;
          
          
          TVector3 initDir(cos(initPhi), sin(initPhi), 1./tan(initTheta));
          initDir = initDir.Unit();
		      TVector3 initMomentum = initMom*initDir;
		      TVector3 initPosition(hitk->GetX(), hitk->GetY(), hitk->GetZ());
		      TVector3 finalMomentum;
          TVector3 finalPosition;
          Double_t stepSize = 1.;
          
          Bool_t isSeed = false;
          Double_t toZ = fECal->GetECZ(type);
		      if (type == kFAEC){
        
            fFieldStepper->PropagationClassicalRK4(initMomentum, initPosition, toZ, 
                                               charge, stepSize, finalMomentum, finalPosition);
            for (UInt_t ec_count=0; ec_count<fCaloHits.size(); ec_count++){
	            if (fCaloHits[ec_count].fECID != kFAEC) continue; //not FAEC hit
	            if (sqrt( pow(finalPosition.X() - fCaloHits[ec_count].fXPos, 2) +  
	                  pow(finalPosition.Y() - fCaloHits[ec_count].fYPos, 2) ) < 0.2 ) isSeed = true;
	          } 
		      }
		      else if (type == kLAEC){
            fFieldStepper->PropagationClassicalRK4(initMomentum, initPosition, toZ, 
                                               charge, stepSize, finalMomentum, finalPosition);
            for (UInt_t ec_count=0; ec_count<fCaloHits.size(); ec_count++){
	            if (fCaloHits[ec_count].fECID != kLAEC) continue; //not FAEC hit
	            if (sqrt( pow(finalPosition.X() - fCaloHits[ec_count].fXPos, 2) +  
	                  pow(finalPosition.Y() - fCaloHits[ec_count].fYPos, 2) ) < 0.06 ) isSeed = true;
	          } 
		      }
		      if (type == kFAEC && !isSeed) continue;
          
          
          fFieldStepper->PropagationClassicalRK4(initMomentum, initPosition, 
                                             fTargetCenter, charge, stepSize, finalMomentum, finalPosition);
      
          double tx = finalMomentum.X()/finalMomentum.Z();
          double ty = finalMomentum.Y()/finalMomentum.Z();
		      double ReconZ = fTargetCenter + (1./(pow(tx,2) + pow(ty,2)))*
                       (tx*(fBPMX-finalPosition.X()) + ty*(fBPMY-finalPosition.Y()) );
      
          if (type == kFAEC && (ReconZ > fTargetCenter + 0.5 || ReconZ < fTargetCenter - 0.5) ) continue;
          if (type == kLAEC && (ReconZ > fTargetCenter + 0.4 || ReconZ < fTargetCenter - 0.4) ) continue;
          
          
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
          
		      if (dynamic_cast<SoLIDMCGEMHit*>(hitj)->IsSignalHit() && 
		        dynamic_cast<SoLIDMCGEMHit*>(hitk)->IsSignalHit() ){
		          fSeedEfficiency = true;
		      }
		      
          
        }
      }
      
    }
  }
  
}
//___________________________________________________________________________________________________________________
void SoLKalTrackFinder::MergeSeed()
{ 
  //here we will merge the doublet seed into a triplet seed, for which the three type of doublet seed must match at the 
  //common plane. Once a triplet seed is form, its corresponding doublet seeds will be deactivated
  
  for (unsigned int i=0; i<fSeedPool[kMidBack].size(); i++){
    for (unsigned int j=0; j<fSeedPool[kFrontMid].size(); j++){
      //see if there is a common point in between the two doublet seeds
      if (fSeedPool[kMidBack].at(i).hita == fSeedPool[kFrontMid].at(j).hitb){
        //now check the third type of seed
        for (unsigned int k=0; k<fSeedPool[kFrontBack].size(); k++){
          if (fSeedPool[kFrontBack].at(k).hita == fSeedPool[kFrontMid].at(j).hita && 
            fSeedPool[kFrontBack].at(k).hitb == fSeedPool[kMidBack].at(i).hitb){
              
            fSeedPool[kMidBack].at(i).Deactive();
            fSeedPool[kFrontMid].at(j).Deactive();
            fSeedPool[kFrontBack].at(k).Deactive();
            
            SoLKalTrackSite & initSite =  SiteInitWithSeed(&(fSeedPool[kMidBack].at(i)));
            SoLKalTrackSystem *thisSystem = new ((*fCoarseTracks)[fNSeeds++]) SoLKalTrackSystem();
            thisSystem->SetMass(kElectronMass);
            thisSystem->SetCharge(fSeedPool[kMidBack].at(i).charge);
            thisSystem->SetElectron(kTRUE);
            thisSystem->SetAngleFlag(fSeedPool[kMidBack].at(i).flag);
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
//___________________________________________________________________________________________________________________
void SoLKalTrackFinder::TrackFollow()
{
  //this function is responsible for propagating the seed track toward the next tracker, find suitable hits
  //the process stop until the track reach the first tracker upstream (track searching always go backward)
  
  for (Int_t i=0; i<fCoarseTracks->GetLast()+1; i++){
    SoLKalTrackSystem* thisSystem = (SoLKalTrackSystem*)(fCoarseTracks->At(i));
    thisSystem->CheckTrackStatus();
    if (!thisSystem->GetTrackStatus()) continue; //skip the bad tracks
    //--------------------for test----------------------//
      /*SoLKalTrackSystem *newSystem = (SoLKalTrackSystem*)thisSystem->Clone();
      cout<<"chi2: "<<newSystem->GetChi2()<<" "<<thisSystem->GetChi2()<<endl;
      cout<<"mass: "<<newSystem->GetMass()<<" "<<thisSystem->GetMass()<<endl;
      cout<<"state x: "<<(newSystem->GetCurSite()).GetCurState()(0, 0)<<" "<<(thisSystem->GetCurSite()).GetCurState()(0, 0)<<endl;
      cout<<"state x: "<<(newSystem->GetCurSite()).GetCurState()(1, 0)<<" "<<(thisSystem->GetCurSite()).GetCurState()(1, 0)<<endl;
      cout<<"state x: "<<(newSystem->GetCurSite()).GetCurState().GetZ0()<<" "<<(thisSystem->GetCurSite()).GetCurState().GetZ0()<<endl;*/
    //--------------------------------------------------//
    
    thisSystem->SetCurInstancePtr(thisSystem);
    Int_t currentTracker = ((thisSystem->GetCurSite()).GetHit())->GetTrackerID();
    
    Int_t lastTracker = 0;
    
    //seed from type kMidBack will skip the front seed plane. We assume for this type of seed, the hit on the
    //front seed plane is missing, (otherwise the seed should be absorbed into the triplet seed)
    
    if (thisSystem->GetSeedType() == kMidBack) currentTracker--;
    
    while (currentTracker > lastTracker){
      currentTracker--;
      
      thisSystem->CheckTrackStatus();    
      if (!thisSystem->GetTrackStatus()) break; //skip the bad tracks
      
      SoLKalTrackState currentState = (thisSystem->GetCurSite()).GetCurState();
      currentState.InitPredictSV();
      SoLKalTrackState *predictState = currentState.PredictSVatNextZ(fGEMTracker[currentTracker]->GetZ());
      
      bool flag = (thisSystem->GetNHits() >= 3);
      
      int size = GetHitsInWindow(currentTracker, (*predictState)(kIdxX0, 0), (predictState->GetCovMat())(kIdxX0, kIdxX0),
                                (*predictState)(kIdxY0, 0), (predictState->GetCovMat())(kIdxY0, kIdxY0), flag); 
                                
      
      if (size <= 0){
        //when there are too many hits in a small window (usually should not happen), or
        //when there is no hit found in the window, if the track has not missed a hit so far
        //we will keep the track, otherwise it is a bad track (miss too many hits)
        
        //the only exception will be the 0th tracker for a forward angle track, for which it 
        //is not necessary to have a hit, and in that case, we don't count missing hit
        if (thisSystem->GetAngleFlag() == kFAEC && currentTracker == 0) continue;
        
        thisSystem->AddMissingHits();
        
        /*if (thisSystem->GetNMissingHits() > 1 ){
          thisSystem->SetTrackStatus(kFALSE);
        }*/
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
    //------------check MC track efficiency--------------//
    bool allMC = true;
    for (Int_t j=1; j!=thisSystem->GetLast()+1;j++){
      SoLIDGEMHit* thisHit = (SoLIDGEMHit*)((SoLKalTrackSite*)thisSystem->At(j))->GetHit();
      if (!dynamic_cast<SoLIDMCGEMHit*>(thisHit)->IsSignalHit()) allMC = false;
      break;
    }
    if (allMC) fMcTrackEfficiency = true;
    //---------------------------------------------------//
  }
}
//___________________________________________________________________________________________________________________
void SoLKalTrackFinder::FindandAddVertex()
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
      
      if (thisSystem->GetAngleFlag() == kFAEC && fabs(vertexz - fTargetCenter) > 0.3){
        thisSystem->SetTrackStatus(kFALSE);
        continue;
      }
      else if (thisSystem->GetAngleFlag() == kLAEC && fabs(vertexz - fTargetCenter) > 0.25){
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
//___________________________________________________________________________________________________________________
void SoLKalTrackFinder::FinalSelection(TClonesArray *theTracks)
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
//___________________________________________________________________________________________________________________
void SoLKalTrackFinder::ECalFinalMatch()
{
  for (Int_t i=0; i<fCoarseTracks->GetLast()+1; i++){
    SoLKalTrackSystem *thisSystem = (SoLKalTrackSystem*)(fCoarseTracks->At(i));
    
    if ( !(thisSystem->GetTrackStatus()) ) continue; //skip bad tracks
    Double_t ecalZ = fECal->GetECZ((ECType)thisSystem->GetAngleFlag());
    
    //using Kalman Filter smoother to smooth the track back to the first measurement site
    //so that we don't need to propagate and fit back again
    thisSystem->SetCurInstancePtr(thisSystem);
    thisSystem->SmoothBackTo(1);
    
    SoLKalTrackState currentState = (thisSystem->GetCurSite()).GetCurState();
    currentState.InitPredictSV();
    
    SoLKalTrackState *predictState = currentState.PredictSVatNextZ(ecalZ);
    
    thisSystem->SetTrackStatus(kFALSE);
    for (UInt_t ec_count=0; ec_count<fCaloHits.size(); ec_count++){
	    if (fCaloHits.at(ec_count).fECID != thisSystem->GetAngleFlag()) continue;
	    if ( fabs(fCaloHits.at(ec_count).fXPos - (*predictState)(kIdxX0, 0)) <5.*0.01 &&
	         fabs(fCaloHits.at(ec_count).fYPos - (*predictState)(kIdxY0, 0)) <5.*0.01  ){
	      
	      //for large angle, require also that the momentum of the track needs to match the 
	      //cluster energy. This is difficult to do for forward angle since we detect both hadron
	      //and electron there and there is a long distance betwee the FAEC and the last GEM, during
	      //which there could be significant energy loss but we don't have other tracking detectors
	      //and field integral to measure it
	      Double_t momentum = thisSystem->GetCharge() / (*predictState)(kIdxQP, 0);
	      if (thisSystem->GetAngleFlag() == kLAEC){
	        if ( fabs( (momentum - fCaloHits.at(ec_count).fEdp) / momentum) > 0.5) continue;
	      }
	      
	      thisSystem->SetTrackStatus(kTRUE);
	      thisSystem->fDeltaECX = fCaloHits.at(ec_count).fXPos - (*predictState)(kIdxX0, 0);
	      thisSystem->fDeltaECY = fCaloHits.at(ec_count).fYPos - (*predictState)(kIdxY0, 0);
	      thisSystem->fDeltaECE = (momentum - fCaloHits.at(ec_count).fEdp)/momentum;    
	    }
	  }
	  
	  thisSystem->SetSitePtrToLastSite();
  }
}
//___________________________________________________________________________________________________________________

//assistant functions below
inline SoLKalTrackSite & SoLKalTrackFinder::SiteInitWithSeed(DoubletSeed* thisSeed)
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
  /*for (int index=0; index<kSdim; index++) {
    C(index,index) = 1.;   // dummy error matrix
  }*/
  Double_t phi = (thisSeed->hitb)->GetPhi();
  Double_t dr = 4.e-4;
  Double_t drphi = 4.5e-5;
  
  Double_t dx = sqrt( pow( cos(phi)*dr, 2) + pow( sin(phi)*drphi, 2) );
  Double_t dy = sqrt( pow( sin(phi)*dr, 2) + pow( cos(phi)*drphi, 2) );
  
  C(kIdxX0, kIdxX0) = 10*pow(dx, 2);
  C(kIdxY0, kIdxY0) = 10*pow(dy, 2);
  C(kIdxTX, kIdxTX) = 0.001; 
  C(kIdxTY, kIdxTY) = 0.001;
  C(kIdxQP, kIdxQP) = 0.01;
      
  SoLKalTrackSite& initSite = *new SoLKalTrackSite(thisSeed->hitb, kMdim, kSdim, kMdim*fChi2PerNDFCut);
      
  initSite.Add(new SoLKalTrackState(svd, C, initSite, SoLKalTrackSite::kPredicted));
  initSite.Add(new SoLKalTrackState(svd, C, initSite, SoLKalTrackSite::kFiltered));
  initSite.SetHitResolution(kGiga, kGiga); //give it a very large resolution (100m) since it is a virtual site
  
  return initSite;
}
//___________________________________________________________________________________________________________________
inline Bool_t SoLKalTrackFinder::TriggerCheck(SoLIDGEMHit *theHit, ECType type)
{
  if (type == kLAEC){
    for (UInt_t ec_count=0; ec_count<fCaloHits.size(); ec_count++){
	    if (fCaloHits[ec_count].fECID != kLAEC) continue; //not LAEC hit
	    Double_t ecHitPhi = TMath::ATan2(fCaloHits[ec_count].fYPos, fCaloHits[ec_count].fXPos);
	    Double_t ecHitR = TMath::Sqrt( TMath::Power(fCaloHits[ec_count].fXPos, 2) + 
	    			  	   TMath::Power(fCaloHits[ec_count].fYPos, 2) );
	    Double_t tmpDeltaPhi = CalDeltaPhi(ecHitPhi, theHit->GetPhi());
	    Double_t tmpDeltaR   = CalDeltaR(ecHitR, theHit->GetR());
	    if (theHit->GetTrackerID()==2 && (tmpDeltaPhi < 0.15 && tmpDeltaPhi > -0.05) && (tmpDeltaR < 0.286 && tmpDeltaR > 0.095)){
        return kTRUE;
     }else if (theHit->GetTrackerID()==3 && (tmpDeltaPhi<0.06 && tmpDeltaPhi>-0.06) && (tmpDeltaR<0.054 && tmpDeltaR > -0.039) ){
       return kTRUE;
     }else{
	     return kFALSE;
	   }
	  }
  }
  else if (type == kFAEC){
    for (UInt_t ec_count=0; ec_count<fCaloHits.size(); ec_count++)
	    {
	      if (fCaloHits[ec_count].fECID != kFAEC) continue; //not FAEC hit
	      Double_t ecHitPhi = TMath::ATan2(fCaloHits[ec_count].fYPos, fCaloHits[ec_count].fXPos);
	      Double_t ecHitR = TMath::Sqrt( TMath::Power(fCaloHits[ec_count].fXPos, 2) + 
				  	   TMath::Power(fCaloHits[ec_count].fYPos, 2) );
	      Double_t tmpDeltaPhi = CalDeltaPhi(ecHitPhi, theHit->GetPhi());
	      Double_t tmpDeltaR   = CalDeltaR(ecHitR, theHit->GetR());
	      if (theHit->GetTrackerID()==4 && (tmpDeltaPhi < 1.1 && tmpDeltaPhi > 0.04) && (tmpDeltaR < 1.11 && tmpDeltaR > 0.425)){
	        return kTRUE;
	      }else if (theHit->GetTrackerID()==5 && (tmpDeltaPhi<0.9 && tmpDeltaPhi>0.02) && (tmpDeltaR<0.88 && tmpDeltaR > 0.31) ){
	        return kTRUE;
	      }else{
	        return kFALSE;
	      }
	    }
  }
  return kFALSE;
}
//___________________________________________________________________________________________________________________
inline void SoLKalTrackFinder::CalCircle(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x3,
                                         Double_t y3, Double_t* R,Double_t* Xc, Double_t* Yc){
  if (y1!=y2&& (x1-x2)*(y2-y3)!=(x2-x3)*(y1-y2)){
    *Xc = ((y2-y3)*(x1*x1-x2*x2+y1*y1-y2*y2)/2.-(y1-y2)*(x2*x2-x3*x3+y2*y2-y3*y3)/2.)/((x1-x2)*(y2-y3)-(x2-x3)*(y1-y2));
    *Yc = ((x1*x1-x2*x2+y1*y1-y2*y2)/2.-(*Xc)*(x1-x2))/(y1-y2);
    *R = sqrt((*Xc-x1)*(*Xc-x1)+(*Yc-y1)*(*Yc-y1));
    
  }else if (y2!=y3&&(x2-x3)*(y3-y1)!=(x3-x1)*(y2-y3)){
    *Xc = ((y3-y1)*(x2*x2-x3*x3+y2*y2-y3*y3)/2.-(y2-y3)*(x3*x3-x1*x1+y3*y3-y1*y1)/2.)/((x2-x3)*(y3-y1)-(x3-x1)*(y2-y3));
    *Yc = ((x2*x2-x3*x3+y2*y2-y3*y3)/2.-(*Xc)*(x2-x3))/(y2-y3);
    *R = sqrt((*Xc-x2)*(*Xc-x2)+(*Yc-y2)*(*Yc-y2));
   
  }else if (y3!=y1&&(x3-x1)*(y1-y2)!=(x1-x2)*(y3-y1)){
    *Xc = ((y1-y2)*(x3*x3-x1*x1+y3*y3-y1*y1)/2.-(y3-y1)*(x1*x1-x2*x2+y1*y1-y2*y2)/2.)/((x3-x1)*(y1-y2)-(x1-x2)*(y3-y1));
    *Yc = ((x3*x3-x1*x1+y3*y3-y1*y1)/2.-(*Xc)*(x3-x1))/(y3-y1);
    *R = sqrt((*Xc-x3)*(*Xc-x3)+(*Yc-y3)*(*Yc-y3));
    
  }else{
    *Xc = 0;
    *Yc = 0;
    *R = 0;
  }
}
//___________________________________________________________________________________________________________________
inline double SoLKalTrackFinder::CalDeltaPhi(double phi1, double phi2)
{
  double deltaPhi = phi1 - phi2;
  deltaPhi = TVector2::Phi_mpi_pi(deltaPhi);
  return deltaPhi;
}
//__________________________________________________________________________________________________________________
inline double SoLKalTrackFinder::CalDeltaR(double r1, double r2)
{
  return r1 - r2;
}
//___________________________________________________________________________________________________________________
inline SoLIDGEMHit* SoLKalTrackFinder::FindCloestHitInWindow(double &x, double &y){
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
//___________________________________________________________________________________________________________________
inline int SoLKalTrackFinder::GetHitsInWindow(int plane, double x, double wx, double y, double wy, bool flag)
{
  assert(plane >= 0);
  fWindowHits.clear();
  vector<Int_t> ChamberList;
  
  double thisR = sqrt(x*x + y*y);
  
  
  GetHitChamberList(ChamberList, GetChamIDFromPos(x, y, plane), 1);
  for (int i=0; i<(int)ChamberList.size(); i++){
    TSeqCollection* HitArray = fGEMTracker[plane]->GetChamber(ChamberList.at(i))->GetHits();
    for (int nhit = 0; nhit < HitArray->GetLast()+1; nhit++){ 
      SoLIDGEMHit *hit = (SoLIDGEMHit*)HitArray->At(nhit);
      
      if (hit->IsUsed()) continue;
      if (hit->GetR() < thisR - 0.03) continue;
      if (hit->GetR() > thisR + 0.03) break;     
      
      bool condition;
      if (!flag) condition = sqrt( pow(hit->GetX() - x, 2) + pow(hit->GetY() - y, 2) ) < 0.015 ;
      else condition = ( fabs(hit->GetX() - x) < 10.*sqrt(wx) && fabs(hit->GetY() - y) < 10.*sqrt(wy) ); 
      
      
      
      if (condition){
        fWindowHits.push_back(hit);
        if (fWindowHits.size() > MAXWINDOWHIT) return -1; //too many hits to be considered       
      }        
    }      
  }
  
  return fWindowHits.size();
}
//____________________________________________________________________________________________________________________
inline Double_t SoLKalTrackFinder::FindVertexZ(SoLKalTrackState* thisState)
{
  assert(fabs(thisState->GetZ0() - fTargetCenter) < 0.01);
  Double_t vertexz = fTargetCenter+ (1./(pow((*thisState)(kIdxTX, 0),2) + pow((*thisState)(kIdxTY, 0),2)))*
  ((*thisState)(kIdxTX,0)*(fBPMX-(*thisState)(kIdxX0,0)) + (*thisState)(kIdxTY,0)*(fBPMY-(*thisState)(kIdxY0, 0)));
  
  return vertexz;
}
//____________________________________________________________________________________________________________________
inline void SoLKalTrackFinder::CopyTrack(SoLIDTrack* soltrack, SoLKalTrackSystem* kaltrack)
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
//______________________________________________________________________________________________________________________
inline Bool_t SoLKalTrackFinder::CheckChargeAsy(SoLKalTrackSystem* theSystem)
{
  Int_t countCharge = 0;
  
  for (Int_t i=1; i<theSystem->GetLast()+1; i++){
    SoLIDGEMHit* theHit = (SoLIDMCGEMHit*)(static_cast<SoLKalTrackSite*>(theSystem->At(i))->GetHit());
    if (fabs( (theHit->GetQU() - theHit->GetQV())/(theHit->GetQU() + theHit->GetQV()) ) < 0.4) countCharge++;
  }
  if (countCharge >= 3){
    return kTRUE;
  }
  else{
    return kFALSE;
  }
}
//_______________________________________________________________________________________________________________________
inline void SoLKalTrackFinder::GetHitChamberList(vector<Int_t> &theList, Int_t thisChamber, Int_t size)
{
  //choose 7 chambers around the current hit chamber, in order to reduce the search region 
  theList.resize(2*size+1);
  thisChamber -= size;
  for (unsigned int i=0; i<theList.size(); i++){
    theList[i] = thisChamber;
    thisChamber++;
    if (theList[i] < 0) theList[i] += 30;
    if (theList[i] >29) theList[i] -= 30;
  }
}
//_______________________________________________________________________________________________________________________
inline Int_t SoLKalTrackFinder::GetChamIDFromPos(Double_t &x, Double_t &y, Int_t TrackerID)
{
  double phi = atan2(y, x);
  for (int i=0; i<fGEMTracker[TrackerID]->GetNChamber(); i++){
    SoLIDGEMChamber* thisChamber = fGEMTracker[TrackerID]->GetChamber(i);
    double dphi = phi - thisChamber->GetPhi();
    dphi = TVector2::Phi_mpi_pi(dphi);
    if ( dphi < thisChamber->GetPhiCover()/2. && dphi > -1*thisChamber->GetPhiCover()/2.) 
    return i;
  }
  return -1; 
}
//_______________________________________________________________________________________________________________________
void SoLKalTrackFinder::SetBPM(Double_t x, Double_t y)
{
  fBPMX = x + gRandom->Gaus(0., 3e-4);
  fBPMY = y + gRandom->Gaus(0., 3e-4);
}
//_______________________________________________________________________________________________________________________
inline double SoLKalTrackFinder::PredictR(Int_t &plane, SoLIDGEMHit* hit1, SoLIDGEMHit* hit2)
{
  double z = fGEMTracker[plane]->GetZ();
  return (hit1->GetR() - hit2->GetR())/(hit1->GetZ() - hit2->GetZ())*(z - hit1->GetZ()) 
         + hit1->GetR();
}
//________________________________________________________________________________________________________________________
inline Bool_t SoLKalTrackFinder::CalInitParForPair(SoLIDGEMHit* hita, SoLIDGEMHit* hitb, Double_t &charge,
                                                 Double_t& mom, Double_t& theta, Double_t& phi, ECType& type)
{
  double deltaR = hitb->GetR() - hita->GetR();
  double deltaPhi = TVector2::Phi_mpi_pi(hitb->GetPhi() - hita->GetPhi());
  deltaPhi = fabs(deltaPhi);
  
  double deltaZ = hitb->GetZ() - hita->GetZ();
  double deltaY = hitb->GetY() - hita->GetY();
  double deltaX = hitb->GetX() - hita->GetX();
  
  theta = atan(deltaR/deltaZ);
  phi = atan2(deltaY, deltaX);
  
  if (type == kFAEC && hita->GetTrackerID() == 4 && hitb->GetTrackerID() == 5){
    if (deltaPhi < -0.000260251 + 1e-6) return kFALSE;
    mom = 0.184566/(deltaPhi - -0.000260251);
    theta += 0.00117211 + -0.0519514*deltaPhi + 2.25303*deltaPhi*deltaPhi;
    phi +=(-1.*charge)*( 0.000617591 + 0.827572*deltaPhi + 0.53233*deltaPhi*deltaPhi );
  }
  else if (type == kFAEC && hita->GetTrackerID() == 3 && hitb->GetTrackerID() == 4){
    if (deltaPhi < -0.000456619 + 1e-6) return kFALSE;
    mom = 0.156196/(deltaPhi - -0.000456619);
    theta += 0.000614522 + -0.0317596*deltaPhi + 2.14209*deltaPhi*deltaPhi;
    phi += (-1*charge)*(-0.000117883 + 1.0794*deltaPhi + -0.178718*deltaPhi*deltaPhi);
  }
  else if (type == kFAEC && hita->GetTrackerID() == 3 && hitb->GetTrackerID() == 5){
    if (deltaPhi < -0.00068094 + 1e-6) return kFALSE;
    mom = 0.340493/(deltaPhi - -0.00068094);
    theta += 0.000992593 + -0.0239743*deltaPhi + 0.549664*deltaPhi*deltaPhi;
    phi += (-1*charge)*(0.000672179 + 0.908331*deltaPhi + 0.166936*deltaPhi*deltaPhi);
  }
  else if (type == kLAEC && hita->GetTrackerID() == 2 && hitb->GetTrackerID() == 3){
    if (deltaPhi < -0.000396581 + 1e-6) return kFALSE;
    mom = 0.109309/(deltaPhi - -0.000396581);
    theta += 3.61953e-05 + -1.06941e-05*deltaPhi + 4.23056*deltaPhi*deltaPhi;
    phi += (-1*charge)*(-0.000283361 + 1.2515*deltaPhi + -0.637376*deltaPhi*deltaPhi);
  }
  else if (type == kLAEC && hita->GetTrackerID() == 1 && hitb->GetTrackerID() == 2){
    if (deltaPhi < -0.000371282 + 1e-6) return kFALSE;
    mom = 0.0610608/(deltaPhi - -0.000371282);
    theta += 5.30086e-05 + -0.00477822*deltaPhi + 8.29537*deltaPhi*deltaPhi;
    phi += (-1*charge)*(-0.00018356 + 1.39156*deltaPhi + -1.29865*deltaPhi*deltaPhi);
  }
  else if (type == kLAEC && hita->GetTrackerID() == 1 && hitb->GetTrackerID() == 3){
    if (deltaPhi < -0.000743132 + 1e-6) return kFALSE;
    mom = 0.17019/(deltaPhi - -0.000743132);
    theta += -5.87987e-06 + 0.00195822*deltaPhi + 1.58111*deltaPhi*deltaPhi;
    phi += (-1*charge)*(-0.000508897 + 1.30589*deltaPhi + -0.468172*deltaPhi*deltaPhi);
  }
  
  phi = TVector2::Phi_mpi_pi(phi);
  return kTRUE;
  
}
//________________________________________________________________________________________________________________________
inline Int_t SoLKalTrackFinder::BinarySearchForR(TSeqCollection* array, Double_t &lowr)
{
  //search for the first index that is above lowr, if not found, return the length of the array
  //the array needs to be sorted in increasing r order before use
  
  Int_t low = 0;
  Int_t high = array->GetEntries();
  
  while (low != high){
    Int_t mid = (low + high)/2;
    SoLIDGEMHit* hit = (SoLIDGEMHit*)array->At(mid);
    if ( hit->GetR() <= lowr ){
      
      low = mid + 1;
    }
    else{
      high = mid;
    }
  }
  
  assert(low == high);
  return low;
}



