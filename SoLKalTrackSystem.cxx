//SoLIDTracking
#include "SoLKalTrackSystem.h"

ClassImp(SoLKalTrackSystem)

SoLKalTrackSystem *SoLKalTrackSystem::fgCurInstancePtr = 0;
//__________________________________________________________________
SoLKalTrackSystem::SoLKalTrackSystem(Int_t n)
            :TObjArray(n),
             fCurSitePtr(0),
             fChi2(0.), fIsGood(kTRUE), fNMissingHits(0), 
             fNHits(-1), fNDF(0), fSeedType(kTriplet)
{
}
//___________________________________________________________________
SoLKalTrackSystem::~SoLKalTrackSystem()
{
   if (this == fgCurInstancePtr) fgCurInstancePtr = 0;
   SetOwner(kTRUE);
   Delete();
}
//___________________________________________________________________
Bool_t SoLKalTrackSystem::AddAndFilter(SoLKalTrackSite &next)
{
   SetCurInstancePtr(this);

   //
   // Propagate current state to the next site
   //

   GetState(SoLKalTrackSite::kFiltered).Propagate(next);

   //
   // Calculate new pull and gain matrix
   //

   if (next.Filter()) {
      //
      // Add this to the system if accepted.
      //

      Add(&next);
      fChi2 += next.GetDeltaChi2();
      return kTRUE;
   } else {
      return kFALSE;
   }
}
//____________________________________________________________________
Int_t SoLKalTrackSystem::GetNDF(Bool_t self) const
{
   Int_t ndf    = 0;
   Int_t nsites = GetEntries();
   for (Int_t isite=1; isite<nsites; isite++) {
       SoLKalTrackSite &site = *static_cast<SoLKalTrackSite *>(At(isite));
       /*if (!site.IsVirtualSite())*/
        ndf += site.GetDimension();
   }
   if (self) ndf -= kSdim;
   return ndf;
}
//_____________________________________________________________________
void SoLKalTrackSystem::SmoothBackTo(Int_t k)
{
   TIter previous(this,kIterBackward);
   TIter cur     (this,kIterBackward);

   SoLKalTrackSite  *prePtr;
   SoLKalTrackSite  *curPtr = static_cast<SoLKalTrackSite *>(cur());
   SoLKalTrackState &cura   = curPtr->GetState(SoLKalTrackSite::kFiltered);
   SoLKalTrackState &scura  = curPtr->GetState(SoLKalTrackSite::kSmoothed);
   if (!&scura) {
      curPtr->Add(&curPtr->CreateState(cura, cura.GetCovMat(),
                                       SoLKalTrackSite::kSmoothed));
   }

   while ((curPtr = static_cast<SoLKalTrackSite *>(cur())) &&
          (prePtr = static_cast<SoLKalTrackSite *>(previous()))) {
      curPtr->Smooth(*prePtr);
      fCurSitePtr = curPtr;
      if (IndexOf(curPtr) == k) break;
   }
}
//_____________________________________________________________________
void SoLKalTrackSystem::SmoothAll()
{
   SmoothBackTo(0);
}
//_____________________________________________________________________
void SoLKalTrackSystem::InvFilter(Int_t k)
{
   using namespace std;
   //
   // Check if site k exists
   //
   SoLKalTrackSite  *curPtr = static_cast<SoLKalTrackSite *>(At(k));
   if (!curPtr) {
      cerr << "::::: ERROR in SoLKalTrackSystem::InvFilter(k=" << k << ")"  << endl
           << "  Site " << k << " nonexistent! Abort!"
           << endl;
      ::abort();
   }
   //
   // Check if site k already smoothed
   //
   if (!&curPtr->GetState(SoLKalTrackSite::kSmoothed)) SmoothBackTo(k);
   //
   // Inverse filter site k
   //
   fCurSitePtr = curPtr;
   curPtr->InvFilter();
}
//_____________________________________________________________________
inline void SoLKalTrackSystem::Add(TObject *obj)
{
   TObjArray::Add(obj);
   fCurSitePtr = static_cast<SoLKalTrackSite *>(obj);
   fNHits++;
   fNDF = fNHits*kMdim - kSdim;
}
//_____________________________________________________________________
void SoLKalTrackSystem::CheckTrackStatus()
{
  //monitor the position, momentum and direction of the track every time before 
  //propagation
  
  if (!fIsGood) return; //no need to check a bad track
  
  SoLKalTrackState &a = GetCurSite().GetCurState();
  //Double_t zpos = a.GetZ0();
  Double_t theta = atan(sqrt(pow(a(kIdxTX, 0),2) + pow(a(kIdxTY,0), 2)))*180./TMath::Pi();

  //position
  Double_t r = TMath::Sqrt( pow(a(0,0), 2) + pow(a(1,0), 2) );
  Double_t mom = fCharge/a(4,0);
  if (r >= 3.)
    {
      fIsGood = kFALSE;
    }
    else if ( mom > 100. || mom < 0. ){
      fIsGood = kFALSE;
    }
    else if (theta > 40. ){
      fIsGood = kFALSE;
    }
    else if (fNMissingHits > 1){
      fIsGood = kFALSE;
    }
}
//__________________________________________________________________________
Int_t SoLKalTrackSystem::Compare( const TObject* obj ) const
{
   // Used for sorting tracks in a TSeqCollection (e.g. TList, TClonesArray).
   //a track is considered better than the other if it has more hits 
   //and has a smaller chi2 per degree of freedom
   assert( dynamic_cast<const SoLKalTrackSystem*>(obj) );
   const SoLKalTrackSystem* rhs = static_cast<const SoLKalTrackSystem*>(obj);

  if (fAngleFlag < rhs->GetAngleFlag()){
    return -1;
  }else if (fAngleFlag > rhs->GetAngleFlag()){
    return 1;
  }else{
    if (fCharge < rhs->GetCharge()){
      return -1;
    }else if (fCharge > rhs->GetCharge()){
      return 1;
    }else{
      if (GetNHits() < rhs->GetNHits()){
        return 1;
      }else if (GetNHits() > rhs->GetNHits()){
        return -1;
      }else{
        if (GetChi2perNDF() > rhs->GetChi2perNDF()){
          return 1;
        }else{
          return -1;
        }
      }
    }
  }

}
//_____________________________________________________________________________














