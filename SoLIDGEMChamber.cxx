//c++
#include <cassert>
#include <iostream>
#include <sstream>
//ROOT
#include "TMath.h" 
#include "TVector2.h"
//SoLID Tracking
#include "SoLIDGEMChamber.h" 
#include "SoLIDGEMTracker.h"
#include "SoLIDTrackerSystem.h"
#include "SoLIDUtility.h"
ClassImp(SoLIDGEMChamber)

using namespace std;
SoLIDGEMChamber::SoLIDGEMChamber(Int_t ichamber, const char* name, const char* description,
                                 THaDetectorBase* parent)
  : THaSubDetector(name,description,parent), fChamberID(ichamber)
{
  static const char* const here = "SoLIDGEMChamber";
  assert( name && parent );
  fParentTrackerID = dynamic_cast<SoLIDGEMTracker*>(GetParent())->GetTrackerID();
  assert( dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector()) );
  try {
#ifdef MCDATA
    if( dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->TestBit(SoLIDTrackerSystem::kMCData) )
      fHits = new TClonesArray("SoLIDMCGEMHit", 4096);
    else
#endif
      fHits = new TClonesArray("SoLIDGEMHit", 4096);
  }
  catch( std::bad_alloc ) {
    Error( Here(here), "Out of memory allocating hit array for readout "
           "plane %s. Call expert.", name );
    MakeZombie();
    return;
  }

}
//_________________________________________________________________________________________
SoLIDGEMChamber::~SoLIDGEMChamber()
{
  if( fIsSetup )
    RemoveVariables();
  DeleteContainer(fGEMReadOut);

  delete fHits;
}
//_________________________________________________________________________________________
void SoLIDGEMChamber::Clear( Option_t* opt)
{
  if( !opt or *opt != 'I' ) {
    for (Int_t i=0; i<fNReadOut; ++i){
      fGEMReadOut[i]->Clear(opt);
    }
  }
  fUOccupancy = 0;
  fVOccupancy = 0;
  fUHitOccupancy = 0;
  fVHitOccupancy = 0;
  fHits->Clear(opt);
}
//_________________________________________________________________________________________
Int_t SoLIDGEMChamber::Decode( const THaEvData& evdata)
{

  static const char* const here = "SoLIDGEMChamber::Decode";
  //decode raw data from the two readout plane the chamber has
  for (Int_t i=0; i<fNReadOut; i++){
    fGEMReadOut[i]->Decode(evdata);
  }

  for (Int_t i=0; i<fNReadOut; i++){
    //skip the very noisy event, will slow down the spped a lot. Not worth it.
    if ( fGEMReadOut[i]->IsNoisyEvent() == 1 ){
     cout<<"found a noisy event"<<endl;
     Warning(Here(here), "Noisy event on chamber %d tracker %d",
     fChamberID, dynamic_cast<SoLIDGEMTracker*>(GetParent())->GetTrackerID());
     return -1;
    }
  }
  if (!fDo3DAmCorr) return 0; //not going to do 3d amplitude matching here, if requested

  assert(fNReadOut == 2);//again, require 2D readout

  TSeqCollection* uhits = fGEMReadOut[0]->GetHits();
  TSeqCollection* vhits = fGEMReadOut[1]->GetHits();
  assert(uhits && vhits);
  Int_t nHit = ProcessRawHits(uhits, vhits);

  fUOccupancy = fGEMReadOut[0]->GetOccupancy();
  fVOccupancy = fGEMReadOut[1]->GetOccupancy();
  fUHitOccupancy = fGEMReadOut[0]->GetHitOccupancy();
  fVHitOccupancy = fGEMReadOut[1]->GetHitOccupancy();

  return 1;
}
//_________________________________________________________________________________________
THaAnalysisObject::EStatus SoLIDGEMChamber::Init( const TDatime& date )
{
  EStatus status = THaAnalysisObject::Init(date);
  if (status == kOK){
    //fGEMReadOut = new SoLIDGEMReadOut *[fNReadOut];
    for (Int_t i=0; i<fNReadOut; i++){
      stringstream sn, sd;
      THaDetectorBase *p = GetParent();
      sn <<i;
      sd << "Readout " << i << " of GEM Chamber "<<fChamberID<<
      " on GEM tracker"<<dynamic_cast<SoLIDGEMTracker*>(p)->GetTrackerID();
      fTrackerZ = dynamic_cast<SoLIDGEMTracker*>(p)->GetZ();
      SoLIDGEMReadOut* theReadOut = new SoLIDGEMReadOut(i, sn.str().c_str(), 
                                                      sd.str().c_str(), this);
      fGEMReadOut.push_back(theReadOut);
      status = fGEMReadOut[i]->Init(date);
      if (status) break;
    }
  }
  
  if( status ){
  return fStatus = status;
  }
  
  return fStatus = kOK;
}
//_________________________________________________________________________________________
void SoLIDGEMChamber::Print( Option_t* opt ) const
{

}
//_________________________________________________________________________________________
void SoLIDGEMChamber::PrintDataBase(Int_t level) const
{
  if (level>=0 && level<2){
    if (level == 0){
      Int_t grand_parent_systemID = 
      dynamic_cast<SoLIDTrackerSystem*>((dynamic_cast<SoLIDGEMTracker*>(GetParent()))->GetParent())->GetSystemID();
      Int_t parent_trackerID = dynamic_cast<SoLIDGEMTracker*>(GetParent())->GetTrackerID();
      string out_prefix = Form("solid.trackersystem.%d.%d.", grand_parent_systemID, parent_trackerID);
      cout<<"******parameter from database for chamber "<<fChamberID<<" in tracker "
      <<parent_trackerID<<" in tracker system "<<grand_parent_systemID<<"******"<<endl;
      cout<<out_prefix<<fChamberID<<".do_3d_amcorr = "<<fDo3DAmCorr<<endl;
      cout<<out_prefix<<fChamberID<<".3d_amcorr_cut = "<<f3DAmCorrCut<<endl;
      cout<<out_prefix<<fChamberID<<".nreadout = "<<fNReadOut<<endl;
      cout<<out_prefix<<fChamberID<<".rmin = "<<fRMin<<endl;
      cout<<out_prefix<<fChamberID<<".rmax = "<<fRMax<<endl;
      cout<<out_prefix<<fChamberID<<".dz = "<<fDz<<endl;
      cout<<out_prefix<<fChamberID<<".phi = "<<fPhi<<endl;
      cout<<out_prefix<<fChamberID<<".phi_cover = "<<fPhiCover<<endl;
      cout<<out_prefix<<fChamberID<<".phi_offset = "<<fPhiOffset<<endl;
      cout<<out_prefix<<fChamberID<<".reference = "<<fReference.X()<<" "<<fReference.Y()<<endl;
      cout<<"***************************************************************************"<<endl;
    }else if(level > 0){
      level--;
      for (Int_t i=0; i<fNReadOut; i++){
        fGEMReadOut[i]->PrintDataBase(level);
      }
    }
  }
} 
//_________________________________________________________________________________________
Int_t SoLIDGEMChamber::Begin( THaRunBase* /*r*/ )
{
  return 0;
}
//_________________________________________________________________________________________
Int_t SoLIDGEMChamber::End( THaRunBase* /*r*/ )
{
  return 0;
}
//_________________________________________________________________________________________
Int_t SoLIDGEMChamber::ReadDatabase( const TDatime& date )
{
  fIsInit = kFALSE;
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;
  
   Int_t err = ReadGeometry( file, date );
  if( err ) {
    fclose(file);
    return err;
  }
  
  fNReadOut = -1;
  fDo3DAmCorr = -1;
  f3DAmCorrCut = -1;
  const DBRequest request[] = {
        { "do_3d_amcorr",   &fDo3DAmCorr,      kInt,     0, 1 },
        { "3d_amcorr_cut",  &f3DAmCorrCut,     kDouble,  0, 1 },
        { "nreadout",       &fNReadOut,        kInt,     0, 1 },
        { 0 }
      };

  Int_t status = LoadDB( file, date, request, fPrefix );

  assert(fDo3DAmCorr >=0 && f3DAmCorrCut >=0 && fNReadOut >0 );

  if (status != kOK) { return status; }

  fIsInit = kTRUE;
  return fStatus = kOK;
}
//__________________________________________________________________________________________
Int_t SoLIDGEMChamber::ReadGeometry( FILE* file, const TDatime& date,
                                Bool_t /* required */ )
{

  static const char* const here = "SoLIDGEMChamber::ReadGeometry";
  
  Int_t status = -1;
  fRMin = -1.;
  fRMax = -1.;
  vector<Double_t>* reference = 0;
  try{
    reference = new vector<Double_t>;
    const DBRequest request[] = {
          { "rmin",           &fRMin,            kDouble,  0, 1 },
          { "rmax",           &fRMax,            kDouble,  0, 1 },
          { "dz",             &fDz,              kDouble,  0, 1 },
          { "phi",            &fPhi,             kDouble,  0, 1 },
          { "phi_cover",      &fPhiCover,        kDouble,  0, 1 },
          { "phi_offset",     &fPhiOffset,       kDouble,  0, 1 },
          { "reference",       reference,        kDoubleV},
          { 0 }
      };
    status = LoadDB( file, date, request, fPrefix );
    assert(reference->size() == 2); //a reference is a point in x-y, z is already know from tracker
    //fReference.Set(reference->at(0), reference->at(1));
    delete reference; //it has served its purpose, rest in peace
  }catch(...) {
      delete reference;
      fclose(file);
      throw;
  }
  
   if( fRMin < 0 or fRMax < 0 ) {
    Error( Here(here), "rmin and rmax must be positive. Read %.1lf and %.1lf. "
           "Fix database.", fRMin, fRMax );
    return kInitError;
  }
  if( fRMin >= fRMax ) {
    Error( Here(here), "rmin = %.1lf must be less than rmax = %.1lf. "
           "Fix database.", fRMin, fRMax );
    return kInitError;
  }
  // Limit the opening angle to simplify subsequent calculations
  if( fPhiCover <= 0 or fPhiCover >= 90. ) {
    Error( Here(here), "Illegal value for dphi = %.1lf. Must be > 0 and <= "
           "90 degrees. Fix database.", fPhiCover );
    return kInitError;
  }

  Double_t max_off = 90.0 - 0.5*fPhiCover;
  if( TMath::Abs(fPhiOffset) >= max_off ) {
    Error( Here(here), "Illegal value for phioff = %.1lf. Must be between "
           "%.1lf and %.1lf degrees. Fix database.",
           fPhiOffset, -max_off, max_off );
    return kInitError;
  }

  //change all the angles from deg to rad and keep them in -pi to pi for consistancy

  fPhi = TVector2::Phi_mpi_pi( fPhi*TMath::DegToRad() );
  fPhiCover = TVector2::Phi_mpi_pi( fPhiCover*TMath::DegToRad() );
  fPhiOffset = TVector2::Phi_mpi_pi( fPhiOffset*TMath::DegToRad() );
  //phi offset is something we want to minus, if notice there is a constant shift in phi, this is
  //likely the cause
  fPhiInLab = fPhi - fPhiOffset + dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->GetPhi();
  fPhiInLab = TVector2::Phi_mpi_pi(fPhiInLab);

  Double_t phi2 = 0.5*fPhiCover;
  fPhiMin = -phi2 - fPhiOffset;
  fPhiMax =  phi2 - fPhiOffset;
  assert( fPhiMin < fPhiMax );
  assert( fPhiMin > -TMath::PiOver2() );
  assert( fPhiMax < TMath::PiOver2() );

  // Define the origin of the plane in the same way as in libsolgem: It is
  // the center of the bounding box of the ring segment, calculated WITHOUT
  // any phi offset rotation.
  // Notice that this is no longer true if the chamber is shifted inward in
  // the radial direction. Phi angle coverage is ill-defined know. But, it can
  // still work in this way if we define a reference point, at which the lines
  // that are parallel to the longer edge of the chamber intersect.
  // this info is stored in fReference TVector2. One will need to define this
  // point in the data base, since I don't know how much it will be shifted inward here

  Double_t xmin = fRMin * TMath::Cos(phi2), xmax = fRMax;
  Double_t r_shift = TMath::Sqrt(fReference.X()*fReference.X() +
                                 fReference.Y()*fReference.Y() );
  //TODO: test this later on when there is really a reference, not is (0,0)
  fOrigin.SetXYZ( 0.5*(xmin+xmax) - r_shift, 0.0,
                  dynamic_cast<SoLIDGEMTracker*>(GetParent())->GetZ() + fDz);
  fChamberCenter = fOrigin;

  if( fPhiOffset != 0.0 ) {
  //TODO:think if we want to do all the rotations here once and for all, or we break it down
  //at different level?

    fOrigin.RotateZ(dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->GetPhi() +
                    fPhi-fPhiOffset);
  }

  //just to compare the raw decode result in TreeSearch, we set all the chamber origin relative
  //to the first chamber in the first GEM tracker plane, will delete this later
  //fOrigin.SetXYZ(fOrigin.X(), fOrigin.Y(), fOrigin.Z());
  return kOK;


}
//________________________________________________________________________________________
Int_t SoLIDGEMChamber::DefineVariables(EMode mode)
{
  // initialize global variables
  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  // Register variables in global list
  Int_t ret;
#ifdef MCDATA
  if( !dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->TestBit(SoLIDTrackerSystem::kMCData) ) {
#endif
     // Non-Monte Carlo hit data
    RVarDef nonmcvars[] = {
      { "hit2D.x",       "2D hit x coordinate",          "fHits.SoLIDGEMHit.GetX()"         },
      { "hit2D.y",       "2D hit y coordinate",          "fHits.SoLIDGEMHit.GetY()"         },
      { "hit2D.z",       "2D hit z coordinate",          "fHits.SoLIDGEMHit.GetZ()"         },
      { "hit2D.r",       "2D hit r coordinate",          "fHits.SoLIDGEMHit.GetR()"         },
      { "hit2D.phi",     "2D hit phi coordinate",        "fHits.SoLIDGEMHit.GetPhi()"       },
      { "hit2D.chamber", "2D hit chamber ID",            "fHits.SoLIDGEMHit.GetChamberID()" },
      { 0 }
    };
    ret = DefineVarsFromList( nonmcvars, mode );
  }else{
    //Monte-Carlo hit data
    RVarDef mcvars[] = {
      { "hit2D.x",       "2D hit x coordinate",            "fHits.SoLIDMCGEMHit.GetX()"        },
      { "hit2D.y",       "2D hit y coordinate",            "fHits.SoLIDMCGEMHit.GetY()"        },
      { "hit2D.z",       "2D hit z coordinate",            "fHits.SoLIDMCGEMHit.GetZ()"        },
      { "hit2D.r",       "2D hit r coordinate",            "fHits.SoLIDMCGEMHit.GetR()"        },
      { "hit2D.qu",      "2D hit charge deposition on u",  "fHits.SoLIDMCGEMHit.GetQU()"       },
      { "hit2D.qv",      "2D hit charge deposition on v",  "fHits.SoLIDMCGEMHit.GetQV()"       },
      { "hit2D.phi",     "2D hit phi coordinate",          "fHits.SoLIDMCGEMHit.GetPhi()"      },
      { "hit2D.signal",  "if this hit is signal",          "fHits.SoLIDMCGEMHit.IsSignalHit()" },
      { "hit2D.chamber", "2D hit chamber ID",              "fHits.SoLIDGEMHit.GetChamberID()"  },
      { 0 }
    };
    ret = DefineVarsFromList( mcvars, mode );
  }
  return ret;
}
//_________________________________________________________________________________________
Int_t SoLIDGEMChamber::ProcessRawHits(TSeqCollection* uhits, TSeqCollection* vhits)
{
  assert(fNReadOut == 2); //if not 2D readout, we need to talk
  UInt_t nHit = 0;
  Bool_t mc_data = dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector())->TestBit(SoLIDTrackerSystem::kMCData);
  TIterator* uit = uhits->MakeIterator();
  TIterator* vit = vhits->MakeIterator();

  SoLIDRawHit* auhit = 0;
  SoLIDRawHit* avhit = 0;
  while( (auhit = static_cast<SoLIDRawHit*>(uit->Next()) ) ) {
    if (!auhit->GetStatus()) continue;
    vit -> Reset();
    while( (avhit = static_cast<SoLIDRawHit*>(vit->Next()) ) ){
      if (!avhit->GetStatus()) continue;
      Double_t r = 0;
      Double_t phi = 0;
      UVtoCylinCoor(auhit->GetPos(), avhit->GetPos(), &r, &phi);
      if (Contains( &r, &phi ) && CheckChargeAsymmetry( dynamic_cast<SoLIDRawHit*>(auhit)->GetADCsum(),
          dynamic_cast<SoLIDRawHit*>(avhit)->GetADCsum()) ){
          RotateToLab(&phi);
          if (!mc_data){
            new ( (*fHits)[nHit++]) SoLIDGEMHit( fChamberID, fParentTrackerID, r, phi, GetZ(), auhit, avhit); 
          }
#ifdef MCDATA
          else{
            new ( (*fHits)[nHit++]) SoLIDMCGEMHit( fChamberID, fParentTrackerID, r, phi, GetZ(), auhit, avhit);
          }
#endif
      }
    }
  }
  //sort the hit array now from small radius to large radius
  fHits->Sort();
  return nHit;
}
//____________________________________________________________________________________________
inline void SoLIDGEMChamber::UVtoCylinCoor(Double_t upos, Double_t vpos, Double_t* r, Double_t* phi)
{

  assert(fNReadOut == 2);
  Double_t su = fGEMReadOut[0]->GetSinStripAngle();
  Double_t sv = fGEMReadOut[1]->GetSinStripAngle();
  Double_t cu = fGEMReadOut[0]->GetCosStripAngle();
  Double_t cv = fGEMReadOut[1]->GetCosStripAngle();
  Double_t ba = sv*cu-su*cv;

  Double_t x = (upos*sv - vpos*su)/ba;
  Double_t y = (vpos*cu - upos*cv)/ba;

  *r = TMath::Sqrt(x*x + y*y);
  *phi = TMath::ATan2(y, x); //this is in a frame where center of the symmetric axis of the chamber is
                             //the same as x axis
}
//____________________________________________________________________________________________
inline Bool_t SoLIDGEMChamber::Contains(Double_t* r, Double_t *phi)
{
  if ( (*r <= fRMax && *r >= fRMin ) && ( *phi <= fPhiCover/2.
        && *phi >= -fPhiCover/2. ) ){ return kTRUE; }
  else { return kFALSE; }
}
//____________________________________________________________________________________________
inline Bool_t SoLIDGEMChamber::CheckChargeAsymmetry(Double_t qu, Double_t qv)
{
  if ( TMath::Abs((qu-qv)/(qu+qv)) < f3DAmCorrCut ) return kTRUE;
  else return kFALSE;
}
//____________________________________________________________________________________________
inline void SoLIDGEMChamber::RotateToLab(Double_t *phi)
{
  *phi += fPhiInLab;
  *phi = TVector2::Phi_mpi_pi(*phi);
}














