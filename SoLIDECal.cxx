//c++
#include <cassert>
#include <sstream>
#include <iostream>

//ROOT
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
//SoLIDTracking
#include "SoLIDECal.h"
#include "SoLIDTrackerSystem.h"


//________________________________________________________________________________________
SoLIDECal::SoLIDECal( const char* name, const char* description,
                                 THaDetectorBase* parent)
:THaSubDetector(name,description,parent),
fIsLAECTriggered(kFALSE), fIsFAECTriggered(kFALSE), fLAECEdpCut(0.), fFAECEdpCut(0.)
{}
//________________________________________________________________________________________
SoLIDECal::~SoLIDECal()
{}
//________________________________________________________________________________________
void SoLIDECal::Clear(Option_t* opt)
{
  fNLAECHits = 0;
  fNFAECHits = 0;
  fIsLAECTriggered = kFALSE;
  fIsFAECTriggered = kFALSE;
  fCaloHits.clear();
}
//________________________________________________________________________________________
Int_t SoLIDECal::Decode(const THaEvData& evdata)
{
  vector<Float_t> laecXPos;
  vector<Float_t> laecYPos;
  vector<Float_t> laecEdp;
  vector<Float_t> faecXPos;
  vector<Float_t> faecYPos;
  vector<Float_t> faecEdp;
  
  if( GetNLAECHits() == 0 && GetNFAECHits() == 0) {
  

    for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
      // Decode data
      THaDetMap::Module* d = fDetMap->GetModule(imod);
      // Read active channels of this module
      Int_t nchan = evdata.GetNumChan( d->crate, d->slot );
      for( Int_t ichan = 0; ichan < nchan; ++ichan ) {
        Int_t chan = evdata.GetNextChan( d->crate, d->slot, ichan );
        if( chan < d->lo or chan > d->hi ) continue; // not part of this detector


        Int_t nhit = evdata.GetNumHits( d->crate, d->slot, chan );
        for( Int_t ihit = 0; ihit < nhit; ++ihit ) {
          // The hit's data and raw data words hold the x and y coordinates,
          // respectively
          union FloatIntUnion {
            Float_t f;
            Int_t   i;
          } datx, daty;
          datx.i = evdata.GetData( d->crate, d->slot, chan, ihit );
          daty.i = evdata.GetRawData( d->crate, d->slot, chan, ihit );
          
          
          switch(imod){
            case kLAECPos:
              SmearPosition(&datx.f, &daty.f);
              laecXPos.push_back(datx.f);
              laecYPos.push_back(daty.f);
              break;
            case kLAECEdp:
              SmearEnergy(&datx.f);
              if (datx.f > fLAECEdpCut) fIsLAECTriggered = kTRUE;
              laecEdp.push_back(datx.f);
              break;
            case kFAECPos:
              SmearPosition(&datx.f, &daty.f, kFAECPos);
              faecXPos.push_back(datx.f);
              faecYPos.push_back(daty.f);
              break;
            case kFAECEdp:
              SmearEnergy(&datx.f);
              if (datx.f > fFAECEdpCut) fIsFAECTriggered = kTRUE;
              faecEdp.push_back(datx.f);
              break;
            default:
              break;
          }
          
        }
      }
    } 
  }
  assert(laecXPos.size() == laecEdp.size() && faecXPos.size() == faecEdp.size());
  fNLAECHits = laecXPos.size();
  fNFAECHits = faecXPos.size();
  
  for (unsigned int i=0; i<fNLAECHits; i++) fCaloHits.push_back(SoLIDCaloHit( laecXPos.at(i), laecYPos.at(i), 
								     kLAEC,  laecEdp.at(i)));

  for (unsigned int i=0; i<fNFAECHits; i++) fCaloHits.push_back(SoLIDCaloHit( faecXPos.at(i), faecYPos.at(i), 
								     kFAEC,  faecEdp.at(i)));
  
  return kOK;
}
//________________________________________________________________________________________
THaAnalysisObject::EStatus SoLIDECal::Init( const TDatime& date )
{
  EStatus status = THaAnalysisObject::Init(date);
  if( status ){
  return fStatus = status;
  }
  return fStatus = kOK;
}
//_________________________________________________________________________________________
void SoLIDECal::Print( Option_t* /*opt*/ ) const
{

}
//__________________________________________________________________________________________
void SoLIDECal::PrintDataBase() const
{

}
//__________________________________________________________________________________________
Int_t SoLIDECal::Begin( THaRunBase* /*r*/ )
{
  return 0;
}
//_________________________________________________________________________________________
Int_t SoLIDECal::End( THaRunBase* /*r*/ )
{
  return 0;
}
//________________________________________________________________________________________
Int_t SoLIDECal::ReadDatabase( const TDatime& date )
{
  static const char* const here = "SoLIDECal::ReadDatabase";
  fIsInit = kFALSE;
  FILE* file = OpenFile( date );
  if( !file ) return kFileError;
  
  Int_t status = -1;
  fLAECZ = 0;
  fFAECZ = 0;
  fPosReso = -1.;//m
  fEReso = -1.;
  fMRPCPitchWidth = -1.;
  fMRPCNSectors = -1;
  fMRPCPhiCover = 0.;
  vector<Int_t>* laec_detmap_pos = 0;
  vector<Int_t>* laec_detmap_edp = 0;
  vector<Int_t>* faec_detmap_pos = 0;
  vector<Int_t>* faec_detmap_edp = 0;
  try{
   laec_detmap_pos = new vector<Int_t>;
   laec_detmap_edp = new vector<Int_t>;
   faec_detmap_pos = new vector<Int_t>;
   faec_detmap_edp = new vector<Int_t>;
   const DBRequest request[] = {
         { "laec_detmap_pos",      laec_detmap_pos,   kIntV },
         { "laec_detmap_edp",      laec_detmap_edp,   kIntV },
         { "faec_detmap_pos",      faec_detmap_pos,   kIntV },
         { "faec_detmap_edp",      faec_detmap_edp,   kIntV },
         { "laec_z",               &fLAECZ,           kDouble,  0, 1 },
         { "faec_z",               &fFAECZ,           kDouble,  0, 1 },
         { "ec_pos_reso",          &fPosReso,         kDouble,  0, 1 },
         { "ec_energy_reso",       &fEReso,           kDouble,  0, 1 },
	 { "mrpc_pitch_width",     &fMRPCPitchWidth,  kDouble,  0, 1 },
         { "mrpc_n_sectors",       &fMRPCNSectors,    kInt, 1, 1},
         { "mrpc_phi_reso",        &fMRPCPhiReso,     kDouble,  0, 1 },
	 { "mrpc_rmin",            &fMRPCRmin,        kDouble,  0, 1 },
         { 0 }
       };
   status = LoadDB( file, date, request, fPrefix );
   assert(fPosReso >= 0 && fEReso >= 0);
    if (status == kOK){
      if( FillDetMap( *laec_detmap_pos, THaDetMap::kDoNotClear, here ) <= 0 ) status = kInitError;
      if( FillDetMap( *laec_detmap_edp, THaDetMap::kDoNotClear, here ) <= 0 ) status = kInitError;
      if( FillDetMap( *faec_detmap_pos, THaDetMap::kDoNotClear, here ) <= 0 ) status = kInitError;
      if( FillDetMap( *faec_detmap_edp, THaDetMap::kDoNotClear, here ) <= 0 ) status = kInitError;
    }
    delete laec_detmap_pos;
    delete laec_detmap_edp;
    delete faec_detmap_pos;
    delete faec_detmap_edp;
  }catch(...) {
      delete laec_detmap_pos;
      delete laec_detmap_edp;
      delete faec_detmap_pos;
      delete faec_detmap_edp;
      fclose(file);
      throw;
  }
  
  fclose(file);
  if( status != kOK ) return status;
  
  for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
    THaDetMap::Module* d = fDetMap->GetModule(imod);
    assert( dynamic_cast<SoLIDTrackerSystem*>(GetMainDetector()) );
    SoLIDTrackerSystem *thisSystem = dynamic_cast<SoLIDTrackerSystem*>( GetMainDetector() );
    thisSystem->LoadDAQmodel(d);
    thisSystem->LoadDAQresolution(d);
    //only ADC for now
    d->MakeADC();
    UInt_t nchan = thisSystem->GetDAQnchan(d);
    if( d->hi >= nchan ) {
      Error( Here(here), "Detector map channel out of range for module "
          "cr/sl/lo/hi = %u/%u/%u/%u. Must be < %u. Fix database.",
          d->crate, d->slot, d->lo, d->hi, nchan );
      return kInitError;
    }
  }

  fMRPCPhiCover = 2.*TMath::Pi() / fMRPCNSectors;
  
  fIsInit = kTRUE;
  return kOK;
}
//___________________________________________________________________________________________________
Int_t SoLIDECal::DefineVariables( EMode mode )
{
  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );
  
  return kOK;
}
//___________________________________________________________________________________________________
void SoLIDECal::SmearPosition(Float_t *x, Float_t *y, Int_t mode)
{
#ifdef SIDIS
  //only in SIDIS forward angle, whre ec hit position will be replaced by hit on MRPC
  if (mode == kFAECPos){
    //if FAEC, use hit on MRPC to replace hit on EC
    
    //which MRPC sector
    //assume the first sector is centered at 0 deg
    double phi = atan2(*y, *x) + fMRPCPhiCover/2.;
    phi = TVector2::Phi_0_2pi(phi);
    int sector = (int)(phi / fMRPCPhiCover);
    assert(sector < fMRPCNSectors);
    
    double phiCenter = sector*fMRPCPhiCover;
    
    //rotate to the frame that with x-axis parallel to the symmetric axis of the sector

    double tmpx = cos(-1.*phiCenter)*(*x) + -1.*sin(-1.*phiCenter)*(*y);
    double tmpy = sin(-1.*phiCenter)*(*x) +     cos(-1.*phiCenter)*(*y);

    assert(tmpx > fMRPCRmin - 0.01); //should not happen

    int ipitch = (tmpx - fMRPCRmin)/fMRPCPitchWidth;
    tmpx  = fMRPCRmin + (ipitch + 0.5)*fMRPCPitchWidth; //using the pitch center
    tmpy += gRandom->Gaus(0., fMRPCPhiReso);

    //rotate back to the lab frame
    *x = cos(phiCenter)*tmpx + -1.*sin(phiCenter)*tmpy;
    *y = sin(phiCenter)*tmpx +     cos(phiCenter)*tmpy;
    
    return; //no need to continue
  }
#endif
  
  *x += gRandom->Gaus(0, fPosReso);
  *y += gRandom->Gaus(0, fPosReso);
  
}
//___________________________________________________________________________________________________
inline void SoLIDECal::SmearEnergy(Float_t *energy)
{
  *energy += gRandom->Gaus(0, fEReso/TMath::Sqrt(*energy));
}
























