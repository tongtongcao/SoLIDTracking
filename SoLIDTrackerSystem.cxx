//*-- Author :    Ole Hansen, Jefferson Lab   06-Feb-2012

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SoLID::GEMTracker                                                         //
//                                                                           //
// This is class GEMTracker in namespace SoLID. It inherits from class       //
// GEMTracker in namespace TreeSearch. Although confusing, this is fine.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//c++
#include <sstream>
#include <cassert>
#include <map>
//root
#include "TClonesArray.h"
#include "TMath.h"
//Hall A analyzer
#include "THaAnalysisObject.h"
#ifndef MCDATA
#include "THaEvData.h"
#else
#include "SimDecoder.h"
#endif
//SoLID tracking
#include "SoLIDTrackerSystem.h"
#include "SoLIDSpectrometer.h"
#include "SoLIDUtility.h"
#include "SoLIDGEMHit.h"
#include "SoLIDTrack.h"
using namespace std;
using namespace Podd;

#define ALL(c) (c).begin(), (c).end()

typedef string::size_type ssiz_t;
typedef vector<vector<Int_t> >::iterator vviter_t;

//typedef vector<Plane*>::size_type vrsiz_t;
//_____________________________________________________________________________
SoLIDTrackerSystem::SoLIDTrackerSystem( const char* name, const char* desc, THaApparatus* app)
  :THaTrackingDetector(name,desc,app), fSystemID(-1), 
   fPhi(0), fTracks(0), fCrateMap(0), fTrackFinder(0)
#ifdef MCDATA
  , fMCDecoder(0), fChecked(false)
#endif
{
#ifdef MCDATA
  fTracks = new TClonesArray("SoLIDMCTrack", 10);
#else
  fTracks = new TClonesArray("SoLIDTrack", 10);
#endif

  fFieldMap = SoLIDFieldMap::GetInstance();
}

//_____________________________________________________________________________
SoLIDTrackerSystem::~SoLIDTrackerSystem()
{
  if( fIsSetup )
    RemoveVariables();
  DeleteContainer(fGEMTracker);
  delete fTracks;
  delete fECal;
  delete fTrackFinder;
  
}
//_____________________________________________________________________________
Int_t SoLIDTrackerSystem::ReadDatabase(const TDatime& date)
{
   // Read generic Tracker keys from database

  static const char* const here = "SoLIDTrackerSystem::ReadDatabase";
  fIsInit = kFALSE;
  // Delete existing configuration (in case we are re-initializing)
  //DeleteContainer( fProj );
  //DeleteContainer( fPlanes );
  //fCalibPlanes.clear();

  FILE* file = OpenFile( date );
  if( !file ) return kFileError;

  
  fOrigin.SetXYZ(0,0,0);
  
  Int_t err = ReadGeometry( file, date );
  if( err ) {
    fclose(file);
    return err;
  }

  // Putting this container on the stack may cause strange stack overflows!
  vector<vector<Int_t> > *cmap = new vector<vector<Int_t> >;

#ifdef MCDATA
  Int_t mc_data = 0;
#endif
  fNTracker = -1;
  fChi2Cut  = -1;
  fNMaxMissHit = -1;
  Int_t do_rawdecode = -1, do_coarsetrack = -1, do_finetrack = -1, do_chi2 = -1;
  assert( GetCrateMapDBcols() >= 5 );
  DBRequest request[] = {
    { "cratemap",          cmap,               kIntM,   GetCrateMapDBcols() },
#ifdef MCDATA
    { "MCdata",            &mc_data,           kInt,    0, 1 },
#endif
    { "do_rawdecode",      &do_rawdecode,      kInt,    0, 1 },
    { "do_coarsetrack",    &do_coarsetrack,    kInt,    0, 1 },
    { "do_finetrack",      &do_finetrack,      kInt,    0, 1 },
    { "do_chi2",           &do_chi2,           kInt,    0, 1 },
    { "chi2_cut",          &fChi2Cut,          kDouble, 0, 1 },
    { "max_miss_hit",      &fNMaxMissHit,      kInt,    0, 1 },
    { "ntracker",          &fNTracker,         kInt,    0, 1 },
    { 0 }
  };
  
  Int_t status = kInitError;
  err = LoadDB( file, date, request, fPrefix );
  assert( fNTracker > 0 && fChi2Cut > 0 && fNMaxMissHit > 0);
  fclose(file);
  if( !err ) {
    if( cmap->empty() ) {
      Error(Here(here), "No cratemap defined. Set \"cratemap\" in database.");
    } else {
      // Build the list of crate map elements
      for( vviter_t it = cmap->begin(); it != cmap->end(); ++it ) {
	vector<int>& row = *it;
	assert( row.size() == GetCrateMapDBcols() ); // failure indicates bug in LoadDB
	for( Int_t slot = row[1]; slot <= row[2]; ++slot ) {
	  DAQmodule* m = 0;
	  if( GetCrateMapDBcols() < 6 ){
	    m = new DAQmodule( row[0], slot, row[3], row[4] );
	  }
	  else{
	    m = new DAQmodule( row[0], slot, row[3], row[4], row[5] );
	  }
	  DAQmodule* found = static_cast<DAQmodule*>(fCrateMap->FindObject(m));
	  if( found ) {
	    m->Copy(*found);  // Later entries override earlier ones
	    delete m;
	  }
	  else{
	    fCrateMap->Add(m);
	  }
	}
      }
      status = kOK;
    }
  }
  delete cmap; cmap = 0;
  if( status != kOK )
    return status;

  // Set common analysis control flags
#ifdef MCDATA
  SetBit( kMCData,        mc_data );
#endif
  SetBit( kDoRawDecode,   do_rawdecode );
  SetBit( kDoCoarse,      do_coarsetrack );
  SetBit( kDoFine,        do_coarsetrack && do_finetrack );
  SetBit( kDoChi2,        do_chi2 );

  cout << endl;
  if( fDebug > 0 ) {
#ifdef MCDATA
    Info( Here(here), "Tracker flags mcdata/rawdecode/coarse/fine/chi2 = "
	  "%d/%d/%d/%d/%d", TestBit(kMCData), TestBit(kDoRawDecode),
	  TestBit(kDoCoarse), TestBit(kDoFine), TestBit(kDoChi2) );
#else
    Info( Here(here), "Tracker flags rawdecode/coarse/fine/chi2 = "
	  "%d/%d/%d/%d",  TestBit(kDoRawDecode),
	  TestBit(kDoCoarse), TestBit(kDoFine), TestBit(kDoChi2) );
#endif
  }
 
  fIsInit = kTRUE;
  
  return kOK;
}
//_____________________________________________________________________________
void SoLIDTrackerSystem::Clear(Option_t* opt)
{
  // Clear event-by-event data, including those of the planes and projections
  THaTrackingDetector::Clear(opt);
  
  if( !opt or *opt != 'I' ) {
    for( Int_t i = 0; i < fNTracker; ++i ){
      fGEMTracker[i]->Clear(opt);
    }
    fECal->Clear(opt);
    fTrackFinder->Clear(opt);
  }
  fTracks->Clear(opt);
}
//_____________________________________________________________________________
Int_t SoLIDTrackerSystem::Decode( const THaEvData& evdata)
{
#ifdef MCDATA
const char* const here = "SoLIDTrackerSystem::Decode";
  if( !fChecked ) {
    if( TestBit(kMCData) ) {
      fMCDecoder = dynamic_cast<const Podd::SimDecoder*>(&evdata);
      if( fMCDecoder == 0 ) {
        Error( Here(here), "MCdata flag set, but decoder is not a SimDecoder. "
               "Fix database or replay configuration." );
      }
    }
    fChecked = true;
  }

#endif

  for (Int_t i=0; i<fNTracker; i++){
    fGEMTracker[i]->Decode(evdata);
  }
  fECal->Decode(evdata);
  return kOK;
}
//_____________________________________________________________________________
THaAnalysisObject::EStatus SoLIDTrackerSystem::Init( const TDatime& date )
{
  assert( fCrateMap == 0 );
  fCrateMap = new THashTable(100);
  fCrateMap->SetOwner();
  EStatus status = THaTrackingDetector::Init(date);
  //all GEM trackers
  if (status == kOK){
    //fGEMTracker = new SoLIDGEMTracker *[fNTracker];
    for (Int_t i=0; i<fNTracker; i++){
      stringstream sn, sd;
      sn <<i;
      sd << "SoLID GEM tracker " << i;
      SoLIDGEMTracker *theTracker = new SoLIDGEMTracker(i, sn.str().c_str(), 
                                                      sd.str().c_str(), this);
      fGEMTracker.push_back(theTracker);
      status = fGEMTracker[i]->Init(date);
      if( status ) break;
    }
  }
  //the e-calorimeter
  if (status == kOK){
    stringstream sn, sd;
      sn <<"ecal";
      sd << "SoLID Electromagnetic Calorimeters ";
      fECal = new SoLIDECal(sn.str().c_str(), sd.str().c_str(), this);
  }
  status = fECal->Init(date);
  
  delete fCrateMap; fCrateMap = 0;
  if( status ){
  return fStatus = status;
  }

  fTrackFinder = new SoLKalTrackFinder(TestBit(kMCData), fNTracker);
  fTrackFinder->SetGEMDetector(fGEMTracker);
  fTrackFinder->SetECalDetector(fECal);
  
  return fStatus = kOK;
}
//_____________________________________________________________________________
Int_t SoLIDTrackerSystem::CoarseTrack( TClonesArray& /*tracks*/ )
{
  /*This is the main function where pattern recognition (track finder ) 
  is done (identify the hits that belong to a high energy signal track)
  1. collect all the hit array from all trackers and all chambers
  2. put the hit arrays into the pattern recognition class, what ever that is
  3. collect the output from the pattern recognition class, and fill the SoLIDTrack
     object (optional, may be done in the pattern recognition class) 
  */
  //static const char* const here = "SoLIDTrackerSystem::CoarseTrack";
  
  //c++ map for storing pointers to the hit object array
  if ( TestBit(kDoCoarse) ){  
  fGoodSignalFlag = 0;
  Int_t signalCount[6] = {0}; //Temperory for checking
  
  //check to see if there is any high energy hit on the calorimeters
  /*if ( fECal->IsLAECTriggered() ){
    map<Int_t, vector<Float_t> > * thisECHitMap = fECal->GetLAECHits();
    vector<Float_t> * ecalXPos = &(thisECHitMap->find(kECalXPos)->second);
    vector<Float_t> * ecalYPos = &(thisECHitMap->find(kECalYPos)->second);
    vector<Float_t> * ecalEdp  = &(thisECHitMap->find(kECalEdp)->second);
    assert(ecalXPos->size() == ecalYPos->size() && ecalXPos->size() == ecalEdp->size() );
    for (UInt_t j=0; j<ecalXPos->size(); j++){
      fTrackFinder->SetCaloHit(ecalXPos->at(j), ecalYPos->at(j), 0, ecalEdp->at(j));
    }
  }else if( fECal->IsFAECTriggered() ){
    map<Int_t, vector<Float_t> > * thisECHitMap = fECal->GetFAECHits();
    vector<Float_t> * ecalXPos = &(thisECHitMap->find(kECalXPos)->second);
    vector<Float_t> * ecalYPos = &(thisECHitMap->find(kECalYPos)->second);
    vector<Float_t> * ecalEdp  = &(thisECHitMap->find(kECalEdp)->second);
    assert(ecalXPos->size() == ecalYPos->size() && ecalXPos->size() == ecalEdp->size() );
    for (UInt_t j=0; j<ecalXPos->size(); j++){
      fTrackFinder->SetCaloHit(ecalXPos->at(j), ecalYPos->at(j), 1, ecalEdp->at(j));
    }
     
  }*/
  
  //getting the Beam spot, using MC info for now
#ifdef MCDATA
  if( TestBit(kMCData) ) {
    assert( dynamic_cast<MCTrack*>(fMCDecoder->GetMCTrack(0)) );
    MCTrack* trk = static_cast<MCTrack*>( fMCDecoder->GetMCTrack(0) );
    assert(trk);
    fTrackFinder->SetBPM(trk->VX(), trk->VY());
  }
#endif
  
  //this is where the actual pattern recognition done
  fTrackFinder->ProcessHits(fTracks);
  }
  return kOK;
}
//_____________________________________________________________________________
Int_t SoLIDTrackerSystem::FineTrack( TClonesArray& /*tracks*/ )
{
  if (TestBit(kDoFine)) {}
  return 1;
}
//_____________________________________________________________________________
Int_t SoLIDTrackerSystem::DefineVariables( EMode mode )
{
  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );
  
  // Register variables in global list
  Int_t ret;

  if( !TestBit(kMCData) ) {

     // Non-Monte Carlo hit data
    RVarDef nonmcvars[] = {
      { "track.0.x",       "x coordinate of the 1st hit",    "fTracks.SoLIDTrack.GetHitX0()"         },
      { "track.1.x",       "x coordinate of the 2nd hit",    "fTracks.SoLIDTrack.GetHitX1()"         },
      { "track.2.x",       "x coordinate of the 3rd hit",    "fTracks.SoLIDTrack.GetHitX2()"         },
      { "track.3.x",       "x coordinate of the 4th hit",    "fTracks.SoLIDTrack.GetHitX3()"         },
      { "track.4.x",       "x coordinate of the 5th hit",    "fTracks.SoLIDTrack.GetHitX4()"         },
      { "track.5.x",       "x coordinate of the 6th hit",    "fTracks.SoLIDTrack.GetHitX5()"         },
      { "track.0.y",       "y coordinate of the 1st hit",    "fTracks.SoLIDTrack.GetHitY0()"         },
      { "track.1.y",       "y coordinate of the 2nd hit",    "fTracks.SoLIDTrack.GetHitY1()"         },
      { "track.2.y",       "y coordinate of the 3rd hit",    "fTracks.SoLIDTrack.GetHitY2()"         },
      { "track.3.y",       "y coordinate of the 4th hit",    "fTracks.SoLIDTrack.GetHitY3()"         },
      { "track.4.y",       "y coordinate of the 5th hit",    "fTracks.SoLIDTrack.GetHitY4()"         },
      { "track.5.y",       "y coordinate of the 6th hit",    "fTracks.SoLIDTrack.GetHitY5()"         },
      { "track.0.predx",       "x coordinate of the 1st hit",    "fTracks.SoLIDTrack.GetPredHitX0()"         },
      { "track.1.predx",       "x coordinate of the 2nd hit",    "fTracks.SoLIDTrack.GetPredHitX1()"         },
      { "track.2.predx",       "x coordinate of the 3rd hit",    "fTracks.SoLIDTrack.GetPredHitX2()"         },
      { "track.3.predx",       "x coordinate of the 4th hit",    "fTracks.SoLIDTrack.GetPredHitX3()"         },
      { "track.4.predx",       "x coordinate of the 5th hit",    "fTracks.SoLIDTrack.GetPredHitX4()"         },
      { "track.5.predx",       "x coordinate of the 6th hit",    "fTracks.SoLIDTrack.GetPredHitX5()"         },
      { "track.0.predy",       "y coordinate of the 1st hit",    "fTracks.SoLIDTrack.GetPredHitY0()"         },
      { "track.1.predy",       "y coordinate of the 2nd hit",    "fTracks.SoLIDTrack.GetPredHitY1()"         },
      { "track.2.predy",       "y coordinate of the 3rd hit",    "fTracks.SoLIDTrack.GetPredHitY2()"         },
      { "track.3.predy",       "y coordinate of the 4th hit",    "fTracks.SoLIDTrack.GetPredHitY3()"         },
      { "track.4.predy",       "y coordinate of the 5th hit",    "fTracks.SoLIDTrack.GetPredHitY4()"         },
      { "track.5.predy",       "y coordinate of the 6th hit",    "fTracks.SoLIDTrack.GetPredHitY5()"         },
      { "track.0.predex",       "x coordinate of the 1st hit",    "fTracks.SoLIDTrack.GetPredHiteX0()"         },
      { "track.1.predex",       "x coordinate of the 2nd hit",    "fTracks.SoLIDTrack.GetPredHiteX1()"         },
      { "track.2.predex",       "x coordinate of the 3rd hit",    "fTracks.SoLIDTrack.GetPredHiteX2()"         },
      { "track.3.predex",       "x coordinate of the 4th hit",    "fTracks.SoLIDTrack.GetPredHiteX3()"         },
      { "track.4.predex",       "x coordinate of the 5th hit",    "fTracks.SoLIDTrack.GetPredHiteX4()"         },
      { "track.5.predex",       "x coordinate of the 6th hit",    "fTracks.SoLIDTrack.GetPredHiteX5()"         },
      { "track.0.predey",       "y coordinate of the 1st hit",    "fTracks.SoLIDTrack.GetPredHiteY0()"         },
      { "track.1.predey",       "y coordinate of the 2nd hit",    "fTracks.SoLIDTrack.GetPredHiteY1()"         },
      { "track.2.predey",       "y coordinate of the 3rd hit",    "fTracks.SoLIDTrack.GetPredHiteY2()"         },
      { "track.3.predey",       "y coordinate of the 4th hit",    "fTracks.SoLIDTrack.GetPredHiteY3()"         },
      { "track.4.predey",       "y coordinate of the 5th hit",    "fTracks.SoLIDTrack.GetPredHiteY4()"         },
      { "track.5.predey",       "y coordinate of the 6th hit",    "fTracks.SoLIDTrack.GetPredHiteY5()"         },
      { "track.0.tracker",      "tracker ID of the 1st hit",      "fTracks.SoLIDTrack.GetHitTracker0()"   },
      { "track.1.tracker",      "tracker ID of the 2nd hit",      "fTracks.SoLIDTrack.GetHitTracker1()"   },
      { "track.2.tracker",      "tracker ID of the 3rd hit",      "fTracks.SoLIDTrack.GetHitTracker2()"   },
      { "track.3.tracker",      "tracker ID of the 4th hit",      "fTracks.SoLIDTrack.GetHitTracker3()"   },
      { "track.4.tracker",      "tracker ID of the 5th hit",      "fTracks.SoLIDTrack.GetHitTracker4()"   },
      { "track.5.tracker",      "tracker ID of the 6th hit",      "fTracks.SoLIDTrack.GetHitTracker5()"   },
      { "track.ntrack",         "number of tracks for the event", "GetNTracks()"                          },
      { "track.p",              "momentum of the track",          "fTracks.SoLIDTrack.GetMomentum()"},
      { "track.theta",          "polar angle of the track",       "fTracks.SoLIDTrack.GetTheta()"},
      { "track.phi",            "azimuthal angle of the track",   "fTracks.SoLIDTrack.GetPhi()"},
      { "track.vertexz",        "vertex z of the track",          "fTracks.SoLIDTrack.GetVertexZ()"},
      { 0 }   
    };
    ret = DefineVarsFromList( nonmcvars, mode );
  }else{
#ifdef MCDATA
    //Monte-Carlo hit data
    RVarDef mcvars[] = {
      { "track.0.x",       "x coordinate of the 1st hit",    "fTracks.SoLIDMCTrack.GetHitX0()"         },
      { "track.1.x",       "x coordinate of the 2nd hit",    "fTracks.SoLIDMCTrack.GetHitX1()"         },
      { "track.2.x",       "x coordinate of the 3rd hit",    "fTracks.SoLIDMCTrack.GetHitX2()"         },
      { "track.3.x",       "x coordinate of the 4th hit",    "fTracks.SoLIDMCTrack.GetHitX3()"         },
      { "track.4.x",       "x coordinate of the 5th hit",    "fTracks.SoLIDMCTrack.GetHitX4()"         },
      { "track.5.x",       "x coordinate of the 6th hit",    "fTracks.SoLIDMCTrack.GetHitX5()"         },
      { "track.0.y",       "y coordinate of the 1st hit",    "fTracks.SoLIDMCTrack.GetHitY0()"         },
      { "track.1.y",       "y coordinate of the 2nd hit",    "fTracks.SoLIDMCTrack.GetHitY1()"         },
      { "track.2.y",       "y coordinate of the 3rd hit",    "fTracks.SoLIDMCTrack.GetHitY2()"         },
      { "track.3.y",       "y coordinate of the 4th hit",    "fTracks.SoLIDMCTrack.GetHitY3()"         },
      { "track.4.y",       "y coordinate of the 5th hit",    "fTracks.SoLIDMCTrack.GetHitY4()"         },
      { "track.5.y",       "y coordinate of the 6th hit",    "fTracks.SoLIDMCTrack.GetHitY5()"         },
      { "track.0.px",       "x coordinate of the 1st hit",    "fTracks.SoLIDMCTrack.GetHitPX0()"         },
      { "track.1.px",       "x coordinate of the 2nd hit",    "fTracks.SoLIDMCTrack.GetHitPX1()"         },
      { "track.2.px",       "x coordinate of the 3rd hit",    "fTracks.SoLIDMCTrack.GetHitPX2()"         },
      { "track.3.px",       "x coordinate of the 4th hit",    "fTracks.SoLIDMCTrack.GetHitPX3()"         },
      { "track.4.px",       "x coordinate of the 5th hit",    "fTracks.SoLIDMCTrack.GetHitPX4()"         },
      { "track.5.px",       "x coordinate of the 6th hit",    "fTracks.SoLIDMCTrack.GetHitPX5()"         },
      { "track.0.py",       "y coordinate of the 1st hit",    "fTracks.SoLIDMCTrack.GetHitPY0()"         },
      { "track.1.py",       "y coordinate of the 2nd hit",    "fTracks.SoLIDMCTrack.GetHitPY1()"         },
      { "track.2.py",       "y coordinate of the 3rd hit",    "fTracks.SoLIDMCTrack.GetHitPY2()"         },
      { "track.3.py",       "y coordinate of the 4th hit",    "fTracks.SoLIDMCTrack.GetHitPY3()"         },
      { "track.4.py",       "y coordinate of the 5th hit",    "fTracks.SoLIDMCTrack.GetHitPY4()"         },
      { "track.5.py",       "y coordinate of the 6th hit",    "fTracks.SoLIDMCTrack.GetHitPY5()"         },
      { "track.0.pz",       "z coordinate of the 1st hit",    "fTracks.SoLIDMCTrack.GetHitPZ0()"         },
      { "track.1.pz",       "z coordinate of the 2nd hit",    "fTracks.SoLIDMCTrack.GetHitPZ1()"         },
      { "track.2.pz",       "z coordinate of the 3rd hit",    "fTracks.SoLIDMCTrack.GetHitPZ2()"         },
      { "track.3.pz",       "z coordinate of the 4th hit",    "fTracks.SoLIDMCTrack.GetHitPZ3()"         },
      { "track.4.pz",       "z coordinate of the 5th hit",    "fTracks.SoLIDMCTrack.GetHitPZ4()"         },
      { "track.5.pz",       "z coordinate of the 6th hit",    "fTracks.SoLIDMCTrack.GetHitPZ5()"         },
      { "track.0.predx",       "x coordinate of the 1st hit",    "fTracks.SoLIDTrack.GetPredHitX0()"         },
      { "track.1.predx",       "x coordinate of the 2nd hit",    "fTracks.SoLIDTrack.GetPredHitX1()"         },
      { "track.2.predx",       "x coordinate of the 3rd hit",    "fTracks.SoLIDTrack.GetPredHitX2()"         },
      { "track.3.predx",       "x coordinate of the 4th hit",    "fTracks.SoLIDTrack.GetPredHitX3()"         },
      { "track.4.predx",       "x coordinate of the 5th hit",    "fTracks.SoLIDTrack.GetPredHitX4()"         },
      { "track.5.predx",       "x coordinate of the 6th hit",    "fTracks.SoLIDTrack.GetPredHitX5()"         },
      { "track.0.predy",       "y coordinate of the 1st hit",    "fTracks.SoLIDTrack.GetPredHitY0()"         },
      { "track.1.predy",       "y coordinate of the 2nd hit",    "fTracks.SoLIDTrack.GetPredHitY1()"         },
      { "track.2.predy",       "y coordinate of the 3rd hit",    "fTracks.SoLIDTrack.GetPredHitY2()"         },
      { "track.3.predy",       "y coordinate of the 4th hit",    "fTracks.SoLIDTrack.GetPredHitY3()"         },
      { "track.4.predy",       "y coordinate of the 5th hit",    "fTracks.SoLIDTrack.GetPredHitY4()"         },
      { "track.5.predy",       "y coordinate of the 6th hit",    "fTracks.SoLIDTrack.GetPredHitY5()"         },
      { "track.0.predex",       "x coordinate of the 1st hit",    "fTracks.SoLIDTrack.GetPredHiteX0()"         },
      { "track.1.predex",       "x coordinate of the 2nd hit",    "fTracks.SoLIDTrack.GetPredHiteX1()"         },
      { "track.2.predex",       "x coordinate of the 3rd hit",    "fTracks.SoLIDTrack.GetPredHiteX2()"         },
      { "track.3.predex",       "x coordinate of the 4th hit",    "fTracks.SoLIDTrack.GetPredHiteX3()"         },
      { "track.4.predex",       "x coordinate of the 5th hit",    "fTracks.SoLIDTrack.GetPredHiteX4()"         },
      { "track.5.predex",       "x coordinate of the 6th hit",    "fTracks.SoLIDTrack.GetPredHiteX5()"         },
      { "track.0.predey",       "y coordinate of the 1st hit",    "fTracks.SoLIDTrack.GetPredHiteY0()"         },
      { "track.1.predey",       "y coordinate of the 2nd hit",    "fTracks.SoLIDTrack.GetPredHiteY1()"         },
      { "track.2.predey",       "y coordinate of the 3rd hit",    "fTracks.SoLIDTrack.GetPredHiteY2()"         },
      { "track.3.predey",       "y coordinate of the 4th hit",    "fTracks.SoLIDTrack.GetPredHiteY3()"         },
      { "track.4.predey",       "y coordinate of the 5th hit",    "fTracks.SoLIDTrack.GetPredHiteY4()"         },
      { "track.5.predey",       "y coordinate of the 6th hit",    "fTracks.SoLIDTrack.GetPredHiteY5()"         },
      { "track.0.tracker", "tracker ID of the 1st hit",      "fTracks.SoLIDMCTrack.GetHitTracker0()"   },
      { "track.1.tracker", "tracker ID of the 2nd hit",      "fTracks.SoLIDMCTrack.GetHitTracker1()"   },
      { "track.2.tracker", "tracker ID of the 3rd hit",      "fTracks.SoLIDMCTrack.GetHitTracker2()"   },
      { "track.3.tracker", "tracker ID of the 4th hit",      "fTracks.SoLIDMCTrack.GetHitTracker3()"   },
      { "track.4.tracker", "tracker ID of the 5th hit",      "fTracks.SoLIDMCTrack.GetHitTracker4()"   },
      { "track.5.tracker", "tracker ID of the 6th hit",      "fTracks.SoLIDMCTrack.GetHitTracker5()"   },
      { "track.0.signal",  "1 if the hit is a signal",       "fTracks.SoLIDMCTrack.IsSignalHit0()"   },
      { "track.1.signal",  "1 if the hit is a signal",       "fTracks.SoLIDMCTrack.IsSignalHit1()"   },
      { "track.2.signal",  "1 if the hit is a signal",       "fTracks.SoLIDMCTrack.IsSignalHit2()"   },
      { "track.3.signal",  "1 if the hit is a signal",       "fTracks.SoLIDMCTrack.IsSignalHit3()"   },
      { "track.4.signal",  "1 if the hit is a signal",       "fTracks.SoLIDMCTrack.IsSignalHit4()"   },
      { "track.5.signal",  "1 if the hit is a signal",       "fTracks.SoLIDMCTrack.IsSignalHit5()"   },
      { "track.nMC",       "signal hit number of the track", "fTracks.SoLIDMCTrack.GetNMCHits()"       },
      { "track.nhits",     "number of hits of the track",    "fTracks.SoLIDMCTrack.GetNHits()"       },
      { "track.ntrack",    "number of tracks for the event", "GetNTracks()"                          },
      { "track.goodsignal","1 if all GEMs have one hit from the signal", "fGoodSignalFlag"           },
      { "track.nseeds",    "number of seeds for tracking",   "GetNSeeds()"},
      { "track.seedeffi1",  "efficiency of first seeding",          "GetFirstSeedEfficiency()"},
      { "track.trackeffi1",  "efficiency of first seeding",          "GetFirstMCTrackEfficiency()"},
      { "track.seedeffi2",  "efficiency of second seeding",          "GetSecondSeedEfficiency()"},
      { "track.trackeffi2",  "efficiency of second seeding",          "GetSecondMCTrackEfficiency()"},
      { "track.chi2perndf", "chi2 per ndf",                  "fTracks.SoLIDMCTrack.GetChi2()"},
      { "track.p",              "momentum of the track",          "fTracks.SoLIDTrack.GetMomentum()"},
      { "track.theta",          "polar angle of the track",       "fTracks.SoLIDTrack.GetTheta()"},
      { "track.phi",            "azimuthal angle of the track",   "fTracks.SoLIDTrack.GetPhi()"},
      { "track.vertexz",        "vertex z of the track",          "fTracks.SoLIDTrack.GetVertexZ()"},
      { "track.deltaecx",        "delta x ec",          "fTracks.SoLIDTrack.GetMomMax()"},
      { "track.deltaecy",        "delta y ec",          "fTracks.SoLIDTrack.GetMomMin()"},
      { "track.deltaece",        "delta E ec in %",          "fTracks.SoLIDTrack.GetThetaMin()"}, 
      { 0 }
    };
    ret = DefineVarsFromList( mcvars, mode );
#endif
  }
  return ret;
}
//_____________________________________________________________________________
void SoLIDTrackerSystem::Print(const Option_t* opt) const
{
  THaTrackingDetector::Print(opt);

  Int_t verbose = 0;
  if( opt ) {
    TString opt_s(opt);
    opt_s.ToLower();
    verbose = opt_s.CountChar('v');
  }
}
//_____________________________________________________________________________
void SoLIDTrackerSystem::PrintDataBase(Int_t level) const
{
  if (level>=0 && level <4){
    if (level == 0){
      string out_prefix = "solid.trackersystem.";
      cout<<"******parameter from database for tracker system "<<fSystemID<<"******"<<endl;
      cout<<out_prefix<<fSystemID<<".MCdata = "<<TestBit(kMCData)<<endl;
      cout<<out_prefix<<fSystemID<<".ntracker = "<<fNTracker<<endl;
      cout<<out_prefix<<fSystemID<<".do_rawdecode = "<<TestBit(kDoRawDecode)<<endl;
      cout<<out_prefix<<fSystemID<<".do_coarsetrack = "<<TestBit(kDoCoarse)<<endl;
      cout<<out_prefix<<fSystemID<<".do_finetrack = "<<TestBit(kDoFine)<<endl;
      cout<<out_prefix<<fSystemID<<".do_chi2 = "<<TestBit(kDoChi2)<<endl;
      cout<<out_prefix<<fSystemID<<".chi2_cut = "<<fChi2Cut<<endl;
      cout<<out_prefix<<fSystemID<<".max_miss_hit = "<<fNMaxMissHit<<endl;
      cout<<"**********************************************************"<<endl;
    }else if (level > 0){
      level--;
      for (Int_t i=0; i<fNTracker; i++){
        fGEMTracker[i]->PrintDataBase(level);
      }
    }
  }
}
//_____________________________________________________________________________
void SoLIDTrackerSystem::SetDebug( Int_t level )
{
  
}
//_____________________________________________________________________________
Int_t SoLIDTrackerSystem::Begin( THaRunBase* r)
{
  return 0;
}
//_____________________________________________________________________________
Int_t SoLIDTrackerSystem::End( THaRunBase* r )
{
  return 0;
}
//_____________________________________________________________________________
const char* SoLIDTrackerSystem::GetDBFileName() const
{
  // Return database file name prefix. For SoLID trackers, this is the detector
  // name with any trailing sector number removed. In this way, all trackers
  // use the same database file.

  // fDBPrefix is set in MakePrefix()
  return fDBPrefix.Data();
}

//_____________________________________________________________________________
void SoLIDTrackerSystem::MakePrefix()
{
  // Set up name prefixes for global variables and database.
  // Global variables and database keys get the standard prefix,
  // e.g. "solid.tracker.3."
  // The database file name gets a different, shorter prefix, allowing
  // the trackers to share a single database file. For the example above,
  // "solid.tracker."

  THaDetector::MakePrefix();
  TString prefix( GetPrefix() );
  //prefix.Append("{allsectors}.");
  assert( prefix.EndsWith(".") );
  Int_t ndot = prefix.CountChar('.');
  if( ndot > 1 ) {
    prefix.Chop();
    if( !prefix.EndsWith(".") ){
      prefix.Remove( prefix.Last('.')+1 );
    }
    else{
      Warning( Here("MakePrefix"), "Double dot in detector prefix = "
	       "\"%s\"?", GetPrefix() );
	  }
  }
  fDBPrefix = prefix;
}


//_____________________________________________________________________________
Int_t SoLIDTrackerSystem::ReadGeometry( FILE* file, const TDatime& date,
				Bool_t /* required */ )
{
  DBRequest request[] = {
    { "phi",  &fPhi,  kDouble, 0, 0 },
    { 0 }
  };
  Int_t err = LoadDB( file, date, request, fPrefix );
  if( err )
    return err;

  // Keep phi in rad and normalize to -pi..pi
  fPhi = TVector2::Phi_mpi_pi( fPhi*TMath::DegToRad() );

  return kOK;
}

#ifdef MCDATA
//_____________________________________________________________________________
Int_t SoLIDTrackerSystem::FitMCPoints( Podd::MCTrack* mctrk ) const
{
  return 1;
}
#endif // MCDATA
//_____________________________________________________________________________
UInt_t SoLIDTrackerSystem::GetCrateMapDBcols() const
{
  // Return number of columns for the detector crate map in the database.
  // 5 columns means that the no resolution value is needed to analyze the
  // data, as is typical for ADCs

  return 5;
}
//_____________________________________________________________________________
inline
static DAQmodule* FindDAQmodule( UShort_t crate, UShort_t slot,
                                 const THashTable* table )
{
  assert(table);
  CSpair m( crate, slot );
  return static_cast<DAQmodule*>( table->FindObject(&m) );
}

//_____________________________________________________________________________
UInt_t SoLIDTrackerSystem::LoadDAQmodel( THaDetMap::Module* mod ) const
{
  // Update detector map module 'mod' with the model number from the cratemap
  DAQmodule* found = FindDAQmodule( mod->crate, mod->slot, fCrateMap );
  UInt_t num = found ? found->fModel : 0;
  mod->SetModel( num );
  return num;
}

//_____________________________________________________________________________
Double_t SoLIDTrackerSystem::LoadDAQresolution( THaDetMap::Module* mod ) const
{
  // Update detector map module 'mod' with the resolution from the cratemap
  DAQmodule* found = FindDAQmodule( mod->crate, mod->slot, fCrateMap );
  Double_t res = found ? found->fResolution : 0.0;
  mod->SetResolution( res );
  return res;
}

//_____________________________________________________________________________
UInt_t SoLIDTrackerSystem::GetDAQnchan( THaDetMap::Module* mod ) const
{
  // Return number of channels for detector map module 'mod' from cratemap
  DAQmodule* found = FindDAQmodule( mod->crate, mod->slot, fCrateMap );
  return found ? found->fNchan : 0;
}

///////////////////////////////////////////////////////////////////////////////

ClassImp(SoLIDTrackerSystem)

