//*-- Author:  Ole Hansen<mailto:ole@jlab.org>; Jefferson Lab; 06-Feb-2012
//*-- Modified by: Weizhi Xiong<weizhi.xiong@duke.edu>; 28-Mar-2016
#ifndef ROOT_SoLID_Tracker_System
#define ROOT_SoLID_Tracker_System

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SoLIDTrackSystem                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//c++
#include <vector>
#include <utility>
#include <set>
#include <list>
#include <cassert>
#include <exception>
//root
#include "TMatrixDSym.h"
#include "TRotation.h"
#include "THashTable.h"
//Hall A analyzer
#include "THaTrackingDetector.h"
#include "THaDetMap.h"
#include "THaAnalysisObject.h"
#ifndef MCDATA
#include "THaEvData.h"
#else
#include "SimDecoder.h"
#endif
//SoLID tracking
#include "SoLIDGEMTracker.h"
#include "SoLIDECal.h"
#include "SoLIDFieldMap.h"
#include "SoLKalTrackFinder.h"

class SoLIDTrackerSystem : public THaTrackingDetector {
  public:
    SoLIDTrackerSystem( const char* name, const char* description = "", THaApparatus* app = 0 );
    SoLIDTrackerSystem() : fECal(0), fTracks(0), fTrackFinder(0) {}
    virtual ~SoLIDTrackerSystem();
    
    virtual Int_t   ReadDatabase( const TDatime& date );
    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    virtual EStatus Init( const TDatime& date );
    virtual Int_t   CoarseTrack( TClonesArray& tracks );
    virtual Int_t   FineTrack( TClonesArray& tracks );
    virtual Int_t   DefineVariables( EMode mode = kDefine );
    virtual void    Print(const Option_t* opt) const;
    virtual void    PrintDataBase(Int_t level) const;
    virtual void    SetDebug( Int_t level );
    virtual Int_t   Begin( THaRunBase* r=0 );
    virtual Int_t   End( THaRunBase* r=0 );
    Double_t GetPhi() const { return fPhi; }//in rad -pi to pi
    
    // Helper functions for getting DAQ module parameters - used by Init
    UInt_t    LoadDAQmodel( THaDetMap::Module* m ) const;
    Double_t  LoadDAQresolution( THaDetMap::Module* m ) const;
    UInt_t    GetDAQnchan( THaDetMap::Module* m ) const;
    
    void    SetSystemID( Int_t i ) { fSystemID = i; }
    Int_t   GetSystemID() const    { return fSystemID; }
    Int_t   GetNTracks()  const    { return fTracks->GetLast() + 1; }
    Int_t   GetNSeeds()   const    { return fTrackFinder->GetNSeeds(); }
    bool    GetFirstSeedEfficiency() const { return fTrackFinder->GetSeedEfficiency(0);} 
    bool    GetFirstMCTrackEfficiency() const { return fTrackFinder->GetMCTrackEfficiency(0);}
    bool    GetSecondSeedEfficiency() const { return fTrackFinder->GetSeedEfficiency(1);}
    bool    GetSecondMCTrackEfficiency() const { return fTrackFinder->GetMCTrackEfficiency(1);}  
 
    SoLIDGEMTracker * GetTracker(Int_t i) const {
      if (i >= 0 && i<fNTracker ) { return fGEMTracker[i]; }
      else { return 0; }
    }
     
     // Analysis control flags. Set via database.
    enum {
#ifdef MCDATA
      kMCData        = BIT(16), // Assume input is Monte Carlo data
#endif
      kDoRawDecode   = BIT(17),
      kDoCoarse      = BIT(18), // Do coarse tracking (if unset, decode only)
      kDoFine        = BIT(19), // Do fine tracking (implies kDoCoarse)
      kDoChi2        = BIT(20), // Apply chi2 cut to 3D tracks
    };


  protected:

		virtual UInt_t GetCrateMapDBcols() const;

     // Podd interface
    virtual const char* GetDBFileName() const;
    virtual void  MakePrefix();
    virtual Int_t ReadGeometry( FILE* file, const TDatime& date,
				Bool_t required = kTRUE );

    // Geometry
    Int_t          fSystemID;       // track system ID, 0 for SIDIS and J/Psi, 0~29 for PVDIS
    Double_t       fPhi;            // Phi rotation of this tracker system
    TRotation      fRotation;       // Rotation Tracker -> lab frame
    TRotation      fInvRot;         // Rotation lab frame -> Tracker
    

    std::vector<SoLIDGEMTracker*> fGEMTracker;
    SoLIDECal*     fECal;
    SoLIDFieldMap* fFieldMap;
    
    TClonesArray*  fTracks;         // array for storing SoLIDTrack(SoLIDMCTrack) objects
    // Configuration
    TString        fDBPrefix;       // Safe storage for database file name prefix
    THashTable*    fCrateMap;
    Int_t          fNTracker;       //total number of GEM detectors in this system SIDIS:6, PVDIS:5
    Double_t       fChi2Cut;        //chi2 cut after fitting the track
    Int_t          fNMaxMissHit;    //maximum number of hits that is allowed in the coarse tracking
    
    
    SoLKalTrackFinder* fTrackFinder; 
#ifdef MCDATA
    const Podd::SimDecoder* fMCDecoder; //! MC data decoder (if kMCdata)
    Bool_t fChecked;
    virtual Int_t FitMCPoints( Podd::MCTrack* mctrk ) const;
    Int_t         fGoodSignalFlag;
#endif

    ClassDef(SoLIDTrackerSystem,0)   // Collection of GEM trackers in one SoLID sector
};



///////////////////////////////////////////////////////////////////////////////

#endif
