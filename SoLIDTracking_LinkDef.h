#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
//#pragma link C++ namespace SoLID;

#pragma link C++ class SoLIDSpectrometer+;
#pragma link C++ class SolTrackInfo+;
#pragma link C++ class SoLIDTrackerSystem+;
#pragma link C++ class SoLIDGEMTracker+;
#pragma link C++ class SoLIDGEMChamber+;
#pragma link C++ class SoLIDGEMReadOut+;
#pragma link C++ class SoLIDECal+;
#pragma link C++ class Hit+;
#pragma link C++ class SoLIDRawHit+;
#pragma link C++ class SoLIDGEMHit+;
#pragma link C++ class SoLIDTrack+;
#pragma link C++ class SoLIDFieldMap+;
#pragma link C++ class SoLKalTrackFinder+;
#pragma link C++ class SoLKalMatrix+;
#pragma link C++ class SoLKalTrackSystem+;
#pragma link C++ class SoLKalTrackState+;
#pragma link C++ class SoLKalTrackSite+;
#pragma link C++ class SoLKalFieldStepper+;
#ifdef MCDATA
#pragma link C++ class SoLIDMCRawHit+;
#pragma link C++ class SoLIDMCGEMHit+;
#pragma link C++ class SoLIDMCTrack+;
#endif
#endif

