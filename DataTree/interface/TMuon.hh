#ifndef HIGGSANA_NTUPLER_TMUON_HH
#define HIGGSANA_NTUPLER_TMUON_HH

#include <TObject.h>

namespace higgsana 
{
  class TMuon : public TObject
  {
    public:
      TMuon(){}
      ~TMuon(){} 
  
      Float_t pt, eta, phi;           // kinematics
      Float_t pfPt, pfEta, pfPhi;     // kinematics
      Float_t pterr;                  // pt error
      Float_t staPt, staEta, staPhi;  // standalone muon measurements
      Float_t trkIso03;	              // track isolation
      Float_t emIso03;	              // ECAL-based isolation
      Float_t hadIso03;	              // HCAL-based isolation
      Float_t hoIso03;	              // HO-based isolation
      Float_t trkIso05;	              // track isolation
      Float_t emIso05;	              // ECAL-based isolation
      Float_t hadIso05;	              // HCAL-based isolation
      Float_t hoIso05;	              // HO-based isolation
      
      Float_t d0, d0Err, dz;          // impact parameter
      Float_t tkNchi2;	              // track chi^2/ndf 
      Float_t muNchi2;	              // global muon chi^2/ndf
      Int_t q;		              // charge
      Int_t nValidHits;	              // number of valid hits in muon system
      UInt_t qualityBits;             // bits for various muon quality criteria
      UInt_t typeBits;	              // global muon, tracker muon, or standalone muon
      UInt_t nTkHits;	              // number of inner tracker hits
      UInt_t nPixHits;	              // number of pixel hits
      UInt_t nSeg;  	              // number of muon segments
      UInt_t nMatch;                  // number of muon chambers matched to segments
      ULong_t hltMatchBits;           // bits for matching with HLT primitives 
      Int_t  isMCReal;	              // used in a dimuon combination?     
      Float_t ip3d, ip3dSig;          // IP3D PV

      Float_t SegmentCompatibility;
      Float_t CaloCompatilibity;

      Float_t TrkKink;
      Float_t GlobalKink;
      Float_t HadEnergy;
      Float_t HadS9Energy;
      Float_t HoEnergy;
      Float_t HoS9Energy;
      Float_t EmEnergy;
      Float_t EmS9Energy;
      Bool_t  PassPFId;

      Float_t PFMuonEEcal;
      Float_t PFMuonEtEcal;
      Float_t trkLayers;

//       Int_t   Station0NSegments;
//       Float_t Station0TrackDist;
//       Float_t Station0TrackDistErr;
//       Float_t Station0dX;
//       Float_t Station0dY;
//       Float_t Station0PullX;
//       Float_t Station0PullY;

//       Int_t   Station1NSegments;
//       Float_t Station1TrackDist;
//       Float_t Station1TrackDistErr;
//       Float_t Station1dX;
//       Float_t Station1dY;
//       Float_t Station1PullX;
//       Float_t Station1PullY;

//       Int_t   Station2NSegments;
//       Float_t Station2TrackDist;
//       Float_t Station2TrackDistErr;
//       Float_t Station2dX;
//       Float_t Station2dY;
//       Float_t Station2PullX;
//       Float_t Station2PullY;

//       Int_t   Station3NSegments;
//       Float_t Station3TrackDist;
//       Float_t Station3TrackDistErr;
//       Float_t Station3dX;
//       Float_t Station3dY;
//       Float_t Station3PullX;
//       Float_t Station3PullY;

//       Int_t   Station4NSegments;
//       Float_t Station4TrackDist;
//       Float_t Station4TrackDistErr;
//       Float_t Station4dX;
//       Float_t Station4dY;
//       Float_t Station4PullX;
//       Float_t Station4PullY;

//       Int_t   Station5NSegments;
//       Float_t Station5TrackDist;
//       Float_t Station5TrackDistErr;
//       Float_t Station5dX;
//       Float_t Station5dY;
//       Float_t Station5PullX;
//       Float_t Station5PullY;

//       Int_t   Station6NSegments;
//       Float_t Station6TrackDist;
//       Float_t Station6TrackDistErr;
//       Float_t Station6dX;
//       Float_t Station6dY;
//       Float_t Station6PullX;
//       Float_t Station6PullY;

//       Int_t   Station7NSegments;
//       Float_t Station7TrackDist;
//       Float_t Station7TrackDistErr;
//       Float_t Station7dX;
//       Float_t Station7dY;
//       Float_t Station7PullX;
//       Float_t Station7PullY;


      ClassDef(TMuon,3)
        };  
}
#endif
