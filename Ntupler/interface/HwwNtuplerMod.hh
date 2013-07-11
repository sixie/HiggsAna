#ifndef HIGGSANA_NTUPLER_HWWNTUPLERMOD_H
#define HIGGSANA_NTUPLER_HWWNTUPLERMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/TriggerMask.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"

#include "HiggsAnaDefs.hh"
#include "TEventInfo.hh"
#include "TGenInfo.hh"
#include <TClonesArray.h>
#include "TElectron.hh"
#include "TJet.hh"
#include "TMuon.hh"
#include "TPhoton.hh"
#include "TGenParticle.hh"
#include "TPFCandidate.hh"

#include <vector>

class TTree;
class TFile;
class TString;

namespace mithep
{
  class HwwNtuplerMod : public BaseMod
  {    
    public:
      HwwNtuplerMod(const char *name="HwwNtuplerMod", const char *title="BAMBU to ntuple");
      ~HwwNtuplerMod();	

      void SetOutputName(const char *f)  { fOutputName = f; }
      
      void SetUseGen(Bool_t flag)        { fUseGen = flag; } 
      void SetSkipIfHLTFail(Bool_t flag) { fSkipIfHLTFail = flag; } 
      void SetFakeRateSkim(Bool_t flag)  { fFakeRateSkim = flag; }
          
      void SetElePtMin(Double_t pt)      { fElePtMin = pt; }
      void SetElePtMax(Double_t pt)      { fElePtMax = pt; }
      void SetEleEtaMin(Double_t eta)    { fEleEtaMin = eta; }
      void SetEleEtaMax(Double_t eta)    { fEleEtaMax = eta; }
      void SetMuonPtMin(Double_t pt)     { fMuonPtMin = pt; }
      void SetMuonPtMax(Double_t pt)     { fMuonPtMax = pt; }
      void SetMuonEtaMin(Double_t eta)   { fMuonEtaMin = eta; }
      void SetMuonEtaMax(Double_t eta)   { fMuonEtaMax = eta; }
      void SetMassMin(Double_t m)        { fMassMin = m; }
      void SetMassMax(Double_t m)        { fMassMax = m; }
      void SetConversionName(const char *f)  { fConversionName = f; }
      void SetCleanJetsName(const char *f)  { fCleanJetsName = f; }
      void SetCleanJetsNoPtCutName(const char *name) { fCleanJetsNoPtCutName  = name; }
      void SetJetPtMin(Double_t et)      { fJetPtMin = et; }
      void SetPhotonEtMin(Double_t et)   { fPhotonEtMin = et; }
      void SetPrintHLT(Bool_t flag)      { fPrintTable = flag; }
      void SetFSRMode(Int_t mode)        { fFSRMode = mode; }
      void SetFillGenOnly(Bool_t b)      { fFillGenOnly = b; }
      void SetPDFName(const char *s)     { fPDFName            = s;    }
      void SetComputePDFWeights(Bool_t b) { fComputePDFWeights = b;    }
      void SetReadPileupInfo(Bool_t b)   { fReadPileupInfo = b;    }

      void AddTrigger(const char* name,  ULong_t id, 
                      ULong_t firstObjectId = 0, const char* firstObjectModuleName = "", 
                      ULong_t secondObjectId = 0, const char* secondObjectModuleName = "" ) {
	fTriggerNamesv.push_back(name);
	fTriggerIdsv.push_back(id);
	fFirstTriggerObjectModuleNamesv.push_back(firstObjectModuleName);
        fFirstTriggerObjectIdsv.push_back(firstObjectId);
	fSecondTriggerObjectModuleNamesv.push_back(secondObjectModuleName);
        fSecondTriggerObjectIdsv.push_back(secondObjectId);
 //         cout << "Add Trigger: " << name << " " << id << " : " << firstObjectModuleName << " " << firstObjectId << " : " 
//               << secondObjectModuleName << " " << secondObjectId << endl;
      }
      void AddL1Trigger(const char* name, ULong_t id) {
	fL1TriggerNamesv.push_back(name);
	fL1TriggerIdsv.push_back(id);
      }            
      void AddL1SeedModule(const char* name, ULong_t id) {
	fL1SeedModuleNamesv.push_back(name);
	fL1SeedModuleIdsv.push_back(id);
      }            
            
    protected:
      void Begin();
      void BeginRun();
      void EndRun();
      void SlaveBegin();
      void SlaveTerminate();
      void Terminate();
      void Process();

      // Fill generator info data object
      void FillGenInfo(const MCParticleCol *GenLeptons, const MCParticleCol *GenNeutrinos, 
                       const MCParticleCol *GenBosons, const MCParticleCol *GenPhotons);
      
      // Fill electron data object
      void FillElectron(const Electron *ele, Int_t isMCMatched, const Vertex *primaryVertices);
      
      // Fill muon data object
      void FillMuon(const Muon *mu, Int_t isMCMatched, const Vertex *primaryVertices );
      void FillMuon(const Track *mu, const Vertex* primaryVertex);

      // Fill jet data object
      void FillJet(const CaloJet *jet);
      void FillJet(const TrackJet *jet);
      void FillJet(const PFJet *jet, const Vertex *primaryVertices);
      
      // Fill photon data object
      void FillPhoton(const Photon *pho);
      
      // Fill PF candidates
      void FillPFCandidate(const PFCandidate *pf, 
                           UInt_t matchedObjectType, UInt_t matchedObjectIndex,
			   const Vertex* primaryVertex );


     // Match electron to HLT primitive
      ULong_t MatchL1(const Double_t eta, const Double_t phi);

      // Match electron to HLT primitive
      ULong_t MatchHLT(const Double_t eta, const Double_t phi, Double_t runNum = 0, Double_t evtNum = 0);

      // Match muon to HLT primitive
      ULong_t MatchHLTMuon(const Double_t pt, const Double_t eta, const Double_t phi, Double_t runNum = 0, Double_t evtNum = 0);
      
      TFile                  *fOutputFile;      // output file handle
      TString                 fOutputName;      // output file name
      
      TString                 fPDFName;         //PDF name
      TString                 fPartName;        // MC particle collection name
      TString                 fMCEvtInfoName;   // MC event info name
      TString                 fElectronName;    // electron collection name
      TString                 fMuonName;        // muon collection name
      TString                 fTrackName;       // track collection name
      TString                 fPrimVtxName;     // primary vertex collection name
      TString                 fBeamSpotName;    // pointer to beam spot branch
      TString                 fCaloJetName;     //name of jet collection used in b-tagging
      TString                 fCleanJetsName;   // clean jets collection
      TString                 fCleanJetsNoPtCutName;    //name of clean all jets collection
      TString                 fTrigMaskName;    // trigger mask name
      TString                 fCaloMetName;     // calorimeter MET collection name
      TString                 fTCMetName;       // track-corrected MET collection name
      TString                 fPFMetName;       // particle flow MET collection name
      TString                 fConversionName;  // conversion collection name         


      const MCParticleCol    *fParticles;       // MC particle collection handle
      const MCEventInfo      *fMCEvtInfo;       // MC event info handle
      const GenJetCol        *fGenJets;         // MC particle collection handle
      const ElectronCol      *fElectrons;       // electron collection handle
      const MuonCol          *fMuons;           // muon collection handle
      const TrackCol         *fTracks;          // track collection handle
      const VertexCol        *fPrimVerts;       // primary vertex collection handle
      const BeamSpotCol      *fBeamSpot;        // pointer to beam spot branch
      const CaloJetCol       *fCaloJets;        // calorimeter jet collection handle
      const TrackJetCol      *fTrackJets;       // track jet collection handle
      const PFJetCol         *fPFJets;          // particle flow jet collection handle
      const PFCandidateCol   *fPFCandidates;    // PF Candidates
      const PhotonCol        *fPhotons;         // photon collection handle
      const TriggerMask      *fTrigMask;        // trigger mask handle
      const L1TriggerMask    *fL1TrigMask;      // L1 trigger mask handle
      const CaloMetCol       *fCaloMet;         // calorimeter MET handle
      const MetCol           *fTCMet;           // track-corrected MET handle
      const PFMetCol         *fPFMet;           // particle flow MET handle
      const DecayParticleCol *fConversions;     // conversion collection handle 
      const PileupInfoCol    *fPileupInfo;
      const PileupEnergyDensityCol *fPileupEnergyDensity;
      const VertexCol        *fGoodPrimaryVertices;  // good primary vertex collection handle

      Bool_t                  fPrintDebug;
      Bool_t                  fUseGen;          // flag whether to look at generator info
      Bool_t                  fPrintTable;      // flag whether to print out HLT table
      Bool_t                  fSkipIfHLTFail;   // flag whether to skip event processing if HLT does not accept
      Bool_t                  fFakeRateSkim;    // flag to keep all electrons in ntuple or just those from selected dielectrons
      Bool_t                  fFillGenOnly;     // flag to skip reco objects       
      Bool_t                  fComputePDFWeights; //compute weights for different pdfs?
      Bool_t                  fReadPileupInfo;  

      Double_t                fElePtMin;        // minimum supercluster ET
      Double_t                fElePtMax;        // maximum supercluster ET
      Double_t                fEleEtaMin;       // minimum supercluster eta
      Double_t                fEleEtaMax;       // maximum supercluster eta
      Double_t                fMuonPtMin;       // minimum reco muon pT
      Double_t                fMuonPtMax;       // maximum reco muon pT
      Double_t                fMuonEtaMin;      // minimum reco muon eta
      Double_t                fMuonEtaMax;      // maximum reco muon eta
      Double_t                fMassMin;         // minimum reco dielectron mass
      Double_t                fMassMax;         // maximum reco dielectron mass
      Double_t                fJetPtMin;        // minimum jet ET
      Double_t                fPhotonEtMin;     // minimum photon ET
      
      TTree*                  fEventTree;       // event tree
                
      MuonTools               *fMuonTools;      // interface to tools for muon ID
            
      TEventInfo              fEventInfo;       // general event information
      TGenInfo                fGenInfo;         // generator information
      TClonesArray           *fElectronArr;     // electron array
      TClonesArray           *fMuonArr;         // muon array
      TClonesArray           *fPFJetArr;        // particle flow jet array
      TClonesArray           *fPhotonArr;       // photon array
      TClonesArray           *fGenParticleArr;  // genparticle array
      TClonesArray           *fPFCandidateArr;  // PFCandidate array
      TClonesArray           *fTestArr;  // PFCandidate array

      Int_t                   fNPDFMembers;     // Count How many PDF Members we have
      Float_t                 fPDFWeights[100]; // Save weights for each of the PDFs

      vector<TString>         fTriggerNamesv;       // names of triggers we're interested in 
      vector<ULong_t>         fTriggerIdsv;     // corresponding ETriggerBit value
      vector<TString>         fFirstTriggerObjectModuleNamesv; // 
      vector<TString>         fSecondTriggerObjectModuleNamesv; // 
      vector<ULong_t>         fFirstTriggerObjectIdsv;   // ETriggerObjectBit
      vector<ULong_t>         fSecondTriggerObjectIdsv;  // ETriggerObjectBit
      vector<TString>         fL1TriggerNamesv; // names of L1 triggers we're interested in 
      vector<ULong_t>         fL1TriggerIdsv;   // corresponding L1TriggerBit value
      vector<TString>         fL1SeedModuleNamesv; // names of L1 Seed Modules
      vector<ULong_t>         fL1SeedModuleIdsv;   // corresponding L1TriggerBit value

      Int_t                   fFSRMode;         // flag to indicate how to parse GEN-level tree (depends on MC tool)
                                                // 0 -> PYTHIA: traverse down tree from status=3 boson
					        // 1 -> HORACE: look at immediate daughters of status=3 boson
					        // 2 -> PHOTOS: look at immediate daughters of status=2 boson       
    ClassDef(HwwNtuplerMod,1)
  };
}

#endif    
