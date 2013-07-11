#include "HiggsAna/BAMBUNtupler/interface/HwwNtuplerMod.hh"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/Track.h"
#include "MitAna/DataTree/interface/Vertex.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/CaloJet.h"
#include "MitAna/DataTree/interface/TrackJet.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/L1TriggerMask.h"
#include "MitAna/DataTree/interface/TriggerTable.h"
#include "MitAna/DataTree/interface/TriggerObjectsTable.h"
#include "MitAna/DataTree/interface/TriggerName.h"
#include "MitAna/DataTree/interface/CaloMet.h"
#include "MitAna/DataTree/interface/Met.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/DecayParticle.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/GeneratorTools.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include "LHAPDF/LHAPDF.h"
#include <vector>
#include <TSystem.h>

//#define __DATA_36X__

using namespace mithep;
using namespace higgsana;

ClassImp(mithep::HwwNtuplerMod)

HwwNtuplerMod::HwwNtuplerMod(const char *name, const char *title):
  BaseMod        (name,title),
  fOutputFile    (0),
  fOutputName    ("ntuple.root"),
  fPDFName       ("cteq66.LHgrid"),
  fPartName      (Names::gkMCPartBrn),
  fMCEvtInfoName (Names::gkMCEvtInfoBrn),
  fElectronName  (Names::gkElectronBrn),
  fMuonName     (Names::gkMuonBrn),
  fTrackName     (Names::gkTrackBrn),
  fPrimVtxName   (Names::gkPVBeamSpotBrn),
  fBeamSpotName  ("BeamSpot"),
  fCaloJetName   ("AKt5Jets"),
  fCleanJetsName ("NoDefaultNameSet"),
  fCleanJetsNoPtCutName("NoDefaultNameSet"),
  fTrigMaskName  (Names::gkHltBitBrn),
  fCaloMetName   (Names::gkCaloMetBrn),
  fTCMetName     ("TCMet"),
  fPFMetName     ("PFMet"),
  fConversionName(Names::gkMvfConversionBrn),
  fParticles     (0),
  fMCEvtInfo     (0),
  fGenJets       (0),
  fElectrons     (0),
  fTracks        (0),
  fPrimVerts     (0),
  fBeamSpot      (0),
  fCaloJets      (0),
  fTrackJets     (0),
  fPFJets        (0),
  fPFCandidates  (0),
  fPhotons       (0),
  fTrigMask      (0),
  fCaloMet       (0),
  fTCMet         (0),
  fPFMet         (0),
  fConversions   (0),  
  fGoodPrimaryVertices (0),
  fUseGen        (kTRUE),
  fPrintTable    (kFALSE),
  fSkipIfHLTFail (kFALSE),
  fFakeRateSkim  (kFALSE),
  fFillGenOnly   (kFALSE),
  fComputePDFWeights(kFALSE),
  fReadPileupInfo(kTRUE),
  fElePtMin       (5),
  fElePtMax       (14000),
  fEleEtaMin      (-3),
  fEleEtaMax      (3),
  fMuonPtMin     (5),
  fMuonPtMax     (14000),
  fMuonEtaMin    (-3),
  fMuonEtaMax    (3),
  fMassMin       (0),
  fMassMax       (1000),
  fJetPtMin      (10),
  fPhotonEtMin   (10), 
  fEventTree     (0),
  fFSRMode       (0)
{
  // Constructor

  // Don't write TObject part of the objects
  higgsana::TEventInfo::Class()->IgnoreTObjectStreamer();
  higgsana::TGenInfo::Class()->IgnoreTObjectStreamer();
  higgsana::TElectron::Class()->IgnoreTObjectStreamer();
  higgsana::TMuon::Class()->IgnoreTObjectStreamer();
  higgsana::TJet::Class()->IgnoreTObjectStreamer();
  higgsana::TPhoton::Class()->IgnoreTObjectStreamer();
}

//--------------------------------------------------------------------------------------------------
HwwNtuplerMod::~HwwNtuplerMod()
{
  // Destructor
}
	//--------------------------------------------------------------------------------------------------      
void HwwNtuplerMod::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::BeginRun()
{
//   if(HasHLTInfo() && fPrintTable) { GetHLTTable()->Print(); GetL1AlgoTable()->Print(); }
}

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::SlaveBegin()
{
  //
  // Request BAMBU branches
  //
  ReqBranch(fPartName,      fParticles); 
  ReqBranch(Names::gkPVBrn, fPrimVerts); 
  ReqBranch(fMCEvtInfoName, fMCEvtInfo);
//   ReqBranch(Names::gkGenJetBrn , fGenJets); 
  ReqBranch(fElectronName,  fElectrons);
  ReqBranch(fMuonName,      fMuons);
  ReqBranch(fTrackName,     fTracks);
  ReqBranch(Names::gkPFCandidatesBrn,     fPFCandidates);
  ReqBranch(fBeamSpotName,  fBeamSpot);
  ReqBranch(fTrigMaskName,  fTrigMask);
  ReqBranch("L1AlgoBitsBeforeMask",  fL1TrigMask);
  ReqBranch(fTCMetName,     fTCMet);
  ReqBranch(fPFMetName,     fPFMet);
  ReqBranch(fConversionName,fConversions);
  ReqBranch(fCaloJetName,   fCaloJets);
  ReqBranch(Names::gkPhotonBrn,    fPhotons);
  ReqBranch(Names::gkPFCandidatesBrn, fPFCandidates);
  cout << "ReadPileup: " << fReadPileupInfo << endl;
  cout << "UseGen : " << fUseGen << endl;
  if (fReadPileupInfo) {
    cout << "useGen: " << fUseGen << endl;
    if (Bool_t(fUseGen)) ReqBranch(Names::gkPileupInfoBrn,   fPileupInfo);
    ReqBranch(Names::gkPileupEnergyDensityBrn, fPileupEnergyDensity);
  }

  //
  // Set up arrays
  //
  fElectronArr   = new TClonesArray("higgsana::TElectron");   assert(fElectronArr);
  fMuonArr       = new TClonesArray("higgsana::TMuon");       assert(fMuonArr);
  fPFJetArr      = new TClonesArray("higgsana::TJet");        assert(fPFJetArr);
  fCaloJetArr    = new TClonesArray("higgsana::TCaloJet");    assert(fCaloJetArr);
  fPhotonArr     = new TClonesArray("higgsana::TPhoton");     assert(fPhotonArr);
  fGenParticleArr = new TClonesArray("higgsana::TGenParticle");     assert(fGenParticleArr);
  fPFCandidateArr = new TClonesArray("higgsana::TPFCandidate",2000);     assert(fPFCandidateArr);


  //
  // Create output file
  //
  fOutputFile = new TFile(fOutputName, "RECREATE");

  //
  // Initialize data trees and structs
  // 
  fEventTree = new TTree("Events","Events");

  fEventTree->Branch("Info",&fEventInfo);
  if(fUseGen) {
    fEventTree->Branch("Gen",&fGenInfo);
    fEventTree->Branch("GenParticle",&fGenParticleArr);  
  }

  fEventTree->Branch("Electron",   &fElectronArr);
  fEventTree->Branch("Muon",       &fMuonArr);
  fEventTree->Branch("PFJet",      &fPFJetArr);
  fEventTree->Branch("CaloJet",    &fCaloJetArr);
  fEventTree->Branch("Photon",     &fPhotonArr);
  fEventTree->Branch("PFCandidate",&fPFCandidateArr);

  if (fComputePDFWeights) {
 //    gSystem->Setenv("LHAPATH","/afs/cern.ch/cms/sw/slc4_ia32_gcc345/external/lhapdf/5.6.0-cms4/share/lhapdf/PDFsets/");
    LHAPDF::setVerbosity(LHAPDF::SILENT);
    LHAPDF::initPDFSet(fPDFName.Data());
    LHAPDF::getDescription();
    
    //Branches for PDFset weights
    fEventTree->Branch("NPDFMembers",&fNPDFMembers,"NPDFMembers/I");  
    fEventTree->Branch("PDFWeights",fPDFWeights,"PDFWeights[NPDFMembers]/F");
  }

//   AddOutput(fEventTree);

  fMuonTools = new MuonTools();
  
  //photon regression initialization
  egcorHgg2011 = new EGEnergyCorrector;
  egcorHgg2012 = new EGEnergyCorrector;
  egcorHgg2011->Initialize("4_2",gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFixGRPV22.dat"),gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/gbrv2ph_52x.root"));
  egcorHgg2012->Initialize("4_2",gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/PhotonFixGRPV22.dat"),gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/gbrv3ph_52x.root"));

}

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::SlaveTerminate()
{
  //
  // Save to ROOT file
  //
  fEventTree->Print();

  fOutputFile->cd();
  fOutputFile->Write();
  fOutputFile->Close();
  
  delete fElectronArr;
  delete fMuonArr;
  delete fPFJetArr;
  delete fCaloJetArr;
  delete fPhotonArr;
  delete fGenParticleArr;
  delete fPFCandidateArr;
  delete egcorHgg2011;
  delete egcorHgg2012;

}  

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::Process()
{
  //
  // Load branches
  //
  if(fUseGen) LoadBranch(fPartName);
  if(fUseGen) LoadBranch(fMCEvtInfoName);
//   if(fUseGen) LoadBranch(Names::gkGenJetBrn);
  if (!fFillGenOnly) {
    LoadBranch(Names::gkPVBrn);
    LoadBranch(fElectronName);
    LoadBranch(fMuonName);
    LoadBranch(fTrackName);
    LoadBranch(Names::gkPFCandidatesBrn);
    LoadBranch(fBeamSpotName);
    LoadBranch(fTrigMaskName);
    LoadBranch("L1AlgoBitsBeforeMask");
    LoadBranch(fTCMetName);
    LoadBranch(fPFMetName); 
    LoadBranch(fConversionName);
    LoadBranch(fCaloJetName);
    LoadBranch(Names::gkPhotonBrn);
    LoadBranch(Names::gkPFCandidatesBrn);
    if (fReadPileupInfo) {
      if (Bool_t(fUseGen)) LoadBranch(Names::gkPileupInfoBrn);
      LoadBranch(Names::gkPileupEnergyDensityBrn);
    }
  }


 //***********************************************************************************************
  //Some MC Specific Requirements
  //***********************************************************************************************
  const MCParticleCol *GenLeptons = 0;
  GenLeptons = GetObjThisEvt<MCParticleCol>(ModNames::gkMCLeptonsName,0);
  const MCParticleCol *GenNeutrinos = GetObjThisEvt<MCParticleCol>(ModNames::gkMCNeutrinosName,0);
  const MCParticleCol *GenTaus = GetObjThisEvt<MCParticleCol>(ModNames::gkMCTausName,0);
  const MCParticleCol *GenPhotons = GetObjThisEvt<MCParticleCol>(ModNames::gkMCPhotonsName,0);  
  const MCParticleCol *GenBosons = GetObjThisEvt<MCParticleCol>(ModNames::gkMCBosonsName,0);

  //exclude Wgamma , Zgamma
  Bool_t FoundZ = false;
  Bool_t FoundZGamma = false;
  Bool_t FoundWGammaEvent = false;
  Bool_t FoundW = false;


  if (fUseGen) {

    //clear gen particle array
    fGenParticleArr->Clear();

    for (UInt_t i=0; i<GenBosons->GetEntries(); i++) {  
      if (GenBosons->At(i)->AbsPdgId() == 23) {
        FoundZ = true;
      }
    }
    
    for (UInt_t i=0; i<GenPhotons->GetEntries(); i++) {  
      if (GenPhotons->At(i)->Pt() > 10.0) {
        
        //ISR Photon
        if ( (GenPhotons->At(i)->Mother() && GenPhotons->At(i)->Mother()->IsParton())
            || (GenPhotons->At(i)->Mother()->AbsPdgId() == 22 
                && GenPhotons->At(i)->Mother()->Status() ==3  
                && GenPhotons->At(i)->Mother()->Mother() 
                && GenPhotons->At(i)->Mother()->Mother()->IsParton())
          ) {
          FoundZGamma = true;    
        }
        
        //Pythia FSR
        if (GenPhotons->At(i)->Mother() && (GenPhotons->At(i)->Mother()->Status() == 3 || GenPhotons->At(i)->Mother()->Status() == 2)
            && (GenPhotons->At(i)->Mother()->AbsPdgId() == 11 
                || GenPhotons->At(i)->Mother()->AbsPdgId() == 13
                || GenPhotons->At(i)->Mother()->AbsPdgId() == 15)
          ) {
          CompositeParticle *object = new CompositeParticle();
          object->AddDaughter(GenPhotons->At(i));
          object->AddDaughter(GenPhotons->At(i)->Mother());
          if(object->Mass() > 1.0) FoundZGamma = true;
          delete object;          
        }
      }
    }
    
    //For WJets sample: Remove FSR W+gamma Events.
   
    for (UInt_t i=0; i<GenBosons->GetEntries(); i++) {  
      if (GenBosons->At(i)->AbsPdgId() == 24) {
        FoundW = true;
      }
    }

    for (UInt_t i=0; i<GenPhotons->GetEntries(); i++) {  
      if (GenPhotons->At(i)->Pt() > 10.0) {
        
        //ISR Photon
        if ((GenPhotons->At(i)->Mother() && GenPhotons->At(i)->Mother()->IsParton())
            || (GenPhotons->At(i)->Mother()->AbsPdgId() == 22 
                && GenPhotons->At(i)->Mother()->Status() ==3  
                && GenPhotons->At(i)->Mother()->Mother() 
                && GenPhotons->At(i)->Mother()->Mother()->IsParton())
          ) {
          FoundWGammaEvent = true;    
        }
        
        //WWgamma vertex
        if ((GenPhotons->At(i)->Mother() && GenPhotons->At(i)->Mother()->AbsPdgId() == 24) 
            || 
            (GenPhotons->At(i)->Mother()->AbsPdgId() == 22 
             && GenPhotons->At(i)->Mother()->Status() == 3
             && GenPhotons->At(i)->Mother()->Mother()
             && GenPhotons->At(i)->Mother()->Mother()->AbsPdgId() == 24
              )
          ) {
          FoundWGammaEvent = true;
        }
        
        //Pythia FSR
        if (GenPhotons->At(i)->Mother() && (GenPhotons->At(i)->Mother()->Status() == 3 || GenPhotons->At(i)->Mother()->Status() == 2)
            && (GenPhotons->At(i)->Mother()->AbsPdgId() == 11 
                || GenPhotons->At(i)->Mother()->AbsPdgId() == 13
                || GenPhotons->At(i)->Mother()->AbsPdgId() == 15)
          ) {
          CompositeParticle *object = new CompositeParticle();
          object->AddDaughter(GenPhotons->At(i));
          object->AddDaughter(GenPhotons->At(i)->Mother());
          if(object->Mass() > 1.0) FoundWGammaEvent = true;
          delete object;          
        }
      }
    } //for all gen photons   
   
    //Fill Gen Info
    FillGenInfo(GenLeptons, GenNeutrinos, GenBosons, GenPhotons);  // fill the data structure



  }

  fPrintDebug = kFALSE;
  if ( (GetEventHeader()->RunNum() == 162926 && GetEventHeader()->LumiSec() == 647 && GetEventHeader()->EvtNum() == 389987568 )
       || 
       (GetEventHeader()->RunNum() == 162926 && GetEventHeader()->LumiSec() == 647 && GetEventHeader()->EvtNum() == 390151968 )
    ) {
    fPrintDebug = kTRUE;    
  }
  if (fPrintDebug) {
    cout << "Event " << GetEventHeader()->RunNum() << " " << GetEventHeader()->LumiSec() << " " << GetEventHeader()->EvtNum() << "\n";
  }

  //************************************************************************************************
  //Fill PDF information
  //************************************************************************************************
  if (fComputePDFWeights && fUseGen) {
    //cout << "Event " << fMCEvtInfo->Weight() << " " << fMCEvtInfo->Id1() << " " << fMCEvtInfo->Id2() << " " << fMCEvtInfo->X1() << " " << fMCEvtInfo->X2() << endl;    
    Double_t Q    = fMCEvtInfo->Scale();
    Int_t    id1  = fMCEvtInfo->Id1();
    Double_t x1   = fMCEvtInfo->X1();
    Double_t pdf1 = fMCEvtInfo->Pdf1();
    Int_t    id2  = fMCEvtInfo->Id2();
    Double_t x2   = fMCEvtInfo->X2();
    Double_t pdf2 = fMCEvtInfo->Pdf2();
    
    UInt_t nmembers = LHAPDF::numberPDF() + 1;
    fNPDFMembers = nmembers;
    //    cout << "PDFset has " << nmembers << " members\n";
      
    //don't use default pdf numbers. use values of pdf member 0 as the default
    fPDFWeights[0] = 1;
    LHAPDF::usePDFMember(0);
    pdf1 =  LHAPDF::xfx(x1, Q, id1)/x1;
    pdf2 = LHAPDF::xfx(x2, Q, id2)/x2;   
      
    for (UInt_t i=1; i<nmembers; ++i) {
      LHAPDF::usePDFMember(i);
      Double_t newpdf1 = LHAPDF::xfx(x1, Q, id1)/x1;
      Double_t newpdf2 = LHAPDF::xfx(x2, Q, id2)/x2;
      Double_t TheWeight = newpdf1/pdf1*newpdf2/pdf2;
        
//             cout << i << " --> " << newpdf1 << " " << newpdf2 << " | " 
//                  << pdf1 << " "   << pdf2 << " | "
//                  << x1   << " "   << x2   << " | "
//                  << id1  << " "   << id2  << " | "
//                  << Q    << " : " <<  TheWeight << endl;
        
      fPDFWeights[i] = TheWeight;
    } 
  } // end if comp

  //If fill gen only, then skip the rest of the event.
  if (fFillGenOnly){ 
    fEventInfo.eventweight  = fMCEvtInfo->Weight();
    fEventInfo.runNum       = GetEventHeader()->RunNum();
    fEventInfo.evtNum       = GetEventHeader()->EvtNum();
    fEventInfo.lumiSec      = GetEventHeader()->LumiSec();
    fEventTree->Fill();

    return;
  }


  //************************************************************************************************
  //Obtain all the good objects from the event cleaning module
  //************************************************************************************************
//   ObjArray<Muon> *CleanMuons = dynamic_cast<ObjArray<Muon>* >(FindObjThisEvt(ModNames::gkCleanMuonsName));
//   ObjArray<Electron> *CleanElectrons = dynamic_cast<ObjArray<Electron>* >(FindObjThisEvt(ModNames::gkCleanElectronsName));
  ObjArray<Photon> *CleanPhotons = dynamic_cast<ObjArray<Photon>* >(FindObjThisEvt(ModNames::gkCleanPhotonsName));

//   JetOArr *CleanJets    = GetObjThisEvt<JetOArr>(fCleanJetsName);
  ObjArray<Jet> *CleanJetsNoPtCut = dynamic_cast<ObjArray<Jet>* >
    (FindObjThisEvt(fCleanJetsNoPtCutName.Data()));



  //************************************************************************************************
  //Fake Rate Skim
  //Keep only events with one denominator electron (v2,v3,or v4) and a jet with pt > 15
  //************************************************************************************************
  if (fFakeRateSkim) {
    Bool_t passSkim = kFALSE;
    ElectronOArr *electronsIsolated    = GetObjThisEvt<ElectronOArr>("ElectronsIsolated");
    ElectronOArr *electronsVBTF80NonIsolated    = GetObjThisEvt<ElectronOArr>("ElectronsVBTF80NonIsolated");
    ElectronOArr *electronsVBTF90PartiallyIsolated    = GetObjThisEvt<ElectronOArr>("ElectronsVBTF90PartiallyIsolated");
    MuonOArr *muonsDenominator    = GetObjThisEvt<MuonOArr>("DenominatorMuons");

    Bool_t foundJet = kFALSE;
    if (CleanJetsNoPtCut) {
      for (UInt_t k=0; k<CleanJetsNoPtCut->GetEntries() ; ++k) {

        if (CleanJetsNoPtCut->At(k)->Pt() < 15.0) continue;

        Bool_t hasOverlap = kFALSE;
        if (electronsIsolated) {
          for (UInt_t l = 0; l < electronsIsolated->GetEntries() ; ++l) {
            if (MathUtils::DeltaR(*electronsIsolated->At(l),*CleanJetsNoPtCut->At(k)) < 0.3)
              hasOverlap = kTRUE;
          }
        }
        if (electronsVBTF80NonIsolated) {
          for (UInt_t l = 0; l < electronsVBTF80NonIsolated->GetEntries() ; ++l) {
            if (MathUtils::DeltaR(*electronsVBTF80NonIsolated->At(l),*CleanJetsNoPtCut->At(k)) < 0.3)
              hasOverlap = kTRUE;
          }
        }
        if (electronsVBTF90PartiallyIsolated) {
          for (UInt_t l = 0; l < electronsVBTF90PartiallyIsolated->GetEntries() ; ++l) {
            if (MathUtils::DeltaR(*electronsVBTF90PartiallyIsolated->At(l),*CleanJetsNoPtCut->At(k)) < 0.3)
              hasOverlap = kTRUE;
          }
        }
        if (!hasOverlap) {
          foundJet = kTRUE;
          break;
        }
      }
    }

    Int_t NRecoElectrons = 0;
    for(UInt_t k=0; k<fElectrons->GetEntries(); ++k) {
      if (fElectrons->At(k)->Pt() > 10.0) NRecoElectrons++;
    } 

    //if we didn't find an electron and a jet then we don't save the event
    if ( NRecoElectrons + electronsIsolated->GetEntries() + 
         electronsVBTF80NonIsolated->GetEntries() + 
         electronsVBTF90PartiallyIsolated->GetEntries() >= 1 && foundJet         
      ) {
      passSkim = kTRUE;
    }


    if (CleanJetsNoPtCut) {
      for (UInt_t k=0; k<CleanJetsNoPtCut->GetEntries() ; ++k) {
        
        if (CleanJetsNoPtCut->At(k)->Pt() < 15.0) continue;
        
        Bool_t hasOverlap = kFALSE;
        if (muonsDenominator) {
          for (UInt_t l = 0; l < muonsDenominator->GetEntries() ; ++l) {
            if (MathUtils::DeltaR(*muonsDenominator->At(l),*CleanJetsNoPtCut->At(k)) < 0.3)
              hasOverlap = kTRUE;
          }
        }        
        if (!hasOverlap) {
          foundJet = kTRUE;
          break;
        }
      }
    }
    
    //if we didn't find a muon and a jet then we don't save the event
    if (muonsDenominator->GetEntries() >= 1 && foundJet) {
      passSkim = kTRUE;
    }

    if (!passSkim) {
      return;
    }

  }




  //
  // Get HLT info. Trigger objects can be matched by name to the corresponding trigger that passed.
  //
  ULong_t trigbits=0;
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();

    assert(hltTable);
    for(UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
       
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      if(fTrigMask->At(trigname->Id())) { trigbits |= fTriggerIdsv[itrig]; }

      if (fPrintDebug)
      {
        cout << "Trig : " << trigname->GetName() << " " << fTrigMask->At(trigname->Id())  << " "
	     << fTriggerIdsv[itrig] << " ::: " << trigbits << " " 
	     << endl;
      }
  
    }  
  }
  if(fSkipIfHLTFail && (trigbits==0))
    return;
  
  ULong_t l1trigbits = 0;
  if (fL1TrigMask) {
    const TriggerTable *l1Table = GetL1AlgoTable();
    assert(l1Table);
    for(UInt_t itrig=0; itrig<fL1TriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = l1Table->Get(fL1TriggerNamesv[itrig].Data());
      if(!trigname) continue;


      if (fPrintDebug)
      {
        cout << "Event " << GetEventHeader()->RunNum() << " " << GetEventHeader()->LumiSec() << " " << GetEventHeader()->EvtNum() << "\n";
        cout << "L1Trig : " << trigname->GetName() << " " << trigname->Id() << " ";
        cout << fL1TrigMask->At(trigname->Id()) << endl;
      }
      
      if(fL1TrigMask->At(trigname->Id())) { l1trigbits |= fL1TriggerIdsv[itrig]; }
    }  
  }


  //*****************************
  //Clear Arrays
  //*****************************
  IncNEventsProcessed();  
  fElectronArr->Clear();
  fMuonArr->Clear();
  fPFJetArr->Clear();
  fCaloJetArr->Clear();
  fPhotonArr->Clear();
  fPFCandidateArr->Clear();


  //
  // Get beam spot. If no beam spot information is available, default the coordinates to 99999
  //
  Double_t bsx=99999, bsy=99999, bsz=99999;
  if(fBeamSpot) {
    if(fBeamSpot->GetEntries() > 1) 
      std::cout << "********** More than 1 beam spot! **********" << std::endl;
    const BeamSpot *bs = fBeamSpot->At(0);
    bsx = bs->X();
    bsy = bs->Y();
    bsz = bs->Z();
  }


  //
  // Get primary vertex (for corrected d0)
  // Take the first primary vertex listed (should be ordered by sum-pT as in CMSSW)
  // NOTE: if no PV is found from fitting tracks, the beamspot is used
  //
  fGoodPrimaryVertices = GetObjThisEvt<VertexOArr>(ModNames::gkGoodVertexesName);
  const Vertex* primaryVertex = 0;
  Bool_t hasGoodPV = kFALSE;   
  if (fGoodPrimaryVertices->At(0)) {
    primaryVertex = fGoodPrimaryVertices->At(0);
    hasGoodPV = kTRUE;
  } 
  
  //-----------------------------------------
  //
  // Loop through electrons.
  //
  //-----------------------------------------
  vector<const Electron*> elev;  // array of pointers to preselected electrons ... 
  ElectronTools eleTools;        // helper class for electron ID decisions

  assert(fElectrons);                
  for(UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *ele = fElectrons->At(i);  
    if (ele->Pt()  < fElePtMin  || ele->Pt()  > fElePtMax ) continue;
    if (ele->SCluster()->Eta() < fEleEtaMin || ele->SCluster()->Eta() > fEleEtaMax) continue;
    
    elev.push_back(ele);       
  }
  
  for(UInt_t i=0; i<elev.size(); i++) {  

    //********************************************************************************
    //Match MC Truth
    //********************************************************************************
    Int_t isMCMatched = 0;
    Bool_t GenEleMatch = kFALSE;
    Bool_t GenMuMatch = kFALSE;
    Bool_t GenTauMatch = kFALSE;
    Bool_t GenPhotonMatch = kFALSE;
    if (GenLeptons) {
      for (UInt_t l=0; l<GenLeptons->GetEntries(); ++l) {
        if (MathUtils::DeltaR(*GenLeptons->At(l), *(elev[i])) < 0.3  ) {
          if (GenLeptons->At(l)->AbsPdgId() == 11) {
            GenEleMatch = kTRUE;
          } 
          if (GenLeptons->At(l)->AbsPdgId() == 13) {
            GenMuMatch = kTRUE;
          }           
        }
      }
    }
    if (GenTaus) {
      for (UInt_t l=0; l<GenTaus->GetEntries(); ++l) {
        if (MathUtils::DeltaR(*GenTaus->At(l), *(elev[i])) < 0.5 ) {
          GenTauMatch = kTRUE;
          break;
        }
      }
    }
    if (GenPhotons) {
      for (UInt_t l=0; l < GenPhotons->GetEntries(); l++) {  
        if (GenPhotons->At(l)->Pt() > 10.0 &&  MathUtils::DeltaR(*GenPhotons->At(l), *(elev[i])) < 0.3 ) {
        
          //ISR Photon
          if ( (GenPhotons->At(l)->Mother() && GenPhotons->At(l)->Mother()->IsParton())
               || (GenPhotons->At(l)->Mother()->AbsPdgId() == 22 
                   && GenPhotons->At(l)->Mother()->Status() ==3  
                   && GenPhotons->At(l)->Mother()->Mother() 
                   && GenPhotons->At(l)->Mother()->Mother()->IsParton())
            ) {
            GenPhotonMatch = kTRUE;    
          }
        
          //WWgamma vertex
          if ((GenPhotons->At(l)->Mother() && GenPhotons->At(l)->Mother()->AbsPdgId() == 24) 
              || 
              (GenPhotons->At(l)->Mother()->AbsPdgId() == 22 
               && GenPhotons->At(l)->Mother()->Status() == 3
               && GenPhotons->At(l)->Mother()->Mother()
               && GenPhotons->At(l)->Mother()->Mother()->AbsPdgId() == 24
                )
            ) {
            GenPhotonMatch = true;
          }

          //Pythia FSR
          if (GenPhotons->At(l)->Mother() && (GenPhotons->At(l)->Mother()->Status() == 3 || GenPhotons->At(l)->Mother()->Status() == 2)
              && (GenPhotons->At(l)->Mother()->AbsPdgId() == 11 
                  || GenPhotons->At(l)->Mother()->AbsPdgId() == 13
                  || GenPhotons->At(l)->Mother()->AbsPdgId() == 15)
            ) {
            CompositeParticle *object = new CompositeParticle();
            object->AddDaughter(GenPhotons->At(l));
            object->AddDaughter(GenPhotons->At(l)->Mother());
            if(object->Mass() > 1.0) GenPhotonMatch = kTRUE;
            delete object;          
          }      
        }
      }
    }


    if (GenEleMatch) isMCMatched += 1;
    if (GenMuMatch) isMCMatched += 2;
    if (GenTauMatch) isMCMatched += 4;
    if (GenPhotonMatch) isMCMatched += 8;

    
    //***********************************************************************************************
    //Match Electron to a genJet.
    //***********************************************************************************************  
    Double_t minDR = 5000.0;
    const GenJet *matchedGenJet = NULL;
//     if (fUseGen) {
//       for(UInt_t j=0;j<fGenJets->GetEntries();j++) {
// 	Double_t DR = MathUtils::DeltaR(elev[i]->Phi(), elev[i]->Eta(), fGenJets->At(j)->Phi(),
// 					fGenJets->At(j)->Eta());
// 	if (DR < minDR && DR < 0.5) {
// 	  minDR = DR;
// 	  matchedGenJet = fGenJets->At(j);
// 	}
//       }
//     }

    if (matchedGenJet) {
      if (abs(matchedGenJet->MatchedMCFlavor()) == 1) isMCMatched += 16;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 2) isMCMatched += 32;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 3) isMCMatched += 64;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 4) isMCMatched += 128;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 5) isMCMatched += 256;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 6) isMCMatched += 512;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 21) isMCMatched += 1024;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 0) {

	//Look explicitly for parton within cone around lepton
	Int_t tmpMatchedFlavor = -1;
	for (UInt_t j=0; j<fParticles->GetEntries(); ++j) {
	  //consider only partons
	  if (!(abs(fParticles->At(j)->PdgId()) <= 6 || fParticles->At(j)->PdgId() == 21)) continue;
	    
	  Double_t DR = MathUtils::DeltaR(elev[i]->Phi(), elev[i]->Eta(), fParticles->At(j)->Phi(),
					  fParticles->At(j)->Eta());
	  if (DR < 0.5) {
	    tmpMatchedFlavor = abs(fParticles->At(j)->PdgId());
	    break;
	  }
	}

	if (tmpMatchedFlavor == 1) isMCMatched += 16;
	else if (tmpMatchedFlavor == 2) isMCMatched += 32;
	else if (tmpMatchedFlavor == 3) isMCMatched += 64;
	else if (tmpMatchedFlavor == 4) isMCMatched += 128;
	else if (tmpMatchedFlavor == 5) isMCMatched += 256;
	else if (tmpMatchedFlavor == 6) isMCMatched += 512;
	else if (tmpMatchedFlavor == 21) isMCMatched += 1024;
	else if (tmpMatchedFlavor == -1) isMCMatched += 2048;
      }
    }
    
    //*****************************************************************
    //Debug MC Matching
    //*****************************************************************
//     if ((isMCMatched&2048) == 2048 && !(GenEleMatch || GenMuMatch || GenTauMatch || GenPhotonMatch)) {
//       cout << "weird electron : " << elev[i]->Pt() << " " << elev[i]->Eta() << " " << elev[i]->Phi() << " "
// 	   << endl;

//       cout << "GenTaus: " << GenTaus->GetEntries() << endl;
//       for (UInt_t l=0; l<GenTaus->GetEntries(); ++l) {        
// 	  cout << "GenTau " << l << " : " << GenTaus->At(l)->Pt() << " " << GenTaus->At(l)->Eta() << " " << GenTaus->At(l)->Phi() << " : " << MathUtils::DeltaR(*GenTaus->At(l), *(elev[i])) << endl;         
//       }

//       for(UInt_t j=0;j<fGenJets->GetEntries();j++) {
// 	Double_t DR = MathUtils::DeltaR(elev[i]->Phi(), elev[i]->Eta(), fGenJets->At(j)->Phi(),
// 					fGenJets->At(j)->Eta());
// 	cout << "GenJet " << j << " : " << fGenJets->At(j)->Pt() << " " << fGenJets->At(j)->Eta() << " " << fGenJets->At(j)->Phi() << " : " << DR << endl;
//       }
//       GeneratorTools::PrintNearbyParticles(fParticles, elev[i]->Eta(), elev[i]->Phi(), 0.5);
//       GeneratorTools::PrintHepMCTable(fParticles,kTRUE,-1);
//     }


    //********************************************************************************
    //Fill Electron
    //********************************************************************************
    FillElectron(elev[i], isMCMatched, primaryVertex);  // fill electron data object    

//      if (elev[i]->Pt() > 10 && !isMCMatched && (elev[i]->TrackIsolationDr03() + elev[i]->EcalRecHitIsoDr03() + elev[i]->HcalTowerSumEtDr03()) / elev[i]->Pt() < 0.1) {
//        cout << "ELE : " << elev[i]->Pt() << " " << elev[i]->Eta() << " " << elev[i]->Phi() << endl;
//        cout << elev[i]->TrackIsolationDr03() << " " << elev[i]->EcalRecHitIsoDr03() << " " << elev[i]->HcalTowerSumEtDr03() << " ";
//        cout << elev[i]->BestTrk()->D0Corrected(primaryVertex) << " " << elev[i]->BestTrk()->DzCorrected(primaryVertex);
//        cout << endl;
       
//        GeneratorTools::PrintNearbyParticles(fParticles, elev[i]->Eta(), elev[i]->Phi(), 0.5);
//        cout << "FULL GEN TABLE\n";
//        GeneratorTools::PrintHepMCTable(fParticles,kTRUE,-1);
//      }

  }
  
  
  //-----------------------------------------
  //
  // Loop through muons (and general tracks if desired).
  // If desired, tracks not matched to muon system information but passing 
  // kinematic preselection will be included in the muon TTree. 
  // Tracks from tracker info only will have typeBits=0
  //
  //-----------------------------------------
  vector<const Muon*> muonv;    // array of pointers to preselected muons ... 
  assert(fMuons);
  for(UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i); 
    if(!mu->HasTrk()) continue; 
    
    // Use tracker tracks for kinematics when available
    const Track *muTrk=0;
    if(mu->HasTrackerTrk())         { muTrk = mu->TrackerTrk(); }
    else if(mu->HasStandaloneTrk()) { muTrk = mu->StandaloneTrk(); } 
          
    if((muTrk->Pt()  < fMuonPtMin)  || (muTrk->Pt()  > fMuonPtMax))  continue;  // muon pT cut
    if((muTrk->Eta() < fMuonEtaMin) || (muTrk->Eta() > fMuonEtaMax)) continue;  // muon eta cut
    
    muonv.push_back(mu);        
  }

  for(UInt_t i=0; i<muonv.size(); i++) {    

    //********************************************************************************
    //Match MC Truth
    //********************************************************************************
    Int_t isMCMatched = 0;
    Bool_t GenEleMatch = kFALSE;
    Bool_t GenMuMatch = kFALSE;
    Bool_t GenTauMatch = kFALSE;
    if (GenLeptons) {
      for (UInt_t l=0; l<GenLeptons->GetEntries(); ++l) {
        if (MathUtils::DeltaR(*GenLeptons->At(l), *(muonv[i])) < 0.3 ) {
          if (GenLeptons->At(l)->AbsPdgId() == 11) {
            GenEleMatch = kTRUE;
          } 
          if (GenLeptons->At(l)->AbsPdgId() == 13) {
            GenMuMatch = kTRUE;
          }           
        }
      }
    }
    if (GenTaus) {
      for (UInt_t l=0; l<GenTaus->GetEntries(); ++l) {
        if (MathUtils::DeltaR(*GenTaus->At(l), *(muonv[i])) < 0.5 ) {
          GenTauMatch = kTRUE;
          break;
        }
      }
    }
    if (GenEleMatch) isMCMatched += 1;
    if (GenMuMatch) isMCMatched += 2;
    if (GenTauMatch) isMCMatched += 4;


    //***********************************************************************************************
    //Match Muon to a genJet.
    //***********************************************************************************************  
    Double_t minDR = 5000.0;
    const GenJet *matchedGenJet = NULL;
//     if (fUseGen) {
//       for(UInt_t j=0;j<fGenJets->GetEntries();j++) {
// 	Double_t DR = MathUtils::DeltaR(muonv[i]->Phi(), muonv[i]->Eta(), fGenJets->At(j)->Phi(),
// 					fGenJets->At(j)->Eta());
// 	if (DR < minDR && DR < 0.5) {
// 	  minDR = DR;
// 	  matchedGenJet = fGenJets->At(j);
// 	}
//       }
//     }

    if (matchedGenJet) {
      if (abs(matchedGenJet->MatchedMCFlavor()) == 1) isMCMatched += 16;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 2) isMCMatched += 32;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 3) isMCMatched += 64;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 4) isMCMatched += 128;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 5) isMCMatched += 256;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 6) isMCMatched += 512;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 21) isMCMatched += 1024;
      else if (abs(matchedGenJet->MatchedMCFlavor()) == 0) {

	//Look explicitly for parton within cone around lepton
	Int_t tmpMatchedFlavor = -1;
	for (UInt_t j=0; j<fParticles->GetEntries(); ++j) {
	  //consider only partons
	  if (!(abs(fParticles->At(j)->PdgId()) <= 6 || fParticles->At(j)->PdgId() == 21)) continue;
	    
	  Double_t DR = MathUtils::DeltaR(muonv[i]->Phi(), muonv[i]->Eta(), fParticles->At(j)->Phi(),
					  fParticles->At(j)->Eta());
	  if (DR < 0.5) {
	    tmpMatchedFlavor = abs(fParticles->At(j)->PdgId());
	    break;
	  }
	}

	if (tmpMatchedFlavor == 1) isMCMatched += 16;
	else if (tmpMatchedFlavor == 2) isMCMatched += 32;
	else if (tmpMatchedFlavor == 3) isMCMatched += 64;
	else if (tmpMatchedFlavor == 4) isMCMatched += 128;
	else if (tmpMatchedFlavor == 5) isMCMatched += 256;
	else if (tmpMatchedFlavor == 6) isMCMatched += 512;
	else if (tmpMatchedFlavor == 21) isMCMatched += 1024;
	else if (tmpMatchedFlavor == -1) isMCMatched += 2048;
      } 
    }

    //*****************************************************************
    //Debug MC Matching
    //*****************************************************************
//     if ((isMCMatched&2048) == 2048 && !(GenEleMatch || GenMuMatch || GenTauMatch )) {
//       cout << "weird muon: " << muonv[i]->Pt() << " " << muonv[i]->Eta() << " " << muonv[i]->Phi() << " "
// 	   << endl;

//       cout << "GenTaus: " << GenTaus->GetEntries() << endl;
//       for (UInt_t l=0; l<GenTaus->GetEntries(); ++l) {        
// 	  cout << "GenTau " << l << " : " << GenTaus->At(l)->Pt() << " " << GenTaus->At(l)->Eta() << " " << GenTaus->At(l)->Phi() << " : " << MathUtils::DeltaR(*GenTaus->At(l), *(muonv[i])) << endl;         
//       }

//       for(UInt_t j=0;j<fGenJets->GetEntries();j++) {
// 	Double_t DR = MathUtils::DeltaR(muonv[i]->Phi(), muonv[i]->Eta(), fGenJets->At(j)->Phi(),
// 					fGenJets->At(j)->Eta());
// 	cout << "GenJet " << j << " : " << fGenJets->At(j)->Pt() << " " << fGenJets->At(j)->Eta() << " " << fGenJets->At(j)->Phi() << " : " << DR << endl;
//       }
//       GeneratorTools::PrintNearbyParticles(fParticles, muonv[i]->Eta(), muonv[i]->Phi(), 0.5);
//       GeneratorTools::PrintHepMCTable(fParticles,kTRUE,-1);
//     }


    //********************************************************************************
    // fill muon data object
    //********************************************************************************
    FillMuon(muonv[i], isMCMatched, primaryVertex);

    if (fPrintDebug) {
      cout << "MUON : " << muonv[i]->Pt() << " " << muonv[i]->Eta() << " " << muonv[i]->Phi() << endl;
      cout << muonv[i]->IsoR03SumPt() << " " << muonv[i]->IsoR03EmEt() << " " << muonv[i]->IsoR03HadEt() << " ";
      
      const Track *muTrk=0;
      if(muonv[i]->HasTrackerTrk())         { muTrk = muonv[i]->TrackerTrk(); }
      else if(muonv[i]->HasStandaloneTrk()) { muTrk = muonv[i]->StandaloneTrk(); } 
      cout << muTrk->D0Corrected(*primaryVertex) << " " << muTrk->DzCorrected(*primaryVertex);
      cout << endl;
    }

//       GeneratorTools::PrintNearbyParticles(fParticles, muonv[i]->Eta(), muonv[i]->Phi(), 0.5);
//       cout << "FULL GEN TABLE\n";
//       GeneratorTools::PrintHepMCTable(fParticles,kTRUE,-1);
//     }

  }

  //-----------------------------------------
  //
  //Fill tracks for muon tag and probe
  //
  //-----------------------------------------
  assert(fTracks);
  for(UInt_t i=0; i<fTracks->GetEntries(); ++i) {
    const Track *track = fTracks->At(i);

    if((track->Pt()  < 9.0)  || (track->Pt()  > 7000.0))  continue;  // pT cut
    if((track->Eta() < fMuonEtaMin) || (track->Eta() > fMuonEtaMax)) continue;  // eta cut
      
    // Check that the track is not associated with a muon.
    // If it is, skip to next track...
    Bool_t isMuon = kFALSE;
    for(UInt_t j=0; j<fMuons->GetEntries(); ++j) {
      if(track == (fMuons->At(j)->TrackerTrk())) isMuon = kTRUE;
    }
    if(isMuon) continue;
  
    FillMuon(track, primaryVertex);
  }



  //-----------------------------------------
  //
  //Fill jets
  //
  //-----------------------------------------
  assert(CleanJetsNoPtCut);
  for(UInt_t i=0; i<CleanJetsNoPtCut->GetEntries(); ++i) {
    const PFJet *jet = (PFJet*)CleanJetsNoPtCut->At(i);

    if (fPrintDebug) {
      cout << "Jet " << i << " : " <<  CleanJetsNoPtCut->At(i)->Pt() << " : " << jet->Et() << " " << jet->Pt() << " " << jet->Eta() << " " << jet->Phi() << " : " 
           << jet->L1OffsetCorrectionScale() << " " << jet->L2RelativeCorrectionScale() << " " << jet->L3AbsoluteCorrectionScale() << " "  
           << jet->L4EMFCorrectionScale() << " " << jet->L5FlavorCorrectionScale() << " " << jet->L6LSBCorrectionScale() << " " << jet->L7PartonCorrectionScale() << " " 
           << jet->CombinedCorrectionScale() << " " 
           << " === " << jet->RawMom().Pt() << endl;
    }

    if(CleanJetsNoPtCut->At(i)->RawMom().Pt() > fJetPtMin ) { 
      FillJet(jet, primaryVertex);
    }
  }

  //-----------------------------------------
  //
  //Fill calojets
  //
  //-----------------------------------------
  assert(fCaloJets);
  for(UInt_t i=0; i<fCaloJets->GetEntries(); ++i) {
    const CaloJet *jet = (CaloJet*)fCaloJets->At(i);

    if (fPrintDebug) {
      cout << "CaloJet " << i << " : " <<  fCaloJets->At(i)->Pt() << " : " << jet->Et() << " " << jet->Pt() << " " << jet->Eta() << " " << jet->Phi() << " : " 
           << jet->L1OffsetCorrectionScale() << " " << jet->L2RelativeCorrectionScale() << " " << jet->L3AbsoluteCorrectionScale() << " "  
           << jet->L4EMFCorrectionScale() << " " << jet->L5FlavorCorrectionScale() << " " << jet->L6LSBCorrectionScale() << " " << jet->L7PartonCorrectionScale() << " " 
           << jet->CombinedCorrectionScale() << " " 
           << " === " << jet->RawMom().Pt() << endl;
    }

    if(fCaloJets->At(i)->Pt() > 30 ) { 
      FillCaloJet(jet);
    }
  }


  //-----------------------------------------
  //
  //Fill Photons
  //
  //-----------------------------------------
  for(UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
    const Photon *pho= fPhotons->At(i);
    //if(pho->HasPixelSeed()) continue;  // require NOT pixel seeded
    if(pho->Pt() >  10 && fabs(pho->Eta()) < 2.5) { FillPhoton(pho); }
  }


  // Find PFPhotons that are within dR < 0.5 of any filled electrons or muons
  vector<const PFCandidate*> FSRPhotons; 

  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    const PFCandidate *pf= fPFCandidates->At(i);
    Bool_t isFSRCandidate = kFALSE;

    //consider only pf photons
    if (!(pf->PFType() == mithep::PFCandidate::eGamma)) continue;

    for (UInt_t j=0; j<muonv.size(); ++j) {           
      if (MathUtils::DeltaR(*fPFCandidates->At(i), *muonv[j]) < 0.5) {
        isFSRCandidate = kTRUE;
      }
    }

    for (UInt_t j=0; j<elev.size(); ++j) {       
      if (MathUtils::DeltaR(*fPFCandidates->At(i), *elev[j]) < 0.5) {
        isFSRCandidate = kTRUE;
      }     
    }

    if (isFSRCandidate) FSRPhotons.push_back(pf);
  }



  //-----------------------------------------
  //
  //Fill PFCandidates within dR = 0.5 of filled electrons, muons
  //
  //-----------------------------------------

  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    const PFCandidate *pf= fPFCandidates->At(i);
    Bool_t doFillPFCandidate = kFALSE;
    //Bool_t doFillPFCandidate = kTRUE;
    UInt_t matchedObjectType = 0;
    UInt_t matchedObjectIndex = 0;

    for (UInt_t j=0; j<muonv.size(); ++j) {           
      if (MathUtils::DeltaR(*fPFCandidates->At(i), *muonv[j]) < 0.5) {
        doFillPFCandidate = kTRUE;
      }

      if(pf->TrackerTrk() && muonv[j]->TrackerTrk() &&
         pf->TrackerTrk() == muonv[j]->TrackerTrk()) {
        matchedObjectType = 13;
        matchedObjectIndex = j;
      }
    }

    if (matchedObjectType == 0) {
      for (UInt_t j=0; j<elev.size(); ++j) {       
        if (MathUtils::DeltaR(*fPFCandidates->At(i), *elev[j]) < 0.5) {
          doFillPFCandidate = kTRUE;
        }
        
        if(
          (pf->GsfTrk() && elev[j]->GsfTrk() &&
           pf->GsfTrk() == elev[j]->GsfTrk()) 
          ||
          (pf->TrackerTrk() && elev[j]->TrackerTrk() &&
           pf->TrackerTrk() == elev[j]->TrackerTrk())
          ||
          (pf->SCluster() && elev[j]->SCluster() && 
           pf->SCluster() == elev[j]->SCluster())
          ) {
          matchedObjectType = 11;
          matchedObjectIndex = j;
        }

      }
    }

    for (UInt_t j=0; j<FSRPhotons.size(); ++j) {       
      if (MathUtils::DeltaR(*fPFCandidates->At(i), *FSRPhotons[j]) < 0.5) {
        doFillPFCandidate = kTRUE;
      }
    }

    
    for (UInt_t j=0; j<fPhotons->GetEntries(); ++j) {       
      if (MathUtils::DeltaR(*fPFCandidates->At(i), *fPhotons->At(j)) < 0.5) {
        doFillPFCandidate = kTRUE;
      }
    }
    

    if (pf->PFType() == mithep::PFCandidate::eElectron || pf->PFType() == mithep::PFCandidate::eMuon) {
      doFillPFCandidate = kTRUE;
    }

    if (pf->Pt() > 5) {
      doFillPFCandidate = kTRUE;
    }
    if (doFillPFCandidate) {
      FillPFCandidate(pf, matchedObjectType, matchedObjectIndex, primaryVertex );
    }
  }

//debug
//   for (UInt_t j=0; j<fPhotons->GetEntries(); ++j) {       
//     if (fPhotons->At(j)->Pt() < 10) continue;
//     cout << "reco photon " << j << endl;
//     for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
//       if (MathUtils::DeltaR(*fPFCandidates->At(i), *fPhotons->At(j)) < 0.5) {
//         cout << "pfcand " << i << " : " << fPhotons->At(j)->Pt() << " " << fPhotons->At(j)->Eta() << " " << fPhotons->At(j)->Phi() << " : " << fPFCandidates->At(i)->Pt() << " " << fPFCandidates->At(i)->Eta() << " " << fPFCandidates->At(i)->Phi() << " : " << MathUtils::DeltaR(*fPFCandidates->At(i), *fPhotons->At(j)) << " : " << fPFCandidates->At(i)->PFType() <<  endl;
//       }
//     }
//   }



  //************************************************************************************************
  // Compute variations of MET
  //************************************************************************************************
  double MET_trk_X     = 0. ;
  double MET_trk_Y     = 0. ;
  double MET_trk_sumPt = 0. ;
  double MET_trkplusneu_X = 0. ;
  double MET_trkplusneu_Y = 0. ;
  double MET_trkplusneu_sumPt = 0. ;
  double MET_neu_noFwd_X = 0. ;
  double MET_neu_noFwd_Y = 0. ;
  double MET_neu_noFwd_sumPt = 0. ;

  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    // charged

    if ((fPFCandidates->At(i)->HasTrackerTrk() && fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*primaryVertex)) < 0.1) ||
        (fPFCandidates->At(i)->HasGsfTrk()     && fabs(fPFCandidates->At(i)->GsfTrk()->DzCorrected(*primaryVertex)    ) < 0.1)) {      
      MET_trk_X -= fPFCandidates->At(i)->Px();
      MET_trk_Y -= fPFCandidates->At(i)->Py();
      MET_trk_sumPt += fPFCandidates->At(i)->Pt();
      MET_trkplusneu_X -= fPFCandidates->At(i)->Px();
      MET_trkplusneu_Y -= fPFCandidates->At(i)->Py();
      MET_trkplusneu_sumPt += fPFCandidates->At(i)->Pt();
    }

    //neutral
    if (fPFCandidates->At(i)->PFType()== PFCandidate::eNeutralHadron || fPFCandidates->At(i)->PFType()== PFCandidate::eGamma) {
      //Not forward
      if (fabs(fPFCandidates->At(i)->Eta()) < 3.0 && fPFCandidates->At(i)->Pt() > 2.0) {
        MET_neu_noFwd_X -= fPFCandidates->At(i)->Px();
        MET_neu_noFwd_Y -= fPFCandidates->At(i)->Py();
        MET_neu_noFwd_sumPt += fPFCandidates->At(i)->Pt();
      }
      
      //trk+neutral
      if (fabs(fPFCandidates->At(i)->Eta()) < 5.0 && fPFCandidates->At(i)->Pt() > 2.0 ) {
        MET_trkplusneu_X -= fPFCandidates->At(i)->Px();
        MET_trkplusneu_Y -= fPFCandidates->At(i)->Py();
        MET_trkplusneu_sumPt += fPFCandidates->At(i)->Pt();
      }
    }

  }


//   cout << "OLD Met calculation: " << endl;
//   MetTools metTools(CleanMuons, CleanElectrons, fPFCandidates, primaryVertex, 0.1, 8.0, 5.0);
//   MetTools metToolsNoFwd(CleanMuons, CleanElectrons, fPFCandidates, primaryVertex, 0.1, 2.0, 3.0);


  //
  // Fill event info tree
  //
  fEventInfo.runNum       = GetEventHeader()->RunNum();
  fEventInfo.evtNum       = GetEventHeader()->EvtNum();
  fEventInfo.lumiSec      = GetEventHeader()->LumiSec();
  fEventInfo.eventweight  = 1.0;

  if (fReadPileupInfo) {
    if (fUseGen) {
      Int_t NPU = 0;
      Int_t NPU_PlusOne = 0;
      Int_t NPU_MinusOne = 0;
      for (UInt_t k=0; k < fPileupInfo->GetEntries() ; ++k) {
        if (fPileupInfo->At(k)->GetBunchCrossing() == 0)  NPU          = fPileupInfo->At(k)->GetPU_NumInteractions();
        if (fPileupInfo->At(k)->GetBunchCrossing() == 1)  NPU_PlusOne  = fPileupInfo->At(k)->GetPU_NumInteractions();
        if (fPileupInfo->At(k)->GetBunchCrossing() == -1) NPU_MinusOne = fPileupInfo->At(k)->GetPU_NumInteractions();
      }
      fEventInfo.nPUEvents   = NPU;
      fEventInfo.nPUMinusOne = NPU_MinusOne;
      fEventInfo.nPUPlusOne  = NPU_PlusOne;
    }
    fEventInfo.RhoKt6PFJets = fPileupEnergyDensity->At(0)->RhoKt6PFJets();
    fEventInfo.RhoKt6PFJetsForIso25 = fPileupEnergyDensity->At(0)->RhoKt6PFJetsForIso25();
    fEventInfo.RhoKt6PFJetsCentralNeutral = fPileupEnergyDensity->At(0)->RhoKt6PFJetsCentralNeutral();
    fEventInfo.RhoDeterministic = fPileupEnergyDensity->At(0)->Rho();
    fEventInfo.RhoDeterminaticLowEta = fPileupEnergyDensity->At(0)->RhoLowEta();
    fEventInfo.RhoRandom = fPileupEnergyDensity->At(0)->RhoRandom();
    fEventInfo.RhoRandomLowEta = fPileupEnergyDensity->At(0)->RhoRandomLowEta();
  } else {
    fEventInfo.nPUEvents = 0;
    fEventInfo.RhoKt6PFJets = 0;
    fEventInfo.RhoKt6PFJetsForIso25 = 0;
    fEventInfo.RhoKt6PFJetsCentralNeutral = 0;
    fEventInfo.RhoDeterministic = 0;
    fEventInfo.RhoDeterminaticLowEta = 0;
    fEventInfo.RhoRandom = 0;
    fEventInfo.RhoRandomLowEta = 0;
  }

  fEventInfo.triggerBits  = trigbits;
  fEventInfo.l1triggerBits  = l1trigbits;
  fEventInfo.hasGoodPV    = hasGoodPV;
  fEventInfo.nPV0         = fGoodPrimaryVertices->GetEntries();
  fEventInfo.pvx          = primaryVertex->X();
  fEventInfo.pvy          = primaryVertex->Y();
  fEventInfo.pvz          = primaryVertex->Z();
  fEventInfo.bsx          = bsx;
  fEventInfo.bsy          = bsy;
  fEventInfo.bsz          = bsz;
  fEventInfo.tcMEx        = fTCMet->At(0)->Mex();
  fEventInfo.tcMEy        = fTCMet->At(0)->Mey();
  fEventInfo.tcSumET      = fTCMet->At(0)->SumEt();
  fEventInfo.pfMEx        = fPFMet->At(0)->Mex();
  fEventInfo.pfMEy        = fPFMet->At(0)->Mey();
  fEventInfo.pfSumET      = fPFMet->At(0)->SumEt();
  fEventInfo.pfTrackMEx   = MET_trk_X;
  fEventInfo.pfTrackMEy   = MET_trk_Y;
  fEventInfo.pfTrackSumET = MET_trk_sumPt;
  fEventInfo.pfNeutralMEx   = MET_trkplusneu_X;
  fEventInfo.pfNeutralMEy   = MET_trkplusneu_Y;
  fEventInfo.pfNeutralSumET = MET_trkplusneu_sumPt;
  fEventInfo.pfNeutralNoFwdMEx   = MET_neu_noFwd_X;
  fEventInfo.pfNeutralNoFwdMEy   = MET_neu_noFwd_Y;
  fEventInfo.pfNeutralNoFwdSumET = MET_neu_noFwd_sumPt;

  fEventTree->Fill();


  //********************************************************************************************
  //Debug
  //********************************************************************************************

  if ( fPrintDebug   
//        ||
//       (GetEventHeader()->RunNum() == 1 && GetEventHeader()->LumiSec() == 127 && GetEventHeader()->EvtNum() == 90831)
//        ||
//     (GetEventHeader()->RunNum() == 148031 && GetEventHeader()->EvtNum() == 507254451)
//              ||
//     (GetEventHeader()->RunNum() == 148822 && GetEventHeader()->EvtNum() == 215107719)
//              ||
//     (GetEventHeader()->RunNum() == 148862 && GetEventHeader()->EvtNum() == 655863681)
//              ||
//     (GetEventHeader()->RunNum() == 149181 && GetEventHeader()->EvtNum() == 469493528)
//       (GetEventHeader()->RunNum() == 142933 && GetEventHeader()->LumiSec() == 926 && GetEventHeader()->EvtNum() == 565024903)
//       ||
//       (GetEventHeader()->RunNum() == 142933 && GetEventHeader()->LumiSec() == 939 && GetEventHeader()->EvtNum() == 572606305) 




    ) {

    cout << "DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    cout << "Event " << GetEventHeader()->RunNum() << " " << GetEventHeader()->LumiSec() << " " << GetEventHeader()->EvtNum() << "\n";
    
    cout << "Met : " << fPFMet->At(0)->Pt() << " " << fPFMet->At(0)->Mex() << " " << fPFMet->At(0)->Mey() << " : " << endl;
    cout << "TrackMet : " << TMath::Sqrt(pow(MET_trk_X,2) + pow(MET_trk_Y,2))  << " " << MET_trk_X << " " << MET_trk_Y << endl;
   

    for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {

//       const Electron *ele = fElectrons->At(i);
      cout << "All Electron" << i << " : " << fElectrons->At(i)->Pt() << " " << fElectrons->At(i)->Eta() << " " 
           << fElectrons->At(i)->Phi() << " "  
           << fElectrons->At(i)->IsEB() << " "  
           << fElectrons->At(i)->Charge() << " " 
           << " : " 
           << fElectrons->At(i)->CoviEtaiEta() << " " 
           << fElectrons->At(i)->HadronicOverEm() << " " 
           << fElectrons->At(i)->EcalRecHitIsoDr03() << " " 
           << fElectrons->At(i)->TrackIsolationDr03() << " " 
           << fElectrons->At(i)->HcalTowerSumEtDr03() << " " 
           << IsolationTools::PFElectronIsolation(fElectrons->At(i), fPFCandidates, primaryVertex, 0.1, 1.0, 0.4, 0.0) << " " 
           << fElectrons->At(i)->BestTrk()->NExpectedHitsInner() << " "
           << fElectrons->At(i)->SCluster()->Seed()->EMax() << " " 
           << fElectrons->At(i)->SCluster()->Seed()->E3x3() << " "
           << fElectrons->At(i)->FBrem() << " "
           << fElectrons->At(i)->ESuperClusterOverP() << " "
           << fElectrons->At(i)->SCluster()->Eta() << " : "
           << fElectrons->At(i)->BestTrk()->D0Corrected(*primaryVertex) << " "
           << fElectrons->At(i)->BestTrk()->DzCorrected(*primaryVertex) << " "
           << " : "
           << ElectronTools::PassCustomID(fElectrons->At(i), ElectronTools::kVBTFWorkingPointLowPtId) << " " 
           << endl;  
      
      
      for (UInt_t p=0; p<fPFCandidates->GetEntries(); ++p) {
        if (fPFCandidates->At(i)->HasTrackerTrk() || fPFCandidates->At(i)->HasGsfTrk()) {
          if ( (fElectrons->At(i)->TrackerTrk() == fPFCandidates->At(p)->TrackerTrk()) or
               (fElectrons->At(i)->HasGsfTrk() and fElectrons->At(i)->GsfTrk() == fPFCandidates->At(p)->GsfTrk()) ) {
            cout << "Matched PFCandidate: " << fPFCandidates->At(p)->Pt() << " " << fPFCandidates->At(p)->Eta() << " " 
                 << fPFCandidates->At(p)->Phi() << " : " << fPFCandidates->At(i)->HasTrackerTrk() << " " << fPFCandidates->At(i)->HasGsfTrk() 
                 << endl;
          }
        }
      }
      
    }

    for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {

//       const Muon *mu = fMuons->At(i);
      cout << "Muon" << i << " : " << fMuons->At(i)->Pt() << " " << fMuons->At(i)->Eta() << " " << fMuons->At(i)->Phi() << " "
	   << fMuons->At(i)->Charge() << " " 
	   << fMuons->At(i)->BestTrk()->NHits() << " " 
	   << fMuons->At(i)->NMatches() << " " 
	   << fMuons->At(i)->NSegments() << " " 
	   << fMuons->At(i)->BestTrk()->NPixelHits() << " " 
	   << fMuons->At(i)->Quality().Quality(MuonQuality::GlobalMuonPromptTight) << " " 
	   << fMuons->At(i)->BestTrk()->PtErr() << " " ;
      if (fMuons->At(i)->GlobalTrk()) {
        cout << fMuons->At(i)->GlobalTrk()->Chi2()/fMuons->At(i)->GlobalTrk()->Ndof() << " ";
      } 
      cout   << fMuons->At(i)->IsoR03SumPt() << " " 
             << fMuons->At(i)->IsoR03EmEt() << " " 
             << fMuons->At(i)->IsoR03HadEt() << " " 	  
             << fMuons->At(i)->BestTrk()->D0Corrected(*primaryVertex) << " "
             << fMuons->At(i)->BestTrk()->DzCorrected(*primaryVertex) << " "
             << endl;   

      for (UInt_t p=0; p<fPFCandidates->GetEntries(); ++p) {
        if (fPFCandidates->At(i)->HasTrackerTrk() || fPFCandidates->At(i)->HasGsfTrk()) {
          if ( (fMuons->At(i)->TrackerTrk() == fPFCandidates->At(p)->TrackerTrk()) ) {
            cout << "Matched PFCandidate: " << fPFCandidates->At(p)->Pt() << " " << fPFCandidates->At(p)->Eta() << " " 
                 << fPFCandidates->At(p)->Phi() << " : " << fPFCandidates->At(i)->HasTrackerTrk() << " " << fPFCandidates->At(i)->HasGsfTrk() 
                 << endl;
          }
        }
      }

     
    }
  }


}

//--------------------------------------------------------------------------------------------------
  void HwwNtuplerMod::FillGenInfo(const MCParticleCol *GenLeptons, const MCParticleCol *GenNeutrinos, 
                                  const MCParticleCol *GenBosons, const MCParticleCol *GenPhotons )
{
  // N.B. TGenInfo::nGenPart and TGenInfo::nGenCh is not filled here
  
//   cout << "Event " << GetEventHeader()->RunNum() << " " << GetEventHeader()->LumiSec() << " " << GetEventHeader()->EvtNum() << "\n";
//   GeneratorTools::PrintHepMCTable(fParticles, kTRUE, -1);

  fGenInfo.id_1     = fMCEvtInfo->Id1();
  fGenInfo.id_2     = fMCEvtInfo->Id2();
  fGenInfo.x_1      = fMCEvtInfo->X1();
  fGenInfo.x_2      = fMCEvtInfo->X2();
  fGenInfo.pdf_1    = fMCEvtInfo->Pdf1();
  fGenInfo.pdf_2    = fMCEvtInfo->Pdf2();
  fGenInfo.scalePdf = fMCEvtInfo->ScalePdf();
  fGenInfo.scale    = fMCEvtInfo->Scale();
  fGenInfo.weight   = fMCEvtInfo->Weight();
  
  fGenInfo.mass = 0;
  fGenInfo.pt = 0;
  fGenInfo.y = 0;
  fGenInfo.phi = 0;
  fGenInfo.pt_1 = 0;
  fGenInfo.eta_1 = 0;
  fGenInfo.phi_1 = 0;
  fGenInfo.pt_2 = 0;
  fGenInfo.eta_2 = 0;
  fGenInfo.phi_2 = 0;
  fGenInfo.met = 0;
  fGenInfo.vmass_1 = 0;
  fGenInfo.vpt_1 = 0;
  fGenInfo.vy_1 = 0;
  fGenInfo.veta_1 = 0;
  fGenInfo.vphi_1 = 0;
  fGenInfo.vmass_2 = 0;
  fGenInfo.vpt_2 = 0;
  fGenInfo.vy_2 = 0;
  fGenInfo.veta_2 = 0;
  fGenInfo.vphi_2 = 0;
  fGenInfo.ptBosonSystem = 0;

  if (GenLeptons->GetEntries() > 0) {
    fGenInfo.pt_1 = GenLeptons->At(0)->Pt();
    fGenInfo.eta_1 = GenLeptons->At(0)->Eta();
    fGenInfo.phi_1 = GenLeptons->At(0)->Phi();
   }
  if (GenLeptons->GetEntries() > 1) {
    fGenInfo.pt_2 = GenLeptons->At(1)->Pt();
    fGenInfo.eta_2 = GenLeptons->At(1)->Eta();
    fGenInfo.phi_2 = GenLeptons->At(1)->Phi();
  }

  if (GenNeutrinos->GetEntries() >= 2) {
    CompositeParticle *neutrinoSystem = new CompositeParticle;
    neutrinoSystem->AddDaughter(GenNeutrinos->At(0));
    neutrinoSystem->AddDaughter(GenNeutrinos->At(1));
    fGenInfo.met = neutrinoSystem->Pt();
    fGenInfo.metPhi = neutrinoSystem->Phi();
    delete neutrinoSystem;
  }

  const MCParticle *GenWBoson1 = 0;
  const MCParticle *GenWBoson2 = 0;
  const MCParticle *GenHiggsBoson = 0;
  for (UInt_t i=0; i< GenBosons->GetEntries(); ++i) {
    if (GenBosons->At(i)->PdgId() == 24) {
      if (!GenWBoson1) {
        GenWBoson1 = GenBosons->At(i);       
      } 
    }
    if (GenBosons->At(i)->PdgId() == -24) {
      if (!GenWBoson2) {
        GenWBoson2 = GenBosons->At(i);       
      }
    }
    if (GenBosons->At(i)->PdgId() == 25) {
      GenHiggsBoson = GenBosons->At(i);
    }
  }

//   if (GenWBoson1) {
//     cout << "W1 : " << GenWBoson1->Pt() << " " << GenWBoson1->Eta() << " " << GenWBoson1->Phi() << " : " << GenWBoson1->PdgId() << endl;
//   } 
//   if (GenWBoson2) {
//     cout << "W2 : " << GenWBoson2->Pt() << " " << GenWBoson2->Eta() << " " << GenWBoson2->Phi() << " : " << GenWBoson2->PdgId() << endl;
//   } 
//   if (GenHiggsBoson) {
//     cout << "Higgs : " << GenHiggsBoson->Pt() << " " << GenHiggsBoson->Eta() << " " << GenHiggsBoson->Phi() << " : " << GenHiggsBoson->PdgId() << endl;
//   } 


  if (GenWBoson1) {
    fGenInfo.vmass_1 = GenWBoson1->Mass();
    fGenInfo.vpt_1 = GenWBoson1->Pt();
    fGenInfo.vy_1 = GenWBoson1->Rapidity();
    fGenInfo.veta_1 = GenWBoson1->Eta();
    fGenInfo.vphi_1 = GenWBoson1->Phi();
  }
  if (GenWBoson2) {
    fGenInfo.vmass_2 = GenWBoson2->Mass();
    fGenInfo.vpt_2 = GenWBoson2->Pt();
    fGenInfo.vy_2 = GenWBoson2->Rapidity();
    fGenInfo.veta_2 = GenWBoson2->Eta();
    fGenInfo.vphi_2 = GenWBoson2->Phi();

    CompositeParticle *BosonSystem = new CompositeParticle;
    BosonSystem->AddDaughter(GenWBoson1);
    BosonSystem->AddDaughter(GenWBoson2);
//     cout << "WW System: " << BosonSystem->Pt() << " " << BosonSystem->Eta() << " " << BosonSystem->Phi() << endl;
//     fGenInfo.ptBosonSystem = BosonSystem->Pt();
//     fGenInfo.mass = BosonSystem->Mass();
//     fGenInfo.pt = BosonSystem->Pt();
//     fGenInfo.y = BosonSystem->Rapidity();
//     fGenInfo.phi = BosonSystem->Phi();
    fGenInfo.ptBosonSystem = BosonSystem->Pt();
    delete BosonSystem;
  }

  //take higgs pt from the higgs particle directly
  if (GenHiggsBoson) {
    fGenInfo.mass = GenHiggsBoson->Mass();
    fGenInfo.pt = GenHiggsBoson->Pt();
    fGenInfo.y = GenHiggsBoson->Rapidity();
    fGenInfo.phi = GenHiggsBoson->Phi();
  }


//   if (GenLeptons->GetEntries() > 0) 
//     cout << "Lepton1: " << fGenInfo.pt_1 << " " << fGenInfo.eta_1 << " " << fGenInfo.phi_1 << endl;
//   if (GenLeptons->GetEntries() > 1)
//     cout << "Lepton2: " << fGenInfo.pt_2 << " " << fGenInfo.eta_2 << " " << fGenInfo.phi_2 << endl;
//   cout << "BosonSystemPt: " << fGenInfo.ptBosonSystem  << endl;


  if (fGenJets) {
    Int_t nGenJets = 0;
    Int_t genJetsIndex = 0;
    fGenInfo.jetpt_1 = 0;
    fGenInfo.jetpt_2 = 0;
    fGenInfo.jetpt_3 = 0;
    fGenInfo.jetpt_4 = 0;
    for (UInt_t i=0; i<fGenJets->GetEntries(); ++i) {
      
      if (fabs(fGenJets->At(i)->Eta()) < 5.0 && fGenJets->At(i)->Pt() > 10.0  
          && (GenLeptons->GetEntries() < 1 || MathUtils::DeltaR(*fGenJets->At(i), *GenLeptons->At(0)) > 0.3)
          && (GenLeptons->GetEntries() < 2 || MathUtils::DeltaR(*fGenJets->At(i), *GenLeptons->At(1)) > 0.3)
        ) {
//              cout << "GenJet " << i << " : " << fGenJets->At(i)->Pt() << " " << fGenJets->At(i)->Eta() << " " << fGenJets->At(i)->Phi() << endl;

        if (genJetsIndex == 0) {
          fGenInfo.jetpt_1 = fGenJets->At(i)->Pt();
 //          for (UInt_t k=0; k < fParticles->GetEntries() ; ++k) {
//             if (MathUtils::DeltaR(*fGenJets->At(i), *fParticles->At(k)) < 0.5) {
//               cout << "GenJet Particle: " << k << " : " << fParticles->At(k)->PdgId() << " " << fParticles->At(k)->Status() << " " << fParticles->At(k)->Pt() << " " << fParticles->At(k)->Eta() << " " << fParticles->At(k)->Phi() << endl;
//             }
//           }
        }
        if (genJetsIndex == 1) fGenInfo.jetpt_2 = fGenJets->At(i)->Pt();
        if (genJetsIndex == 2) fGenInfo.jetpt_3 = fGenJets->At(i)->Pt();
        if (genJetsIndex == 3) fGenInfo.jetpt_4 = fGenJets->At(i)->Pt();
        if (fGenJets->At(i)->Pt() > 25.0) nGenJets++;
        genJetsIndex++;
      }
    }
    fGenInfo.nGenJets = nGenJets;

  }

  //***************************************************************
  //Fill Gen Particles
  //***************************************************************


  //GenLeptons
  for (UInt_t i=0; i<GenLeptons->GetEntries(); ++i) {
    TClonesArray &rGenParticleArr = *fGenParticleArr;
    assert(rGenParticleArr.GetEntries() < rGenParticleArr.GetSize());
    const Int_t index = rGenParticleArr.GetEntries();  
    new(rGenParticleArr[index]) TGenParticle();
    TGenParticle *pGenParticle = (TGenParticle*)rGenParticleArr[index];

    pGenParticle->pt	  = GenLeptons->At(i)->Pt(); 
    pGenParticle->eta  	  = GenLeptons->At(i)->Eta();
    pGenParticle->phi  	  = GenLeptons->At(i)->Phi();
    pGenParticle->mass 	  = GenLeptons->At(i)->Mass();
    pGenParticle->status  = GenLeptons->At(i)->Status();
    pGenParticle->pdgid	  = GenLeptons->At(i)->PdgId();
    pGenParticle->motherPdgID  = 0;    
    if (GenLeptons->At(i)->Mother()) {
      pGenParticle->motherPdgID  = GenLeptons->At(i)->Mother()->PdgId();
    }
  }

  //Status == 3 leptons
  for (UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    TClonesArray &rGenParticleArr = *fGenParticleArr;

    if (fParticles->At(i)->IsLepton() && fParticles->At(i)->Status() == 3) {
      assert(rGenParticleArr.GetEntries() < rGenParticleArr.GetSize());
      const Int_t index = rGenParticleArr.GetEntries();  
      new(rGenParticleArr[index]) TGenParticle();
      TGenParticle *pGenParticle = (TGenParticle*)rGenParticleArr[index];

      pGenParticle->pt	          = fParticles->At(i)->Pt(); 
      pGenParticle->eta  	  = fParticles->At(i)->Eta();
      pGenParticle->phi  	  = fParticles->At(i)->Phi();
      pGenParticle->mass 	  = fParticles->At(i)->Mass();
      pGenParticle->status        = fParticles->At(i)->Status();
      pGenParticle->pdgid	  = fParticles->At(i)->PdgId();
      pGenParticle->motherPdgID   = 0;    
      if (fParticles->At(i)->Mother()) {
        pGenParticle->motherPdgID  = fParticles->At(i)->Mother()->PdgId();
      }
    }

  }


  //GenNeutrinos
  for (UInt_t i=0; i<GenNeutrinos->GetEntries(); ++i) {
    TClonesArray &rGenParticleArr = *fGenParticleArr;
    assert(rGenParticleArr.GetEntries() < rGenParticleArr.GetSize());
    const Int_t index = rGenParticleArr.GetEntries();  
    new(rGenParticleArr[index]) TGenParticle();
    TGenParticle *pGenParticle = (TGenParticle*)rGenParticleArr[index];
    
    pGenParticle->pt	  = GenNeutrinos->At(i)->Et(); 
    pGenParticle->eta  	  = GenNeutrinos->At(i)->Eta();
    pGenParticle->phi  	  = GenNeutrinos->At(i)->Phi();
    pGenParticle->mass 	  = GenNeutrinos->At(i)->Mass();
    pGenParticle->status  = GenNeutrinos->At(i)->Status();
    pGenParticle->pdgid	  = GenNeutrinos->At(i)->PdgId();
    pGenParticle->motherPdgID  = 0;
    if (GenNeutrinos->At(i)->Mother()) {
      pGenParticle->motherPdgID  = GenNeutrinos->At(i)->Mother()->PdgId();
    }
  }

  //GenPhotons
  for (UInt_t i=0; i<GenPhotons->GetEntries(); ++i) {

    Bool_t FillThisGenPhoton = kFALSE;
    if (GenPhotons->At(i)->Pt() > 5.0 && GenPhotons->At(i)->AbsEta() < 2.7) FillThisGenPhoton = kTRUE;
    if( GenPhotons->At(i)->Is(MCParticle::kGamma) && GenPhotons->At(i)->Status() == 1 && 
        GenPhotons->At(i)->DistinctMother() && GenPhotons->At(i)->DistinctMother()->Status() == 3 &&
        (GenPhotons->At(i)->DistinctMother()->Is(MCParticle::kEl)  || GenPhotons->At(i)->DistinctMother()->Is(MCParticle::kMu) ||
         GenPhotons->At(i)->DistinctMother()->Is(MCParticle::kTau) || GenPhotons->At(i)->DistinctMother()->Is(MCParticle::kW))
      ) {
      FillThisGenPhoton = kTRUE;
    }
    if (!FillThisGenPhoton) continue;    

    TClonesArray &rGenParticleArr = *fGenParticleArr;
    assert(rGenParticleArr.GetEntries() < rGenParticleArr.GetSize());
    const Int_t index = rGenParticleArr.GetEntries();  
    new(rGenParticleArr[index]) TGenParticle();
    TGenParticle *pGenParticle = (TGenParticle*)rGenParticleArr[index];
    
    pGenParticle->pt	  = GenPhotons->At(i)->Et(); 
    pGenParticle->eta  	  = GenPhotons->At(i)->Eta();
    pGenParticle->phi  	  = GenPhotons->At(i)->Phi();
    pGenParticle->mass 	  = GenPhotons->At(i)->Mass();
    pGenParticle->status  = GenPhotons->At(i)->Status();
    pGenParticle->pdgid	  = GenPhotons->At(i)->PdgId();
    pGenParticle->motherPdgID  = 0;
    if (GenPhotons->At(i)->DistinctMother()) {
      pGenParticle->motherPdgID  = GenPhotons->At(i)->DistinctMother()->PdgId();
    }
  }

  //GenBosons
  for (UInt_t i=0; i<GenBosons->GetEntries(); ++i) {
    TClonesArray &rGenParticleArr = *fGenParticleArr;
    assert(rGenParticleArr.GetEntries() < rGenParticleArr.GetSize());
    const Int_t index = rGenParticleArr.GetEntries();  
    new(rGenParticleArr[index]) TGenParticle();
    TGenParticle *pGenParticle = (TGenParticle*)rGenParticleArr[index];
    
    pGenParticle->pt	  = GenBosons->At(i)->Et(); 
    pGenParticle->eta  	  = GenBosons->At(i)->Eta();
    pGenParticle->phi  	  = GenBosons->At(i)->Phi();
    pGenParticle->mass 	  = GenBosons->At(i)->Mass();
    pGenParticle->status  = GenBosons->At(i)->Status();
    pGenParticle->pdgid	  = GenBosons->At(i)->PdgId();
    pGenParticle->motherPdgID  = 0;
    if (GenBosons->At(i)->Mother()) {
      pGenParticle->motherPdgID  = GenBosons->At(i)->Mother()->PdgId();
    }
  }
  
}


//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::FillElectron(const Electron *ele, Int_t isMCMatched, const Vertex* primaryVertex
  )
{
  assert(ele);
 
  const PFCandidate *pfEleMatch = 0;

  for (UInt_t p=0; p<fPFCandidates->GetEntries();p++) {   
    const PFCandidate *pf = fPFCandidates->At(p);      
    //find matching PFElectron
    if (!pfEleMatch) {
      if(
        (pf->GsfTrk() && ele->GsfTrk() &&
         pf->GsfTrk() == ele->GsfTrk())
        ||
        (pf->TrackerTrk() && ele->TrackerTrk() &&
         pf->TrackerTrk() == ele->TrackerTrk())
        ) {
        pfEleMatch = pf;
      }
    }
  }

  TClonesArray &rElectronArr = *fElectronArr;
  assert(rElectronArr.GetEntries() < rElectronArr.GetSize());
  const Int_t index = rElectronArr.GetEntries();  
  new(rElectronArr[index]) TElectron();
  TElectron *pElectron = (TElectron*)rElectronArr[index];

  pElectron->pt              = ele->Pt();
  pElectron->eta             = ele->Eta();
  pElectron->phi             = ele->Phi();
  pElectron->p               = ele->BestTrk()->P();
  if (pfEleMatch) {
    pElectron->pfPt              = pfEleMatch->Pt();
    pElectron->pfEta             = pfEleMatch->Eta();
    pElectron->pfPhi             = pfEleMatch->Phi();
  } else {
    pElectron->pfPt              = -999;
    pElectron->pfEta             = 0;
    pElectron->pfPhi             = 0;
  }
  pElectron->scEt            = ele->SCluster()->Et();
  pElectron->scEta           = ele->SCluster()->Eta();
  pElectron->scPhi           = ele->SCluster()->Phi();
  pElectron->q               = ele->Charge();
  pElectron->isMCReal        = isMCMatched;
  pElectron->hltMatchBits    = MatchHLT(ele->SCluster()->Eta(),ele->SCluster()->Phi(), GetEventHeader()->RunNum(), GetEventHeader()->EvtNum());  
  pElectron->l1TriggerMatchBits = MatchL1(ele->SCluster()->Eta(),ele->SCluster()->Phi());
  pElectron->isEcalDriven    = ele->IsEcalDriven();
  pElectron->isTrackerDriven = ele->IsTrackerDriven();
  pElectron->isEB            = ele->IsEB();

  pElectron->d0              = ele->BestTrk()->D0Corrected(*primaryVertex);
  pElectron->d0Err           = ele->BestTrk()->D0Err();
  pElectron->dz              = ele->BestTrk()->DzCorrected(*primaryVertex);
  pElectron->ip3d            = ele->Ip3dPV();
  pElectron->ip3dSig         = ele->Ip3dPVSignificance();

  //***************************************************************
  //Conversion Veto
  //***************************************************************
  Int_t passConvVeto = 0;
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 1, 1e-6, 2.0, kFALSE, kTRUE)) {
    passConvVeto += 1;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 1, 1e-6, 2.0, kFALSE, kFALSE)) {
    passConvVeto += 2;
  }
  //attempts matching to all conversions without regard to double-counting)
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 1, 1e-6, 2.0, kTRUE, kFALSE)) {
  //also attempt matching through ckf track in addition to gsf track
    passConvVeto += 4;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 1, 1e-6, -1.0, kTRUE, kFALSE)) {
  //relax lxy cut
    passConvVeto += 8;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 1, 1e-6, -999.9, kTRUE, kFALSE)) {
  //remove lxy cut
    passConvVeto += 16;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 1, 1e-10, -999.9, kTRUE, kFALSE)) {
  //relax probability cut
    passConvVeto += 32;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 2, 1e-10, -999.9, kTRUE, kFALSE)) {
  //relax hits before vtx cut
    passConvVeto += 64;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 999, 1e-10, -999.9, kTRUE, kFALSE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 128;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 0, 1e-6, 2.0, kFALSE, kTRUE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 256;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 0, 1e-6, 2.0, kFALSE, kFALSE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 512;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 0, 1e-6, 2.0, kTRUE, kFALSE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 1024;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 0, 1e-6, -1.0, kTRUE, kFALSE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 2048;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 0, 1e-6, -999.9, kTRUE, kFALSE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 4096;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 0, 1e-10, -999.9, kTRUE, kFALSE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 8192;
  }
  pElectron->isConv          = passConvVeto;

  pElectron->nExpHitsInner   = ele->CorrectedNExpectedHitsInner();
  pElectron->partnerDeltaCot = ele->ConvPartnerDCotTheta();
  pElectron->partnerDist     = ele->ConvPartnerDist();
  pElectron->partnerRadius   = ele->ConvPartnerRadius();  

  pElectron->nBrem           = ele->NumberOfClusters() - 1;
  pElectron->fBrem           = ele->FBrem();
  pElectron->EOverP          = ele->ESuperClusterOverP();
  pElectron->pIn             = ele->PIn();
  pElectron->ESeedClusterOverPIn  = ele->ESeedClusterOverPIn();
  pElectron->ESeedClusterOverPout = ele->ESeedClusterOverPout();
  pElectron->EEleClusterOverPout  = ele->EEleClusterOverPout();
  pElectron->EcalEnergy      = ele->EcalEnergy();
  pElectron->EcalEnergyError = ele->EcalEnergyError();

  pElectron->deltaEtaIn      = ele->DeltaEtaSuperClusterTrackAtVtx();
  pElectron->deltaPhiIn      = ele->DeltaPhiSuperClusterTrackAtVtx();
  pElectron->dEtaCalo        = ele->DeltaEtaSeedClusterTrackAtCalo();
  pElectron->dPhiCalo        = ele->DeltaPhiSeedClusterTrackAtCalo();
  pElectron->sigiEtaiEta     = ele->CoviEtaiEta();
  pElectron->sigiPhiiPhi     = sqrt(ele->SCluster()->Seed()->CoviPhiiPhi());
  if (isnan(pElectron->sigiPhiiPhi))  pElectron->sigiPhiiPhi = 0.0;

  pElectron->CovIEtaIPhi = ele->SCluster()->Seed()->CoviEtaiPhi();
  pElectron->SCEtaWidth = ele->SCluster()->EtaWidth();
  pElectron->SCPhiWidth = ele->SCluster()->PhiWidth();
  pElectron->R9 = ele->SCluster()->R9();
  pElectron->PreShowerOverRaw = ele->SCluster()->PreshowerEnergy() / ele->SCluster()->RawEnergy();

  pElectron->HoverE          = ele->HadronicOverEm();
  pElectron->GsfTrackChi2OverNdof = ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof();
  if (ele->TrackerTrk()) {
    pElectron->KFTrackChi2OverNdof = ele->TrackerTrk()->RChi2();
    pElectron->KFTrackNHits = ele->TrackerTrk()->NHits();
    pElectron->KFTrackNLayersWithMeasurement = ele->CTFTrkNLayersWithMeasurement();
  } else {
    pElectron->KFTrackChi2OverNdof = 0;
    pElectron->KFTrackNHits = -1;
    pElectron->KFTrackNLayersWithMeasurement = -1;    
  }

  if (ele->SCluster()->Seed()->E5x5() > 0) {
    pElectron->SeedE1x5OverE5x5 = ele->SCluster()->Seed()->E1x5() / ele->SCluster()->Seed()->E5x5();
  } else {
    pElectron->SeedE1x5OverE5x5 = -1.0;
  }

  pElectron->likelihood      = ele->IDLikelihood();
  pElectron->mva             = ele->Mva();

  pElectron->trkIso03        = ele->TrackIsolationDr03();
  pElectron->emIso03         = ele->EcalRecHitIsoDr03();
  pElectron->hadIso03        = ele->HcalTowerSumEtDr03();
  pElectron->trkIso04        = ele->TrackIsolationDr04();
  pElectron->emIso04         = ele->EcalRecHitIsoDr04();
  pElectron->hadIso04        = ele->HcalTowerSumEtDr04();

  //Fill Regression Variables
  pElectron->isEE = ele->IsEE();
  
  pElectron->SCRawEnergy = ele->SCluster()->RawEnergy();
  pElectron->E5x5 = ele->E55();
  pElectron->EtaSeed = ele->SCluster()->Seed()->Eta();
  pElectron->PhiSeed = ele->SCluster()->Seed()->Phi();
  pElectron->ESeed   = ele->SCluster()->Seed()->Energy();
  pElectron->E3x3Seed = ele->SCluster()->Seed()->E3x3();
  pElectron->E5x5Seed = ele->SCluster()->Seed()->E5x5();
  pElectron->EMaxSeed = ele->SCluster()->Seed()->EMax();
  pElectron->E2ndSeed = ele->SCluster()->Seed()->E2nd();
  pElectron->ETopSeed= ele->SCluster()->Seed()->ETop();
  pElectron->EBottomSeed = ele->SCluster()->Seed()->EBottom();
  pElectron->ELeftSeed = ele->SCluster()->Seed()->ELeft();
  pElectron->ERightSeed = ele->SCluster()->Seed()->ERight();
  pElectron->E2x5MaxSeed = ele->SCluster()->Seed()->E2x5Max();
  pElectron->E2x5TopSeed = ele->SCluster()->Seed()->E2x5Top();
  pElectron->E2x5BottomSeed = ele->SCluster()->Seed()->E2x5Bottom();
  pElectron->E2x5LeftSeed = ele->SCluster()->Seed()->E2x5Left();
  pElectron->E2x5RightSeed = ele->SCluster()->Seed()->E2x5Right();
  pElectron->IEtaSeed = ele->SCluster()->Seed()->IEta();
  pElectron->IPhiSeed = ele->SCluster()->Seed()->IPhi();
  pElectron->EtaCrySeed = ele->SCluster()->Seed()->EtaCry();
  pElectron->PhiCrySeed = ele->SCluster()->Seed()->PhiCry();
  pElectron->TrackMomentumError = ele->TrackMomentumError();
  pElectron->Classification = ele->Classification();

}
    
  //--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::FillMuon(const Muon *mu, Int_t isMCMatched, const Vertex* primaryVertex)
{
  // N.B. TMuon::isUsed is not to be filled here
  
  assert(mu);

  const PFCandidate *pfMuonMatched = 0;
  for (UInt_t p=0; p<fPFCandidates->GetEntries();p++) {   
    const PFCandidate *pf = fPFCandidates->At(p);
      
    if (!pfMuonMatched) {
      if(pf->TrackerTrk() && mu->TrackerTrk() &&
         pf->TrackerTrk() == mu->TrackerTrk()) {
        pfMuonMatched = pf;
      }           
    }
  }

  TClonesArray &rMuonArr = *fMuonArr;
  assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
  const Int_t index = rMuonArr.GetEntries();  
  new(rMuonArr[index]) TMuon();
  TMuon *pMuon = (TMuon*)rMuonArr[index];
  
  // Use tracker track when available
  const Track *muTrk=0;
  if(mu->HasTrackerTrk())         { muTrk = mu->TrackerTrk();    }
  else if(mu->HasStandaloneTrk()) { muTrk = mu->StandaloneTrk(); } 
  assert(muTrk);                 
  
  pMuon->pt        = muTrk->Pt();
  pMuon->eta       = muTrk->Eta();
  pMuon->phi       = muTrk->Phi();
  if (pfMuonMatched) {
    pMuon->pfPt        = pfMuonMatched->Pt();
    pMuon->pfEta       = pfMuonMatched->Eta();
    pMuon->pfPhi       = pfMuonMatched->Phi();  
  } else {
    pMuon->pfPt        = -999;
    pMuon->pfEta       = 0;
    pMuon->pfPhi       = 0;
  }
  pMuon->pterr     = muTrk->PtErr();
  pMuon->trkIso03  = mu->IsoR03SumPt();
  pMuon->emIso03   = mu->IsoR03EmEt();
  pMuon->hadIso03  = mu->IsoR03HadEt();
  pMuon->hoIso03   = mu->IsoR03HoEt();
  pMuon->trkIso05  = mu->IsoR05SumPt();
  pMuon->emIso05   = mu->IsoR05EmEt();
  pMuon->hadIso05  = mu->IsoR05HadEt();
  pMuon->hoIso05   = mu->IsoR05HoEt();

  pMuon->d0        = muTrk->D0Corrected(*primaryVertex);
  pMuon->d0Err     = muTrk->D0Err();
  pMuon->dz        = muTrk->DzCorrected(*primaryVertex);
  pMuon->tkNchi2   = muTrk->RChi2();
  pMuon->trkLayers = mu->NTrkLayersHit();

  if(mu->HasGlobalTrk())          { pMuon->muNchi2 = mu->GlobalTrk()->RChi2();     }
  else if(mu->HasStandaloneTrk()) { pMuon->muNchi2 = mu->StandaloneTrk()->RChi2(); }
  else if(mu->HasTrackerTrk())    { pMuon->muNchi2 = mu->TrackerTrk()->RChi2();    }
  
  pMuon->q          = muTrk->Charge();
  pMuon->nValidHits = mu->NValidHits();
  
  pMuon->qualityBits = mu->Quality().QualityMask().Mask();

  //
  // NOTE:
  // It is possible for a muon to be TK+SA. The muon reco associates a TK with a SA if
  // chamber matches for the TK and hits for the SA share DetIDs
  //   (see hypernews thread: https://hypernews.cern.ch/HyperNews/CMS/get/csa08-muons/57/2/1/1/1.html)
  //	      
  pMuon->typeBits = 0;
  if(mu->IsGlobalMuon())     { pMuon->typeBits |= kGlobal; }
  if(mu->IsTrackerMuon())    { pMuon->typeBits |= kTracker; }
  if(mu->IsStandaloneMuon()) { pMuon->typeBits |= kStandalone; }
  
  pMuon->nTkHits      = muTrk->NHits();
  pMuon->nPixHits     = muTrk->NPixelHits();
  pMuon->nSeg         = mu->NSegments();
  pMuon->nMatch       = mu->NMatches();
  pMuon->hltMatchBits = MatchHLTMuon(muTrk->Pt(),muTrk->Eta(),muTrk->Phi(),GetEventHeader()->RunNum(), GetEventHeader()->EvtNum() );
  
  pMuon->staPt  = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Pt()  : -1;  
  pMuon->staEta = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Eta() : -999;
  pMuon->staPhi = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Phi() : -999;

  pMuon->isMCReal = isMCMatched;
  pMuon->ip3d              = mu->Ip3dPV();
  pMuon->ip3dSig           = mu->Ip3dPVSignificance();

  pMuon->SegmentCompatibility     = fMuonTools->GetSegmentCompatability(mu);
  pMuon->CaloCompatilibity        = fMuonTools->GetCaloCompatability(mu, kTRUE, kTRUE);

  pMuon->TrkKink           = mu->TrkKink();
  pMuon->GlobalKink        = mu->GlbKink();
  pMuon->HadEnergy         = mu->HadEnergy();
  pMuon->HadS9Energy       = mu->HadS9Energy();
  pMuon->HoEnergy          = mu->HoEnergy();
  pMuon->HoS9Energy        = mu->HoS9Energy();
  pMuon->EmEnergy          = mu->EmEnergy();
  pMuon->EmS9Energy        = mu->EmS9Energy();


  //***************************************
  // Check for pf ID
  //***************************************
  Double_t PFMuonEEcal = 0;
  Double_t PFMuonEtEcal = 0;
  pMuon->PassPFId = kFALSE;
  for (UInt_t p=0; p<fPFCandidates->GetEntries();p++) {   
    const PFCandidate *pf = fPFCandidates->At(p);
    if ( (pf->TrackerTrk() == mu->TrackerTrk()) && abs(pf->PFType()) == mithep::PFCandidate::eMuon ) {
      pMuon->PassPFId = kTRUE;

      PFMuonEEcal = pf->EECal();
      //compute EtEcal for Type2 FSR recovery
      if (pf->E() > 0) {
        PFMuonEtEcal = pf->EECal()*pf->Pt()/pf->E();
      } else {
        PFMuonEtEcal = pf->EECal()/TMath::CosH(pf->Eta());
      }
      break;
    }
  }
  pMuon->PFMuonEEcal = PFMuonEEcal;
  pMuon->PFMuonEtEcal = PFMuonEtEcal;


  //***************************************
  // Turn this stuff off for now
  //***************************************
//   pMuon->Station0NSegments    = mu->GetNSegments(0);
//   pMuon->Station0TrackDist    = mu->GetTrackDist(0);
//   pMuon->Station0TrackDistErr = mu->GetTrackDistErr(0);
//   pMuon->Station0dX           = mu->GetDX(0);
//   pMuon->Station0dY           = mu->GetDY(0);
//   pMuon->Station0PullX        = mu->GetPullX(0);
//   pMuon->Station0PullY        = mu->GetPullY(0);

//   pMuon->Station1NSegments    = mu->GetNSegments(1);
//   pMuon->Station1TrackDist    = mu->GetTrackDist(1);
//   pMuon->Station1TrackDistErr = mu->GetTrackDistErr(1);
//   pMuon->Station1dX           = mu->GetDX(1);
//   pMuon->Station1dY           = mu->GetDY(1);
//   pMuon->Station1PullX        = mu->GetPullX(1);
//   pMuon->Station1PullY        = mu->GetPullY(1);

//   pMuon->Station2NSegments    = mu->GetNSegments(2);
//   pMuon->Station2TrackDist    = mu->GetTrackDist(2);
//   pMuon->Station2TrackDistErr = mu->GetTrackDistErr(2);
//   pMuon->Station2dX           = mu->GetDX(2);
//   pMuon->Station2dY           = mu->GetDY(2);
//   pMuon->Station2PullX        = mu->GetPullX(2);
//   pMuon->Station2PullY        = mu->GetPullY(2);

//   pMuon->Station3NSegments    = mu->GetNSegments(3);
//   pMuon->Station3TrackDist    = mu->GetTrackDist(3);
//   pMuon->Station3TrackDistErr = mu->GetTrackDistErr(3);
//   pMuon->Station3dX           = mu->GetDX(3);
//   pMuon->Station3dY           = mu->GetDY(3);
//   pMuon->Station3PullX        = mu->GetPullX(3);
//   pMuon->Station3PullY        = mu->GetPullY(3);

//   pMuon->Station4NSegments    = mu->GetNSegments(4);
//   pMuon->Station4TrackDist    = mu->GetTrackDist(4);
//   pMuon->Station4TrackDistErr = mu->GetTrackDistErr(4);
//   pMuon->Station4dX           = mu->GetDX(4);
//   pMuon->Station4dY           = mu->GetDY(4);
//   pMuon->Station4PullX        = mu->GetPullX(4);
//   pMuon->Station4PullY        = mu->GetPullY(4);

//   pMuon->Station5NSegments    = mu->GetNSegments(5);
//   pMuon->Station5TrackDist    = mu->GetTrackDist(5);
//   pMuon->Station5TrackDistErr = mu->GetTrackDistErr(5);
//   pMuon->Station5dX           = mu->GetDX(5);
//   pMuon->Station5dY           = mu->GetDY(5);
//   pMuon->Station5PullX        = mu->GetPullX(5);
//   pMuon->Station5PullY        = mu->GetPullY(5);

//   pMuon->Station6NSegments    = mu->GetNSegments(6);
//   pMuon->Station6TrackDist    = mu->GetTrackDist(6);
//   pMuon->Station6TrackDistErr = mu->GetTrackDistErr(6);
//   pMuon->Station6dX           = mu->GetDX(6);
//   pMuon->Station6dY           = mu->GetDY(6);
//   pMuon->Station6PullX        = mu->GetPullX(6);
//   pMuon->Station6PullY        = mu->GetPullY(6);

//   pMuon->Station7NSegments    = mu->GetNSegments(7);
//   pMuon->Station7TrackDist    = mu->GetTrackDist(7);
//   pMuon->Station7TrackDistErr = mu->GetTrackDistErr(7);
//   pMuon->Station7dX           = mu->GetDX(7);
//   pMuon->Station7dY           = mu->GetDY(7);
//   pMuon->Station7PullX        = mu->GetPullX(7);
//   pMuon->Station7PullY        = mu->GetPullY(7);


}




void HwwNtuplerMod::FillMuon(const Track *mu, const Vertex* primaryVertex)
{
  assert(mu);

  TClonesArray &rMuonArr = *fMuonArr;
  assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
  const Int_t index = rMuonArr.GetEntries();  
  new(rMuonArr[index]) TMuon();
  TMuon *pMuon = (TMuon*)rMuonArr[index];                 
  
  pMuon->pt           = mu->Pt();
  pMuon->eta          = mu->Eta();
  pMuon->phi          = mu->Phi();
  pMuon->pfPt         = -999;
  pMuon->pfEta        = 0;
  pMuon->pfPhi        = 0;
  pMuon->pterr        = mu->PtErr();
  pMuon->trkIso03     = 0;
  pMuon->emIso03      = 0;
  pMuon->hadIso03     = 0;
  pMuon->hoIso03      = 0;
  pMuon->trkIso05     = 0;
  pMuon->emIso05      = 0;
  pMuon->hadIso05     = 0;
  pMuon->hoIso05      = 0;

  pMuon->d0           = mu->D0Corrected(*primaryVertex);
  pMuon->d0Err        = mu->D0Err();
  pMuon->dz           = mu->DzCorrected(*primaryVertex);
  pMuon->tkNchi2      = mu->RChi2();
  pMuon->muNchi2      = mu->RChi2();
  pMuon->q            = mu->Charge();
  pMuon->nValidHits   = 0;  
  pMuon->qualityBits  = 0;      
  pMuon->typeBits     = 0; 
  pMuon->nTkHits      = mu->NHits();
  pMuon->nPixHits     = mu->NPixelHits();
  pMuon->nSeg         = 0;
  pMuon->nMatch       = 0;
  pMuon->hltMatchBits = MatchHLT(mu->Pt(),mu->Eta(),mu->Phi());
  
  pMuon->staPt  = -1;  
  pMuon->staEta = -999;
  pMuon->staPhi = -999;
  

  pMuon->isMCReal = 0;
  pMuon->ip3d              = 0;
  pMuon->ip3dSig           = 0;

  pMuon->SegmentCompatibility     = 0;
  pMuon->CaloCompatilibity        = 0;

  pMuon->TrkKink           = 0;
  pMuon->GlobalKink        = 0;
  pMuon->HadEnergy         = 0;
  pMuon->HadS9Energy       = 0;
  pMuon->HoEnergy          = 0;
  pMuon->HoS9Energy        = 0;
  pMuon->EmEnergy          = 0;
  pMuon->EmS9Energy        = 0;

//   pMuon->Station0NSegments    = 0;
//   pMuon->Station0TrackDist    = 0;
//   pMuon->Station0TrackDistErr = 0;
//   pMuon->Station0dX           = 0;
//   pMuon->Station0dY           = 0;
//   pMuon->Station0PullX        = 0;
//   pMuon->Station0PullY        = 0;

//   pMuon->Station1NSegments    = 0;
//   pMuon->Station1TrackDist    = 0;
//   pMuon->Station1TrackDistErr = 0;
//   pMuon->Station1dX           = 0;
//   pMuon->Station1dY           = 0;
//   pMuon->Station1PullX        = 0;
//   pMuon->Station1PullY        = 0;

//   pMuon->Station2NSegments    = 0;
//   pMuon->Station2TrackDist    = 0;
//   pMuon->Station2TrackDistErr = 0;
//   pMuon->Station2dX           = 0;
//   pMuon->Station2dY           = 0;
//   pMuon->Station2PullX        = 0;
//   pMuon->Station2PullY        = 0;

//   pMuon->Station3NSegments    = 0;
//   pMuon->Station3TrackDist    = 0;
//   pMuon->Station3TrackDistErr = 0;
//   pMuon->Station3dX           = 0;
//   pMuon->Station3dY           = 0;
//   pMuon->Station3PullX        = 0;
//   pMuon->Station3PullY        = 0;

//   pMuon->Station4NSegments    = 0;
//   pMuon->Station4TrackDist    = 0;
//   pMuon->Station4TrackDistErr = 0;
//   pMuon->Station4dX           = 0;
//   pMuon->Station4dY           = 0;
//   pMuon->Station4PullX        = 0;
//   pMuon->Station4PullY        = 0;

//   pMuon->Station5NSegments    = 0;
//   pMuon->Station5TrackDist    = 0;
//   pMuon->Station5TrackDistErr = 0;
//   pMuon->Station5dX           = 0;
//   pMuon->Station5dY           = 0;
//   pMuon->Station5PullX        = 0;
//   pMuon->Station5PullY        = 0;

//   pMuon->Station6NSegments    = 0;
//   pMuon->Station6TrackDist    = 0;
//   pMuon->Station6TrackDistErr = 0;
//   pMuon->Station6dX           = 0;
//   pMuon->Station6dY           = 0;
//   pMuon->Station6PullX        = 0;
//   pMuon->Station6PullY        = 0;

//   pMuon->Station7NSegments    = 0;
//   pMuon->Station7TrackDist    = 0;
//   pMuon->Station7TrackDistErr = 0;
//   pMuon->Station7dX           = 0;
//   pMuon->Station7dY           = 0;
//   pMuon->Station7PullX        = 0;
//   pMuon->Station7PullY        = 0;


}




void HwwNtuplerMod::FillJet(const PFJet *jet, const Vertex* primaryVertex)
{
  TClonesArray &rPFJetArr = *fPFJetArr;
  assert(rPFJetArr.GetEntries() < rPFJetArr.GetSize());
  const Int_t index = rPFJetArr.GetEntries();  
  new(rPFJetArr[index]) higgsana::TJet();
  higgsana::TJet *pPFJet = (higgsana::TJet*)rPFJetArr[index]; 
   
  pPFJet->pt           = jet->Pt();
  pPFJet->eta          = jet->Eta();
  pPFJet->phi          = jet->Phi();
  pPFJet->mass         = jet->Mass();
  pPFJet->hltMatchBits = MatchHLT(jet->Eta(),jet->Phi());  
  pPFJet->TrackCountingHighEffBJetTagsDisc = jet->TrackCountingHighEffBJetTagsDisc();
  pPFJet->TrackCountingHighPurBJetTagsDisc  = jet->TrackCountingHighPurBJetTagsDisc();
  pPFJet->GhostTrackBJetTagsDisc = jet->GhostTrackBJetTagsDisc();
  pPFJet->SoftElectronByPtBJetTagsDisc = jet->SoftElectronByPtBJetTagsDisc();
  pPFJet->SoftElectronByIP3dBJetTagsDisc = jet->SoftElectronByIP3dBJetTagsDisc();
  pPFJet->SoftMuonByPtBJetTagsDisc = jet->SoftMuonByPtBJetTagsDisc();
  pPFJet->SoftMuonByIP3dBJetTagsDisc = jet->SoftMuonByIP3dBJetTagsDisc();
  pPFJet->SoftMuonBJetTagsDisc = jet->SoftMuonBJetTagsDisc();
  pPFJet->SimpleSecondaryVertexHighPurBJetTagsDisc = jet->SimpleSecondaryVertexHighPurBJetTagsDisc();
  pPFJet->SimpleSecondaryVertexHighEffBJetTagsDisc = jet->SimpleSecondaryVertexHighEffBJetTagsDisc();
  pPFJet->SimpleSecondaryVertexBJetTagsDisc = jet->SimpleSecondaryVertexBJetTagsDisc();
  pPFJet->CombinedSecondaryVertexBJetTagsDisc = jet->CombinedSecondaryVertexBJetTagsDisc();
  pPFJet->CombinedSecondaryVertexMVABJetTagsDisc = jet->CombinedSecondaryVertexMVABJetTagsDisc();
  pPFJet->PassBetaVertexAssociationCut = JetTools::PassBetaVertexAssociationCut(jet, fGoodPrimaryVertices->At(0), fGoodPrimaryVertices,0.2);
  pPFJet->NConstituents = jet->NConstituents();
  pPFJet->NeutralHadronFraction = jet->NeutralHadronEnergy() / jet->E();
  pPFJet->NeutralEMFraction = jet->NeutralEmEnergy() / jet->E();
  pPFJet->ChargedHadronFraction = jet->ChargedHadronEnergy() / jet->E();
  pPFJet->ChargedEMFraction = jet->ChargedEmEnergy() / jet->E();
  pPFJet->ChargedMultiplicity = jet->ChargedMultiplicity();
  pPFJet->JetProbabilityBJetTagsDisc = jet->JetProbabilityBJetTagsDisc();
  pPFJet->JetBProbabilityBJetTagsDisc = jet->JetBProbabilityBJetTagsDisc();

  pPFJet->JetArea = jet->JetArea();

  Double_t ChargedPtSqrSum = 0;
  Double_t DzSum = 0;
  Double_t DzPtSqrWeightedSum = 0;
  Int_t NChargedPFCandidates = 0;
  for (UInt_t i = 0; i < jet->NPFCands(); ++i) {
    //Charged PF Cand
    if (jet->PFCand(i)->HasTrackerTrk() || jet->PFCand(i)->HasGsfTrk()) {
      Double_t dz = 0;
      ChargedPtSqrSum += pow(jet->PFCand(i)->Pt(),2);
      NChargedPFCandidates++;
      if (jet->PFCand(i)->HasTrackerTrk()) {
        dz = jet->PFCand(i)->TrackerTrk()->DzCorrected(*primaryVertex);
      } else {
        dz = jet->PFCand(i)->GsfTrk()->DzCorrected(*primaryVertex);
      }

      DzSum += dz;
      DzPtSqrWeightedSum += dz * pow(jet->PFCand(i)->Pt(),2);
    }
  }
  
  pPFJet->DzAvg = DzSum / NChargedPFCandidates;
  pPFJet->DzPtSqrWeightedAvg = DzPtSqrWeightedSum / ChargedPtSqrSum;

  pPFJet->rawPt = jet->RawMom().Pt();

}
   

void HwwNtuplerMod::FillCaloJet(const CaloJet *jet)
{
  TClonesArray &rCaloJetArr = *fCaloJetArr;
  assert(rCaloJetArr.GetEntries() < rCaloJetArr.GetSize());
  const Int_t index = rCaloJetArr.GetEntries();  
  new(rCaloJetArr[index]) higgsana::TCaloJet();
  higgsana::TCaloJet *pCaloJet = (higgsana::TCaloJet*)rCaloJetArr[index]; 
   
  pCaloJet->pt           = jet->Pt();
  pCaloJet->eta          = jet->Eta();
  pCaloJet->phi          = jet->Phi();
  pCaloJet->mass         = jet->Mass();
  pCaloJet->hltMatchBits = MatchHLT(jet->Eta(),jet->Phi());  
  pCaloJet->TrackCountingHighEffBJetTagsDisc = jet->TrackCountingHighEffBJetTagsDisc();
  pCaloJet->TrackCountingHighPurBJetTagsDisc  = jet->TrackCountingHighPurBJetTagsDisc();
  pCaloJet->GhostTrackBJetTagsDisc = jet->GhostTrackBJetTagsDisc();
  pCaloJet->SoftElectronByPtBJetTagsDisc = jet->SoftElectronByPtBJetTagsDisc();
  pCaloJet->SoftElectronByIP3dBJetTagsDisc = jet->SoftElectronByIP3dBJetTagsDisc();
  pCaloJet->SoftMuonByPtBJetTagsDisc = jet->SoftMuonByPtBJetTagsDisc();
  pCaloJet->SoftMuonByIP3dBJetTagsDisc = jet->SoftMuonByIP3dBJetTagsDisc();
  pCaloJet->SoftMuonBJetTagsDisc = jet->SoftMuonBJetTagsDisc();
  pCaloJet->SimpleSecondaryVertexHighPurBJetTagsDisc = jet->SimpleSecondaryVertexHighPurBJetTagsDisc();
  pCaloJet->SimpleSecondaryVertexHighEffBJetTagsDisc = jet->SimpleSecondaryVertexHighEffBJetTagsDisc();
  pCaloJet->SimpleSecondaryVertexBJetTagsDisc = jet->SimpleSecondaryVertexBJetTagsDisc();
  pCaloJet->CombinedSecondaryVertexBJetTagsDisc = jet->CombinedSecondaryVertexBJetTagsDisc();
  pCaloJet->CombinedSecondaryVertexMVABJetTagsDisc = jet->CombinedSecondaryVertexMVABJetTagsDisc();
  pCaloJet->JetProbabilityBJetTagsDisc = jet->JetProbabilityBJetTagsDisc();
  pCaloJet->JetBProbabilityBJetTagsDisc = jet->JetBProbabilityBJetTagsDisc();

  pCaloJet->NConstituents = jet->NConstituents();

  pCaloJet->EmEnergyInEB = jet->EmEnergyInEB();
  pCaloJet->EmEnergyInEE = jet->EmEnergyInEE();
  pCaloJet->EmEnergyInHF = jet->EmEnergyInHF();
  pCaloJet->HadEnergyInHO = jet->HadEnergyInHO();
  pCaloJet->HadEnergyInHB = jet->HadEnergyInHB();
  pCaloJet->HadEnergyInHF = jet->HadEnergyInHF();
  pCaloJet->HadEnergyInHE = jet->HadEnergyInHE();
  pCaloJet->EnergyFractionH = jet->EnergyFractionH();
  pCaloJet->EnergyFractionEm = jet->EnergyFractionEm();

  pCaloJet->JetArea = jet->JetArea();
  
  pCaloJet->rawPt = jet->RawMom().Pt();

}
   


//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::FillPhoton(const Photon *pho)
{
  TClonesArray &rPhotonArr = *fPhotonArr;
  assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
  const Int_t index = rPhotonArr.GetEntries();  
  new(rPhotonArr[index]) TPhoton();
  TPhoton *pPhoton = (TPhoton*)rPhotonArr[index];
  
  pPhoton->et		  = pho->Et(); 
  pPhoton->eta  	  = pho->Eta();
  pPhoton->phi  	  = pho->Phi();
  pPhoton->scEt		  = pho->SCluster()->Et(); 
  pPhoton->scEta  	  = pho->SCluster()->Eta();
  pPhoton->scPhi  	  = pho->SCluster()->Phi();
  pPhoton->trkIso03Hollow = pho->HollowConeTrkIsoDr03(); 
  pPhoton->trkIso03Solid  = pho->SolidConeTrkIsoDr03();
  pPhoton->emIso03        = pho->EcalRecHitIsoDr03();
  pPhoton->hadIso03	  = pho->HcalTowerSumEtDr03(); 
  pPhoton->HoverE	  = pho->HadOverEm();
  pPhoton->R9		  = pho->R9();
  pPhoton->sigiEtaiEta    = pho->CoviEtaiEta();
  pPhoton->hltMatchBits   = MatchHLT(pho->SCluster()->Eta(),pho->SCluster()->Phi());
  pPhoton->scID           = pho->SCluster()->GetUniqueID();
  pPhoton->hasPixelSeed   = pho->HasPixelSeed();
  pPhoton->HoESingleTower = pho->HadOverEmTow();

  const BaseVertex *bsp = NULL;
  if(fBeamSpot->GetEntries() > 0) {
    bsp = fBeamSpot->At(0);
  }
  pPhoton->passConversionSafeEleVeto = PhotonTools::PassElectronVetoConvRecovery(pho,fElectrons,fConversions, bsp );

  //compute regression energy
  pPhoton->energyRegressionHgg2011 = 0;
  pPhoton->energyRegressionHgg2012 = 0;

  Photon *outph2011 = new Photon(*pho);
  Photon *outph2012 = new Photon(*pho);
  if (!egcorHgg2011->IsInitialized()) {cout << "Error: EGCor2011 not initialized.\n"; return;}
  if (!egcorHgg2012->IsInitialized()) {cout << "Error: EGCor2012 not initialized.\n"; return;}
  egcorHgg2012->CorrectEnergyWithError(outph2012,fGoodPrimaryVertices,fPileupEnergyDensity->At(0)->RhoKt6PFJets(),3, fUseGen);
  egcorHgg2011->CorrectEnergyWithError(outph2011,fGoodPrimaryVertices,fPileupEnergyDensity->At(0)->RhoKt6PFJets(),2, fUseGen);
  
  pPhoton->energyRegressionHgg2011 = outph2011->E();  
  pPhoton->energyRegressionHgg2012 = outph2012->E();
  //cout << "photon regression: " << pho->E() << " : " << outph2011->E() << " , " << outph2012->E() << endl;

  delete outph2011;
  delete outph2012;

}
 

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::FillPFCandidate(const PFCandidate *pf, 
                     UInt_t matchedObjectType, UInt_t matchedObjectIndex,
                     const Vertex* primaryVertex ) {

  TClonesArray &rPFCandidateArr = *fPFCandidateArr;

  if (!(rPFCandidateArr.GetEntries() < rPFCandidateArr.GetSize()) || rPFCandidateArr.GetEntries() >= 2000) {
    cout << "PFCandidate Array too big: " << rPFCandidateArr.GetEntries() << " " << rPFCandidateArr.GetSize() << endl;
  } 

  assert(rPFCandidateArr.GetEntries() < rPFCandidateArr.GetSize());
  const Int_t index = rPFCandidateArr.GetEntries();  
  new(rPFCandidateArr[index]) TPFCandidate();
  TPFCandidate *pPFCandidate = (TPFCandidate*)rPFCandidateArr[index]; 
  
  pPFCandidate->pt	  = pf->Pt(); 
  pPFCandidate->eta  	  = pf->Eta();
  pPFCandidate->phi  	  = pf->Phi();
  pPFCandidate->e  	  = pf->E();  

  if (pf->BestTrk()) {
    pPFCandidate->dz       = pf->BestTrk()->DzCorrected(*primaryVertex);
    pPFCandidate->q        = pf->BestTrk()->Charge();
  } else {
    pPFCandidate->dz       = 0.0;
    pPFCandidate->q        = 0;
  }


  //***********************************************************************
  //Figure out if PFCandidate is PFNoPU
  //***********************************************************************
  Bool_t isPFNoPU = kFALSE;
  if(pf->PFType() == PFCandidate::eHadron) {
    if(pf->HasTrackerTrk() && 
       fPrimVerts->At(0)->HasTrack(pf->TrackerTrk()) &&
       fPrimVerts->At(0)->TrackWeight(pf->TrackerTrk()) > 0) {
      isPFNoPU = kTRUE;
    } else { 

      Bool_t vertexFound = kFALSE;
      const Vertex *closestVtx = 0;
      Double_t dzmin = 10000;

      // loop over vertices
      for(UInt_t j = 0; j < fPrimVerts->GetEntries(); j++) {
	const Vertex *vtx = fPrimVerts->At(j);
	assert(vtx);
          
	if(pf->HasTrackerTrk() && 
	   vtx->HasTrack(pf->TrackerTrk()) &&
	   vtx->TrackWeight(pf->TrackerTrk()) > 0) {
	  vertexFound = kTRUE;
	  closestVtx = vtx;
	  break;
	}
	Double_t dz = fabs(pf->SourceVertex().Z() - vtx->Z());
	if(dz < dzmin) {
	  closestVtx = vtx;
	  dzmin = dz;
	}
      }

      Bool_t fCheckClosestZVertex = kTRUE; //we use option 1
      if(fCheckClosestZVertex) {
	// Fallback: if track is not associated with any vertex,
	// associate it with the vertex closest in z
	if(vertexFound || closestVtx != primaryVertex) {
	  isPFNoPU = kFALSE;
	} else {
	  isPFNoPU = kTRUE;
	}
      } else {
	if(vertexFound && closestVtx != primaryVertex) {
	  isPFNoPU = kFALSE;
	} else {
	  isPFNoPU = kTRUE;
	}
      }
    } //hadron & trk stuff
  } else { // hadron
    //
    isPFNoPU = kTRUE;
  }  
  pPFCandidate->IsPFNoPU = isPFNoPU;

  pPFCandidate->pfType    = pf->PFType();
  pPFCandidate->matchedObjectType  = matchedObjectType;
  pPFCandidate->matchedObjectIndex  = matchedObjectIndex;

//   ULong_t tmpFlags = 0;
//   if (pf->Flag(PFCandidate::eNormal)) tmpFlags += 1;
//   if (pf->Flag(PFCandidate::eEMPhiSModules)) tmpFlags += 2;
//   if (pf->Flag(PFCandidate::eEMEta0)) tmpFlags += 4;
//   if (pf->Flag(PFCandidate::eEMEtaModules)) tmpFlags += 8;
//   if (pf->Flag(PFCandidate::eEMBarrelEndcap)) tmpFlags += 16;
//   if (pf->Flag(PFCandidate::eEMPreshowerEdge)) tmpFlags += 32;
//   if (pf->Flag(PFCandidate::eEMPreshower)) tmpFlags += 64;
//   if (pf->Flag(PFCandidate::eEMEndCapEdge)) tmpFlags += 128;
//   if (pf->Flag(PFCandidate::eHEta0)) tmpFlags += 256;
//   if (pf->Flag(PFCandidate::eHBarrelEndcap)) tmpFlags += 512;
//   if (pf->Flag(PFCandidate::eHEndcapVFCal)) tmpFlags += 1024;
//   if (pf->Flag(PFCandidate::eHVFCalEdge)) tmpFlags += 2048;
//   if (pf->Flag(PFCandidate::eToDispVtx)) tmpFlags += 4096;
//   if (pf->Flag(PFCandidate::eFromDispVtx)) tmpFlags += 8192;
//   if (pf->Flag(PFCandidate::eFromV0)) tmpFlags += 16384;
//   if (pf->Flag(PFCandidate::eFromGammaConv)) tmpFlags += 32768;
//   if (pf->Flag(PFCandidate::eToConversion)) tmpFlags += 65536;
//   pPFCandidate->pfFlags  = tmpFlags;



}

//--------------------------------------------------------------------------------------------------
ULong_t HwwNtuplerMod::MatchL1(const Double_t eta, const Double_t phi)
{

  ULong_t bits = 0;
  
  const Double_t L1MatchDR = 0.5;
  
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);

    for(ULong_t iseed=0; iseed<fL1SeedModuleNamesv.size(); iseed++) {

      for(ULong_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {

        const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
        if(!trigname) continue;
   
        const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
        if(!list) continue;
        TIter iter(list->MakeIterator());
        const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
    
        while(to) {         
        
          
        
          if(to->IsL1()) {
            Bool_t match = true;
            
            //match L1 seed module name
            if (!(string(to->ModuleName()) == string(fL1SeedModuleNamesv[iseed].Data()) )) { match = false; }

//             if (match) cout << to->Pt() << " " << to->Eta() << " " << to->Phi() << " : " << to->IsHLT() << " " << to->IsL1() << " " << to->TrigName() << " " << to->FilterName() << " " << to->ModuleName() << " " << to->TagName() << " " << endl;
        


            // eta-phi matching
            if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > L1MatchDR) { match = false; }

            // set appropriate bits
            if(match) { 
              bits |= fL1SeedModuleIdsv[iseed]; 
            }
          }
          to = dynamic_cast<const TriggerObject*>(iter.Next());
        }
      }
    }
  }
  
  return bits;
}


ULong_t HwwNtuplerMod::MatchHLT(const Double_t eta, const Double_t phi, Double_t runNum, Double_t evtNum)
{



//   if (fPrintDebug) {
//     cout << "Match To HLT Object : " << eta << " " << phi << endl;
//   }

  ULong_t bits = 0;
  
  const Double_t hltMatchR = 0.1;
  
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(ULong_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      
      bool doDebug = kFALSE;

      if (doDebug) cout << "DEBUG: Match To HLT Object : " << eta << " " << phi << endl;

      const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if(!list) continue;
      TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
    
      while(to) {         
        
        if (doDebug || fPrintDebug) {         
          cout << to->Pt() << " " << to->Eta() << " " << to->Phi() << " : " << to->IsHLT() << " " << to->IsL1() << " " << to->TrigName() << " " << to->FilterName() << " " << to->ModuleName() << " " << to->TagName() << " : " 
               << MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) << " --> " << hltMatchR << " : " 
               << fFirstTriggerObjectIdsv[itrig] << " " << fFirstTriggerObjectModuleNamesv[itrig].Data() << " "  
               << (fFirstTriggerObjectIdsv[itrig] != 0) << " " 
               << (string(to->ModuleName()) == string(fFirstTriggerObjectModuleNamesv[itrig].Data()) || string(fFirstTriggerObjectModuleNamesv[itrig].Data()) == "") << " " 
               << endl;          
        }
        
        if(to->IsHLT()) {
          
          //Match First Trigger Object
          if (fFirstTriggerObjectIdsv[itrig] != 0 && 
              (string(to->ModuleName()) == string(fFirstTriggerObjectModuleNamesv[itrig].Data()) || string(fFirstTriggerObjectModuleNamesv[itrig].Data()) == "")
            ) {
            Bool_t match = true;
            
            // eta-phi matching
            if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) { match = false; }

            // set appropriate bits
            if(match) { 
              bits |= fFirstTriggerObjectIdsv[itrig]; 
              if (doDebug ||fPrintDebug) {
                cout << "Matched First : " << fFirstTriggerObjectModuleNamesv[itrig].Data() << " : " << endl; 
              }
            }
          }
          
          //Match Second Trigger Object
          if (fSecondTriggerObjectIdsv[itrig] != 0 && 
              (string(to->ModuleName()) == string(fSecondTriggerObjectModuleNamesv[itrig].Data()) || string(fSecondTriggerObjectModuleNamesv[itrig].Data()) == "")
            ) {
            Bool_t match = true;
            
            // eta-phi matching
            if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) { match = false; }
            
            // set appropriate bits
            if(match) { 
              bits |= fSecondTriggerObjectIdsv[itrig];               
              if (doDebug ||fPrintDebug) {
                cout << "Matched Second : " << fFirstTriggerObjectModuleNamesv[itrig].Data() << " : " << endl; 
              }
            }
          }
        }
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }
    }
  }
  
  return bits;
}

ULong_t HwwNtuplerMod::MatchHLTMuon(const Double_t pt, const Double_t eta, const Double_t phi, Double_t runNum, Double_t evtNum)
{
  ULong_t bits = 0;
  
//   if (fPrintDebug) {
//     cout << "Match HLT Mu : " << pt << " " << eta << " " << phi << endl;
//   }


  const Double_t hltMatchR = 0.2;
//   const Double_t hltMatchPtFrac = 10;
  
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(ULong_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      
      const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if(!list) continue;
      TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
    
      while(to) {  

//         if (fPrintDebug) {
    
//           cout << to->Pt() << " " << to->Eta() << " " << to->Phi() << " : " << to->IsHLT() << " " << to->IsL1() << " " << to->TrigName() << " " << to->FilterName() << " " << to->ModuleName() << " " << to->TagName() << " " << endl;
          
//         }
        
        if(to->IsHLT()) {
          
          //Match First Trigger Object
          if (fFirstTriggerObjectIdsv[itrig] != 0 && 
              (string(to->ModuleName()) == string(fFirstTriggerObjectModuleNamesv[itrig].Data()) || string(fFirstTriggerObjectModuleNamesv[itrig].Data()) == "")
            ) {
            Bool_t match = true;
            
            // eta-phi matching
            if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) { match = false; }
            
//             // pT matching
//             if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt())) { match = false; }
            
            // set appropriate bits
            if(match) { 
              bits |= fFirstTriggerObjectIdsv[itrig]; 
//               if (fPrintDebug) {
//                 cout << "Matched First : " << fFirstTriggerObjectModuleNamesv[itrig].Data() << " : " << endl; 
//               }              
            }
          }
          
          //Match Second Trigger Object
          if (fSecondTriggerObjectIdsv[itrig] != 0 && 
              (string(to->ModuleName()) == string(fSecondTriggerObjectModuleNamesv[itrig].Data()) || string(fSecondTriggerObjectModuleNamesv[itrig].Data()) == "")
            ) {
            Bool_t match = true;
            
            // eta-phi matching
            if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) { match = false; }
            
//             // pT matching
//             if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt())) { match = false; }
            
            // set appropriate bits
            if(match) { 
              bits |= fSecondTriggerObjectIdsv[itrig]; 
//               if (fPrintDebug) {
//                 cout << "Matched Second : " << fFirstTriggerObjectModuleNamesv[itrig].Data() << " : " << endl; 
//               }
            }
          }

        }
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }    
    }
  }
  
  return bits;
}


