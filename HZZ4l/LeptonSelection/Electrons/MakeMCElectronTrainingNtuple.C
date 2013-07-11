//================================================================================================
//
// HWW selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TGraphAsymmErrors.h>     
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

// define structures to read in ntuple
#include "HiggsAna/DataTree/interface/HiggsAnaDefs.hh"
#include "HiggsAna/DataTree/interface/TEventInfo.hh"
#include "HiggsAna/DataTree/interface/TElectron.hh"
#include "HiggsAna/DataTree/interface/TPhoton.hh"
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/TJet.hh"
#include "HiggsAna/DataTree/interface/TPFCandidate.hh"
#include "HiggsAna/DataTree/interface/TGenParticle.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "HiggsAna/CommonData/interface/ElectronTree.h"
#include "HiggsAna/Utils/LeptonTools.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"


#endif

//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
void MakeNtuple(const string inputFilename,  const string outputFilename, const string PUReweightFile, Bool_t SelectRealElectrons, Int_t Option = 0);

//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = TFile::Open(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwNtuplerMod");
    if (!dir) {
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}


//=== MAIN MACRO =================================================================================================
void MakeMCElectronTrainingNtuple(Int_t Sample = 0) {
 
  //HZZ 
  if (Sample == 1) {
    MakeNtuple("/data/smurf/sixie/hist/HZZ4lNtuples_old/mc/AllNtuple_HZZ4lNtuple_s12-h125zz4l-gf-v9_noskim_0004.root","ElectronSelectionTraining.s11HZZ125.Training.root", "/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root", kTRUE);
  }

  //HWW
  if (Sample == 2) {
    
  }

} 


void MakeMCElectronTrainingNtuple(const string inputFilename, const string outputFilename, const string PUReweightFile, Bool_t SelectRealElectrons, Int_t Option = 0) {
  MakeNtuple(inputFilename,outputFilename,PUReweightFile,SelectRealElectrons, Option);
}



void MakeNtuple(const string inputFilename, const string outputFilename, const string PUReweightFile, Bool_t SelectRealElectrons, Int_t Option )
{  
  gBenchmark->Start("WWTemplate");

  TFile *fPUFile = 0;
  if (PUReweightFile != "") {
    fPUFile = TFile::Open(PUReweightFile.c_str());
  }
  TH1D *fhDPU = 0;
  if (fPUFile) fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  if(fhDPU) fhDPU->SetDirectory(0);
  if (fPUFile) delete fPUFile;

  //*****************************************************************************************
  //Define Data Era
  //*****************************************************************************************
  UInt_t DataEra = kDataEra_NONE;
  if (PUReweightFile == "/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root") {
    DataEra = kDataEra_2012_MC;
  } else if (PUReweightFile == "/data/smurf/sixie/Pileup/weights/PileupReweighting.Fall11_To_Full2011.root") {
    DataEra = kDataEra_2012_MC;
  } else if (PUReweightFile == "/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/puWeights_Summer12_5000ipb_71mb.root") {
    DataEra = kDataEra_2012_MC;
  }

  if (Option == 1) {
    DataEra = kDataEra_NONE;    
  }

  cout << "using DataEra = " << DataEra << endl;


  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  TFile *outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  ElectronTree *eleTree = new ElectronTree;
  eleTree->CreateTree();
  eleTree->tree_->SetAutoFlush(0);

  UInt_t NElectronsFilled = 0;
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  higgsana::TEventInfo *info    = new higgsana::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("higgsana::TElectron");
  TClonesArray *pfcandidateArr = new TClonesArray("higgsana::TPFCandidate");
  TClonesArray *genparticleArr = new TClonesArray("higgsana::TGenParticle");

  Int_t NEvents = 0;



  vector<string> inputfiles;
  if (inputFilename == "LIST") {
    inputfiles.push_back("/data/blue/sixie/ntuples/HWW/mc/HwwAnalysis_s11-h115ww2l-gf-v11-pu_noskim_normalized.root");
  } else if (inputFilename == "ZPlusFake") {
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zmmm1020-powheg-v14b-pu_noskim.MCFakeLeptonSkimmed.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim.MCFakeLeptonSkimmed.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zmmm20-powheg-v14b-pu_noskim.MCFakeLeptonSkimmed.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim.MCFakeLeptonSkimmed.root");
  } else if (inputFilename == "TTBAR") {
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-tt2l-powheg-v14b-pu_noskim_0000.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-tt2l-powheg-v14b-pu_noskim_0002.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-tt2l-powheg-v14b-pu_noskim_0003.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-tt2l-powheg-v14b-pu_noskim_0004.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-tt2l-powheg-v14b-pu_noskim_0005.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-tt2l-powheg-v14b-pu_noskim_0006.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-tt2l-powheg-v14b-pu_noskim_0007.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-tt2l-powheg-v14b-pu_noskim_0008.root");
  } else if (inputFilename == "Fall11ZeeM1020") {
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0000.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0001.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0002.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0003.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0004.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0005.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0006.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0007.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0008.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0009.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0010.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0011.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0012.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0013.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0014.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0015.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0016.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0017.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0018.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0019.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0020.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0021.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0022.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0023.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0024.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0025.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0026.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0027.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0028.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0029.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0030.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0031.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0032.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0033.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem1020-powheg-v14b-pu_noskim_0034.root");
  } else if (inputFilename == "Fall11ZeeM20") {
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0000.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0001.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0002.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0003.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0004.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0005.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0006.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0007.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0008.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0009.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0010.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0011.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0012.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0013.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0014.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0015.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0016.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0017.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0018.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0019.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0020.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0021.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0022.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0023.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0024.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0025.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0026.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0027.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0028.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0029.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0030.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0031.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0032.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0033.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0034.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0035.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0036.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0037.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0038.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0039.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0040.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0041.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0042.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0043.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0044.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0045.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0046.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0047.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0048.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0049.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0050.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0051.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0052.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0053.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0054.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0055.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0056.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0057.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0058.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0059.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0060.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0061.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0062.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0063.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0064.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0065.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0066.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0067.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0068.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0069.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0070.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0071.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0072.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0073.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0074.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0075.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0076.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0077.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0078.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0079.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0080.root");
//     inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-zeem20-powheg-v14b-pu_noskim_0081.root");
  } else if (inputFilename == "Summer12ZJets_51X") {
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/mc/AllNtuple_HZZ4lNtuple_s12-zjets-m50-v15_noskim_0000.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/mc/AllNtuple_HZZ4lNtuple_s12-zjets-m50-v15_noskim_0001.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/mc/AllNtuple_HZZ4lNtuple_s12-zjets-m50-v15_noskim_0002.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/mc/AllNtuple_HZZ4lNtuple_s12-zjets-m50-v15_noskim_0003.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/mc/AllNtuple_HZZ4lNtuple_s12-zjets-m50-v15_noskim_0004.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/mc/AllNtuple_HZZ4lNtuple_s12-zjets-m50-v15_noskim_0005.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/mc/AllNtuple_HZZ4lNtuple_s12-zjets-m50-v15_noskim_0006.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/mc/AllNtuple_HZZ4lNtuple_s12-zjets-m50-v15_noskim_0007.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/mc/AllNtuple_HZZ4lNtuple_s12-zjets-m50-v15_noskim_0008.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/mc/AllNtuple_HZZ4lNtuple_s12-zjets-m50-v15_noskim_0009.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/mc/AllNtuple_HZZ4lNtuple_s12-zjets-m50-v15_noskim_0010.root");
  } else {
    inputfiles.push_back(inputFilename);
  }

  for (UInt_t f = 0; f < inputfiles.size(); ++f) {

    //********************************************************
    // Get Tree
    //********************************************************
    eventTree = getTreeFromFile(inputfiles[f].c_str(),"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *pfcandidateBr;
    TBranch *genparticleBr;


    //*****************************************************************************************
    //Loop over Data Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr); pfcandidateBr = eventTree->GetBranch("PFCandidate");
    eventTree->SetBranchAddress("GenParticle", &genparticleArr); genparticleBr = eventTree->GetBranch("GenParticle");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 10000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;

      NEvents++;

      //********************************************************
      // Load the branches
      //********************************************************
      electronArr->Clear(); 
      pfcandidateArr->Clear(); 
      genparticleArr->Clear(); 
      electronBr->GetEntry(ientry);
      pfcandidateBr->GetEntry(ientry);
      genparticleBr->GetEntry(ientry);

      //********************************************************
      double npuWeight = 1;
      if (fhDPU) {
        double mynpu = TMath::Min((double)info->nPUEvents,34.999);
        Int_t npuxbin = fhDPU->GetXaxis()->FindBin(mynpu);
        npuWeight = fhDPU->GetBinContent(npuxbin);
      }
      //********************************************************

      //********************************************************
      // Pileup Energy Density
      //********************************************************
      Double_t rhoEleIso = 0;
      UInt_t EleEAEra = 0;

      if (DataEra == kDataEra_2011_MC) {     
        if (!(isnan(info->RhoKt6PFJetsForIso25) || 
              isinf(info->RhoKt6PFJetsForIso25))) {
          rhoEleIso = info->RhoKt6PFJetsForIso25;
        }
        EleEAEra = kDataEra_2011_Data;
      } else if (DataEra == kDataEra_2012_MC) {
        if (!(isnan(info->RhoKt6PFJets) || 
              isinf(info->RhoKt6PFJets))) {
          rhoEleIso = info->RhoKt6PFJets;
        }
        EleEAEra = kDataEra_2012_Data;
      }


      //********************************************************
      // Loop Over Electrons
      //********************************************************
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);

        //********************************************************
        // Select MC Truth Electrons
        //********************************************************
              
        //use only real electrons
        if (SelectRealElectrons) {
          if (!((UInt_t(abs(max(0,ele->isMCReal))) & 1) == 1)) continue;
        } else {
          //select fake electrons only
          if (!(ele->isMCReal == 0)) continue;          
        }

        //protect against pathologies
        if (TMath::IsNaN(ele->sigiPhiiPhi)) continue;
        
        //Fill These Electrons
        FillElectronTree( eleTree, ele, pfcandidateArr, rhoEleIso, EleEAEra, 
                          info->nPV0, info->runNum, info->lumiSec, info->evtNum, npuWeight);

      } //loop over electrons

    }

    cout << "Total Electrons: " << NElectronsFilled << endl;

  } //end loop over files

  delete info;
  delete electronArr;
  delete pfcandidateArr;
  delete genparticleArr;

  cout << "Total Electrons: " << NElectronsFilled << endl;

  outputFile->Write();
  outputFile->Close();
  
  delete outputFile;
  if (eleTree) delete eleTree;

  gBenchmark->Show("WWTemplate");       
} 


