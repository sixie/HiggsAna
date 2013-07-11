//****************************************************************************************************
//
// foreach f ( `ls /data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu*Tight* | sed 's/\/data\/smurf\/sixie\/hist\/HZZ4lNtuples\/data\/unprocessed\///' ` )
// root -l -b -q HiggsAna/HZZ4l/LeptonSelection/Electrons/MakeFakeElectronTrainingNtupleFromWPlusJetSample.C+\(\"/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/$f\",\"/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet.${f}.root\"\)
// end
//
//****************************************************************************************************
//
//root -l -b -q EWKAna/Hww/FakeRate/ComputeElectronFakeRate_Data.C+\(\)
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
#include "HiggsAna/DataTree/interface/Types.h"
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "HiggsAna/Utils/LeptonTools.hh"
#include "HiggsAna/Utils/CommonTools.hh"
#include "HiggsAna/CommonData/interface/ElectronTree.h"
#include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"

#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
// #include "MitCommon/DataFormats/interface/Types.h"
// #include "MitCommon/MathTools/interface/MathUtils.h"


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
void MakeNtuple(const string inputFilename,  const string outputFilename);

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
void MakeFakeElectronTrainingNtupleFromWPlusJetSample(Int_t Option = 0) {

  // MakeNtuple("LIST","ElectronSelectionTraining.Fake_WPlusJet.Training.root");
  // MakeNtuple("LIST","ElectronSelectionTraining.Fake_WPlusJet.Testing.root");
//   MakeNtuple("WToMuNuPlusFake","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet.root");
  //MakeNtuple("WToENuPlusFake","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToENuPlusJet.root");

  if (Option == 3) {
    MakeNtuple("WToMuNuPlusFake2012A","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet_2012A.root");
  }
  if (Option == 4) {
    MakeNtuple("WToMuNuPlusFake2012B","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToMuNuPlusJet_2012B.root");
  }

  if (Option == 5) {
    MakeNtuple("WToENuPlusFake2012A","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToENuPlusJet_2012A.root");
  }
  if (Option == 6) {
    MakeNtuple("WToENuPlusFake2012B","/data/blue/sixie/HZZ4l/LeptonSelection/Electrons/ElectronSelectionTraining.Fake_WToENuPlusJet_2012B.root");
  }

  if (Option == -1) {
    MakeNtuple("Test","ElectronSelectionTraining.Fake_WToENuPlusJet_Test.root");
  }

}

//=== MAIN MACRO =================================================================================================
void MakeFakeElectronTrainingNtupleFromWPlusJetSample(string inputFilename, string outputFilename) {

  MakeNtuple(inputFilename, outputFilename);
  
}




void MakeNtuple(const string inputFilename, const string outputFilename)
{  
  gBenchmark->Start("WWTemplate");

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
  TClonesArray *muonArr = new TClonesArray("higgsana::TMuon");
  TClonesArray *jetArr = new TClonesArray("higgsana::TJet");
  TClonesArray *photonArr = new TClonesArray("higgsana::TPhoton");
  TClonesArray *pfcandidateArr = new TClonesArray("higgsana::TPFCandidate");

  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("/data/smurf/data/Winter11_4700ipb/auxiliar/hww.Full2011.json"); 
  rlrm.AddJSONFile("/data/blue/sixie/HZZ4l/auxiliar/2012/Cert_190456-196531_8TeV_PromptReco_Collisions12_JSON.txt");

  Int_t NEvents = 0;

  UInt_t DataEra = kDataEra_2012_MC;

  enum { eMC, eMuEl, eDiMu, eMu, eDiEl, eEl };  // dataset type

  vector<string> inputfiles;
  if (inputFilename == "LIST") {
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-smu-m10-v1.TightPlusReco.root"); 
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-smu-pr-v4.TightPlusReco.root");	
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-smu-a05-v1.TightPlusReco.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-smu-o03-v1.TightPlusReco.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11b-smu-pr-v1.TightPlusReco.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-sel-a05-v1.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-sel-m10-v1.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-sel-o03-v1.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-sel-pr-v4.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11b-sel-pr-v1.TightPlusRecoSkimmed.root");
  }
  else if (inputFilename == "Test") {
//     inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0010.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0010.root.TightPlusRecoSkimmed.root");
    DataEra = kDataEra_2012_MC;
  }
  else if (inputFilename == "WToMuNuPlusFake") {
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-smu-a05-v1_noskim.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-smu-m10-v1_noskim.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-smu-o03-v1_noskim.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-smu-pr-v4_noskim.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11b-smu-pr-v1_noskim.TightPlusRecoSkimmed.root");
  }
  else if (inputFilename == "WToENuPlusFake") {
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-sel-a05-v1.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-sel-m10-v1.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-sel-o03-v1.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11a-sel-pr-v4.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/blue/sixie/hist/HZZ4lNtuples_NEW/data/AllNtuple_HZZ4lNtuple_r11b-sel-pr-v1.TightPlusRecoSkimmed.root");
  }
  else if (inputFilename == "WToMuNuPlusFake2012A") {
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0000.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0001.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0002.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0003.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0004.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0005.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0006.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0007.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0008.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0009.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0010.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0011.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0012.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0013.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0014.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0015.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0016.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0017.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0018.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0019.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0020.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0021.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0022.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0023.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0024.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0025.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0026.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0027.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0028.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0029.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0030.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0031.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0032.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0033.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0034.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0035.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0036.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0037.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0038.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0039.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0040.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0041.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0042.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0043.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0044.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0045.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0046.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0047.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0048.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0049.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0050.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0051.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0052.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0053.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0054.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0055.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-smu-pr-v1_noskim_0056.root.TightPlusRecoSkimmed.root");
    DataEra = kDataEra_2012_MC;
  }
  else if (inputFilename == "WToMuNuPlusFake2012B") {
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0000.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0001.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0002.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0003.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0004.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0005.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0006.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0007.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0008.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0009.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0010.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0011.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0012.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0013.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0014.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0015.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0016.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0017.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0018.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0019.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0020.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0021.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0022.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0023.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0024.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0025.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0026.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0027.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0028.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0029.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0030.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0031.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0032.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0033.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0034.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0035.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0036.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0037.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0038.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0039.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0040.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0041.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0042.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0043.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0044.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0045.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0046.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0047.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0048.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0049.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0050.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0051.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0052.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0053.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0054.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0055.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0056.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0057.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0058.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0059.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0060.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0061.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0062.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0063.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0064.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0065.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0066.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0067.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0068.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0069.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0070.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0071.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0072.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0073.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0074.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0075.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0076.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0077.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0078.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0079.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0080.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0081.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0082.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0083.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0084.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0085.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0086.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0087.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0088.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0089.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0090.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0091.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0092.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0093.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0094.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0095.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0096.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0097.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0098.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0099.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0100.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0101.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0102.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0103.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0104.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0105.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0106.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0107.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0108.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0109.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0110.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0111.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0112.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0113.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0114.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0115.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0116.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0117.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0118.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0119.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0120.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0121.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0122.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0123.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0124.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0125.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0126.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0127.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0128.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0129.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0130.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0131.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0132.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0133.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0134.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0135.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0136.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0137.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0138.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0139.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0140.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0141.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0142.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0143.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0144.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0145.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0146.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0147.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0148.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0149.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0150.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0151.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0152.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0153.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0154.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0155.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0156.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0157.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0158.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0159.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0160.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0161.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0162.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0163.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0164.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0165.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0166.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0167.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0168.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0169.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0170.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0171.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0172.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0173.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0174.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0175.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0176.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0177.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0178.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0179.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0180.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0181.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0182.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0183.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0184.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0185.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0186.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0187.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0188.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0189.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0190.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0191.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0192.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0193.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0194.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0195.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0196.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0197.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0198.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0199.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0200.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0201.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0202.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0203.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0204.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0205.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0206.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0207.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0208.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0209.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0210.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0211.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0212.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0213.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0214.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0215.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0216.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0217.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0218.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0219.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0220.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0221.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0222.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0223.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0224.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0225.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0226.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0227.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0228.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0229.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0230.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0231.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0232.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0233.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0234.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0235.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0236.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0237.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0238.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0239.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0240.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0241.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0242.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0243.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0244.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0245.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0246.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0247.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0248.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0249.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0250.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0251.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0252.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0253.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0254.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0255.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0256.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0257.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0258.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0259.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0260.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0261.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0262.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0263.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-smu-pr-v1_noskim_0264.root.TightPlusRecoSkimmed.root");
    DataEra = kDataEra_2012_MC;
  }
  else if (inputFilename == "WToENuPlusFake2012A") {
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0000.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0001.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0002.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0003.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0004.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0006.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0007.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0008.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0009.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0010.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0011.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0012.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0013.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0014.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0015.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0016.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0017.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0018.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0019.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0020.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0021.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0022.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0023.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0024.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0025.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0027.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0028.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0029.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0030.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0031.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0032.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0033.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0034.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0035.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0036.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0037.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0038.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0039.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0040.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0041.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0042.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0043.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0044.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0045.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0046.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0047.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0048.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0049.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0050.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0051.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0052.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0053.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0054.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0055.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0056.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0057.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0058.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0059.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0060.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0061.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0062.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0063.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0064.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0065.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0066.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0067.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0068.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0069.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0070.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0071.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12a-sel-pr-v1_noskim_0072.root.TightPlusRecoSkimmed.root");
    DataEra = kDataEra_2012_MC;
  }
  else if (inputFilename == "WToENuPlusFake2012B") {
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0000.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0001.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0002.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0003.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0004.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0005.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0006.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0007.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0008.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0009.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0010.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0011.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0012.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0013.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0014.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0015.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0016.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0017.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0018.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0019.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0020.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0021.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0022.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0023.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0024.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0025.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0026.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0027.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0028.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0029.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0030.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0031.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0032.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0033.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0034.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0035.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0036.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0037.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0038.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0039.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0040.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0041.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0042.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0043.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0044.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0045.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0046.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0047.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0048.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0049.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0050.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0051.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0052.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0053.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0054.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0055.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0056.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0057.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0058.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0059.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0060.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0061.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0062.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0063.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0064.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0065.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0066.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0067.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0068.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0069.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0070.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0071.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0072.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0073.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0074.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0075.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0076.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0077.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0078.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0079.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0080.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0081.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0082.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0083.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0084.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0085.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0086.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0087.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0088.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0089.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0090.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0091.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0092.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0093.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0094.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0095.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0096.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0097.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0098.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0099.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0100.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0101.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0102.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0103.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0104.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0106.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0107.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0108.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0109.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0110.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0111.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0112.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0113.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0114.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0115.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0116.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0117.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0118.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0119.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0120.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0121.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0122.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0123.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0124.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0125.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0126.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0127.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0128.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0129.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0130.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0131.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0132.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0133.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0134.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0135.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0136.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0137.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0138.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0139.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0140.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0141.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0142.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0143.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0144.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0145.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0146.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0147.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0148.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0149.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0150.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0151.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0152.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0153.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0154.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0155.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0156.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0157.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0158.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0159.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0160.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0161.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0162.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0163.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0164.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0165.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0166.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0167.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0168.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0169.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0170.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0171.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0172.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0173.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0174.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0175.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0176.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0177.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0178.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0179.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0180.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0181.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0182.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0183.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0184.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0185.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0186.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0187.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0188.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0189.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0190.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0191.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0192.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0193.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0194.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0195.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0196.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0197.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0198.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0199.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0200.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0201.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0202.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0203.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0204.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0205.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0206.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0208.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0209.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0210.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0211.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0212.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0213.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0214.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0215.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0216.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0217.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0218.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0219.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0220.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0221.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0222.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0223.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0224.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0225.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0226.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0227.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0228.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0229.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0230.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0231.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0232.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0233.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0234.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0235.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0236.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0237.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0238.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0239.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0240.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0241.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0242.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0243.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0244.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0245.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0246.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0247.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0248.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0249.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0250.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0251.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0252.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0253.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0254.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0255.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0256.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0257.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0258.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0259.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0260.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0261.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0262.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0263.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0264.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0265.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0266.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0267.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0268.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0269.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0270.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0271.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0272.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0273.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0274.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0275.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0276.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0277.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0278.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0279.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0280.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0281.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0282.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0283.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0284.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0285.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0286.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0287.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0288.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0289.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0290.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0291.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0292.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0293.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0294.root.TightPlusRecoSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/unprocessed/AllNtuple_HZZ4lNtuple_r12b-sel-pr-v1_noskim_0295.root.TightPlusRecoSkimmed.root");
    DataEra = kDataEra_2012_MC;
  }
  else {
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
    TBranch *jetBr;
    TBranch *photonBr;
    TBranch *pfcandidateBr;


    //*****************************************************************************************
    //Loop over Data Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);      infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr     = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Photon", &photonArr);     photonBr   = eventTree->GetBranch("Photon");
    eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr      = eventTree->GetBranch("PFJet");
    eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr); pfcandidateBr = eventTree->GetBranch("PFCandidate");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if(ientry % 100000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;

      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

      NEvents++;

      //********************************************************
      // Load the branches
      //********************************************************
      electronArr->Clear(); 
      muonArr->Clear(); 
      photonArr->Clear(); 
      jetArr->Clear(); 
      pfcandidateArr->Clear(); 
      electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);
      photonBr->GetEntry(ientry);
      jetBr->GetEntry(ientry);
      pfcandidateBr->GetEntry(ientry);

      //********************************************************
      // Pileup Energy Density
      //********************************************************
      Double_t rhoMuonIso = 0;
      Double_t rhoEleIso = 0;
      UInt_t MuonEAEra = 0;
      UInt_t EleEAEra = 0;

      if (DataEra == kDataEra_2011_MC) {     
        if (!(isnan(info->RhoKt6PFJetsForIso25) || 
              isinf(info->RhoKt6PFJetsForIso25))) {
          rhoMuonIso = info->RhoKt6PFJetsForIso25;
          rhoEleIso = info->RhoKt6PFJetsForIso25;
        }
        MuonEAEra = kDataEra_2011_Data;
        EleEAEra = kDataEra_2011_Data;
      } else if (DataEra == kDataEra_2012_MC) {
        if (!(isnan(info->RhoKt6PFJetsCentralNeutral) || 
              isinf(info->RhoKt6PFJetsCentralNeutral))) {
          rhoMuonIso = info->RhoKt6PFJetsCentralNeutral;
        }
        if (!(isnan(info->RhoKt6PFJets) || 
              isinf(info->RhoKt6PFJets))) {
          rhoEleIso = info->RhoKt6PFJets;
        }
        MuonEAEra = kDataEra_2012_Data;
        EleEAEra = kDataEra_2012_Data;
      }


      //********************************************************
      // TcMet
      //********************************************************
      TVector3 pfMet;        
      if(info->pfMEx!=0 || info->pfMEy!=0) {       
        pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
      }
      Double_t met = pfMet.Pt();


      //********************************************************
      //Some  Events Selection : suppress ttbar
      //********************************************************
      int njet=0;
      float maxtche=0;
      for(Int_t ijet=0; ijet<jetArr->GetEntries(); ijet++) {
        const higgsana::TJet *jet = (higgsana::TJet*)((*jetArr)[ijet]);
	if(jet->pt > 30) njet++;
	if(jet->pt > 15 && jet->TrackCountingHighEffBJetTagsDisc > maxtche) maxtche = jet->TrackCountingHighEffBJetTagsDisc;
      }
      // jet veto
      if(njet>1) continue;
      // btag veto
      if(maxtche >= 2.1) continue;

      //********************************************************
      // W-Selection
      //********************************************************
      vector<const higgsana::TMuon*> goodishMuons,goodMuons;
      vector<const higgsana::TElectron*> goodishElectrons,goodElectrons;
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const higgsana::TMuon *mu = (higgsana::TMuon*)((*muonArr)[i]);
        
        goodishMuons.push_back(mu);

        if(!(mu->pt > 35)) continue;
        if (ComputeMuonPFIso04(mu,pfcandidateArr,rhoMuonIso,MuonEAEra)/mu->pt > 0.15) continue;

 	ULong_t bits = info->triggerBits;
 	ULong_t mubits = mu->hltMatchBits;

//         cout << "Trigger: " << info->triggerBits << " " << bits << " : " << (bits & kHLT_IsoMu24) << " " << Bool_t(bits & kHLT_IsoMu24) << endl;
//         cout << "Object: " << mu->hltMatchBits << " " << mubits << " : " 
//              << (mubits & kHLTObject_IsoMu24 ) << " " << Bool_t(mubits & kHLTObject_IsoMu24) << " : " 
//              << (mubits & kHLTObject_IsoMu30 ) << " " << Bool_t(mubits & kHLTObject_IsoMu30) << " : " 
//              << (mubits & kHLTObject_IsoMu40 ) << " " << Bool_t(mubits & kHLTObject_IsoMu40) << " : " 
//              << endl;
        
        //turn off trigger checking for 2012B samples temporarily - there's a bug with the trigger bits

        //if (std::string::npos != inputFilename.find("r12b"))

        if (!(inputFilename == "WToMuNuPlusFake2012B" || std::string::npos != inputFilename.find("r12b-smu"))) {          
          if(!(
               (bits & kHLT_Mu24	   &&    mubits & kHLTObject_Mu24        &&				mu->pt > 24)    ||
               (bits & kHLT_Mu30	   &&    mubits & kHLTObject_Mu30        &&				mu->pt > 30)    ||
               (bits & kHLT_Mu15	   &&    mubits & kHLTObject_Mu15        &&				mu->pt > 15)    ||
               (bits & kHLT_IsoMu12    &&    mubits & kHLTObject_IsoMu12     &&				mu->pt > 12)    ||
               (bits & kHLT_IsoMu17    &&    mubits & kHLTObject_IsoMu17     &&				mu->pt > 17)    ||
               (bits & kHLT_IsoMu24    &&    mubits & kHLTObject_IsoMu24     &&				mu->pt > 24)    ||
               (bits & kHLT_IsoMu30    &&    mubits & (kHLTObject_IsoMu30 | kHLTObject_IsoMu34)    &&	mu->pt > 30)    ||
               (bits & kHLT_IsoMu24    &&    mubits & kHLTObject_IsoMu30     &&				mu->pt > 30)    ||
               (bits & kHLT_IsoMu24    &&    mubits & kHLTObject_IsoMu34     &&				mu->pt > 34)    ||
               (bits & kHLT_IsoMu24    &&    mubits & kHLTObject_IsoMu40     &&				mu->pt > 40)    
               ))
            continue;
        }

	goodMuons.push_back(mu);
      }
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);

        goodishElectrons.push_back(ele);
	ULong_t bits = info->triggerBits;
	ULong_t elebits = ele->hltMatchBits;

	if(!(
	     (bits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT     &&     elebits & kHLTObject_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT     &&    ele->pt > 27) ||
	     (bits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT     &&     elebits & kHLTObject_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT     &&    ele->pt > 32) ||
	     (bits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT     &&     elebits & kHLTObject_Ele32_WP70     &&    ele->pt > 32) ||
	     (bits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT     &&     elebits & kHLTObject_Ele32_WP80     &&    ele->pt > 32) ||
	     (bits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT     &&     elebits & kHLTObject_Ele27_WP80     &&    ele->pt > 27) 
	     ))
	  continue;
        if(!(ele->pt > 35)) continue;
        if(passCutBasedEleID(ele,ComputeElePFIso04(ele,pfcandidateArr,rhoEleIso,EleEAEra))) goodElectrons.push_back(ele);
      }

      TLorentzVector lepton;
      if(goodMuons.size()==1)          lepton.SetPtEtaPhiM(goodMuons[0]->pt,goodMuons[0]->eta,goodMuons[0]->phi,0.105658);
      else if(goodElectrons.size()==1) lepton.SetPtEtaPhiM(goodElectrons[0]->pt,goodElectrons[0]->eta,goodElectrons[0]->phi,0.000511);
      else continue;
      double dphimetlep = higgsana::deltaPhi(pfMet.Phi(),lepton.Phi());
      double mt = sqrt((2*lepton.Pt()*pfMet.Pt()*(1.0 - cos(dphimetlep))));
      if(mt<50) continue;


      //require also met
      //make tighter cut
      //if(met<35) continue;
      

      //******************************************************************************
      // z veto
      //******************************************************************************
      TLorentzVector looselepton;
      bool veto=false;
      if(goodMuons.size()==1) {
	for(unsigned imu=0; imu<goodishMuons.size(); imu++) {
	  const higgsana::TMuon *mu = goodishMuons[imu];
	  if(goodMuons[0] == mu)       continue;
	  if(goodMuons[0]->q == mu->q) continue;
	  looselepton.SetPtEtaPhiM(mu->pt,mu->eta,mu->phi,0.10568);
	  TLorentzVector dilepton(lepton+looselepton);
	  if(dilepton.M() > 75.0 && dilepton.M() < 105.0) veto=true;
	}
      } else if(goodElectrons.size()==1) {
	for(unsigned iele=0; iele<goodishElectrons.size(); iele++) {
	  const higgsana::TElectron *ele = goodishElectrons[iele];
	  if(goodElectrons[0] == ele)       continue;
	  if(goodElectrons[0]->q == ele->q) continue;
	  looselepton.SetPtEtaPhiM(ele->pt,ele->eta,ele->phi,0.000511);
	  TLorentzVector dilepton(lepton+looselepton);
	  if(dilepton.M() > 75.0 && dilepton.M() < 105.0) veto=true;
	}
      } else continue;
      if(veto) continue;



      //******************************************************************************
      //loop over electrons 
      //******************************************************************************
      bool used=false;
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);

	// Z veto
	TLorentzVector elev;
	elev.SetPtEtaPhiM(ele->pt,ele->eta,ele->phi,0.51099892e-3);
	TLorentzVector dilepton(lepton+elev);
	if(dilepton.M() > 75.0 && dilepton.M() < 105.0) continue;
	//if(dilepton.M() < 30) continue;

        //make cut on dz
        if(fabs(ele->dz) > 0.1) continue;
	// make sure it's not close to the W candidate
	if(elev.DeltaR(lepton) < 0.5) continue;
	
	used = true;
	
        NElectronsFilled++;
        FillElectronTree( eleTree, ele, pfcandidateArr, rhoEleIso, EleEAEra, 
                          info->nPV0, info->runNum, info->lumiSec, info->evtNum);
         
      }
    }

    cout << "Total Electrons: " << NElectronsFilled << endl;

  } //end loop over files

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;
  delete photonArr;
  delete pfcandidateArr;


  cout << "Total Electrons: " << NElectronsFilled << endl;
  outputFile->Write();
  outputFile->Close();
  delete outputFile;
  if (eleTree) delete eleTree;
  

  gBenchmark->Show("WWTemplate");

} 
