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
#include "HiggsAna/DataTree/interface/TMuon.hh"
#include "HiggsAna/DataTree/interface/TJet.hh"
#include "HiggsAna/DataTree/interface/TGenInfo.hh"
#include "HiggsAna/DataTree/interface/TGenParticle.hh"
#include "HiggsAna/DataTree/interface/TPFCandidate.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// Helper functions for Electron ID selection
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"
#include "HiggsAna/Utils/CommonDefs.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"
#include "HiggsAna/Utils/EfficiencyUtils.hh"
#include "HiggsAna/HWW/Utils/LeptonSelection.hh"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

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
Bool_t passElectronDenominatorCuts(ULong_t triggerBits, const higgsana::TElectron *ele, Int_t DenominatorType, string SampleType);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);
void DoComputeElectronFakeRate(const string inputFilename,
                               const string label, 
                               const string outputFilename,  
                               const string smurfOutputFilename, Int_t Option = 0);

//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
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
void ComputeElectronFakeRate_Data() {


//************************************
//MVA 
//************************************

//------------------
//Fully Inclusive
//------------------
    DoComputeElectronFakeRate("2012Data","ElectronFakeRate","ElectronFakeRate.HWWICHEP2012WP.root", "FakeRates_Electron_HWWICHEP2012", 20);
//   DoComputeElectronFakeRate("2012Data","ElectronFakeRate","ElectronFakeRate.HWWIDIsoMVAV3WP.root", "FakeRates_Electron_HWWIDIsoMVAV3WP", 21);
//     DoComputeElectronFakeRate("2012Data","ElectronFakeRate","ElectronFakeRate.HWWIDIsoMVAV4WP.root", "FakeRates_Electron_HWWIDIsoMVAV4WP", 22);
 
 
}



void DoComputeElectronFakeRate(const string inputFilename,
                               const string label, 
                               const string outputFilename, 
                               const string smurfOutputFilename,
                               Int_t Option)
{  
  gBenchmark->Start("WWTemplate");

  //*****************************************************************************************
  //Setup
  //***************************************************************************************** 
  EGammaMvaEleEstimator *eleIDMVA = new EGammaMvaEleEstimator();
  vector<string> weightFiles;
  
  if (Option == 20) {
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_TrigV0_Cat1.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_TrigV0_Cat2.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_TrigV0_Cat3.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_TrigV0_Cat4.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_TrigV0_Cat5.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/MitPhysics/data/ElectronMVAWeights/Electrons_BDTG_TrigV0_Cat6.weights.xml")).c_str());  
    eleIDMVA->initialize( "BDT", EGammaMvaEleEstimator::kTrig, kTRUE, weightFiles);
  }
  if (Option == 21) {
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HWW/data/ElectronIDMVAWeights/ElectronIDMVA_Trig_V3_EtaBin0LowPt_V3_BDTG.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HWW/data/ElectronIDMVAWeights/ElectronIDMVA_Trig_V3_EtaBin1LowPt_V3_BDTG.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HWW/data/ElectronIDMVAWeights/ElectronIDMVA_Trig_V3_EtaBin2LowPt_V3_BDTG.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HWW/data/ElectronIDMVAWeights/ElectronIDMVA_Trig_V3_EtaBin0HighPt_V3_BDTG.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HWW/data/ElectronIDMVAWeights/ElectronIDMVA_Trig_V3_EtaBin1HighPt_V3_BDTG.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HWW/data/ElectronIDMVAWeights/ElectronIDMVA_Trig_V3_EtaBin2HighPt_V3_BDTG.weights.xml")).c_str()); 
    eleIDMVA->initialize( "BDT", EGammaMvaEleEstimator::kTrigIDIsoCombinedPUCorrected, kTRUE, weightFiles);
  }
  if (Option == 22) {
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HWW/data/ElectronIDMVAWeights/ElectronIDMVA_Trig_V4_EtaBin0LowPt_V4_BDTG.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HWW/data/ElectronIDMVAWeights/ElectronIDMVA_Trig_V4_EtaBin1LowPt_V4_BDTG.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HWW/data/ElectronIDMVAWeights/ElectronIDMVA_Trig_V4_EtaBin2LowPt_V4_BDTG.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HWW/data/ElectronIDMVAWeights/ElectronIDMVA_Trig_V4_EtaBin0HighPt_V4_BDTG.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HWW/data/ElectronIDMVAWeights/ElectronIDMVA_Trig_V4_EtaBin1HighPt_V4_BDTG.weights.xml")).c_str());
    weightFiles.push_back((string(getenv("CMSSW_BASE"))+string("/src/HiggsAna/HWW/data/ElectronIDMVAWeights/ElectronIDMVA_Trig_V4_EtaBin2HighPt_V4_BDTG.weights.xml")).c_str());
    eleIDMVA->initialize( "BDT", EGammaMvaEleEstimator::kTrigIDIsoCombined, kTRUE, weightFiles);
  }

  //*****************************************************************************************
  //Define Pt bins
  //*****************************************************************************************
  vector<double> ptbins;
  ptbins.push_back(10);  
  ptbins.push_back(12.5);  
  ptbins.push_back(15);  
  ptbins.push_back(17.5);  
  ptbins.push_back(20);  
  ptbins.push_back(22.5);  
  ptbins.push_back(25);  
  ptbins.push_back(27.5);  
  ptbins.push_back(30);  
  ptbins.push_back(32.5);  
  ptbins.push_back(35);  
  ptbins.push_back(37.5);  
  ptbins.push_back(40);  
  ptbins.push_back(50);  
  ptbins.push_back(60);  
  ptbins.push_back(70);  
  ptbins.push_back(80);  
  ptbins.push_back(100);  


  vector<double> etabins;
  etabins.push_back(0.0);
  etabins.push_back(0.25);
  etabins.push_back(0.5);
  etabins.push_back(0.75);
  etabins.push_back(1.0);
  etabins.push_back(1.25);
  etabins.push_back(1.479);
  etabins.push_back(1.75);
  etabins.push_back(2.0);
  etabins.push_back(2.25);
  etabins.push_back(2.5);
  etabins.push_back(3.0);


  vector<double> phibins;
  phibins.push_back(-3.25);
  phibins.push_back(-2.75);
  phibins.push_back(-2.25);
  phibins.push_back(-1.75);
  phibins.push_back(-1.25);
  phibins.push_back(-0.75);
  phibins.push_back(-0.25);
  phibins.push_back(0.25);
  phibins.push_back(0.75);
  phibins.push_back(1.25);
  phibins.push_back(1.75);
  phibins.push_back(2.25);
  phibins.push_back(2.75);
  phibins.push_back(3.25);

  vector<double> nvtxbins;
  nvtxbins.push_back(0);
  nvtxbins.push_back(2);
  nvtxbins.push_back(4);
  nvtxbins.push_back(6);
  nvtxbins.push_back(8);
  nvtxbins.push_back(10);
  nvtxbins.push_back(12);
  nvtxbins.push_back(14);
  nvtxbins.push_back(16);
  nvtxbins.push_back(18);
  nvtxbins.push_back(20);
  nvtxbins.push_back(22);
  nvtxbins.push_back(24);
  nvtxbins.push_back(26);



  vector<double> rhobins;
  rhobins.push_back(0);
  rhobins.push_back(2);
  rhobins.push_back(4);
  rhobins.push_back(6);
  rhobins.push_back(8);
  rhobins.push_back(10);
  rhobins.push_back(12);
  rhobins.push_back(14);
  rhobins.push_back(16);
  rhobins.push_back(18);
  rhobins.push_back(20);
  rhobins.push_back(22);
  rhobins.push_back(24);
  rhobins.push_back(26);



  vector<double> ptbins2D;
  ptbins2D.push_back(10);  
  ptbins2D.push_back(15);  
  ptbins2D.push_back(20);  
  ptbins2D.push_back(25);  
  ptbins2D.push_back(30);  
  ptbins2D.push_back(35);  
  vector<double> etabins2D;
  etabins2D.push_back(0.0);
  etabins2D.push_back(1.0);
  etabins2D.push_back(1.479);
  etabins2D.push_back(2.0);
  etabins2D.push_back(2.5);


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *histLeadingJetPt = new TH1F("leadingJetPt" , "; p_{T} [GeV/c] ; Number of Events ",  200, 0 , 200);

  //3D array, indices give: [denominatorType][SampleType][ptThreshold]

  vector<vector<vector<TH1F*> > > DenominatorVector_Pt;
  vector<vector<vector<TH1F*> > > DenominatorVector_Eta;
  vector<vector<vector<TH1F*> > > DenominatorVector_Phi;
  vector<vector<vector<TH1F*> > > DenominatorVector_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Rho;
  vector<vector<vector<TH2F*> > > DenominatorVector_PtEta;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt;
  vector<vector<vector<TH1F*> > > NumeratorVector_Eta;
  vector<vector<vector<TH1F*> > > NumeratorVector_Phi;
  vector<vector<vector<TH1F*> > > NumeratorVector_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Rho;
  vector<vector<vector<TH2F*> > > NumeratorVector_PtEta;
  vector<vector<vector<TH1F*> > > LeptonJetPt;
  vector<vector<vector<TH1F*> > > DenominatorIsolation;

  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Barrel_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Endcap_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Barrel_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Endcap_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Barrel_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Endcap_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Barrel_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Endcap_Rho;


  vector<Double_t> denominatorType;
  denominatorType.push_back(4);
  vector<string> sampleLabel;
  sampleLabel.push_back("QCDTightTrigCombinedSample");
  sampleLabel.push_back("QCDTightPlusJetTrigCombinedSample");
  sampleLabel.push_back("QCDCombinedSample");
  vector<Double_t> ptThreshold;
  ptThreshold.push_back(0);
  ptThreshold.push_back(15);
  ptThreshold.push_back(20);
  ptThreshold.push_back(25);
  ptThreshold.push_back(30);
  ptThreshold.push_back(35);
  ptThreshold.push_back(40);
  ptThreshold.push_back(45);
  ptThreshold.push_back(50);
  
  for (UInt_t denominatorTypeIndex = 0; denominatorTypeIndex < denominatorType.size(); ++denominatorTypeIndex) {
    vector<vector<TH1F*> > tmpDenominatorVector_Pt;
    vector<vector<TH1F*> > tmpDenominatorVector_Eta;
    vector<vector<TH1F*> > tmpDenominatorVector_Phi;
    vector<vector<TH1F*> > tmpDenominatorVector_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Rho;
    vector<vector<TH2F*> > tmpDenominatorVector_PtEta;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt;
    vector<vector<TH1F*> > tmpNumeratorVector_Eta;
    vector<vector<TH1F*> > tmpNumeratorVector_Phi;
    vector<vector<TH1F*> > tmpNumeratorVector_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Rho;
    vector<vector<TH2F*> > tmpNumeratorVector_PtEta;
    vector<vector<TH1F*> > tmpLeptonJetPt;
    vector<vector<TH1F*> > tmpDenominatorIsolation;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Barrel_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Barrel_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Endcap_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Endcap_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Barrel_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Barrel_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Endcap_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Endcap_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Barrel_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Barrel_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Endcap_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Endcap_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Barrel_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Barrel_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Endcap_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Endcap_Rho;


    for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
      vector<TH1F*>  tmptmpDenominatorVector_Pt;
      vector<TH1F*>  tmptmpDenominatorVector_Eta;
      vector<TH1F*>  tmptmpDenominatorVector_Phi;
      vector<TH1F*>  tmptmpDenominatorVector_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Rho;
      vector<TH2F*>  tmptmpDenominatorVector_PtEta;
      vector<TH1F*>  tmptmpNumeratorVector_Pt;
      vector<TH1F*>  tmptmpNumeratorVector_Eta;
      vector<TH1F*>  tmptmpNumeratorVector_Phi;
      vector<TH1F*>  tmptmpNumeratorVector_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Rho;
      vector<TH2F*>  tmptmpNumeratorVector_PtEta;
      vector<TH1F*>  tmptmpLeptonJetPt;
      vector<TH1F*>  tmptmpDenominatorIsolation;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Barrel_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Barrel_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Endcap_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Endcap_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Barrel_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Barrel_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Endcap_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Endcap_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho;
      for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {
        TH1F *histDenominator_Pt = new TH1F(("histDenominator_Pt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str(), "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histDenominator_Eta = new TH1F(("histDenominator_Eta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, 0 , 3.0);
        TH1F *histDenominator_Phi = new TH1F(("histDenominator_Phi_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH1F *histDenominator_NVtx = new TH1F(("histDenominator_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Rho = new TH1F(("histDenominator_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH2F *histDenominator_PtEta = new TH2F(("histDenominator_PtEta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, 0, 3.0);
        TH1F *histNumerator_Pt = new TH1F(("histNumerator_Pt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histNumerator_Eta = new TH1F(("histNumerator_Eta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, 0.0 , 3.0);
        TH1F *histNumerator_Phi = new TH1F(("histNumerator_Phi_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH1F *histNumerator_NVtx = new TH1F(("histNumerator_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Rho = new TH1F(("histNumerator_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH2F *histNumerator_PtEta = new TH2F(("histNumerator_PtEta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, 0, 3.0);        
        TH1F *histLeptonJetPt = new TH1F(("histLeptonJetPt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histDenominatorIsolation = new TH1F(("histDenominatorIsolation_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; PF RelIso ; Number of Events ",  100, 0 , 1.0);

        TH1F *histDenominator_Pt10To20_Barrel_NVtx = new TH1F(("histDenominator_Pt10To20_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt10To20_Barrel_Rho = new TH1F(("histDenominator_Pt10To20_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt10To20_Endcap_NVtx = new TH1F(("histDenominator_Pt10To20_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt10To20_Endcap_Rho = new TH1F(("histDenominator_Pt10To20_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt20ToInf_Barrel_NVtx = new TH1F(("histDenominator_Pt20ToInf_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt20ToInf_Barrel_Rho = new TH1F(("histDenominator_Pt20ToInf_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt20ToInf_Endcap_NVtx = new TH1F(("histDenominator_Pt20ToInf_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt20ToInf_Endcap_Rho = new TH1F(("histDenominator_Pt20ToInf_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt10To20_Barrel_NVtx = new TH1F(("histNumerator_Pt10To20_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt10To20_Barrel_Rho = new TH1F(("histNumerator_Pt10To20_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt10To20_Endcap_NVtx = new TH1F(("histNumerator_Pt10To20_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt10To20_Endcap_Rho = new TH1F(("histNumerator_Pt10To20_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt20ToInf_Barrel_NVtx = new TH1F(("histNumerator_Pt20ToInf_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt20ToInf_Barrel_Rho = new TH1F(("histNumerator_Pt20ToInf_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt20ToInf_Endcap_NVtx = new TH1F(("histNumerator_Pt20ToInf_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt20ToInf_Endcap_Rho = new TH1F(("histNumerator_Pt20ToInf_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);



        tmptmpDenominatorVector_Pt.push_back(histDenominator_Pt);
        tmptmpDenominatorVector_Eta.push_back(histDenominator_Eta);
        tmptmpDenominatorVector_Phi.push_back(histDenominator_Phi);
        tmptmpDenominatorVector_NVtx.push_back(histDenominator_NVtx);
        tmptmpDenominatorVector_Rho.push_back(histDenominator_Rho);
        tmptmpDenominatorVector_PtEta.push_back(histDenominator_PtEta);
        tmptmpNumeratorVector_Pt.push_back(histNumerator_Pt);
        tmptmpNumeratorVector_Eta.push_back(histNumerator_Eta);
        tmptmpNumeratorVector_Phi.push_back(histNumerator_Phi);
        tmptmpNumeratorVector_NVtx.push_back(histNumerator_NVtx);
        tmptmpNumeratorVector_Rho.push_back(histNumerator_Rho);
        tmptmpNumeratorVector_PtEta.push_back(histNumerator_PtEta);
        tmptmpLeptonJetPt.push_back(histLeptonJetPt);
        tmptmpDenominatorIsolation.push_back(histDenominatorIsolation);
        tmptmpDenominatorVector_Pt10To20_Barrel_NVtx.push_back(histDenominator_Pt10To20_Barrel_NVtx);
        tmptmpDenominatorVector_Pt10To20_Barrel_Rho.push_back(histDenominator_Pt10To20_Barrel_Rho);
        tmptmpDenominatorVector_Pt10To20_Endcap_NVtx.push_back(histDenominator_Pt10To20_Endcap_NVtx);
        tmptmpDenominatorVector_Pt10To20_Endcap_Rho.push_back(histDenominator_Pt10To20_Endcap_Rho);
        tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(histDenominator_Pt20ToInf_Barrel_NVtx);
        tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho.push_back(histDenominator_Pt20ToInf_Barrel_Rho);
        tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(histDenominator_Pt20ToInf_Endcap_NVtx);
        tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho.push_back(histDenominator_Pt20ToInf_Endcap_Rho);
        tmptmpNumeratorVector_Pt10To20_Barrel_NVtx.push_back(histNumerator_Pt10To20_Barrel_NVtx);
        tmptmpNumeratorVector_Pt10To20_Barrel_Rho.push_back(histNumerator_Pt10To20_Barrel_Rho);
        tmptmpNumeratorVector_Pt10To20_Endcap_NVtx.push_back(histNumerator_Pt10To20_Endcap_NVtx);
        tmptmpNumeratorVector_Pt10To20_Endcap_Rho.push_back(histNumerator_Pt10To20_Endcap_Rho);
        tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(histNumerator_Pt20ToInf_Barrel_NVtx);
        tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho.push_back(histNumerator_Pt20ToInf_Barrel_Rho);
        tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(histNumerator_Pt20ToInf_Endcap_NVtx);
        tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho.push_back(histNumerator_Pt20ToInf_Endcap_Rho);

      }
      tmpDenominatorVector_Pt.push_back(tmptmpDenominatorVector_Pt);
      tmpDenominatorVector_Eta.push_back(tmptmpDenominatorVector_Eta);
      tmpDenominatorVector_Phi.push_back(tmptmpDenominatorVector_Phi);
      tmpDenominatorVector_NVtx.push_back(tmptmpDenominatorVector_NVtx);
      tmpDenominatorVector_Rho.push_back(tmptmpDenominatorVector_Rho);
      tmpDenominatorVector_PtEta.push_back(tmptmpDenominatorVector_PtEta);
      tmpNumeratorVector_Pt.push_back(tmptmpNumeratorVector_Pt);
      tmpNumeratorVector_Eta.push_back(tmptmpNumeratorVector_Eta);
      tmpNumeratorVector_Phi.push_back(tmptmpNumeratorVector_Phi);
      tmpNumeratorVector_NVtx.push_back(tmptmpNumeratorVector_NVtx);
      tmpNumeratorVector_Rho.push_back(tmptmpNumeratorVector_Rho);
      tmpNumeratorVector_PtEta.push_back(tmptmpNumeratorVector_PtEta);
      tmpLeptonJetPt.push_back(tmptmpLeptonJetPt);
      tmpDenominatorIsolation.push_back(tmptmpDenominatorIsolation);
      tmpDenominatorVector_Pt10To20_Barrel_NVtx.push_back(tmptmpDenominatorVector_Pt10To20_Barrel_NVtx);
      tmpDenominatorVector_Pt10To20_Barrel_Rho.push_back(tmptmpDenominatorVector_Pt10To20_Barrel_Rho);
      tmpDenominatorVector_Pt10To20_Endcap_NVtx.push_back(tmptmpDenominatorVector_Pt10To20_Endcap_NVtx);
      tmpDenominatorVector_Pt10To20_Endcap_Rho.push_back(tmptmpDenominatorVector_Pt10To20_Endcap_Rho);
      tmpDenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx);
      tmpDenominatorVector_Pt20ToInf_Barrel_Rho.push_back(tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho);
      tmpDenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx);
      tmpDenominatorVector_Pt20ToInf_Endcap_Rho.push_back(tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho);
      tmpNumeratorVector_Pt10To20_Barrel_NVtx.push_back(tmptmpNumeratorVector_Pt10To20_Barrel_NVtx);
      tmpNumeratorVector_Pt10To20_Barrel_Rho.push_back(tmptmpNumeratorVector_Pt10To20_Barrel_Rho);
      tmpNumeratorVector_Pt10To20_Endcap_NVtx.push_back(tmptmpNumeratorVector_Pt10To20_Endcap_NVtx);
      tmpNumeratorVector_Pt10To20_Endcap_Rho.push_back(tmptmpNumeratorVector_Pt10To20_Endcap_Rho);
      tmpNumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx);
      tmpNumeratorVector_Pt20ToInf_Barrel_Rho.push_back(tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho);
      tmpNumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx);
      tmpNumeratorVector_Pt20ToInf_Endcap_Rho.push_back(tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho);
    }
    DenominatorVector_Pt.push_back(tmpDenominatorVector_Pt);
    DenominatorVector_Eta.push_back(tmpDenominatorVector_Eta);
    DenominatorVector_Phi.push_back(tmpDenominatorVector_Phi);
    DenominatorVector_NVtx.push_back(tmpDenominatorVector_NVtx);
    DenominatorVector_Rho.push_back(tmpDenominatorVector_Rho);
    DenominatorVector_PtEta.push_back(tmpDenominatorVector_PtEta);
    NumeratorVector_Pt.push_back(tmpNumeratorVector_Pt);
    NumeratorVector_Eta.push_back(tmpNumeratorVector_Eta);
    NumeratorVector_Phi.push_back(tmpNumeratorVector_Phi);
    NumeratorVector_NVtx.push_back(tmpNumeratorVector_NVtx);
    NumeratorVector_Rho.push_back(tmpNumeratorVector_Rho);
    NumeratorVector_PtEta.push_back(tmpNumeratorVector_PtEta);
    LeptonJetPt.push_back(tmpLeptonJetPt);
    DenominatorIsolation.push_back(tmpDenominatorIsolation);
    DenominatorVector_Pt10To20_Barrel_NVtx.push_back(tmpDenominatorVector_Pt10To20_Barrel_NVtx);
    DenominatorVector_Pt10To20_Barrel_Rho.push_back(tmpDenominatorVector_Pt10To20_Barrel_Rho);
    DenominatorVector_Pt10To20_Endcap_NVtx.push_back(tmpDenominatorVector_Pt10To20_Endcap_NVtx);
    DenominatorVector_Pt10To20_Endcap_Rho.push_back(tmpDenominatorVector_Pt10To20_Endcap_Rho);
    DenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(tmpDenominatorVector_Pt20ToInf_Barrel_NVtx);
    DenominatorVector_Pt20ToInf_Barrel_Rho.push_back(tmpDenominatorVector_Pt20ToInf_Barrel_Rho);
    DenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(tmpDenominatorVector_Pt20ToInf_Endcap_NVtx);
    DenominatorVector_Pt20ToInf_Endcap_Rho.push_back(tmpDenominatorVector_Pt20ToInf_Endcap_Rho);
    NumeratorVector_Pt10To20_Barrel_NVtx.push_back(tmpNumeratorVector_Pt10To20_Barrel_NVtx);
    NumeratorVector_Pt10To20_Barrel_Rho.push_back(tmpNumeratorVector_Pt10To20_Barrel_Rho);
    NumeratorVector_Pt10To20_Endcap_NVtx.push_back(tmpNumeratorVector_Pt10To20_Endcap_NVtx);
    NumeratorVector_Pt10To20_Endcap_Rho.push_back(tmpNumeratorVector_Pt10To20_Endcap_Rho);
    NumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(tmpNumeratorVector_Pt20ToInf_Barrel_NVtx);
    NumeratorVector_Pt20ToInf_Barrel_Rho.push_back(tmpNumeratorVector_Pt20ToInf_Barrel_Rho);
    NumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(tmpNumeratorVector_Pt20ToInf_Endcap_NVtx);
    NumeratorVector_Pt20ToInf_Endcap_Rho.push_back(tmpNumeratorVector_Pt20ToInf_Endcap_Rho);
  }

  ofstream eventListFile("fakeEventList.ComputeFakeRate.txt");

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
  rlrm.AddJSONFile("/afs/cern.ch/work/s/sixie/public/HZZ4l/auxiliar/2012/Cert_190456-196531_8TeV_PromptReco_Collisions12_JSON.txt"); 

  Int_t NEvents = 0;


  UInt_t DataEra = kDataEra_NONE;


  vector<string> inputfiles;
  if (inputFilename == "2012Data") {
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r12a-del-pr-v1_FakeRateTriggerAndDenominatorSkimmed.root");
    inputfiles.push_back("/data/smurf/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r12b-del-pr-v1_FakeRateTriggerAndDenominatorSkimmed.root");
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
    eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Photon", &photonArr);     photonBr = eventTree->GetBranch("Photon");
    eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");
    eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr); pfcandidateBr = eventTree->GetBranch("PFCandidate"); 

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
      
      //For MVA, only use odd event numbers because even event numbers were used for training
      if (Option >= 20 && Option < 30) {
        if (info->evtNum % 2 == 0) continue;
      }

      Double_t eventweight = info->eventweight;

      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

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
      // TcMet
      //********************************************************
      TVector3 pfMet;        
      if(info->pfMEx!=0 || info->pfMEy!=0) {       
        pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
      }
      Double_t met = pfMet.Pt();

      Int_t NElectrons = electronArr->GetEntries();
      
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const higgsana::TElectron *ele = (higgsana::TElectron*)((*electronArr)[i]);

        Double_t leadingJetPt = -1;
        Double_t leptonJetPt = -1;
        //pass event selection     
        for(Int_t j=0; j<jetArr->GetEntries(); j++) {
          const higgsana::TJet *jet = (higgsana::TJet*)((*jetArr)[j]);        
          if (jet->pt > leadingJetPt &&
              higgsana::deltaR(jet->eta, jet->phi, ele->eta, ele->phi) > 1.0) {
            leadingJetPt = jet->pt;          
          }
          if (higgsana::deltaR(jet->eta, jet->phi, ele->eta, ele->phi) < 0.5) {
            leptonJetPt = jet->pt;
          }
        }
      
        //if there's no jet ( == isolated lepton?) then take pt of lepton
        if (leptonJetPt < 0) {
          leptonJetPt = ele->pt;
        }


        for ( UInt_t denominatorTypeIndex = 0 ; denominatorTypeIndex < denominatorType.size() ; ++denominatorTypeIndex ) {
          for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
            for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {
          
              //********************************************************
              // Event Selection Cuts
              //********************************************************

              if (NElectrons > 1) continue;

              if (sampleLabel[sampleTypeIndex] == "QCDCombinedSample" 
                  || sampleLabel[sampleTypeIndex] == "QCDTightTrigCombinedSample"
                  || sampleLabel[sampleTypeIndex] == "QCDTightPlusJetTrigCombinedSample"
                ) {
                if (met > 20) continue;

                Bool_t passJetSelection = kFALSE;
                if (ptThreshold[ptThresholdIndex] == 0) passJetSelection = kTRUE;
                if (leadingJetPt > ptThreshold[ptThresholdIndex]) {
                  passJetSelection = kTRUE;               
                }
                if (!passJetSelection) continue;
              }

              if (!(ele->pt > 10 && ele->pt < 35.0)) continue;

              if (passElectronDenominatorCuts(info->triggerBits, ele, denominatorType[denominatorTypeIndex], sampleLabel[sampleTypeIndex])) {
              
                if (denominatorTypeIndex == 0 && sampleLabel[sampleTypeIndex] == "QCDTightTrigCombinedSample" && ptThresholdIndex == 0 && Option == 20) {
                  eventListFile << info->runNum << " " << info->lumiSec << " " << info->evtNum << " : " << ele->pt << " " << ele->eta << " : "
                                << Bool_t(PassEleHWWICHEP2012IDMVA( ele, eleIDMVA)) << " "
                                << ComputeElePFIso04(ele, pfcandidateArr, rhoEleIso,EleEAEra) / ele->pt << " " 
                                << endl;

//                   cout << "run lumi event: " << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
//                   PassEleHWWICHEP2012IDMVA( ele, eleIDMVA, kTRUE);

                }

                
                DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt,eventweight);
                DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(ele->eta),eventweight);
                DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->phi,eventweight);
                DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rhoEleIso,eventweight);
                DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt, fabs(ele->eta), eventweight);
               
                if (ele->pt < 20 && fabs(ele->eta) < 1.479) {
                  DenominatorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rhoEleIso,eventweight);
                }
                if (ele->pt < 20 && fabs(ele->eta) >= 1.479) {
                  DenominatorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rhoEleIso,eventweight);
                } 
                if (ele->pt >= 20 && fabs(ele->eta) < 1.479) {
                  DenominatorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rhoEleIso,eventweight);
                } 
                if (ele->pt >= 20 && fabs(ele->eta) >= 1.479) {
                  DenominatorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rhoEleIso,eventweight);
                } 
              
                if (
                  (Option == 20 
                   && 
                   PassEleHWWICHEP2012IDMVA( ele, eleIDMVA) 
                   &&
                   (ComputeElePFIso04(ele, pfcandidateArr, rhoEleIso,EleEAEra) / ele->pt < 0.15)
                    )
                  ||
                  (Option == 21 
                   && 
                   PassEleHWWIDIsoMVAV3( ele,pfcandidateArr,rhoEleIso,EleEAEra,eleIDMVA)
                    )
                  ||
                  (Option == 22  
                   && 
                   PassEleHWWIDIsoMVAV4( ele,pfcandidateArr,rhoEleIso,kDataEra_NONE,eleIDMVA)
                    )
                  ) {

                  NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt, eventweight);
                  NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(ele->eta), eventweight);
                  NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->phi, eventweight);
                  NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rhoEleIso,eventweight);
                  NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt, fabs(ele->eta) , eventweight);   
        
                  if (ele->pt < 20 && fabs(ele->eta) < 1.479) {
                    NumeratorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rhoEleIso,eventweight);
                  }
                  if (ele->pt < 20 && fabs(ele->eta) >= 1.479) {
                    NumeratorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rhoEleIso,eventweight);
                  } 
                  if (ele->pt >= 20 && fabs(ele->eta) < 1.479) {
                    NumeratorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rhoEleIso,eventweight);
                  } 
                  if (ele->pt >= 20 && fabs(ele->eta) >= 1.479) {
                    NumeratorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rhoEleIso,eventweight);
                  }
                }
              }

            } //loop over denominator types
          } //loop over sample types
        } //loop over ptThresholds

      } //loop over electrons

    } //end loop over data    
    cout << "done " << inputfiles[f] << endl;

  } //end loop over files

  eventListFile.close();

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;
  
  
  //*****************************************************************************************
  //Make Efficiency Plots
  //*****************************************************************************************
  
  for (UInt_t denominatorTypeIndex = 0; denominatorTypeIndex < denominatorType.size(); ++denominatorTypeIndex) {
    for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
      for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {

        Bool_t printDebug = kFALSE;
        
        TGraphAsymmErrors *efficiency_pt = higgsana::createEfficiencyGraph(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt", ptbins,  -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_eta = higgsana::createEfficiencyGraph(NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Eta", etabins,  -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_phi = higgsana::createEfficiencyGraph(NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Phi", phibins,  -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_nvtx = higgsana::createEfficiencyGraph(NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_NVtx", nvtxbins,  -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_rho = higgsana::createEfficiencyGraph(NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Rho", rhobins,  -99, -99, 0, 1);
        TH2F *efficiency_PtEta = 
          higgsana::createEfficiencyHist2D(NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                           DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_PtEta", ptbins2D, etabins2D);
        
        
        TGraphAsymmErrors *efficiency_Pt10To20_Barrel_nvtx = higgsana::createEfficiencyGraph(NumeratorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Barrel_NVtx", nvtxbins,  -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Barrel_rho = higgsana::createEfficiencyGraph(NumeratorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Barrel_Rho", rhobins,  -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Endcap_nvtx = higgsana::createEfficiencyGraph(NumeratorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Endcap_NVtx", nvtxbins,  -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Endcap_rho = higgsana::createEfficiencyGraph(NumeratorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Endcap_Rho", rhobins,  -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Barrel_nvtx = higgsana::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Barrel_NVtx", nvtxbins,  -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Barrel_rho = higgsana::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Barrel_Rho", rhobins,  -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Endcap_nvtx = higgsana::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Endcap_NVtx", nvtxbins,  -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Endcap_rho = higgsana::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Endcap_Rho", rhobins, -99, -99, 0, 1);



        TFile *file = TFile::Open(outputFilename.c_str(), "UPDATE");
        file->cd();
        
        file->WriteTObject(efficiency_pt, efficiency_pt->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_eta, efficiency_eta->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_phi, efficiency_phi->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_nvtx, efficiency_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_rho, efficiency_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_PtEta, efficiency_PtEta->GetName(), "WriteDelete");

        cout << denominatorTypeIndex << " " << sampleTypeIndex << " " << ptThresholdIndex << " : " << endl;
        cout << Bool_t(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]) << endl;

//         file->WriteTObject(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        
        file->WriteTObject(efficiency_Pt10To20_Barrel_nvtx, efficiency_Pt10To20_Barrel_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Barrel_rho, efficiency_Pt10To20_Barrel_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Endcap_nvtx, efficiency_Pt10To20_Endcap_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Endcap_rho, efficiency_Pt10To20_Endcap_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Barrel_nvtx, efficiency_Pt20ToInf_Barrel_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Barrel_rho, efficiency_Pt20ToInf_Barrel_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Endcap_nvtx, efficiency_Pt20ToInf_Endcap_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Endcap_rho, efficiency_Pt20ToInf_Endcap_rho->GetName(), "WriteDelete");

        file->Close();
        
        //*****************************************************
        // Combined Fake Rate
        //*****************************************************
        if (denominatorType[denominatorTypeIndex] == 4 && 
            (sampleLabel[sampleTypeIndex] == "QCDTightPlusJetTrigCombinedSample"
              ) 
          ) {
          
          TH2F *eff = 0;
          
          file = TFile::Open((smurfOutputFilename + ".root").c_str(), "UPDATE");

          cout << Bool_t(NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]) << endl;
          cout << NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName() << endl;
          cout << Bool_t(DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]) << endl;
          cout << DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName() << endl;
          eff = higgsana::createEfficiencyHist2D(NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                 DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                 "ElectronFakeRate_V"+IntToString(denominatorType[denominatorTypeIndex])+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_PtEta", 
                                                          ptbins2D, etabins2D);
          file->WriteTObject(eff, eff->GetName(), "WriteDelete");
          delete file;
        }

      }
    }
  }


  cout << "Total Events: " << NEvents << endl;


  gBenchmark->Show("WWTemplate");       
} 



Bool_t passElectronDenominatorCuts(ULong_t triggerBits, const higgsana::TElectron *ele, Int_t DenominatorType, string SampleType) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(ele->pt > 10 && fabs(ele->eta) < 2.5)) pass = kFALSE;

  //match to HLT

  if (SampleType == "QCDTightTrigCombinedSample") {
    if (!(
          ((triggerBits & kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL)
           && ele->hltMatchBits & kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            )
          ||
          ((triggerBits & kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30)
           && ele->hltMatchBits & kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            )
          ||
          ((triggerBits & kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL)
           && ele->hltMatchBits & kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            )
          ||
          ((triggerBits & kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30)
           && ele->hltMatchBits & kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            )
          )
      ) {
      pass = kFALSE;
    }
  }

  if (SampleType == "QCDTightPlusJetTrigCombinedSample") {
    if (!(
          ((triggerBits & kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30)
           && ele->hltMatchBits & kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            )
          ||
          ((triggerBits & kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30)
           && ele->hltMatchBits & kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            )
          )
      ) {
      pass = kFALSE;
    }
  }

  if (SampleType == "QCDCombinedSample") {
    if (!(
          ((triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL)
           && ele->hltMatchBits & kHLTObject_Ele8_CaloIdL_CaloIsoVL
            )
          ||
          ((triggerBits & kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL)
           && ele->hltMatchBits & kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            )
          ||
          ((triggerBits & kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30)
           && ele->hltMatchBits & kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            )
          ||
          ((triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL)
           && ele->hltMatchBits & kHLTObject_Ele17_CaloIdL_CaloIsoVL
            )
          ||
          ((triggerBits & kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL)
           && ele->hltMatchBits & kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            )
          ||
          ((triggerBits & kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30)
           && ele->hltMatchBits & kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            )
          )
      ) {
      pass = kFALSE;
    }
  }

  if (DenominatorType == 4) {
    if (!passElectronTriggerDenominator(ele)) pass = kFALSE;
  }

  return pass;
}


//--------------------------------------------------------------------------------------------------
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2)
{
  ofs <<   runNum << " " ;
  ofs <<  lumiSec << " ";
  ofs << evtNum<< " ";
  ofs << mass<< " ";

//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
//   ofs << "    pt    |    eta    |    phi    |   iso    |    d0      | ntk | npx | nseg | nval | chi^2/ndf | TM | HLT" << endl;
//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
  ofs << " " ;
  ofs << setw(9) << pt1 << " |";
  ofs << setw(10) << eta1 << " |";
  ofs << setw(10) << phi1 << " |";
  ofs << setw(10) << leptonCharge1 << " |";
  ofs << setw(9) << pt2 << " |";
  ofs << setw(10) << eta2 << " |";
  ofs << setw(10) << phi2 << " |";
  ofs << setw(10) << leptonCharge2 << " |";
  ofs << endl;
  
}
