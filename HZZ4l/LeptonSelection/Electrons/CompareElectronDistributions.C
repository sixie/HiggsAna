//root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt",0)'
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
#include <MyStyle.h>
#include "TLegend.h"
#include "TEfficiency.h"

#include "HiggsAna/CommonData/interface/ElectronTree.h"

// lumi section selection with JSON files
// #include "MitCommon/DataFormats/interface/Types.h"
// #include "MitAna/DataCont/interface/RunLumiRangeMap.h"
// #include "MitCommon/MathTools/interface/MathUtils.h"
// #include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
// #include "MitHiggs/Utils/interface/EfficiencyUtils.h"
// #include "MitHiggs/Utils/interface/PlotUtils.h"
// #include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"
// #include "HiggsAna/Ntupler/interface/HiggsAnaDefs.hh"

#endif


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
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("eleIDdir");
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


//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}





//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void CompareElectronDistributions(string RealElectronFile, string FakeElectronFile, 
                                  string RealElectronLabel, string FakeElectronLabel,
                                  string Label, Int_t Option)
{  

  string label = "";
  if (Label != "") label = "_" + Label;


  //*****************************************************************************************
  //Plotting Setup
  //*****************************************************************************************


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *Ele_Pt_Real = new TH1F(("Ele_Pt_Real"+label).c_str(), "; Pt ; Number of Events ",  100, 0 , 50);
  TH1F *Ele_Eta_Real = new TH1F(("Ele_Eta_Real"+label).c_str(), "; Eta ; Number of Events ",  100, -2.5 , 2.5);
  TH1F *Ele_Rho_Real = new TH1F(("Ele_Rho_Real"+label).c_str(), "; Rho ; Number of Events ",  100, 0 , 50);
  TH1F *Ele_D0_Real = new TH1F(("Ele_D0_Real"+label).c_str(), "; D0 ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_DZ_Real = new TH1F(("Ele_DZ_Real"+label).c_str(), "; DZ ; Number of Events ",  100, 0 , 0.2);
  TH1F *Ele_IP3D_Real = new TH1F(("Ele_IP3D_Real"+label).c_str(), "; IP3D ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_NBrem_Real = new TH1F(("Ele_NBrem_Real"+label).c_str(), "; NBrem ; Number of Events ",  100, 0 , 10);
  TH1F *Ele_FBrem_Real = new TH1F(("Ele_FBrem_Real"+label).c_str(), "; FBrem ; Number of Events ",  100, 0 , 1);
  TH1F *Ele_EoP_Real = new TH1F(("Ele_EoP_Real"+label).c_str(), "; EoP ; Number of Events ",  100, 0 , 2);
  TH1F *Ele_EEleoPout_Real = new TH1F(("Ele_EEleoPout_Real"+label).c_str(), "; EEleoPout ; Number of Events ",  100, 0 , 2);
  TH1F *Ele_IoEmIoP_Real = new TH1F(("Ele_IoEmIoP_Real"+label).c_str(), "; IoEmIoP ; Number of Events ",  100, -0.2 , 0.3);
  TH1F *Ele_DEta_Real = new TH1F(("Ele_DEta_Real"+label).c_str(), "; DEta ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_DPhi_Real = new TH1F(("Ele_DPhi_Real"+label).c_str(), "; DPhi ; Number of Events ",  100, -1 , 1);
  TH1F *Ele_DEtaCalo_Real = new TH1F(("Ele_DEtaCalo_Real"+label).c_str(), "; DEtaCalo ; Number of Events ",  100, -0.2 , 0.2);
  TH1F *Ele_DPhiCalo_Real = new TH1F(("Ele_DPhiCalo_Real"+label).c_str(), "; DPhiCalo ; Number of Events ",  100, -1 , 1);
  TH1F *Ele_SigmaIEtaIEta_Real = new TH1F(("Ele_SigmaIEtaIEta_Real"+label).c_str(), "; SigmaIEtaIEta ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_SigmaIPhiIPhi_Real = new TH1F(("Ele_SigmaIPhiIPhi_Real"+label).c_str(), "; SigmaIPhiIPhi ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_EtaWidth_Real = new TH1F(("Ele_EtaWidth_Real"+label).c_str(), "; EtaWidth ; Number of Events ",  100, 0 , 0.5);
  TH1F *Ele_PhiWidth_Real = new TH1F(("Ele_PhiWidth_Real"+label).c_str(), "; PhiWidth ; Number of Events ",  100, 0 , 0.5);
  TH1F *Ele_R9_Real = new TH1F(("Ele_R9_Real"+label).c_str(), "; R9 ; Number of Events ",  100, 0 , 2);
  TH1F *Ele_PreShowerOverRaw_Real = new TH1F(("Ele_PreShowerOverRaw_Real"+label).c_str(), "; PreShowerOverRaw ; Number of Events ",  100, 0 , 2);
  TH1F *Ele_HoE_Real = new TH1F(("Ele_HoE_Real"+label).c_str(), "; HoE ; Number of Events ",  100, 0 , 2);
  TH1F *Ele_GsfChi2_Real = new TH1F(("Ele_GsfChi2_Real"+label).c_str(), "; GsfChi2 ; Number of Events ",  100, 0 , 200);
  TH1F *Ele_KFChi2_Real = new TH1F(("Ele_KFChi2_Real"+label).c_str(), "; KFChi2 ; Number of Events ",  100, 0 , 10);
  TH1F *Ele_KFLayers_Real = new TH1F(("Ele_KFLayers_Real"+label).c_str(), "; KFLayers ; Number of Events ",  21, -0.5 , 20.5);
  TH1F *Ele_OneMinusSeedE1x5OverE5x5_Real = new TH1F(("Ele_OneMinusSeedE1x5OverE5x5_Real"+label).c_str(), "; OneMinusSeedE1x5OverE5x5 ; Number of Events ",  100, -1.5 , 1.5);
  TH1F *Ele_PFIso04OverPt_Real = new TH1F(("Ele_PFIso04OverPt_Real"+label).c_str(), "; PFIso04OverPt ; Number of Events ",  100, 0 , 5);

  TH1F *Ele_Pt_Fake = new TH1F(("Ele_Pt_Fake"+label).c_str(), "; Pt ; Number of Events ",  100, 0 , 50);
  TH1F *Ele_Eta_Fake = new TH1F(("Ele_Eta_Fake"+label).c_str(), "; Eta ; Number of Events ",  100, -2.5 , 2.5);
  TH1F *Ele_Rho_Fake = new TH1F(("Ele_Rho_Fake"+label).c_str(), "; Rho ; Number of Events ",  100, 0 , 50);
  TH1F *Ele_D0_Fake = new TH1F(("Ele_D0_Fake"+label).c_str(), "; D0 ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_DZ_Fake = new TH1F(("Ele_DZ_Fake"+label).c_str(), "; DZ ; Number of Events ",  100, 0 , 0.2);
  TH1F *Ele_IP3D_Fake = new TH1F(("Ele_IP3D_Fake"+label).c_str(), "; IP3D ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_NBrem_Fake = new TH1F(("Ele_NBrem_Fake"+label).c_str(), "; NBrem ; Number of Events ",  100, 0 , 10);
  TH1F *Ele_FBrem_Fake = new TH1F(("Ele_FBrem_Fake"+label).c_str(), "; FBrem ; Number of Events ",  100, 0 , 1);
  TH1F *Ele_EoP_Fake = new TH1F(("Ele_EoP_Fake"+label).c_str(), "; EoP ; Number of Events ",  100, 0 , 2);
  TH1F *Ele_EEleoPout_Fake = new TH1F(("Ele_EEleoPout_Fake"+label).c_str(), "; EEleoPout ; Number of Events ",  100, 0 , 2);
  TH1F *Ele_IoEmIoP_Fake = new TH1F(("Ele_IoEmIoP_Fake"+label).c_str(), "; IoEmIoP ; Number of Events ",  100, -0.2 , 0.3);
  TH1F *Ele_DEta_Fake = new TH1F(("Ele_DEta_Fake"+label).c_str(), "; DEta ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_DPhi_Fake = new TH1F(("Ele_DPhi_Fake"+label).c_str(), "; DPhi ; Number of Events ",  100, -1 , 1);
  TH1F *Ele_DEtaCalo_Fake = new TH1F(("Ele_DEtaCalo_Fake"+label).c_str(), "; DEtaCalo ; Number of Events ",  100, -0.2 , 0.2);
  TH1F *Ele_DPhiCalo_Fake = new TH1F(("Ele_DPhiCalo_Fake"+label).c_str(), "; DPhiCalo ; Number of Events ",  100, -1 , 1);
  TH1F *Ele_SigmaIEtaIEta_Fake = new TH1F(("Ele_SigmaIEtaIEta_Fake"+label).c_str(), "; SigmaIEtaIEta ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_SigmaIPhiIPhi_Fake = new TH1F(("Ele_SigmaIPhiIPhi_Fake"+label).c_str(), "; SigmaIPhiIPhi ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_EtaWidth_Fake = new TH1F(("Ele_EtaWidth_Fake"+label).c_str(), "; EtaWidth ; Number of Events ",  100, 0 , 0.5);
  TH1F *Ele_PhiWidth_Fake = new TH1F(("Ele_PhiWidth_Fake"+label).c_str(), "; PhiWidth ; Number of Events ",  100, 0 , 0.5);
  TH1F *Ele_R9_Fake = new TH1F(("Ele_R9_Fake"+label).c_str(), "; R9 ; Number of Events ",  100, 0 , 2);
  TH1F *Ele_PreShowerOverRaw_Fake = new TH1F(("Ele_PreShowerOverRaw_Fake"+label).c_str(), "; PreShowerOverRaw ; Number of Events ",  100, 0 , 2);
  TH1F *Ele_HoE_Fake = new TH1F(("Ele_HoE_Fake"+label).c_str(), "; HoE ; Number of Events ",  100, 0 , 2);
  TH1F *Ele_GsfChi2_Fake = new TH1F(("Ele_GsfChi2_Fake"+label).c_str(), "; GsfChi2 ; Number of Events ",  100, 0 , 200);
  TH1F *Ele_KFChi2_Fake = new TH1F(("Ele_KFChi2_Fake"+label).c_str(), "; KFChi2 ; Number of Events ",  100, 0 , 10);
  TH1F *Ele_KFLayers_Fake = new TH1F(("Ele_KFLayers_Fake"+label).c_str(), "; KFLayers ; Number of Events ",  21, -0.5 , 20.5);
  TH1F *Ele_OneMinusSeedE1x5OverE5x5_Fake = new TH1F(("Ele_OneMinusSeedE1x5OverE5x5_Fake"+label).c_str(), "; OneMinusSeedE1x5OverE5x5 ; Number of Events ",  100, -1.5 , 1.5);
  TH1F *Ele_PFIso04OverPt_Fake = new TH1F(("Ele_PFIso04OverPt_Fake"+label).c_str(), "; PFIso04OverPt ; Number of Events ",  100, 0 , 5);



  //*****************************************************************************************
  //Load MVA Branches
  //*****************************************************************************************
  Float_t       fEleBDTGV0;
  Float_t       fEleBDTGV1;
  Float_t       fEleBDTGV2;
  Float_t       fEleBDTGV3;
  Float_t       fEleBDTGV4;
  Float_t       fEleBDTGV5;
  Float_t       fEleBDTGV6;
  Float_t       fEleBDTGV7;
  Float_t       fEleBDTGV8;
  Float_t       fEleBDTGV9;
  Float_t       fEleBDTGV10;
  Float_t       fEleBDTGV11;
  Float_t       fEleBDTGV12;
  Float_t       fEleBDTGV13;
  Float_t       fEleBDTGV14;
  Float_t       fEleBDTGV15;
  Float_t       fEleBDTGV16;
  Float_t       fEleBDTGV17;
  Float_t       fEleBDTGV18;

  //*****************************************************************************************
  //RealEleTree
  //*****************************************************************************************
  ElectronTree RealEleTree;
  RealEleTree.LoadTree(RealElectronFile.c_str());
  RealEleTree.InitTree();
  RealEleTree.tree_->SetBranchAddress( "EleIDMVA_BDTG_V0", &fEleBDTGV0); 
  RealEleTree.tree_->SetBranchAddress( "EleIDMVA_BDTG_V1", &fEleBDTGV1); 
  RealEleTree.tree_->SetBranchAddress( "EleIDMVA_BDTG_V2", &fEleBDTGV2); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV3", &fEleBDTGV3); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV4", &fEleBDTGV4); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV5", &fEleBDTGV5); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV6", &fEleBDTGV6); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV7", &fEleBDTGV7); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV8", &fEleBDTGV8); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV9", &fEleBDTGV9); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV10", &fEleBDTGV10); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV11", &fEleBDTGV11); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV12", &fEleBDTGV12); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV13", &fEleBDTGV13); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV14", &fEleBDTGV14); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV15", &fEleBDTGV15); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV16", &fEleBDTGV16); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV17", &fEleBDTGV17); 
  RealEleTree.tree_->SetBranchAddress( "BDTGV18", &fEleBDTGV18); 

  for(UInt_t ientry=0; ientry < RealEleTree.tree_->GetEntries(); ientry++) {       	
    RealEleTree.tree_->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //don't evaluate performance using training events
    if (RealEleTree.fElePt < 5) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(RealEleTree.fEleEta) < 0.8) subdet = 0;
    else if (fabs(RealEleTree.fEleEta) < 1.485) subdet = 1;
    else subdet = 2;

    Int_t ptBin = 0;
    if (RealEleTree.fElePt > 10.0) ptBin = 1;
    if (RealEleTree.fElePt > 20.0) ptBin = 2;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);    
    if (Option == 6) passCuts = (subdet == 0 && ptBin == 2);
    if (Option == 7) passCuts = (subdet == 1 && ptBin == 2);
    if (Option == 8) passCuts = (subdet == 2 && ptBin == 2);    
    if (Option == 10) passCuts = (ptBin == 0 );

    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

//     //apply denominator cuts
//     if (!(RealEleTree.fEleDenomFake==1 && fabs(RealEleTree.fEleDZ) < 0.1 && fabs(RealEleTree.fEleD0) < 0.02)) continue;
//     if (RealEleTree.fEleMatchedConversion) continue;
    if (!(RealEleTree.fEleNMissHits <= 1)) continue;
    if( RealEleTree.fEleIP3dSig > 4 ) continue;

    
    //Fill Histograms
    Ele_SigmaIEtaIEta_Real->Fill(RealEleTree.fEleSigmaIEtaIEta, RealEleTree.fWeight);
    Ele_Pt_Real->Fill(RealEleTree.fElePt, RealEleTree.fWeight);
    Ele_Eta_Real->Fill(RealEleTree.fEleEta, RealEleTree.fWeight);
    Ele_Rho_Real->Fill(RealEleTree.fRho, RealEleTree.fWeight);
    Ele_D0_Real->Fill(RealEleTree.fEleD0, RealEleTree.fWeight);
    Ele_DZ_Real->Fill(RealEleTree.fEleDZ, RealEleTree.fWeight);
    Ele_IP3D_Real->Fill(RealEleTree.fEleIP3d, RealEleTree.fWeight);
    Ele_NBrem_Real->Fill(RealEleTree.fEleNBrem, RealEleTree.fWeight);
    Ele_FBrem_Real->Fill(RealEleTree.fEleFBrem, RealEleTree.fWeight);
    Ele_EoP_Real->Fill(RealEleTree.fEleEOverP, RealEleTree.fWeight);
    Ele_EEleoPout_Real->Fill(RealEleTree.fEleEEleClusterOverPout, RealEleTree.fWeight);
    Ele_IoEmIoP_Real->Fill(RealEleTree.fEleOneOverEMinusOneOverP, RealEleTree.fWeight);
    Ele_DEta_Real->Fill(RealEleTree.fEleDEtaIn, RealEleTree.fWeight);
    Ele_DPhi_Real->Fill(RealEleTree.fEleDPhiIn, RealEleTree.fWeight);
    Ele_DEtaCalo_Real->Fill(RealEleTree.fEledEtaCalo, RealEleTree.fWeight);
    Ele_DPhiCalo_Real->Fill(RealEleTree.fEledPhiCalo, RealEleTree.fWeight);
    Ele_SigmaIEtaIEta_Real->Fill(RealEleTree.fEleSigmaIEtaIEta, RealEleTree.fWeight);
    Ele_SigmaIPhiIPhi_Real->Fill(RealEleTree.fEleSigmaIPhiIPhi, RealEleTree.fWeight);
    Ele_EtaWidth_Real->Fill(RealEleTree.fEleSCEtaWidth, RealEleTree.fWeight);
    Ele_PhiWidth_Real->Fill(RealEleTree.fEleSCPhiWidth, RealEleTree.fWeight);
    Ele_R9_Real->Fill(RealEleTree.fEleR9, RealEleTree.fWeight);
    Ele_PreShowerOverRaw_Real->Fill(RealEleTree.fElePreShowerOverRaw, RealEleTree.fWeight);
    Ele_HoE_Real->Fill(RealEleTree.fEleHoverE, RealEleTree.fWeight);
    Ele_GsfChi2_Real->Fill(RealEleTree.fEleGsfTrackChi2OverNdof, RealEleTree.fWeight);
    Ele_KFChi2_Real->Fill(RealEleTree.fEleKFTrackChi2OverNDoF, RealEleTree.fWeight);
    Ele_KFLayers_Real->Fill(RealEleTree.fEleKFTrackNLayersWithMeasurement, RealEleTree.fWeight);
    Ele_OneMinusSeedE1x5OverE5x5_Real->Fill(RealEleTree.fEleOneMinusSeedE1x5OverE5x5, RealEleTree.fWeight);
    Ele_PFIso04OverPt_Real->Fill(RealEleTree.fElePFIso04 / RealEleTree.fElePt, RealEleTree.fWeight);


  } 
  





  //*****************************************************************************************
  //FakeEleTree
  //*****************************************************************************************
  ElectronTree FakeEleTree;
  FakeEleTree.LoadTree(FakeElectronFile.c_str());
  FakeEleTree.InitTree();

  FakeEleTree.tree_->SetBranchAddress( "EleIDMVA_BDTG_V0", &fEleBDTGV0); 
  FakeEleTree.tree_->SetBranchAddress( "EleIDMVA_BDTG_V1", &fEleBDTGV1); 
  FakeEleTree.tree_->SetBranchAddress( "EleIDMVA_BDTG_V2", &fEleBDTGV2); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV3", &fEleBDTGV3); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV4", &fEleBDTGV4); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV5", &fEleBDTGV5); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV6", &fEleBDTGV6); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV7", &fEleBDTGV7); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV8", &fEleBDTGV8); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV9", &fEleBDTGV9); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV10", &fEleBDTGV10); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV11", &fEleBDTGV11); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV12", &fEleBDTGV12); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV13", &fEleBDTGV13); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV14", &fEleBDTGV14); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV15", &fEleBDTGV15); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV16", &fEleBDTGV16); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV17", &fEleBDTGV17); 
  FakeEleTree.tree_->SetBranchAddress( "BDTGV18", &fEleBDTGV18); 

  for(UInt_t ientry=0; ientry < FakeEleTree.tree_->GetEntries(); ientry++) {       	
    FakeEleTree.tree_->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //don't evaluate performance using training events
    if (FakeEleTree.fElePt < 5) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(FakeEleTree.fEleEta) < 0.8) subdet = 0;
    else if (fabs(FakeEleTree.fEleEta) < 1.485) subdet = 1;
    else subdet = 2;

    Int_t ptBin = 0;
    if (FakeEleTree.fElePt > 10.0) ptBin = 1;
    if (FakeEleTree.fElePt > 20.0) ptBin = 2;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);    
    if (Option == 6) passCuts = (subdet == 0 && ptBin == 2);
    if (Option == 7) passCuts = (subdet == 1 && ptBin == 2);
    if (Option == 8) passCuts = (subdet == 2 && ptBin == 2);    
    if (Option == 10) passCuts = (ptBin == 0 );
    if (!passCuts) continue;    

//     //apply denominator cuts
//     if (!(FakeEleTree.fEleDenomFake==1 && fabs(FakeEleTree.fEleDZ) < 0.1 && fabs(FakeEleTree.fEleD0) < 0.02)) continue;
//     if (FakeEleTree.fEleMatchedConversion) continue;
    if (!(FakeEleTree.fEleNMissHits <= 1)) continue;
    if( FakeEleTree.fEleIP3dSig > 4 ) continue;

    if( FakeEleTree.fElePFIso04 /FakeEleTree.fElePt  < 1.0 ) continue;

    //Fill Histograms
    Ele_SigmaIEtaIEta_Fake->Fill(FakeEleTree.fEleSigmaIEtaIEta, FakeEleTree.fWeight);
    Ele_Pt_Fake->Fill(FakeEleTree.fElePt, FakeEleTree.fWeight);
    Ele_Eta_Fake->Fill(FakeEleTree.fEleEta, FakeEleTree.fWeight);
    Ele_Rho_Fake->Fill(FakeEleTree.fRho, FakeEleTree.fWeight);
    Ele_D0_Fake->Fill(FakeEleTree.fEleD0, FakeEleTree.fWeight);
    Ele_DZ_Fake->Fill(FakeEleTree.fEleDZ, FakeEleTree.fWeight);
    Ele_IP3D_Fake->Fill(FakeEleTree.fEleIP3d, FakeEleTree.fWeight);
    Ele_NBrem_Fake->Fill(FakeEleTree.fEleNBrem, FakeEleTree.fWeight);
    Ele_FBrem_Fake->Fill(FakeEleTree.fEleFBrem, FakeEleTree.fWeight);
    Ele_EoP_Fake->Fill(FakeEleTree.fEleEOverP, FakeEleTree.fWeight);
    Ele_EEleoPout_Fake->Fill(FakeEleTree.fEleEEleClusterOverPout, FakeEleTree.fWeight);
    Ele_IoEmIoP_Fake->Fill(FakeEleTree.fEleOneOverEMinusOneOverP, FakeEleTree.fWeight);
    Ele_DEta_Fake->Fill(FakeEleTree.fEleDEtaIn, FakeEleTree.fWeight);
    Ele_DPhi_Fake->Fill(FakeEleTree.fEleDPhiIn, FakeEleTree.fWeight);
    Ele_DEtaCalo_Fake->Fill(FakeEleTree.fEledEtaCalo, FakeEleTree.fWeight);
    Ele_DPhiCalo_Fake->Fill(FakeEleTree.fEledPhiCalo, FakeEleTree.fWeight);
    Ele_SigmaIEtaIEta_Fake->Fill(FakeEleTree.fEleSigmaIEtaIEta, FakeEleTree.fWeight);
    Ele_SigmaIPhiIPhi_Fake->Fill(FakeEleTree.fEleSigmaIPhiIPhi, FakeEleTree.fWeight);
    Ele_EtaWidth_Fake->Fill(FakeEleTree.fEleSCEtaWidth, FakeEleTree.fWeight);
    Ele_PhiWidth_Fake->Fill(FakeEleTree.fEleSCPhiWidth, FakeEleTree.fWeight);
    Ele_R9_Fake->Fill(FakeEleTree.fEleR9, FakeEleTree.fWeight);
    Ele_PreShowerOverRaw_Fake->Fill(FakeEleTree.fElePreShowerOverRaw, FakeEleTree.fWeight);
    Ele_HoE_Fake->Fill(FakeEleTree.fEleHoverE, FakeEleTree.fWeight);
    Ele_GsfChi2_Fake->Fill(FakeEleTree.fEleGsfTrackChi2OverNdof, FakeEleTree.fWeight);
    Ele_KFChi2_Fake->Fill(FakeEleTree.fEleKFTrackChi2OverNDoF, FakeEleTree.fWeight);
    Ele_KFLayers_Fake->Fill(FakeEleTree.fEleKFTrackNLayersWithMeasurement, FakeEleTree.fWeight);
    Ele_OneMinusSeedE1x5OverE5x5_Fake->Fill(FakeEleTree.fEleOneMinusSeedE1x5OverE5x5, FakeEleTree.fWeight);
    Ele_PFIso04OverPt_Fake->Fill(FakeEleTree.fElePFIso04 / FakeEleTree.fElePt, FakeEleTree.fWeight);


  } //loop over electrons
  



  
  //*****************************************************************************************
  //MakePlots
  //*****************************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_Pt_Real,RealElectronLabel.c_str(), "L");
  legend->AddEntry(Ele_Pt_Fake,FakeElectronLabel.c_str(), "L");
  Ele_Pt_Real->SetLineColor(kRed);
  Ele_Pt_Fake->SetLineColor(kBlue);
  NormalizeHist(Ele_Pt_Real);
  NormalizeHist(Ele_Pt_Fake);
  Ele_Pt_Real->Draw("hist");
  Ele_Pt_Fake->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_Pt" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_SigmaIEtaIEta_Real,RealElectronLabel.c_str(), "L");
  legend->AddEntry(Ele_SigmaIEtaIEta_Fake,FakeElectronLabel.c_str(), "L");
  Ele_SigmaIEtaIEta_Real->SetLineColor(kRed);
  Ele_SigmaIEtaIEta_Fake->SetLineColor(kBlue);
  NormalizeHist(Ele_SigmaIEtaIEta_Real);
  NormalizeHist(Ele_SigmaIEtaIEta_Fake);
  Ele_SigmaIEtaIEta_Real->Draw("hist");
  Ele_SigmaIEtaIEta_Fake->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_SigmaIEtaIEta" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_PFIso04OverPt_Real,RealElectronLabel.c_str(), "L");
  legend->AddEntry(Ele_PFIso04OverPt_Fake,FakeElectronLabel.c_str(), "L");
  Ele_PFIso04OverPt_Real->SetLineColor(kRed);
  Ele_PFIso04OverPt_Fake->SetLineColor(kBlue);
  NormalizeHist(Ele_PFIso04OverPt_Real);
  NormalizeHist(Ele_PFIso04OverPt_Fake);
  Ele_PFIso04OverPt_Real->Draw("hist");
  Ele_PFIso04OverPt_Fake->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_PFIso04OverPt" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_PhiWidth_Real,RealElectronLabel.c_str(), "L");
  legend->AddEntry(Ele_PhiWidth_Fake,FakeElectronLabel.c_str(), "L");
  Ele_PhiWidth_Real->SetLineColor(kRed);
  Ele_PhiWidth_Fake->SetLineColor(kBlue);
  NormalizeHist(Ele_PhiWidth_Real);
  NormalizeHist(Ele_PhiWidth_Fake);
  Ele_PhiWidth_Real->Draw("hist");
  Ele_PhiWidth_Fake->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_PhiWidth" + label + ".gif").c_str());



  gBenchmark->Show("WWTemplate");       
} 

