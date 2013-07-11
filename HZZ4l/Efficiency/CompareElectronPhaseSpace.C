//root -l HiggsAna/HZZ4l/Efficiency/CompareElectronPhaseSpace.C+'("/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.HZZ125.PtBelow20.root","/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12DY_53X.ZeeTP.PtBelow20.root","/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12DY_53X.ZeeGamma.PtBelow20.root","Bin0",0)'
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
#include "MyStyle.C"
#include "TLegend.h"
#include "TEfficiency.h"

#include "CITCommon/CommonData/interface/ElectronTree.h"

 
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
void CompareElectronPhaseSpace( string signalFile, string zeeFile, string zeeGammaFile,
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
  TH1F *Ele_Pt_Signal = new TH1F(("Ele_Pt_Signal"+label).c_str(), "; Pt ; Number of Events ",  40, 0 , 20);
  TH1F *Ele_Eta_Signal = new TH1F(("Ele_Eta_Signal"+label).c_str(), "; Eta ; Number of Events ",  40, -2.5 , 2.5);
  TH1F *Ele_Rho_Signal = new TH1F(("Ele_Rho_Signal"+label).c_str(), "; Rho ; Number of Events ",  40, 0 , 50);
  TH1F *Ele_D0_Signal = new TH1F(("Ele_D0_Signal"+label).c_str(), "; D0 ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_DZ_Signal = new TH1F(("Ele_DZ_Signal"+label).c_str(), "; DZ ; Number of Events ",  40, 0 , 0.2);
  TH1F *Ele_IP3D_Signal = new TH1F(("Ele_IP3D_Signal"+label).c_str(), "; IP3D ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_NBrem_Signal = new TH1F(("Ele_NBrem_Signal"+label).c_str(), "; NBrem ; Number of Events ",  40, 0 , 10);
  TH1F *Ele_FBrem_Signal = new TH1F(("Ele_FBrem_Signal"+label).c_str(), "; FBrem ; Number of Events ",  20, 0 , 1);
  TH1F *Ele_EoP_Signal = new TH1F(("Ele_EoP_Signal"+label).c_str(), "; EoP ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_EEleoPout_Signal = new TH1F(("Ele_EEleoPout_Signal"+label).c_str(), "; EEleoPout ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_IoEmIoP_Signal = new TH1F(("Ele_IoEmIoP_Signal"+label).c_str(), "; IoEmIoP ; Number of Events ",  40, -0.2 , 0.3);
  TH1F *Ele_DEta_Signal = new TH1F(("Ele_DEta_Signal"+label).c_str(), "; DEta ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_DPhi_Signal = new TH1F(("Ele_DPhi_Signal"+label).c_str(), "; DPhi ; Number of Events ",  40, -1 , 1);
  TH1F *Ele_DEtaCalo_Signal = new TH1F(("Ele_DEtaCalo_Signal"+label).c_str(), "; DEtaCalo ; Number of Events ",  40, -0.2 , 0.2);
  TH1F *Ele_DPhiCalo_Signal = new TH1F(("Ele_DPhiCalo_Signal"+label).c_str(), "; DPhiCalo ; Number of Events ",  40, -1 , 1);
  TH1F *Ele_SigmaIEtaIEta_Signal = new TH1F(("Ele_SigmaIEtaIEta_Signal"+label).c_str(), "; SigmaIEtaIEta ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_SigmaIPhiIPhi_Signal = new TH1F(("Ele_SigmaIPhiIPhi_Signal"+label).c_str(), "; SigmaIPhiIPhi ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_EtaWidth_Signal = new TH1F(("Ele_EtaWidth_Signal"+label).c_str(), "; EtaWidth ; Number of Events ",  40, 0 , 0.5);
  TH1F *Ele_PhiWidth_Signal = new TH1F(("Ele_PhiWidth_Signal"+label).c_str(), "; PhiWidth ; Number of Events ",  40, 0 , 0.3);
  TH1F *Ele_R9_Signal = new TH1F(("Ele_R9_Signal"+label).c_str(), "; R9 ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_PreShowerOverRaw_Signal = new TH1F(("Ele_PreShowerOverRaw_Signal"+label).c_str(), "; PreShowerOverRaw ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_HoE_Signal = new TH1F(("Ele_HoE_Signal"+label).c_str(), "; HoE ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_GsfChi2_Signal = new TH1F(("Ele_GsfChi2_Signal"+label).c_str(), "; GsfChi2 ; Number of Events ",  40, 0 , 200);
  TH1F *Ele_KFChi2_Signal = new TH1F(("Ele_KFChi2_Signal"+label).c_str(), "; KFChi2 ; Number of Events ",  40, 0 , 10);
  TH1F *Ele_KFLayers_Signal = new TH1F(("Ele_KFLayers_Signal"+label).c_str(), "; KFLayers ; Number of Events ",  21, -0.5 , 20.5);
  TH1F *Ele_OneMinusSeedE1x5OverE5x5_Signal = new TH1F(("Ele_OneMinusSeedE1x5OverE5x5_Signal"+label).c_str(), "; OneMinusSeedE1x5OverE5x5 ; Number of Events ",  40, -1.5 , 1.5);
  TH1F *Ele_PFIso04OverPt_Signal = new TH1F(("Ele_PFIso04OverPt_Signal"+label).c_str(), "; PFIso04OverPt ; Number of Events ",  40, 0 , 5);
  TH1F *Ele_EGammaNonTriggeringMVA_Signal = new TH1F(("Ele_EGammaNonTriggeringMVA_Signal"+label).c_str(), "; EGammaNonTriggeringMVA ; Number of Events ",  40, -1 , 1);

  TH1F *Ele_Pt_Zee = new TH1F(("Ele_Pt_Zee"+label).c_str(), "; Pt ; Number of Events ",  40, 0 , 20);
  TH1F *Ele_Eta_Zee = new TH1F(("Ele_Eta_Zee"+label).c_str(), "; Eta ; Number of Events ",  40, -2.5 , 2.5);
  TH1F *Ele_Rho_Zee = new TH1F(("Ele_Rho_Zee"+label).c_str(), "; Rho ; Number of Events ",  40, 0 , 50);
  TH1F *Ele_D0_Zee = new TH1F(("Ele_D0_Zee"+label).c_str(), "; D0 ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_DZ_Zee = new TH1F(("Ele_DZ_Zee"+label).c_str(), "; DZ ; Number of Events ",  40, 0 , 0.2);
  TH1F *Ele_IP3D_Zee = new TH1F(("Ele_IP3D_Zee"+label).c_str(), "; IP3D ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_NBrem_Zee = new TH1F(("Ele_NBrem_Zee"+label).c_str(), "; NBrem ; Number of Events ",  40, 0 , 10);
  TH1F *Ele_FBrem_Zee = new TH1F(("Ele_FBrem_Zee"+label).c_str(), "; FBrem ; Number of Events ",  20, 0 , 1);
  TH1F *Ele_EoP_Zee = new TH1F(("Ele_EoP_Zee"+label).c_str(), "; EoP ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_EEleoPout_Zee = new TH1F(("Ele_EEleoPout_Zee"+label).c_str(), "; EEleoPout ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_IoEmIoP_Zee = new TH1F(("Ele_IoEmIoP_Zee"+label).c_str(), "; IoEmIoP ; Number of Events ",  40, -0.2 , 0.3);
  TH1F *Ele_DEta_Zee = new TH1F(("Ele_DEta_Zee"+label).c_str(), "; DEta ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_DPhi_Zee = new TH1F(("Ele_DPhi_Zee"+label).c_str(), "; DPhi ; Number of Events ",  40, -1 , 1);
  TH1F *Ele_DEtaCalo_Zee = new TH1F(("Ele_DEtaCalo_Zee"+label).c_str(), "; DEtaCalo ; Number of Events ",  40, -0.2 , 0.2);
  TH1F *Ele_DPhiCalo_Zee = new TH1F(("Ele_DPhiCalo_Zee"+label).c_str(), "; DPhiCalo ; Number of Events ",  40, -1 , 1);
  TH1F *Ele_SigmaIEtaIEta_Zee = new TH1F(("Ele_SigmaIEtaIEta_Zee"+label).c_str(), "; SigmaIEtaIEta ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_SigmaIPhiIPhi_Zee = new TH1F(("Ele_SigmaIPhiIPhi_Zee"+label).c_str(), "; SigmaIPhiIPhi ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_EtaWidth_Zee = new TH1F(("Ele_EtaWidth_Zee"+label).c_str(), "; EtaWidth ; Number of Events ",  40, 0 , 0.5);
  TH1F *Ele_PhiWidth_Zee = new TH1F(("Ele_PhiWidth_Zee"+label).c_str(), "; PhiWidth ; Number of Events ",  40, 0 , 0.3);
  TH1F *Ele_R9_Zee = new TH1F(("Ele_R9_Zee"+label).c_str(), "; R9 ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_PreShowerOverRaw_Zee = new TH1F(("Ele_PreShowerOverRaw_Zee"+label).c_str(), "; PreShowerOverRaw ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_HoE_Zee = new TH1F(("Ele_HoE_Zee"+label).c_str(), "; HoE ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_GsfChi2_Zee = new TH1F(("Ele_GsfChi2_Zee"+label).c_str(), "; GsfChi2 ; Number of Events ",  40, 0 , 200);
  TH1F *Ele_KFChi2_Zee = new TH1F(("Ele_KFChi2_Zee"+label).c_str(), "; KFChi2 ; Number of Events ",  40, 0 , 10);
  TH1F *Ele_KFLayers_Zee = new TH1F(("Ele_KFLayers_Zee"+label).c_str(), "; KFLayers ; Number of Events ",  21, -0.5 , 20.5);
  TH1F *Ele_OneMinusSeedE1x5OverE5x5_Zee = new TH1F(("Ele_OneMinusSeedE1x5OverE5x5_Zee"+label).c_str(), "; OneMinusSeedE1x5OverE5x5 ; Number of Events ",  40, -1.5 , 1.5);
  TH1F *Ele_PFIso04OverPt_Zee = new TH1F(("Ele_PFIso04OverPt_Zee"+label).c_str(), "; PFIso04OverPt ; Number of Events ",  40, 0 , 5);
  TH1F *Ele_EGammaNonTriggeringMVA_Zee = new TH1F(("Ele_EGammaNonTriggeringMVA_Zee"+label).c_str(), "; EGammaNonTriggeringMVA ; Number of Events ",  40, -1 , 1);

  TH1F *Ele_Pt_ZeeGamma = new TH1F(("Ele_Pt_ZeeGamma"+label).c_str(), "; Pt ; Number of Events ",  40, 0 , 20);
  TH1F *Ele_Eta_ZeeGamma = new TH1F(("Ele_Eta_ZeeGamma"+label).c_str(), "; Eta ; Number of Events ",  40, -2.5 , 2.5);
  TH1F *Ele_Rho_ZeeGamma = new TH1F(("Ele_Rho_ZeeGamma"+label).c_str(), "; Rho ; Number of Events ",  40, 0 , 50);
  TH1F *Ele_D0_ZeeGamma = new TH1F(("Ele_D0_ZeeGamma"+label).c_str(), "; D0 ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_DZ_ZeeGamma = new TH1F(("Ele_DZ_ZeeGamma"+label).c_str(), "; DZ ; Number of Events ",  40, 0 , 0.2);
  TH1F *Ele_IP3D_ZeeGamma = new TH1F(("Ele_IP3D_ZeeGamma"+label).c_str(), "; IP3D ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_NBrem_ZeeGamma = new TH1F(("Ele_NBrem_ZeeGamma"+label).c_str(), "; NBrem ; Number of Events ",  40, 0 , 10);
  TH1F *Ele_FBrem_ZeeGamma = new TH1F(("Ele_FBrem_ZeeGamma"+label).c_str(), "; FBrem ; Number of Events ",  20, 0 , 1);
  TH1F *Ele_EoP_ZeeGamma = new TH1F(("Ele_EoP_ZeeGamma"+label).c_str(), "; EoP ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_EEleoPout_ZeeGamma = new TH1F(("Ele_EEleoPout_ZeeGamma"+label).c_str(), "; EEleoPout ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_IoEmIoP_ZeeGamma = new TH1F(("Ele_IoEmIoP_ZeeGamma"+label).c_str(), "; IoEmIoP ; Number of Events ",  40, -0.2 , 0.3);
  TH1F *Ele_DEta_ZeeGamma = new TH1F(("Ele_DEta_ZeeGamma"+label).c_str(), "; DEta ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_DPhi_ZeeGamma = new TH1F(("Ele_DPhi_ZeeGamma"+label).c_str(), "; DPhi ; Number of Events ",  40, -1 , 1);
  TH1F *Ele_DEtaCalo_ZeeGamma = new TH1F(("Ele_DEtaCalo_ZeeGamma"+label).c_str(), "; DEtaCalo ; Number of Events ",  40, -0.2 , 0.2);
  TH1F *Ele_DPhiCalo_ZeeGamma = new TH1F(("Ele_DPhiCalo_ZeeGamma"+label).c_str(), "; DPhiCalo ; Number of Events ",  40, -1 , 1);
  TH1F *Ele_SigmaIEtaIEta_ZeeGamma = new TH1F(("Ele_SigmaIEtaIEta_ZeeGamma"+label).c_str(), "; SigmaIEtaIEta ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_SigmaIPhiIPhi_ZeeGamma = new TH1F(("Ele_SigmaIPhiIPhi_ZeeGamma"+label).c_str(), "; SigmaIPhiIPhi ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_EtaWidth_ZeeGamma = new TH1F(("Ele_EtaWidth_ZeeGamma"+label).c_str(), "; EtaWidth ; Number of Events ",  40, 0 , 0.5);
  TH1F *Ele_PhiWidth_ZeeGamma = new TH1F(("Ele_PhiWidth_ZeeGamma"+label).c_str(), "; PhiWidth ; Number of Events ",  40, 0 , 0.3);
  TH1F *Ele_R9_ZeeGamma = new TH1F(("Ele_R9_ZeeGamma"+label).c_str(), "; R9 ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_PreShowerOverRaw_ZeeGamma = new TH1F(("Ele_PreShowerOverRaw_ZeeGamma"+label).c_str(), "; PreShowerOverRaw ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_HoE_ZeeGamma = new TH1F(("Ele_HoE_ZeeGamma"+label).c_str(), "; HoE ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_GsfChi2_ZeeGamma = new TH1F(("Ele_GsfChi2_ZeeGamma"+label).c_str(), "; GsfChi2 ; Number of Events ",  40, 0 , 200);
  TH1F *Ele_KFChi2_ZeeGamma = new TH1F(("Ele_KFChi2_ZeeGamma"+label).c_str(), "; KFChi2 ; Number of Events ",  40, 0 , 10);
  TH1F *Ele_KFLayers_ZeeGamma = new TH1F(("Ele_KFLayers_ZeeGamma"+label).c_str(), "; KFLayers ; Number of Events ",  21, -0.5 , 20.5);
  TH1F *Ele_OneMinusSeedE1x5OverE5x5_ZeeGamma = new TH1F(("Ele_OneMinusSeedE1x5OverE5x5_ZeeGamma"+label).c_str(), "; OneMinusSeedE1x5OverE5x5 ; Number of Events ",  40, -1.5 , 1.5);
  TH1F *Ele_PFIso04OverPt_ZeeGamma = new TH1F(("Ele_PFIso04OverPt_ZeeGamma"+label).c_str(), "; PFIso04OverPt ; Number of Events ",  40, 0 , 5);
  TH1F *Ele_EGammaNonTriggeringMVA_ZeeGamma = new TH1F(("Ele_EGammaNonTriggeringMVA_ZeeGamma"+label).c_str(), "; EGammaNonTriggeringMVA ; Number of Events ",  40, -1 , 1);

  //*****************************************************************************************
  //SignalEleTree
  //*****************************************************************************************
  citana::ElectronTree SignalEleTree;
  SignalEleTree.LoadTree(signalFile.c_str());
  SignalEleTree.InitTree();

  for(UInt_t ientry=0; ientry < SignalEleTree.tree_->GetEntries(); ientry++) {       	
    SignalEleTree.tree_->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //don't evaluate performance using training events
    if (SignalEleTree.fElePt < 7) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(SignalEleTree.fEleEta) < 1.485) subdet = 0;
    else subdet = 1;

    Int_t ptBin = 0;
    if (SignalEleTree.fElePt > 15.0) ptBin = 1;
    if (SignalEleTree.fElePt > 20.0) ptBin = 2;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 3) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 0 && ptBin == 2);
    if (Option == 5) passCuts = (subdet == 1 && ptBin == 2);
    if (Option == 10) passCuts = (ptBin == 0 );

    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

    if (!(SignalEleTree.fEleNMissHits <= 1)) continue;
    if( SignalEleTree.fEleIP3dSig > 4 ) continue;

    
    //Fill Histograms
    Ele_SigmaIEtaIEta_Signal->Fill(SignalEleTree.fEleSigmaIEtaIEta, SignalEleTree.fWeight);
    Ele_Pt_Signal->Fill(SignalEleTree.fElePt, SignalEleTree.fWeight);
    Ele_Eta_Signal->Fill(SignalEleTree.fEleEta, SignalEleTree.fWeight);
    Ele_Rho_Signal->Fill(SignalEleTree.fRho, SignalEleTree.fWeight);
    Ele_D0_Signal->Fill(SignalEleTree.fEleD0, SignalEleTree.fWeight);
    Ele_DZ_Signal->Fill(SignalEleTree.fEleDZ, SignalEleTree.fWeight);
    Ele_IP3D_Signal->Fill(SignalEleTree.fEleIP3d, SignalEleTree.fWeight);
    Ele_NBrem_Signal->Fill(SignalEleTree.fEleNBrem, SignalEleTree.fWeight);
    Ele_FBrem_Signal->Fill(SignalEleTree.fEleFBrem, SignalEleTree.fWeight);
    Ele_EoP_Signal->Fill(SignalEleTree.fEleEOverP, SignalEleTree.fWeight);
    Ele_EEleoPout_Signal->Fill(SignalEleTree.fEleEEleClusterOverPout, SignalEleTree.fWeight);
    Ele_IoEmIoP_Signal->Fill(SignalEleTree.fEleOneOverEMinusOneOverP, SignalEleTree.fWeight);
    Ele_DEta_Signal->Fill(SignalEleTree.fEleDEtaIn, SignalEleTree.fWeight);
    Ele_DPhi_Signal->Fill(SignalEleTree.fEleDPhiIn, SignalEleTree.fWeight);
    Ele_DEtaCalo_Signal->Fill(SignalEleTree.fEledEtaCalo, SignalEleTree.fWeight);
    Ele_DPhiCalo_Signal->Fill(SignalEleTree.fEledPhiCalo, SignalEleTree.fWeight);
    Ele_SigmaIEtaIEta_Signal->Fill(SignalEleTree.fEleSigmaIEtaIEta, SignalEleTree.fWeight);
    Ele_SigmaIPhiIPhi_Signal->Fill(SignalEleTree.fEleSigmaIPhiIPhi, SignalEleTree.fWeight);
    Ele_EtaWidth_Signal->Fill(SignalEleTree.fEleSCEtaWidth, SignalEleTree.fWeight);
    Ele_PhiWidth_Signal->Fill(SignalEleTree.fEleSCPhiWidth, SignalEleTree.fWeight);
    Ele_R9_Signal->Fill(SignalEleTree.fEleR9, SignalEleTree.fWeight);
    Ele_PreShowerOverRaw_Signal->Fill(SignalEleTree.fElePreShowerOverRaw, SignalEleTree.fWeight);
    Ele_HoE_Signal->Fill(SignalEleTree.fEleHoverE, SignalEleTree.fWeight);
    Ele_GsfChi2_Signal->Fill(SignalEleTree.fEleGsfTrackChi2OverNdof, SignalEleTree.fWeight);
    Ele_KFChi2_Signal->Fill(SignalEleTree.fEleKFTrackChi2OverNDoF, SignalEleTree.fWeight);
    Ele_KFLayers_Signal->Fill(SignalEleTree.fEleKFTrackNLayersWithMeasurement, SignalEleTree.fWeight);
    Ele_OneMinusSeedE1x5OverE5x5_Signal->Fill(SignalEleTree.fEleOneMinusSeedE1x5OverE5x5, SignalEleTree.fWeight);
    Ele_PFIso04OverPt_Signal->Fill(SignalEleTree.fElePFIso04 / SignalEleTree.fElePt, SignalEleTree.fWeight);
    Ele_EGammaNonTriggeringMVA_Signal->Fill(SignalEleTree.fEGammaNonTriggeringMVA , SignalEleTree.fWeight);


  } 
  




  //*****************************************************************************************
  //ZeeEleTree
  //*****************************************************************************************
  citana::ElectronTree ZeeEleTree;
  ZeeEleTree.LoadTree(zeeFile.c_str());
  ZeeEleTree.InitTree();

  for(UInt_t ientry=0; ientry < ZeeEleTree.tree_->GetEntries(); ientry++) {       	
    ZeeEleTree.tree_->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //don't evaluate performance using training events
    if (ZeeEleTree.fElePt < 7) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(ZeeEleTree.fEleEta) < 1.485) subdet = 0;
    else subdet = 1;

    Int_t ptBin = 0;
    if (ZeeEleTree.fElePt > 15.0) ptBin = 1;
    if (ZeeEleTree.fElePt > 20.0) ptBin = 2;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 3) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 0 && ptBin == 2);
    if (Option == 5) passCuts = (subdet == 1 && ptBin == 2);
    if (Option == 10) passCuts = (ptBin == 0 );

    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

    if (!(ZeeEleTree.fEleNMissHits <= 1)) continue;
    if( ZeeEleTree.fEleIP3dSig > 4 ) continue;

    
    //Fill Histograms
    Ele_SigmaIEtaIEta_Zee->Fill(ZeeEleTree.fEleSigmaIEtaIEta, ZeeEleTree.fWeight);
    Ele_Pt_Zee->Fill(ZeeEleTree.fElePt, ZeeEleTree.fWeight);
    Ele_Eta_Zee->Fill(ZeeEleTree.fEleEta, ZeeEleTree.fWeight);
    Ele_Rho_Zee->Fill(ZeeEleTree.fRho, ZeeEleTree.fWeight);
    Ele_D0_Zee->Fill(ZeeEleTree.fEleD0, ZeeEleTree.fWeight);
    Ele_DZ_Zee->Fill(ZeeEleTree.fEleDZ, ZeeEleTree.fWeight);
    Ele_IP3D_Zee->Fill(ZeeEleTree.fEleIP3d, ZeeEleTree.fWeight);
    Ele_NBrem_Zee->Fill(ZeeEleTree.fEleNBrem, ZeeEleTree.fWeight);
    Ele_FBrem_Zee->Fill(ZeeEleTree.fEleFBrem, ZeeEleTree.fWeight);
    Ele_EoP_Zee->Fill(ZeeEleTree.fEleEOverP, ZeeEleTree.fWeight);
    Ele_EEleoPout_Zee->Fill(ZeeEleTree.fEleEEleClusterOverPout, ZeeEleTree.fWeight);
    Ele_IoEmIoP_Zee->Fill(ZeeEleTree.fEleOneOverEMinusOneOverP, ZeeEleTree.fWeight);
    Ele_DEta_Zee->Fill(ZeeEleTree.fEleDEtaIn, ZeeEleTree.fWeight);
    Ele_DPhi_Zee->Fill(ZeeEleTree.fEleDPhiIn, ZeeEleTree.fWeight);
    Ele_DEtaCalo_Zee->Fill(ZeeEleTree.fEledEtaCalo, ZeeEleTree.fWeight);
    Ele_DPhiCalo_Zee->Fill(ZeeEleTree.fEledPhiCalo, ZeeEleTree.fWeight);
    Ele_SigmaIEtaIEta_Zee->Fill(ZeeEleTree.fEleSigmaIEtaIEta, ZeeEleTree.fWeight);
    Ele_SigmaIPhiIPhi_Zee->Fill(ZeeEleTree.fEleSigmaIPhiIPhi, ZeeEleTree.fWeight);
    Ele_EtaWidth_Zee->Fill(ZeeEleTree.fEleSCEtaWidth, ZeeEleTree.fWeight);
    Ele_PhiWidth_Zee->Fill(ZeeEleTree.fEleSCPhiWidth, ZeeEleTree.fWeight);
    Ele_R9_Zee->Fill(ZeeEleTree.fEleR9, ZeeEleTree.fWeight);
    Ele_PreShowerOverRaw_Zee->Fill(ZeeEleTree.fElePreShowerOverRaw, ZeeEleTree.fWeight);
    Ele_HoE_Zee->Fill(ZeeEleTree.fEleHoverE, ZeeEleTree.fWeight);
    Ele_GsfChi2_Zee->Fill(ZeeEleTree.fEleGsfTrackChi2OverNdof, ZeeEleTree.fWeight);
    Ele_KFChi2_Zee->Fill(ZeeEleTree.fEleKFTrackChi2OverNDoF, ZeeEleTree.fWeight);
    Ele_KFLayers_Zee->Fill(ZeeEleTree.fEleKFTrackNLayersWithMeasurement, ZeeEleTree.fWeight);
    Ele_OneMinusSeedE1x5OverE5x5_Zee->Fill(ZeeEleTree.fEleOneMinusSeedE1x5OverE5x5, ZeeEleTree.fWeight);
    Ele_PFIso04OverPt_Zee->Fill(ZeeEleTree.fElePFIso04 / ZeeEleTree.fElePt, ZeeEleTree.fWeight);
    Ele_EGammaNonTriggeringMVA_Zee->Fill(ZeeEleTree.fEGammaNonTriggeringMVA , ZeeEleTree.fWeight);


  } 
  


  //*****************************************************************************************
  //ZeeGammaEleTree
  //*****************************************************************************************
  citana::ElectronTree ZeeGammaEleTree;
  ZeeGammaEleTree.LoadTree(zeeGammaFile.c_str());
  ZeeGammaEleTree.InitTree();

  for(UInt_t ientry=0; ientry < ZeeGammaEleTree.tree_->GetEntries(); ientry++) {       	
    ZeeGammaEleTree.tree_->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //don't evaluate performance using training events
    if (ZeeGammaEleTree.fElePt < 7) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(ZeeGammaEleTree.fEleEta) < 1.485) subdet = 0;
    else subdet = 1;

    Int_t ptBin = 0;
    if (ZeeGammaEleTree.fElePt > 15.0) ptBin = 1;
    if (ZeeGammaEleTree.fElePt > 20.0) ptBin = 2;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 3) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 0 && ptBin == 2);
    if (Option == 5) passCuts = (subdet == 1 && ptBin == 2);
    if (Option == 10) passCuts = (ptBin == 0 );

    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

    if (!(ZeeGammaEleTree.fEleNMissHits <= 1)) continue;
    if( ZeeGammaEleTree.fEleIP3dSig > 4 ) continue;

    
    //Fill Histograms
    Ele_SigmaIEtaIEta_ZeeGamma->Fill(ZeeGammaEleTree.fEleSigmaIEtaIEta, ZeeGammaEleTree.fWeight);
    Ele_Pt_ZeeGamma->Fill(ZeeGammaEleTree.fElePt, ZeeGammaEleTree.fWeight);
    Ele_Eta_ZeeGamma->Fill(ZeeGammaEleTree.fEleEta, ZeeGammaEleTree.fWeight);
    Ele_Rho_ZeeGamma->Fill(ZeeGammaEleTree.fRho, ZeeGammaEleTree.fWeight);
    Ele_D0_ZeeGamma->Fill(ZeeGammaEleTree.fEleD0, ZeeGammaEleTree.fWeight);
    Ele_DZ_ZeeGamma->Fill(ZeeGammaEleTree.fEleDZ, ZeeGammaEleTree.fWeight);
    Ele_IP3D_ZeeGamma->Fill(ZeeGammaEleTree.fEleIP3d, ZeeGammaEleTree.fWeight);
    Ele_NBrem_ZeeGamma->Fill(ZeeGammaEleTree.fEleNBrem, ZeeGammaEleTree.fWeight);
    Ele_FBrem_ZeeGamma->Fill(ZeeGammaEleTree.fEleFBrem, ZeeGammaEleTree.fWeight);
    Ele_EoP_ZeeGamma->Fill(ZeeGammaEleTree.fEleEOverP, ZeeGammaEleTree.fWeight);
    Ele_EEleoPout_ZeeGamma->Fill(ZeeGammaEleTree.fEleEEleClusterOverPout, ZeeGammaEleTree.fWeight);
    Ele_IoEmIoP_ZeeGamma->Fill(ZeeGammaEleTree.fEleOneOverEMinusOneOverP, ZeeGammaEleTree.fWeight);
    Ele_DEta_ZeeGamma->Fill(ZeeGammaEleTree.fEleDEtaIn, ZeeGammaEleTree.fWeight);
    Ele_DPhi_ZeeGamma->Fill(ZeeGammaEleTree.fEleDPhiIn, ZeeGammaEleTree.fWeight);
    Ele_DEtaCalo_ZeeGamma->Fill(ZeeGammaEleTree.fEledEtaCalo, ZeeGammaEleTree.fWeight);
    Ele_DPhiCalo_ZeeGamma->Fill(ZeeGammaEleTree.fEledPhiCalo, ZeeGammaEleTree.fWeight);
    Ele_SigmaIEtaIEta_ZeeGamma->Fill(ZeeGammaEleTree.fEleSigmaIEtaIEta, ZeeGammaEleTree.fWeight);
    Ele_SigmaIPhiIPhi_ZeeGamma->Fill(ZeeGammaEleTree.fEleSigmaIPhiIPhi, ZeeGammaEleTree.fWeight);
    Ele_EtaWidth_ZeeGamma->Fill(ZeeGammaEleTree.fEleSCEtaWidth, ZeeGammaEleTree.fWeight);
    Ele_PhiWidth_ZeeGamma->Fill(ZeeGammaEleTree.fEleSCPhiWidth, ZeeGammaEleTree.fWeight);
    Ele_R9_ZeeGamma->Fill(ZeeGammaEleTree.fEleR9, ZeeGammaEleTree.fWeight);
    Ele_PreShowerOverRaw_ZeeGamma->Fill(ZeeGammaEleTree.fElePreShowerOverRaw, ZeeGammaEleTree.fWeight);
    Ele_HoE_ZeeGamma->Fill(ZeeGammaEleTree.fEleHoverE, ZeeGammaEleTree.fWeight);
    Ele_GsfChi2_ZeeGamma->Fill(ZeeGammaEleTree.fEleGsfTrackChi2OverNdof, ZeeGammaEleTree.fWeight);
    Ele_KFChi2_ZeeGamma->Fill(ZeeGammaEleTree.fEleKFTrackChi2OverNDoF, ZeeGammaEleTree.fWeight);
    Ele_KFLayers_ZeeGamma->Fill(ZeeGammaEleTree.fEleKFTrackNLayersWithMeasurement, ZeeGammaEleTree.fWeight);
    Ele_OneMinusSeedE1x5OverE5x5_ZeeGamma->Fill(ZeeGammaEleTree.fEleOneMinusSeedE1x5OverE5x5, ZeeGammaEleTree.fWeight);
    Ele_PFIso04OverPt_ZeeGamma->Fill(ZeeGammaEleTree.fElePFIso04 / ZeeGammaEleTree.fElePt, ZeeGammaEleTree.fWeight);
    Ele_EGammaNonTriggeringMVA_ZeeGamma->Fill(ZeeGammaEleTree.fEGammaNonTriggeringMVA , ZeeGammaEleTree.fWeight);


  } 
  





  
  //*****************************************************************************************
  //MakePlots
  //*****************************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;
//   gStyle->SetStatOpt(0);


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.60,0.50,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_Pt_Signal, "HZZ Signal", "L");
  legend->AddEntry(Ele_Pt_Zee, "Zee", "L");
  legend->AddEntry(Ele_Pt_ZeeGamma, "ZeeGamma", "L");
  Ele_Pt_Signal->SetLineColor(kRed);
  Ele_Pt_Zee->SetLineColor(kBlue);
  Ele_Pt_ZeeGamma->SetLineColor(kMagenta);
  NormalizeHist(Ele_Pt_Signal);
  NormalizeHist(Ele_Pt_Zee);
  NormalizeHist(Ele_Pt_ZeeGamma);
  Ele_Pt_Signal->SetMaximum( 1.2*max(max(Ele_Pt_Zee->GetMaximum(),Ele_Pt_ZeeGamma->GetMaximum()),Ele_Pt_Signal->GetMaximum()) );
  Ele_Pt_Signal->Draw("hist");
  Ele_Pt_Zee->Draw("hist,same");
  Ele_Pt_ZeeGamma->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_Pt" + label + ".gif").c_str());


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_Eta_Signal, "HZZ Signal", "L");
  legend->AddEntry(Ele_Eta_Zee, "Zee", "L");
  legend->AddEntry(Ele_Eta_ZeeGamma, "ZeeGamma", "L");
  Ele_Eta_Signal->SetLineColor(kRed);
  Ele_Eta_Zee->SetLineColor(kBlue);
  Ele_Eta_ZeeGamma->SetLineColor(kMagenta);
  NormalizeHist(Ele_Eta_Signal);
  NormalizeHist(Ele_Eta_Zee);
  NormalizeHist(Ele_Eta_ZeeGamma);
  Ele_Eta_Signal->SetMaximum( 1.2*max(max(Ele_Eta_Zee->GetMaximum(),Ele_Eta_ZeeGamma->GetMaximum()),Ele_Eta_Signal->GetMaximum()) );
  Ele_Eta_Signal->Draw("hist");
  Ele_Eta_Zee->Draw("hist,same");
  Ele_Eta_ZeeGamma->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_Eta" + label + ".gif").c_str());
  return;

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.70,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_SigmaIEtaIEta_Signal, "HZZ Signal", "L");
  legend->AddEntry(Ele_SigmaIEtaIEta_Zee, "Zee", "L");
  legend->AddEntry(Ele_SigmaIEtaIEta_ZeeGamma, "ZeeGamma", "L");
  Ele_SigmaIEtaIEta_Signal->SetLineColor(kRed);
  Ele_SigmaIEtaIEta_Zee->SetLineColor(kBlue);
  Ele_SigmaIEtaIEta_ZeeGamma->SetLineColor(kMagenta);
  NormalizeHist(Ele_SigmaIEtaIEta_Signal);
  NormalizeHist(Ele_SigmaIEtaIEta_Zee);
  NormalizeHist(Ele_SigmaIEtaIEta_ZeeGamma);
  Ele_SigmaIEtaIEta_Signal->SetMaximum( 1.2*max(max(Ele_SigmaIEtaIEta_Zee->GetMaximum(),Ele_SigmaIEtaIEta_ZeeGamma->GetMaximum()),Ele_SigmaIEtaIEta_Signal->GetMaximum()) );
  Ele_SigmaIEtaIEta_Signal->GetXaxis()->SetRangeUser(0,0.04);
  Ele_SigmaIEtaIEta_Signal->Draw("hist");
  Ele_SigmaIEtaIEta_Zee->Draw("hist,same");
  Ele_SigmaIEtaIEta_ZeeGamma->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_SigmaIEtaIEta" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.60,0.94,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_FBrem_Signal, "HZZ Signal", "L");
  legend->AddEntry(Ele_FBrem_Zee, "Zee", "L");
  legend->AddEntry(Ele_FBrem_ZeeGamma, "ZeeGamma", "L");
  Ele_FBrem_Signal->SetLineColor(kRed);
  Ele_FBrem_Zee->SetLineColor(kBlue);
  Ele_FBrem_ZeeGamma->SetLineColor(kMagenta);
  NormalizeHist(Ele_FBrem_Signal);
  NormalizeHist(Ele_FBrem_Zee);
  NormalizeHist(Ele_FBrem_ZeeGamma);
  Ele_FBrem_Signal->SetMaximum( 1.2*max(max(Ele_FBrem_Zee->GetMaximum(),Ele_FBrem_ZeeGamma->GetMaximum()),Ele_FBrem_Signal->GetMaximum()) );
  Ele_FBrem_Signal->Draw("hist");
  Ele_FBrem_Zee->Draw("hist,same");
  Ele_FBrem_ZeeGamma->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_FBrem" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.50,0.50,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_EoP_Signal, "HZZ Signal", "L");
  legend->AddEntry(Ele_EoP_Zee, "Zee", "L");
  legend->AddEntry(Ele_EoP_ZeeGamma, "ZeeGamma", "L");
  Ele_EoP_Signal->SetLineColor(kRed);
  Ele_EoP_Zee->SetLineColor(kBlue);
  Ele_EoP_ZeeGamma->SetLineColor(kMagenta);
  NormalizeHist(Ele_EoP_Signal);
  NormalizeHist(Ele_EoP_Zee);
  NormalizeHist(Ele_EoP_ZeeGamma);
  Ele_EoP_Signal->SetMaximum( 1.2*max(max(Ele_EoP_Zee->GetMaximum(),Ele_EoP_ZeeGamma->GetMaximum()),Ele_EoP_Signal->GetMaximum()) );
  Ele_EoP_Signal->Draw("hist");
  Ele_EoP_Zee->Draw("hist,same");
  Ele_EoP_ZeeGamma->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_EoP" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.70,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_R9_Signal, "HZZ Signal", "L");
  legend->AddEntry(Ele_R9_Zee, "Zee", "L");
  legend->AddEntry(Ele_R9_ZeeGamma, "ZeeGamma", "L");
  Ele_R9_Signal->SetLineColor(kRed);
  Ele_R9_Zee->SetLineColor(kBlue);
  Ele_R9_ZeeGamma->SetLineColor(kMagenta);
  NormalizeHist(Ele_R9_Signal);
  NormalizeHist(Ele_R9_Zee);
  NormalizeHist(Ele_R9_ZeeGamma);
  Ele_R9_Signal->SetMaximum( 1.2*max(max(Ele_R9_Zee->GetMaximum(),Ele_R9_ZeeGamma->GetMaximum()),Ele_R9_Signal->GetMaximum()) );
  Ele_R9_Signal->Draw("hist");
  Ele_R9_Zee->Draw("hist,same");
  Ele_R9_ZeeGamma->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_R9" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_PhiWidth_Signal, "HZZ Signal", "L");
  legend->AddEntry(Ele_PhiWidth_Zee, "Zee", "L");
  legend->AddEntry(Ele_PhiWidth_ZeeGamma, "ZeeGamma", "L");
  Ele_PhiWidth_Signal->SetLineColor(kRed);
  Ele_PhiWidth_Zee->SetLineColor(kBlue);
  Ele_PhiWidth_ZeeGamma->SetLineColor(kMagenta);
  NormalizeHist(Ele_PhiWidth_Signal);
  NormalizeHist(Ele_PhiWidth_Zee);
  NormalizeHist(Ele_PhiWidth_ZeeGamma);
  Ele_PhiWidth_Signal->SetMaximum( 1.2*max(max(Ele_PhiWidth_Zee->GetMaximum(),Ele_PhiWidth_ZeeGamma->GetMaximum()),Ele_PhiWidth_Signal->GetMaximum()) );
  Ele_PhiWidth_Signal->Draw("hist");
  Ele_PhiWidth_Zee->Draw("hist,same");
  Ele_PhiWidth_ZeeGamma->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_PhiWidth" + label + ".gif").c_str());


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.50,0.50,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_EGammaNonTriggeringMVA_Signal, "HZZ Signal", "L");
  legend->AddEntry(Ele_EGammaNonTriggeringMVA_Zee, "Zee", "L");
  legend->AddEntry(Ele_EGammaNonTriggeringMVA_ZeeGamma, "ZeeGamma", "L");
  Ele_EGammaNonTriggeringMVA_Signal->SetLineColor(kRed);
  Ele_EGammaNonTriggeringMVA_Zee->SetLineColor(kBlue);
  Ele_EGammaNonTriggeringMVA_ZeeGamma->SetLineColor(kMagenta);
  NormalizeHist(Ele_EGammaNonTriggeringMVA_Signal);
  NormalizeHist(Ele_EGammaNonTriggeringMVA_Zee);
  NormalizeHist(Ele_EGammaNonTriggeringMVA_ZeeGamma);
  Ele_EGammaNonTriggeringMVA_Signal->SetMaximum( 1.2*max(max(Ele_EGammaNonTriggeringMVA_Zee->GetMaximum(),Ele_EGammaNonTriggeringMVA_ZeeGamma->GetMaximum()),Ele_EGammaNonTriggeringMVA_Signal->GetMaximum()) );
  Ele_EGammaNonTriggeringMVA_Signal->Draw("hist");
  Ele_EGammaNonTriggeringMVA_Zee->Draw("hist,same");
  Ele_EGammaNonTriggeringMVA_ZeeGamma->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_EGammaNonTriggeringMVA" + label + ".gif").c_str());


  gBenchmark->Show("WWTemplate");       
} 

