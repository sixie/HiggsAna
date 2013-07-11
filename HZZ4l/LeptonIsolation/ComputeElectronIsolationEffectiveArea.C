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
#include <MitStyle.h>
#include "TLegend.h"
#include "TProfile.h"
#include "TF1.h"

// define structures to read in ntuple
#include "HiggsAna/Ntupler/interface/HiggsAnaDefs.hh"
#include "HiggsAna/Ntupler/interface/TEventInfo.hh"
#include "HiggsAna/Ntupler/interface/TElectron.hh"
#include "HiggsAna/Ntupler/interface/TPhoton.hh"
#include "HiggsAna/Ntupler/interface/TMuon.hh"
#include "HiggsAna/Ntupler/interface/TJet.hh"
#include "HiggsAna/Ntupler/interface/TPFCandidate.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

// helper functions for lepton ID selection
#include "HiggsAna/Utils/LeptonIDCuts.hh"

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
//Make Projection Graph
//*************************************************************************************************
TGraphAsymmErrors *MakeMeanVsNPUGraph(TH2F *hist2D, string name, string axisLabel = "") {

  const UInt_t nPoints = hist2D->GetYaxis()->GetNbins();
  double NPU[nPoints];
  double NPUErr[nPoints];
  double MeanValue[nPoints];
  double MeanErr[nPoints];

//   cout << "Npoints: " << nPoints << endl;

  TFile *file = new TFile("Test.root", "UPDATE");
  
  for(UInt_t b=0; b < nPoints; ++b) {
    TH1F *Projection = (TH1F*)hist2D->ProjectionX("_px", b,b+1);
    NPU[b] = (hist2D->GetYaxis()->GetBinUpEdge(b) + hist2D->GetYaxis()->GetBinLowEdge(b)) / 2;
    NPUErr[b] = 0.1;
    MeanValue[b] = Projection->GetMean();
    MeanErr[b] = Projection->GetMeanError();

//     cout << b << " : " << NPU[b] << " : " << MeanValue[b] << " " << MeanErr[b] << endl;
    char buffer[200];
    sprintf(buffer,"%d",b);
    if (name == "Ele_HoverE_VS_NPU_Fall11ZeeMC_Graph" || name == "Ele_HcalDepth1OverEcal_VS_NPU_Fall11ZeeMC_Graph" ||name == "Ele_HcalDepth2OverEcal_VS_NPU_Fall11ZeeMC_Graph") {
      file->WriteTObject(Projection , (name + "_Proj" + string(buffer)).c_str(), "WriteDelete");  
    }
  }
  file->Close();  

  TGraphAsymmErrors *Graph = new TGraphAsymmErrors (nPoints,  NPU, MeanValue, NPUErr, NPUErr, MeanErr, MeanErr);
  Graph->SetName(name.c_str());
  Graph->SetTitle("");
  Graph->SetMarkerColor(kBlack);
  Graph->GetXaxis()->SetTitleOffset(1.02);
  Graph->GetXaxis()->SetTitle("Number of Pileup Events");
  Graph->GetYaxis()->SetTitleOffset(1.05);
  Graph->GetYaxis()->SetTitle(axisLabel.c_str());

  return Graph;

}


//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void FillIsolationHistograms(string ElectronFile, string Label, Int_t EtaBin = -1)
{  

  string label = "";
  if (Label != "") label = "_" + Label;

  if (EtaBin == 0) {
    label = label + "_Eta0To1";
  }
  if (EtaBin == 1) {
    label = label + "_Eta1To1p479";
  }
  if (EtaBin == 2) {
    label = label + "_Eta1p479To2";
  }
  if (EtaBin == 3) {
    label = label + "_Eta2To2p2";
  }
  if (EtaBin == 4) {
    label = label + "_Eta2p2To2p3";
  }
  if (EtaBin == 5) {
    label = label + "_Eta2p3To2p4";
  }
  if (EtaBin == 6) {
    label = label + "_Eta2p4To2p5";
  }

  //*****************************************************************************************
  //Plotting Setup
  //*****************************************************************************************

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH2F *Ele_ChargedIso03_VS_NPU = new TH2F(("Ele_ChargedIso03_VS_NPU"+label).c_str(), "; ChargedIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_NeutralHadronIso03_VS_NPU = new TH2F(("Ele_NeutralHadronIso03_VS_NPU"+label).c_str(), "; NeutralHadronIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_GammaIso03_VS_NPU = new TH2F(("Ele_GammaIso03_VS_NPU"+label).c_str(), "; GammaIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_GammaIsoVetoEtaStrip03_VS_NPU = new TH2F(("Ele_GammaIsoVetoEtaStrip03_VS_NPU"+label).c_str(), "; GammaIsoVetoEtaStrip03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_ChargedIso04_VS_NPU = new TH2F(("Ele_ChargedIso04_VS_NPU"+label).c_str(), "; ChargedIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_NeutralHadronIso04_VS_NPU = new TH2F(("Ele_NeutralHadronIso04_VS_NPU"+label).c_str(), "; NeutralHadronIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_GammaIso04_VS_NPU = new TH2F(("Ele_GammaIso04_VS_NPU"+label).c_str(), "; GammaIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_GammaIsoVetoEtaStrip04_VS_NPU = new TH2F(("Ele_GammaIsoVetoEtaStrip04_VS_NPU"+label).c_str(), "; GammaIsoVetoEtaStrip04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_NeutralHadronIso007_VS_NPU = new TH2F(("Ele_NeutralHadronIso007_VS_NPU"+label).c_str(), "; NeutralHadronIso007 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 100, 50, -0.5, 49.5);
  TH2F *Ele_HoverE_VS_NPU = new TH2F(("Ele_HoverE_VS_NPU"+label).c_str(), "; HoverE ; Number of Pileup Events ; Number of Events ",  1000, 0 , 20, 50, -0.5, 49.5);
  TH2F *Ele_HcalDepth1OverEcal_VS_NPU = new TH2F(("Ele_HcalDepth1OverEcal_VS_NPU"+label).c_str(), "; HcalDepth1OverEcal ; Number of Pileup Events ; Number of Events ",  1000, 0 , 20, 50, -0.5, 49.5);
  TH2F *Ele_HcalDepth2OverEcal_VS_NPU = new TH2F(("Ele_HcalDepth2OverEcal_VS_NPU"+label).c_str(), "; HcalDepth2OverEcal ; Number of Pileup Events ; Number of Events ",  1000, 0 , 20, 50, -0.5, 49.5);
  TH2F *Ele_TrkIso03_VS_NPU = new TH2F(("Ele_TrkIso03_VS_NPU"+label).c_str(), "; TrkIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_EMIso03_VS_NPU = new TH2F(("Ele_EMIso03_VS_NPU"+label).c_str(), "; EMIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_HadIso03_VS_NPU = new TH2F(("Ele_HadIso03_VS_NPU"+label).c_str(), "; HadIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_TrkIso04_VS_NPU = new TH2F(("Ele_TrkIso04_VS_NPU"+label).c_str(), "; TrkIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_EMIso04_VS_NPU = new TH2F(("Ele_EMIso04_VS_NPU"+label).c_str(), "; EMIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_HadIso04_VS_NPU = new TH2F(("Ele_HadIso04_VS_NPU"+label).c_str(), "; HadIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_Rho_VS_NPU = new TH2F(("Ele_Rho_VS_NPU"+label).c_str(), "; Rho ; Number of Pileup Events ; Number of Events ",  1000, 0 , 200, 50, -0.5, 49.5);

  TH2F *Ele_ChargedIso03_VS_NVtx = new TH2F(("Ele_ChargedIso03_VS_NVtx"+label).c_str(), "; ChargedIso03 ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_NeutralHadronIso03_VS_NVtx = new TH2F(("Ele_NeutralHadronIso03_VS_NVtx"+label).c_str(), "; NeutralHadronIso03 ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_GammaIso03_VS_NVtx = new TH2F(("Ele_GammaIso03_VS_NVtx"+label).c_str(), "; GammaIso03 ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_GammaIsoVetoEtaStrip03_VS_NVtx = new TH2F(("Ele_GammaIsoVetoEtaStrip03_VS_NVtx"+label).c_str(), "; GammaIsoVetoEtaStrip03 ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_ChargedIso04_VS_NVtx = new TH2F(("Ele_ChargedIso04_VS_NVtx"+label).c_str(), "; ChargedIso04 ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_NeutralHadronIso04_VS_NVtx = new TH2F(("Ele_NeutralHadronIso04_VS_NVtx"+label).c_str(), "; NeutralHadronIso04 ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_GammaIso04_VS_NVtx = new TH2F(("Ele_GammaIso04_VS_NVtx"+label).c_str(), "; GammaIso04 ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_GammaIsoVetoEtaStrip04_VS_NVtx = new TH2F(("Ele_GammaIsoVetoEtaStrip04_VS_NVtx"+label).c_str(), "; GammaIsoVetoEtaStrip04 ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_NeutralHadronIso007_VS_NVtx = new TH2F(("Ele_NeutralHadronIso007_VS_NVtx"+label).c_str(), "; NeutralHadronIso007 ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 100, 50, -0.5, 49.5);
  TH2F *Ele_HoverE_VS_NVtx = new TH2F(("Ele_HoverE_VS_NVtx"+label).c_str(), "; HoverE ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 20, 50, -0.5, 49.5);
  TH2F *Ele_HcalDepth1OverEcal_VS_NVtx = new TH2F(("Ele_HcalDepth1OverEcal_VS_NVtx"+label).c_str(), "; HcalDepth1OverEcal ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 20, 50, -0.5, 49.5);
  TH2F *Ele_HcalDepth2OverEcal_VS_NVtx = new TH2F(("Ele_HcalDepth2OverEcal_VS_NVtx"+label).c_str(), "; HcalDepth2OverEcal ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 20, 50, -0.5, 49.5);
  TH2F *Ele_TrkIso03_VS_NVtx = new TH2F(("Ele_TrkIso03_VS_NVtx"+label).c_str(), "; TrkIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_EMIso03_VS_NVtx = new TH2F(("Ele_EMIso03_VS_NVtx"+label).c_str(), "; EMIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_HadIso03_VS_NVtx = new TH2F(("Ele_HadIso03_VS_NVtx"+label).c_str(), "; HadIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_TrkIso04_VS_NVtx = new TH2F(("Ele_TrkIso04_VS_NVtx"+label).c_str(), "; TrkIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_EMIso04_VS_NVtx = new TH2F(("Ele_EMIso04_VS_NVtx"+label).c_str(), "; EMIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_HadIso04_VS_NVtx = new TH2F(("Ele_HadIso04_VS_NVtx"+label).c_str(), "; HadIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Ele_Rho_VS_NVtx = new TH2F(("Ele_Rho_VS_NVtx"+label).c_str(), "; Rho ; Number of Reconstructed Primary Vertices ; Number of Events ",  1000, 0 , 200, 50, -0.5, 49.5);

  TH2F* Ele_ChargedIso_DR0p0To0p1_VS_NVtx = new TH2F(("Ele_ChargedIso_DR0p0To0p1_VS_NVtx"+Label).c_str(), ";ChargedIso_DR0p0To0p1;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_ChargedIso_DR0p1To0p2_VS_NVtx = new TH2F(("Ele_ChargedIso_DR0p1To0p2_VS_NVtx"+Label).c_str(), ";ChargedIso_DR0p1To0p2;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_ChargedIso_DR0p2To0p3_VS_NVtx = new TH2F(("Ele_ChargedIso_DR0p2To0p3_VS_NVtx"+Label).c_str(), ";ChargedIso_DR0p2To0p3;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_ChargedIso_DR0p3To0p4_VS_NVtx = new TH2F(("Ele_ChargedIso_DR0p3To0p4_VS_NVtx"+Label).c_str(), ";ChargedIso_DR0p3To0p4;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_ChargedIso_DR0p4To0p5_VS_NVtx = new TH2F(("Ele_ChargedIso_DR0p4To0p5_VS_NVtx"+Label).c_str(), ";ChargedIso_DR0p4To0p5;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_ChargedIso_DR0p5To0p7_VS_NVtx = new TH2F(("Ele_ChargedIso_DR0p5To0p7_VS_NVtx"+Label).c_str(), ";ChargedIso_DR0p5To0p7;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_ChargedIso_DR0p7To1p0_VS_NVtx = new TH2F(("Ele_ChargedIso_DR0p7To1p0_VS_NVtx"+Label).c_str(), ";ChargedIso_DR0p7To1p0;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_GammaIso_DR0p0To0p1_VS_NVtx = new TH2F(("Ele_GammaIso_DR0p0To0p1_VS_NVtx"+Label).c_str(), ";GammaIso_DR0p0To0p1;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_GammaIso_DR0p1To0p2_VS_NVtx = new TH2F(("Ele_GammaIso_DR0p1To0p2_VS_NVtx"+Label).c_str(), ";GammaIso_DR0p1To0p2;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_GammaIso_DR0p2To0p3_VS_NVtx = new TH2F(("Ele_GammaIso_DR0p2To0p3_VS_NVtx"+Label).c_str(), ";GammaIso_DR0p2To0p3;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_GammaIso_DR0p3To0p4_VS_NVtx = new TH2F(("Ele_GammaIso_DR0p3To0p4_VS_NVtx"+Label).c_str(), ";GammaIso_DR0p3To0p4;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_GammaIso_DR0p4To0p5_VS_NVtx = new TH2F(("Ele_GammaIso_DR0p4To0p5_VS_NVtx"+Label).c_str(), ";GammaIso_DR0p4To0p5;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_GammaIso_DR0p5To0p7_VS_NVtx = new TH2F(("Ele_GammaIso_DR0p5To0p7_VS_NVtx"+Label).c_str(), ";GammaIso_DR0p5To0p7;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_GammaIso_DR0p7To1p0_VS_NVtx = new TH2F(("Ele_GammaIso_DR0p7To1p0_VS_NVtx"+Label).c_str(), ";GammaIso_DR0p7To1p0;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_NeutralHadronIso_DR0p0To0p1_VS_NVtx = new TH2F(("Ele_NeutralHadronIso_DR0p0To0p1_VS_NVtx"+Label).c_str(), ";NeutralHadronIso_DR0p0To0p1;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_NeutralHadronIso_DR0p1To0p2_VS_NVtx = new TH2F(("Ele_NeutralHadronIso_DR0p1To0p2_VS_NVtx"+Label).c_str(), ";NeutralHadronIso_DR0p1To0p2;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_NeutralHadronIso_DR0p2To0p3_VS_NVtx = new TH2F(("Ele_NeutralHadronIso_DR0p2To0p3_VS_NVtx"+Label).c_str(), ";NeutralHadronIso_DR0p2To0p3;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_NeutralHadronIso_DR0p3To0p4_VS_NVtx = new TH2F(("Ele_NeutralHadronIso_DR0p3To0p4_VS_NVtx"+Label).c_str(), ";NeutralHadronIso_DR0p3To0p4;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_NeutralHadronIso_DR0p4To0p5_VS_NVtx = new TH2F(("Ele_NeutralHadronIso_DR0p4To0p5_VS_NVtx"+Label).c_str(), ";NeutralHadronIso_DR0p4To0p5;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_NeutralHadronIso_DR0p5To0p7_VS_NVtx = new TH2F(("Ele_NeutralHadronIso_DR0p5To0p7_VS_NVtx"+Label).c_str(), ";NeutralHadronIso_DR0p5To0p7;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_NeutralHadronIso_DR0p7To1p0_VS_NVtx = new TH2F(("Ele_NeutralHadronIso_DR0p7To1p0_VS_NVtx"+Label).c_str(), ";NeutralHadronIso_DR0p7To1p0;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);

  TH2F* Ele_ChargedIso_DR0p0To0p3_VS_NVtx = new TH2F(("Ele_ChargedIso_DR0p0To0p3_VS_NVtx"+Label).c_str(), ";ChargedIso_DR0p0To0p3;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_ChargedIso_DR0p0To0p4_VS_NVtx = new TH2F(("Ele_ChargedIso_DR0p0To0p4_VS_NVtx"+Label).c_str(), ";ChargedIso_DR0p0To0p4;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_GammaIso_DR0p0To0p3_VS_NVtx = new TH2F(("Ele_GammaIso_DR0p0To0p3_VS_NVtx"+Label).c_str(), ";GammaIso_DR0p0To0p3;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_GammaIso_DR0p0To0p4_VS_NVtx = new TH2F(("Ele_GammaIso_DR0p0To0p4_VS_NVtx"+Label).c_str(), ";GammaIso_DR0p0To0p4;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_NeutralHadronIso_DR0p0To0p3_VS_NVtx = new TH2F(("Ele_NeutralHadronIso_DR0p0To0p3_VS_NVtx"+Label).c_str(), ";NeutralHadronIso_DR0p0To0p3;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F* Ele_NeutralHadronIso_DR0p0To0p4_VS_NVtx = new TH2F(("Ele_NeutralHadronIso_DR0p0To0p4_VS_NVtx"+Label).c_str(), ";NeutralHadronIso_DR0p0To0p4;Number of Events",  1000, 0 , 1000, 50, -0.5, 49.5);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  TClonesArray *photonArr = new TClonesArray("mithep::TPhoton");
  TClonesArray *pfcandidateArr = new TClonesArray("mithep::TPFCandidate");

  Int_t NEvents = 0;
  Bool_t isMC = kTRUE;

  vector<string> inputfiles;
  if (ElectronFile == "Fall11MC") {
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-h115zz4l-gf-v14b-pu_noskim_0000.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-h120zz4l-gf-v14b-pu_noskim_0000.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-h130zz4l-gf-v14b-pu_noskim_0000.root");
  } else if (ElectronFile == "Summer11MC") {
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_s11-h115zz4l-gf-v1g1-pu_noskim_0000.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_s11-h120zz4l-gf-v1g1-pu_noskim_0000.root");
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_s11-h130zz4l-gf-v1g1-pu_noskim_0000.root");
  }else if (ElectronFile == "Data2011") {
    isMC = kFALSE;
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-m10-v1.TightPlusReco.root"); 
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-pr-v4.TightPlusReco.root"); 
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-a05-v1.TightPlusReco.root"); 
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11a-del-o03-v1.TightPlusReco.root"); 
    inputfiles.push_back("/home/sixie/hist/HZZ4lNtuples/data/AllNtuple_HZZ4lNtuple_r11b-del-pr-v1.TightPlusReco.root");    
  } else {
    inputfiles.push_back(ElectronFile);
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
    eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);         pfcandidateBr = eventTree->GetBranch("PFCandidate");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;

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

      if (isMC) {
        for(Int_t i=0; i<electronArr->GetEntries(); i++) {
          const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
          
          //********************************************************
          // Select MC Truth Muons
          //********************************************************
          
          if (!((UInt_t(abs(max(0,ele->isMCReal))) & 1) == 1)) continue;

          if (EtaBin == 0) {
            if (!(fabs(ele->scEta) < 1.0)) continue;
          }
          if (EtaBin == 1) {
            if (!(fabs(ele->scEta) >= 1.0 && fabs(ele->scEta) < 1.479)) continue;
          }
          if (EtaBin == 2) {
            if (!(fabs(ele->scEta) >= 1.479 && fabs(ele->scEta) < 2.0)) continue;
          }
          if (EtaBin == 3) {
            if (!(fabs(ele->scEta) >= 2.0 && fabs(ele->scEta) < 2.2)) continue;
          }
          if (EtaBin == 4) {
            if (!(fabs(ele->scEta) >= 2.2 && fabs(ele->scEta) < 2.3)) continue;
          }
          if (EtaBin == 5) {
            if (!(fabs(ele->scEta) >= 2.3 && fabs(ele->scEta) < 2.4)) continue;
          }
          if (EtaBin == 6) {
            if (!(fabs(ele->scEta) >= 2.4 && fabs(ele->scEta) < 2.5)) continue;
          }

          Ele_ChargedIso03_VS_NPU->Fill(ele->ChargedIso03,info->nPUEvents);
          Ele_NeutralHadronIso03_VS_NPU->Fill(ele->NeutralHadronIso03_05Threshold,info->nPUEvents); 
          Ele_GammaIso03_VS_NPU->Fill(ele->GammaIso03_05Threshold,info->nPUEvents);
          Ele_GammaIsoVetoEtaStrip03_VS_NPU->Fill(ele->GammaIsoVetoEtaStrip03_05Threshold,info->nPUEvents);
          Ele_ChargedIso04_VS_NPU->Fill(ele->ChargedIso04,info->nPUEvents);
          Ele_NeutralHadronIso04_VS_NPU->Fill(ele->NeutralHadronIso04_05Threshold,info->nPUEvents); 
          Ele_GammaIso04_VS_NPU->Fill(ele->GammaIso04_05Threshold,info->nPUEvents); 
          Ele_GammaIsoVetoEtaStrip04_VS_NPU->Fill(ele->GammaIsoVetoEtaStrip04_05Threshold,info->nPUEvents);
          Ele_NeutralHadronIso007_VS_NPU->Fill(ele->NeutralHadronIso007_05Threshold,info->nPUEvents);
          Ele_HoverE_VS_NPU->Fill(ele->HoverE,info->nPUEvents);
          Ele_HcalDepth1OverEcal_VS_NPU->Fill(ele->HcalDepth1OverEcal,info->nPUEvents);
          Ele_HcalDepth2OverEcal_VS_NPU->Fill(ele->HcalDepth2OverEcal,info->nPUEvents);
          Ele_TrkIso03_VS_NPU->Fill(ele->trkIso03,info->nPUEvents);
          Ele_EMIso03_VS_NPU->Fill(ele->emIso03,info->nPUEvents);
          Ele_HadIso03_VS_NPU->Fill(ele->hadIso03,info->nPUEvents);
          Ele_TrkIso04_VS_NPU->Fill(ele->trkIso04,info->nPUEvents);
          Ele_EMIso04_VS_NPU->Fill(ele->emIso04,info->nPUEvents);
          Ele_HadIso04_VS_NPU->Fill(ele->hadIso04,info->nPUEvents);

          Ele_Rho_VS_NPU->Fill(info->PileupEnergyDensity,info->nPUEvents);
  
          Ele_ChargedIso03_VS_NVtx->Fill(ele->ChargedIso03,info->nPV0);
          Ele_NeutralHadronIso03_VS_NVtx->Fill(ele->NeutralHadronIso03_05Threshold,info->nPV0); 
          Ele_GammaIso03_VS_NVtx->Fill(ele->GammaIso03_05Threshold,info->nPV0); 
          Ele_GammaIsoVetoEtaStrip03_VS_NVtx->Fill(ele->GammaIsoVetoEtaStrip03_05Threshold,info->nPV0);
          Ele_ChargedIso04_VS_NVtx->Fill(ele->ChargedIso04,info->nPV0);
          Ele_NeutralHadronIso04_VS_NVtx->Fill(ele->NeutralHadronIso04_05Threshold,info->nPV0); 
          Ele_GammaIso04_VS_NVtx->Fill(ele->GammaIso04_05Threshold,info->nPV0); 
          Ele_GammaIsoVetoEtaStrip04_VS_NVtx->Fill(ele->GammaIsoVetoEtaStrip04_05Threshold,info->nPV0);
          Ele_NeutralHadronIso007_VS_NVtx->Fill(ele->NeutralHadronIso007_05Threshold,info->nPV0);
          Ele_HoverE_VS_NVtx->Fill(ele->HoverE,info->nPV0);
          Ele_HcalDepth1OverEcal_VS_NVtx->Fill(ele->HcalDepth1OverEcal,info->nPV0);
          Ele_HcalDepth2OverEcal_VS_NVtx->Fill(ele->HcalDepth2OverEcal,info->nPV0);
          Ele_TrkIso03_VS_NVtx->Fill(ele->trkIso03,info->nPV0);
          Ele_EMIso03_VS_NVtx->Fill(ele->emIso03,info->nPV0);
          Ele_HadIso03_VS_NVtx->Fill(ele->hadIso03,info->nPV0);
          Ele_TrkIso04_VS_NVtx->Fill(ele->trkIso04,info->nPV0);
          Ele_EMIso04_VS_NVtx->Fill(ele->emIso04,info->nPV0);
          Ele_HadIso04_VS_NVtx->Fill(ele->hadIso04,info->nPV0);
          Ele_Rho_VS_NVtx->Fill(info->PileupEnergyDensity,info->nPV0);

            //*************************************************
            //Isolation STUFF
            //*************************************************
            Double_t tmpChargedIso_DR0p0To0p1  = 0;
            Double_t tmpChargedIso_DR0p1To0p2  = 0;
            Double_t tmpChargedIso_DR0p2To0p3  = 0;
            Double_t tmpChargedIso_DR0p3To0p4  = 0;
            Double_t tmpChargedIso_DR0p4To0p5  = 0;
            Double_t tmpChargedIso_DR0p5To0p7  = 0; 
            Double_t tmpChargedIso_DR0p7To1p0  = 0; 
            Double_t tmpGammaIso_DR0p0To0p1  = 0;
            Double_t tmpGammaIso_DR0p1To0p2  = 0;
            Double_t tmpGammaIso_DR0p2To0p3  = 0;
            Double_t tmpGammaIso_DR0p3To0p4  = 0;
            Double_t tmpGammaIso_DR0p4To0p5  = 0;
            Double_t tmpGammaIso_DR0p5To0p7  = 0; 
            Double_t tmpGammaIso_DR0p7To1p0  = 0; 
            Double_t tmpNeutralHadronIso_DR0p0To0p1  = 0;
            Double_t tmpNeutralHadronIso_DR0p1To0p2  = 0;
            Double_t tmpNeutralHadronIso_DR0p2To0p3  = 0;
            Double_t tmpNeutralHadronIso_DR0p3To0p4  = 0;
            Double_t tmpNeutralHadronIso_DR0p4To0p5  = 0;
            Double_t tmpNeutralHadronIso_DR0p5To0p7  = 0;
            Double_t tmpNeutralHadronIso_DR0p7To1p0  = 0;



            for(Int_t k=0; k<pfcandidateArr->GetEntries(); ++k) {
              const mithep::TPFCandidate *pf = (mithep::TPFCandidate*)((*pfcandidateArr)[k]);
              if (pf->matchedObjectType == 11 && pf->matchedObjectIndex == i) continue;
              Double_t deta = (ele->eta - pf->eta);
              Double_t dphi = mithep::MathUtils::DeltaPhi(Double_t(ele->phi),Double_t(pf->phi));
              Double_t dr = fabs(mithep::MathUtils::DeltaR(ele->phi,ele->eta, pf->phi, pf->eta));
              if (dr > 1.0) continue;

              //charged
              if (pf->q != 0) {
                if (abs(pf->dz) > 0.1) continue;
                //************************************************************
                // Veto any PFmuon, or PFEle
                if (pf->pfType == eElectron || pf->pfType == eMuon) continue;
                //************************************************************
                //************************************************************
                // Footprint Veto only for electrons that failed PF ele mva cut
                if (fabs(ele->scEta) > 1.479) {  
                  if (dr < 0.015) continue;
                }
                //************************************************************

                if (dr < 0.1) tmpChargedIso_DR0p0To0p1 += pf->pt;
                if (dr >= 0.1 && dr < 0.2) tmpChargedIso_DR0p1To0p2 += pf->pt;
                if (dr >= 0.2 && dr < 0.3) tmpChargedIso_DR0p2To0p3 += pf->pt;
                if (dr >= 0.3 && dr < 0.4) tmpChargedIso_DR0p3To0p4 += pf->pt;
                if (dr >= 0.4 && dr < 0.5) tmpChargedIso_DR0p4To0p5 += pf->pt;
                if (dr >= 0.5 && dr < 0.7) tmpChargedIso_DR0p5To0p7 += pf->pt;
                if (dr >= 0.7 && dr < 1.0) tmpChargedIso_DR0p7To1p0 += pf->pt;
             }
              else if (pf->pfType == eGamma) {
                //************************************************************
                // Footprint Veto
                if (fabs(ele->scEta) > 1.479) {                  
                  if (dr < 0.08) continue;
                }
                //************************************************************                
                if (dr < 0.1) tmpGammaIso_DR0p0To0p1 += pf->pt;
                if (dr >= 0.1 && dr < 0.2) tmpGammaIso_DR0p1To0p2 += pf->pt;
                if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += pf->pt;
                if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += pf->pt;
                if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += pf->pt;
                if (dr >= 0.5 && dr < 0.7) tmpGammaIso_DR0p5To0p7 += pf->pt;
                if (dr >= 0.7 && dr < 1.0) tmpGammaIso_DR0p7To1p0 += pf->pt;
              }
              else {
                if (dr < 0.1) tmpNeutralHadronIso_DR0p0To0p1 += pf->pt;
                if (dr >= 0.1 && dr < 0.2) tmpNeutralHadronIso_DR0p1To0p2 += pf->pt;
                if (dr >= 0.2 && dr < 0.3) tmpNeutralHadronIso_DR0p2To0p3 += pf->pt;
                if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += pf->pt;
                if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += pf->pt;
                if (dr >= 0.5 && dr < 0.7) tmpNeutralHadronIso_DR0p5To0p7 += pf->pt;
                if (dr >= 0.7 && dr < 1.0) tmpNeutralHadronIso_DR0p7To1p0 += pf->pt;                
              }
            }

            //***********************
            //Fill isolation rings
            //***********************
            
            Ele_ChargedIso_DR0p0To0p1_VS_NVtx->Fill(tmpChargedIso_DR0p0To0p1,info->nPV0);
            Ele_ChargedIso_DR0p1To0p2_VS_NVtx->Fill(tmpChargedIso_DR0p1To0p2,info->nPV0);
            Ele_ChargedIso_DR0p2To0p3_VS_NVtx->Fill(tmpChargedIso_DR0p2To0p3,info->nPV0);
            Ele_ChargedIso_DR0p3To0p4_VS_NVtx->Fill(tmpChargedIso_DR0p3To0p4,info->nPV0);
            Ele_ChargedIso_DR0p4To0p5_VS_NVtx->Fill(tmpChargedIso_DR0p4To0p5,info->nPV0);
            Ele_ChargedIso_DR0p5To0p7_VS_NVtx->Fill(tmpChargedIso_DR0p5To0p7,info->nPV0);
            Ele_ChargedIso_DR0p7To1p0_VS_NVtx->Fill(tmpChargedIso_DR0p7To1p0,info->nPV0);
            Ele_GammaIso_DR0p0To0p1_VS_NVtx->Fill(tmpGammaIso_DR0p0To0p1,info->nPV0);
            Ele_GammaIso_DR0p1To0p2_VS_NVtx->Fill(tmpGammaIso_DR0p1To0p2,info->nPV0);
            Ele_GammaIso_DR0p2To0p3_VS_NVtx->Fill(tmpGammaIso_DR0p2To0p3,info->nPV0);
            Ele_GammaIso_DR0p3To0p4_VS_NVtx->Fill(tmpGammaIso_DR0p3To0p4,info->nPV0);
            Ele_GammaIso_DR0p4To0p5_VS_NVtx->Fill(tmpGammaIso_DR0p4To0p5,info->nPV0);
            Ele_GammaIso_DR0p5To0p7_VS_NVtx->Fill(tmpGammaIso_DR0p5To0p7,info->nPV0);
            Ele_GammaIso_DR0p7To1p0_VS_NVtx->Fill(tmpGammaIso_DR0p7To1p0,info->nPV0);
            Ele_NeutralHadronIso_DR0p0To0p1_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p0To0p1,info->nPV0);
            Ele_NeutralHadronIso_DR0p1To0p2_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p1To0p2,info->nPV0);
            Ele_NeutralHadronIso_DR0p2To0p3_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p2To0p3,info->nPV0);
            Ele_NeutralHadronIso_DR0p3To0p4_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p3To0p4,info->nPV0);
            Ele_NeutralHadronIso_DR0p4To0p5_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p4To0p5,info->nPV0);
            Ele_NeutralHadronIso_DR0p5To0p7_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p5To0p7,info->nPV0);
            Ele_NeutralHadronIso_DR0p7To1p0_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p7To1p0,info->nPV0);

            Ele_ChargedIso_DR0p0To0p3_VS_NVtx->Fill(tmpChargedIso_DR0p0To0p1 + tmpChargedIso_DR0p1To0p2 + tmpChargedIso_DR0p2To0p3,info->nPV0);
            Ele_ChargedIso_DR0p0To0p4_VS_NVtx->Fill(tmpChargedIso_DR0p0To0p1 + tmpChargedIso_DR0p1To0p2 + tmpChargedIso_DR0p2To0p3 + tmpChargedIso_DR0p3To0p4,info->nPV0);
            Ele_GammaIso_DR0p0To0p3_VS_NVtx->Fill(tmpGammaIso_DR0p0To0p1 + tmpGammaIso_DR0p1To0p2 + tmpGammaIso_DR0p2To0p3,info->nPV0);
            Ele_GammaIso_DR0p0To0p4_VS_NVtx->Fill(tmpGammaIso_DR0p0To0p1 + tmpGammaIso_DR0p1To0p2 + tmpGammaIso_DR0p2To0p3 + tmpGammaIso_DR0p3To0p4,info->nPV0);
            Ele_NeutralHadronIso_DR0p0To0p3_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p0To0p1 + tmpNeutralHadronIso_DR0p1To0p2 + tmpNeutralHadronIso_DR0p2To0p3,info->nPV0);
            Ele_NeutralHadronIso_DR0p0To0p4_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p0To0p1 + tmpNeutralHadronIso_DR0p1To0p2 + tmpNeutralHadronIso_DR0p2To0p3 + tmpNeutralHadronIso_DR0p3To0p4,info->nPV0);


        } //loop over electrons
      } else {

        Bool_t FilledRho = kFALSE;

        //Require single mu triggers
        ULong_t trigger = kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 |
          kHLT_Ele32_CaloIdL_CaloIsoVL_SC17 |
          kHLT_Ele17_CaloIdL_CaloIsoVL;
        
        if(!(info->triggerBits & trigger)) continue;      
    


        for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
          const mithep::TElectron *tag = (mithep::TElectron*)((*electronArr)[i]);
          
          if(tag->pt          < 20)  continue;
          if(fabs(tag->scEta) > 2.5) continue;
          if (!passCutBasedTightEleID(tag)) continue;          
          
          //Match Tag with Tight Leg of T&P Trigger
          if(
            !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (tag->hltMatchBits & kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT)) &&
            !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (tag->hltMatchBits & kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT)) &&
            !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdL_CaloIsoVL)) &&
            !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT)) &&
            !((info->triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL) && (tag->hltMatchBits & kHLTObject_Ele17_CaloIdL_CaloIsoVL))
            )
            continue;
          
          const Double_t m = 0.000511;
          TLorentzVector vtag;
          vtag.SetPtEtaPhiM(tag->pt, tag->eta, tag->phi, m);
          
          for(Int_t j=0; j<electronArr->GetEntriesFast(); j++) {
            if(i==j) continue;
            
            const mithep::TElectron *probe = (mithep::TElectron*)((*electronArr)[j]);
            if(probe->q == tag->q) continue;
            
            //match probe with loose leg of T&P trigger
            if(
              !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (probe->hltMatchBits & kHLTObject_SC8)) &&
              !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (probe->hltMatchBits & kHLTObject_Ele8)) &&
              !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (probe->hltMatchBits & kHLTObject_SC17))
              )
              continue;
            
            if(fabs(probe->scEta) > 2.5) continue;

            if (EtaBin == 0) {
              if (!(fabs(probe->scEta) < 1.0)) continue;
            }
            if (EtaBin == 1) {
              if (!(fabs(probe->scEta) >= 1.0 && fabs(probe->scEta) < 1.479)) continue;
            }
            if (EtaBin == 2) {
              if (!(fabs(probe->scEta) >= 1.479 && fabs(probe->scEta) < 2.0)) continue;
            }
            if (EtaBin == 3) {
              if (!(fabs(probe->scEta) >= 2.0 && fabs(probe->scEta) < 2.2)) continue;
            }
            if (EtaBin == 4) {
              if (!(fabs(probe->scEta) >= 2.2 && fabs(probe->scEta) < 2.3)) continue;
            }
            if (EtaBin == 5) {
              if (!(fabs(probe->scEta) >= 2.3 && fabs(probe->scEta) < 2.4)) continue;
            }
            if (EtaBin == 6) {
              if (!(fabs(probe->scEta) >= 2.4 && fabs(probe->scEta) < 2.5)) continue;
            }
            
            TLorentzVector vprobe;
            vprobe.SetPtEtaPhiM(probe->pt, probe->eta, probe->phi, m);
            
            TLorentzVector vdielectron = vtag + vprobe;
            if((vdielectron.M()<75) || (vdielectron.M()>105)) continue;
            
            //probe is ok
            Ele_ChargedIso03_VS_NVtx->Fill(probe->ChargedIso03,info->nPV0);
            Ele_NeutralHadronIso03_VS_NVtx->Fill(probe->NeutralHadronIso03_05Threshold,info->nPV0); 
            Ele_GammaIso03_VS_NVtx->Fill(probe->GammaIso03_05Threshold,info->nPV0); 
            Ele_GammaIsoVetoEtaStrip03_VS_NVtx->Fill(probe->GammaIsoVetoEtaStrip03_05Threshold,info->nPV0); 
            Ele_ChargedIso04_VS_NVtx->Fill(probe->ChargedIso04,info->nPV0);
            Ele_NeutralHadronIso04_VS_NVtx->Fill(probe->NeutralHadronIso04_05Threshold,info->nPV0); 
            Ele_GammaIso04_VS_NVtx->Fill(probe->GammaIso04_05Threshold,info->nPV0); 
            Ele_GammaIsoVetoEtaStrip04_VS_NVtx->Fill(probe->GammaIsoVetoEtaStrip04_05Threshold,info->nPV0); 
            Ele_NeutralHadronIso007_VS_NVtx->Fill(probe->NeutralHadronIso007_05Threshold,info->nPV0); 
            Ele_HoverE_VS_NVtx->Fill(probe->HoverE,info->nPV0);
            Ele_HcalDepth1OverEcal_VS_NVtx->Fill(probe->HcalDepth1OverEcal,info->nPV0);
            Ele_HcalDepth2OverEcal_VS_NVtx->Fill(probe->HcalDepth2OverEcal,info->nPV0);
            Ele_TrkIso03_VS_NVtx->Fill(probe->trkIso03,info->nPV0);
            Ele_EMIso03_VS_NVtx->Fill(probe->emIso03,info->nPV0);
            Ele_HadIso03_VS_NVtx->Fill(probe->hadIso03,info->nPV0);
            Ele_TrkIso04_VS_NVtx->Fill(probe->trkIso04,info->nPV0);
            Ele_EMIso04_VS_NVtx->Fill(probe->emIso04,info->nPV0);
            Ele_HadIso04_VS_NVtx->Fill(probe->hadIso04,info->nPV0);
            
            //*************************************************
            //Isolation STUFF
            //*************************************************
            Double_t tmpChargedIso_DR0p0To0p1  = 0;
            Double_t tmpChargedIso_DR0p1To0p2  = 0;
            Double_t tmpChargedIso_DR0p2To0p3  = 0;
            Double_t tmpChargedIso_DR0p3To0p4  = 0;
            Double_t tmpChargedIso_DR0p4To0p5  = 0;
            Double_t tmpChargedIso_DR0p5To0p7  = 0; 
            Double_t tmpChargedIso_DR0p7To1p0  = 0; 
            Double_t tmpGammaIso_DR0p0To0p1  = 0;
            Double_t tmpGammaIso_DR0p1To0p2  = 0;
            Double_t tmpGammaIso_DR0p2To0p3  = 0;
            Double_t tmpGammaIso_DR0p3To0p4  = 0;
            Double_t tmpGammaIso_DR0p4To0p5  = 0;
            Double_t tmpGammaIso_DR0p5To0p7  = 0; 
            Double_t tmpGammaIso_DR0p7To1p0  = 0; 
            Double_t tmpNeutralHadronIso_DR0p0To0p1  = 0;
            Double_t tmpNeutralHadronIso_DR0p1To0p2  = 0;
            Double_t tmpNeutralHadronIso_DR0p2To0p3  = 0;
            Double_t tmpNeutralHadronIso_DR0p3To0p4  = 0;
            Double_t tmpNeutralHadronIso_DR0p4To0p5  = 0;
            Double_t tmpNeutralHadronIso_DR0p5To0p7  = 0;
            Double_t tmpNeutralHadronIso_DR0p7To1p0  = 0;



            for(Int_t k=0; k<pfcandidateArr->GetEntries(); ++k) {
              const mithep::TPFCandidate *pf = (mithep::TPFCandidate*)((*pfcandidateArr)[k]);
              if (pf->matchedObjectType == 11 && pf->matchedObjectIndex == j) continue;
              Double_t deta = (probe->eta - pf->eta);
              Double_t dphi = mithep::MathUtils::DeltaPhi(Double_t(probe->phi),Double_t(pf->phi));
              Double_t dr = fabs(mithep::MathUtils::DeltaR(probe->phi,probe->eta, pf->phi, pf->eta));
              if (dr > 1.0) continue;

              //charged
              if (pf->q != 0) {
                if (abs(pf->dz) > 0.1) continue;
                //************************************************************
                // Veto any PFmuon, or PFEle
                if (pf->pfType == eElectron || pf->pfType == eMuon) continue;
                //************************************************************
                //************************************************************
                // Footprint Veto
                if (fabs(probe->scEta) > 1.479) {  
                  if (dr < 0.015) continue;
                }
                //************************************************************
                if (dr < 0.1) tmpChargedIso_DR0p0To0p1 += pf->pt;
                if (dr >= 0.1 && dr < 0.2) tmpChargedIso_DR0p1To0p2 += pf->pt;
                if (dr >= 0.2 && dr < 0.3) tmpChargedIso_DR0p2To0p3 += pf->pt;
                if (dr >= 0.3 && dr < 0.4) tmpChargedIso_DR0p3To0p4 += pf->pt;
                if (dr >= 0.4 && dr < 0.5) tmpChargedIso_DR0p4To0p5 += pf->pt;
                if (dr >= 0.5 && dr < 0.7) tmpChargedIso_DR0p5To0p7 += pf->pt;
                if (dr >= 0.7 && dr < 1.0) tmpChargedIso_DR0p7To1p0 += pf->pt;
             }
              else if (pf->pfType == eGamma) {
                //************************************************************
                // Footprint Veto
                if (fabs(probe->scEta) > 1.479) {                  
                  if (dr < 0.08) continue;
                }
                //************************************************************
                if (dr < 0.1) tmpGammaIso_DR0p0To0p1 += pf->pt;
                if (dr >= 0.1 && dr < 0.2) tmpGammaIso_DR0p1To0p2 += pf->pt;
                if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += pf->pt;
                if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += pf->pt;
                if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += pf->pt;
                if (dr >= 0.5 && dr < 0.7) tmpGammaIso_DR0p5To0p7 += pf->pt;
                if (dr >= 0.7 && dr < 1.0) tmpGammaIso_DR0p7To1p0 += pf->pt;
              }
              else {
                if (dr < 0.1) tmpNeutralHadronIso_DR0p0To0p1 += pf->pt;
                if (dr >= 0.1 && dr < 0.2) tmpNeutralHadronIso_DR0p1To0p2 += pf->pt;
                if (dr >= 0.2 && dr < 0.3) tmpNeutralHadronIso_DR0p2To0p3 += pf->pt;
                if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += pf->pt;
                if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += pf->pt;
                if (dr >= 0.5 && dr < 0.7) tmpNeutralHadronIso_DR0p5To0p7 += pf->pt;
                if (dr >= 0.7 && dr < 1.0) tmpNeutralHadronIso_DR0p7To1p0 += pf->pt;
              }
            }

            //***********************
            //Fill isolation rings
            //***********************
            Ele_ChargedIso_DR0p0To0p1_VS_NVtx->Fill(tmpChargedIso_DR0p0To0p1,info->nPV0);
            Ele_ChargedIso_DR0p1To0p2_VS_NVtx->Fill(tmpChargedIso_DR0p1To0p2,info->nPV0);
            Ele_ChargedIso_DR0p2To0p3_VS_NVtx->Fill(tmpChargedIso_DR0p2To0p3,info->nPV0);
            Ele_ChargedIso_DR0p3To0p4_VS_NVtx->Fill(tmpChargedIso_DR0p3To0p4,info->nPV0);
            Ele_ChargedIso_DR0p4To0p5_VS_NVtx->Fill(tmpChargedIso_DR0p4To0p5,info->nPV0);
            Ele_ChargedIso_DR0p5To0p7_VS_NVtx->Fill(tmpChargedIso_DR0p5To0p7,info->nPV0);
            Ele_ChargedIso_DR0p7To1p0_VS_NVtx->Fill(tmpChargedIso_DR0p7To1p0,info->nPV0);
            Ele_GammaIso_DR0p0To0p1_VS_NVtx->Fill(tmpGammaIso_DR0p0To0p1,info->nPV0);
            Ele_GammaIso_DR0p1To0p2_VS_NVtx->Fill(tmpGammaIso_DR0p1To0p2,info->nPV0);
            Ele_GammaIso_DR0p2To0p3_VS_NVtx->Fill(tmpGammaIso_DR0p2To0p3,info->nPV0);
            Ele_GammaIso_DR0p3To0p4_VS_NVtx->Fill(tmpGammaIso_DR0p3To0p4,info->nPV0);
            Ele_GammaIso_DR0p4To0p5_VS_NVtx->Fill(tmpGammaIso_DR0p4To0p5,info->nPV0);
            Ele_GammaIso_DR0p5To0p7_VS_NVtx->Fill(tmpGammaIso_DR0p5To0p7,info->nPV0);
            Ele_GammaIso_DR0p7To1p0_VS_NVtx->Fill(tmpGammaIso_DR0p7To1p0,info->nPV0);
            Ele_NeutralHadronIso_DR0p0To0p1_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p0To0p1,info->nPV0);
            Ele_NeutralHadronIso_DR0p1To0p2_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p1To0p2,info->nPV0);
            Ele_NeutralHadronIso_DR0p2To0p3_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p2To0p3,info->nPV0);
            Ele_NeutralHadronIso_DR0p3To0p4_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p3To0p4,info->nPV0);
            Ele_NeutralHadronIso_DR0p4To0p5_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p4To0p5,info->nPV0);
            Ele_NeutralHadronIso_DR0p5To0p7_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p5To0p7,info->nPV0);
            Ele_NeutralHadronIso_DR0p7To1p0_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p7To1p0,info->nPV0);

            Ele_ChargedIso_DR0p0To0p3_VS_NVtx->Fill(tmpChargedIso_DR0p0To0p1 + tmpChargedIso_DR0p1To0p2 + tmpChargedIso_DR0p2To0p3,info->nPV0);
            Ele_ChargedIso_DR0p0To0p4_VS_NVtx->Fill(tmpChargedIso_DR0p0To0p1 + tmpChargedIso_DR0p1To0p2 + tmpChargedIso_DR0p2To0p3 + tmpChargedIso_DR0p3To0p4,info->nPV0);
            Ele_GammaIso_DR0p0To0p3_VS_NVtx->Fill(tmpGammaIso_DR0p0To0p1 + tmpGammaIso_DR0p1To0p2 + tmpGammaIso_DR0p2To0p3,info->nPV0);
            Ele_GammaIso_DR0p0To0p4_VS_NVtx->Fill(tmpGammaIso_DR0p0To0p1 + tmpGammaIso_DR0p1To0p2 + tmpGammaIso_DR0p2To0p3 + tmpGammaIso_DR0p3To0p4,info->nPV0);
            Ele_NeutralHadronIso_DR0p0To0p3_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p0To0p1 + tmpNeutralHadronIso_DR0p1To0p2 + tmpNeutralHadronIso_DR0p2To0p3,info->nPV0);
            Ele_NeutralHadronIso_DR0p0To0p4_VS_NVtx->Fill(tmpNeutralHadronIso_DR0p0To0p1 + tmpNeutralHadronIso_DR0p1To0p2 + tmpNeutralHadronIso_DR0p2To0p3 + tmpNeutralHadronIso_DR0p3To0p4,info->nPV0);

            if (!FilledRho) {
              Ele_Rho_VS_NVtx->Fill(info->PileupEnergyDensity,info->nPV0);
              FilledRho = kTRUE;
            }
            
          }//loop over probes
        } //loop over tags        
      } // if data

    } //loop over events

  } //end loop over files

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  gBenchmark->Show("WWTemplate");       



  TGraphAsymmErrors *Ele_ChargedIso03_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_NeutralHadronIso03_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_GammaIso03_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_ChargedIso04_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_NeutralHadronIso04_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_GammaIso04_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_NeutralHadronIso007_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_HoverE_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_HcalDepth1OverEcal_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_HcalDepth2OverEcal_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_TrkIso03_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_EMIso03_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_HadIso03_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_TrkIso04_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_EMIso04_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_HadIso04_Vs_NPU_Graph = 0;
  TGraphAsymmErrors *Ele_Rho_Vs_NPU_Graph = 0;
  if(isMC) {
    Ele_ChargedIso03_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_ChargedIso03_VS_NPU,("Ele_ChargedIso03_VS_NPU"+label+"_Graph").c_str(),"ChargedIso03 [GeV]");
    Ele_NeutralHadronIso03_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_NeutralHadronIso03_VS_NPU,("Ele_NeutralHadronIso03_VS_NPU"+label+"_Graph").c_str(),"NeutralHadronIso03 [GeV]");
    Ele_GammaIso03_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_GammaIso03_VS_NPU,("Ele_GammaIso03_VS_NPU"+label+"_Graph").c_str(),"GammaIso03 [GeV]");
    Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_GammaIsoVetoEtaStrip03_VS_NPU,("Ele_GammaIsoVetoEtaStrip03_VS_NPU"+label+"_Graph").c_str(),"GammaIsoVetoEtaStrip03 [GeV]");
    Ele_ChargedIso04_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_ChargedIso04_VS_NPU,("Ele_ChargedIso04_VS_NPU"+label+"_Graph").c_str(),"ChargedIso04 [GeV]");
    Ele_NeutralHadronIso04_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_NeutralHadronIso04_VS_NPU,("Ele_NeutralHadronIso04_VS_NPU"+label+"_Graph").c_str(),"NeutralHadronIso04 [GeV]");
    Ele_GammaIso04_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_GammaIso04_VS_NPU,("Ele_GammaIso04_VS_NPU"+label+"_Graph").c_str(),"GammaIso04 [GeV]");
    Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_GammaIsoVetoEtaStrip04_VS_NPU,("Ele_GammaIsoVetoEtaStrip04_VS_NPU"+label+"_Graph").c_str(),"GammaIsoVetoEtaStrip04 [GeV]");
    Ele_NeutralHadronIso007_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_NeutralHadronIso007_VS_NPU,("Ele_NeutralHadronIso007_VS_NPU"+label+"_Graph").c_str(),"NeutralHadronIso007 [GeV]");
    Ele_HoverE_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_HoverE_VS_NPU,("Ele_HoverE_VS_NPU"+label+"_Graph").c_str(),"HoverE");
    Ele_HcalDepth1OverEcal_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_HcalDepth1OverEcal_VS_NPU,("Ele_HcalDepth1OverEcal_VS_NPU"+label+"_Graph").c_str(),"HcalDepth1OverEcal");
    Ele_HcalDepth2OverEcal_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_HcalDepth2OverEcal_VS_NPU,("Ele_HcalDepth2OverEcal_VS_NPU"+label+"_Graph").c_str(),"HcalDepth2OverEcal");
    Ele_TrkIso03_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_TrkIso03_VS_NPU,("Ele_TrkIso03_VS_NPU"+label+"_Graph").c_str(),"TrkIso03 [GeV]");
    Ele_EMIso03_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_EMIso03_VS_NPU,("Ele_EMIso03_VS_NPU"+label+"_Graph").c_str(),"EMIso03 [GeV]");
    Ele_HadIso03_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_HadIso03_VS_NPU,("Ele_HadIso03_VS_NPU"+label+"_Graph").c_str(),"HadIso03 [GeV]");
    Ele_TrkIso04_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_TrkIso04_VS_NPU,("Ele_TrkIso04_VS_NPU"+label+"_Graph").c_str(),"TrkIso04 [GeV]");
    Ele_EMIso04_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_EMIso04_VS_NPU,("Ele_EMIso04_VS_NPU"+label+"_Graph").c_str(),"EMIso04 [GeV]");
    Ele_HadIso04_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_HadIso04_VS_NPU,("Ele_HadIso04_VS_NPU"+label+"_Graph").c_str(),"HadIso04 [GeV]");
    Ele_Rho_Vs_NPU_Graph = MakeMeanVsNPUGraph(Ele_Rho_VS_NPU,("Ele_Rho_VS_NPU"+label+"_Graph").c_str(),"#rho [GeV]");
  }

  TGraphAsymmErrors *Ele_ChargedIso03_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_ChargedIso03_VS_NVtx,("Ele_ChargedIso03_VS_NVtx"+label+"_Graph").c_str(),"ChargedIso03 [GeV]");
  TGraphAsymmErrors *Ele_NeutralHadronIso03_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_NeutralHadronIso03_VS_NVtx,("Ele_NeutralHadronIso03_VS_NVtx"+label+"_Graph").c_str(),"NeutralHadronIso03 [GeV]");
  TGraphAsymmErrors *Ele_GammaIso03_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_GammaIso03_VS_NVtx,("Ele_GammaIso03_VS_NVtx"+label+"_Graph").c_str(),"GammaIso03 [GeV]");
  TGraphAsymmErrors *Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_GammaIsoVetoEtaStrip03_VS_NVtx,("Ele_GammaIsoVetoEtaStrip03_VS_NVtx"+label+"_Graph").c_str(),"GammaIsoVetoEtaStrip03 [GeV]");
  TGraphAsymmErrors *Ele_ChargedIso04_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_ChargedIso04_VS_NVtx,("Ele_ChargedIso04_VS_NVtx"+label+"_Graph").c_str(),"ChargedIso04 [GeV]");
  TGraphAsymmErrors *Ele_NeutralHadronIso04_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_NeutralHadronIso04_VS_NVtx,("Ele_NeutralHadronIso04_VS_NVtx"+label+"_Graph").c_str(),"NeutralHadronIso04 [GeV]");
  TGraphAsymmErrors *Ele_GammaIso04_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_GammaIso04_VS_NVtx,("Ele_GammaIso04_VS_NVtx"+label+"_Graph").c_str(),"GammaIso04 [GeV]");
  TGraphAsymmErrors *Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_GammaIsoVetoEtaStrip04_VS_NVtx,("Ele_GammaIsoVetoEtaStrip04_VS_NVtx"+label+"_Graph").c_str(),"GammaIsoVetoEtaStrip04 [GeV]");
  TGraphAsymmErrors *Ele_NeutralHadronIso007_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_NeutralHadronIso007_VS_NVtx,("Ele_NeutralHadronIso007_VS_NVtx"+label+"_Graph").c_str(),"NeutralHadronIso007 [GeV]");
  TGraphAsymmErrors *Ele_HoverE_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_HoverE_VS_NVtx,("Ele_HoverE_VS_NVtx"+label+"_Graph").c_str(),"HoverE");
  TGraphAsymmErrors *Ele_HcalDepth1OverEcal_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_HcalDepth1OverEcal_VS_NVtx,("Ele_HcalDepth1OverEcal_VS_NVtx"+label+"_Graph").c_str(),"HcalDepth1OverEcal");
  TGraphAsymmErrors *Ele_HcalDepth2OverEcal_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_HcalDepth2OverEcal_VS_NVtx,("Ele_HcalDepth2OverEcal_VS_NVtx"+label+"_Graph").c_str(),"HcalDepth2OverEcal");
  TGraphAsymmErrors *Ele_TrkIso03_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_TrkIso03_VS_NVtx,("Ele_TrkIso03_VS_NVtx"+label+"_Graph").c_str(),"TrkIso03 [GeV]");
  TGraphAsymmErrors *Ele_EMIso03_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_EMIso03_VS_NVtx,("Ele_EMIso03_VS_NVtx"+label+"_Graph").c_str(),"EMIso03 [GeV]");
  TGraphAsymmErrors *Ele_HadIso03_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_HadIso03_VS_NVtx,("Ele_HadIso03_VS_NVtx"+label+"_Graph").c_str(),"HadIso03 [GeV]");
  TGraphAsymmErrors *Ele_TrkIso04_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_TrkIso04_VS_NVtx,("Ele_TrkIso04_VS_NVtx"+label+"_Graph").c_str(),"TrkIso04 [GeV]");
  TGraphAsymmErrors *Ele_EMIso04_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_EMIso04_VS_NVtx,("Ele_EMIso04_VS_NVtx"+label+"_Graph").c_str(),"EMIso04 [GeV]");
  TGraphAsymmErrors *Ele_HadIso04_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_HadIso04_VS_NVtx,("Ele_HadIso04_VS_NVtx"+label+"_Graph").c_str(),"HadIso04 [GeV]");
  TGraphAsymmErrors *Ele_Rho_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Ele_Rho_VS_NVtx,("Ele_Rho_VS_NVtx"+label+"_Graph").c_str(),"#rho [GeV]");
 
  TGraphAsymmErrors* Ele_ChargedIso_DR0p0To0p1_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_ChargedIso_DR0p0To0p1_VS_NVtx, ("Ele_ChargedIso_DR0p0To0p1_VS_NVtx"+label+"_Graph").c_str(), "ChargedIso_DR0p0To0p1 [GeV]");
  TGraphAsymmErrors* Ele_ChargedIso_DR0p1To0p2_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_ChargedIso_DR0p1To0p2_VS_NVtx, ("Ele_ChargedIso_DR0p1To0p2_VS_NVtx"+label+"_Graph").c_str(), "ChargedIso_DR0p1To0p2 [GeV]");
  TGraphAsymmErrors* Ele_ChargedIso_DR0p2To0p3_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_ChargedIso_DR0p2To0p3_VS_NVtx, ("Ele_ChargedIso_DR0p2To0p3_VS_NVtx"+label+"_Graph").c_str(), "ChargedIso_DR0p2To0p3 [GeV]");
  TGraphAsymmErrors* Ele_ChargedIso_DR0p3To0p4_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_ChargedIso_DR0p3To0p4_VS_NVtx, ("Ele_ChargedIso_DR0p3To0p4_VS_NVtx"+label+"_Graph").c_str(), "ChargedIso_DR0p3To0p4 [GeV]");
  TGraphAsymmErrors* Ele_ChargedIso_DR0p4To0p5_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_ChargedIso_DR0p4To0p5_VS_NVtx, ("Ele_ChargedIso_DR0p4To0p5_VS_NVtx"+label+"_Graph").c_str(), "ChargedIso_DR0p4To0p5 [GeV]");
  TGraphAsymmErrors* Ele_ChargedIso_DR0p5To0p7_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_ChargedIso_DR0p5To0p7_VS_NVtx, ("Ele_ChargedIso_DR0p5To0p7_VS_NVtx"+label+"_Graph").c_str(), "ChargedIso_DR0p5To0p7 [GeV]");
  TGraphAsymmErrors* Ele_ChargedIso_DR0p7To1p0_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_ChargedIso_DR0p7To1p0_VS_NVtx, ("Ele_ChargedIso_DR0p7To1p0_VS_NVtx"+label+"_Graph").c_str(), "ChargedIso_DR0p7To1p0 [GeV]");
  TGraphAsymmErrors* Ele_GammaIso_DR0p0To0p1_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_GammaIso_DR0p0To0p1_VS_NVtx, ("Ele_GammaIso_DR0p0To0p1_VS_NVtx"+label+"_Graph").c_str(), "GammaIso_DR0p0To0p1 [GeV]");
  TGraphAsymmErrors* Ele_GammaIso_DR0p1To0p2_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_GammaIso_DR0p1To0p2_VS_NVtx, ("Ele_GammaIso_DR0p1To0p2_VS_NVtx"+label+"_Graph").c_str(), "GammaIso_DR0p1To0p2 [GeV]");
  TGraphAsymmErrors* Ele_GammaIso_DR0p2To0p3_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_GammaIso_DR0p2To0p3_VS_NVtx, ("Ele_GammaIso_DR0p2To0p3_VS_NVtx"+label+"_Graph").c_str(), "GammaIso_DR0p2To0p3 [GeV]");
  TGraphAsymmErrors* Ele_GammaIso_DR0p3To0p4_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_GammaIso_DR0p3To0p4_VS_NVtx, ("Ele_GammaIso_DR0p3To0p4_VS_NVtx"+label+"_Graph").c_str(), "GammaIso_DR0p3To0p4 [GeV]");
  TGraphAsymmErrors* Ele_GammaIso_DR0p4To0p5_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_GammaIso_DR0p4To0p5_VS_NVtx, ("Ele_GammaIso_DR0p4To0p5_VS_NVtx"+label+"_Graph").c_str(), "GammaIso_DR0p4To0p5 [GeV]");
  TGraphAsymmErrors* Ele_GammaIso_DR0p5To0p7_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_GammaIso_DR0p5To0p7_VS_NVtx, ("Ele_GammaIso_DR0p5To0p7_VS_NVtx"+label+"_Graph").c_str(), "GammaIso_DR0p5To0p7 [GeV]");
  TGraphAsymmErrors* Ele_GammaIso_DR0p7To1p0_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_GammaIso_DR0p7To1p0_VS_NVtx, ("Ele_GammaIso_DR0p7To1p0_VS_NVtx"+label+"_Graph").c_str(), "GammaIso_DR0p7To1p0 [GeV]");
  TGraphAsymmErrors* Ele_NeutralHadronIso_DR0p0To0p1_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_NeutralHadronIso_DR0p0To0p1_VS_NVtx, ("Ele_NeutralHadronIso_DR0p0To0p1_VS_NVtx"+label+"_Graph").c_str(), "NeutralHadronIso_DR0p0To0p1 [GeV]");
  TGraphAsymmErrors* Ele_NeutralHadronIso_DR0p1To0p2_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_NeutralHadronIso_DR0p1To0p2_VS_NVtx, ("Ele_NeutralHadronIso_DR0p1To0p2_VS_NVtx"+label+"_Graph").c_str(), "NeutralHadronIso_DR0p1To0p2 [GeV]");
  TGraphAsymmErrors* Ele_NeutralHadronIso_DR0p2To0p3_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_NeutralHadronIso_DR0p2To0p3_VS_NVtx, ("Ele_NeutralHadronIso_DR0p2To0p3_VS_NVtx"+label+"_Graph").c_str(), "NeutralHadronIso_DR0p2To0p3 [GeV]");
  TGraphAsymmErrors* Ele_NeutralHadronIso_DR0p3To0p4_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_NeutralHadronIso_DR0p3To0p4_VS_NVtx, ("Ele_NeutralHadronIso_DR0p3To0p4_VS_NVtx"+label+"_Graph").c_str(), "NeutralHadronIso_DR0p3To0p4 [GeV]");
  TGraphAsymmErrors* Ele_NeutralHadronIso_DR0p4To0p5_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_NeutralHadronIso_DR0p4To0p5_VS_NVtx, ("Ele_NeutralHadronIso_DR0p4To0p5_VS_NVtx"+label+"_Graph").c_str(), "NeutralHadronIso_DR0p4To0p5 [GeV]");
  TGraphAsymmErrors* Ele_NeutralHadronIso_DR0p5To0p7_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_NeutralHadronIso_DR0p5To0p7_VS_NVtx, ("Ele_NeutralHadronIso_DR0p5To0p7_VS_NVtx"+label+"_Graph").c_str(), "NeutralHadronIso_DR0p5To0p7 [GeV]");
  TGraphAsymmErrors* Ele_NeutralHadronIso_DR0p7To1p0_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_NeutralHadronIso_DR0p7To1p0_VS_NVtx, ("Ele_NeutralHadronIso_DR0p7To1p0_VS_NVtx"+label+"_Graph").c_str(), "NeutralHadronIso_DR0p7To1p0 [GeV]");

  TGraphAsymmErrors* Ele_ChargedIso_DR0p0To0p3_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_ChargedIso_DR0p0To0p3_VS_NVtx, ("Ele_ChargedIso_DR0p0To0p3_VS_NVtx"+label+"_Graph").c_str(), "ChargedIso_DR0p0To0p3 [GeV]");
  TGraphAsymmErrors* Ele_ChargedIso_DR0p0To0p4_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_ChargedIso_DR0p0To0p4_VS_NVtx, ("Ele_ChargedIso_DR0p0To0p4_VS_NVtx"+label+"_Graph").c_str(), "ChargedIso_DR0p0To0p4 [GeV]");
  TGraphAsymmErrors* Ele_GammaIso_DR0p0To0p3_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_GammaIso_DR0p0To0p3_VS_NVtx, ("Ele_GammaIso_DR0p0To0p3_VS_NVtx"+label+"_Graph").c_str(), "GammaIso_DR0p0To0p3 [GeV]");
  TGraphAsymmErrors* Ele_GammaIso_DR0p0To0p4_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_GammaIso_DR0p0To0p4_VS_NVtx, ("Ele_GammaIso_DR0p0To0p4_VS_NVtx"+label+"_Graph").c_str(), "GammaIso_DR0p0To0p4 [GeV]");
  TGraphAsymmErrors* Ele_NeutralHadronIso_DR0p0To0p3_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_NeutralHadronIso_DR0p0To0p3_VS_NVtx, ("Ele_NeutralHadronIso_DR0p0To0p3_VS_NVtx"+label+"_Graph").c_str(), "NeutralHadronIso_DR0p0To0p3 [GeV]");
  TGraphAsymmErrors* Ele_NeutralHadronIso_DR0p0To0p4_VS_NVtx_Graph = MakeMeanVsNPUGraph( Ele_NeutralHadronIso_DR0p0To0p4_VS_NVtx, ("Ele_NeutralHadronIso_DR0p0To0p4_VS_NVtx"+label+"_Graph").c_str(), "NeutralHadronIso_DR0p0To0p4 [GeV]");


  //*****************************************************************************************
  //Save Histograms in file
  //*****************************************************************************************
  TFile *file = new TFile("ElectronEffectiveArea.root", "UPDATE");
  file->cd();
  if(isMC) {
    file->WriteTObject(Ele_ChargedIso03_Vs_NPU_Graph , Ele_ChargedIso03_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_NeutralHadronIso03_Vs_NPU_Graph , Ele_NeutralHadronIso03_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_GammaIso03_Vs_NPU_Graph , Ele_GammaIso03_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Graph , Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_ChargedIso04_Vs_NPU_Graph , Ele_ChargedIso04_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_NeutralHadronIso04_Vs_NPU_Graph , Ele_NeutralHadronIso04_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_GammaIso04_Vs_NPU_Graph , Ele_GammaIso04_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Graph , Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_NeutralHadronIso007_Vs_NPU_Graph , Ele_NeutralHadronIso007_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_HoverE_Vs_NPU_Graph , Ele_HoverE_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_HcalDepth1OverEcal_Vs_NPU_Graph , Ele_HcalDepth1OverEcal_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_HcalDepth2OverEcal_Vs_NPU_Graph , Ele_HcalDepth2OverEcal_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_TrkIso03_Vs_NPU_Graph , Ele_TrkIso03_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_EMIso03_Vs_NPU_Graph , Ele_EMIso03_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_HadIso03_Vs_NPU_Graph , Ele_HadIso03_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_TrkIso04_Vs_NPU_Graph , Ele_TrkIso04_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_EMIso04_Vs_NPU_Graph , Ele_EMIso04_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_HadIso04_Vs_NPU_Graph , Ele_HadIso04_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Ele_Rho_Vs_NPU_Graph , Ele_Rho_Vs_NPU_Graph->GetName(), "WriteDelete");  
  }

  file->WriteTObject(Ele_ChargedIso03_Vs_NVtx_Graph , Ele_ChargedIso03_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_NeutralHadronIso03_Vs_NVtx_Graph , Ele_NeutralHadronIso03_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIso03_Vs_NVtx_Graph , Ele_GammaIso03_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Graph , Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_ChargedIso04_Vs_NVtx_Graph , Ele_ChargedIso04_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_NeutralHadronIso04_Vs_NVtx_Graph , Ele_NeutralHadronIso04_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIso04_Vs_NVtx_Graph , Ele_GammaIso04_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Graph , Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_NeutralHadronIso007_Vs_NVtx_Graph , Ele_NeutralHadronIso007_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_HoverE_Vs_NVtx_Graph , Ele_HoverE_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_HcalDepth1OverEcal_Vs_NVtx_Graph , Ele_HcalDepth1OverEcal_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_HcalDepth2OverEcal_Vs_NVtx_Graph , Ele_HcalDepth2OverEcal_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_TrkIso03_Vs_NVtx_Graph , Ele_TrkIso03_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_EMIso03_Vs_NVtx_Graph , Ele_EMIso03_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_HadIso03_Vs_NVtx_Graph , Ele_HadIso03_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_TrkIso04_Vs_NVtx_Graph , Ele_TrkIso04_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_EMIso04_Vs_NVtx_Graph , Ele_EMIso04_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_HadIso04_Vs_NVtx_Graph , Ele_HadIso04_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_ChargedIso_DR0p0To0p1_VS_NVtx_Graph , Ele_ChargedIso_DR0p0To0p1_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_ChargedIso_DR0p1To0p2_VS_NVtx_Graph , Ele_ChargedIso_DR0p1To0p2_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_ChargedIso_DR0p2To0p3_VS_NVtx_Graph , Ele_ChargedIso_DR0p2To0p3_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_ChargedIso_DR0p3To0p4_VS_NVtx_Graph , Ele_ChargedIso_DR0p3To0p4_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_ChargedIso_DR0p4To0p5_VS_NVtx_Graph , Ele_ChargedIso_DR0p4To0p5_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_ChargedIso_DR0p5To0p7_VS_NVtx_Graph , Ele_ChargedIso_DR0p5To0p7_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_ChargedIso_DR0p7To1p0_VS_NVtx_Graph , Ele_ChargedIso_DR0p7To1p0_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIso_DR0p0To0p1_VS_NVtx_Graph , Ele_GammaIso_DR0p0To0p1_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIso_DR0p1To0p2_VS_NVtx_Graph , Ele_GammaIso_DR0p1To0p2_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIso_DR0p2To0p3_VS_NVtx_Graph , Ele_GammaIso_DR0p2To0p3_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIso_DR0p3To0p4_VS_NVtx_Graph , Ele_GammaIso_DR0p3To0p4_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIso_DR0p4To0p5_VS_NVtx_Graph , Ele_GammaIso_DR0p4To0p5_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIso_DR0p5To0p7_VS_NVtx_Graph , Ele_GammaIso_DR0p5To0p7_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIso_DR0p7To1p0_VS_NVtx_Graph , Ele_GammaIso_DR0p7To1p0_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_NeutralHadronIso_DR0p0To0p1_VS_NVtx_Graph , Ele_NeutralHadronIso_DR0p0To0p1_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_NeutralHadronIso_DR0p1To0p2_VS_NVtx_Graph , Ele_NeutralHadronIso_DR0p1To0p2_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_NeutralHadronIso_DR0p2To0p3_VS_NVtx_Graph , Ele_NeutralHadronIso_DR0p2To0p3_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_NeutralHadronIso_DR0p3To0p4_VS_NVtx_Graph , Ele_NeutralHadronIso_DR0p3To0p4_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_NeutralHadronIso_DR0p4To0p5_VS_NVtx_Graph , Ele_NeutralHadronIso_DR0p4To0p5_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_NeutralHadronIso_DR0p5To0p7_VS_NVtx_Graph , Ele_NeutralHadronIso_DR0p5To0p7_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_NeutralHadronIso_DR0p7To1p0_VS_NVtx_Graph , Ele_NeutralHadronIso_DR0p7To1p0_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_Rho_Vs_NVtx_Graph , Ele_Rho_Vs_NVtx_Graph->GetName(), "WriteDelete");   

  file->WriteTObject(Ele_ChargedIso_DR0p0To0p3_VS_NVtx_Graph , Ele_ChargedIso_DR0p0To0p3_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_ChargedIso_DR0p0To0p4_VS_NVtx_Graph , Ele_ChargedIso_DR0p0To0p4_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIso_DR0p0To0p3_VS_NVtx_Graph , Ele_GammaIso_DR0p0To0p3_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_GammaIso_DR0p0To0p4_VS_NVtx_Graph , Ele_GammaIso_DR0p0To0p4_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_NeutralHadronIso_DR0p0To0p3_VS_NVtx_Graph , Ele_NeutralHadronIso_DR0p0To0p3_VS_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Ele_NeutralHadronIso_DR0p0To0p4_VS_NVtx_Graph , Ele_NeutralHadronIso_DR0p0To0p4_VS_NVtx_Graph->GetName(), "WriteDelete");  

  file->Close();
  delete file;

  gBenchmark->Show("WWTemplate");       
} 



void ComputeEffectiveArea(string Label, Bool_t isMC = kTRUE, Int_t EtaBin = -1) {

  string label = "";
  if (Label != "") label = "_" + Label;
  
  if (EtaBin == 0) {
    label = label + "_Eta0To1";
  }
  if (EtaBin == 1) {
    label = label + "_Eta1To1p479";
  }
  if (EtaBin == 2) {
    label = label + "_Eta1p479To2";
  }
  if (EtaBin == 3) {
    label = label + "_Eta2To2p2";
  }
  if (EtaBin == 4) {
    label = label + "_Eta2p2To2p3";
  }
  if (EtaBin == 5) {
    label = label + "_Eta2p3To2p4";
  }
  if (EtaBin == 6) {
    label = label + "_Eta2p4To2p5";
  }

  TFile *file = new TFile("ElectronEffectiveArea.root", "READ");
    
  TGraphAsymmErrors *Ele_ChargedIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_NeutralHadronIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_GammaIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_ChargedIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_NeutralHadronIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_GammaIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_NeutralHadronIso007_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_HoverE_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_HcalDepth1OverEcal_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_HcalDepth2OverEcal_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_TrkIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_EMIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_HadIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_TrkIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_EMIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_HadIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_Rho_Vs_NPU_Graph;
  if (isMC) {
    Ele_ChargedIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso03_VS_NPU"+label+"_Graph").c_str());
    Ele_NeutralHadronIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso03_VS_NPU"+label+"_Graph").c_str());
    Ele_GammaIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso03_VS_NPU"+label+"_Graph").c_str());
    Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIsoVetoEtaStrip03_VS_NPU"+label+"_Graph").c_str());
    Ele_ChargedIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso04_VS_NPU"+label+"_Graph").c_str());
    Ele_NeutralHadronIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso04_VS_NPU"+label+"_Graph").c_str());
    Ele_GammaIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso04_VS_NPU"+label+"_Graph").c_str());
    Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIsoVetoEtaStrip04_VS_NPU"+label+"_Graph").c_str());
    Ele_NeutralHadronIso007_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso007_VS_NPU"+label+"_Graph").c_str());
    Ele_HoverE_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HoverE_VS_NPU"+label+"_Graph").c_str());
    Ele_HcalDepth1OverEcal_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HcalDepth1OverEcal_VS_NPU"+label+"_Graph").c_str());
    Ele_HcalDepth2OverEcal_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HcalDepth2OverEcal_VS_NPU"+label+"_Graph").c_str());
    Ele_TrkIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_TrkIso03_VS_NPU"+label+"_Graph").c_str());
    Ele_EMIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_EMIso03_VS_NPU"+label+"_Graph").c_str());
    Ele_HadIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HadIso03_VS_NPU"+label+"_Graph").c_str());
    Ele_TrkIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_TrkIso04_VS_NPU"+label+"_Graph").c_str());
    Ele_EMIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_EMIso04_VS_NPU"+label+"_Graph").c_str());
    Ele_HadIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HadIso04_VS_NPU"+label+"_Graph").c_str());
    Ele_Rho_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_Rho_VS_NPU"+label+"_Graph").c_str());
  }

  TGraphAsymmErrors *Ele_ChargedIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIsoVetoEtaStrip03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_ChargedIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIsoVetoEtaStrip04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso007_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso007_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_HoverE_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HoverE_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_HcalDepth1OverEcal_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HcalDepth1OverEcal_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_HcalDepth2OverEcal_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HcalDepth2OverEcal_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_TrkIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_TrkIso03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_EMIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_EMIso03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_HadIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HadIso03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_TrkIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_TrkIso04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_EMIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_EMIso04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_HadIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HadIso04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso_DR0p0To0p1_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso_DR0p1To0p2_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso_DR0p2To0p3_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso_DR0p3To0p4_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso_DR0p4To0p5_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso_DR0p5To0p7_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso_DR0p7To1p0_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso_DR0p0To0p1_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso_DR0p0To0p1_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso_DR0p1To0p2_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso_DR0p1To0p2_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso_DR0p2To0p3_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso_DR0p2To0p3_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso_DR0p3To0p4_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso_DR0p3To0p4_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso_DR0p4To0p5_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso_DR0p4To0p5_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso_DR0p5To0p7_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso_DR0p5To0p7_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso_DR0p7To1p0_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso_DR0p7To1p0_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso_DR0p0To0p1_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso_DR0p1To0p2_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso_DR0p2To0p3_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso_DR0p3To0p4_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso_DR0p4To0p5_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso_DR0p5To0p7_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso_DR0p7To1p0_VS_NVtx"+label+"_Graph").c_str());


  TGraphAsymmErrors *Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso_DR0p0To0p3_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso_DR0p0To0p4_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso_DR0p0To0p3_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso_DR0p0To0p3_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso_DR0p0To0p4_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso_DR0p0To0p4_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso_DR0p0To0p3_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso_DR0p0To0p4_VS_NVtx"+label+"_Graph").c_str());

  TGraphAsymmErrors *Ele_Rho_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_Rho_VS_NVtx"+label+"_Graph").c_str());
   
  assert(Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_Graph);
  assert(Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_Graph);
  assert(Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_Graph);
  assert(Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_Graph);
  assert(Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_Graph);
  assert(Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_Graph);
  assert(Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_Graph);
  assert(Ele_GammaIso_DR0p0To0p1_Vs_NVtx_Graph);
  assert(Ele_GammaIso_DR0p1To0p2_Vs_NVtx_Graph);
  assert(Ele_GammaIso_DR0p2To0p3_Vs_NVtx_Graph);
  assert(Ele_GammaIso_DR0p3To0p4_Vs_NVtx_Graph);
  assert(Ele_GammaIso_DR0p4To0p5_Vs_NVtx_Graph);
  assert(Ele_GammaIso_DR0p5To0p7_Vs_NVtx_Graph);
  assert(Ele_GammaIso_DR0p7To1p0_Vs_NVtx_Graph);
  assert(Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_Graph);
  assert(Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_Graph);
  assert(Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_Graph);
  assert(Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_Graph);
  assert(Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_Graph);
  assert(Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_Graph);
  assert(Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_Graph);
  assert(Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_Graph);
  assert(Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_Graph);
  assert(Ele_GammaIso_DR0p0To0p3_Vs_NVtx_Graph);
  assert(Ele_GammaIso_DR0p0To0p4_Vs_NVtx_Graph);
  assert(Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_Graph);
  assert(Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_Graph);
  assert(Ele_Rho_Vs_NVtx_Graph);

  Double_t Ele_ChargedIso03_Vs_NPU_EffArea = 0;
  Double_t Ele_ChargedIso03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_NeutralHadronIso03_Vs_NPU_EffArea = 0;
  Double_t Ele_NeutralHadronIso03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_GammaIso03_Vs_NPU_EffArea = 0;
  Double_t Ele_GammaIso03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffArea = 0;
  Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_ChargedIso04_Vs_NPU_EffArea = 0;
  Double_t Ele_ChargedIso04_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_NeutralHadronIso04_Vs_NPU_EffArea = 0;
  Double_t Ele_NeutralHadronIso04_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_GammaIso04_Vs_NPU_EffArea = 0;
  Double_t Ele_GammaIso04_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffArea = 0;
  Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_NeutralHadronIso007_Vs_NPU_EffArea = 0;
  Double_t Ele_NeutralHadronIso007_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_HoverE_Vs_NPU_EffArea = 0;
  Double_t Ele_HoverE_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_HcalDepth1OverEcal_Vs_NPU_EffArea = 0;
  Double_t Ele_HcalDepth1OverEcal_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_HcalDepth2OverEcal_Vs_NPU_EffArea = 0;
  Double_t Ele_HcalDepth2OverEcal_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_TrkIso03_Vs_NPU_EffArea = 0;
  Double_t Ele_TrkIso03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_EMIso03_Vs_NPU_EffArea = 0;
  Double_t Ele_EMIso03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_HadIso03_Vs_NPU_EffArea = 0;
  Double_t Ele_HadIso03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_TrkIso04_Vs_NPU_EffArea = 0;
  Double_t Ele_TrkIso04_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_EMIso04_Vs_NPU_EffArea = 0;
  Double_t Ele_EMIso04_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_HadIso04_Vs_NPU_EffArea = 0;
  Double_t Ele_HadIso04_Vs_NPU_EffAreaErr = 0;

  TF1 *f1= 0;

  if (isMC) {
    //*****************************************************************************************************
    //*****************************************************************************************************
    //Vs NPU
    //*****************************************************************************************************
    //*****************************************************************************************************

    //----------------------
    // Rho
    //----------------------
    f1 = new TF1(("Ele_Rho_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_Rho_Vs_NPU_Graph->Fit(("Ele_Rho_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_Rho_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_Rho_Vs_NPU_SlopeErr = f1->GetParError(1);

    //----------------------
    // Energy Observables
    //----------------------
    f1 = new TF1(("Ele_ChargedIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_ChargedIso03_Vs_NPU_Graph->Fit(("Ele_ChargedIso03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_ChargedIso03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_ChargedIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_ChargedIso03_Vs_NPU_EffArea = Ele_ChargedIso03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_ChargedIso03_Vs_NPU_EffAreaErr = Ele_ChargedIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_ChargedIso03_Vs_NPU_SlopeErr/Ele_ChargedIso03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_NeutralHadronIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_NeutralHadronIso03_Vs_NPU_Graph->Fit(("Ele_NeutralHadronIso03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_NeutralHadronIso03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_NeutralHadronIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_NeutralHadronIso03_Vs_NPU_EffArea = Ele_NeutralHadronIso03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_NeutralHadronIso03_Vs_NPU_EffAreaErr = Ele_NeutralHadronIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_NeutralHadronIso03_Vs_NPU_SlopeErr/Ele_NeutralHadronIso03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_GammaIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_GammaIso03_Vs_NPU_Graph->Fit(("Ele_GammaIso03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_GammaIso03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_GammaIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_GammaIso03_Vs_NPU_EffArea = Ele_GammaIso03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_GammaIso03_Vs_NPU_EffAreaErr = Ele_GammaIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_GammaIso03_Vs_NPU_SlopeErr/Ele_GammaIso03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_GammaIsoVetoEtaStrip03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Graph->Fit(("Ele_GammaIsoVetoEtaStrip03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffArea = Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffAreaErr = Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_GammaIsoVetoEtaStrip03_Vs_NPU_SlopeErr/Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_ChargedIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_ChargedIso04_Vs_NPU_Graph->Fit(("Ele_ChargedIso04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_ChargedIso04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_ChargedIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_ChargedIso04_Vs_NPU_EffArea = Ele_ChargedIso04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_ChargedIso04_Vs_NPU_EffAreaErr = Ele_ChargedIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_ChargedIso04_Vs_NPU_SlopeErr/Ele_ChargedIso04_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_NeutralHadronIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_NeutralHadronIso04_Vs_NPU_Graph->Fit(("Ele_NeutralHadronIso04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_NeutralHadronIso04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_NeutralHadronIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_NeutralHadronIso04_Vs_NPU_EffArea = Ele_NeutralHadronIso04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_NeutralHadronIso04_Vs_NPU_EffAreaErr = Ele_NeutralHadronIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_NeutralHadronIso04_Vs_NPU_SlopeErr/Ele_NeutralHadronIso04_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_GammaIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_GammaIso04_Vs_NPU_Graph->Fit(("Ele_GammaIso04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_GammaIso04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_GammaIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_GammaIso04_Vs_NPU_EffArea = Ele_GammaIso04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_GammaIso04_Vs_NPU_EffAreaErr = Ele_GammaIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_GammaIso04_Vs_NPU_SlopeErr/Ele_GammaIso04_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_GammaIsoVetoEtaStrip04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Graph->Fit(("Ele_GammaIsoVetoEtaStrip04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffArea = Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffAreaErr = Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_GammaIsoVetoEtaStrip04_Vs_NPU_SlopeErr/Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_NeutralHadronIso007_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_NeutralHadronIso007_Vs_NPU_Graph->Fit(("Ele_NeutralHadronIso007_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_NeutralHadronIso007_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_NeutralHadronIso007_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_NeutralHadronIso007_Vs_NPU_EffArea = Ele_NeutralHadronIso007_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_NeutralHadronIso007_Vs_NPU_EffAreaErr = Ele_NeutralHadronIso007_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_NeutralHadronIso007_Vs_NPU_SlopeErr/Ele_NeutralHadronIso007_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_HoverE_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_HoverE_Vs_NPU_Graph->Fit(("Ele_HoverE_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_HoverE_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_HoverE_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_HoverE_Vs_NPU_EffArea = Ele_HoverE_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_HoverE_Vs_NPU_EffAreaErr = Ele_HoverE_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_HoverE_Vs_NPU_SlopeErr/Ele_HoverE_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_HcalDepth1OverEcal_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_HcalDepth1OverEcal_Vs_NPU_Graph->Fit(("Ele_HcalDepth1OverEcal_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_HcalDepth1OverEcal_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_HcalDepth1OverEcal_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_HcalDepth1OverEcal_Vs_NPU_EffArea = Ele_HcalDepth1OverEcal_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_HcalDepth1OverEcal_Vs_NPU_EffAreaErr = Ele_HcalDepth1OverEcal_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_HcalDepth1OverEcal_Vs_NPU_SlopeErr/Ele_HcalDepth1OverEcal_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_HcalDepth2OverEcal_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_HcalDepth2OverEcal_Vs_NPU_Graph->Fit(("Ele_HcalDepth2OverEcal_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_HcalDepth2OverEcal_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_HcalDepth2OverEcal_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_HcalDepth2OverEcal_Vs_NPU_EffArea = Ele_HcalDepth2OverEcal_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_HcalDepth2OverEcal_Vs_NPU_EffAreaErr = Ele_HcalDepth2OverEcal_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_HcalDepth2OverEcal_Vs_NPU_SlopeErr/Ele_HcalDepth2OverEcal_Vs_NPU_Slope,2));


    f1 = new TF1(("Ele_TrkIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_TrkIso03_Vs_NPU_Graph->Fit(("Ele_TrkIso03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_TrkIso03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_TrkIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_TrkIso03_Vs_NPU_EffArea = Ele_TrkIso03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_TrkIso03_Vs_NPU_EffAreaErr = Ele_TrkIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_TrkIso03_Vs_NPU_SlopeErr/Ele_TrkIso03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_EMIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_EMIso03_Vs_NPU_Graph->Fit(("Ele_EMIso03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_EMIso03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_EMIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_EMIso03_Vs_NPU_EffArea = Ele_EMIso03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_EMIso03_Vs_NPU_EffAreaErr = Ele_EMIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_EMIso03_Vs_NPU_SlopeErr/Ele_EMIso03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_HadIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_HadIso03_Vs_NPU_Graph->Fit(("Ele_HadIso03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_HadIso03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_HadIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_HadIso03_Vs_NPU_EffArea = Ele_HadIso03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_HadIso03_Vs_NPU_EffAreaErr = Ele_HadIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_HadIso03_Vs_NPU_SlopeErr/Ele_HadIso03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_TrkIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_TrkIso04_Vs_NPU_Graph->Fit(("Ele_TrkIso04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_TrkIso04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_TrkIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_TrkIso04_Vs_NPU_EffArea = Ele_TrkIso04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_TrkIso04_Vs_NPU_EffAreaErr = Ele_TrkIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_TrkIso04_Vs_NPU_SlopeErr/Ele_TrkIso04_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_EMIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_EMIso04_Vs_NPU_Graph->Fit(("Ele_EMIso04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_EMIso04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_EMIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_EMIso04_Vs_NPU_EffArea = Ele_EMIso04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_EMIso04_Vs_NPU_EffAreaErr = Ele_EMIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_EMIso04_Vs_NPU_SlopeErr/Ele_EMIso04_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_HadIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_HadIso04_Vs_NPU_Graph->Fit(("Ele_HadIso04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_HadIso04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_HadIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_HadIso04_Vs_NPU_EffArea = Ele_HadIso04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_HadIso04_Vs_NPU_EffAreaErr = Ele_HadIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_HadIso04_Vs_NPU_SlopeErr/Ele_HadIso04_Vs_NPU_Slope,2));


  }

  //*****************************************************************************************************
  //*****************************************************************************************************
  //Vs NVtx
  //*****************************************************************************************************
  //*****************************************************************************************************

  //----------------------
  // Rho
  //----------------------
  f1 = new TF1(("Ele_Rho_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_Rho_Vs_NVtx_Graph->Fit(("Ele_Rho_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_Rho_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_Rho_Vs_NVtx_SlopeErr = f1->GetParError(1);

  //----------------------
  // Energy Observables
  //----------------------
  f1 = new TF1(("Ele_ChargedIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_ChargedIso03_Vs_NVtx_Graph->Fit(("Ele_ChargedIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso03_Vs_NVtx_EffArea = Ele_ChargedIso03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso03_Vs_NVtx_EffAreaErr = Ele_ChargedIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso03_Vs_NVtx_SlopeErr/Ele_ChargedIso03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_NeutralHadronIso03_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso03_Vs_NVtx_EffArea = Ele_NeutralHadronIso03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso03_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso03_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_GammaIso03_Vs_NVtx_Graph->Fit(("Ele_GammaIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso03_Vs_NVtx_EffArea = Ele_GammaIso03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso03_Vs_NVtx_EffAreaErr = Ele_GammaIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso03_Vs_NVtx_SlopeErr/Ele_GammaIso03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIsoVetoEtaStrip03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Graph->Fit(("Ele_GammaIsoVetoEtaStrip03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffArea = Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffAreaErr = Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_SlopeErr/Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_ChargedIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_ChargedIso04_Vs_NVtx_Graph->Fit(("Ele_ChargedIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso04_Vs_NVtx_EffArea = Ele_ChargedIso04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso04_Vs_NVtx_EffAreaErr = Ele_ChargedIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso04_Vs_NVtx_SlopeErr/Ele_ChargedIso04_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_NeutralHadronIso04_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso04_Vs_NVtx_EffArea = Ele_NeutralHadronIso04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso04_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso04_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso04_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_GammaIso04_Vs_NVtx_Graph->Fit(("Ele_GammaIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso04_Vs_NVtx_EffArea = Ele_GammaIso04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso04_Vs_NVtx_EffAreaErr = Ele_GammaIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso04_Vs_NVtx_SlopeErr/Ele_GammaIso04_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIsoVetoEtaStrip04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Graph->Fit(("Ele_GammaIsoVetoEtaStrip04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffArea = Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffAreaErr = Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_SlopeErr/Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso007_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_NeutralHadronIso007_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso007_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso007_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso007_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso007_Vs_NVtx_EffArea = Ele_NeutralHadronIso007_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso007_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso007_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso007_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso007_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_HoverE_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_HoverE_Vs_NVtx_Graph->Fit(("Ele_HoverE_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_HoverE_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_HoverE_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_HoverE_Vs_NVtx_EffArea = Ele_HoverE_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_HoverE_Vs_NVtx_EffAreaErr = Ele_HoverE_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_HoverE_Vs_NVtx_SlopeErr/Ele_HoverE_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_HcalDepth1OverEcal_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_HcalDepth1OverEcal_Vs_NVtx_Graph->Fit(("Ele_HcalDepth1OverEcal_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_HcalDepth1OverEcal_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_HcalDepth1OverEcal_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_HcalDepth1OverEcal_Vs_NVtx_EffArea = Ele_HcalDepth1OverEcal_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_HcalDepth1OverEcal_Vs_NVtx_EffAreaErr = Ele_HcalDepth1OverEcal_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_HcalDepth1OverEcal_Vs_NVtx_SlopeErr/Ele_HcalDepth1OverEcal_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_HcalDepth2OverEcal_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_HcalDepth2OverEcal_Vs_NVtx_Graph->Fit(("Ele_HcalDepth2OverEcal_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_HcalDepth2OverEcal_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_HcalDepth2OverEcal_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_HcalDepth2OverEcal_Vs_NVtx_EffArea = Ele_HcalDepth2OverEcal_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_HcalDepth2OverEcal_Vs_NVtx_EffAreaErr = Ele_HcalDepth2OverEcal_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_HcalDepth2OverEcal_Vs_NVtx_SlopeErr/Ele_HcalDepth2OverEcal_Vs_NVtx_Slope,2));


  f1 = new TF1(("Ele_TrkIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_TrkIso03_Vs_NVtx_Graph->Fit(("Ele_TrkIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_TrkIso03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_TrkIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_TrkIso03_Vs_NVtx_EffArea = Ele_TrkIso03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_TrkIso03_Vs_NVtx_EffAreaErr = Ele_TrkIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_TrkIso03_Vs_NVtx_SlopeErr/Ele_TrkIso03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_EMIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_EMIso03_Vs_NVtx_Graph->Fit(("Ele_EMIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_EMIso03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_EMIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_EMIso03_Vs_NVtx_EffArea = Ele_EMIso03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_EMIso03_Vs_NVtx_EffAreaErr = Ele_EMIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_EMIso03_Vs_NVtx_SlopeErr/Ele_EMIso03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_HadIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_HadIso03_Vs_NVtx_Graph->Fit(("Ele_HadIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_HadIso03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_HadIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_HadIso03_Vs_NVtx_EffArea = Ele_HadIso03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_HadIso03_Vs_NVtx_EffAreaErr = Ele_HadIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_HadIso03_Vs_NVtx_SlopeErr/Ele_HadIso03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_TrkIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_TrkIso04_Vs_NVtx_Graph->Fit(("Ele_TrkIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_TrkIso04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_TrkIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_TrkIso04_Vs_NVtx_EffArea = Ele_TrkIso04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_TrkIso04_Vs_NVtx_EffAreaErr = Ele_TrkIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_TrkIso04_Vs_NVtx_SlopeErr/Ele_TrkIso04_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_EMIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_EMIso04_Vs_NVtx_Graph->Fit(("Ele_EMIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_EMIso04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_EMIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_EMIso04_Vs_NVtx_EffArea = Ele_EMIso04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_EMIso04_Vs_NVtx_EffAreaErr = Ele_EMIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_EMIso04_Vs_NVtx_SlopeErr/Ele_EMIso04_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_HadIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_HadIso04_Vs_NVtx_Graph->Fit(("Ele_HadIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_HadIso04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_HadIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_HadIso04_Vs_NVtx_EffArea = Ele_HadIso04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_HadIso04_Vs_NVtx_EffAreaErr = Ele_HadIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_HadIso04_Vs_NVtx_SlopeErr/Ele_HadIso04_Vs_NVtx_Slope,2));


  f1 = new TF1(("Ele_ChargedIso_DR0p0To0p1_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_Graph->Fit(("Ele_ChargedIso_DR0p0To0p1_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_EffArea = Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_EffAreaErr = Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_SlopeErr/Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_ChargedIso_DR0p1To0p2_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_Graph->Fit(("Ele_ChargedIso_DR0p1To0p2_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_EffArea = Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_EffAreaErr = Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_SlopeErr/Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_Slope,2));


  f1 = new TF1(("Ele_ChargedIso_DR0p2To0p3_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_Graph->Fit(("Ele_ChargedIso_DR0p2To0p3_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_EffArea = Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_EffAreaErr = Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_SlopeErr/Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_ChargedIso_DR0p3To0p4_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_Graph->Fit(("Ele_ChargedIso_DR0p3To0p4_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_EffArea = Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_EffAreaErr = Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_SlopeErr/Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_ChargedIso_DR0p4To0p5_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_Graph->Fit(("Ele_ChargedIso_DR0p4To0p5_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_EffArea = Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_EffAreaErr = Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_SlopeErr/Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_ChargedIso_DR0p5To0p7_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_Graph->Fit(("Ele_ChargedIso_DR0p5To0p7_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_EffArea = Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_EffAreaErr = Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_SlopeErr/Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_ChargedIso_DR0p7To1p0_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_Graph->Fit(("Ele_ChargedIso_DR0p7To1p0_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_EffArea = Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_EffAreaErr = Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_SlopeErr/Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso_DR0p0To0p1_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_GammaIso_DR0p0To0p1_Vs_NVtx_Graph->Fit(("Ele_GammaIso_DR0p0To0p1_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso_DR0p0To0p1_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso_DR0p0To0p1_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso_DR0p0To0p1_Vs_NVtx_EffArea = Ele_GammaIso_DR0p0To0p1_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso_DR0p0To0p1_Vs_NVtx_EffAreaErr = Ele_GammaIso_DR0p0To0p1_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso_DR0p0To0p1_Vs_NVtx_SlopeErr/Ele_GammaIso_DR0p0To0p1_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso_DR0p1To0p2_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_GammaIso_DR0p1To0p2_Vs_NVtx_Graph->Fit(("Ele_GammaIso_DR0p1To0p2_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso_DR0p1To0p2_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso_DR0p1To0p2_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso_DR0p1To0p2_Vs_NVtx_EffArea = Ele_GammaIso_DR0p1To0p2_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso_DR0p1To0p2_Vs_NVtx_EffAreaErr = Ele_GammaIso_DR0p1To0p2_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso_DR0p1To0p2_Vs_NVtx_SlopeErr/Ele_GammaIso_DR0p1To0p2_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso_DR0p2To0p3_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_GammaIso_DR0p2To0p3_Vs_NVtx_Graph->Fit(("Ele_GammaIso_DR0p2To0p3_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso_DR0p2To0p3_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso_DR0p2To0p3_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso_DR0p2To0p3_Vs_NVtx_EffArea = Ele_GammaIso_DR0p2To0p3_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso_DR0p2To0p3_Vs_NVtx_EffAreaErr = Ele_GammaIso_DR0p2To0p3_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso_DR0p2To0p3_Vs_NVtx_SlopeErr/Ele_GammaIso_DR0p2To0p3_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso_DR0p3To0p4_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_GammaIso_DR0p3To0p4_Vs_NVtx_Graph->Fit(("Ele_GammaIso_DR0p3To0p4_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso_DR0p3To0p4_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso_DR0p3To0p4_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso_DR0p3To0p4_Vs_NVtx_EffArea = Ele_GammaIso_DR0p3To0p4_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso_DR0p3To0p4_Vs_NVtx_EffAreaErr = Ele_GammaIso_DR0p3To0p4_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso_DR0p3To0p4_Vs_NVtx_SlopeErr/Ele_GammaIso_DR0p3To0p4_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso_DR0p4To0p5_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_GammaIso_DR0p4To0p5_Vs_NVtx_Graph->Fit(("Ele_GammaIso_DR0p4To0p5_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso_DR0p4To0p5_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso_DR0p4To0p5_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso_DR0p4To0p5_Vs_NVtx_EffArea = Ele_GammaIso_DR0p4To0p5_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso_DR0p4To0p5_Vs_NVtx_EffAreaErr = Ele_GammaIso_DR0p4To0p5_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso_DR0p4To0p5_Vs_NVtx_SlopeErr/Ele_GammaIso_DR0p4To0p5_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso_DR0p5To0p7_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_GammaIso_DR0p5To0p7_Vs_NVtx_Graph->Fit(("Ele_GammaIso_DR0p5To0p7_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso_DR0p5To0p7_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso_DR0p5To0p7_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso_DR0p5To0p7_Vs_NVtx_EffArea = Ele_GammaIso_DR0p5To0p7_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso_DR0p5To0p7_Vs_NVtx_EffAreaErr = Ele_GammaIso_DR0p5To0p7_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso_DR0p5To0p7_Vs_NVtx_SlopeErr/Ele_GammaIso_DR0p5To0p7_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso_DR0p7To1p0_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_GammaIso_DR0p7To1p0_Vs_NVtx_Graph->Fit(("Ele_GammaIso_DR0p7To1p0_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso_DR0p7To1p0_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso_DR0p7To1p0_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso_DR0p7To1p0_Vs_NVtx_EffArea = Ele_GammaIso_DR0p7To1p0_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso_DR0p7To1p0_Vs_NVtx_EffAreaErr = Ele_GammaIso_DR0p7To1p0_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso_DR0p7To1p0_Vs_NVtx_SlopeErr/Ele_GammaIso_DR0p7To1p0_Vs_NVtx_Slope,2));





  f1 = new TF1(("Ele_NeutralHadronIso_DR0p0To0p1_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso_DR0p0To0p1_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_EffArea = Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso_DR0p1To0p2_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso_DR0p1To0p2_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_EffArea = Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso_DR0p2To0p3_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso_DR0p2To0p3_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_EffArea = Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso_DR0p3To0p4_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso_DR0p3To0p4_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_EffArea = Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso_DR0p4To0p5_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso_DR0p4To0p5_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_EffArea = Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso_DR0p5To0p7_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso_DR0p5To0p7_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_EffArea = Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso_DR0p7To1p0_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso_DR0p7To1p0_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_EffArea = Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_Slope,2));


  f1 = new TF1(("Ele_ChargedIso_DR0p0To0p3_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_Graph->Fit(("Ele_ChargedIso_DR0p0To0p3_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_EffArea = Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_EffAreaErr = Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_SlopeErr/Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_ChargedIso_DR0p0To0p4_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_Graph->Fit(("Ele_ChargedIso_DR0p0To0p4_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_EffArea = Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_EffAreaErr = Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_SlopeErr/Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso_DR0p0To0p3_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_GammaIso_DR0p0To0p3_Vs_NVtx_Graph->Fit(("Ele_GammaIso_DR0p0To0p3_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso_DR0p0To0p3_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso_DR0p0To0p3_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso_DR0p0To0p3_Vs_NVtx_EffArea = Ele_GammaIso_DR0p0To0p3_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso_DR0p0To0p3_Vs_NVtx_EffAreaErr = Ele_GammaIso_DR0p0To0p3_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso_DR0p0To0p3_Vs_NVtx_SlopeErr/Ele_GammaIso_DR0p0To0p3_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso_DR0p0To0p4_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_GammaIso_DR0p0To0p4_Vs_NVtx_Graph->Fit(("Ele_GammaIso_DR0p0To0p4_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso_DR0p0To0p4_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso_DR0p0To0p4_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso_DR0p0To0p4_Vs_NVtx_EffArea = Ele_GammaIso_DR0p0To0p4_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso_DR0p0To0p4_Vs_NVtx_EffAreaErr = Ele_GammaIso_DR0p0To0p4_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso_DR0p0To0p4_Vs_NVtx_SlopeErr/Ele_GammaIso_DR0p0To0p4_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso_DR0p0To0p3_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso_DR0p0To0p3_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_EffArea = Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso_DR0p0To0p4_VS_NVtx"+label+"_Fit").c_str(), "pol1", 2.5, 15.5);
  f1->SetLineColor(kRed);
  Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso_DR0p0To0p4_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_EffArea = Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_Slope,2));




  cout << endl << endl 
       << "***********************************************************************************"
       << endl << endl ;

  char buffer[200];
  char buffer2[200];

  if (isMC) {
    sprintf(buffer,"%.3f +/- %.3f",Ele_ChargedIso03_Vs_NPU_EffArea,Ele_ChargedIso03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso03_Vs_NVtx_EffArea,Ele_ChargedIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_NeutralHadronIso03_Vs_NPU_EffArea,Ele_NeutralHadronIso03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso03_Vs_NVtx_EffArea,Ele_NeutralHadronIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_GammaIso03_Vs_NPU_EffArea,Ele_GammaIso03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso03_Vs_NVtx_EffArea,Ele_GammaIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffArea,Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffArea,Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_ChargedIso04_Vs_NPU_EffArea,Ele_ChargedIso04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso04_Vs_NVtx_EffArea,Ele_ChargedIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_NeutralHadronIso04_Vs_NPU_EffArea,Ele_NeutralHadronIso04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso04_Vs_NVtx_EffArea,Ele_NeutralHadronIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_GammaIso04_Vs_NPU_EffArea,Ele_GammaIso04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso04_Vs_NVtx_EffArea,Ele_GammaIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffArea,Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffArea,Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_NeutralHadronIso007_Vs_NPU_EffArea,Ele_NeutralHadronIso007_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso007_Vs_NVtx_EffArea,Ele_NeutralHadronIso007_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso007_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.5f +/- %.5f",Ele_HoverE_Vs_NPU_EffArea,Ele_HoverE_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.5f +/- %.5f",Ele_HoverE_Vs_NVtx_EffArea,Ele_HoverE_Vs_NVtx_EffAreaErr);
    cout << "Ele_HoverE_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.5f +/- %.5f",Ele_HcalDepth1OverEcal_Vs_NPU_EffArea,Ele_HcalDepth1OverEcal_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.5f +/- %.5f",Ele_HcalDepth1OverEcal_Vs_NVtx_EffArea,Ele_HcalDepth1OverEcal_Vs_NVtx_EffAreaErr);
    cout << "Ele_HcalDepth1OverEcal_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.5f +/- %.5f",Ele_HcalDepth2OverEcal_Vs_NPU_EffArea,Ele_HcalDepth2OverEcal_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.5f +/- %.5f",Ele_HcalDepth2OverEcal_Vs_NVtx_EffArea,Ele_HcalDepth2OverEcal_Vs_NVtx_EffAreaErr);
    cout << "Ele_HcalDepth2OverEcal_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_TrkIso03_Vs_NPU_EffArea,Ele_TrkIso03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_TrkIso03_Vs_NVtx_EffArea,Ele_TrkIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_TrkIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_EMIso03_Vs_NPU_EffArea,Ele_EMIso03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_EMIso03_Vs_NVtx_EffArea,Ele_EMIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_EMIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_HadIso03_Vs_NPU_EffArea,Ele_HadIso03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_HadIso03_Vs_NVtx_EffArea,Ele_HadIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_HadIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_TrkIso04_Vs_NPU_EffArea,Ele_TrkIso04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_TrkIso04_Vs_NVtx_EffArea,Ele_TrkIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_TrkIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_EMIso04_Vs_NPU_EffArea,Ele_EMIso04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_EMIso04_Vs_NVtx_EffArea,Ele_EMIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_EMIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_HadIso04_Vs_NPU_EffArea,Ele_HadIso04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_HadIso04_Vs_NVtx_EffArea,Ele_HadIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_HadIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

  } else {

    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso03_Vs_NVtx_EffArea,Ele_ChargedIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso03_Vs_NVtx_EffArea,Ele_NeutralHadronIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso03_Vs_NVtx_EffArea,Ele_GammaIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffArea,Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso04_Vs_NVtx_EffArea,Ele_ChargedIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso04_Vs_NVtx_EffArea,Ele_NeutralHadronIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso04_Vs_NVtx_EffArea,Ele_GammaIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffArea,Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso007_Vs_NVtx_EffArea,Ele_NeutralHadronIso007_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso007_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.5f +/- %.5f",Ele_HoverE_Vs_NVtx_EffArea,Ele_HoverE_Vs_NVtx_EffAreaErr);
    cout << "Ele_HoverE_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.5f +/- %.5f",Ele_HcalDepth1OverEcal_Vs_NVtx_EffArea,Ele_HcalDepth1OverEcal_Vs_NVtx_EffAreaErr);
    cout << "Ele_HcalDepth1OverEcal_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.5f +/- %.5f",Ele_HcalDepth2OverEcal_Vs_NVtx_EffArea,Ele_HcalDepth2OverEcal_Vs_NVtx_EffAreaErr);
    cout << "Ele_HcalDepth2OverEcal_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_TrkIso03_Vs_NVtx_EffArea,Ele_TrkIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_TrkIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_EMIso03_Vs_NVtx_EffArea,Ele_EMIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_EMIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_HadIso03_Vs_NVtx_EffArea,Ele_HadIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_HadIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_TrkIso04_Vs_NVtx_EffArea,Ele_TrkIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_TrkIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_EMIso04_Vs_NVtx_EffArea,Ele_EMIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_EMIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_HadIso04_Vs_NVtx_EffArea,Ele_HadIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_HadIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
  
    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_EffArea,Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso_DR0p0To0p1_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_EffArea,Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso_DR0p1To0p2_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_EffArea,Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso_DR0p2To0p3_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_EffArea,Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso_DR0p3To0p4_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_EffArea,Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso_DR0p4To0p5_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_EffArea,Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso_DR0p5To0p7_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_EffArea,Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso_DR0p7To1p0_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso_DR0p0To0p1_Vs_NVtx_EffArea,Ele_GammaIso_DR0p0To0p1_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso_DR0p0To0p1_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso_DR0p1To0p2_Vs_NVtx_EffArea,Ele_GammaIso_DR0p1To0p2_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso_DR0p1To0p2_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso_DR0p2To0p3_Vs_NVtx_EffArea,Ele_GammaIso_DR0p2To0p3_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso_DR0p2To0p3_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso_DR0p3To0p4_Vs_NVtx_EffArea,Ele_GammaIso_DR0p3To0p4_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso_DR0p3To0p4_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso_DR0p4To0p5_Vs_NVtx_EffArea,Ele_GammaIso_DR0p4To0p5_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso_DR0p4To0p5_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso_DR0p5To0p7_Vs_NVtx_EffArea,Ele_GammaIso_DR0p5To0p7_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso_DR0p5To0p7_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso_DR0p7To1p0_Vs_NVtx_EffArea,Ele_GammaIso_DR0p7To1p0_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso_DR0p7To1p0_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_EffArea,Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso_DR0p0To0p1_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_EffArea,Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso_DR0p1To0p2_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_EffArea,Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso_DR0p2To0p3_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_EffArea,Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso_DR0p3To0p4_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_EffArea,Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso_DR0p4To0p5_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_EffArea,Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso_DR0p5To0p7_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_EffArea,Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso_DR0p7To1p0_Vs_NVtx_EffArea : " << buffer2 << endl;

    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_EffArea,Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso_DR0p0To0p3_Vs_NVtx_EffArea : " << buffer2 << endl;
    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_EffArea,Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso_DR0p0To0p4_Vs_NVtx_EffArea : " << buffer2 << endl;
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso_DR0p0To0p3_Vs_NVtx_EffArea,Ele_GammaIso_DR0p0To0p3_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso_DR0p0To0p3_Vs_NVtx_EffArea : " << buffer2 << endl;
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso_DR0p0To0p4_Vs_NVtx_EffArea,Ele_GammaIso_DR0p0To0p4_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso_DR0p0To0p4_Vs_NVtx_EffArea : " << buffer2 << endl;
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_EffArea,Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso_DR0p0To0p3_Vs_NVtx_EffArea : " << buffer2 << endl;
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_EffArea,Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso_DR0p0To0p4_Vs_NVtx_EffArea : " << buffer2 << endl;



  }
   
}

void ComputeElectronIsolationEffectiveArea() {
  
//   FillIsolationHistograms("/data/blue/sixie/ntuples/HWW/mc/HwwAnalysis_s11-h115ww2l-gf-v11-pu_noskim_normalized.root", "HWW115");


  //Fall11 Z->mumu
//   FillIsolationHistograms("Fall11ZeeMC", "Fall11ZeeMC", -1);
//    ComputeEffectiveArea("Fall11ZeeMC",kTRUE, -1);

//    FillIsolationHistograms("Fall11ZeeMC", "Fall11ZeeMC", 0);
//   FillIsolationHistograms("Fall11ZeeMC", "Fall11ZeeMC", 1);
//   FillIsolationHistograms("Fall11ZeeMC", "Fall11ZeeMC", 2);
//   FillIsolationHistograms("Fall11ZeeMC", "Fall11ZeeMC", 3);
//   FillIsolationHistograms("Fall11ZeeMC", "Fall11ZeeMC", 4);
  
//    ComputeEffectiveArea("Fall11ZeeMC",kTRUE, 0);
//   ComputeEffectiveArea("Fall11ZeeMC",kTRUE, 1);
//   ComputeEffectiveArea("Fall11ZeeMC",kTRUE, 2);
//   ComputeEffectiveArea("Fall11ZeeMC",kTRUE, 3);
//   ComputeEffectiveArea("Fall11ZeeMC",kTRUE, 4);

//   Summer11MC
//   FillIsolationHistograms("Summer11MC", "Summer11MC", 0);
//   FillIsolationHistograms("Summer11MC", "Summer11MC", 1);
//   FillIsolationHistograms("Summer11MC", "Summer11MC", 2);
//   FillIsolationHistograms("Summer11MC", "Summer11MC", 3);
//   FillIsolationHistograms("Summer11MC", "Summer11MC", 4);
//   FillIsolationHistograms("Summer11MC", "Summer11MC", 5);
//   FillIsolationHistograms("Summer11MC", "Summer11MC", 6);
  
//   ComputeEffectiveArea("Summer11MC",kFALSE, 0);
//   ComputeEffectiveArea("Summer11MC",kFALSE, 1);
//   ComputeEffectiveArea("Summer11MC",kFALSE, 2);
//   ComputeEffectiveArea("Summer11MC",kFALSE, 3);
//   ComputeEffectiveArea("Summer11MC",kFALSE, 4);
//   ComputeEffectiveArea("Summer11MC",kFALSE, 5);
//   ComputeEffectiveArea("Summer11MC",kFALSE, 6);

//   Fall11MC
//   FillIsolationHistograms("Fall11MC", "Fall11MC", 0);
//   FillIsolationHistograms("Fall11MC", "Fall11MC", 1);
//   FillIsolationHistograms("Fall11MC", "Fall11MC", 2);
//   FillIsolationHistograms("Fall11MC", "Fall11MC", 3);
//   FillIsolationHistograms("Fall11MC", "Fall11MC", 4);
//   FillIsolationHistograms("Fall11MC", "Fall11MC", 5);
//   FillIsolationHistograms("Fall11MC", "Fall11MC", 6);
  
//   ComputeEffectiveArea("Fall11MC",kFALSE, 0);
//   ComputeEffectiveArea("Fall11MC",kFALSE, 1);
//   ComputeEffectiveArea("Fall11MC",kFALSE, 2);
//   ComputeEffectiveArea("Fall11MC",kFALSE, 3);
//   ComputeEffectiveArea("Fall11MC",kFALSE, 4);
//   ComputeEffectiveArea("Fall11MC",kFALSE, 5);
//   ComputeEffectiveArea("Fall11MC",kFALSE, 6);



  //Data
//   FillIsolationHistograms("Data2011", "Data2011", -1);
//   ComputeEffectiveArea("Data2011",kFALSE, -1);

//   FillIsolationHistograms("Data2011", "Data2011", 0);
//   FillIsolationHistograms("Data2011", "Data2011", 1);
//   FillIsolationHistograms("Data2011", "Data2011", 2);
//   FillIsolationHistograms("Data2011", "Data2011", 3);
//   FillIsolationHistograms("Data2011", "Data2011", 4);
//   FillIsolationHistograms("Data2011", "Data2011", 5);
//   FillIsolationHistograms("Data2011", "Data2011", 6);
  
   ComputeEffectiveArea("Data2011",kFALSE, 0);
  ComputeEffectiveArea("Data2011",kFALSE, 1);
  ComputeEffectiveArea("Data2011",kFALSE, 2);
  ComputeEffectiveArea("Data2011",kFALSE, 3);
  ComputeEffectiveArea("Data2011",kFALSE, 4);
  ComputeEffectiveArea("Data2011",kFALSE, 5);
  ComputeEffectiveArea("Data2011",kFALSE, 6);

}
