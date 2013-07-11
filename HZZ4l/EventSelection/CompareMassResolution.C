//root -l EWKAna/HZZ4l/Selection/HZZ4lSelection.C+\(\"\"\)

//================================================================================================
//
// HZZ4l selection macro
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
#include "HiggsAna/Ntupler/interface/HiggsAnaDefs.hh"
#include "HiggsAna/Ntupler/interface/TEventInfo.hh"
#include "HiggsAna/Ntupler/interface/TElectron.hh"
#include "HiggsAna/Ntupler/interface/TMuon.hh"
#include "HiggsAna/Ntupler/interface/TJet.hh"
#include "HiggsAna/Ntupler/interface/TGenParticle.hh"
#include "HiggsAna/Ntupler/interface/TPFCandidate.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"
#include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"

// output data structs
#include "HiggsAna/HZZ4l/interface/HZZKinematics.hh"
#include "HiggsAna/HZZ4l/interface/HZZGenInfo.hh"
#include "HiggsAna/HZZ4l/interface/HZZ4lDefs.hh"

#include "TMVAGui.C"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TLegend.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================

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

void CompareMassResolution(const string Label = "") 
{  
  gBenchmark->Start("HZZTemplate");
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================



  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *Mass4l_4e_IDOnly = new TH1F("Mass_4e_IDOnly", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_4m_IDOnly = new TH1F("Mass_4m_IDOnly", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_2e2m_IDOnly = new TH1F("Mass_2e2m_IDOnly", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_4l_IDOnly = new TH1F("Mass_4l_IDOnly", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);

  TH1F *Mass4l_4e_PairwiseIso = new TH1F("Mass_4e_PairwiseIso", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_4m_PairwiseIso = new TH1F("Mass_4m_PairwiseIso", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_2e2m_PairwiseIso = new TH1F("Mass_2e2m_PairwiseIso", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_4l_PairwiseIso = new TH1F("Mass_4l_PairwiseIso", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);

  TH1F *Mass4l_4e_MVAIsoLoose = new TH1F("Mass_4e_MVAIsoLoose", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_4m_MVAIsoLoose = new TH1F("Mass_4m_MVAIsoLoose", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_2e2m_MVAIsoLoose = new TH1F("Mass_2e2m_MVAIsoLoose", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_4l_MVAIsoLoose = new TH1F("Mass_4l_MVAIsoLoose", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);

  TH1F *Mass4l_4e_MVAIsoTight = new TH1F("Mass_4e_MVAIsoTight", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_4m_MVAIsoTight = new TH1F("Mass_4m_MVAIsoTight", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_2e2m_MVAIsoTight = new TH1F("Mass_2e2m_MVAIsoTight", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_4l_MVAIsoTight = new TH1F("Mass_4l_MVAIsoTight", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);

  TH1F *Mass4l_4e_MVAIsoLooseFailPairwiseIso = new TH1F("Mass_4e_MVAIsoLooseFailPairwiseIso", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_4m_MVAIsoLooseFailPairwiseIso = new TH1F("Mass_4m_MVAIsoLooseFailPairwiseIso", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_2e2m_MVAIsoLooseFailPairwiseIso = new TH1F("Mass_2e2m_MVAIsoLooseFailPairwiseIso", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);
  TH1F *Mass4l_4l_MVAIsoLooseFailPairwiseIso = new TH1F("Mass_4l_MVAIsoLooseFailPairwiseIso", ";m_{4l} [GeV/c^{2}] ; Number of Events", 80, 100, 140);


  
  //--------------------------------------------------------------------------------------------------------------
  // input ntuples
  //==============================================================================================================  
  HZZKinematics *kinematics = new HZZKinematics();
  HZZGenInfo    *geninfo = new HZZGenInfo();
  TBranch *kinematicsBr;
  TBranch *geninfoBr;

  //--------------------------------------------------------------------------------------------------------------
  // ID-Only Loop
  //==============================================================================================================  

  TFile *fFile_IDOnly = new TFile("HZZOutput_IDOnly.root", "READ");
  TTree *fTree_IDOnly = (TTree*)fFile_IDOnly->Get("hzz4l");
  assert(fTree_IDOnly);
  fTree_IDOnly->SetBranchAddress("geninfo",&geninfo);          geninfoBr = fTree_IDOnly->GetBranch("geninfo");
  fTree_IDOnly->SetBranchAddress("kinematics",&kinematics);    kinematicsBr = fTree_IDOnly->GetBranch("kinematics");

  
  for(UInt_t ientry=0; ientry<fTree_IDOnly->GetEntries(); ientry++) {       	
    kinematicsBr->GetEntry(ientry);
    if (kinematics->channel == kFourEle) Mass4l_4e_IDOnly->Fill(kinematics->m4l);
    if (kinematics->channel == kFourMu) Mass4l_4m_IDOnly->Fill(kinematics->m4l);
    if (kinematics->channel == kTwoEleTwoMu || kinematics->channel == kTwoMuTwoEle) Mass4l_2e2m_IDOnly->Fill(kinematics->m4l);
    Mass4l_4l_IDOnly->Fill(kinematics->m4l);
  }

  //--------------------------------------------------------------------------------------------------------------
  // Pairwise Iso Loop
  //==============================================================================================================  

  TFile *fFile_PairwiseIso = new TFile("HZZOutput_PairwiseIso.root", "READ");
  TTree *fTree_PairwiseIso = (TTree*)fFile_PairwiseIso->Get("hzz4l");
  assert(fTree_PairwiseIso);
  fTree_PairwiseIso->SetBranchAddress("geninfo",&geninfo);          geninfoBr = fTree_PairwiseIso->GetBranch("geninfo");
  fTree_PairwiseIso->SetBranchAddress("kinematics",&kinematics);    kinematicsBr = fTree_PairwiseIso->GetBranch("kinematics");

  
  for(UInt_t ientry=0; ientry<fTree_PairwiseIso->GetEntries(); ientry++) {       	
    kinematicsBr->GetEntry(ientry);
    if (kinematics->channel == kFourEle) Mass4l_4e_PairwiseIso->Fill(kinematics->m4l);
    if (kinematics->channel == kFourMu) Mass4l_4m_PairwiseIso->Fill(kinematics->m4l);
    if (kinematics->channel == kTwoEleTwoMu || kinematics->channel == kTwoMuTwoEle) Mass4l_2e2m_PairwiseIso->Fill(kinematics->m4l);
    Mass4l_4l_PairwiseIso->Fill(kinematics->m4l);
  }

  //--------------------------------------------------------------------------------------------------------------
  // MVALooseIso Loop
  //==============================================================================================================  

  TFile *fFile_MVAIsoLoose = new TFile("HZZOutput_MVAIsoLoose.root", "READ");
  TTree *fTree_MVAIsoLoose = (TTree*)fFile_MVAIsoLoose->Get("hzz4l");
  assert(fTree_MVAIsoLoose);
  fTree_MVAIsoLoose->SetBranchAddress("geninfo",&geninfo);          geninfoBr = fTree_MVAIsoLoose->GetBranch("geninfo");
  fTree_MVAIsoLoose->SetBranchAddress("kinematics",&kinematics);    kinematicsBr = fTree_MVAIsoLoose->GetBranch("kinematics");

  
  for(UInt_t ientry=0; ientry<fTree_MVAIsoLoose->GetEntries(); ientry++) {       	
    kinematicsBr->GetEntry(ientry);
    if (kinematics->channel == kFourEle) Mass4l_4e_MVAIsoLoose->Fill(kinematics->m4l);
    if (kinematics->channel == kFourMu) Mass4l_4m_MVAIsoLoose->Fill(kinematics->m4l);
    if (kinematics->channel == kTwoEleTwoMu || kinematics->channel == kTwoMuTwoEle) Mass4l_2e2m_MVAIsoLoose->Fill(kinematics->m4l);
    Mass4l_4l_MVAIsoLoose->Fill(kinematics->m4l);
  }

  //--------------------------------------------------------------------------------------------------------------
  // MVALooseIso Loop
  //==============================================================================================================  

  TFile *fFile_MVAIsoTight = new TFile("HZZOutput_MVAIsoTight.root", "READ");
  TTree *fTree_MVAIsoTight = (TTree*)fFile_MVAIsoTight->Get("hzz4l");
  assert(fTree_MVAIsoTight);
  fTree_MVAIsoTight->SetBranchAddress("geninfo",&geninfo);          geninfoBr = fTree_MVAIsoTight->GetBranch("geninfo");
  fTree_MVAIsoTight->SetBranchAddress("kinematics",&kinematics);    kinematicsBr = fTree_MVAIsoTight->GetBranch("kinematics");

  
  for(UInt_t ientry=0; ientry<fTree_MVAIsoTight->GetEntries(); ientry++) {       	
    kinematicsBr->GetEntry(ientry);
    if (kinematics->channel == kFourEle) Mass4l_4e_MVAIsoTight->Fill(kinematics->m4l);
    if (kinematics->channel == kFourMu) Mass4l_4m_MVAIsoTight->Fill(kinematics->m4l);
    if (kinematics->channel == kTwoEleTwoMu || kinematics->channel == kTwoMuTwoEle) Mass4l_2e2m_MVAIsoTight->Fill(kinematics->m4l);
    Mass4l_4l_MVAIsoTight->Fill(kinematics->m4l);
  }

  //--------------------------------------------------------------------------------------------------------------
  // ID-Only Loop
  //==============================================================================================================  

  TFile *fFile_MVAIsoLooseFailPairwiseIso = new TFile("HZZOutput_MVAIsoLooseFailPairwiseIso.root", "READ");
  TTree *fTree_MVAIsoLooseFailPairwiseIso = (TTree*)fFile_MVAIsoLooseFailPairwiseIso->Get("hzz4l");
  assert(fTree_MVAIsoLooseFailPairwiseIso);
  fTree_MVAIsoLooseFailPairwiseIso->SetBranchAddress("geninfo",&geninfo);          geninfoBr = fTree_MVAIsoLooseFailPairwiseIso->GetBranch("geninfo");
  fTree_MVAIsoLooseFailPairwiseIso->SetBranchAddress("kinematics",&kinematics);    kinematicsBr = fTree_MVAIsoLooseFailPairwiseIso->GetBranch("kinematics");

  
  for(UInt_t ientry=0; ientry<fTree_MVAIsoLooseFailPairwiseIso->GetEntries(); ientry++) {       	
    kinematicsBr->GetEntry(ientry);
    if (kinematics->channel == kFourEle) Mass4l_4e_MVAIsoLooseFailPairwiseIso->Fill(kinematics->m4l);
    if (kinematics->channel == kFourMu) Mass4l_4m_MVAIsoLooseFailPairwiseIso->Fill(kinematics->m4l);
    if (kinematics->channel == kTwoEleTwoMu || kinematics->channel == kTwoMuTwoEle) Mass4l_2e2m_MVAIsoLooseFailPairwiseIso->Fill(kinematics->m4l);
    Mass4l_4l_MVAIsoLooseFailPairwiseIso->Fill(kinematics->m4l);
  }




  delete kinematics;
  delete geninfo;



  //*************************************************************************
  //Set Colors
  //*************************************************************************
  Mass4l_4l_IDOnly->SetLineColor(kBlue);
  Mass4l_4e_IDOnly->SetLineColor(kBlue);
  Mass4l_4m_IDOnly->SetLineColor(kBlue);
  Mass4l_2e2m_IDOnly->SetLineColor(kBlue);
  Mass4l_4l_PairwiseIso->SetLineColor(kRed);
  Mass4l_4e_PairwiseIso->SetLineColor(kRed);
  Mass4l_4m_PairwiseIso->SetLineColor(kRed);
  Mass4l_2e2m_PairwiseIso->SetLineColor(kRed);
  Mass4l_4l_MVAIsoLoose->SetLineColor(kGreen);
  Mass4l_4e_MVAIsoLoose->SetLineColor(kGreen);
  Mass4l_4m_MVAIsoLoose->SetLineColor(kGreen);
  Mass4l_2e2m_MVAIsoLoose->SetLineColor(kGreen);
  Mass4l_4l_MVAIsoTight->SetLineColor(kMagenta);
  Mass4l_4e_MVAIsoTight->SetLineColor(kMagenta);
  Mass4l_4m_MVAIsoTight->SetLineColor(kMagenta);
  Mass4l_2e2m_MVAIsoTight->SetLineColor(kMagenta);
  Mass4l_4l_MVAIsoLooseFailPairwiseIso->SetLineColor(kCyan);
  Mass4l_4e_MVAIsoLooseFailPairwiseIso->SetLineColor(kCyan);
  Mass4l_4m_MVAIsoLooseFailPairwiseIso->SetLineColor(kCyan);
  Mass4l_2e2m_MVAIsoLooseFailPairwiseIso->SetLineColor(kCyan);


  //*************************************************************************
  //Make Plot
  //*************************************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;
  
  //*************************************************************************
  //Compare Failing and passing
  //*************************************************************************
  cv = new TCanvas("cv","cv",800,600);
  legend = new TLegend(0.15,0.75,0.45,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  legend->AddEntry(Mass4l_4l_MVAIsoLoose,"MVAIsoLoose", "L");
  legend->AddEntry(Mass4l_4l_PairwiseIso,"PairwiseIso", "L");
  legend->AddEntry(Mass4l_4l_MVAIsoLooseFailPairwiseIso,"MVAIsoLoose && !PairwiseIso", "L");
  
  Mass4l_4l_PairwiseIso->GetXaxis()->SetRangeUser(110,130);
  Mass4l_4l_PairwiseIso->SetMaximum(1900);

  Mass4l_4l_PairwiseIso->Draw("hist");
  Mass4l_4l_MVAIsoLoose->Draw("hist,same");
  Mass4l_4l_MVAIsoLooseFailPairwiseIso->Draw("hist,same");
  legend->Draw();  
  cv->SaveAs("Mass4l_4l.AdditionalSignal.gif");

  return;

  //*************************************************************************
  //Normalize Mass spectra
  //*************************************************************************
  NormalizeHist(Mass4l_4e_IDOnly);
  NormalizeHist(Mass4l_4m_IDOnly);
  NormalizeHist(Mass4l_2e2m_IDOnly);
  NormalizeHist(Mass4l_4l_IDOnly);

  Mass4l_4l_IDOnly->SetLineColor(kBlue);
  Mass4l_4e_IDOnly->SetLineColor(kBlue);
  Mass4l_4m_IDOnly->SetLineColor(kBlue);
  Mass4l_2e2m_IDOnly->SetLineColor(kBlue);

  NormalizeHist(Mass4l_4e_PairwiseIso);
  NormalizeHist(Mass4l_4m_PairwiseIso);
  NormalizeHist(Mass4l_2e2m_PairwiseIso);
  NormalizeHist(Mass4l_4l_PairwiseIso);

  Mass4l_4l_PairwiseIso->SetLineColor(kRed);
  Mass4l_4e_PairwiseIso->SetLineColor(kRed);
  Mass4l_4m_PairwiseIso->SetLineColor(kRed);
  Mass4l_2e2m_PairwiseIso->SetLineColor(kRed);

  NormalizeHist(Mass4l_4e_MVAIsoLoose);
  NormalizeHist(Mass4l_4m_MVAIsoLoose);
  NormalizeHist(Mass4l_2e2m_MVAIsoLoose);
  NormalizeHist(Mass4l_4l_MVAIsoLoose);

  Mass4l_4l_MVAIsoLoose->SetLineColor(kGreen);
  Mass4l_4e_MVAIsoLoose->SetLineColor(kGreen);
  Mass4l_4m_MVAIsoLoose->SetLineColor(kGreen);
  Mass4l_2e2m_MVAIsoLoose->SetLineColor(kGreen);


  NormalizeHist(Mass4l_4e_MVAIsoTight);
  NormalizeHist(Mass4l_4m_MVAIsoTight);
  NormalizeHist(Mass4l_2e2m_MVAIsoTight);
  NormalizeHist(Mass4l_4l_MVAIsoTight);

  Mass4l_4l_MVAIsoTight->SetLineColor(kMagenta);
  Mass4l_4e_MVAIsoTight->SetLineColor(kMagenta);
  Mass4l_4m_MVAIsoTight->SetLineColor(kMagenta);
  Mass4l_2e2m_MVAIsoTight->SetLineColor(kMagenta);


  NormalizeHist(Mass4l_4e_MVAIsoLooseFailPairwiseIso);
  NormalizeHist(Mass4l_4m_MVAIsoLooseFailPairwiseIso);
  NormalizeHist(Mass4l_2e2m_MVAIsoLooseFailPairwiseIso);
  NormalizeHist(Mass4l_4l_MVAIsoLooseFailPairwiseIso);

  Mass4l_4l_MVAIsoLooseFailPairwiseIso->SetLineColor(kCyan);
  Mass4l_4e_MVAIsoLooseFailPairwiseIso->SetLineColor(kCyan);
  Mass4l_4m_MVAIsoLooseFailPairwiseIso->SetLineColor(kCyan);
  Mass4l_2e2m_MVAIsoLooseFailPairwiseIso->SetLineColor(kCyan);



  //*************************************************************************
  //Make normalized Plot
  //*************************************************************************
  cv = new TCanvas("cv","cv",800,600);
  legend = new TLegend(0.20,0.75,0.45,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Mass4l_4l_IDOnly,"ID-Only", "LP");
  legend->AddEntry(Mass4l_4l_PairwiseIso,"PairwiseIso", "L");
  legend->AddEntry(Mass4l_4l_MVAIsoLoose,"MVAIsoLoose", "L");
  legend->AddEntry(Mass4l_4l_MVAIsoTight,"MVAIsoTight", "L");
  
  Mass4l_4l_IDOnly->GetXaxis()->SetRangeUser(110,130);
  Mass4l_4l_IDOnly->SetMaximum(0.12);
  Mass4l_4l_IDOnly->Draw("hist");
  Mass4l_4l_PairwiseIso->Draw("hist,same");
  Mass4l_4l_MVAIsoLoose->Draw("hist,same");
  Mass4l_4l_MVAIsoTight->Draw("hist,same");
//   Mass4l_4l_MVAIsoLooseFailPairwiseIso->Draw("hist,same");
  legend->Draw();  
  cv->SaveAs("Mass4l_4l.gif");




  gBenchmark->Show("WWTemplate");       
} 



