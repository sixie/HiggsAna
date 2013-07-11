//root -l CITHZZ/scripts/rootlogon.C HiggsAna/HZZ4l/Efficiency/CompareElectronDistributions.C+'("/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.ZeeTP.DataRun2012.PtBelow20.root","/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12DY_53X.ZeeTP.PtBelow20.root","ZeeTP_Bin0",0,0)'
//root -l CITHZZ/scripts/rootlogon.C HiggsAna/HZZ4l/Efficiency/CompareElectronDistributions.C+'("/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.ZeeGamma.DataRun2012.PtBelow20.root","/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12DY_53X.ZeeGamma.PtBelow20.root","ZeeGammaTP_Bin0",0,1)'

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
void CompareElectronDistributions( string dataFile, string mcFile,
                                   string Label, Int_t Option, Int_t SelectionType)
{  

  string label = "";
  if (Label != "") label = "_" + Label;


  //*****************************************************************************************
  //Plotting Setup
  //*****************************************************************************************


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *Ele_Pt_Data = new TH1F(("Ele_Pt_Data"+label).c_str(), "; Pt ; Number of Events ",  40, 0 , 20);
  TH1F *Ele_Eta_Data = new TH1F(("Ele_Eta_Data"+label).c_str(), "; Eta ; Number of Events ",  40, -2.5 , 2.5);
  TH1F *Ele_Rho_Data = new TH1F(("Ele_Rho_Data"+label).c_str(), "; Rho ; Number of Events ",  40, 0 , 50);
  TH1F *Ele_D0_Data = new TH1F(("Ele_D0_Data"+label).c_str(), "; D0 ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_DZ_Data = new TH1F(("Ele_DZ_Data"+label).c_str(), "; DZ ; Number of Events ",  40, 0 , 0.2);
  TH1F *Ele_IP3D_Data = new TH1F(("Ele_IP3D_Data"+label).c_str(), "; IP3D ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_IP3DSig_Data = new TH1F(("Ele_IP3DSig_Data"+label).c_str(), "; IP3DSig ; Number of Events ",  40, 0 , 5);
  TH1F *Ele_NBrem_Data = new TH1F(("Ele_NBrem_Data"+label).c_str(), "; NBrem ; Number of Events ",  40, 0 , 10);
  TH1F *Ele_FBrem_Data = new TH1F(("Ele_FBrem_Data"+label).c_str(), "; FBrem ; Number of Events ",  20, 0 , 1);
  TH1F *Ele_EoP_Data = new TH1F(("Ele_EoP_Data"+label).c_str(), "; EoP ; Number of Events ",  120, 0 , 2);
  TH1F *Ele_EEleoPout_Data = new TH1F(("Ele_EEleoPout_Data"+label).c_str(), "; EEleoPout ; Number of Events ",  120, 0 , 2);
  TH1F *Ele_IoEmIoP_Data = new TH1F(("Ele_IoEmIoP_Data"+label).c_str(), "; IoEmIoP ; Number of Events ",  120, -0.2 , 0.3);
  TH1F *Ele_DEta_Data = new TH1F(("Ele_DEta_Data"+label).c_str(), "; DEta ; Number of Events ",  40, 0 , 0.02);
  TH1F *Ele_DPhi_Data = new TH1F(("Ele_DPhi_Data"+label).c_str(), "; DPhi ; Number of Events ",  40, -0.3 , 0.3);
  TH1F *Ele_DEtaCalo_Data = new TH1F(("Ele_DEtaCalo_Data"+label).c_str(), "; DEtaCalo ; Number of Events ",  40, -0.2 , 0.2);
  TH1F *Ele_DPhiCalo_Data = new TH1F(("Ele_DPhiCalo_Data"+label).c_str(), "; DPhiCalo ; Number of Events ",  40, -1 , 1);
  TH1F *Ele_SigmaIEtaIEta_Data = new TH1F(("Ele_SigmaIEtaIEta_Data"+label).c_str(), "; SigmaIEtaIEta ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_SigmaIPhiIPhi_Data = new TH1F(("Ele_SigmaIPhiIPhi_Data"+label).c_str(), "; SigmaIPhiIPhi ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_EtaWidth_Data = new TH1F(("Ele_EtaWidth_Data"+label).c_str(), "; EtaWidth ; Number of Events ",  40, 0 , 0.5);
  TH1F *Ele_PhiWidth_Data = new TH1F(("Ele_PhiWidth_Data"+label).c_str(), "; PhiWidth ; Number of Events ",  40, 0 , 0.3);
  TH1F *Ele_R9_Data = new TH1F(("Ele_R9_Data"+label).c_str(), "; R9 ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_PreShowerOverRaw_Data = new TH1F(("Ele_PreShowerOverRaw_Data"+label).c_str(), "; PreShowerOverRaw ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_HoE_Data = new TH1F(("Ele_HoE_Data"+label).c_str(), "; HoE ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_GsfChi2_Data = new TH1F(("Ele_GsfChi2_Data"+label).c_str(), "; GsfChi2 ; Number of Events ",  40, 0 , 200);
  TH1F *Ele_KFChi2_Data = new TH1F(("Ele_KFChi2_Data"+label).c_str(), "; KFChi2 ; Number of Events ",  40, 0 , 10);
  TH1F *Ele_KFLayers_Data = new TH1F(("Ele_KFLayers_Data"+label).c_str(), "; KFLayers ; Number of Events ",  21, -0.5 , 20.5);
  TH1F *Ele_OneMinusSeedE1x5OverE5x5_Data = new TH1F(("Ele_OneMinusSeedE1x5OverE5x5_Data"+label).c_str(), "; OneMinusSeedE1x5OverE5x5 ; Number of Events ",  40, -1.5 , 1.5);
  TH1F *Ele_PFIso04OverPt_Data = new TH1F(("Ele_PFIso04OverPt_Data"+label).c_str(), "; PFIso04OverPt ; Number of Events ",  40, 0 , 5);
  TH1F *Ele_EGammaNonTriggeringMVA_Data = new TH1F(("Ele_EGammaNonTriggeringMVA_Data"+label).c_str(), "; EGammaNonTriggeringMVA ; Number of Events ",  40, -1 , 1);

  TH1F *Ele_Pt_MC = new TH1F(("Ele_Pt_MC"+label).c_str(), "; Pt ; Number of Events ",  40, 0 , 20);
  TH1F *Ele_Eta_MC = new TH1F(("Ele_Eta_MC"+label).c_str(), "; Eta ; Number of Events ",  40, -2.5 , 2.5);
  TH1F *Ele_Rho_MC = new TH1F(("Ele_Rho_MC"+label).c_str(), "; Rho ; Number of Events ",  40, 0 , 50);
  TH1F *Ele_D0_MC = new TH1F(("Ele_D0_MC"+label).c_str(), "; D0 ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_DZ_MC = new TH1F(("Ele_DZ_MC"+label).c_str(), "; DZ ; Number of Events ",  40, 0 , 0.2);
  TH1F *Ele_IP3D_MC = new TH1F(("Ele_IP3D_MC"+label).c_str(), "; IP3D ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_IP3DSig_MC = new TH1F(("Ele_IP3DSig_MC"+label).c_str(), "; IP3DSig ; Number of Events ",  40, 0 , 5);
  TH1F *Ele_NBrem_MC = new TH1F(("Ele_NBrem_MC"+label).c_str(), "; NBrem ; Number of Events ",  40, 0 , 10);
  TH1F *Ele_FBrem_MC = new TH1F(("Ele_FBrem_MC"+label).c_str(), "; FBrem ; Number of Events ",  20, 0 , 1);
  TH1F *Ele_EoP_MC = new TH1F(("Ele_EoP_MC"+label).c_str(), "; EoP ; Number of Events ",  120, 0 , 2);
  TH1F *Ele_EEleoPout_MC = new TH1F(("Ele_EEleoPout_MC"+label).c_str(), "; EEleoPout ; Number of Events ",  120, 0 , 2);
  TH1F *Ele_IoEmIoP_MC = new TH1F(("Ele_IoEmIoP_MC"+label).c_str(), "; IoEmIoP ; Number of Events ",  120, -0.2 , 0.3);
  TH1F *Ele_DEta_MC = new TH1F(("Ele_DEta_MC"+label).c_str(), "; DEta ; Number of Events ",  40, 0 , 0.02);
  TH1F *Ele_DPhi_MC = new TH1F(("Ele_DPhi_MC"+label).c_str(), "; DPhi ; Number of Events ",  40, -0.3 , 0.3);
  TH1F *Ele_DEtaCalo_MC = new TH1F(("Ele_DEtaCalo_MC"+label).c_str(), "; DEtaCalo ; Number of Events ",  40, -0.2 , 0.2);
  TH1F *Ele_DPhiCalo_MC = new TH1F(("Ele_DPhiCalo_MC"+label).c_str(), "; DPhiCalo ; Number of Events ",  40, -1 , 1);
  TH1F *Ele_SigmaIEtaIEta_MC = new TH1F(("Ele_SigmaIEtaIEta_MC"+label).c_str(), "; SigmaIEtaIEta ; Number of Events ",  100, 0 , 0.1);
  TH1F *Ele_SigmaIPhiIPhi_MC = new TH1F(("Ele_SigmaIPhiIPhi_MC"+label).c_str(), "; SigmaIPhiIPhi ; Number of Events ",  40, 0 , 0.1);
  TH1F *Ele_EtaWidth_MC = new TH1F(("Ele_EtaWidth_MC"+label).c_str(), "; EtaWidth ; Number of Events ",  40, 0 , 0.5);
  TH1F *Ele_PhiWidth_MC = new TH1F(("Ele_PhiWidth_MC"+label).c_str(), "; PhiWidth ; Number of Events ",  40, 0 , 0.3);
  TH1F *Ele_R9_MC = new TH1F(("Ele_R9_MC"+label).c_str(), "; R9 ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_PreShowerOverRaw_MC = new TH1F(("Ele_PreShowerOverRaw_MC"+label).c_str(), "; PreShowerOverRaw ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_HoE_MC = new TH1F(("Ele_HoE_MC"+label).c_str(), "; HoE ; Number of Events ",  40, 0 , 2);
  TH1F *Ele_GsfChi2_MC = new TH1F(("Ele_GsfChi2_MC"+label).c_str(), "; GsfChi2 ; Number of Events ",  40, 0 , 200);
  TH1F *Ele_KFChi2_MC = new TH1F(("Ele_KFChi2_MC"+label).c_str(), "; KFChi2 ; Number of Events ",  40, 0 , 10);
  TH1F *Ele_KFLayers_MC = new TH1F(("Ele_KFLayers_MC"+label).c_str(), "; KFLayers ; Number of Events ",  21, -0.5 , 20.5);
  TH1F *Ele_OneMinusSeedE1x5OverE5x5_MC = new TH1F(("Ele_OneMinusSeedE1x5OverE5x5_MC"+label).c_str(), "; OneMinusSeedE1x5OverE5x5 ; Number of Events ",  40, -1.5 , 1.5);
  TH1F *Ele_PFIso04OverPt_MC = new TH1F(("Ele_PFIso04OverPt_MC"+label).c_str(), "; PFIso04OverPt ; Number of Events ",  40, 0 , 5);
  TH1F *Ele_EGammaNonTriggeringMVA_MC = new TH1F(("Ele_EGammaNonTriggeringMVA_MC"+label).c_str(), "; EGammaNonTriggeringMVA ; Number of Events ",  40, -1 , 1);


  //*****************************************************************************************
  //DataEleTree
  //*****************************************************************************************
  citana::ElectronTree DataEleTree;
  DataEleTree.LoadTree(dataFile.c_str());
  DataEleTree.InitTree();

  for(UInt_t ientry=0; ientry < DataEleTree.tree_->GetEntries(); ientry++) {       	
    DataEleTree.tree_->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //don't evaluate performance using training events
    if (DataEleTree.fElePt < 7) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(DataEleTree.fEleEta) < 0.8) subdet = 0;
    else if (fabs(DataEleTree.fEleEta) < 1.479) subdet = 1;
    else subdet = 2;
    Int_t ptBin = 0;
    if (DataEleTree.fElePt > 10) ptBin = 1;


    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);
    if (Option == 10) passCuts = (ptBin == 0 );

    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

//     if (!(DataEleTree.fEleNMissHits <= 1)) continue;
//     if( DataEleTree.fEleIP3dSig > 4 ) continue;

    //make isolation cut

    if( DataEleTree.fElePFIso04 / DataEleTree.fElePt > 0.4) continue;
    if (SelectionType == 0) {
      if(!(DataEleTree.fWeight > 85 && DataEleTree.fWeight < 95)) continue;
    } else if (SelectionType == 1) {
      if(!(DataEleTree.fWeight > 80 && DataEleTree.fWeight < 95)) continue;
    }
    
    cout << DataEleTree.fRunNumber << " " << DataEleTree.fLumiSectionNumber  << " " << DataEleTree.fEventNumber  << endl;

    //Fill Histograms
    Ele_SigmaIEtaIEta_Data->Fill(DataEleTree.fEleSigmaIEtaIEta);
    Ele_Pt_Data->Fill(DataEleTree.fElePt);
    Ele_Eta_Data->Fill(DataEleTree.fEleEta);
    Ele_Rho_Data->Fill(DataEleTree.fRho);
    Ele_D0_Data->Fill(DataEleTree.fEleD0);
    Ele_DZ_Data->Fill(DataEleTree.fEleDZ);
    Ele_IP3D_Data->Fill(DataEleTree.fEleIP3d);
    Ele_IP3DSig_Data->Fill(DataEleTree.fEleIP3dSig);
    Ele_NBrem_Data->Fill(DataEleTree.fEleNBrem);
    Ele_FBrem_Data->Fill(DataEleTree.fEleFBrem);
//     Ele_EoP_Data->Fill(DataEleTree.fEleEOverP * 1.1); //correct by 10% in data
//     Ele_EEleoPout_Data->Fill(DataEleTree.fEleEEleClusterOverPout * 1.025);
//     Ele_IoEmIoP_Data->Fill(DataEleTree.fEleOneOverEMinusOneOverP - 0.005);
    Ele_EoP_Data->Fill(DataEleTree.fEleEOverP ); 
    Ele_EEleoPout_Data->Fill(DataEleTree.fEleEEleClusterOverPout );
    Ele_IoEmIoP_Data->Fill(DataEleTree.fEleOneOverEMinusOneOverP );
    Ele_DEta_Data->Fill(DataEleTree.fEleDEtaIn);
    Ele_DPhi_Data->Fill(DataEleTree.fEleDPhiIn);
    Ele_DEtaCalo_Data->Fill(DataEleTree.fEledEtaCalo);
    Ele_DPhiCalo_Data->Fill(DataEleTree.fEledPhiCalo);
    Ele_SigmaIEtaIEta_Data->Fill(DataEleTree.fEleSigmaIEtaIEta);
    Ele_SigmaIPhiIPhi_Data->Fill(DataEleTree.fEleSigmaIPhiIPhi);
    Ele_EtaWidth_Data->Fill(DataEleTree.fEleSCEtaWidth);
    Ele_PhiWidth_Data->Fill(DataEleTree.fEleSCPhiWidth);
    Ele_R9_Data->Fill(DataEleTree.fEleR9);
    Ele_PreShowerOverRaw_Data->Fill(DataEleTree.fElePreShowerOverRaw);
    Ele_HoE_Data->Fill(DataEleTree.fEleHoverE);
    Ele_GsfChi2_Data->Fill(DataEleTree.fEleGsfTrackChi2OverNdof);
    Ele_KFChi2_Data->Fill(DataEleTree.fEleKFTrackChi2OverNDoF);
    Ele_KFLayers_Data->Fill(DataEleTree.fEleKFTrackNLayersWithMeasurement);
    Ele_OneMinusSeedE1x5OverE5x5_Data->Fill(DataEleTree.fEleOneMinusSeedE1x5OverE5x5);
    Ele_PFIso04OverPt_Data->Fill(DataEleTree.fElePFIso04 / DataEleTree.fElePt);
    Ele_EGammaNonTriggeringMVA_Data->Fill(DataEleTree.fEGammaNonTriggeringMVA);


  } 
  




  //*****************************************************************************************
  //MCEleTree
  //*****************************************************************************************
  citana::ElectronTree MCEleTree;
  MCEleTree.LoadTree(mcFile.c_str());
  MCEleTree.InitTree();

  for(UInt_t ientry=0; ientry < MCEleTree.tree_->GetEntries(); ientry++) {       	
    MCEleTree.tree_->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //don't evaluate performance using training events
    if (MCEleTree.fElePt < 7) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(MCEleTree.fEleEta) < 0.8) subdet = 0;
    else if (fabs(MCEleTree.fEleEta) < 1.479) subdet = 1;
    else subdet = 2;
    Int_t ptBin = 0;
    if (MCEleTree.fElePt > 10) ptBin = 1;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);
    if (Option == 10) passCuts = (ptBin == 0 );

    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

//     if (!(MCEleTree.fEleNMissHits <= 1)) continue;
//     if( MCEleTree.fEleIP3dSig > 4 ) continue;

    if( MCEleTree.fElePFIso04 / MCEleTree.fElePt > 0.4) continue;
    if (SelectionType == 0) {
      if(!(MCEleTree.fWeight > 85 && MCEleTree.fWeight < 95)) continue;
    } else if (SelectionType == 1) {
      if(!(MCEleTree.fWeight > 80 && MCEleTree.fWeight < 95)) continue;
    }
    
    //Fill Histograms
    Ele_SigmaIEtaIEta_MC->Fill(MCEleTree.fEleSigmaIEtaIEta);
    Ele_Pt_MC->Fill(MCEleTree.fElePt);
    Ele_Eta_MC->Fill(MCEleTree.fEleEta);
    Ele_Rho_MC->Fill(MCEleTree.fRho);
    Ele_D0_MC->Fill(MCEleTree.fEleD0);
    Ele_DZ_MC->Fill(MCEleTree.fEleDZ);
    Ele_IP3D_MC->Fill(MCEleTree.fEleIP3d);
    Ele_IP3DSig_MC->Fill(MCEleTree.fEleIP3dSig);
    Ele_NBrem_MC->Fill(MCEleTree.fEleNBrem);
    Ele_FBrem_MC->Fill(MCEleTree.fEleFBrem);
    Ele_EoP_MC->Fill(MCEleTree.fEleEOverP);
    Ele_EEleoPout_MC->Fill(MCEleTree.fEleEEleClusterOverPout);
    Ele_IoEmIoP_MC->Fill(MCEleTree.fEleOneOverEMinusOneOverP);
    Ele_DEta_MC->Fill(MCEleTree.fEleDEtaIn);
    Ele_DPhi_MC->Fill(MCEleTree.fEleDPhiIn);
    Ele_DEtaCalo_MC->Fill(MCEleTree.fEledEtaCalo);
    Ele_DPhiCalo_MC->Fill(MCEleTree.fEledPhiCalo);
    Ele_SigmaIEtaIEta_MC->Fill(MCEleTree.fEleSigmaIEtaIEta);
    Ele_SigmaIPhiIPhi_MC->Fill(MCEleTree.fEleSigmaIPhiIPhi);
    Ele_EtaWidth_MC->Fill(MCEleTree.fEleSCEtaWidth);
    Ele_PhiWidth_MC->Fill(MCEleTree.fEleSCPhiWidth);
    Ele_R9_MC->Fill(MCEleTree.fEleR9);
    Ele_PreShowerOverRaw_MC->Fill(MCEleTree.fElePreShowerOverRaw);
    Ele_HoE_MC->Fill(MCEleTree.fEleHoverE);
    Ele_GsfChi2_MC->Fill(MCEleTree.fEleGsfTrackChi2OverNdof);
    Ele_KFChi2_MC->Fill(MCEleTree.fEleKFTrackChi2OverNDoF);
    Ele_KFLayers_MC->Fill(MCEleTree.fEleKFTrackNLayersWithMeasurement);
    Ele_OneMinusSeedE1x5OverE5x5_MC->Fill(MCEleTree.fEleOneMinusSeedE1x5OverE5x5);
    Ele_PFIso04OverPt_MC->Fill(MCEleTree.fElePFIso04 / MCEleTree.fElePt);
    Ele_EGammaNonTriggeringMVA_MC->Fill(MCEleTree.fEGammaNonTriggeringMVA );


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
  legend->AddEntry(Ele_Pt_Data, "Data", "L");
  legend->AddEntry(Ele_Pt_MC, "MC", "L");
  Ele_Pt_Data->SetLineColor(kRed);
  Ele_Pt_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_Pt_Data);
  NormalizeHist(Ele_Pt_MC);
  Ele_Pt_Data->SetMaximum( 1.2*max(Ele_Pt_MC->GetMaximum(),Ele_Pt_Data->GetMaximum()) );
  Ele_Pt_Data->Draw("hist");
  Ele_Pt_MC->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_Pt" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_Eta_Data, "Data", "L");
  legend->AddEntry(Ele_Eta_MC, "MC", "L");
  Ele_Eta_Data->SetLineColor(kRed);
  Ele_Eta_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_Eta_Data);
  NormalizeHist(Ele_Eta_MC);
  Ele_Eta_Data->SetMaximum( 1.2*max(Ele_Eta_MC->GetMaximum(),Ele_Eta_Data->GetMaximum()) );
  Ele_Eta_Data->Draw("hist");
  Ele_Eta_MC->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_Eta" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.70,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_SigmaIEtaIEta_Data, "Data", "L");
  legend->AddEntry(Ele_SigmaIEtaIEta_MC, "MC", "L");
  Ele_SigmaIEtaIEta_Data->SetLineColor(kRed);
  Ele_SigmaIEtaIEta_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_SigmaIEtaIEta_Data);
  NormalizeHist(Ele_SigmaIEtaIEta_MC);
  Ele_SigmaIEtaIEta_Data->SetMaximum( 1.2*max(Ele_SigmaIEtaIEta_MC->GetMaximum(),Ele_SigmaIEtaIEta_Data->GetMaximum()) );
  Ele_SigmaIEtaIEta_Data->GetXaxis()->SetRangeUser(0,0.04);
  Ele_SigmaIEtaIEta_Data->Draw("hist");
  Ele_SigmaIEtaIEta_MC->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_SigmaIEtaIEta" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.60,0.94,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_FBrem_Data, "Data", "L");
  legend->AddEntry(Ele_FBrem_MC, "MC", "L");
  Ele_FBrem_Data->SetLineColor(kRed);
  Ele_FBrem_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_FBrem_Data);
  NormalizeHist(Ele_FBrem_MC);
  Ele_FBrem_Data->SetMaximum( 1.2*max(Ele_FBrem_MC->GetMaximum(),Ele_FBrem_Data->GetMaximum()) );
  Ele_FBrem_Data->Draw("hist");
  Ele_FBrem_MC->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_FBrem" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.50,0.50,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_EoP_Data, "Data", "L");
  legend->AddEntry(Ele_EoP_MC, "MC", "L");
  Ele_EoP_Data->SetLineColor(kRed);
  Ele_EoP_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_EoP_Data);
  NormalizeHist(Ele_EoP_MC);
  Ele_EoP_Data->SetMaximum( 1.2*max(Ele_EoP_MC->GetMaximum(),Ele_EoP_Data->GetMaximum()) );
  Ele_EoP_Data->Draw("hist");
  Ele_EoP_MC->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_EoP" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.50,0.50,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_EEleoPout_Data, "Data", "L");
  legend->AddEntry(Ele_EEleoPout_MC, "MC", "L");
  Ele_EEleoPout_Data->SetLineColor(kRed);
  Ele_EEleoPout_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_EEleoPout_Data);
  NormalizeHist(Ele_EEleoPout_MC);
  Ele_EEleoPout_Data->SetMaximum( 1.2*max(Ele_EEleoPout_MC->GetMaximum(),Ele_EEleoPout_Data->GetMaximum()) );
  Ele_EEleoPout_Data->Draw("hist");
  Ele_EEleoPout_MC->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_EEleoPout" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.50,0.50,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_IoEmIoP_Data, "Data", "L");
  legend->AddEntry(Ele_IoEmIoP_MC, "MC", "L");
  Ele_IoEmIoP_Data->SetLineColor(kRed);
  Ele_IoEmIoP_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_IoEmIoP_Data);
  NormalizeHist(Ele_IoEmIoP_MC);
  Ele_IoEmIoP_Data->SetMaximum( 1.2*max(Ele_IoEmIoP_MC->GetMaximum(),Ele_IoEmIoP_Data->GetMaximum()) );
  Ele_IoEmIoP_Data->Draw("hist");
  Ele_IoEmIoP_MC->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_IoEmIoP" + label + ".gif").c_str());




  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.70,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_R9_Data, "Data", "L");
  legend->AddEntry(Ele_R9_MC, "MC", "L");
  Ele_R9_Data->SetLineColor(kRed);
  Ele_R9_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_R9_Data);
  NormalizeHist(Ele_R9_MC);
  Ele_R9_Data->SetMaximum( 1.2*max(Ele_R9_MC->GetMaximum(),Ele_R9_Data->GetMaximum()) );
  Ele_R9_Data->Draw("hist");
  Ele_R9_MC->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_R9" + label + ".gif").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_PhiWidth_Data, "Data", "L");
  legend->AddEntry(Ele_PhiWidth_MC, "MC", "L");
  Ele_PhiWidth_Data->SetLineColor(kRed);
  Ele_PhiWidth_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_PhiWidth_Data);
  NormalizeHist(Ele_PhiWidth_MC);
  Ele_PhiWidth_Data->SetMaximum( 1.2*max(Ele_PhiWidth_MC->GetMaximum(),Ele_PhiWidth_Data->GetMaximum()) );
  Ele_PhiWidth_Data->Draw("hist");
  Ele_PhiWidth_MC->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_PhiWidth" + label + ".gif").c_str());


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_IP3DSig_Data, "Data", "L");
  legend->AddEntry(Ele_IP3DSig_MC, "MC", "L");
  Ele_IP3DSig_Data->SetLineColor(kRed);
  Ele_IP3DSig_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_IP3DSig_Data);
  NormalizeHist(Ele_IP3DSig_MC);
  Ele_IP3DSig_Data->SetMaximum( 1.2*max(Ele_IP3DSig_MC->GetMaximum(),Ele_IP3DSig_Data->GetMaximum()) );
  Ele_IP3DSig_Data->Draw("hist");
  Ele_IP3DSig_MC->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_IP3DSig" + label + ".gif").c_str());


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_DEta_Data, "Data", "L");
  legend->AddEntry(Ele_DEta_MC, "MC", "L");
  Ele_DEta_Data->SetLineColor(kRed);
  Ele_DEta_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_DEta_Data);
  NormalizeHist(Ele_DEta_MC);
  Ele_DEta_Data->SetMaximum( 1.2*max(Ele_DEta_MC->GetMaximum(),Ele_DEta_Data->GetMaximum()) );
  Ele_DEta_Data->Draw("hist");
  Ele_DEta_MC->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_DEta" + label + ".gif").c_str());


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_DPhi_Data, "Data", "L");
  legend->AddEntry(Ele_DPhi_MC, "MC", "L");
  Ele_DPhi_Data->SetLineColor(kRed);
  Ele_DPhi_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_DPhi_Data);
  NormalizeHist(Ele_DPhi_MC);
  Ele_DPhi_Data->SetMaximum( 1.2*max(Ele_DPhi_MC->GetMaximum(),Ele_DPhi_Data->GetMaximum()) );
  Ele_DPhi_Data->Draw("hist");
  Ele_DPhi_MC->Draw("hist,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_DPhi" + label + ".gif").c_str());



  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.50,0.50,0.80);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Ele_EGammaNonTriggeringMVA_Data, "Data", "L");
  legend->AddEntry(Ele_EGammaNonTriggeringMVA_MC, "MC", "L");
  Ele_EGammaNonTriggeringMVA_Data->SetLineColor(kRed);
  Ele_EGammaNonTriggeringMVA_MC->SetLineColor(kBlue);
  NormalizeHist(Ele_EGammaNonTriggeringMVA_Data);
  NormalizeHist(Ele_EGammaNonTriggeringMVA_MC);
  Ele_EGammaNonTriggeringMVA_Data->SetMaximum( 1.2*max(Ele_EGammaNonTriggeringMVA_MC->GetMaximum(),Ele_EGammaNonTriggeringMVA_Data->GetMaximum()) );
  Ele_EGammaNonTriggeringMVA_Data->Draw("hist");
  Ele_EGammaNonTriggeringMVA_MC->Draw("e1,same");
  legend->Draw();
  cv->SaveAs(("EleDistribution_DataVsMC_EGammaNonTriggeringMVA" + label + ".gif").c_str());


  gBenchmark->Show("WWTemplate");       
} 

