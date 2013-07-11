//root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/MuonMVA/output/MuonNtuple.Real.Subdet0LowPt.root","/home/sixie/CMSSW_analysis/src/MuonMVA/output/MuonNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt",0)'
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
#include "TLegend.h"

// define structures to read in ntuple
#include "HiggsAna/Ntupler/interface/HiggsAnaDefs.hh"
#include "HiggsAna/Ntupler/interface/TEventInfo.hh"
#include "HiggsAna/Ntupler/interface/TElectron.hh"
#include "HiggsAna/Ntupler/interface/TPhoton.hh"
#include "HiggsAna/Ntupler/interface/TMuon.hh"
#include "HiggsAna/Ntupler/interface/TJet.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
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
//
//*************************************************************************************************
TGraphAsymmErrors* MakeSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, string name, Bool_t GreaterThan = kTRUE ) {
  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double SigEff[nPoints];
  double BkgEff[nPoints];
  double SigEffErrLow[nPoints];
  double SigEffErrHigh[nPoints];
  double BkgEffErrLow[nPoints];
  double BkgEffErrHigh[nPoints];
  double NSigTotal = 0;
  double NBkgTotal = 0;
  
  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
    NBkgTotal += bkgHist->GetBinContent(q);
  }

  if (GreaterThan) {
    for(UInt_t b=0; b < nPoints; ++b) {
      Double_t nsig = 0;
      Double_t nbkg = 0;
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nsig += signalHist->GetBinContent(q);
        nbkg += bkgHist->GetBinContent(q);
      }

      Double_t ratio;
      Double_t errLow;
      Double_t errHigh;     
      Double_t n1 = 0;
      Double_t n2 = 0;

      n1 = TMath::Nint(nsig);
      n2 = TMath::Nint(NSigTotal);
      mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
      SigEff[b] = ratio;
      SigEffErrLow[b] = 0;
      SigEffErrHigh[b] = 0;
//     SigEffErrLow[b] = errLow;
//     SigEffErrHigh[b] = errHigh;

      n1 = TMath::Nint(nbkg);
      n2 = TMath::Nint(NBkgTotal);
      mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
      BkgEff[b] = ratio;
      BkgEffErrLow[b] = 0;
      BkgEffErrHigh[b] = 0;
//     BkgEffErrLow[b] = errLow;
//     BkgEffErrHigh[b] = errHigh;
    }  
  } else {
    for(UInt_t b=0; b < nPoints; ++b) {
      Double_t nsig = 0;
      Double_t nbkg = 0;
      for (UInt_t q=0; q <= b; ++q) {
        nsig += signalHist->GetBinContent(q);
        nbkg += bkgHist->GetBinContent(q);
      }
    
      Double_t ratio;
      Double_t errLow;
      Double_t errHigh;     
      Double_t n1 = 0;
      Double_t n2 = 0;

      n1 = TMath::Nint(nsig);
      n2 = TMath::Nint(NSigTotal);
      mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
      SigEff[b] = ratio;
      SigEffErrLow[b] = 0;
      SigEffErrHigh[b] = 0;
//     SigEffErrLow[b] = errLow;
//     SigEffErrHigh[b] = errHigh;

      n1 = TMath::Nint(nbkg);
      n2 = TMath::Nint(NBkgTotal);
      mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
      BkgEff[b] = ratio;
      BkgEffErrLow[b] = 0;
      BkgEffErrHigh[b] = 0;
//     BkgEffErrLow[b] = errLow;
//     BkgEffErrHigh[b] = errHigh;
    }
  }

  TGraphAsymmErrors *tmpSigEffVsBkgEff = new TGraphAsymmErrors (nPoints, BkgEff, SigEff, BkgEffErrLow, BkgEffErrHigh, SigEffErrLow, SigEffErrHigh );
  tmpSigEffVsBkgEff->SetName(name.c_str());
  tmpSigEffVsBkgEff->SetTitle("");
  tmpSigEffVsBkgEff->GetXaxis()->SetTitle("Bkg Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitle("Signal Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsBkgEff->GetXaxis()->SetTitleOffset(1.05);
  tmpSigEffVsBkgEff->SetMarkerSize(0.5);
  tmpSigEffVsBkgEff->SetMarkerStyle(20);

  return tmpSigEffVsBkgEff;
}

//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeCurrentWPSigEffVsBkgEffGraph(Double_t signalEff, Double_t bkgEff, string name ) {
  //Make Met Plots
  double SigEff[1];
  double BkgEff[1];
  double SigEffErrLow[1];
  double SigEffErrHigh[1];
  double BkgEffErrLow[1];
  double BkgEffErrHigh[1];
  double NSigTotal = 0;
  double NBkgTotal = 0;
  double cutValue;

  SigEff[0] = signalEff;
  SigEffErrLow[0] = 0;
  SigEffErrHigh[0] = 0;
  BkgEff[0] = bkgEff;
  BkgEffErrLow[0] = 0;
  BkgEffErrHigh[0] = 0;

  TGraphAsymmErrors *tmpSigEffVsBkgEff = new TGraphAsymmErrors (1, BkgEff, SigEff, BkgEffErrLow, BkgEffErrHigh , SigEffErrLow, SigEffErrHigh );
  tmpSigEffVsBkgEff->SetName(name.c_str());
  tmpSigEffVsBkgEff->SetTitle("");
  tmpSigEffVsBkgEff->GetXaxis()->SetTitle("Bkg Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitle("Signal Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsBkgEff->GetXaxis()->SetTitleOffset(1.05);
  tmpSigEffVsBkgEff->SetMarkerColor(kBlack);
  tmpSigEffVsBkgEff->SetLineColor(kBlack);
  tmpSigEffVsBkgEff->SetMarkerSize(1.5);

  return tmpSigEffVsBkgEff;
}


//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeSigEffVsCutValueGraph(TH1F* signalHist, string name ) {

  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double cutValue[nPoints];
  double cutValueErr[nPoints];
  double SigEff[nPoints];
  double SigEffErrLow[nPoints];
  double SigEffErrHigh[nPoints];
  double NSigTotal = 0;
  
  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }

  for(UInt_t b=0; b < nPoints; ++b) {
    cutValue[b] = signalHist->GetXaxis()->GetBinCenter(b);
    cutValueErr[b] = 0;
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    SigEff[b] = ratio;
    SigEffErrLow[b] = 0;
    SigEffErrHigh[b] = 0;
//     SigEffErrLow[b] = errLow;
//     SigEffErrHigh[b] = errHigh;

  }

  TGraphAsymmErrors *tmpSigEffVsCut = new TGraphAsymmErrors (nPoints, cutValue, SigEff, cutValueErr, cutValueErr, SigEffErrLow, SigEffErrHigh  );
  tmpSigEffVsCut->SetName(name.c_str());
  tmpSigEffVsCut->SetTitle("");
  tmpSigEffVsCut->GetXaxis()->SetTitle("Cut Value");
  tmpSigEffVsCut->GetYaxis()->SetTitle("Efficiency");
  tmpSigEffVsCut->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsCut->GetXaxis()->SetTitleOffset(1.05);

  return tmpSigEffVsCut;
}


//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeCurrentWPSigEffVsCutValueGraph(TH1F* signalHist, string name, Double_t myCutValue ) {
  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double cutValue[1] = {0};
  double cutValueErr[1] = {0};
  double SigEff[1] = {0};
  double SigEffErrLow[1] = {0};
  double SigEffErrHigh[1] = {0};
  double NSigTotal = 0;
  
  Double_t effDiff = 9999;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }

  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    
      cout << myCutValue << " : " << signalHist->GetXaxis()->GetBinCenter(b) << " , " << cutValue[0] << endl;
    if (fabs(myCutValue - signalHist->GetXaxis()->GetBinCenter(b)) < fabs(myCutValue - cutValue[0])) {
      SigEff[0] = ratio;
      SigEffErrLow[0] = 0;
      SigEffErrHigh[0] = 0;
//     SigEffErrLow[0] = errLow;
//     SigEffErrHigh[0] = errHigh;
      cutValue[0] = signalHist->GetXaxis()->GetBinCenter(b);
      cutValueErr[0] = 0;
    }
  }

//   cout << "Final: " << cutValue[0] << " , " << SigEff[0] << endl;

  TGraphAsymmErrors *tmpSigEffVsCut = new TGraphAsymmErrors (1, cutValue, SigEff, cutValueErr, cutValueErr, SigEffErrLow, SigEffErrHigh  );
  tmpSigEffVsCut->SetName(name.c_str());
  tmpSigEffVsCut->SetTitle("");
  tmpSigEffVsCut->GetXaxis()->SetTitle("Cut Value");
  tmpSigEffVsCut->GetYaxis()->SetTitle("Efficiency");
  tmpSigEffVsCut->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsCut->GetXaxis()->SetTitleOffset(1.05);
  tmpSigEffVsCut->SetMarkerColor(kBlack);
  tmpSigEffVsCut->SetLineColor(kBlack);
  tmpSigEffVsCut->SetMarkerSize(1.5);

  return tmpSigEffVsCut;
}



//*************************************************************************************************
//
//*************************************************************************************************
Double_t FindCutValueAtFixedEfficiency(TH1F* hist, Double_t targetEff ) {
  //Make Met Plots


  Double_t targetCutValue = -9999;
  Double_t bestCurrentEff = 0;
  const UInt_t nPoints = hist->GetXaxis()->GetNbins();
  double NSigTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += hist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += hist->GetBinContent(q);
    }

    Double_t ratio = nsig / NSigTotal;
//     cout << targetEff << " : " << ratio << " , " << Hist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetEff - ratio) < fabs(targetEff - bestCurrentEff)) {
      targetCutValue = hist->GetXaxis()->GetBinCenter(b);
      bestCurrentEff = ratio;
    }
  }

  return targetCutValue;
}





//*************************************************************************************************
//
//*************************************************************************************************
Double_t FindBkgEffAtFixedSignalEfficiency(TH1F* signalHist, TH1F* bkgHist, Double_t targetSignalEff, Bool_t GreaterThan = kTRUE ) {
  //Make Met Plots


  Double_t targetBkgEff = 0;
  Double_t bestCurrentSignalEff = 0;
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double NSigTotal = 0;
  double NBkgTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
    NBkgTotal += bkgHist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    Double_t nbkg = 0;

    if (GreaterThan) {
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nbkg += bkgHist->GetBinContent(q);
      }
    } else {
      for (UInt_t q=0; q <= b; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
      for (UInt_t q=0; q <= b; ++q) {
        nbkg += bkgHist->GetBinContent(q);
      }
    }

    Double_t ratio = nsig / NSigTotal;
    Double_t bkgEff = nbkg / NBkgTotal;
//     cout << targetSignalEff << " : " << ratio << " , " << bkgEff << " : " << signalHist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetSignalEff - ratio) < fabs(targetSignalEff - bestCurrentSignalEff)) {
      bestCurrentSignalEff = ratio;
      targetBkgEff = bkgEff;
    }
  }

  return targetBkgEff;
}

//*************************************************************************************************
//
//*************************************************************************************************
Double_t FindSigEffAtFixedBkgEfficiency(TH1F* signalHist, TH1F* bkgHist, Double_t targetBkgEff, Bool_t GreaterThan = kTRUE ) {
  //Make Met Plots

  Double_t targetSignalEff = 0;
  Double_t bestCurrentBkgEff = 0;
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double NSigTotal = 0;
  double NBkgTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
    NBkgTotal += bkgHist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {

    Double_t nsig = 0;
    Double_t nbkg = 0;

    if (GreaterThan) {
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nbkg += bkgHist->GetBinContent(q);
      }
    } else {
      for (UInt_t q=0; q <= b; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
      for (UInt_t q=0; q <= b; ++q) {
        nbkg += bkgHist->GetBinContent(q);
      }
    }

    Double_t sigEff = nsig / NSigTotal;
    Double_t bkgEff = nbkg / NBkgTotal;
//     cout << targetSignalEff << " : " << ratio << " , " << bkgEff << " : " << signalHist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetBkgEff - bkgEff) < fabs(targetBkgEff - bestCurrentBkgEff)) {
      bestCurrentBkgEff = bkgEff;
      targetSignalEff = sigEff;
    }
  }

  return targetSignalEff;
}



Bool_t passEleIsoMVASameAsDetIso025( Double_t fMuPt, Double_t fMuEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fMuEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fMuPt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;

  Double_t MVACut = -999;
  //Cut Values for MVA-V10 (Detector Based Iso)
  if (MVABin == 0) MVACut = -0.616;
  if (MVABin == 1) MVACut = -0.740;
  if (MVABin == 2) MVACut = -0.045;
  if (MVABin == 3) MVACut = 0.205;

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}

Bool_t passEleIsoMVASameAsDetIso015( Double_t fMuPt, Double_t fMuEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fMuEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fMuPt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;

  Double_t MVACut = -999;
  //Cut Values for MVA-V10 (Detector Based Iso)
  if (MVABin == 0) MVACut = -0.277;
  if (MVABin == 1) MVACut = -0.566;
  if (MVABin == 2) MVACut = 0.525;
  if (MVABin == 3) MVACut = 0.600;

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}






//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void MakeElectronIsoMVAPerformancePlots(string RealElectronFile, string FakeElectronFile, string Label, Int_t Option, Int_t NVtxBin = -1)
{  

  string label = "";
  if (Label != "") label = "_" + Label;


  //*****************************************************************************************
  //Plotting Setup
  //*****************************************************************************************


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  

  TH1F *EleDetIso03_Real = new TH1F(("EleDetIso03_Real"+label).c_str(), "; DetIso03 ; Number of Events ",  10000, 0 , 10);
  TH1F *ElePFIsoHWW_Real = new TH1F(("ElePFIsoHWW_Real"+label).c_str(), "; PFIsoHWW ; Number of Events ",  10000, 0 , 10);
  TH1F *ElePFIso03_Real = new TH1F(("ElePFIso03_Real"+label).c_str(), "; PFIso03 ; Number of Events ",  10000, 0 , 10);
  TH1F *ElePFIso04_Real = new TH1F(("ElePFIso04_Real"+label).c_str(), "; PFIso04 ; Number of Events ",  10000, 0 , 10);
  TH1F *EleIsoBDT_V0_Real = new TH1F(("EleIsoBDT_V0_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDT_V1_Real = new TH1F(("EleIsoBDT_V1_Real"+label).c_str(), "; BDTG V1 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDT_V2_Real = new TH1F(("EleIsoBDT_V2_Real"+label).c_str(), "; BDTG V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDT_V3_Real = new TH1F(("EleIsoBDT_V3_Real"+label).c_str(), "; BDTG V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDT_V4_Real = new TH1F(("EleIsoBDT_V4_Real"+label).c_str(), "; BDTG V4 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDT_V5_Real = new TH1F(("EleIsoBDT_V5_Real"+label).c_str(), "; BDTG V5 ; Number of Events ",  10000, -2 , 2);

  TH1F *EleDetIso03_Fake = new TH1F(("EleDetIso03_Fake"+label).c_str(), "; DetIso03 ; Number of Events ",  10000, 0 , 10);
  TH1F *ElePFIsoHWW_Fake = new TH1F(("ElePFIsoHWW_Fake"+label).c_str(), "; PFIsoHWW ; Number of Events ",  10000, 0 , 10);
  TH1F *ElePFIso03_Fake = new TH1F(("ElePFIso03_Fake"+label).c_str(), "; PFIso03 ; Number of Events ",  10000, 0 , 10);
  TH1F *ElePFIso04_Fake = new TH1F(("ElePFIso04_Fake"+label).c_str(), "; PFIso04 ; Number of Events ",  10000, 0 , 10);
  TH1F *EleIsoBDT_V0_Fake = new TH1F(("EleIsoBDT_V0_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDT_V1_Fake = new TH1F(("EleIsoBDT_V1_Fake"+label).c_str(), "; BDTG V1 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDT_V2_Fake = new TH1F(("EleIsoBDT_V2_Fake"+label).c_str(), "; BDTG V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDT_V3_Fake = new TH1F(("EleIsoBDT_V3_Fake"+label).c_str(), "; BDTG V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDT_V4_Fake = new TH1F(("EleIsoBDT_V4_Fake"+label).c_str(), "; BDTG V4 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDT_V5_Fake = new TH1F(("EleIsoBDT_V5_Fake"+label).c_str(), "; BDTG V5 ; Number of Events ",  10000, -2 , 2);

  Double_t RealElectrons = 0;
  Double_t FakeElectrons = 0;
  Double_t RealElectronPassDetIso025 = 0;
  Double_t FakeElectronPassDetIso025 = 0;
  Double_t RealElectronPassDetIso015 = 0;
  Double_t FakeElectronPassDetIso015 = 0;
  Double_t RealElectronPassDetIso010 = 0;
  Double_t FakeElectronPassDetIso010 = 0;
  Double_t RealElectronPassDetIso008 = 0;
  Double_t FakeElectronPassDetIso008 = 0;
  Double_t RealElectronPassDetIso005 = 0;
  Double_t FakeElectronPassDetIso005 = 0;

  Double_t RealElectronPassMVAIsoLoose = 0;
  Double_t FakeElectronPassMVAIsoLoose = 0;
  Double_t RealElectronPassMVAIsoTight = 0;
  Double_t FakeElectronPassMVAIsoTight = 0;

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

  //Variables
  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Bool_t                  fEleEventNumberParity;
  Float_t                 fElePt; 
  Float_t                 fEleEta; 
  Float_t                 fElePhi; 
  Float_t                 fEleSCEt; 
  Float_t                 fEleSCEta; 
  Float_t                 fEleSCPhi; 
  Bool_t                  fEleIsEcalDriven;
  Float_t                 fElePFIso; 

  //Conversion Variables
  Bool_t                  fEleMatchedConversion;
  Float_t                 fEleConvDist;
  Float_t                 fEleConvDCot;
  UInt_t                  fEleNMissHits;

  //CutBased Variables
  Float_t                 fEleSigmaIEtaIEta; 
  Float_t                 fEleDEtaIn; 
  Float_t                 fEleDPhiIn; 
  Float_t                 fEleHoverE; 
  Float_t                 fEleD0; 
  Float_t                 fEleDZ; 
  Float_t                 fEleFBrem; 
  Float_t                 fEleEOverP; 

  //Additional Vars used in Likelihood
  Float_t                 fEleESeedClusterOverPout; 
  Float_t                 fEleSigmaIPhiIPhi; 
  Float_t                 fEleNBrem; 
  Float_t                 fEleOneOverEMinusOneOverP; 
  Float_t                 fEleESeedClusterOverPIn; 
  Float_t                 fEleIP3d; 
  Float_t                 fEleIP3dSig; 

  Float_t                 fEleHcalDepth1OverEcal;
  Float_t                 fEleHcalDepth2OverEcal;
  Float_t                 fEledEtaCalo;
  Float_t                 fEledPhiCalo;
  Float_t                 fElePreShowerOverRaw;
  Float_t                 fEleSigmaIEtaIPhi;
  Float_t                 fEleSCEtaWidth;
  Float_t                 fEleSCPhiWidth;
  Float_t                 fEleGsfTrackChi2OverNdof;
  UInt_t                  fEleKFTrackNHits;
  Float_t                 fEleKFTrackChi2OverNDoF;
  Float_t                 fEleR9;

  Float_t                 fEleSeedEMaxOverE;
  Float_t                 fEleSeedETopOverE;
  Float_t                 fEleSeedEBottomOverE;
  Float_t                 fEleSeedELeftOverE;
  Float_t                 fEleSeedERightOverE;
  Float_t                 fEleSeedE2ndOverE;
  Float_t                 fEleSeedE2x5RightOverE;
  Float_t                 fEleSeedE2x5LeftOverE;
  Float_t                 fEleSeedE2x5TopOverE;
  Float_t                 fEleSeedE2x5BottomOverE;
  Float_t                 fEleSeedE2x5MaxOverE;
  Float_t                 fEleSeedE1x3OverE;
  Float_t                 fEleSeedE3x1OverE;
  Float_t                 fEleSeedE1x5OverE;
  Float_t                 fEleSeedE2x2OverE;
  Float_t                 fEleSeedE3x2OverE;
  Float_t                 fEleSeedE3x3OverE;
  Float_t                 fEleSeedE4x4OverE;
  Float_t                 fEleSeedE5x5OverE;
  Float_t                 fEleE1x5OverE5x5;
  Float_t                 fEleE3x3OverE5x5;
  Float_t                 fEleE2x5MaxOverE5x5;


  //Isolation Variables
  Float_t                 fEleStandardLikelihood;
  Float_t                 fElePFMVA;
  Float_t                 fEleChargedIso03; 
  Float_t                 fEleNeutralHadronIso03; 
  Float_t                 fEleGammaIso03; 
  Float_t                 fEleChargedIso04; 
  Float_t                 fEleNeutralHadronIso04; 
  Float_t                 fEleGammaIso04; 
  Float_t                 fEleChargedIso04FromOtherVertices; 
  Float_t                 fEleNeutralHadronIso04_10Threshold; 
  Float_t                 fEleGammaIso04_10Threshold; 
  Float_t                 fEleTrkIso03; 
  Float_t                 fEleEMIso03; 
  Float_t                 fEleHadIso03; 
  Float_t                 fEleTrkIso04; 
  Float_t                 fEleEMIso04; 
  Float_t                 fEleHadIso04; 
  Float_t                 fRho; 
  UInt_t                  fNVertices; 

  UInt_t                  fEleTriggerBit;
  Bool_t                  fElePassDenominator;
  Bool_t                  fElePassDenominatorSmurf;

  //More Isolation Variables
  Float_t fChargedIso_DR0p0To0p1;
  Float_t fChargedIso_DR0p1To0p2;
  Float_t fChargedIso_DR0p2To0p3;
  Float_t fChargedIso_DR0p3To0p4;
  Float_t fChargedIso_DR0p4To0p5;
  Float_t fGammaIso_DR0p0To0p1;
  Float_t fGammaIso_DR0p1To0p2;
  Float_t fGammaIso_DR0p2To0p3;
  Float_t fGammaIso_DR0p3To0p4;
  Float_t fGammaIso_DR0p4To0p5;
  Float_t fNeutralHadronIso_DR0p0To0p1;
  Float_t fNeutralHadronIso_DR0p1To0p2;
  Float_t fNeutralHadronIso_DR0p2To0p3;
  Float_t fNeutralHadronIso_DR0p3To0p4;
  Float_t fNeutralHadronIso_DR0p4To0p5;

  Float_t                 fEleIsoMVA_BDTG_V0; 
  Float_t                 fEleIsoMVA_BDTG_V1; 
  Float_t                 fEleIsoMVA_BDTG_V2; 
  Float_t                 fEleIsoMVA_BDTG_V3; 
  Float_t                 fEleIsoMVA_BDTG_V4; 
  Float_t                 fEleIsoMVA_BDTG_V5; 



  //*****************************************************************************************
  //RealEleTree
  //*****************************************************************************************
  TFile *RealEleFile = new TFile(RealElectronFile.c_str(), "READ");
  TTree *RealEleTree = (TTree*)RealEleFile->Get("Electrons");

  RealEleTree->SetBranchAddress( "weight", &fWeight);
  RealEleTree->SetBranchAddress( "run", &fRunNumber);
  RealEleTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  RealEleTree->SetBranchAddress( "event", &fEventNumber);
  RealEleTree->SetBranchAddress( "pt", &fElePt); 
  RealEleTree->SetBranchAddress( "eta", &fEleEta); 
  RealEleTree->SetBranchAddress( "phi", &fElePhi); 
  RealEleTree->SetBranchAddress( "scet", &fEleSCEt); 
  RealEleTree->SetBranchAddress( "sceta", &fEleSCEta); 
  RealEleTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  RealEleTree->SetBranchAddress( "combPFIsoHWW", &fElePFIso); 
  RealEleTree->SetBranchAddress( "see", &fEleSigmaIEtaIEta); 
  RealEleTree->SetBranchAddress( "deta", &fEleDEtaIn); 
  RealEleTree->SetBranchAddress( "dphi", &fEleDPhiIn); 
  RealEleTree->SetBranchAddress( "HoE", &fEleHoverE); 
  RealEleTree->SetBranchAddress( "d0", &fEleD0); 
  RealEleTree->SetBranchAddress( "dZ", &fEleDZ); 
  RealEleTree->SetBranchAddress( "fbrem", &fEleFBrem); 
  RealEleTree->SetBranchAddress( "EoP", &fEleEOverP); 
  RealEleTree->SetBranchAddress( "EoPout", &fEleESeedClusterOverPout); 
  RealEleTree->SetBranchAddress( "spp", &fEleSigmaIPhiIPhi); 
  RealEleTree->SetBranchAddress( "nbrem", &fEleNBrem); 
  RealEleTree->SetBranchAddress( "IoEmIoP", &fEleOneOverEMinusOneOverP); 
  RealEleTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  RealEleTree->SetBranchAddress( "ip3d", &fEleIP3d); 
  RealEleTree->SetBranchAddress( "ip3ds", &fEleIP3dSig); 
  RealEleTree->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  RealEleTree->SetBranchAddress( "chaPFIso", &fEleChargedIso04); 
  RealEleTree->SetBranchAddress( "neuPFIso", &fEleNeutralHadronIso04); 
  RealEleTree->SetBranchAddress( "phoPFIso", &fEleGammaIso04); 
  RealEleTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  RealEleTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  RealEleTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  RealEleTree->SetBranchAddress( "rho", &fRho); 
  RealEleTree->SetBranchAddress( "vertices", &fNVertices); 

  RealEleTree->SetBranchAddress("trkIso03",&fEleTrkIso03);
  RealEleTree->SetBranchAddress("ecalIso03",&fEleEMIso03);
  RealEleTree->SetBranchAddress("hcalIso03",&fEleHadIso03);
  RealEleTree->SetBranchAddress("ChargedIso_DR0p0To0p1",&fChargedIso_DR0p0To0p1);
  RealEleTree->SetBranchAddress("ChargedIso_DR0p1To0p2",&fChargedIso_DR0p1To0p2);
  RealEleTree->SetBranchAddress("ChargedIso_DR0p2To0p3",&fChargedIso_DR0p2To0p3);
  RealEleTree->SetBranchAddress("ChargedIso_DR0p3To0p4",&fChargedIso_DR0p3To0p4);
  RealEleTree->SetBranchAddress("ChargedIso_DR0p4To0p5",&fChargedIso_DR0p4To0p5);
  RealEleTree->SetBranchAddress("GammaIso_DR0p0To0p1",&fGammaIso_DR0p0To0p1);
  RealEleTree->SetBranchAddress("GammaIso_DR0p1To0p2",&fGammaIso_DR0p1To0p2);
  RealEleTree->SetBranchAddress("GammaIso_DR0p2To0p3",&fGammaIso_DR0p2To0p3);
  RealEleTree->SetBranchAddress("GammaIso_DR0p3To0p4",&fGammaIso_DR0p3To0p4);
  RealEleTree->SetBranchAddress("GammaIso_DR0p4To0p5",&fGammaIso_DR0p4To0p5);
  RealEleTree->SetBranchAddress("NeutralHadronIso_DR0p0To0p1",&fNeutralHadronIso_DR0p0To0p1);
  RealEleTree->SetBranchAddress("NeutralHadronIso_DR0p1To0p2",&fNeutralHadronIso_DR0p1To0p2);
  RealEleTree->SetBranchAddress("NeutralHadronIso_DR0p2To0p3",&fNeutralHadronIso_DR0p2To0p3);
  RealEleTree->SetBranchAddress("NeutralHadronIso_DR0p3To0p4",&fNeutralHadronIso_DR0p3To0p4);
  RealEleTree->SetBranchAddress("NeutralHadronIso_DR0p4To0p5",&fNeutralHadronIso_DR0p4To0p5);

  RealEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V0", &fEleIsoMVA_BDTG_V0); 
  RealEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V1", &fEleIsoMVA_BDTG_V1); 
  RealEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V2", &fEleIsoMVA_BDTG_V2); 
  RealEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V3", &fEleIsoMVA_BDTG_V3); 
  RealEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V4", &fEleIsoMVA_BDTG_V4); 
  RealEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V5", &fEleIsoMVA_BDTG_V5); 


  for(UInt_t ientry=0; ientry < RealEleTree->GetEntries(); ientry++) {       	
    RealEleTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
        
    Double_t rho = 0;
    if (!(TMath::IsNaN(fRho) || isinf(fRho))) rho = fRho;

    //don't evaluate performance using training events
//     if (fEventNumber % 2 == 0) continue;

     if (fElePt < 5) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fEleEta) < 1.479) subdet = 0;
    else subdet = 1;
    Int_t ptBin = 0;
    if (fElePt > 10.0) ptBin = 1;
//    if (fElePt > 20.0) ptBin = 2;
//     if (fElePt > 35.0) ptBin = 3;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 3) passCuts = (subdet == 1 && ptBin == 1);
    if (!passCuts) continue;    

    if (NVtxBin == 0) if (!(fNVertices >= 0 && fNVertices <= 5)) continue;
    if (NVtxBin == 1) if (!(fNVertices >= 6 && fNVertices <= 10)) continue;
    if (NVtxBin == 2) if (!(fNVertices >= 11 && fNVertices <= 15)) continue;
    if (NVtxBin == 3) if (!(fNVertices >= 16 && fNVertices <= 20)) continue;
    if (NVtxBin == 4) if (!(fNVertices >= 21 && fNVertices <= 25)) continue;

    //apply some electron preselection
//     if( fabs(fEleSCEta) < 1.479 ) {
//       if (! (fEleSigmaIEtaIEta < 0.02 ) ) continue;
//     } else {
//       if (! (fEleSigmaIEtaIEta < 0.04 ) ) continue;
//     }    
    if( fEleIP3dSig > 4 ) continue;
    if( fEleTrkIso03/fElePt > 0.7 ) continue;

    RealElectrons += fWeight;
  
    Double_t DetIso03PUCorrection = 0;
    if (fabs(fEleEta) < 1.5) DetIso03PUCorrection = rho * (0.078 + 0.026);
    else DetIso03PUCorrection = rho * (0.046 + 0.072);
    Double_t PFIso03PUCorrection = rho * (ElectronEffectiveArea(kEleGammaIsoDR0p0To0p2,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleGammaIsoDR0p2To0p3,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleNeutralHadronIsoDR0p0To0p2,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleNeutralHadronIsoDR0p2To0p3,fEleSCEta,kFALSE));
    Double_t PFIso04PUCorrection = rho * (ElectronEffectiveArea(kEleGammaIsoDR0p0To0p2,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleGammaIsoDR0p2To0p3,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleGammaIsoDR0p3To0p4,fEleSCEta,kFALSE)+ ElectronEffectiveArea(kEleNeutralHadronIsoDR0p0To0p2,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleNeutralHadronIsoDR0p2To0p3,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleNeutralHadronIsoDR0p3To0p4,fEleSCEta,kFALSE));

    if ( (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.25)
      RealElectronPassDetIso025 += fWeight;
    if ( (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.15)
      RealElectronPassDetIso015 += fWeight;
    if ( (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.10)
      RealElectronPassDetIso010 += fWeight;
    if ( (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.08)
      RealElectronPassDetIso008 += fWeight;
    if ( (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.05)
      RealElectronPassDetIso005 += fWeight;
    if (passEleIsoMVASameAsDetIso025( fElePt, fEleEta, fEleIsoMVA_BDTG_V0))
      RealElectronPassMVAIsoLoose += fWeight;
    if (passEleIsoMVASameAsDetIso015( fElePt, fEleEta, fEleIsoMVA_BDTG_V0))
      RealElectronPassMVAIsoTight += fWeight;

    //Fill Histograms
    EleDetIso03_Real->Fill((TMath::Max(fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection,0.0))/fElePt,fWeight);
   ElePFIsoHWW_Real->Fill(fElePFIso,fWeight);
    ElePFIso03_Real->Fill(fChargedIso_DR0p0To0p1 + fChargedIso_DR0p1To0p2 + fChargedIso_DR0p2To0p3 + fGammaIso_DR0p0To0p1 + fGammaIso_DR0p1To0p2 + fGammaIso_DR0p2To0p3 + fNeutralHadronIso_DR0p0To0p1 + fNeutralHadronIso_DR0p1To0p2+ fNeutralHadronIso_DR0p2To0p3,fWeight);
    ElePFIso04_Real->Fill(fChargedIso_DR0p0To0p1 + fChargedIso_DR0p1To0p2 + fChargedIso_DR0p2To0p3 + fChargedIso_DR0p3To0p4 + fGammaIso_DR0p0To0p1 + fGammaIso_DR0p1To0p2 + fGammaIso_DR0p2To0p3 + fGammaIso_DR0p3To0p4 + fNeutralHadronIso_DR0p0To0p1 + fNeutralHadronIso_DR0p1To0p2 + fNeutralHadronIso_DR0p2To0p3 + fNeutralHadronIso_DR0p3To0p4,fWeight);


    EleIsoBDT_V0_Real->Fill(fEleIsoMVA_BDTG_V0,fWeight);
    EleIsoBDT_V1_Real->Fill(fEleIsoMVA_BDTG_V1,fWeight);
    EleIsoBDT_V2_Real->Fill(fEleIsoMVA_BDTG_V2,fWeight);
    EleIsoBDT_V3_Real->Fill(fEleIsoMVA_BDTG_V3,fWeight);
    EleIsoBDT_V4_Real->Fill(fEleIsoMVA_BDTG_V4,fWeight);
    EleIsoBDT_V5_Real->Fill(fEleIsoMVA_BDTG_V5,fWeight);
    
  }
  





  //*****************************************************************************************
  //FakeEleTree
  //*****************************************************************************************
  TFile *FakeEleFile = new TFile(FakeElectronFile.c_str(), "READ");
  TTree *FakeEleTree = (TTree*)FakeEleFile->Get("Electrons");
  FakeEleTree->SetBranchAddress( "weight", &fWeight);
  FakeEleTree->SetBranchAddress( "run", &fRunNumber);
  FakeEleTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  FakeEleTree->SetBranchAddress( "event", &fEventNumber);
  FakeEleTree->SetBranchAddress( "pt", &fElePt); 
  FakeEleTree->SetBranchAddress( "eta", &fEleEta); 
  FakeEleTree->SetBranchAddress( "phi", &fElePhi); 
  FakeEleTree->SetBranchAddress( "scet", &fEleSCEt); 
  FakeEleTree->SetBranchAddress( "sceta", &fEleSCEta); 
  FakeEleTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  FakeEleTree->SetBranchAddress( "combPFIsoHWW", &fElePFIso); 
  FakeEleTree->SetBranchAddress( "see", &fEleSigmaIEtaIEta); 
  FakeEleTree->SetBranchAddress( "deta", &fEleDEtaIn); 
  FakeEleTree->SetBranchAddress( "dphi", &fEleDPhiIn); 
  FakeEleTree->SetBranchAddress( "HoE", &fEleHoverE); 
  FakeEleTree->SetBranchAddress( "d0", &fEleD0); 
  FakeEleTree->SetBranchAddress( "dZ", &fEleDZ); 
  FakeEleTree->SetBranchAddress( "fbrem", &fEleFBrem); 
  FakeEleTree->SetBranchAddress( "EoP", &fEleEOverP); 
  FakeEleTree->SetBranchAddress( "EoPout", &fEleESeedClusterOverPout); 
  FakeEleTree->SetBranchAddress( "spp", &fEleSigmaIPhiIPhi); 
  FakeEleTree->SetBranchAddress( "nbrem", &fEleNBrem); 
  FakeEleTree->SetBranchAddress( "IoEmIoP", &fEleOneOverEMinusOneOverP); 
  FakeEleTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  FakeEleTree->SetBranchAddress( "ip3d", &fEleIP3d); 
  FakeEleTree->SetBranchAddress( "ip3ds", &fEleIP3dSig); 
  FakeEleTree->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  FakeEleTree->SetBranchAddress( "chaPFIso", &fEleChargedIso04); 
  FakeEleTree->SetBranchAddress( "neuPFIso", &fEleNeutralHadronIso04); 
  FakeEleTree->SetBranchAddress( "phoPFIso", &fEleGammaIso04); 
  FakeEleTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  FakeEleTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  FakeEleTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  FakeEleTree->SetBranchAddress( "rho", &fRho); 
  FakeEleTree->SetBranchAddress( "vertices", &fNVertices); 

  FakeEleTree->SetBranchAddress("trkIso03",&fEleTrkIso03);
  FakeEleTree->SetBranchAddress("ecalIso03",&fEleEMIso03);
  FakeEleTree->SetBranchAddress("hcalIso03",&fEleHadIso03);
  FakeEleTree->SetBranchAddress("ChargedIso_DR0p0To0p1",&fChargedIso_DR0p0To0p1);
  FakeEleTree->SetBranchAddress("ChargedIso_DR0p1To0p2",&fChargedIso_DR0p1To0p2);
  FakeEleTree->SetBranchAddress("ChargedIso_DR0p2To0p3",&fChargedIso_DR0p2To0p3);
  FakeEleTree->SetBranchAddress("ChargedIso_DR0p3To0p4",&fChargedIso_DR0p3To0p4);
  FakeEleTree->SetBranchAddress("ChargedIso_DR0p4To0p5",&fChargedIso_DR0p4To0p5);
  FakeEleTree->SetBranchAddress("GammaIso_DR0p0To0p1",&fGammaIso_DR0p0To0p1);
  FakeEleTree->SetBranchAddress("GammaIso_DR0p1To0p2",&fGammaIso_DR0p1To0p2);
  FakeEleTree->SetBranchAddress("GammaIso_DR0p2To0p3",&fGammaIso_DR0p2To0p3);
  FakeEleTree->SetBranchAddress("GammaIso_DR0p3To0p4",&fGammaIso_DR0p3To0p4);
  FakeEleTree->SetBranchAddress("GammaIso_DR0p4To0p5",&fGammaIso_DR0p4To0p5);
  FakeEleTree->SetBranchAddress("NeutralHadronIso_DR0p0To0p1",&fNeutralHadronIso_DR0p0To0p1);
  FakeEleTree->SetBranchAddress("NeutralHadronIso_DR0p1To0p2",&fNeutralHadronIso_DR0p1To0p2);
  FakeEleTree->SetBranchAddress("NeutralHadronIso_DR0p2To0p3",&fNeutralHadronIso_DR0p2To0p3);
  FakeEleTree->SetBranchAddress("NeutralHadronIso_DR0p3To0p4",&fNeutralHadronIso_DR0p3To0p4);
  FakeEleTree->SetBranchAddress("NeutralHadronIso_DR0p4To0p5",&fNeutralHadronIso_DR0p4To0p5);

  FakeEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V0", &fEleIsoMVA_BDTG_V0); 
  FakeEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V1", &fEleIsoMVA_BDTG_V1); 
  FakeEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V2", &fEleIsoMVA_BDTG_V2); 
  FakeEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V3", &fEleIsoMVA_BDTG_V3); 
  FakeEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V4", &fEleIsoMVA_BDTG_V4); 
  FakeEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V5", &fEleIsoMVA_BDTG_V5); 


  for(UInt_t ientry=0; ientry < FakeEleTree->GetEntries(); ientry++) {       	
    FakeEleTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    Double_t rho = 0;
    if (!(TMath::IsNaN(fRho) || isinf(fRho))) rho = fRho;

    //don't evaluate performance using training events
    if (fEventNumber % 2 == 0) continue;

     if (fElePt < 5) continue;

   //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fEleEta) < 1.479) subdet = 0;
    else subdet = 1;
    Int_t ptBin = 0;
    if (fElePt > 10.0) ptBin = 1;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 3) passCuts = (subdet == 1 && ptBin == 1);
    if (!passCuts) continue;    

    if (NVtxBin == 0) if (!(fNVertices >= 0 && fNVertices <= 5)) continue;
    if (NVtxBin == 1) if (!(fNVertices >= 6 && fNVertices <= 10)) continue;
    if (NVtxBin == 2) if (!(fNVertices >= 11 && fNVertices <= 15)) continue;
    if (NVtxBin == 3) if (!(fNVertices >= 16 && fNVertices <= 20)) continue;
    if (NVtxBin == 4) if (!(fNVertices >= 21 && fNVertices <= 25)) continue;

    //apply some electron preselection
//     if( fabs(fEleSCEta) < 1.479 ) {
//       if (! (fEleSigmaIEtaIEta < 0.02 ) ) continue;
//     } else {
//       if (! (fEleSigmaIEtaIEta < 0.04 ) ) continue;
//     }    
    if( fEleIP3dSig > 4 ) continue;
    if( fEleTrkIso03/fElePt > 0.7 ) continue;


    FakeElectrons += fWeight;

    Double_t DetIso03PUCorrection = 0;
    if (fabs(fEleEta) < 1.5) DetIso03PUCorrection = rho * (0.078 + 0.026);
    else DetIso03PUCorrection = rho * (0.046 + 0.072);
    Double_t PFIso03PUCorrection = rho * (ElectronEffectiveArea(kEleGammaIsoDR0p0To0p2,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleGammaIsoDR0p2To0p3,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleNeutralHadronIsoDR0p0To0p2,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleNeutralHadronIsoDR0p2To0p3,fEleSCEta,kFALSE));
    Double_t PFIso04PUCorrection = rho * (ElectronEffectiveArea(kEleGammaIsoDR0p0To0p2,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleGammaIsoDR0p2To0p3,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleGammaIsoDR0p3To0p4,fEleSCEta,kFALSE)+ ElectronEffectiveArea(kEleNeutralHadronIsoDR0p0To0p2,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleNeutralHadronIsoDR0p2To0p3,fEleSCEta,kFALSE) + ElectronEffectiveArea(kEleNeutralHadronIsoDR0p3To0p4,fEleSCEta,kFALSE));

    if ( (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.25)
      FakeElectronPassDetIso025 += fWeight;
    if ( (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.15)
      FakeElectronPassDetIso015 += fWeight;
    if ( (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.10)
      FakeElectronPassDetIso010 += fWeight;
    if ( (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.08)
      FakeElectronPassDetIso008 += fWeight;
    if ( (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.05)
      FakeElectronPassDetIso005 += fWeight;
    if (passEleIsoMVASameAsDetIso025( fElePt, fEleEta, fEleIsoMVA_BDTG_V0))
      FakeElectronPassMVAIsoLoose += fWeight;
    if (passEleIsoMVASameAsDetIso015( fElePt, fEleEta, fEleIsoMVA_BDTG_V0))
      FakeElectronPassMVAIsoTight += fWeight;

    //Fill Histograms
    EleDetIso03_Fake->Fill((TMath::Max(fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection,0.0))/fElePt,fWeight);
    ElePFIsoHWW_Fake->Fill(fElePFIso,fWeight);
    ElePFIso03_Fake->Fill(fChargedIso_DR0p0To0p1 + fChargedIso_DR0p1To0p2 + fChargedIso_DR0p2To0p3 + fGammaIso_DR0p0To0p1 + fGammaIso_DR0p1To0p2 + fGammaIso_DR0p2To0p3 + fNeutralHadronIso_DR0p0To0p1 + fNeutralHadronIso_DR0p1To0p2+ fNeutralHadronIso_DR0p2To0p3,fWeight);
    ElePFIso04_Fake->Fill(fChargedIso_DR0p0To0p1 + fChargedIso_DR0p1To0p2 + fChargedIso_DR0p2To0p3 + fChargedIso_DR0p3To0p4 + fGammaIso_DR0p0To0p1 + fGammaIso_DR0p1To0p2 + fGammaIso_DR0p2To0p3 + fGammaIso_DR0p3To0p4 + fNeutralHadronIso_DR0p0To0p1 + fNeutralHadronIso_DR0p1To0p2 + fNeutralHadronIso_DR0p2To0p3 + fNeutralHadronIso_DR0p3To0p4,fWeight);
 
    //Fill Histograms
    EleIsoBDT_V0_Fake->Fill(fEleIsoMVA_BDTG_V0,fWeight);
    EleIsoBDT_V1_Fake->Fill(fEleIsoMVA_BDTG_V1,fWeight);
    EleIsoBDT_V2_Fake->Fill(fEleIsoMVA_BDTG_V2,fWeight);
    EleIsoBDT_V3_Fake->Fill(fEleIsoMVA_BDTG_V3,fWeight);
    EleIsoBDT_V4_Fake->Fill(fEleIsoMVA_BDTG_V4,fWeight);
    EleIsoBDT_V5_Fake->Fill(fEleIsoMVA_BDTG_V5,fWeight);

  } //loop over electrons
  



  
  //*****************************************************************************************
  //Current Working Points
  //*****************************************************************************************
  cout << "DetIso < 0.25 : Real Electron Efficiency : " << RealElectronPassDetIso025 << " / " << RealElectrons << " = " << RealElectronPassDetIso025/RealElectrons << endl;
  cout << "DetIso < 0.25 : Fake Electron Efficiency : " << FakeElectronPassDetIso025 << " / " << FakeElectrons << " = " << FakeElectronPassDetIso025/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_DetIso025 = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassDetIso025/RealElectrons , FakeElectronPassDetIso025/FakeElectrons, "ROC_DetIso025"+label);

  cout << "DetIso < 0.15 : Real Electron Efficiency : " << RealElectronPassDetIso015 << " / " << RealElectrons << " = " << RealElectronPassDetIso015/RealElectrons << endl;
  cout << "DetIso < 0.15 : Fake Electron Efficiency : " << FakeElectronPassDetIso015 << " / " << FakeElectrons << " = " << FakeElectronPassDetIso015/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_DetIso015 = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassDetIso015/RealElectrons , FakeElectronPassDetIso015/FakeElectrons, "ROC_DetIso015"+label);

  cout << "DetIso < 0.10 : Real Electron Efficiency : " << RealElectronPassDetIso010 << " / " << RealElectrons << " = " << RealElectronPassDetIso010/RealElectrons << endl;
  cout << "DetIso < 0.10 : Fake Electron Efficiency : " << FakeElectronPassDetIso010 << " / " << FakeElectrons << " = " << FakeElectronPassDetIso010/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_DetIso010 = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassDetIso010/RealElectrons , FakeElectronPassDetIso010/FakeElectrons, "ROC_DetIso010"+label);

  cout << "DetIso < 0.08 : Real Electron Efficiency : " << RealElectronPassDetIso008 << " / " << RealElectrons << " = " << RealElectronPassDetIso008/RealElectrons << endl;
  cout << "DetIso < 0.08 : Fake Electron Efficiency : " << FakeElectronPassDetIso008 << " / " << FakeElectrons << " = " << FakeElectronPassDetIso008/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_DetIso008 = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassDetIso008/RealElectrons , FakeElectronPassDetIso008/FakeElectrons, "ROC_DetIso008"+label);

  cout << "DetIso < 0.05 : Real Electron Efficiency : " << RealElectronPassDetIso005 << " / " << RealElectrons << " = " << RealElectronPassDetIso005/RealElectrons << endl;
  cout << "DetIso < 0.05 : Fake Electron Efficiency : " << FakeElectronPassDetIso005 << " / " << FakeElectrons << " = " << FakeElectronPassDetIso005/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_DetIso005 = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassDetIso005/RealElectrons , FakeElectronPassDetIso005/FakeElectrons, "ROC_DetIso005"+label);

  cout << "MVA Loose : Real Electron Efficiency : " << RealElectronPassMVAIsoLoose << " / " << RealElectrons << " = " << RealElectronPassMVAIsoLoose/RealElectrons << endl;
  cout << "MVA Tight : Fake Electron Efficiency : " << FakeElectronPassMVAIsoTight << " / " << FakeElectrons << " = " << FakeElectronPassMVAIsoTight/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_MVAIsoLoose = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassMVAIsoLoose/RealElectrons , FakeElectronPassMVAIsoLoose/FakeElectrons, "ROC_MVALoose"+label);
  TGraphAsymmErrors* ROC_MVAIsoTight = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassMVAIsoTight/RealElectrons , FakeElectronPassMVAIsoTight/FakeElectrons, "ROC_MVATight"+label);

  cout << "**********************\n";
  cout << "Bkg At DetIso < 0.15 Signal Eff\n";

  Double_t BkgEffCutBased = FakeElectronPassDetIso025/FakeElectrons;
  Double_t SigEffCutBased = RealElectronPassDetIso025/RealElectrons;
  Double_t SigEffPFIso03_SameBkg = FindSigEffAtFixedBkgEfficiency(ElePFIso03_Real, ElePFIso03_Fake, BkgEffCutBased, kFALSE);
  Double_t SigEffPFIso04_SameBkg = FindSigEffAtFixedBkgEfficiency(ElePFIso04_Real, ElePFIso04_Fake, BkgEffCutBased, kFALSE);
  Double_t SigEffBDT_V0_SameBkg = FindSigEffAtFixedBkgEfficiency(EleIsoBDT_V0_Real, EleIsoBDT_V0_Fake, BkgEffCutBased);
  Double_t SigEffBDT_V1_SameBkg = FindSigEffAtFixedBkgEfficiency(EleIsoBDT_V1_Real, EleIsoBDT_V1_Fake, BkgEffCutBased);
  Double_t SigEffBDT_V2_SameBkg = FindSigEffAtFixedBkgEfficiency(EleIsoBDT_V2_Real, EleIsoBDT_V2_Fake, BkgEffCutBased);
  Double_t SigEffBDT_V3_SameBkg = FindSigEffAtFixedBkgEfficiency(EleIsoBDT_V3_Real, EleIsoBDT_V3_Fake, BkgEffCutBased);
  Double_t SigEffBDT_V4_SameBkg = FindSigEffAtFixedBkgEfficiency(EleIsoBDT_V4_Real, EleIsoBDT_V4_Fake, BkgEffCutBased);
  Double_t SigEffBDT_V5_SameBkg = FindSigEffAtFixedBkgEfficiency(EleIsoBDT_V5_Real, EleIsoBDT_V5_Fake, BkgEffCutBased);

  cout << "Signal Efficiency (wrt Cut-based) for : same bkg \n";
  cout << "PFIso03 : " << SigEffPFIso03_SameBkg/SigEffCutBased << endl;
  cout << "PFIso04 : " << SigEffPFIso04_SameBkg/SigEffCutBased << endl;
  cout << "BDT V0 : " << SigEffBDT_V0_SameBkg/SigEffCutBased << endl;
  cout << "BDT V1 : " << SigEffBDT_V1_SameBkg/SigEffCutBased << endl;
  cout << "BDT V2 : " << SigEffBDT_V2_SameBkg/SigEffCutBased << endl;
  cout << "BDT V3 : " << SigEffBDT_V3_SameBkg/SigEffCutBased << endl;
  cout << "BDT V4 : " << SigEffBDT_V4_SameBkg/SigEffCutBased << endl;
  cout << "BDT V5 : " << SigEffBDT_V5_SameBkg/SigEffCutBased << endl;
  cout << "**********************\n";

  Double_t BkgEffPFIso03_SameSig = FindBkgEffAtFixedSignalEfficiency(ElePFIso03_Real, ElePFIso03_Fake, SigEffCutBased, kFALSE);
  Double_t BkgEffPFIso04_SameSig = FindBkgEffAtFixedSignalEfficiency(ElePFIso04_Real, ElePFIso04_Fake, SigEffCutBased, kFALSE);
  Double_t BkgEffBDT_V0_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIsoBDT_V0_Real, EleIsoBDT_V0_Fake, SigEffCutBased);
  Double_t BkgEffBDT_V1_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIsoBDT_V1_Real, EleIsoBDT_V1_Fake, SigEffCutBased);
  Double_t BkgEffBDT_V2_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIsoBDT_V2_Real, EleIsoBDT_V2_Fake, SigEffCutBased);
  Double_t BkgEffBDT_V3_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIsoBDT_V3_Real, EleIsoBDT_V3_Fake, SigEffCutBased);
  Double_t BkgEffBDT_V4_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIsoBDT_V4_Real, EleIsoBDT_V4_Fake, SigEffCutBased);
  Double_t BkgEffBDT_V5_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIsoBDT_V5_Real, EleIsoBDT_V5_Fake, SigEffCutBased);

  cout << "Bkg Efficiency (wrt Cut-based) for same sig eff \n";
  cout << "PFIso03 : " << BkgEffPFIso03_SameSig/BkgEffCutBased << endl;
  cout << "PFIso04 : " << BkgEffPFIso04_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV0 : " << BkgEffBDT_V0_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV1 : " << BkgEffBDT_V1_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV2 : " << BkgEffBDT_V2_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV3 : " << BkgEffBDT_V3_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV4 : " << BkgEffBDT_V4_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV5 : " << BkgEffBDT_V5_SameSig/BkgEffCutBased << endl;

  cout << "**********************\n";




  //*****************************************************************************************
  //Make ROC curves
  //*****************************************************************************************
  TGraphAsymmErrors* ROC_DetIso03 = MakeSigEffVsBkgEffGraph(EleDetIso03_Real, EleDetIso03_Fake, "ROC_EleDetIso03"+label , kFALSE );
  TGraphAsymmErrors* ROC_PFIso03 = MakeSigEffVsBkgEffGraph(ElePFIso03_Real, ElePFIso03_Fake, "ROC_ElePFIso03"+label , kFALSE );
  TGraphAsymmErrors* ROC_PFIso04 = MakeSigEffVsBkgEffGraph(ElePFIso04_Real, ElePFIso04_Fake, "ROC_ElePFIso04"+label , kFALSE );
  //TGraphAsymmErrors* ROC_PFIsoHWW = MakeSigEffVsBkgEffGraph(ElePFIsoHWW_Real, ElePFIsoHWW_Fake, "ROC_ElePFIsoHWW"+label , kFALSE );

  TGraphAsymmErrors* ROC_BDTG_V0 = MakeSigEffVsBkgEffGraph(EleIsoBDT_V0_Real, EleIsoBDT_V0_Fake, "ROC_BDTGV0"+label );
  TGraphAsymmErrors* ROC_BDTG_V1 = MakeSigEffVsBkgEffGraph(EleIsoBDT_V1_Real, EleIsoBDT_V1_Fake, "ROC_BDTGV1"+label );
  TGraphAsymmErrors* ROC_BDTG_V2 = MakeSigEffVsBkgEffGraph(EleIsoBDT_V2_Real, EleIsoBDT_V2_Fake, "ROC_BDTGV2"+label );
//   TGraphAsymmErrors* ROC_BDTG_V3 = MakeSigEffVsBkgEffGraph(EleIsoBDT_V3_Real, EleIsoBDT_V3_Fake, "ROC_BDTGV3"+label );
//   TGraphAsymmErrors* ROC_BDTG_V4 = MakeSigEffVsBkgEffGraph(EleIsoBDT_V4_Real, EleIsoBDT_V4_Fake, "ROC_BDTGV4"+label );
//   TGraphAsymmErrors* ROC_BDTG_V5 = MakeSigEffVsBkgEffGraph(EleIsoBDT_V5_Real, EleIsoBDT_V5_Fake, "ROC_BDTGV5"+label );


  //*****************************************************************************************
  //Find Cut with same signal efficiency Make ROC curves
  //*****************************************************************************************
  Double_t CutValue_SameBkgAsDetIso025 = FindCutValueAtFixedEfficiency(EleIsoBDT_V0_Fake, FakeElectronPassDetIso025/FakeElectrons );
  cout << "MVA Iso Cut Value @ Same Bkg as DetIso03 < 0.25: " << CutValue_SameBkgAsDetIso025  << " " << endl;
  Double_t CutValue_SameBkgAsDetIso015 = FindCutValueAtFixedEfficiency(EleIsoBDT_V0_Fake, FakeElectronPassDetIso015/FakeElectrons );
  cout << "MVA Iso Cut Value @ Same Bkg as DetIso03 < 0.15: " << CutValue_SameBkgAsDetIso015 << " " << endl;


//   TFile *canvasFile = new TFile("ElectronIDMVAPerformancePlots.root","UPDATE");
  TLegend* legend;
  TCanvas* cv;
  string plotname;

//   //*****************************************************************************************
//   //Plot Distributions
//   //*****************************************************************************************


//   vector<TH1F*> hists;
//   vector<string> histLabels;




  //*****************************************************************************************
  //Plot ROC Curves
  //*****************************************************************************************
  vector<TGraphAsymmErrors*> ROCGraphs;
  vector<string> GraphLabels;
  vector<Int_t> colors;

//   vector<Int_t> markers;
//   vector<Int_t> colors;
//   colors.push_back(kRed);     markers.push_back(20);
//   colors.push_back(kCyan);    markers.push_back(21);
// //   colors.push_back(kBlue);    markers.push_back(21);
//   colors.push_back(kMagenta); markers.push_back(22);
//   colors.push_back(kCyan);    markers.push_back(34);
//   colors.push_back(kBlack);   markers.push_back(29);
//   colors.push_back(kGreen);   markers.push_back(33);
//   colors.push_back(kRed-2);   markers.push_back(33);
//   colors.push_back(kOrange);   markers.push_back(33);
//   colors.push_back(kBlue-2);   markers.push_back(33);
//   colors.push_back(kGreen-2);   markers.push_back(33);
//   colors.push_back(kMagenta-2);   markers.push_back(33);

  //*****************************************************************************************
  //*****************************************************************************************
  ROCGraphs.clear();
  GraphLabels.clear();
  plotname = "ElectronIDMVA"+label;


  ROCGraphs.push_back(ROC_DetIso03);
  GraphLabels.push_back("DetIso");
  colors.push_back(kBlue);

  ROCGraphs.push_back(ROC_PFIso03);
  GraphLabels.push_back("PFIso");
  colors.push_back(kRed);

//    ROCGraphs.push_back(ROC_PFIso04);
//    GraphLabels.push_back("PFIso04");
//    colors.push_back(kRed);


  ROCGraphs.push_back(ROC_BDTG_V0);
  GraphLabels.push_back("IsoRings");
  colors.push_back(kMagenta);

//   ROCGraphs.push_back(ROC_BDTG_V2);
//   GraphLabels.push_back("W+Jets Trained + cuts");
//   colors.push_back(kOrange);

//   ROCGraphs.push_back(ROC_BDTG_V3);
//   GraphLabels.push_back("IsoRings+CategorizedShapes");
//   colors.push_back(kMagenta);
  
//   ROCGraphs.push_back(ROC_BDTG_V4);
//   GraphLabels.push_back("IsoRings+CatShapes+Directional");
//   colors.push_back(kOrange);
  
//   ROCGraphs.push_back(ROC_BDTG_V5);
//   GraphLabels.push_back("SmallerIsoRings+CatShapes+Directional");
//   colors.push_back(kBlack);
  



  //*****************************************************************************************
  Double_t xmin = 0.0;
  Double_t xmax = 1.0;
  Double_t ymin = 0.0;
  Double_t ymax = 1.0;
// //   if (Option == 0 )                              { xmin = 0.15; xmax = 0.45; ymin = 0.75; ymax = 0.85; }

// //FOr Data
//    if (Option == 0 )                              { xmin = 0.05; xmax = 0.26; ymin = 0.35; ymax = 0.85; }
// //   if (Option == 1 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.00; ymax = 1.00; }
// //   if (Option == 2 )                              { xmin = 0.00; xmax = 0.40; ymin = 0.00; ymax = 1.00; }
// //   if (Option == 3 )                              { xmin = 0.25; xmax = 0.65; ymin = 0.80; ymax = 1.00; }
//   if (Option == 4 )                                 { xmin = 0.02; xmax = 0.28; ymin = 0.75; ymax = 1.00; }
//   if (Option == 5 )                                 { xmin = 0.15; xmax = 0.32; ymin = 0.75; ymax = 1.00; }



//   //For MC
//    if (Option == 0 )                              { xmin = 0.05; xmax = 0.26; ymin = 0.35; ymax = 0.85; }
// //   if (Option == 1 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.00; ymax = 1.00; }
// //   if (Option == 2 )                              { xmin = 0.00; xmax = 0.40; ymin = 0.00; ymax = 1.00; }
// //   if (Option == 3 )                              { xmin = 0.25; xmax = 0.65; ymin = 0.80; ymax = 1.00; }
//   if (Option == 4 )                                 { xmin = 0.02; xmax = 0.23; ymin = 0.75; ymax = 1.00; }
//   if (Option == 5 )                                 { xmin = 0.15; xmax = 0.32; ymin = 0.75; ymax = 1.00; }


//FOr Data
//   if (Option == 0 )                              { xmin = 0.10; xmax = 0.50; ymin = 0.70; ymax = 0.90; }
//   if (Option == 1 )                              { xmin = 0.05; xmax = 0.35; ymin = 0.30; ymax = 0.90; }
//   if (Option == 2 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.30; ymax = 0.80; }
//   if (Option == 3 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.80; ymax = 1.00; }
//   if (Option == 4 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.70; ymax = 1.00; }
//   if (Option == 5 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.70; ymax = 1.00; }

  cv = new TCanvas("cv", "cv", 800, 600);

//    legend = new TLegend(0.45,0.20,0.75,0.50);
  legend = new TLegend(0.54,0.14,0.90,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(ROCGraphs[i],GraphLabels[i].c_str(), "LP");

    ROCGraphs[i]->SetMarkerColor(colors[i]);
    ROCGraphs[i]->SetLineColor(colors[i]);
    ROCGraphs[i]->SetMarkerSize(0.75);
   
    ROCGraphs[i]->GetXaxis()->SetRangeUser(xmin,xmax);    
    ROCGraphs[i]->GetYaxis()->SetRangeUser(ymin,ymax);    
    if (i==0) {
      ROCGraphs[i]->Draw("AP");
    } else {
      ROCGraphs[i]->Draw("Psame");
    }
  }


  legend->AddEntry(ROC_DetIso025, "DetIso < 0.25", "P");
  ROC_DetIso025->SetFillColor(kGreen+1);
  ROC_DetIso025->SetMarkerColor(kGreen+1);
  ROC_DetIso025->SetMarkerStyle(34);
  ROC_DetIso025->SetMarkerSize(2.5);
  ROC_DetIso025->Draw("Psame");

  legend->AddEntry(ROC_DetIso015, "DetIso < 0.15", "P");
  ROC_DetIso015->SetFillColor(kGreen+2);
  ROC_DetIso015->SetMarkerColor(kGreen+2);
  ROC_DetIso015->SetMarkerStyle(34);
  ROC_DetIso015->SetMarkerSize(2.5);
  ROC_DetIso015->Draw("Psame");

//   legend->AddEntry(ROC_DetIso010, "DetIso < 0.1", "P");
//   ROC_DetIso010->SetFillColor(kGreen+3);
//   ROC_DetIso010->SetMarkerColor(kGreen+3);
//   ROC_DetIso010->SetMarkerStyle(34);
//   ROC_DetIso010->SetMarkerSize(2.5);
//   ROC_DetIso010->Draw("Psame");

//   legend->AddEntry(ROC_DetIso008, "DetIso < 0.08", "P");
//   ROC_DetIso008->SetFillColor(kGreen+4);
//   ROC_DetIso008->SetMarkerColor(kGreen+4);
//   ROC_DetIso008->SetMarkerStyle(34);
//   ROC_DetIso008->SetMarkerSize(2.5);
//   ROC_DetIso008->Draw("Psame");

//   legend->AddEntry(ROC_DetIso005, "DetIso < 0.05", "P");
//   ROC_DetIso005->SetFillColor(kGreen);
//   ROC_DetIso005->SetMarkerColor(kGreen);
//   ROC_DetIso005->SetMarkerStyle(34);
//   ROC_DetIso005->SetMarkerSize(2.5);
//   ROC_DetIso005->Draw("Psame");

//   legend->AddEntry(ROC_MVAIsoLoose, "MVA Iso Loose", "P");
//   ROC_MVAIsoLoose->SetFillColor(kBlue+1);
//   ROC_MVAIsoLoose->SetMarkerColor(kBlue+1);
//   ROC_MVAIsoLoose->SetMarkerStyle(34);
//   ROC_MVAIsoLoose->SetMarkerSize(2.5);
//   ROC_MVAIsoLoose->Draw("Psame");

//   legend->AddEntry(ROC_MVAIsoTight, "MVA Iso Tight", "P");
//   ROC_MVAIsoTight->SetFillColor(kBlue+3);
//   ROC_MVAIsoTight->SetMarkerColor(kBlue+3);
//   ROC_MVAIsoTight->SetMarkerStyle(34);
//   ROC_MVAIsoTight->SetMarkerSize(2.5);
//   ROC_MVAIsoTight->Draw("Psame");


  legend->Draw();
  
  cv->SaveAs(("ROCGraphs_" + plotname + ".gif").c_str());
//   canvasFile->WriteTObject(cv,("ROCGraphs_" + plotname).c_str(), "WriteDelete") ;

 





  gBenchmark->Show("WWTemplate");       
} 

