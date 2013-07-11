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
//
//*************************************************************************************************
TGraphAsymmErrors* MakeSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, string name ) {
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

  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    Double_t nbkg = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
      nbkg += bkgHist->GetBinContent(q);
    }

    Double_t ratio;
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    ratio = n1/n2;


    SigEff[b] = ratio;
    SigEffErrLow[b] = 0;
    SigEffErrHigh[b] = 0;

    n1 = TMath::Nint(nbkg);
    n2 = TMath::Nint(NBkgTotal);
    ratio = n1/n2;
    BkgEff[b] = ratio;
    BkgEffErrLow[b] = 0;
    BkgEffErrHigh[b] = 0;
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
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    ratio = n1/n2;
    SigEff[b] = ratio;
    SigEffErrLow[b] = 0;
    SigEffErrHigh[b] = 0;

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
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    ratio = n1/n2;
    
      cout << myCutValue << " : " << signalHist->GetXaxis()->GetBinCenter(b) << " , " << cutValue[0] << endl;
    if (fabs(myCutValue - signalHist->GetXaxis()->GetBinCenter(b)) < fabs(myCutValue - cutValue[0])) {
      SigEff[0] = ratio;
      SigEffErrLow[0] = 0;
      SigEffErrHigh[0] = 0;
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
Double_t FindCutValueAtFixedEfficiency(TH1F* signalHist, Double_t targetSignalEff ) {
  //Make Met Plots


  Double_t targetCutValue = -9999;
  Double_t bestCurrentSignalEff = 0;
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double NSigTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t ratio = nsig / NSigTotal;
//     cout << targetSignalEff << " : " << ratio << " , " << signalHist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetSignalEff - ratio) < fabs(targetSignalEff - bestCurrentSignalEff)) {
      targetCutValue = signalHist->GetXaxis()->GetBinCenter(b);
      bestCurrentSignalEff = ratio;
    }
  }

  return targetCutValue;
}


//*************************************************************************************************
//
//*************************************************************************************************
Double_t FindBkgEffAtFixedSignalEfficiency(TH1F* signalHist, TH1F* bkgHist, Double_t targetSignalEff ) {
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
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t nbkg = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nbkg += bkgHist->GetBinContent(q);
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
Double_t FindSigEffAtFixedBkgEfficiency(TH1F* signalHist, TH1F* bkgHist, Double_t targetBkgEff ) {
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
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t nbkg = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nbkg += bkgHist->GetBinContent(q);
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



Bool_t passICHEP2012( Double_t fElePt, Double_t fEleSCEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fEleSCEta) < 0.8) subdet = 0;
  else if (fabs(fEleSCEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (fElePt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = 0.470;
  if (MVABin == 1) MVACut = 0.004;
  if (MVABin == 2) MVACut = 0.295; 
  if (MVABin == 3) MVACut = 0.500;
  if (MVABin == 4) MVACut = 0.120;
  if (MVABin == 5) MVACut = 0.600;  


  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}


Bool_t passIDMVA( Double_t fElePt, Double_t fEleSCEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fEleSCEta) < 0.8) subdet = 0;
  else if (fabs(fEleSCEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (fElePt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = 0.5094;
  if (MVABin == 1) MVACut = -0.0394;
  if (MVABin == 2) MVACut = 0.1902; 
  if (MVABin == 3) MVACut = 0.839;
  if (MVABin == 4) MVACut = 0.673;
  if (MVABin == 5) MVACut = 0.8078;  


  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}




Bool_t passCutBasedIsoOnly(Double_t fElePt, 
                    Double_t fEleEta, 
                    Double_t fElePFIso ) {

  //apply full isolation cut
  Bool_t passIsoCuts = kFALSE;
  if (fElePt >= 10 && fElePt < 20) passIsoCuts = ( fElePFIso  < 0.09 ); 
  if (fElePt >= 20) passIsoCuts = ( fElePFIso  < 0.13 ); 

  return passIsoCuts;
}



//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void MakeElectronIDMVAPerformancePlots(string RealElectronFile, string FakeElectronFile, string Label, Int_t Option)
{  

  string label = "";
  if (Label != "") label = "_" + Label;


  //*****************************************************************************************
  //Plotting Setup
  //*****************************************************************************************
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


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *EleIDBDTGV0_Real = new TH1F(("EleIDBDTGV0_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV0_Fake = new TH1F(("EleIDBDTGV0_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV1_Real = new TH1F(("EleIDBDTGV1_Real"+label).c_str(), "; BDTG V1 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV1_Fake = new TH1F(("EleIDBDTGV1_Fake"+label).c_str(), "; BDTG V1 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV2_Real = new TH1F(("EleIDBDTGV2_Real"+label).c_str(), "; BDTG V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV2_Fake = new TH1F(("EleIDBDTGV2_Fake"+label).c_str(), "; BDTG V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV3_Real = new TH1F(("EleIDBDTGV3_Real"+label).c_str(), "; BDTG V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV3_Fake = new TH1F(("EleIDBDTGV3_Fake"+label).c_str(), "; BDTG V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV4_Real = new TH1F(("EleIDBDTGV4_Real"+label).c_str(), "; BDTG V4 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV4_Fake = new TH1F(("EleIDBDTGV4_Fake"+label).c_str(), "; BDTG V4 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV5_Real = new TH1F(("EleIDBDTGV5_Real"+label).c_str(), "; BDTG V5 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV5_Fake = new TH1F(("EleIDBDTGV5_Fake"+label).c_str(), "; BDTG V5 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV6_Real = new TH1F(("EleIDBDTGV6_Real"+label).c_str(), "; BDTG V6 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV6_Fake = new TH1F(("EleIDBDTGV6_Fake"+label).c_str(), "; BDTG V6 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV7_Real = new TH1F(("EleIDBDTGV7_Real"+label).c_str(), "; BDTG V7 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV7_Fake = new TH1F(("EleIDBDTGV7_Fake"+label).c_str(), "; BDTG V7 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV8_Real = new TH1F(("EleIDBDTGV8_Real"+label).c_str(), "; BDTG V8 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV8_Fake = new TH1F(("EleIDBDTGV8_Fake"+label).c_str(), "; BDTG V8 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV9_Real = new TH1F(("EleIDBDTGV9_Real"+label).c_str(), "; BDTG V9 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV9_Fake = new TH1F(("EleIDBDTGV9_Fake"+label).c_str(), "; BDTG V9 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV10_Real = new TH1F(("EleIDBDTGV10_Real"+label).c_str(), "; BDTG V10 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV10_Fake = new TH1F(("EleIDBDTGV10_Fake"+label).c_str(), "; BDTG V10 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV11_Real = new TH1F(("EleIDBDTGV11_Real"+label).c_str(), "; BDTG V11 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV11_Fake = new TH1F(("EleIDBDTGV11_Fake"+label).c_str(), "; BDTG V11 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV12_Real = new TH1F(("EleIDBDTGV12_Real"+label).c_str(), "; BDTG V12 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV12_Fake = new TH1F(("EleIDBDTGV12_Fake"+label).c_str(), "; BDTG V12 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV13_Real = new TH1F(("EleIDBDTGV13_Real"+label).c_str(), "; BDTG V13 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV13_Fake = new TH1F(("EleIDBDTGV13_Fake"+label).c_str(), "; BDTG V13 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV14_Real = new TH1F(("EleIDBDTGV14_Real"+label).c_str(), "; BDTG V14 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV14_Fake = new TH1F(("EleIDBDTGV14_Fake"+label).c_str(), "; BDTG V14 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV15_Real = new TH1F(("EleIDBDTGV15_Real"+label).c_str(), "; BDTG V15 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV15_Fake = new TH1F(("EleIDBDTGV15_Fake"+label).c_str(), "; BDTG V15 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV16_Real = new TH1F(("EleIDBDTGV16_Real"+label).c_str(), "; BDTG V16 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV16_Fake = new TH1F(("EleIDBDTGV16_Fake"+label).c_str(), "; BDTG V16 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV17_Real = new TH1F(("EleIDBDTGV17_Real"+label).c_str(), "; BDTG V17 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV17_Fake = new TH1F(("EleIDBDTGV17_Fake"+label).c_str(), "; BDTG V17 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV18_Real = new TH1F(("EleIDBDTGV18_Real"+label).c_str(), "; BDTG V18 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV18_Fake = new TH1F(("EleIDBDTGV18_Fake"+label).c_str(), "; BDTG V18 ; Number of Events ",  10000, -2 , 2);

  Double_t RealElectrons = 0;
  Double_t FakeElectrons = 0;
  Double_t RealElectronPassIDMVA = 0;
  Double_t FakeElectronPassIDMVA = 0;
  Double_t RealElectronPassICHEP2012 = 0;
  Double_t FakeElectronPassICHEP2012 = 0;

  Double_t RealElectronPassBDTGV1_SameCutBasedSig = 0;
  Double_t FakeElectronPassBDTGV1_SameCutBasedSig = 0;
  Double_t RealElectronPassBDTGV1_HalfCutBasedBkg = 0;
  Double_t FakeElectronPassBDTGV1_HalfCutBasedBkg = 0;
  Double_t RealElectronPassBDTGV1_OneThirdCutBasedBkg = 0;
  Double_t FakeElectronPassBDTGV1_OneThirdCutBasedBkg = 0;

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
    if (RealEleTree.fEventNumber % 2 == 0) continue;
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

    RealElectrons += RealEleTree.fWeight;

    if (passICHEP2012( RealEleTree.fElePt, RealEleTree.fEleEta, fEleBDTGV0))
      RealElectronPassICHEP2012 += RealEleTree.fWeight;
    if (passIDMVA( RealEleTree.fElePt, RealEleTree.fEleEta, fEleBDTGV0))
      RealElectronPassIDMVA += RealEleTree.fWeight;
    
    EleIDBDTGV0_Real->Fill(fEleBDTGV0,RealEleTree.fWeight);
    EleIDBDTGV1_Real->Fill(fEleBDTGV1,RealEleTree.fWeight);
    EleIDBDTGV2_Real->Fill(fEleBDTGV2,RealEleTree.fWeight);
    EleIDBDTGV3_Real->Fill(fEleBDTGV3,RealEleTree.fWeight);
    EleIDBDTGV4_Real->Fill(fEleBDTGV4,RealEleTree.fWeight);
    EleIDBDTGV5_Real->Fill(fEleBDTGV5,RealEleTree.fWeight);
    EleIDBDTGV6_Real->Fill(fEleBDTGV6,RealEleTree.fWeight);
    EleIDBDTGV7_Real->Fill(fEleBDTGV7,RealEleTree.fWeight);
    EleIDBDTGV8_Real->Fill(fEleBDTGV8,RealEleTree.fWeight);
    EleIDBDTGV9_Real->Fill(fEleBDTGV9,RealEleTree.fWeight);      
    EleIDBDTGV10_Real->Fill(fEleBDTGV10,RealEleTree.fWeight);
    EleIDBDTGV11_Real->Fill(fEleBDTGV11,RealEleTree.fWeight);
    EleIDBDTGV12_Real->Fill(fEleBDTGV12,RealEleTree.fWeight);
    EleIDBDTGV13_Real->Fill(fEleBDTGV13,RealEleTree.fWeight);
    EleIDBDTGV14_Real->Fill(fEleBDTGV14,RealEleTree.fWeight);
    EleIDBDTGV15_Real->Fill(fEleBDTGV15,RealEleTree.fWeight);
    EleIDBDTGV16_Real->Fill(fEleBDTGV16,RealEleTree.fWeight);
    EleIDBDTGV17_Real->Fill(fEleBDTGV17,RealEleTree.fWeight);
    EleIDBDTGV18_Real->Fill(fEleBDTGV18,RealEleTree.fWeight);


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
    if (FakeEleTree.fEventNumber % 2 == 0) continue;
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

    FakeElectrons += FakeEleTree.fWeight;

    if (passICHEP2012( FakeEleTree.fElePt, FakeEleTree.fEleEta, fEleBDTGV0))
      FakeElectronPassICHEP2012 += FakeEleTree.fWeight;
    if (passIDMVA( FakeEleTree.fElePt, FakeEleTree.fEleEta, fEleBDTGV0) ) 
      FakeElectronPassIDMVA += FakeEleTree.fWeight;

    EleIDBDTGV0_Fake->Fill(fEleBDTGV0,FakeEleTree.fWeight);
    EleIDBDTGV1_Fake->Fill(fEleBDTGV1,FakeEleTree.fWeight);
    EleIDBDTGV2_Fake->Fill(fEleBDTGV2,FakeEleTree.fWeight);
    EleIDBDTGV3_Fake->Fill(fEleBDTGV3,FakeEleTree.fWeight);
    EleIDBDTGV4_Fake->Fill(fEleBDTGV4,FakeEleTree.fWeight);
    EleIDBDTGV5_Fake->Fill(fEleBDTGV5,FakeEleTree.fWeight);
    EleIDBDTGV6_Fake->Fill(fEleBDTGV6,FakeEleTree.fWeight);
    EleIDBDTGV7_Fake->Fill(fEleBDTGV7,FakeEleTree.fWeight);
    EleIDBDTGV8_Fake->Fill(fEleBDTGV8,FakeEleTree.fWeight);
    EleIDBDTGV9_Fake->Fill(fEleBDTGV9,FakeEleTree.fWeight);
    EleIDBDTGV10_Fake->Fill(fEleBDTGV10,FakeEleTree.fWeight);
    EleIDBDTGV11_Fake->Fill(fEleBDTGV11,FakeEleTree.fWeight);
    EleIDBDTGV12_Fake->Fill(fEleBDTGV12,FakeEleTree.fWeight);
    EleIDBDTGV13_Fake->Fill(fEleBDTGV13,FakeEleTree.fWeight);
    EleIDBDTGV14_Fake->Fill(fEleBDTGV14,FakeEleTree.fWeight);
    EleIDBDTGV15_Fake->Fill(fEleBDTGV15,FakeEleTree.fWeight);
    EleIDBDTGV16_Fake->Fill(fEleBDTGV16,FakeEleTree.fWeight);
    EleIDBDTGV17_Fake->Fill(fEleBDTGV17,FakeEleTree.fWeight);
    EleIDBDTGV18_Fake->Fill(fEleBDTGV18,FakeEleTree.fWeight);

  } //loop over electrons
  



  
  //*****************************************************************************************
  //Current Working Points
  //*****************************************************************************************
  cout << "ICHEP2012 Real Electron Efficiency : " << RealElectronPassICHEP2012 << " / " << RealElectrons << " = " << RealElectronPassICHEP2012/RealElectrons << endl;
  cout << "ICHEP2012 Fake Electron Efficiency : " << FakeElectronPassICHEP2012 << " / " << FakeElectrons << " = " << FakeElectronPassICHEP2012/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_ICHEP2012WP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassICHEP2012/RealElectrons , FakeElectronPassICHEP2012/FakeElectrons, "ROC_ICHEP2012WP"+label);

  cout << "IDMVA Real Electron Efficiency : " << RealElectronPassIDMVA << " / " << RealElectrons << " = " << RealElectronPassIDMVA/RealElectrons << endl;
  cout << "IDMVA Fake Electron Efficiency : " << FakeElectronPassIDMVA << " / " << FakeElectrons << " = " << FakeElectronPassIDMVA/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_IDMVAWP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassIDMVA/RealElectrons , FakeElectronPassIDMVA/FakeElectrons, "ROC_IDMVAWP"+label);



  cout << "BDTGV1 SameCutBasedSig Real Electron Efficiency : " << RealElectronPassBDTGV1_SameCutBasedSig << " / " << RealElectrons << " = " << RealElectronPassBDTGV1_SameCutBasedSig/RealElectrons << endl;
  cout << "BDTGV1 SameCutBasedSig Fake Electron Efficiency : " << FakeElectronPassBDTGV1_SameCutBasedSig << " / " << FakeElectrons << " = " << FakeElectronPassBDTGV1_SameCutBasedSig/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_BDTGV1_SameCutBasedSig_WP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassBDTGV1_SameCutBasedSig/RealElectrons , FakeElectronPassBDTGV1_SameCutBasedSig/FakeElectrons, "ROC_BDTGV1SameCutBasedSigWP"+label);

  cout << "BDTGV1 HalfCutBasedBkg Real Electron Efficiency : " << RealElectronPassBDTGV1_HalfCutBasedBkg << " / " << RealElectrons << " = " << RealElectronPassBDTGV1_HalfCutBasedBkg/RealElectrons << endl;
  cout << "BDTGV1 HalfCutBasedBkg Fake Electron Efficiency : " << FakeElectronPassBDTGV1_HalfCutBasedBkg << " / " << FakeElectrons << " = " << FakeElectronPassBDTGV1_HalfCutBasedBkg/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_BDTGV1_HalfCutBasedBkg_WP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassBDTGV1_HalfCutBasedBkg/RealElectrons , FakeElectronPassBDTGV1_HalfCutBasedBkg/FakeElectrons, "ROC_BDTGV1HalfCutBasedBkgWP"+label);

  cout << "BDTGV1 OneThirdCutBasedBkg Real Electron Efficiency : " << RealElectronPassBDTGV1_OneThirdCutBasedBkg << " / " << RealElectrons << " = " << RealElectronPassBDTGV1_OneThirdCutBasedBkg/RealElectrons << endl;
  cout << "BDTGV1 OneThirdCutBasedBkg Fake Electron Efficiency : " << FakeElectronPassBDTGV1_OneThirdCutBasedBkg << " / " << FakeElectrons << " = " << FakeElectronPassBDTGV1_OneThirdCutBasedBkg/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_BDTGV1_OneThirdCutBasedBkg_WP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassBDTGV1_OneThirdCutBasedBkg/RealElectrons , FakeElectronPassBDTGV1_OneThirdCutBasedBkg/FakeElectrons, "ROC_BDTGV1OneThirdCutBasedBkgWP"+label);

  cout << "**********************\n";
  cout << "Bkg At LHTight Signal Eff\n";

  Double_t BkgEffICHEP2012 = FakeElectronPassICHEP2012/FakeElectrons;
  Double_t SigEffICHEP2012 = RealElectronPassICHEP2012/RealElectrons;
  Double_t SigEffBDTGV0_SameBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV0_Real, EleIDBDTGV0_Fake, BkgEffICHEP2012);
  Double_t SigEffBDTGV1_SameBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV1_Real, EleIDBDTGV1_Fake, BkgEffICHEP2012);
  Double_t SigEffBDTGV2_SameBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV2_Real, EleIDBDTGV2_Fake, BkgEffICHEP2012);

  cout << "Signal Efficiency (wrt Cut-based) for : same bkg \n";
  cout << "BDTGV0 : " << SigEffBDTGV0_SameBkg/SigEffICHEP2012 <<  endl;
  cout << "BDTGV1 : " << SigEffBDTGV1_SameBkg/SigEffICHEP2012 <<  endl;
  cout << "BDTGV2 : " << SigEffBDTGV2_SameBkg/SigEffICHEP2012 <<  endl;


  cout << "**********************\n";

  Double_t BkgEffBDTGV0_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV0_Real, EleIDBDTGV0_Fake, SigEffICHEP2012);
  Double_t BkgEffBDTGV1_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV1_Real, EleIDBDTGV1_Fake, SigEffICHEP2012);
  Double_t BkgEffBDTGV2_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV2_Real, EleIDBDTGV2_Fake, SigEffICHEP2012);
  cout << "Bkg Efficiency (wrt Cut-based) for same sig eff \n";
  cout << "BDTGV0 : " << BkgEffBDTGV0_SameSig/BkgEffICHEP2012 << endl;
  cout << "BDTGV1 : " << BkgEffBDTGV1_SameSig/BkgEffICHEP2012 << endl;
  cout << "BDTGV2 : " << BkgEffBDTGV2_SameSig/BkgEffICHEP2012 << endl;

  cout << "**********************\n";



  //*****************************************************************************************
  //Make ROC curves
  //*****************************************************************************************
//   TGraphAsymmErrors* ROC_StandardLikelihood = MakeSigEffVsBkgEffGraph(EleIDStandardLikelihood_Real, EleIDStandardLikelihood_Fake, "ROC_StandardLikelihood"+label );
//   TGraphAsymmErrors* ROC_TMVALikelihood = MakeSigEffVsBkgEffGraph(EleIDTMVALikelihood_Real, EleIDTMVALikelihood_Fake, "ROC_TMVALikelihood"+label );
//   TGraphAsymmErrors* ROC_KNN = MakeSigEffVsBkgEffGraph(EleIDKNN_Real, EleIDKNN_Fake, "ROC_KNN"+label );
//   TGraphAsymmErrors* ROC_MLP = MakeSigEffVsBkgEffGraph(EleIDMLP_Real, EleIDMLP_Fake, "ROC_MLP"+label );
//   TGraphAsymmErrors* ROC_MLPBNN = MakeSigEffVsBkgEffGraph(EleIDMLPBNN_Real, EleIDMLPBNN_Fake, "ROC_MLPBNN"+label );
//   TGraphAsymmErrors* ROC_BDT = MakeSigEffVsBkgEffGraph(EleIDBDT_Real, EleIDBDT_Fake, "ROC_BDT"+label );
//   TGraphAsymmErrors* ROC_BDTG = MakeSigEffVsBkgEffGraph(EleIDBDTG_Real, EleIDBDTG_Fake, "ROC_BDTG"+label );
//   TGraphAsymmErrors* ROC_CombinedMVA = MakeSigEffVsBkgEffGraph(EleIDCombinedMVA_Real, EleIDCombinedMVA_Fake, "ROC_CombinedMVA"+label );

//   TGraphAsymmErrors* ROC_BDTG_C = MakeSigEffVsBkgEffGraph(EleIDBDTG_C_Real, EleIDBDTG_C_Fake, "ROC_BDTG_C"+label );
//   TGraphAsymmErrors* ROC_CombinedMVA_C = MakeSigEffVsBkgEffGraph(EleIDCombinedMVA_C_Real, EleIDCombinedMVA_C_Fake, "ROC_CombinedMVA_C"+label );


  TGraphAsymmErrors* ROC_BDTGV0 = MakeSigEffVsBkgEffGraph(EleIDBDTGV0_Real, EleIDBDTGV0_Fake, "ROC_BDTGV0"+label );
  TGraphAsymmErrors* ROC_BDTGV1 = MakeSigEffVsBkgEffGraph(EleIDBDTGV1_Real, EleIDBDTGV1_Fake, "ROC_BDTGV1"+label );
  TGraphAsymmErrors* ROC_BDTGV2 = MakeSigEffVsBkgEffGraph(EleIDBDTGV2_Real, EleIDBDTGV2_Fake, "ROC_BDTGV2"+label );
  TGraphAsymmErrors* ROC_BDTGV3 = MakeSigEffVsBkgEffGraph(EleIDBDTGV3_Real, EleIDBDTGV3_Fake, "ROC_BDTGV3"+label );
  TGraphAsymmErrors* ROC_BDTGV4 = MakeSigEffVsBkgEffGraph(EleIDBDTGV4_Real, EleIDBDTGV4_Fake, "ROC_BDTGV4"+label );
//   TGraphAsymmErrors* ROC_BDTGV5 = MakeSigEffVsBkgEffGraph(EleIDBDTGV5_Real, EleIDBDTGV5_Fake, "ROC_BDTGV5"+label );
//   TGraphAsymmErrors* ROC_BDTGV6 = MakeSigEffVsBkgEffGraph(EleIDBDTGV6_Real, EleIDBDTGV6_Fake, "ROC_BDTGV6"+label );
//   TGraphAsymmErrors* ROC_BDTGV7 = MakeSigEffVsBkgEffGraph(EleIDBDTGV7_Real, EleIDBDTGV7_Fake, "ROC_BDTGV7"+label );
//   TGraphAsymmErrors* ROC_BDTGV8 = MakeSigEffVsBkgEffGraph(EleIDBDTGV8_Real, EleIDBDTGV8_Fake, "ROC_BDTGV8"+label );
//   TGraphAsymmErrors* ROC_BDTGV9 = MakeSigEffVsBkgEffGraph(EleIDBDTGV9_Real, EleIDBDTGV9_Fake, "ROC_BDTGV9"+label );
//   TGraphAsymmErrors* ROC_BDTGV10 = MakeSigEffVsBkgEffGraph(EleIDBDTGV10_Real, EleIDBDTGV10_Fake, "ROC_BDTGV10"+label );
//   TGraphAsymmErrors* ROC_BDTGV11 = MakeSigEffVsBkgEffGraph(EleIDBDTGV11_Real, EleIDBDTGV11_Fake, "ROC_BDTGV11"+label );
//   TGraphAsymmErrors* ROC_BDTGV12 = MakeSigEffVsBkgEffGraph(EleIDBDTGV12_Real, EleIDBDTGV12_Fake, "ROC_BDTGV12"+label );
//   TGraphAsymmErrors* ROC_BDTGV13 = MakeSigEffVsBkgEffGraph(EleIDBDTGV13_Real, EleIDBDTGV13_Fake, "ROC_BDTGV13"+label );
//   TGraphAsymmErrors* ROC_BDTGV14 = MakeSigEffVsBkgEffGraph(EleIDBDTGV14_Real, EleIDBDTGV14_Fake, "ROC_BDTGV14"+label );
//   TGraphAsymmErrors* ROC_BDTGV15 = MakeSigEffVsBkgEffGraph(EleIDBDTGV15_Real, EleIDBDTGV15_Fake, "ROC_BDTGV15"+label );
//   TGraphAsymmErrors* ROC_BDTGV16 = MakeSigEffVsBkgEffGraph(EleIDBDTGV16_Real, EleIDBDTGV16_Fake, "ROC_BDTGV16"+label );
//   TGraphAsymmErrors* ROC_BDTGV17 = MakeSigEffVsBkgEffGraph(EleIDBDTGV17_Real, EleIDBDTGV17_Fake, "ROC_BDTGV17"+label );
//   TGraphAsymmErrors* ROC_BDTGV18 = MakeSigEffVsBkgEffGraph(EleIDBDTGV18_Real, EleIDBDTGV18_Fake, "ROC_BDTGV18"+label );

  //*****************************************************************************************
  //Find Cut with same signal efficiency Make ROC curves
  //*****************************************************************************************
  Double_t CutValue_BDTGV2_SameSig = FindCutValueAtFixedEfficiency(EleIDBDTGV2_Real, SigEffICHEP2012 );
  Double_t CutValue_BDTGV2_SameBkg = FindCutValueAtFixedEfficiency(EleIDBDTGV2_Fake, BkgEffICHEP2012 );
  cout << "BDTG V2 Cut Value @ Same Cut-Based Sig: " << CutValue_BDTGV2_SameSig << endl;
  cout << "BDTG V2 Cut Value @ Same Cut-Based Bkg: " << CutValue_BDTGV2_SameBkg << endl;

  //TFile *canvasFile = new TFile("ElectronIDMVAPerformancePlots.root","UPDATE");
  TLegend* legend;
  TCanvas* cv;
  string plotname;

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
//   colors.push_back(kBlue);    markers.push_back(21);
//   colors.push_back(kMagenta); markers.push_back(22);
//   colors.push_back(kCyan);    markers.push_back(34);
//   colors.push_back(kBlack);   markers.push_back(29);
//   colors.push_back(kGreen);   markers.push_back(33);
//   colors.push_back(kRed-2);   markers.push_back(33);
//   colors.push_back(kOrange);  markers.push_back(33);
//   colors.push_back(kBlue-2);  markers.push_back(33);
//   colors.push_back(kGreen-2);  markers.push_back(33);
//   colors.push_back(kMagenta-2);   markers.push_back(33);

  //*****************************************************************************************
  //*****************************************************************************************
  ROCGraphs.clear();
  GraphLabels.clear();
  plotname = "ElectronIDMVA"+label;


//   ROCGraphs.push_back(ROC_BDTGV8);
//   GraphLabels.push_back("2011MVA (V0)");
//   colors.push_back(kRed);
  
  ROCGraphs.push_back(ROC_BDTGV0);
  GraphLabels.push_back("ICHEP2012");
  colors.push_back(kBlue);
  
  ROCGraphs.push_back(ROC_BDTGV1);
  GraphLabels.push_back("V1");
  colors.push_back(kGreen+2);
  
//   ROCGraphs.push_back(ROC_BDTGV2);
//   GraphLabels.push_back("V2");
//   colors.push_back(kMagenta);
  
//   ROCGraphs.push_back(ROC_BDTGV3);
//   GraphLabels.push_back("Adaptive NTrees=4000");
//   colors.push_back(kCyan);
  
//   ROCGraphs.push_back(ROC_BDTGV4);
//   GraphLabels.push_back("Gradient");
//   colors.push_back(kOrange);

//   ROCGraphs.push_back(ROC_BDTGV5);
//   GraphLabels.push_back("V0+Single Crystal");
//   colors.push_back(kBlack);

//   ROCGraphs.push_back(ROC_BDTGV6);
//   GraphLabels.push_back("All ID Combined");
//   colors.push_back(kAzure+3);

//   ROCGraphs.push_back(ROC_BDTGV7);
//   GraphLabels.push_back("V7");
//   colors.push_back(kRed-2);

//   ROCGraphs.push_back(ROC_BDTGV9);
//   GraphLabels.push_back("All ID Combined (except H/E)");
//   colors.push_back(kRed-2);

//   ROCGraphs.push_back(ROC_BDTGV10);
//   GraphLabels.push_back("V0+PFIso 0.3");
//   colors.push_back(kMagenta+2);

//   ROCGraphs.push_back(ROC_BDTGV11);
//   GraphLabels.push_back("V0+PFIso 0.3 & 0.4");
//   colors.push_back(kCyan+2);

//   ROCGraphs.push_back(ROC_BDTGV12);
//   GraphLabels.push_back("V12");
//   colors.push_back(kGreen);
  
//   ROCGraphs.push_back(ROC_BDTGV13);
//   GraphLabels.push_back("V13");
//   colors.push_back(kMagenta-2);

//   ROCGraphs.push_back(ROC_BDTGV14);
//   GraphLabels.push_back("V0+AllID+AllIso");
//   colors.push_back(kGreen+2);

//   ROCGraphs.push_back(ROC_BDTGV15);
//   GraphLabels.push_back("V0+AllID+AllRelIso");
//   colors.push_back(kMagenta);

//   ROCGraphs.push_back(ROC_BDTGV16);
//   GraphLabels.push_back("V0+AllID+DetIso03");
//   colors.push_back(kGreen);

//   ROCGraphs.push_back(ROC_BDTGV17);
//   GraphLabels.push_back("V0+AllID+DetIso04");
//   colors.push_back(kCyan+2);

//   ROCGraphs.push_back(ROC_BDTGV18);
//   GraphLabels.push_back("V0+AllID+AllDetIso");
//   colors.push_back(kAzure+3);

  //*****************************************************************************************
  Double_t xmin = 0.0;
  Double_t xmax = 1.0;
  Double_t ymin = 0.0;
  Double_t ymax = 1.0;
// //   if (Option == 0 )                              { xmin = 0.15; xmax = 0.45; ymin = 0.75; ymax = 0.85; }

//   //For MC
//    if (Option == 0 )                                { xmin = 0.055; xmax = 0.090; ymin = 0.54; ymax = 0.66; }
// //   if (Option == 1 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.00; ymax = 1.00; }
//    if (Option == 2 )                              { xmin = 0.02; xmax = 0.04; ymin = 0.30; ymax = 0.45; }
//    if (Option == 3 )                                { xmin = 0.11; xmax = 0.145; ymin = 0.84; ymax = 0.93; }
// //   if (Option == 4 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.75; ymax = 1.00; }
//    if (Option == 5 )                              { xmin = 0.07; xmax = 0.12; ymin = 0.75; ymax = 0.85; }

//FOr Data
//   if (Option == 0 )                              { xmin = 0.10; xmax = 0.50; ymin = 0.70; ymax = 0.90; }
//   if (Option == 1 )                              { xmin = 0.05; xmax = 0.35; ymin = 0.30; ymax = 0.90; }
//   if (Option == 2 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.30; ymax = 0.80; }
//   if (Option == 3 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.80; ymax = 1.00; }
//   if (Option == 4 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.70; ymax = 1.00; }
//   if (Option == 5 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.70; ymax = 1.00; }

//   //**************************************
//   //Data: V0 vs V0+New ID Variables
//   //**************************************
//    if (Option == 0 )                                { xmin = 0.055; xmax = 0.090; ymin = 0.54; ymax = 0.66; }
// //   if (Option == 1 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.00; ymax = 1.00; }
//    if (Option == 2 )                              { xmin = 0.02; xmax = 0.04; ymin = 0.30; ymax = 0.45; }
//    if (Option == 3 )                                { xmin = 0.11; xmax = 0.145; ymin = 0.84; ymax = 0.93; }
// //   if (Option == 4 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.75; ymax = 1.00; }
//    if (Option == 5 )                              { xmin = 0.07; xmax = 0.12; ymin = 0.75; ymax = 0.85; }

//   //**************************************
//   //Data: V0 vs V0+Iso
//   //**************************************
//   if (Option == 0 )                                { xmin = 0.055; xmax = 0.090; ymin = 0.54; ymax = 0.70; }
// //   if (Option == 1 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.00; ymax = 1.00; }
//   if (Option == 2 )                              { xmin = 0.02; xmax = 0.04; ymin = 0.30; ymax = 0.45; }
//   if (Option == 3 )                                { xmin = 0.11; xmax = 0.145; ymin = 0.84; ymax = 0.95; }
// //   if (Option == 4 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.75; ymax = 1.00; }
//   if (Option == 5 )                              { xmin = 0.07; xmax = 0.12; ymin = 0.77; ymax = 0.88; }
  
  //**************************************
  //Data: V0 vs Best MVAs
  //**************************************
//   if (Option == 0 )                                { xmin = 0.040; xmax = 0.090; ymin = 0.40; ymax = 0.70; }
//   if (Option == 1 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.00; ymax = 1.00; }
//   if (Option == 2 )                              { xmin = 0.015; xmax = 0.035; ymin = 0.30; ymax = 0.45; }
//   if (Option == 3 )                                { xmin = 0.09; xmax = 0.15; ymin = 0.75; ymax = 0.95; }
//   if (Option == 4 )                              { xmin = 0.00; xmax = 0.2; ymin = 0.75; ymax = 1.00; }
//   if (Option == 5 )                              { xmin = 0.05; xmax = 0.12; ymin = 0.65; ymax = 0.90; }


//HWW115
//   if (Option == 3 )                                { xmin = 0.07; xmax = 0.16; ymin = 0.75; ymax = 0.95; }
//Zee
//   if (Option == 3 )                                { xmin = 0.07; xmax = 0.16; ymin = 0.75; ymax = 0.95; }
//WJets vs HWW115
//    if (Option == 3 )                                { xmin = 0.009; xmax = 0.045; ymin = 0.75; ymax = 0.96; }



  cv = new TCanvas("cv", "cv", 800, 600);

//    legend = new TLegend(0.45,0.20,0.75,0.50);
  legend = new TLegend(0.54,0.14,0.94,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(ROCGraphs[i],GraphLabels[i].c_str(), "LP");

    ROCGraphs[i]->SetMarkerColor(colors[i]);
    ROCGraphs[i]->SetLineColor(colors[i]);
    ROCGraphs[i]->SetMarkerSize(0.5);
   
    ROCGraphs[i]->GetXaxis()->SetRangeUser(xmin,xmax);    
    ROCGraphs[i]->GetYaxis()->SetRangeUser(ymin,ymax);    
    if (i==0) {
      ROCGraphs[i]->Draw("AP");
    } else {
      ROCGraphs[i]->Draw("Psame");
    }
  }

//   legend->AddEntry(ROC_LHTightWP, "NEW LH WP", "P");
//   ROC_LHTightWP->SetFillColor(kBlue);
//   ROC_LHTightWP->SetMarkerColor(kBlue);
//   ROC_LHTightWP->SetMarkerStyle(34);
//   ROC_LHTightWP->SetMarkerSize(2.5);
//   ROC_LHTightWP->Draw("Psame");
  
//   legend->AddEntry(ROC_BDTGV1WP, "V1 WP", "P");
//   ROC_BDTGV1WP->SetFillColor(kBlack);
//   ROC_BDTGV1WP->SetMarkerColor(kBlack);
//   ROC_BDTGV1WP->SetMarkerStyle(34);
//   ROC_BDTGV1WP->SetMarkerSize(2.5);
//   ROC_BDTGV1WP->Draw("Psame");

//   legend->AddEntry(ROC_BDTGV2WP, "V2 WP", "P");
//   ROC_BDTGV2WP->SetFillColor(kRed);
//   ROC_BDTGV2WP->SetMarkerColor(kRed);
//   ROC_BDTGV2WP->SetMarkerStyle(34);
//   ROC_BDTGV2WP->SetMarkerSize(2.5);
//   ROC_BDTGV2WP->Draw("Psame");
 

  legend->AddEntry(ROC_ICHEP2012WP, "ICHEP2012 WP", "P");
  ROC_ICHEP2012WP->SetFillColor(kGreen+3);
  ROC_ICHEP2012WP->SetMarkerColor(kGreen+3);
  ROC_ICHEP2012WP->SetMarkerStyle(34);
  ROC_ICHEP2012WP->SetMarkerSize(2.5);
  ROC_ICHEP2012WP->Draw("Psame");

  legend->AddEntry(ROC_IDMVAWP, "IDMVA", "P");
  ROC_IDMVAWP->SetFillColor(kBlack);
  ROC_IDMVAWP->SetMarkerColor(kBlack);
  ROC_IDMVAWP->SetMarkerStyle(34);
  ROC_IDMVAWP->SetMarkerSize(2.5);
  ROC_IDMVAWP->Draw("Psame");

//   legend->AddEntry(ROC_BDTGV15_SameCutBasedSig_WP, "BDT @ Same 2011 MVA Sig", "P");
//   ROC_BDTGV15_SameCutBasedSig_WP->SetFillColor(kMagenta);
//   ROC_BDTGV15_SameCutBasedSig_WP->SetMarkerColor(kMagenta);
//   ROC_BDTGV15_SameCutBasedSig_WP->SetMarkerStyle(34);
//   ROC_BDTGV15_SameCutBasedSig_WP->SetMarkerSize(2.5);
//   ROC_BDTGV15_SameCutBasedSig_WP->Draw("Psame");

//   legend->AddEntry(ROC_BDTGV15_HalfCutBasedBkg_WP, "BDT @ 50% 2011 MVA Bkg", "P");
//   ROC_BDTGV15_HalfCutBasedBkg_WP->SetFillColor(kRed);
//   ROC_BDTGV15_HalfCutBasedBkg_WP->SetMarkerColor(kRed);
//   ROC_BDTGV15_HalfCutBasedBkg_WP->SetMarkerStyle(34);
//   ROC_BDTGV15_HalfCutBasedBkg_WP->SetMarkerSize(2.5);
//   ROC_BDTGV15_HalfCutBasedBkg_WP->Draw("Psame");

  legend->Draw();
  
  cv->SaveAs(("ROCGraphs_" + plotname + ".gif").c_str());
//   canvasFile->WriteTObject(cv,("ROCGraphs_" + plotname).c_str(), "WriteDelete") ;


  gBenchmark->Show("WWTemplate");       
} 

