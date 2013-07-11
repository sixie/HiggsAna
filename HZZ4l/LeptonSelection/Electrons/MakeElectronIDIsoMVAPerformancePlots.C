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
#include <MitStyle.h>
#include "TLegend.h"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"

#include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"
#include "HiggsAna/Ntupler/interface/HiggsAnaDefs.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"

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
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
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


//*************************************************************************************************
//
//*************************************************************************************************
void OptimizeIDAndIsoCuts(TH2F* signalHist, TH2F* bkgHist, 
                          Double_t &CutValueX, Double_t &CutValueY,
                          Double_t &SigEffAtTarget,
                          Double_t targetBkgEff ) {

  Double_t targetSignalEff = 0;
  Double_t targetXCutValue = -9999;
  Double_t targetYCutValue = -9999;

  const UInt_t nXPoints = signalHist->GetXaxis()->GetNbins();
  const UInt_t nYPoints = signalHist->GetYaxis()->GetNbins();
 
  //first integrate to get total
  double NSigTotal = 0;
  double NBkgTotal = 0;
  for (UInt_t q=0; q < nXPoints+2; ++q) {
    for (UInt_t r=0; r < nYPoints+2; ++r) {
      NSigTotal += signalHist->GetBinContent(q,r);
      NBkgTotal += bkgHist->GetBinContent(q,r);
    }
  }

  cout << "Total Sig: " << NSigTotal << endl;
  cout << "Total Bkg: " << NBkgTotal << endl;
  

  for(UInt_t a=0; a < nXPoints; ++a) {

    Double_t bestCurrentBkgEff = 0;
    Double_t bestCurrentSigEff = 0;
    Double_t bestCurrentXCutValue = 0;
    Double_t bestCurrentYCutValue = 0;
    for(UInt_t b=0; b < nYPoints; ++b) {
      //cout << "start " << a << " " << b << endl;
      Double_t nsig = 0;
      Double_t nbkg = 0;
      for (UInt_t q=a; q < nXPoints+2; ++q) {
        for (UInt_t r=b; r < nYPoints+2; ++r) {
          nsig += signalHist->GetBinContent(q,r);
          nbkg += bkgHist->GetBinContent(q,r);
        }
      }
      Double_t sigEff = nsig / NSigTotal;
      Double_t bkgEff = nbkg / NBkgTotal;
//     cout << targetSignalEff << " : " << ratio << " , " << bkgEff << " : " << signalHist->GetXaxis()->GetBinCenter(b) << endl;      
      if (fabs(targetBkgEff - bkgEff) < fabs(targetBkgEff - bestCurrentBkgEff)) {
        bestCurrentBkgEff = bkgEff;
        bestCurrentSigEff = sigEff;
        bestCurrentXCutValue = signalHist->GetXaxis()->GetBinCenter(a);
        bestCurrentYCutValue = signalHist->GetYaxis()->GetBinCenter(b);       
      }
      //cout << "done " << a << " " << b << endl;
    }
    cout << "For x bin " << a << ", Best Sig Eff: " << bestCurrentSigEff << " , bkg eff = " << bestCurrentBkgEff << endl;
    //check if it's the best one so far
    if (bestCurrentSigEff > targetSignalEff && fabs(targetBkgEff - bestCurrentBkgEff)/bestCurrentBkgEff < 0.05) {
      targetSignalEff = bestCurrentSigEff;
      targetXCutValue = bestCurrentXCutValue;
      targetYCutValue = bestCurrentYCutValue;
    }
  }

  cout << "Final Best Sig Eff: " << targetSignalEff << endl;
  cout << "Final Cuts: " << targetXCutValue << " " << targetYCutValue << endl;
  CutValueX = targetXCutValue;
  CutValueY = targetYCutValue;
  SigEffAtTarget = targetSignalEff;

}





Bool_t passIDIsoMVA( Double_t fElePt, Double_t fEleSCEta, Double_t MVAValue) {

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
  if (MVABin == 0) MVACut = 0.5094 ;
  if (MVABin == 1) MVACut = -0.0394 ;
  if (MVABin == 2) MVACut = 0.1902 ; 
  if (MVABin == 3) MVACut = 0.839 ;
  if (MVABin == 4) MVACut = 0.673 ;
  if (MVABin == 5) MVACut = 0.8078 ;  


  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}



Bool_t passIDMVA( Double_t fElePt, Double_t fEleSCEta, Double_t MVAValue, Int_t Option = 0) {

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

  Double_t cutAdjustment = Option * (0.1);

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = 0.5094 + cutAdjustment;
  if (MVABin == 1) MVACut = -0.0394 + cutAdjustment;
  if (MVABin == 2) MVACut = 0.1902 + cutAdjustment; 
  if (MVABin == 3) MVACut = 0.839 + cutAdjustment;
  if (MVABin == 4) MVACut = 0.673 + cutAdjustment;
  if (MVABin == 5) MVACut = 0.8078 + cutAdjustment;  


  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}



Bool_t passIsoMVA( Double_t fMuPt, Double_t fMuEta, Double_t MVAValue, Int_t Option = 0) {

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

  Double_t cutAdjustment = Option * (0.1);

  Double_t MVACut = -999;
  //Cut Values for MVA-V10 (Detector Based Iso)
  if (MVABin == 0) MVACut = -0.616 + cutAdjustment;
  if (MVABin == 1) MVACut = -0.740 + cutAdjustment;
  if (MVABin == 2) MVACut = -0.045 + cutAdjustment;
  if (MVABin == 3) MVACut = 0.205 + cutAdjustment;

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}


Bool_t passIDAndIsoMVA_Loose( Double_t fPt, Double_t fEta, Double_t IDMVAValue, Double_t IsoMVAValue) {

  Int_t subdet = 0;
  if (fabs(fEta) < 0.8) subdet = 0;
  else if (fabs(fEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (fPt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t IDMVACut = -999;
  Double_t IsoMVACut = -999;
  //Cut Values for MVA-V10 (Detector Based Iso)
  if (MVABin == 0) {
    IDMVACut =  0.369;
    IsoMVACut = 0.385;
  }
  if (MVABin == 1) {
    IDMVACut = -0.025;
    IsoMVACut = -0.083;
  }
  if (MVABin == 2) {
    IDMVACut = 0.531;
    IsoMVACut = -0.573;
  }
  if (MVABin == 3) {
    IDMVACut = 0.735;
    IsoMVACut = 0.413;
  }
  if (MVABin == 4) {
    IDMVACut = 0.467;
    IsoMVACut = 0.271;
  }
  if (MVABin == 5) {
    IDMVACut = 0.795;
    IsoMVACut = 0.135;
  }

  if (IDMVAValue > IDMVACut && IsoMVAValue > IsoMVACut) return kTRUE;
  return kFALSE;
}

Bool_t passIDAndIsoMVA_Tight( Double_t fPt, Double_t fEta, Double_t IDMVAValue, Double_t IsoMVAValue) {

  Int_t subdet = 0;
  if (fabs(fEta) < 0.8) subdet = 0;
  else if (fabs(fEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (fPt > 10.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t IDMVACut = -999;
  Double_t IsoMVACut = -999;
  //Cut Values for MVA-V10 (Detector Based Iso)
  if (MVABin == 0) {
    IDMVACut =  0.093;
    IsoMVACut = 0.553;
  }
  if (MVABin == 1) {
    IDMVACut = 0.451;
    IsoMVACut = -0.237;
  }
  if (MVABin == 2) {
    IDMVACut = 0.595;
    IsoMVACut = -0.573;
  }
  if (MVABin == 3) {
    IDMVACut = 0.881;
    IsoMVACut = 0.521;
  }
  if (MVABin == 4) {
    IDMVACut = 0.731;
    IsoMVACut = 0.531;
  }
  if (MVABin == 5) {
    IDMVACut = 0.819;
    IsoMVACut = 0.493;
  }

  if (IDMVAValue > IDMVACut && IsoMVAValue > IsoMVACut) return kTRUE;
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
void MakeElectronIDIsoMVAPerformancePlots(string RealElectronFile, string FakeElectronFile, string Label, Int_t Option)
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

  TH2F *Ele_IDMVAV0_Vs_IsoMVA_Real = new TH2F(("Ele_IDMVAV0_Vs_IsoMVA_Real"+label).c_str(), "; ID BDT V0; Iso BDT ; Number of Events ",  100, -1 , 1, 100, -1, 1);
  TH2F *Ele_IDMVAV0_Vs_IsoMVA_Fake = new TH2F(("Ele_IDMVAV0_Vs_IsoMVA_Fake"+label).c_str(), "; ID BDT V0; Iso BDT ; Number of Events ",  100, -1 , 1, 100, -1, 1);


   TH2F *Ele_IDMVAV1_Vs_IsoMVA_Real = new TH2F(("Ele_IDMVAV1_Vs_IsoMVA_Real"+label).c_str(), "; ID BDT V1; Iso BDT ; Number of Events ",  100, -1 , 1, 100, -1, 1);
   TH2F *Ele_IDMVAV1_Vs_IsoMVA_Fake = new TH2F(("Ele_IDMVAV1_Vs_IsoMVA_Fake"+label).c_str(), "; ID BDT V1; Iso BDT ; Number of Events ",  100, -1 , 1, 100, -1, 1);

//   TH2F *Ele_IDMVAV1_Vs_IsoMVA_Real = new TH2F(("Ele_IDMVAV1_Vs_IsoMVA_Real"+label).c_str(), "; ID BDT V1; Iso BDT ; Number of Events ",  100, 0.77 , 0.97, 100, 0.31, 0.51);
//   TH2F *Ele_IDMVAV1_Vs_IsoMVA_Fake = new TH2F(("Ele_IDMVAV1_Vs_IsoMVA_Fake"+label).c_str(), "; ID BDT V1; Iso BDT ; Number of Events ",  100, 0.77 , 0.97, 100, 0.31, 0.51);


  TH2F *Ele_IDMVAV2_Vs_IsoMVA_Real = new TH2F(("Ele_IDMVAV2_Vs_IsoMVA_Real"+label).c_str(), "; ID BDT V2; Iso BDT ; Number of Events ",  100, -1 , 1, 100, -1, 1);
  TH2F *Ele_IDMVAV2_Vs_IsoMVA_Fake = new TH2F(("Ele_IDMVAV2_Vs_IsoMVA_Fake"+label).c_str(), "; ID BDT V2; Iso BDT ; Number of Events ",  100, -1 , 1, 100, -1, 1);



  TH1F *EleIDIsoBDT_V0_Real = new TH1F(("EleIDIsoBDT_V0_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  40000, -2 , 2);
  TH1F *EleIDIsoBDT_V0_Fake = new TH1F(("EleIDIsoBDT_V0_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  40000, -2 , 2);
  TH1F *EleIDIsoBDT_V1_Real = new TH1F(("EleIDIsoBDT_V1_Real"+label).c_str(), "; BDTG V1 ; Number of Events ",  40000, -2 , 2);
  TH1F *EleIDIsoBDT_V1_Fake = new TH1F(("EleIDIsoBDT_V1_Fake"+label).c_str(), "; BDTG V1 ; Number of Events ",  40000, -2 , 2);

  TH1F *EleDetIsoWithCiCTight_Real = new TH1F(("EleDetIsoWithCiCTight_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -1 , 3);
  TH1F *EleDetIsoWithCiCTight_Fake = new TH1F(("EleDetIsoWithCiCTight_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -1 , 3);

  TH1F *EleIDBDTWithIsoBDTCut0_Real = new TH1F(("EleIDBDTWithIsoBDTCut0_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTWithIsoBDTCut1_Real = new TH1F(("EleIDBDTWithIsoBDTCut1_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTWithIsoBDTCut2_Real = new TH1F(("EleIDBDTWithIsoBDTCut2_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTWithIsoBDTCut3_Real = new TH1F(("EleIDBDTWithIsoBDTCut3_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTWithIsoBDTCut4_Real = new TH1F(("EleIDBDTWithIsoBDTCut4_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTWithIsoBDTCut5_Real = new TH1F(("EleIDBDTWithIsoBDTCut5_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTWithIsoBDTCut0_Fake = new TH1F(("EleIDBDTWithIsoBDTCut0_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTWithIsoBDTCut1_Fake = new TH1F(("EleIDBDTWithIsoBDTCut1_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTWithIsoBDTCut2_Fake = new TH1F(("EleIDBDTWithIsoBDTCut2_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTWithIsoBDTCut3_Fake = new TH1F(("EleIDBDTWithIsoBDTCut3_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTWithIsoBDTCut4_Fake = new TH1F(("EleIDBDTWithIsoBDTCut4_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTWithIsoBDTCut5_Fake = new TH1F(("EleIDBDTWithIsoBDTCut5_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);

  TH1F *EleIsoBDTWithIDBDTCut0_Real = new TH1F(("EleIsoBDTWithIDBDTCut0_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDTWithIDBDTCut1_Real = new TH1F(("EleIsoBDTWithIDBDTCut1_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDTWithIDBDTCut2_Real = new TH1F(("EleIsoBDTWithIDBDTCut2_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDTWithIDBDTCut3_Real = new TH1F(("EleIsoBDTWithIDBDTCut3_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDTWithIDBDTCut4_Real = new TH1F(("EleIsoBDTWithIDBDTCut4_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDTWithIDBDTCut5_Real = new TH1F(("EleIsoBDTWithIDBDTCut5_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDTWithIDBDTCut0_Fake = new TH1F(("EleIsoBDTWithIDBDTCut0_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDTWithIDBDTCut1_Fake = new TH1F(("EleIsoBDTWithIDBDTCut1_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDTWithIDBDTCut2_Fake = new TH1F(("EleIsoBDTWithIDBDTCut2_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDTWithIDBDTCut3_Fake = new TH1F(("EleIsoBDTWithIDBDTCut3_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDTWithIDBDTCut4_Fake = new TH1F(("EleIsoBDTWithIDBDTCut4_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIsoBDTWithIDBDTCut5_Fake = new TH1F(("EleIsoBDTWithIDBDTCut5_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);

  Double_t RealElectrons = 0;
  Double_t FakeElectrons = 0;
  Double_t RealElectronPassCiCTightDetIso = 0;
  Double_t FakeElectronPassCiCTightDetIso = 0;
  Double_t RealElectronPassCiCTightDetIso015 = 0;
  Double_t FakeElectronPassCiCTightDetIso015 = 0;
  Double_t RealElectronPassIDMVAAndIsoMVALoose = 0;
  Double_t FakeElectronPassIDMVAAndIsoMVALoose = 0;
  Double_t RealElectronPassIDMVAAndIsoMVATight = 0;
  Double_t FakeElectronPassIDMVAAndIsoMVATight = 0;
  Double_t RealElectronPassIDIsoMVA = 0;
  Double_t FakeElectronPassIDIsoMVA = 0;


  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

  //Variables
  Float_t                 fWeight = 1.0;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Float_t                 fElePt; 
  Float_t                 fEleEta; 
  Float_t                 fElePhi; 
  Float_t                 fEleSCEt; 
  Float_t                 fEleSCEta; 
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
  Float_t                 fEleOneOverEMinusOneOverP; 
  Float_t                 fEleESeedClusterOverPIn; 
  Float_t                 fEleIP3d; 
  Float_t                 fEleIP3dSig; 

  //Isolation Variables
  //Float_t                 fElePFMVA;
//   Float_t                 fEleChargedIso03; 
//   Float_t                 fEleNeutralHadronIso03; 
//   Float_t                 fEleGammaIso03; 
//   Float_t                 fEleChargedIso04; 
//   Float_t                 fEleNeutralHadronIso04; 
//   Float_t                 fEleGammaIso04; 
  Float_t                 fEleTrkIso03; 
  Float_t                 fEleEMIso03; 
  Float_t                 fEleHadIso03; 
//   Float_t                 fEleTrkIso04; 
//   Float_t                 fEleEMIso04; 
//   Float_t                 fEleHadIso04; 
  Float_t                 fRho; 
  UInt_t                  fNVertices; 
//   UInt_t                  fEleTriggerBit;
//   Bool_t                  fElePassDenominator;
//   Bool_t                  fElePassDenominatorSmurf;

  Float_t                 fEleIDBDT_V0;
  Float_t                 fEleIDBDT_V1;
  Float_t                 fEleIDBDT_V2;
  Float_t                 fEleIsoBDT;
  Float_t                 fEleIDIsoBDT_V0;
  Float_t                 fEleIDIsoBDT_V1;


  Bool_t fEleDenomFake;


  //*****************************************************************************************
  //RealEleTree
  //*****************************************************************************************
  TTree *RealEleTree = getTreeFromFile(RealElectronFile.c_str(), "Electrons");
  assert(RealEleTree);
  RealEleTree->SetBranchAddress( "weight", &fWeight);
  RealEleTree->SetBranchAddress( "run", &fRunNumber);
  RealEleTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  RealEleTree->SetBranchAddress( "event", &fEventNumber);
  RealEleTree->SetBranchAddress( "pt", &fElePt); 
  RealEleTree->SetBranchAddress( "eta", &fEleEta); 
  RealEleTree->SetBranchAddress( "phi", &fElePhi); 
  RealEleTree->SetBranchAddress( "scEt", &fEleSCEt); 
  RealEleTree->SetBranchAddress( "scEta", &fEleSCEta); 
  RealEleTree->SetBranchAddress( "combPFIsoHWW", &fElePFIso);
  RealEleTree->SetBranchAddress("matchConv",&fEleMatchedConversion); 
  RealEleTree->SetBranchAddress("dist",&fEleConvDist); 
  RealEleTree->SetBranchAddress("dcot",&fEleConvDCot); 
  RealEleTree->SetBranchAddress("missHits",&fEleNMissHits); 
  RealEleTree->SetBranchAddress( "see", &fEleSigmaIEtaIEta); 
  RealEleTree->SetBranchAddress( "deta", &fEleDEtaIn); 
  RealEleTree->SetBranchAddress( "dphi", &fEleDPhiIn); 
  RealEleTree->SetBranchAddress( "HoE", &fEleHoverE); 
  RealEleTree->SetBranchAddress( "d0", &fEleD0); 
  RealEleTree->SetBranchAddress( "dz", &fEleDZ); 
  RealEleTree->SetBranchAddress( "fbrem", &fEleFBrem); 
  RealEleTree->SetBranchAddress( "EoP", &fEleEOverP); 
  RealEleTree->SetBranchAddress( "EoPout", &fEleESeedClusterOverPout); 
  RealEleTree->SetBranchAddress( "spp", &fEleSigmaIPhiIPhi); 
  RealEleTree->SetBranchAddress( "IoEmIoP", &fEleOneOverEMinusOneOverP); 
  RealEleTree->SetBranchAddress( "EoPin", &fEleESeedClusterOverPIn); 
  RealEleTree->SetBranchAddress( "ip3d", &fEleIP3d); 
  RealEleTree->SetBranchAddress( "ip3ds", &fEleIP3dSig); 
  RealEleTree->SetBranchAddress( "rho", &fRho); 
  RealEleTree->SetBranchAddress( "vertices", &fNVertices); 

  RealEleTree->SetBranchAddress("trkIso03",&fEleTrkIso03);
  RealEleTree->SetBranchAddress("ecalIso03",&fEleEMIso03);
  RealEleTree->SetBranchAddress("hcalIso03",&fEleHadIso03);

  RealEleTree->SetBranchAddress( "EleIDMVA_BDTG_V0", &fEleIDBDT_V0 );
  RealEleTree->SetBranchAddress( "EleIDMVA_BDTG_V1", &fEleIDBDT_V1 );
  RealEleTree->SetBranchAddress( "EleIDMVA_BDTG_V2", &fEleIDBDT_V2 );
  RealEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V0", &fEleIsoBDT); 
  RealEleTree->SetBranchAddress( "EleIDIsoMVA_BDTG_V0", &fEleIDIsoBDT_V0); 
  RealEleTree->SetBranchAddress( "EleIDIsoMVA_BDTG_V1", &fEleIDIsoBDT_V1); 

  RealEleTree->SetBranchAddress( "DenomFake", &fEleDenomFake); 


  for(UInt_t ientry=0; ientry < RealEleTree->GetEntries(); ientry++) {       	
    RealEleTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
        
    Double_t rho = 0;
    if (!(TMath::IsNaN(fRho) || isinf(fRho))) rho = fRho;

   //don't evaluate performance using training events
    if (fElePt < 5) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fEleEta) < 0.8) subdet = 0;
    else if (fabs(fEleEta) < 1.479) subdet = 1;
    else subdet = 2;

    Int_t ptBin = 0;
    if (fElePt > 10.0) ptBin = 1;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 5) passCuts = (subdet == 2 && ptBin == 1); 

    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;

    if( fEleIP3dSig > 4 ) continue;
    if( fEleTrkIso03/fElePt > 0.7 ) continue;
    if( fEleNMissHits >  1 ) continue;


    RealElectrons += fWeight;

    Double_t DetIso03PUCorrection = 0;
    if (fabs(fEleEta) < 1.5) DetIso03PUCorrection = rho * (0.078 + 0.026);
    else DetIso03PUCorrection = rho * (0.046 + 0.072);

    

    //baseline working point
    if (PassCiCID(fElePt, fEleSCEt, fEleEta, 1, Bool_t(fabs(fEleEta) < 1.479), 
                  fEleFBrem, fEleEOverP, fEleHoverE, fEleSigmaIEtaIEta, fEleDPhiIn, fEleDEtaIn, 
                  fEleESeedClusterOverPIn, 0, 
                  0, 0, 0, 
                  0, 0, 
                  fEleD0, fEleIP3dSig, 
                  1, kFALSE)
        &&
        (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.25
      ) {
      RealElectronPassCiCTightDetIso += fWeight;
    }

    if (PassCiCID(fElePt, fEleSCEt, fEleEta, 1, Bool_t(fabs(fEleEta) < 1.479), 
                  fEleFBrem, fEleEOverP, fEleHoverE, fEleSigmaIEtaIEta, fEleDPhiIn, fEleDEtaIn, 
                  fEleESeedClusterOverPIn, 0, 
                  0, 0, 0, 
                  0, 0, 
                  fEleD0, fEleIP3dSig, 
                  1, kFALSE)
        &&
        (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.15
      ) {
      RealElectronPassCiCTightDetIso015 += fWeight;
    }



    //mva working point
    if (passIDIsoMVA( fElePt, fEleEta, fEleIDIsoBDT_V1)) {
      RealElectronPassIDIsoMVA += fWeight;
    }

    //Separated Working point
    if (passIDAndIsoMVA_Loose(fElePt, fEleEta, fEleIDBDT_V1,fEleIsoBDT)) {
      RealElectronPassIDMVAAndIsoMVALoose += fWeight;
    }
    if (passIDAndIsoMVA_Tight(fElePt, fEleEta, fEleIDBDT_V1,fEleIsoBDT)) {
      RealElectronPassIDMVAAndIsoMVATight += fWeight;
    }

    if (PassCiCID(fElePt, fEleSCEt, fEleEta, 1, Bool_t(fabs(fEleEta) < 1.479), 
                  fEleFBrem, fEleEOverP, fEleHoverE, fEleSigmaIEtaIEta, fEleDPhiIn, fEleDEtaIn, 
                  fEleESeedClusterOverPIn, 0, 
                  0, 0, 0, 
                  0, 0, 
                  fEleD0, fEleIP3dSig, 
                  1)) {
 //      cout << "Real DetIso: " << (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt << endl;
      EleDetIsoWithCiCTight_Real->Fill(min((fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt, 2.5), fWeight);
    } else {
      EleDetIsoWithCiCTight_Real->Fill(2.5,fWeight);
    }


    //Scan IDMVA for fixed IsoMVA cut
    if (passIDMVA(fElePt, fEleEta, fEleIDBDT_V2, -3)) {
      EleIsoBDTWithIDBDTCut0_Real->Fill(fEleIsoBDT,fWeight);
    } else {
      EleIsoBDTWithIDBDTCut0_Real->Fill(-0.99,fWeight);
    }
    if (passIDMVA(fElePt, fEleEta, fEleIDBDT_V2, -2)) {
      EleIsoBDTWithIDBDTCut1_Real->Fill(fEleIsoBDT,fWeight);
    } else {
      EleIsoBDTWithIDBDTCut1_Real->Fill(-0.99,fWeight);
    }
    if (passIDMVA(fElePt, fEleEta, fEleIDBDT_V2, -1)) {
      EleIsoBDTWithIDBDTCut2_Real->Fill(fEleIsoBDT,fWeight);
    } else {
      EleIsoBDTWithIDBDTCut2_Real->Fill(-0.99,fWeight);
    }
    if (passIDMVA(fElePt, fEleEta, fEleIDBDT_V2, 0)) {
      EleIsoBDTWithIDBDTCut3_Real->Fill(fEleIsoBDT,fWeight);
    } else {
      EleIsoBDTWithIDBDTCut3_Real->Fill(-0.99,fWeight);
    }
    if (passIDMVA(fElePt, fEleEta, fEleIDBDT_V2, 1)) {
      EleIsoBDTWithIDBDTCut4_Real->Fill(fEleIsoBDT,fWeight);
    } else {
      EleIsoBDTWithIDBDTCut4_Real->Fill(-0.99,fWeight);
    }
    if (passIDMVA(fElePt, fEleEta, fEleIDBDT_V2, 2)) {
      EleIsoBDTWithIDBDTCut5_Real->Fill(fEleIsoBDT,fWeight);
    } else {
      EleIsoBDTWithIDBDTCut5_Real->Fill(-0.99,fWeight);
    }
    //Scan IsoMVA for fixed IDMVA cut    

    if (passIsoMVA(fElePt, fEleEta, fEleIsoBDT, -3)) {
      EleIDBDTWithIsoBDTCut0_Real->Fill(fEleIDBDT_V2,fWeight);
    } else {
      EleIDBDTWithIsoBDTCut0_Real->Fill(-0.99,fWeight);
    }
    if (passIsoMVA(fElePt, fEleEta, fEleIsoBDT, -2)) {
      EleIDBDTWithIsoBDTCut1_Real->Fill(fEleIDBDT_V2,fWeight);
    } else {
      EleIDBDTWithIsoBDTCut1_Real->Fill(-0.99,fWeight);
    }
    if (passIsoMVA(fElePt, fEleEta, fEleIsoBDT, -1)) {
      EleIDBDTWithIsoBDTCut2_Real->Fill(fEleIDBDT_V2,fWeight);
    } else {
      EleIDBDTWithIsoBDTCut2_Real->Fill(-0.99,fWeight);
    }
    if (passIsoMVA(fElePt, fEleEta, fEleIsoBDT, 0)) {
      EleIDBDTWithIsoBDTCut3_Real->Fill(fEleIDBDT_V2,fWeight);
    } else {
      EleIDBDTWithIsoBDTCut3_Real->Fill(-0.99,fWeight);
    }
    if (passIsoMVA(fElePt, fEleEta, fEleIsoBDT, 1)) {
      EleIDBDTWithIsoBDTCut4_Real->Fill(fEleIDBDT_V2,fWeight);
    } else {
      EleIDBDTWithIsoBDTCut4_Real->Fill(-0.99,fWeight);
    }
    if (passIsoMVA(fElePt, fEleEta, fEleIsoBDT, 2)) {
      EleIDBDTWithIsoBDTCut5_Real->Fill(fEleIDBDT_V2,fWeight);
    } else {
      EleIDBDTWithIsoBDTCut5_Real->Fill(-0.99,fWeight);
    }

    Ele_IDMVAV0_Vs_IsoMVA_Real->Fill(fEleIDBDT_V0,fEleIsoBDT,fWeight);
    Ele_IDMVAV1_Vs_IsoMVA_Real->Fill(fEleIDBDT_V1,fEleIsoBDT,fWeight);
    Ele_IDMVAV2_Vs_IsoMVA_Real->Fill(fEleIDBDT_V2,fEleIsoBDT,fWeight);

    //IDIsoCombined
    EleIDIsoBDT_V0_Real->Fill(fEleIDIsoBDT_V0,fWeight);
    EleIDIsoBDT_V1_Real->Fill(fEleIDIsoBDT_V1,fWeight);
//     cout << "Real " << fEleIDIsoBDT << endl;

  } 
  





  //*****************************************************************************************
  //FakeEleTree
  //*****************************************************************************************
  ofstream eventListFile("ROCEventList.txt");
  TTree *FakeEleTree = getTreeFromFile(FakeElectronFile.c_str(),"Electrons");
  assert(FakeEleTree);
  FakeEleTree->SetBranchAddress( "weight", &fWeight);
  FakeEleTree->SetBranchAddress( "run", &fRunNumber);
  FakeEleTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  FakeEleTree->SetBranchAddress( "event", &fEventNumber);
  FakeEleTree->SetBranchAddress( "pt", &fElePt); 
  FakeEleTree->SetBranchAddress( "eta", &fEleEta); 
  FakeEleTree->SetBranchAddress( "phi", &fElePhi); 
  FakeEleTree->SetBranchAddress( "scEt", &fEleSCEt); 
  FakeEleTree->SetBranchAddress( "scEta", &fEleSCEta); 
  FakeEleTree->SetBranchAddress( "combPFIsoHWW", &fElePFIso); 
  FakeEleTree->SetBranchAddress("matchConv",&fEleMatchedConversion); 
  FakeEleTree->SetBranchAddress("dist",&fEleConvDist); 
  FakeEleTree->SetBranchAddress("dcot",&fEleConvDCot); 
  FakeEleTree->SetBranchAddress("missHits",&fEleNMissHits); 
  FakeEleTree->SetBranchAddress( "see", &fEleSigmaIEtaIEta); 
  FakeEleTree->SetBranchAddress( "deta", &fEleDEtaIn); 
  FakeEleTree->SetBranchAddress( "dphi", &fEleDPhiIn); 
  FakeEleTree->SetBranchAddress( "HoE", &fEleHoverE); 
  FakeEleTree->SetBranchAddress( "d0", &fEleD0); 
  FakeEleTree->SetBranchAddress( "dz", &fEleDZ); 
  FakeEleTree->SetBranchAddress( "fbrem", &fEleFBrem); 
  FakeEleTree->SetBranchAddress( "EoP", &fEleEOverP); 
  FakeEleTree->SetBranchAddress( "EoPout", &fEleESeedClusterOverPout); 
  FakeEleTree->SetBranchAddress( "spp", &fEleSigmaIPhiIPhi); 
  FakeEleTree->SetBranchAddress( "IoEmIoP", &fEleOneOverEMinusOneOverP); 
  FakeEleTree->SetBranchAddress( "EoPin", &fEleESeedClusterOverPIn); 
  FakeEleTree->SetBranchAddress( "ip3d", &fEleIP3d); 
  FakeEleTree->SetBranchAddress( "ip3ds", &fEleIP3dSig); 
  FakeEleTree->SetBranchAddress( "rho", &fRho); 
  FakeEleTree->SetBranchAddress( "vertices", &fNVertices); 

  FakeEleTree->SetBranchAddress("trkIso03",&fEleTrkIso03);
  FakeEleTree->SetBranchAddress("ecalIso03",&fEleEMIso03);
  FakeEleTree->SetBranchAddress("hcalIso03",&fEleHadIso03);

  FakeEleTree->SetBranchAddress( "EleIDMVA_BDTG_V0", &fEleIDBDT_V0 );
  FakeEleTree->SetBranchAddress( "EleIDMVA_BDTG_V1", &fEleIDBDT_V1 );
  FakeEleTree->SetBranchAddress( "EleIDMVA_BDTG_V2", &fEleIDBDT_V2 );
  FakeEleTree->SetBranchAddress( "EleIsoMVA_BDTG_V0", &fEleIsoBDT); 
  FakeEleTree->SetBranchAddress( "EleIDIsoMVA_BDTG_V0", &fEleIDIsoBDT_V0); 
  FakeEleTree->SetBranchAddress( "EleIDIsoMVA_BDTG_V1", &fEleIDIsoBDT_V1); 
 
  FakeEleTree->SetBranchAddress( "DenomFake", &fEleDenomFake); 

  for(UInt_t ientry=0; ientry < FakeEleTree->GetEntries(); ientry++) {       	
    FakeEleTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    Double_t rho = 0;
    if (!(TMath::IsNaN(fRho) || isinf(fRho))) rho = fRho;

    //don't evaluate performance using training events
    if (fElePt < 5) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fEleEta) < 0.8) subdet = 0;
    else if (fabs(fEleEta) < 1.485) subdet = 1;
    else subdet = 2;

    Int_t ptBin = 0;
    if (fElePt > 10.0) ptBin = 1;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);    
    if (Option == 10) passCuts = (ptBin == 0 );
    if (!passCuts) continue;    


    if( fEleIP3dSig > 4 ) continue;
    if( fEleTrkIso03/fElePt > 0.7 ) continue;
    if( fEleNMissHits >  1 ) continue;

    eventListFile << fRunNumber << " " << fLumiSectionNumber << " " << fEventNumber << " : "
                  << fElePt << " " << fEleEta << " " << fElePhi << " "
                  << PassCiCID(fElePt, fEleSCEt, fEleEta, 1, Bool_t(fabs(fEleEta) < 1.479), 
                               fEleFBrem, fEleEOverP, fEleHoverE, fEleSigmaIEtaIEta, fEleDPhiIn, fEleDEtaIn, 
                               fEleESeedClusterOverPIn, 0, 
                               0, 0, 0, 
                               0, 0, 
                               fEleD0, fEleIP3dSig, 
                               1) << " "
                  << endl;


    if (fEventNumber == 122200967) {
      cout << "DEBUG\n";
      PassCiCID(fElePt, fEleSCEt, fEleEta, 1, Bool_t(fabs(fEleEta) < 1.479), 
                fEleFBrem, fEleEOverP, fEleHoverE, fEleSigmaIEtaIEta, fEleDPhiIn, fEleDEtaIn, 
                fEleESeedClusterOverPIn, 0, 
                0, 0, 0, 
                0, 0, 
                fEleD0, fEleIP3dSig, 
                1, kTRUE);
      cout << "DEBUG DONE\n";
    }



    FakeElectrons += fWeight;

    Double_t DetIso03PUCorrection = 0;
    if (fabs(fEleEta) < 1.5) DetIso03PUCorrection = rho * (0.078 + 0.026);
    else DetIso03PUCorrection = rho * (0.046 + 0.072);

    //baseline working point
    if (PassCiCID(fElePt, fEleSCEt, fEleEta, 1, Bool_t(fabs(fEleEta) < 1.479), 
                  fEleFBrem, fEleEOverP, fEleHoverE, fEleSigmaIEtaIEta, fEleDPhiIn, fEleDEtaIn, 
                  fEleESeedClusterOverPIn, 0, 
                  0, 0, 0, 
                  0, 0, 
                  fEleD0, fEleIP3dSig, 
                  1)
        &&
        (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.25
      ) {
      FakeElectronPassCiCTightDetIso += fWeight;
    }
    if (PassCiCID(fElePt, fEleSCEt, fEleEta, 1, Bool_t(fabs(fEleEta) < 1.479), 
                  fEleFBrem, fEleEOverP, fEleHoverE, fEleSigmaIEtaIEta, fEleDPhiIn, fEleDEtaIn, 
                  fEleESeedClusterOverPIn, 0, 
                  0, 0, 0, 
                  0, 0, 
                  fEleD0, fEleIP3dSig, 
                  1)
        &&
        (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt < 0.15
      ) {
      FakeElectronPassCiCTightDetIso015 += fWeight;
    }

    //mva working point
    if (passIDIsoMVA( fElePt, fEleEta, fEleIDIsoBDT_V1)) {
      FakeElectronPassIDIsoMVA += fWeight;
    }

    //Separated Working point
    if (passIDAndIsoMVA_Loose(fElePt, fEleEta, fEleIDBDT_V1,fEleIsoBDT)) {
      FakeElectronPassIDMVAAndIsoMVALoose += fWeight;
    }
    if (passIDAndIsoMVA_Tight(fElePt, fEleEta, fEleIDBDT_V1,fEleIsoBDT)) {
      FakeElectronPassIDMVAAndIsoMVATight += fWeight;
    }

    //DetIsolation
    if (PassCiCID(fElePt, fEleSCEt, fEleEta, 1, Bool_t(fabs(fEleEta) < 1.479), 
                  fEleFBrem, fEleEOverP, fEleHoverE, fEleSigmaIEtaIEta, fEleDPhiIn, fEleDEtaIn, 
                  fEleESeedClusterOverPIn, 0, 
                  0, 0, 0, 
                  0, 0, 
                  fEleD0, fEleIP3dSig, 
                  1)) {
//       cout << "Fake DetIso: " << (fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt << endl;
      EleDetIsoWithCiCTight_Fake->Fill(min((fEleTrkIso03 + fEleEMIso03 + fEleHadIso03 - DetIso03PUCorrection)/fElePt, 2.5), fWeight);
    } else {
      EleDetIsoWithCiCTight_Fake->Fill(2.5,fWeight);
    }

   //Scan IDMVA for fixed IsoMVA cut
    if (passIDMVA(fElePt, fEleEta, fEleIDBDT_V2, -3)) {
      EleIsoBDTWithIDBDTCut0_Fake->Fill(fEleIsoBDT,fWeight);
    } else {
      EleIsoBDTWithIDBDTCut0_Fake->Fill(-0.99,fWeight);
    }
    if (passIDMVA(fElePt, fEleEta, fEleIDBDT_V2, -2)) {
      EleIsoBDTWithIDBDTCut1_Fake->Fill(fEleIsoBDT,fWeight);
    } else {
      EleIsoBDTWithIDBDTCut1_Fake->Fill(-0.99,fWeight);
    }
    if (passIDMVA(fElePt, fEleEta, fEleIDBDT_V2, -1)) {
      EleIsoBDTWithIDBDTCut2_Fake->Fill(fEleIsoBDT,fWeight);
    } else {
      EleIsoBDTWithIDBDTCut2_Fake->Fill(-0.99,fWeight);
    }
    if (passIDMVA(fElePt, fEleEta, fEleIDBDT_V2, 0)) {
      EleIsoBDTWithIDBDTCut3_Fake->Fill(fEleIsoBDT,fWeight);
    } else {
      EleIsoBDTWithIDBDTCut3_Fake->Fill(-0.99,fWeight);
    }
    if (passIDMVA(fElePt, fEleEta, fEleIDBDT_V2, 1)) {
      EleIsoBDTWithIDBDTCut4_Fake->Fill(fEleIsoBDT,fWeight);
    } else {
      EleIsoBDTWithIDBDTCut4_Fake->Fill(-0.99,fWeight);
    }
    if (passIDMVA(fElePt, fEleEta, fEleIDBDT_V2, 2)) {
      EleIsoBDTWithIDBDTCut5_Fake->Fill(fEleIsoBDT,fWeight);
    } else {
      EleIsoBDTWithIDBDTCut5_Fake->Fill(-0.99,fWeight);
    }
    //Scan IsoMVA for fixed IDMVA cut    

    if (passIsoMVA(fElePt, fEleEta, fEleIsoBDT, -3)) {
      EleIDBDTWithIsoBDTCut0_Fake->Fill(fEleIDBDT_V2,fWeight);
    } else {
      EleIDBDTWithIsoBDTCut0_Fake->Fill(-0.99,fWeight);
    }
    if (passIsoMVA(fElePt, fEleEta, fEleIsoBDT, -2)) {
      EleIDBDTWithIsoBDTCut1_Fake->Fill(fEleIDBDT_V2,fWeight);
    } else {
      EleIDBDTWithIsoBDTCut1_Fake->Fill(-0.99,fWeight);
    }
    if (passIsoMVA(fElePt, fEleEta, fEleIsoBDT, -1)) {
      EleIDBDTWithIsoBDTCut2_Fake->Fill(fEleIDBDT_V2,fWeight);
    } else {
      EleIDBDTWithIsoBDTCut2_Fake->Fill(-0.99,fWeight);
    }
    if (passIsoMVA(fElePt, fEleEta, fEleIsoBDT, 0)) {
      EleIDBDTWithIsoBDTCut3_Fake->Fill(fEleIDBDT_V2,fWeight);
    } else {
      EleIDBDTWithIsoBDTCut3_Fake->Fill(-0.99,fWeight);
    }
    if (passIsoMVA(fElePt, fEleEta, fEleIsoBDT, 1)) {
      EleIDBDTWithIsoBDTCut4_Fake->Fill(fEleIDBDT_V2,fWeight);
    } else {
      EleIDBDTWithIsoBDTCut4_Fake->Fill(-0.99,fWeight);
    }
    if (passIsoMVA(fElePt, fEleEta, fEleIsoBDT, 2)) {
      EleIDBDTWithIsoBDTCut5_Fake->Fill(fEleIDBDT_V2,fWeight);
    } else {
      EleIDBDTWithIsoBDTCut5_Fake->Fill(-0.99,fWeight);
    }

    Ele_IDMVAV0_Vs_IsoMVA_Fake->Fill(fEleIDBDT_V0,fEleIsoBDT,fWeight);
    Ele_IDMVAV1_Vs_IsoMVA_Fake->Fill(fEleIDBDT_V1,fEleIsoBDT,fWeight);
    Ele_IDMVAV2_Vs_IsoMVA_Fake->Fill(fEleIDBDT_V2,fEleIsoBDT,fWeight);

    //IDIsoCombined
    EleIDIsoBDT_V0_Fake->Fill(fEleIDIsoBDT_V0,fWeight);
    EleIDIsoBDT_V1_Fake->Fill(fEleIDIsoBDT_V1,fWeight);
//     cout << "Fake " << fEleIDIsoBDT << endl;

  } //loop over electrons
  
  eventListFile.close();

  
  //*****************************************************************************************
  //Current Working Points
  //*****************************************************************************************
  cout << "CiC Tight && DetIso < 0.25 : " << RealElectronPassCiCTightDetIso << " / " << RealElectrons << " = " << RealElectronPassCiCTightDetIso/RealElectrons << endl;
  cout << "CiC Tight && DetIso < 0.25 : " << FakeElectronPassCiCTightDetIso << " / " << FakeElectrons << " = " << FakeElectronPassCiCTightDetIso/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_CiCTightDetIso = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassCiCTightDetIso/RealElectrons , FakeElectronPassCiCTightDetIso/FakeElectrons, "ROC_CiCTightDetIso"+label);


  cout << "CiC Tight && DetIso < 0.15 : " << RealElectronPassCiCTightDetIso015 << " / " << RealElectrons << " = " << RealElectronPassCiCTightDetIso015/RealElectrons << endl;
  cout << "CiC Tight && DetIso < 0.15 : " << FakeElectronPassCiCTightDetIso015 << " / " << FakeElectrons << " = " << FakeElectronPassCiCTightDetIso015/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_CiCTightDetIso015 = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassCiCTightDetIso015/RealElectrons , FakeElectronPassCiCTightDetIso015/FakeElectrons, "ROC_CiCTightDetIso015"+label);


  cout << "IDMVA && IsoMVA Loose: " << RealElectronPassIDMVAAndIsoMVALoose << " / " << RealElectrons << " = " << RealElectronPassIDMVAAndIsoMVALoose/RealElectrons << endl;
  cout << "IDMVA && IsoMVA Loose: " << FakeElectronPassIDMVAAndIsoMVALoose << " / " << FakeElectrons << " = " << FakeElectronPassIDMVAAndIsoMVALoose/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_IDMVAAndIsoMVALoose = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassIDMVAAndIsoMVALoose/RealElectrons , FakeElectronPassIDMVAAndIsoMVALoose/FakeElectrons, "ROC_IDMVAAndIsoMVALoose"+label);


  cout << "IDMVA && IsoMVA Tight: " << RealElectronPassIDMVAAndIsoMVATight << " / " << RealElectrons << " = " << RealElectronPassIDMVAAndIsoMVATight/RealElectrons << endl;
  cout << "IDMVA && IsoMVA Tight: " << FakeElectronPassIDMVAAndIsoMVATight << " / " << FakeElectrons << " = " << FakeElectronPassIDMVAAndIsoMVATight/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_IDMVAAndIsoMVATight = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassIDMVAAndIsoMVATight/RealElectrons , FakeElectronPassIDMVAAndIsoMVATight/FakeElectrons, "ROC_IDMVAAndIsoMVATight"+label);


  cout << "2011 MVA Real Electron Efficiency : " << RealElectronPassIDIsoMVA << " / " << RealElectrons << " = " << RealElectronPassIDIsoMVA/RealElectrons << endl;
  cout << "2011 MVA Fake Electron Efficiency : " << FakeElectronPassIDIsoMVA << " / " << FakeElectrons << " = " << FakeElectronPassIDIsoMVA/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_IDMVA = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassIDIsoMVA/RealElectrons , FakeElectronPassIDIsoMVA/FakeElectrons, "ROC_IDMVA"+label);



//   cout << "**********************\n";
//   cout << "Bkg At LHTight Signal Eff\n";

   Double_t SigEffBaseline = RealElectronPassCiCTightDetIso/RealElectrons;
   Double_t BkgEffBaseline = FakeElectronPassCiCTightDetIso/FakeElectrons;
   Double_t SigEffTight = RealElectronPassCiCTightDetIso015/RealElectrons;
   Double_t BkgEffTight = FakeElectronPassCiCTightDetIso015/FakeElectrons;
//   Double_t SigEffBDTGV0_SameBkg = FindSigEffAtFixedBkgEfficiency(EleIDIsoBDT_V0_Real, EleIDIsoBDT_V0_Fake, BkgEffBaseline);
//   Double_t SigEffBDTGV1_SameBkg = FindSigEffAtFixedBkgEfficiency(EleIDIsoBDT_V1_Real, EleIDIsoBDT_V1_Fake, BkgEffBaseline);

//   cout << SigEffBDTGV0_SameBkg << " , " << SigEffBaseline << endl;
//   cout << "Signal Efficiency (wrt Cut-based) for : same bkg \n";
//   cout << "BDTGV0 : " << SigEffBDTGV0_SameBkg/SigEffBaseline <<  endl;
//   cout << "BDTGV1 : " << SigEffBDTGV1_SameBkg/SigEffBaseline <<  endl;


//   cout << "**********************\n";

//   Double_t BkgEffBDTGV0_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDIsoBDT_V0_Real, EleIDIsoBDT_V0_Fake, BkgEffBaseline);
//   Double_t BkgEffBDTGV1_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDIsoBDT_V1_Real, EleIDIsoBDT_V1_Fake, BkgEffBaseline);
//   cout << "Bkg Efficiency (wrt Cut-based) for same sig eff \n";
//   cout << "BDTGV0 : " << BkgEffBDTGV0_SameSig/BkgEffBaseline << endl;
//   cout << "BDTGV1 : " << BkgEffBDTGV1_SameSig/BkgEffBaseline << endl;

//   cout << "**********************\n";



  //*****************************************************************************************
  //Make ROC curves
  //*****************************************************************************************


   TGraphAsymmErrors* ROC_DetIsoWithCiCTight = MakeSigEffVsBkgEffGraph(EleDetIsoWithCiCTight_Real, EleDetIsoWithCiCTight_Fake, "ROC_DetIsoWithCiCTight"+label , kFALSE);
//    TGraphAsymmErrors* ROC_IDIsoMVACombined_V0 = MakeSigEffVsBkgEffGraph(EleIDIsoBDT_V0_Real, EleIDIsoBDT_V0_Fake, "ROC_IDIsoMVACombined_V0"+label );
   TGraphAsymmErrors* ROC_IDIsoMVACombined_V1 = MakeSigEffVsBkgEffGraph(EleIDIsoBDT_V1_Real, EleIDIsoBDT_V1_Fake, "ROC_IDIsoMVACombined_V1"+label );

//   TGraphAsymmErrors* ROC_IDWithIsoCut0 = MakeSigEffVsBkgEffGraph(EleIDBDTWithIsoBDTCut0_Real, EleIDBDTWithIsoBDTCut0_Fake, "ROC_IDWithIsoCut0"+label );
//   TGraphAsymmErrors* ROC_IDWithIsoCut1 = MakeSigEffVsBkgEffGraph(EleIDBDTWithIsoBDTCut1_Real, EleIDBDTWithIsoBDTCut1_Fake, "ROC_IDWithIsoCut1"+label );
//   TGraphAsymmErrors* ROC_IDWithIsoCut2 = MakeSigEffVsBkgEffGraph(EleIDBDTWithIsoBDTCut2_Real, EleIDBDTWithIsoBDTCut2_Fake, "ROC_IDWithIsoCut2"+label );
  TGraphAsymmErrors* ROC_IDWithIsoCut3 = MakeSigEffVsBkgEffGraph(EleIDBDTWithIsoBDTCut3_Real, EleIDBDTWithIsoBDTCut3_Fake, "ROC_IDWithIsoCut3"+label );
//   TGraphAsymmErrors* ROC_IDWithIsoCut4 = MakeSigEffVsBkgEffGraph(EleIDBDTWithIsoBDTCut4_Real, EleIDBDTWithIsoBDTCut4_Fake, "ROC_IDWithIsoCut4"+label );
//   TGraphAsymmErrors* ROC_IDWithIsoCut5 = MakeSigEffVsBkgEffGraph(EleIDBDTWithIsoBDTCut5_Real, EleIDBDTWithIsoBDTCut5_Fake, "ROC_IDWithIsoCut5"+label );

//   TGraphAsymmErrors* ROC_IsoWithIDCut0 = MakeSigEffVsBkgEffGraph(EleIsoBDTWithIDBDTCut0_Real, EleIsoBDTWithIDBDTCut0_Fake, "ROC_IsoWithIDCut0"+label );
//   TGraphAsymmErrors* ROC_IsoWithIDCut1 = MakeSigEffVsBkgEffGraph(EleIsoBDTWithIDBDTCut1_Real, EleIsoBDTWithIDBDTCut1_Fake, "ROC_IsoWithIDCut1"+label );
//   TGraphAsymmErrors* ROC_IsoWithIDCut2 = MakeSigEffVsBkgEffGraph(EleIsoBDTWithIDBDTCut2_Real, EleIsoBDTWithIDBDTCut2_Fake, "ROC_IsoWithIDCut2"+label );
  TGraphAsymmErrors* ROC_IsoWithIDCut3 = MakeSigEffVsBkgEffGraph(EleIsoBDTWithIDBDTCut3_Real, EleIsoBDTWithIDBDTCut3_Fake, "ROC_IsoWithIDCut3"+label );
//   TGraphAsymmErrors* ROC_IsoWithIDCut4 = MakeSigEffVsBkgEffGraph(EleIsoBDTWithIDBDTCut4_Real, EleIsoBDTWithIDBDTCut4_Fake, "ROC_IsoWithIDCut4"+label );
//   TGraphAsymmErrors* ROC_IsoWithIDCut5 = MakeSigEffVsBkgEffGraph(EleIsoBDTWithIDBDTCut5_Real, EleIsoBDTWithIDBDTCut5_Fake, "ROC_IsoWithIDCut5"+label );

//   //*****************************************************************************************
//   //Find Cut with same signal efficiency Make ROC curves
//   //*****************************************************************************************
//   Double_t CutValue_BDTGV0_SameSig = FindCutValueAtFixedEfficiency(EleIDIsoBDT_V0_Real, BkgEffBaseline );
//   Double_t CutValue_BDTGV0_SameBkg = FindCutValueAtFixedEfficiency(EleIDIsoBDT_V0_Fake, BkgEffBaseline );
//   cout << "BDTG V0 Cut Value @ Same Cut-Based Sig: " << CutValue_BDTGV0_SameSig << endl;
//   cout << "BDTG V0 Cut Value @ Same Cut-Based Bkg: " << CutValue_BDTGV0_SameBkg << endl;

  //*****************************************************************************************
  //Find Cut on ID and Iso with same bkg efficiency as Baseline
  //*****************************************************************************************
  Double_t CutValue_IsoMVA_V1_SameBkg = 0;
  Double_t CutValue_IDMVAV1_V1_SameBkg = 0;
  Double_t SigEff_AtSameBkg_V1 = 0;
  OptimizeIDAndIsoCuts(Ele_IDMVAV1_Vs_IsoMVA_Real, Ele_IDMVAV1_Vs_IsoMVA_Fake, CutValue_IDMVAV1_V1_SameBkg, CutValue_IsoMVA_V1_SameBkg, SigEff_AtSameBkg_V1, 0.95*BkgEffBaseline );
  cout << "\n";
  cout << "Target BkgEff: " << 0.95*BkgEffBaseline << endl;
  cout << "CutValues (IDMVA V1, Iso MVA) : " << CutValue_IDMVAV1_V1_SameBkg << " , " << CutValue_IsoMVA_V1_SameBkg << endl;
  cout << "Sig Eff at Same Bkg : " << SigEff_AtSameBkg_V1 << endl;



  Double_t CutValue_IsoMVA_V1_SameBkgDetIso015 = 0;
  Double_t CutValue_IDMVAV1_V1_SameBkgDetIso015 = 0;
  Double_t SigEff_AtSameBkgDetIso015_V1 = 0;
  OptimizeIDAndIsoCuts(Ele_IDMVAV1_Vs_IsoMVA_Real, Ele_IDMVAV1_Vs_IsoMVA_Fake, CutValue_IDMVAV1_V1_SameBkgDetIso015, CutValue_IsoMVA_V1_SameBkgDetIso015, SigEff_AtSameBkgDetIso015_V1, 0.95*BkgEffTight );
  cout << "\n";
  cout << "Tight Target BkgEff: " << 0.95*BkgEffTight << endl;
  cout << "CutValues (IDMVA V1, Iso MVA) : " << CutValue_IDMVAV1_V1_SameBkgDetIso015 << " , " << CutValue_IsoMVA_V1_SameBkgDetIso015 << endl;
  cout << "Sig Eff at Same Bkg : " << SigEff_AtSameBkgDetIso015_V1 << endl;



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


  ROCGraphs.push_back(ROC_DetIsoWithCiCTight);
  GraphLabels.push_back("DetIso With CiCTight");
  colors.push_back(kRed);

//   ROCGraphs.push_back(ROC_IDIsoMVACombined_V0);
//   GraphLabels.push_back("ID+Iso Combined WJets");
//   colors.push_back(kBlue);

//   ROCGraphs.push_back(ROC_IDIsoMVACombined_V1);
//   GraphLabels.push_back("ID+Iso Combined ZJets");
//   colors.push_back(kBlack);
  
//   ROCGraphs.push_back(ROC_IDWithIsoCut0);
//   GraphLabels.push_back("ID+FixedIsoMVACut");
//   colors.push_back(kGreen+2);
//   ROCGraphs.push_back(ROC_IDWithIsoCut1);
//   GraphLabels.push_back("ID+FixedIsoMVACut");
//   colors.push_back(kGreen+2);
//   ROCGraphs.push_back(ROC_IDWithIsoCut2);
//   GraphLabels.push_back("ID+FixedIsoMVACut");
//   colors.push_back(kGreen+2);
//     ROCGraphs.push_back(ROC_IDWithIsoCut3);
//     GraphLabels.push_back("ID+FixedIsoMVACut");
//     colors.push_back(kGreen+2);
//   ROCGraphs.push_back(ROC_IDWithIsoCut4);
//   GraphLabels.push_back("ID+FixedIsoMVACut");
//   colors.push_back(kGreen+2);
//   ROCGraphs.push_back(ROC_IDWithIsoCut5);
//   GraphLabels.push_back("ID+FixedIsoMVACut");
//   colors.push_back(kGreen+2);
  
//   ROCGraphs.push_back(ROC_IsoWithIDCut0);
//   GraphLabels.push_back("Iso+FixedIDMVACut");
//   colors.push_back(kMagenta);
//   ROCGraphs.push_back(ROC_IsoWithIDCut1);
//   GraphLabels.push_back("Iso+FixedIDMVACut");
//   colors.push_back(kMagenta);
//   ROCGraphs.push_back(ROC_IsoWithIDCut2);
//   GraphLabels.push_back("Iso+FixedIDMVACut");
//   colors.push_back(kMagenta);
    ROCGraphs.push_back(ROC_IsoWithIDCut3);
    GraphLabels.push_back("Iso+FixedIDMVACut");
    colors.push_back(kMagenta);
//   ROCGraphs.push_back(ROC_IsoWithIDCut4);
//   GraphLabels.push_back("Iso+FixedIDMVACut");
//   colors.push_back(kMagenta);
//   ROCGraphs.push_back(ROC_IsoWithIDCut5);
//   GraphLabels.push_back("Iso+FixedIDMVACut");
//   colors.push_back(kMagenta);

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
    if (i==0 || i==1 || i==2 || i==3 || i==9) {
      legend->AddEntry(ROCGraphs[i],GraphLabels[i].c_str(), "LP");
    }

    ROCGraphs[i]->SetMarkerColor(colors[i]);
    ROCGraphs[i]->SetLineColor(colors[i]);

    ROCGraphs[i]->SetMarkerSize(0.75);
    if (i>=3) {
      ROCGraphs[i]->SetMarkerSize(0.25);
    }
   
    ROCGraphs[i]->GetXaxis()->SetRangeUser(xmin,xmax);    
    ROCGraphs[i]->GetYaxis()->SetRangeUser(ymin,ymax);    
    if (i==0) {
      ROCGraphs[i]->Draw("AP");
    } else {
      ROCGraphs[i]->Draw("Psame");
    }
  }


  legend->AddEntry(ROC_CiCTightDetIso, "CiCTight+DetIso<0.25", "P");
  ROC_CiCTightDetIso->SetFillColor(kGreen);
  ROC_CiCTightDetIso->SetMarkerColor(kGreen);
  ROC_CiCTightDetIso->SetMarkerStyle(34);
  ROC_CiCTightDetIso->SetMarkerSize(2.5);
  ROC_CiCTightDetIso->Draw("Psame");

  legend->AddEntry(ROC_CiCTightDetIso015, "CiCTight+DetIso<0.15", "P");
  ROC_CiCTightDetIso015->SetFillColor(kGreen+3);
  ROC_CiCTightDetIso015->SetMarkerColor(kGreen+3);
  ROC_CiCTightDetIso015->SetMarkerStyle(34);
  ROC_CiCTightDetIso015->SetMarkerSize(2.5);
  ROC_CiCTightDetIso015->Draw("Psame");

  legend->AddEntry(ROC_IDMVAAndIsoMVALoose, "MVAID&Iso Loose", "P");
  ROC_IDMVAAndIsoMVALoose->SetFillColor(kRed);
  ROC_IDMVAAndIsoMVALoose->SetMarkerColor(kRed);
  ROC_IDMVAAndIsoMVALoose->SetMarkerStyle(34);
  ROC_IDMVAAndIsoMVALoose->SetMarkerSize(2.5);
  ROC_IDMVAAndIsoMVALoose->Draw("Psame");

  legend->AddEntry(ROC_IDMVAAndIsoMVATight, "MVAID&Iso Tight", "P");
  ROC_IDMVAAndIsoMVATight->SetFillColor(kRed+3);
  ROC_IDMVAAndIsoMVATight->SetMarkerColor(kRed+3);
  ROC_IDMVAAndIsoMVATight->SetMarkerStyle(34);
  ROC_IDMVAAndIsoMVATight->SetMarkerSize(2.5);
  ROC_IDMVAAndIsoMVATight->Draw("Psame");


//   legend->AddEntry(ROC_IDMVA, "ID+Iso Combined", "P");
//   ROC_IDMVA->SetFillColor(kBlack);
//   ROC_IDMVA->SetMarkerColor(kBlack);
//   ROC_IDMVA->SetMarkerStyle(34);
//   ROC_IDMVA->SetMarkerSize(2.5);
//   ROC_IDMVA->Draw("Psame");


  legend->Draw();
  
  cv->SaveAs(("ROCGraphs_" + plotname + ".gif").c_str());
//   canvasFile->WriteTObject(cv,("ROCGraphs_" + plotname).c_str(), "WriteDelete") ;


  gBenchmark->Show("WWTemplate");       
} 

