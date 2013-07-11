//root -l MakePileupObserved.C+'("PileupTarget_198049To200601.true.root","PileupTarget_198049To200601.obs.root")'

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for 4-vector calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "TH1D.h"
#include "TCanvas.h"
#include "RooPoisson.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TMath.h"

using namespace RooFit;


#endif


void NormalizeHist(TH1D *hist) {
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


void MakePileupObserved(string inputfilename, string outputfilename) {

  TFile *file = new TFile ( inputfilename.c_str(), "read");
  TH1D *histTruth = (TH1D*)file->Get("pileup");
  NormalizeHist(histTruth);
  TH1F *histObs = new TH1F( "pileup", " ; NPU; Number of Events", 201, -0.5, 200.5);

  vector<double> PoissonMeans;
  vector<double> Weights;

  for (UInt_t b=1; int(b)<histTruth->GetXaxis()->GetNbins()+1; ++b) {
    double tmpmean = histTruth->GetXaxis()->GetBinCenter(b);
    double tmpnorm = histTruth->GetBinContent(b);
    cout << "bin " << b << " : " << tmpmean << " : " << tmpnorm << endl;
    PoissonMeans.push_back(tmpmean);
    Weights.push_back(tmpnorm);

  }

  for (UInt_t b=1; int(b)<histObs->GetXaxis()->GetNbins()+1; ++b) {
    double x = histObs->GetXaxis()->GetBinCenter(b);

    //compute the sum of poissons
    double value = 0;
    for (uint k=0; k < PoissonMeans.size(); ++k) {
      value +=  Weights[k]*TMath::Poisson(x, PoissonMeans[k]);
    }
    histObs->SetBinContent(b,value);
    cout << "Obs: " << x << " : " << value << endl;
  }
  NormalizeHist(histObs);

  TCanvas *cv = new TCanvas("cv","cv", 800, 600);
  histObs->Draw("hist");
  histTruth->SetLineColor(kRed);
  histTruth->Draw("hist,same");
  cv->SaveAs("PileupObserved.gif");

  TFile *fileoutput = new TFile(outputfilename.c_str(), "RECREATE");
  fileoutput->WriteTObject(histObs, "pileup", "WriteDelete");
  fileoutput->Close();
  delete fileoutput;

}
