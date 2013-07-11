//================================================================================================
// root -l -b -q ZeeGammaMassFitSystematicStudyTest.C+'("/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/Data2012_EleHZZICHEP2012WPWithZeeGamma/basic2_76_106_LowPt/plots/FitWorkspaceFile_etapt_0.root",1234,0,3)'

//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1D.h>                   // 1D histograms
#include <TH2D.h>                   // 2D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TEfficiency.h>            // class to handle efficiency calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "HiggsAna/Utils/CPlot.hh"	    // helper class for plots
#include "HiggsAna/Utils/PlotStyle.hh"  // style settings for drawing
#include "HiggsAna/Utils/CEffUser1D.hh"     // class for handling efficiency graphs
#include "HiggsAna/Utils/CEffUser2D.hh"     // class for handling efficiency tables

// structure for output ntuple
#include "CITCommon/CommonData/interface/EffData.hh"

#include "CITCommon/FitModels/ZSignals.hh"
#include "CITCommon/FitModels/ZBackgrounds.hh"
#endif

// RooFit headers
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"

#include "TRandom3.h"
#include "TTree.h"

//=== MAIN MACRO ================================================================================================= 

void ZeeGammaMassFitSystematicStudy(string workspaceFile, const Int_t seed = 1234, 
                                    Int_t Option = 0, Int_t NToys = 1) {


  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================    
  TRandom3 *randomnumber = new TRandom3(seed);
//   RooRealVar m("m","mass",60,130);

  RooCategory sample("sample","");
  sample.defineType("Pass",1);
  sample.defineType("Fail",2);

  //--------------------------------------------------------------------------------------------------------------
  //Load Workspace
  //==============================================================================================================    
  TFile *f = new TFile (workspaceFile.c_str(), "READ");
  RooWorkspace *w = (RooWorkspace*)f->Get("MassFitWorkspace");

  //--------------------------------------------------------------------------------------------------------------
  //Setup output tree
  //==============================================================================================================    
  TFile *outputfile = new TFile (Form("EffToyResults_Option%d_Seed%d.root",Option, seed), "RECREATE");
  float varEff = 0;
  float varEffErrL = 0;
  float varEffErrH = 0;
  TTree *outTree = new TTree("eff","eff");
  outTree->Branch("eff",&varEff, "eff/F");
  outTree->Branch("efferrl",&varEffErrL, "efferrl/F");
  outTree->Branch("efferrh",&varEffErrH, "efferrh/F");
  
  //--------------------------------------------------------------------------------------------------------------
  //Load Model
  //==============================================================================================================    
  RooSimultaneous *totalPdf = (RooSimultaneous*)w->pdf("totalPdf");
  RooRealVar *m_default = (RooRealVar*)w->var("m");
  m_default->setRange("signalRange",85, 95);
  
  //get default models
  RooAddPdf *modelPass_default = (RooAddPdf*)w->pdf("modelPass");
  RooAddPdf *modelFail_default = (RooAddPdf*)w->pdf("modelFail");

  //get variables
  RooRealVar *Nsig = (RooRealVar*)w->var("Nsig");
  RooRealVar *eff = (RooRealVar*)w->var("eff");
  RooRealVar *NbkgFail = (RooRealVar*)w->var("NbkgFail");

  RooFormulaVar NsigPass("NsigPass","eff*Nsig",RooArgList(*eff,*Nsig));	 
  RooFormulaVar NsigFail("NsigFail","(1.0-eff)*Nsig",RooArgList(*eff,*Nsig));

  //get number of expected events
  Double_t npass = 100;
  Double_t nfail = 169;

  //*************************************************************************************
  //make alternative model
  //*************************************************************************************
  RooRealVar *tFail_default = (RooRealVar*)w->var("tFail");
  RooRealVar *fracFail_default = (RooRealVar*)w->var("fracFail");
 

  RooRealVar *meanFail_default = (RooRealVar*)w->var("meanFail");
  RooRealVar *sigmaFail_default = (RooRealVar*)w->var("sigmaFail");
  RooHistPdf *bkgFailTemplate_default = (RooHistPdf*)w->pdf("bkgHistPdfFail");
  RooFFTConvPdf *sigFail_default = (RooFFTConvPdf*)w->pdf("signalFail");
  RooFFTConvPdf *bkgFail_default = (RooFFTConvPdf*)w->pdf("bkgConvPdfFail");
  RooExtendPdf *esignalFail_default = (RooExtendPdf *)w->pdf("esignalFail");
  RooExtendPdf *ebackgroundFail_default = (RooExtendPdf *)w->pdf("ebackgroundFail");
  RooExponential *bkgexpFail_default = (RooExponential*)w->pdf("bkgexpFail");
  RooAddPdf *backgroundFail_default = (RooAddPdf*)w->pdf("backgroundFail");
  RooGaussian *bkggausFail_default = (RooGaussian*)w->pdf("bkggausFail");

  //shifted mean
  RooRealVar *meanFail_shifted = new RooRealVar("meanFail_shifted","meanFail_shifted", 0, -5, 5);
  meanFail_shifted->setVal(meanFail_default->getVal());
  if (Option == 1) meanFail_shifted->setVal(meanFail_default->getVal()-1.0);
  else if (Option == 2) meanFail_shifted->setVal(meanFail_default->getVal()+1.0);  
  else if (Option == 11) meanFail_shifted->setVal(meanFail_default->getVal()-2.0);
  else if (Option == 12) meanFail_shifted->setVal(meanFail_default->getVal()+2.0);  

  RooRealVar *sigmaFail_shifted = new RooRealVar("sigmaFail_shifted","sigmaFail_shifted", 0, -5, 5);
  sigmaFail_shifted->setVal(sigmaFail_default->getVal());
  if (Option == 3) sigmaFail_shifted->setVal(sigmaFail_default->getVal()*1.2);
  else if (Option == 4) sigmaFail_shifted->setVal(sigmaFail_default->getVal()*0.8);

  CMCBkgTemplateConvGaussianPlusExp *bkgFailModel = new CMCBkgTemplateConvGaussianPlusExp(*m_default,bkgFailTemplate_default,false,meanFail_shifted,sigmaFail_shifted, "shifted");
  bkgFailModel->t->setVal(tFail_default->getVal());
  bkgFailModel->frac->setVal(fracFail_default->getVal());

  cout << "mean : " << meanFail_default->getVal() << " - " << meanFail_shifted->getVal() << endl;
  cout << "sigma : " << sigmaFail_default->getVal() << " - " << sigmaFail_shifted->getVal() << endl;
  cout << "t: " << tFail_default->getVal() << " - " << bkgFailModel->t->getVal() << endl;
  cout << "frac: " << fracFail_default->getVal() << " - " << bkgFailModel->frac->getVal() << endl;
  
  cout << "eff: " << eff->getVal() << " : " << NsigPass.getVal() << " / " << (NsigPass.getVal() + NsigFail.getVal()) << endl;
  cout << "NbkgFail: " << NbkgFail->getVal() << endl;


  //make alternative fail model
  RooAddPdf *modelFail=0;
  RooExtendPdf *esignalFail=0, *ebackgroundFail=0;
  ebackgroundFail = new RooExtendPdf("ebackgroundFail_shifted","ebackgroundFail_shifted",*(bkgFailModel->model),*NbkgFail,"signalRange");
  modelFail       = new RooAddPdf("modelFail","Model for FAIL sample", RooArgList(*esignalFail_default,*ebackgroundFail));



  cout << "*************************************\n";
  ebackgroundFail->Print();
  cout << "*************************************\n";
  ebackgroundFail_default->Print();
  cout << "*************************************\n";
  modelFail->Print();
  cout << "*************************************\n";
  modelFail_default->Print();
  cout << "*************************************\n";

  TCanvas *cv = new TCanvas("cv","cv",800,600);

  RooPlot *mframeFail_default = m_default->frame(Bins(Int_t(130-60)/2));
  modelFail_default->plotOn(mframeFail_default);
  modelFail_default->plotOn(mframeFail_default,Components("ebackgroundFail"),LineStyle(kDashed),LineColor(kRed));
  modelFail_default->plotOn(mframeFail_default,Components("bkgexpFail"),LineStyle(kDashed),LineColor(kGreen+2));
  mframeFail_default->Draw();
  cv->SaveAs("DefaultModel.gif");

  RooPlot *mframeFail = m_default->frame(Bins(Int_t(130-60)/2));
  modelFail->plotOn(mframeFail);
  modelFail->plotOn(mframeFail,Components("ebackgroundFail_shifted"),LineStyle(kDashed),LineColor(kRed));
  modelFail->plotOn(mframeFail,Components("bkgexpFail_shifted"),LineStyle(kDashed),LineColor(kGreen+2));
  mframeFail->Draw();
  cv->SaveAs(Form("ShiftedModel_%d.gif",Option));


  //*************************************************************************************
  //Do Toys
  //*************************************************************************************
  for(uint t=0; t < NToys; ++t) {

    RooDataSet *pseudoData_pass    = modelPass_default->generate(*m_default, randomnumber->Poisson(npass));
    RooDataSet *pseudoData_fail  = 0;
    pseudoData_fail    = modelFail->generate(*m_default, randomnumber->Poisson(nfail));
    RooDataSet *pseudoDataCombined = new RooDataSet("pseudoDataCombined","pseudoDataCombined",RooArgList(*m_default),
                                                    RooFit::Index(sample),
                                                    RooFit::Import("Pass",*pseudoData_pass),
                                                    RooFit::Import("Fail",*pseudoData_fail));

    pseudoDataCombined->write(Form("toy%d.txt",t));

    RooFitResult *fitResult=0;
    fitResult = totalPdf->fitTo(*pseudoDataCombined,
                                RooFit::Extended(),
                                RooFit::Strategy(2),
                                //RooFit::Minos(RooArgSet(eff)),
                                RooFit::Save());

    cout << "\n\n";
    cout << "Eff Fit: " << eff->getVal() << " -" << fabs(eff->getErrorLo()) << " +" << eff->getErrorHi() << endl;

    //Fill Tree
    varEff = eff->getVal();
    varEffErrL = fabs(eff->getErrorLo());
    varEffErrH = eff->getErrorHi();
    outTree->Fill();


//   //*************************************************************************************
//   //Plot Toys
//   //*************************************************************************************
//   TCanvas *cv = new TCanvas("cv","cv",800,600);
//   char pname[50];
//   char binlabelx[100];
//   char binlabely[100];
//   char yield[50];
//   char effstr[100];
//   char nsigstr[100];
//   char nbkgstr[100];
//   char chi2str[100];

//   //
//   // Plot passing probes
//   //

//   RooPlot *mframeFail_default = m.frame(Bins(Int_t(130-60)/2));
//   modelFail_default->plotOn(mframeFail_default);
//   modelFail_default->plotOn(mframeFail_default,Components("ebackgroundFail"),LineStyle(kDashed),LineColor(kRed));
//   modelFail_default->plotOn(mframeFail_default,Components("bkgexpFail"),LineStyle(kDashed),LineColor(kGreen+2));
//   mframeFail_default->Draw();
//   cv->SaveAs("DefaultModel.gif");





//   RooPlot *mframeFail = m.frame(Bins(Int_t(130-60)/2));
//   modelFail->plotOn(mframeFail);
//   modelFail->plotOn(mframeFail,Components("ebackgroundFail_shifted"),LineStyle(kDashed),LineColor(kRed));
//   modelFail->plotOn(mframeFail,Components("bkgexpFail_shifted"),LineStyle(kDashed),LineColor(kGreen+2));

//   sprintf(yield,"%u Events",(Int_t)passTree->GetEntries());
//   sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",NsigPass.getVal(),NsigPass.getPropagatedError(*fitResult));
//     plotPass.AddTextBox(yield,0.21,0.76,0.51,0.80,0,kBlack,-1);    
//   plotPass.AddTextBox(effstr,0.70,0.85,0.94,0.90,0,kBlack,-1);
//     plotPass.AddTextBox(0.70,0.73,0.94,0.83,0,kBlack,-1,1,nsigstr);//,chi2str);

//   mframeFail->Draw();
//   cv->SaveAs(Form("ShiftedModel_%d.gif",Option));


 
//   //
//   // Plot failing probes
//   //
//   sprintf(pname,"fail%s_%i",name.Data(),ibin);
//   sprintf(yield,"%u Events",(Int_t)failTree->GetEntries());
//   sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",NsigFail.getVal(),NsigFail.getPropagatedError(*fitResult));
//   sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",NbkgFail.getVal(),NbkgFail.getPropagatedError(*fitResult));
//   sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframePass->chiSquare(nflfail));
//   CPlot plotFail(pname,mframeFail,"Failing probes","tag-probe mass [GeV/c^{2}]","Events / 2.0 GeV/c^{2}");
//   plotFail.AddTextBox(binlabelx,0.21,0.85,0.51,0.90,0,kBlack,-1);
//   if((name.CompareTo("etapt")==0) || (name.CompareTo("etaphi")==0)) {
//     plotFail.AddTextBox(binlabely,0.21,0.80,0.51,0.85,0,kBlack,-1);    
//     plotFail.AddTextBox(yield,0.21,0.76,0.51,0.80,0,kBlack,-1);    
//   } else {
//     plotFail.AddTextBox(yield,0.21,0.81,0.51,0.85,0,kBlack,-1);
//   }
//   plotFail.AddTextBox(effstr,0.70,0.85,0.94,0.90,0,kBlack,-1);  
//   plotFail.AddTextBox(0.70,0.68,0.94,0.83,0,kBlack,-1,2,nsigstr,nbkgstr);//,chi2str);
//   plotFail.Draw(cfail,kTRUE,format);  




  } //for loop over all toys
  



  //*************************************************************************************
  //Save To File
  //*************************************************************************************
  outputfile->WriteTObject(outTree, outTree->GetName(), "WriteDelete");

}

