//root -l Smurf/ReferenceAnalyses/ZllNormalization.C+\(\"/data/smurf/sixie/data/Run2011_Summer11_EPSHZZV0/mitf-alljets/data_2l.goodlumi1092ipb.root\",\"/data/smurf/sixie/data/Run2011_Summer11_EPSHZZV0/mitf-alljets/zll.root\",1.143\)

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <TClonesArray.h>          
#include "TRandom.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include <iomanip>
#include <TMath.h>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Smurf/Core/SmurfTree.h"
#include "Smurf/Core/LeptonScaleLookup.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

#include "CITCommon/CommonData/interface/ZeeEventTree.h"
#include "CITHZZ/CommonCode/CommonDefs.hh"


int    verboseLevel =   0;
const double sigmaB = 0.35;

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


//------------------------------------------------------------------------------
// Fill MC PU Distributions
//------------------------------------------------------------------------------
void FillMCPileupDistribution
(
 string mcInputFile    = "/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/dymm.root",
 string outputLabel    = "SmurfV6_dymm"
 )
{
  string outputlabel = ""; if (outputLabel != "") outputlabel = "_" + outputLabel;
  TH1F *NPU = 0;
  TH1F *NPV = 0;

  TFile *infile = new TFile(mcInputFile.c_str(),"READ");
  if (infile->IsOpen()) {
    NPU = (TH1F*)infile->Get("hNPU"); 
    NPV = (TH1F*)infile->Get("hNVtx"); 
  }

  if (!NPU || !NPV) {
    cout << "The NPU and NVtx histograms could not be found from file " << mcInputFile << endl;
    return;
  }

  NormalizeHist(NPU);
  NormalizeHist(NPV);

  //*****************************************************************************
  // Data Loop
  //*****************************************************************************
  TFile *f = new TFile("PileupSource.root", "UPDATE");
  f->WriteTObject(NPU, ("NPU"+outputlabel).c_str(), "WriteDelete");
  f->WriteTObject(NPV, ("NPV"+outputlabel).c_str(), "WriteDelete");
  f->Close();
  delete f;

  infile->Close();
  delete infile;
}

//------------------------------------------------------------------------------
// Fill MC PU Distributions
//------------------------------------------------------------------------------
void FillMCPileupDistributionNtuple
(
 string mcInputFile    = "PileupNtuple_DYmm.root",
 string outputLabel    = "SmurfV6_dymm"
 )
{
  string outputlabel = ""; if (outputLabel != "") outputlabel = "_" + outputLabel;

  TFile *infile = new TFile(mcInputFile.c_str(),"READ");
  TH1F *NPU = (TH1F*)infile->Get("hNPU"); 
  TH1F *NPV = (TH1F*)infile->Get("hNVtx"); 

  NormalizeHist(NPU);
  NormalizeHist(NPV);

  //*****************************************************************************
  // Data Loop
  //*****************************************************************************
  TFile *f = new TFile("PileupReweighting.root", "UPDATE");
  f->WriteTObject(NPU, ("NPU"+outputlabel).c_str(), "WriteDelete");
  f->WriteTObject(NPV, ("NPV"+outputlabel).c_str(), "WriteDelete");
  f->Close();
  delete f;

  infile->Close();
  delete infile;
}



void ComputeWeights(string TargetFileName, string SourceHistName, 
                    string Label, Int_t lastPileupBin = 35) {

  TH1F *targetNPU = 0;
  TH1F *sourceNPU = 0;

  TFile *targetfile = new TFile(TargetFileName.c_str(), "READ");
  if (targetfile->IsOpen()) targetNPU = (TH1F*)targetfile->Get("pileup");
  TFile *sourcefile = new TFile("/afs/cern.ch/work/s/sixie/public/Pileup/PileupSource.root", "READ");
  if (sourcefile->IsOpen()) sourceNPU = (TH1F*)sourcefile->Get(SourceHistName.c_str());

  if (!targetNPU || !sourceNPU) {
    if (!sourceNPU) {
      cout << "Cannot find source hist : " << SourceHistName << "\n"; 
    }
    if (!targetNPU) {
      cout << "Cannot find target hist (\"pileup\") in file : " << TargetFileName << "\n"; 

    }
    return;
  }

  targetNPU->SetDirectory(0);
  NormalizeHist(targetNPU);
  NormalizeHist(sourceNPU);


  TH1D *ReweightingFactor = (TH1D*)targetNPU->Clone("puWeights");
  ReweightingFactor->SetBinContent(0,1.0);

  assert(lastPileupBin <= ReweightingFactor->GetXaxis()->GetNbins());

  for(UInt_t b=1; b < lastPileupBin; ++b) {
    if (sourceNPU->GetBinContent(b)>0) {
      ReweightingFactor->SetBinContent(b, targetNPU->GetBinContent(b) / sourceNPU->GetBinContent(b));
    } else {
      cout << "Error: Bin " << b << " in source histogram has a 0 entry.\n";
      assert(0);
    }
  }

  Double_t sourceLastBin = 0;
  Double_t targetLastBin = 0;
  for(UInt_t b=lastPileupBin; b < ReweightingFactor->GetXaxis()->GetNbins()+2; ++b) {
    sourceLastBin += sourceNPU->GetBinContent(b);
    targetLastBin += targetNPU->GetBinContent(b);
  }
  if (sourceLastBin>0) {
    for (UInt_t b=lastPileupBin; b < ReweightingFactor->GetXaxis()->GetNbins()+2; ++b) {
      ReweightingFactor->SetBinContent(b, targetLastBin / sourceLastBin);
    }
  } else {
    cout << "Error: Last Bin in source histogram is 0.\n";
    assert(0);
  }

  TFile *outputfile = new TFile(("PileupReweighting." + Label + ".root").c_str(), "UPDATE");
  outputfile->WriteTObject(ReweightingFactor, ReweightingFactor->GetName(), "WriteDelete");
  outputfile->Close();
  delete outputfile;
  sourcefile->Close();
  delete sourcefile;
  targetfile->Close();
  delete targetfile;


  TFile *file = new TFile("PileupTargets.root", "UPDATE");
  file->WriteTObject(targetNPU, ("NPU_Target_" + Label).c_str(), "WriteDelete");
  file->Close();
  delete file;

}












void DrawValidationPlots(string Label) {
  string label = ""; if (Label != "") label = "_" + Label;

  TFile *f = new TFile("PileupReweightingValidation.root", "READ");
 
  TH1F* NVtx_MC = (TH1F*)f->Get(("NPV_MC"+label).c_str());
  TH1F* NVtx_Data = (TH1F*)f->Get(("NPV_Data"+label).c_str());
  TH1F* Rho_MC = (TH1F*)f->Get(("Rho_MC"+label).c_str());
  TH1F* Rho_Data = (TH1F*)f->Get(("Rho_Data"+label).c_str());


  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(0);
  tmpLegend->SetFillStyle(0);


  if (Rho_MC && Rho_Data) {
    tmpLegend->Clear();
    tmpLegend->AddEntry(Rho_MC, "Reweighted MC", "L");  
    tmpLegend->AddEntry(Rho_Data, "Data", "L");
    
    Rho_MC->SetMarkerColor(kRed);
    Rho_MC->SetLineColor(kRed);
    Rho_MC->SetTitle("");
    Rho_MC->SetMaximum(max(double(Rho_MC->GetMaximum()),double(Rho_Data->GetMaximum())) * 1.2);
    Rho_MC->GetYaxis()->SetTitleOffset(1.1);
    Rho_MC->GetXaxis()->SetTitleOffset(1.05);
    
    Rho_MC->GetYaxis()->SetTitle("Fraction of Events");
    Rho_MC->Draw("hist");
    Rho_Data->Draw("hist,same");
    tmpLegend->Draw();
    
    cv->SaveAs(("PileupReweightingValidation_Rho" + label + ".png").c_str());
  }

  if (NVtx_MC && NVtx_Data) {
    tmpLegend->Clear();
    tmpLegend->AddEntry(NVtx_MC, "Reweighted MC", "L");  
    tmpLegend->AddEntry(NVtx_Data, "Data", "L");
    
    NVtx_MC->SetMarkerColor(kRed);
    NVtx_MC->SetLineColor(kRed);
    NVtx_MC->SetTitle("");
    NVtx_MC->SetMaximum(max(double(NVtx_MC->GetMaximum()),double(NVtx_Data->GetMaximum())) * 1.2);
    NVtx_MC->GetYaxis()->SetTitleOffset(1.1);
    NVtx_MC->GetXaxis()->SetTitleOffset(1.05);
    
    NVtx_MC->GetYaxis()->SetTitle("Fraction of Events");
    NVtx_MC->Draw("hist");
    NVtx_Data->Draw("hist,same");
    tmpLegend->Draw();
    
    cv->SaveAs(("PileupReweightingValidation_NVtx" + label + ".png").c_str());
  }



}




void DrawPileupPlots() {

  //*************************************************************************
  //NPU Distributions
  //*************************************************************************
  TFile *sourceFile = new TFile("/afs/cern.ch/work/s/sixie/public/Pileup/PileupSource.root", "READ");
  TFile *targetFile = new TFile("PileupTargets.root", "READ");
 
  TH1F* NPU_MC = (TH1F*)sourceFile->Get("NPU_Fall11DYmm");
  TH1F* NPU_Target_2011A    = (TH1F*)targetFile->Get("NPU_Target_Fall11DYmm_To_Run2011A");
  TH1F* NPU_Target_2011B    = (TH1F*)targetFile->Get("NPU_Target_Fall11DYmm_To_Run2011B");
  TH1F* NPU_Target_Full2011 = (TH1F*)targetFile->Get("NPU_Target_Fall11DYmm_To_Full2011");

  TCanvas *cv = 0;
  TLegend *tmpLegend = 0;

  if (NPU_MC && NPU_Target_2011A && NPU_Target_2011B && NPU_Target_Full2011) {
    cv = new TCanvas("cv","cv", 800,600);
    tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
    tmpLegend->SetTextSize(0.04);
    tmpLegend->SetBorderSize(0);


    tmpLegend->Clear();
    tmpLegend->AddEntry(NPU_MC,              "Fall 11 MC"   , "L");  
    tmpLegend->AddEntry(NPU_Target_2011A,    "Target Run2011A", "L");  
    tmpLegend->AddEntry(NPU_Target_2011B,    "Target Run2011B", "L");  
    tmpLegend->AddEntry(NPU_Target_Full2011, "Target Full2011", "L");  
  
    NPU_MC->SetMarkerColor(kRed);
    NPU_MC->SetLineColor(kRed);
    NPU_MC->SetTitle("");
    NPU_MC->SetMaximum(max(max(max(double(NPU_MC->GetMaximum()),double(NPU_Target_2011A->GetMaximum())),
                               NPU_Target_2011B->GetMaximum()),NPU_Target_Full2011->GetMaximum()) * 1.2);
    NPU_MC->GetYaxis()->SetTitleOffset(1.1);
    NPU_MC->GetXaxis()->SetTitleOffset(1.05);
    NPU_MC->GetXaxis()->SetRangeUser(-0.5,49.5);
 
    NPU_Target_2011A->SetLineColor(kBlue);
    NPU_Target_2011B->SetLineColor(kMagenta);
    NPU_Target_Full2011->SetLineColor(kBlack);
    NPU_MC->SetLineWidth(2);
    NPU_Target_2011A->SetLineWidth(2);
    NPU_Target_2011B->SetLineWidth(2);
    NPU_Target_Full2011->SetLineWidth(2);
 
 
    NPU_MC->Draw("hist");
    NPU_Target_2011A->Draw("hist,same");
    NPU_Target_2011B->Draw("hist,same");
    NPU_Target_Full2011->Draw("hist,same");
    tmpLegend->Draw(); 

    cv->SaveAs("NPUDistributions.png");
  } else {
    cout << "Could not retrieve all the source or target pileup distributions.\n";
  }


  //*************************************************************************
  //Weights
  //*************************************************************************
  TFile *f = 0;
  TH1F* Weights_Run2011A = 0;
  TH1F* Weights_Run2011B = 0;
  TH1F* Weights_Full2011 = 0;

  f = new TFile("PileupReweighting.Summer11_To_Run2011A.root", "READ");
  Weights_Run2011A = (TH1F*)f->Get("puWeights");
  f = new TFile("PileupReweighting.Summer11_To_Run2011B.root", "READ");
  Weights_Run2011B = (TH1F*)f->Get("puWeights");
  f = new TFile("PileupReweighting.Summer11_To_Full2011.root", "READ");
  Weights_Full2011 = (TH1F*)f->Get("puWeights");

  if (Weights_Run2011A && Weights_Run2011B && Weights_Full2011) { 
    tmpLegend = new TLegend(0.7,0.75,0.93,0.90);   
    tmpLegend->SetTextSize(0.04);
    tmpLegend->SetBorderSize(0);

    tmpLegend->Clear();
    tmpLegend->AddEntry(Weights_Run2011A,    "Run2011A", "L");  
    tmpLegend->AddEntry(Weights_Run2011B,    "Run2011B", "L");  
    tmpLegend->AddEntry(Weights_Full2011,    "Full2011", "L");  
  
    Weights_Run2011A->SetMarkerColor(kBlue);
    Weights_Run2011A->SetLineColor(kBlue);
    Weights_Run2011A->SetTitle("");
    Weights_Run2011A->SetMaximum(max(max(double(Weights_Run2011A->GetMaximum()),double(Weights_Run2011B->GetMaximum())),
                                     Weights_Full2011->GetMaximum()) * 1.2);
    Weights_Run2011A->GetYaxis()->SetTitleOffset(1.1);
    Weights_Run2011A->GetXaxis()->SetTitleOffset(1.05);
    Weights_Run2011A->GetXaxis()->SetRangeUser(-0.5,49.5);
 
    Weights_Run2011B->SetLineColor(kMagenta);
    Weights_Full2011->SetLineColor(kBlack);
    Weights_Run2011A->SetLineWidth(2);
    Weights_Run2011B->SetLineWidth(2);
    Weights_Full2011->SetLineWidth(2);

    Weights_Run2011A->Draw("hist");
    Weights_Run2011B->Draw("hist,same");
    Weights_Full2011->Draw("hist,same");
    tmpLegend->Draw();

    cv->SaveAs("ReweightingFactors_Summer11Source.png");
  } else {
    cout << "Could not find all pileup reweighting histograms for Summer11 -> Run2011. \n";
  }


  f = new TFile("/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedHalfBkgWP/auxiliar/PileupReweighting.Fall11DYmm_To_Run2011A.root", "READ");
  Weights_Run2011A = (TH1F*)f->Get("puWeights");
  f = new TFile("/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedHalfBkgWP/auxiliar/PileupReweighting.Fall11DYmm_To_Run2011B.root", "READ");
  Weights_Run2011B = (TH1F*)f->Get("puWeights");
  f = new TFile("/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedHalfBkgWP/auxiliar/PileupReweighting.Fall11DYmm_To_Full2011.root", "READ");
  Weights_Full2011 = (TH1F*)f->Get("puWeights");

  if (Weights_Run2011A && Weights_Run2011B && Weights_Full2011) { 
    tmpLegend = new TLegend(0.7,0.75,0.93,0.90);   
    tmpLegend->SetTextSize(0.04);
    tmpLegend->SetBorderSize(0);


    tmpLegend->Clear();
    tmpLegend->AddEntry(Weights_Run2011A,    "Run2011A", "L");  
    tmpLegend->AddEntry(Weights_Run2011B,    "Run2011B", "L");  
    tmpLegend->AddEntry(Weights_Full2011,    "Full2011", "L");  
  
    Weights_Run2011A->SetMarkerColor(kBlue);
    Weights_Run2011A->SetLineColor(kBlue);
    Weights_Run2011A->SetTitle("");
    Weights_Run2011A->SetMaximum(max(max(double(Weights_Run2011A->GetMaximum()),double(Weights_Run2011B->GetMaximum())),
                                     Weights_Full2011->GetMaximum()) * 1.2);
    Weights_Run2011A->GetYaxis()->SetTitleOffset(1.1);
    Weights_Run2011A->GetXaxis()->SetTitleOffset(1.05);
    Weights_Run2011A->GetXaxis()->SetRangeUser(-0.5,49.5);
 
    Weights_Run2011B->SetLineColor(kMagenta);
    Weights_Full2011->SetLineColor(kBlack);
    Weights_Run2011A->SetLineWidth(2);
    Weights_Run2011B->SetLineWidth(2);
    Weights_Full2011->SetLineWidth(2);

    Weights_Run2011A->Draw("hist");
    Weights_Run2011B->Draw("hist,same");
    Weights_Full2011->Draw("hist,same");
    tmpLegend->Draw();
    cv->SaveAs("ReweightingFactors_Fall11Source.png");
  } else {
    cout << "Could not find all pileup reweighting histograms for Fall11 -> Run2011. \n";
  }

}




void ValidateReweighting( string MCFilename, string DataFilename, string ReweightFilename,
                          string label ) {

  string Label = label;
  if (label != "") Label = "_"+label;

  // ***********************************************************************************************
  // Load Reweight File
  // ***********************************************************************************************
  TFile *PUReweightFile = TFile::Open(ReweightFilename.c_str());
  TH1D *PUReweightHist = (TH1D*)(PUReweightFile->Get("puWeights"));
  assert(PUReweightHist);
  PUReweightHist->SetDirectory(0);
  delete PUReweightFile;


  //*************************************************************
  //Histograms
  //*************************************************************
  TH1F *NPV_MC = new TH1F(("NPV_MC"+Label).c_str(), ";Number of reconstructed primary vertices; NEvents; ", 50, -0.5, 49.5);
  TH1F *NPV_Data = new TH1F(("NPV_Data"+Label).c_str(), ";Number of reconstructed primary vertices; NEvents; ", 50, -0.5, 49.5);
  TH1F *Rho_MC = new TH1F(("Rho_MC"+Label).c_str(), ";#rho (Energy density) [GeV]; NEvents; ", 50, -0.5, 49.5);
  TH1F *Rho_Data = new TH1F(("Rho_Data"+Label).c_str(), ";#rho (Energy density) [GeV]; NEvents; ", 50, -0.5, 49.5);


  // Load Z Tree
  citana::ZeeEventTree *dataZeeTree = new citana::ZeeEventTree();
  dataZeeTree->LoadTree(DataFilename.c_str(), citana::ZeeEventTree::kCITZeeEvent);
  citana::ZeeEventTree *mcZeeTree = new citana::ZeeEventTree();
  mcZeeTree->LoadTree(MCFilename.c_str(), citana::ZeeEventTree::kCITZeeEvent);

  //*************************************************************
  //Loop Over Data
  //*************************************************************
  for (int i = 0; i < dataZeeTree->tree_->GetEntries(); i++) {
    dataZeeTree->tree_->GetEntry(i);

    //*************************************************************************
    //Z Selection
    //*************************************************************************
    if (!(dataZeeTree->fEle1Pt > 20 && dataZeeTree->fEle2Pt > 20)) continue;
    if (!(dataZeeTree->fEle1PassHZZICHEP2012 == 1 && dataZeeTree->fEle2PassHZZICHEP2012 == 1)) continue;
    if (!(dataZeeTree->fMass > 75 && dataZeeTree->fMass < 105 )) continue;
    
    //*************************************************************************
    //Fill Histograms
    //*************************************************************************
    NPV_Data->Fill(dataZeeTree->fNVertices);
    Rho_Data->Fill(dataZeeTree->fRho);

  }


  //*************************************************************
  //Loop Over MC
  //*************************************************************
  for (int i = 0; i < mcZeeTree->tree_->GetEntries(); i++) {
    mcZeeTree->tree_->GetEntry(i);

    //*************************************************************************
    //Z Selection
    //*************************************************************************
    if (!(mcZeeTree->fEle1Pt > 20 && mcZeeTree->fEle2Pt > 20)) continue;
    if (!(mcZeeTree->fEle1PassHZZICHEP2012 == 1 && mcZeeTree->fEle2PassHZZICHEP2012 == 1)) continue;
    if (!(mcZeeTree->fMass > 75 && mcZeeTree->fMass < 105 )) continue;
    
    //Get PU weight
    double mynpu = mcZeeTree->fNPU;
    Int_t npuxbin = PUReweightHist->GetXaxis()->FindBin(mynpu);
    Double_t weight = PUReweightHist->GetBinContent(npuxbin);

    //*************************************************************************
    //Fill Histograms
    //*************************************************************************
    NPV_MC->Fill(mcZeeTree->fNVertices, weight);
    Rho_MC->Fill(mcZeeTree->fRho, weight);

  }

  NormalizeHist(NPV_MC);
  NormalizeHist(NPV_Data);
  NormalizeHist(Rho_MC);
  NormalizeHist(Rho_Data);

  //*****************************************************************************
  // SavePlots
  //*****************************************************************************
  TFile *f = new TFile("PileupReweightingValidation.root", "UPDATE");
  f->WriteTObject(NPV_MC, NPV_MC->GetName(), "WriteDelete");
  f->WriteTObject(NPV_Data, NPV_Data->GetName(), "WriteDelete");
  f->WriteTObject(Rho_MC, Rho_MC->GetName(), "WriteDelete");
  f->WriteTObject(Rho_Data, Rho_Data->GetName(), "WriteDelete");
  f->Close();


}





void PileupReweighting(Int_t Option = -1) {

  //**************************************************************************
  // Get Source Distributions from MC
  // Provide input file with the NPU histogram for ALL generated events
  //**************************************************************************
  if (Option == 0) {
    cout << "**********************************************************" << endl;
    cout << "Retrieving Source Pileup Distributions from Monte Carlo..." << endl;
    FillMCPileupDistribution("/home/ceballos/condor/old_5x/histo_s12-zllm50-2-v9_all_noskim.root","Summer12DY");
  }
  

  //**************************************************************************
  // Load Target Distributions and produce reweighting histograms
  //**************************************************************************
  if (Option == -1 || Option == 1) {

    //Summer12 -> Run2012AB
//     ComputeWeights("/afs/cern.ch/work/s/sixie/public/Pileup/PileupTarget_190456To196531.obs.69400.root","NPU_Summer12DY","Summer12DY_To_Run2012AB", 70);
    ComputeWeights("/afs/cern.ch/work/s/sixie/public/Pileup/PileupTarget_190456To196531.obs.73500.root","NPU_Summer12DY","Summer12DY_To_Run2012AB", 70);

    //Summer12 -> Run2012C
    ComputeWeights("/afs/cern.ch/work/s/sixie/public/Pileup/PileupTarget_198049To200601.obs.69400.root","NPU_Summer12DY","Summer12DY_To_Run2012C", 70);

    //Summer12 -> Run2012ABC
    ComputeWeights("/afs/cern.ch/work/s/sixie/public/Pileup/PileupTarget_190456To200601.obs.69400.root","NPU_Summer12DY","Summer12DY_To_Run2012ABC", 70);

  }



  //**************************************************************************
  // Produce Validation plots of reweighting : From SMURF
  //**************************************************************************
  if (Option == -1 || Option == 2) {
    //Run2012AB
    ValidateReweighting("/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY.root",
                        "/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Data2012AB.root",
                        "PileupReweighting.Summer12DY_To_Run2012AB.root",
                        "Summer12DY_2012AB");
    DrawValidationPlots("Summer12DY_2012AB");    


//     //Run2012C
//     ValidateReweighting("/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY.root",
//                         "/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Data2012C.root",
//                         "PileupReweighting.Summer12DY_To_Run2012C.root",
//                         "Summer12DY_2012C");
//     DrawValidationPlots("Summer12DY_2012C");    


//     //Run2012ABC
//     ValidateReweighting("/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY.root",
//                         "/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Data2012ABC.root",
//                         "PileupReweighting.Summer12DY_To_Run2012ABC.root",
//                         "Summer12DY_2012ABC");
//     DrawValidationPlots("Summer12DY_2012ABC");    


  }
  
//   //**************************************************************************
//   // Plots for Note
//   //**************************************************************************
//   if (Option == -1 || Option == 3) {
//     DrawPileupPlots();
//   }

}
 
