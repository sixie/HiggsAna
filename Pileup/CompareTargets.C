void CompareTargets() {

  TFile *f = new TFile("/data/smurf/sixie/Pileup/PUTarget.Run2011AMay10ReReco.160404-163869.root");
  TH1F* Target_minlumiCut0p1_73p5mb = (TH1F*)f->Get("pileup");
  f = new TFile("/data/smurf/sixie/Pileup/68mb/PUTarget.May10ReReco.160404-163869.68mb.root");
  TH1F* Target_minlumiCut0p005_68mb = (TH1F*)f->Get("pileup");
  f = new TFile("/data/smurf/sixie/Pileup/68mb/PUTarget.May10ReReco.160404-163869.73p5mb_AndreDavidCuts.root");
  TH1F* Target_minlumiCut0p005_73p5mb = (TH1F*)f->Get("pileup");


  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.4,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(1);

  tmpLegend->AddEntry(Target_minlumiCut0p1_73p5mb, "minlumi > 0.1 , 73.5mb" , "L");
  tmpLegend->AddEntry(Target_minlumiCut0p005_68mb, "minlumi > 0.005 , 68mb" , "L");
  tmpLegend->AddEntry(Target_minlumiCut0p005_73p5mb, "minlumi > 0.005 , 73.5mb" , "L");

  Target_minlumiCut0p005_68mb->SetTitle("");
  Target_minlumiCut0p005_68mb->GetXaxis()->SetTitle("Target Number of Pileup Events");
  Target_minlumiCut0p005_68mb->SetLineColor(kRed);
  Target_minlumiCut0p1_73p5mb->SetLineColor(kBlack);
  Target_minlumiCut0p005_73p5mb->SetLineColor(kBlue);
  Target_minlumiCut0p005_68mb->SetLineWidth(2);
  Target_minlumiCut0p1_73p5mb->SetLineWidth(2);
  Target_minlumiCut0p005_73p5mb->SetLineWidth(2);

  Target_minlumiCut0p005_68mb->Draw("hist");
  Target_minlumiCut0p005_73p5mb->Draw("hist,same");
  Target_minlumiCut0p1_73p5mb->Draw("hist,same");

  tmpLegend->Draw();

  cv->SaveAs("May10Rereco.PUTargetComparison.gif");

  }
