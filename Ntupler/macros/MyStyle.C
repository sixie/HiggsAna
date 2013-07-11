#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TStyle.h>
#include <MyStyle.h>

void MyStyle() {
}

TCanvas* MakeCanvas(const char* name, const char *title, int dX, int dY)
{
  // Start with a canvas
  TCanvas *canvas = new TCanvas(name,title,0,0,dX,dY);
  // General overall stuff
  canvas->SetFillColor      (0);
  canvas->SetBorderMode     (0);
  canvas->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canvas->SetLeftMargin     (0.20);
  canvas->SetRightMargin    (0.05);
  canvas->SetTopMargin      (0.08);
  canvas->SetBottomMargin   (0.15);
  // Setup a frame which makes sense
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);

  return canvas;
}

void InitSubPad(TPad* pad, int i)
{
  //printf("Pad: %p, index: %d\n",pad,i);

  pad->cd(i);
  TPad *tmpPad = (TPad*) pad->GetPad(i);
  tmpPad->SetLeftMargin  (0.20);
  tmpPad->SetTopMargin   (0.05);
  tmpPad->SetRightMargin (0.07);
  tmpPad->SetBottomMargin(0.15);
  return;
}

void InitHist(TH1 *hist, const char *xtit, const char *ytit, EColor color)
{
  hist->SetXTitle(xtit);
  hist->SetYTitle(ytit);
  hist->SetLineColor(color);
  hist->SetFillStyle(0);
  hist->SetTitleSize  (0.055,"Y");
  hist->SetTitleOffset(1.400,"Y");
  hist->SetLabelOffset(0.014,"Y");
  hist->SetLabelSize  (0.050,"Y");
  hist->SetLabelFont  (42   ,"Y");
  hist->SetTitleSize  (0.055,"X");
  hist->SetTitleOffset(1.150,"X");
  hist->SetLabelOffset(0.014,"X");
  hist->SetLabelSize  (0.050,"X");
  hist->SetLabelFont  (42   ,"X");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize (0.6);
  // Strangely enough this cannot be set anywhere else??
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->SetTitle("");  
  return;
}

void SetStyle()
{
  TStyle *MYStyle = new TStyle("Style","Style");
  gStyle = MYStyle;

  // Canvas
  MYStyle->SetCanvasColor     (0);
  MYStyle->SetCanvasBorderSize(10);
  MYStyle->SetCanvasBorderMode(0);
  MYStyle->SetCanvasDefH      (700);
  MYStyle->SetCanvasDefW      (700);
  MYStyle->SetCanvasDefX      (100);
  MYStyle->SetCanvasDefY      (100);

  // Pads
  MYStyle->SetPadColor       (0);
  MYStyle->SetPadBorderSize  (10);
  MYStyle->SetPadBorderMode  (0);
  MYStyle->SetPadBottomMargin(0.13);
  MYStyle->SetPadTopMargin   (0.08);
  MYStyle->SetPadLeftMargin  (0.15);
  MYStyle->SetPadRightMargin (0.05);
  MYStyle->SetPadGridX       (0);
  MYStyle->SetPadGridY       (0);
  MYStyle->SetPadTickX       (0);
  MYStyle->SetPadTickY       (0);

  // Frames
  MYStyle->SetFrameFillStyle ( 0);
  MYStyle->SetFrameFillColor ( 0);
  MYStyle->SetFrameLineColor ( 1);
  MYStyle->SetFrameLineStyle ( 0);
  MYStyle->SetFrameLineWidth ( 1);
  MYStyle->SetFrameBorderSize(10);
  MYStyle->SetFrameBorderMode( 0);

  // Histograms
  MYStyle->SetHistFillColor(2);
  MYStyle->SetHistFillStyle(0);
  MYStyle->SetHistLineColor(1);
  MYStyle->SetHistLineStyle(0);
  MYStyle->SetHistLineWidth(2);
  MYStyle->SetNdivisions(505);

  // Functions
  MYStyle->SetFuncColor(1);
  MYStyle->SetFuncStyle(0);
  MYStyle->SetFuncWidth(2);

  // Various
  MYStyle->SetMarkerStyle(20);
  MYStyle->SetMarkerColor(kBlack);
  MYStyle->SetMarkerSize (1.2);

  MYStyle->SetTitleSize  (0.055,"X");
  MYStyle->SetTitleOffset(1.150,"X");
  MYStyle->SetLabelOffset(0.005,"X");
  MYStyle->SetLabelSize  (0.050,"X");
  MYStyle->SetLabelFont  (42   ,"X");

  MYStyle->SetStripDecimals(kFALSE);

  MYStyle->SetTitleSize  (0.055,"Y");
  MYStyle->SetTitleOffset(1.400,"Y");
  MYStyle->SetLabelOffset(0.010,"Y");
  MYStyle->SetLabelSize  (0.050,"Y");
  MYStyle->SetLabelFont  (42   ,"Y");

  MYStyle->SetTextSize   (0.055);
  MYStyle->SetTextFont   (42);

  MYStyle->SetStatFont   (42);
  MYStyle->SetTitleFont  (42);
  MYStyle->SetTitleFont  (42,"X");
  MYStyle->SetTitleFont  (42,"Y");

  MYStyle->SetOptStat    (0);
  return;
}
