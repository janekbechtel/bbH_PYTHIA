#include "TH1.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"

void setTDRStyle(Int_t xlog, Int_t ylog, Int_t zlog) {

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(500); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(2);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(4);
  //  tdrStyle->SetErrorMarker(20);
  //  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

//For the fit/function:
  tdrStyle->SetOptFit(0);
  //  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(1);
  tdrStyle->SetFuncStyle(2);
  tdrStyle->SetFuncWidth(4);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.13);
  tdrStyle->SetPadRightMargin(0.05);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);

  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(0.9);


  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(xlog);
  tdrStyle->SetOptLogy(ylog);
  tdrStyle->SetOptLogz(zlog);

// Postscript options:

//  tdrStyle->SetPaperSize(7.5,7.5);

  tdrStyle->SetPaperSize(15.,15.);

//  tdrStyle->SetPaperSize(20.,20.);

  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();
}

void Draw()
{
  // try RooFit

  gSystem->Load("libRooFit");
  using namespace RooFit;


  TCanvas* c1 = new TCanvas("X","Y",1);

  setTDRStyle(0,0,0);

  TLegend *leg = new TLegend(0.50,0.65,0.90,0.85,NULL,"brNDC");
  leg->SetFillColor(10);

  TFile* file = new TFile("vn30GeV_L_A_histos.root");
  hpt_2ndmu_1stbump->GetXaxis()->SetTitle("sub-leading muon p_{T}, GeV");
  hpt_2ndmu_1stbump->GetXaxis()->SetTitleSize(0.05);
  hpt_2ndmu_1stbump->GetXaxis()->SetLabelSize(0.05);
  hpt_2ndmu_1stbump->GetYaxis()->SetTitle("N events");
  hpt_2ndmu_1stbump->GetYaxis()->SetTitleSize(0.05); 
  hpt_2ndmu_1stbump->GetYaxis()->SetLabelSize(0.05);
  hpt_2ndmu_1stbump->SetMaximum(50.);
  hpt_2ndmu_1stbump->SetLineWidth(4);
  hpt_2ndmu_1stbump->SetLineStyle(1);
  hpt_2ndmu_1stbump->Draw("hist");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{L}, (V-B-q)_{A}","L");

  TFile* file = new TFile("vn30GeV_R_A_histos.root");
  hpt_2ndmu_1stbump->SetLineWidth(2);
  hpt_2ndmu_1stbump->SetLineStyle(2);
  hpt_2ndmu_1stbump->Draw("same");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{R}, (V-B-q)_{A}","L");

  TFile* file = new TFile("vn30GeV_A_A_histos.root");
  hpt_2ndmu_1stbump->SetLineWidth(2);
  hpt_2ndmu_1stbump->SetLineStyle(3);
  hpt_2ndmu_1stbump->Draw("same");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{A}, (V-B-q)_{A}","L");

  TFile* file = new TFile("vn30GeV_V_A_histos.root");
  hpt_2ndmu_1stbump->SetLineWidth(1);
  hpt_2ndmu_1stbump->SetLineStyle(1);
  hpt_2ndmu_1stbump->Draw("same");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{V}, (V-B-q)_{A}","L");

  leg->Draw();
  TLatex *t = new TLatex();
  t->SetTextSize(0.040);
  t->DrawLatex(20,47.0,"1st di-muon excess, VLQ models");
  //
  TLatex *tex = new TLatex(0.15,0.96,"CMS di-muon excess, L=19.8 fb^{-1} (8 TeV)");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.8,0.96,"Simulation");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  c1->SaveAs("pt_2ndmuon_VLQ_1.pdf"); 


  TCanvas* c2 = new TCanvas("X","Y",1);

  setTDRStyle(0,0,0);

  TLegend *leg = new TLegend(0.50,0.65,0.90,0.85,NULL,"brNDC");
  leg->SetFillColor(10);

  TFile* file = new TFile("vn30GeV_L_V_histos.root");
  hpt_2ndmu_1stbump->GetXaxis()->SetTitle("sub-leading muon p_{T}, GeV");
  hpt_2ndmu_1stbump->GetXaxis()->SetTitleSize(0.05);
  hpt_2ndmu_1stbump->GetXaxis()->SetLabelSize(0.05);
  hpt_2ndmu_1stbump->GetYaxis()->SetTitle("N events");
  hpt_2ndmu_1stbump->GetYaxis()->SetTitleSize(0.05); 
  hpt_2ndmu_1stbump->GetYaxis()->SetLabelSize(0.05);
  hpt_2ndmu_1stbump->SetMaximum(50.);
  hpt_2ndmu_1stbump->SetLineWidth(4);
  hpt_2ndmu_1stbump->SetLineStyle(1);
  hpt_2ndmu_1stbump->Draw("hist");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{L}, (V-B-q)_{V}","L");

  TFile* file = new TFile("vn30GeV_R_V_histos.root");
  hpt_2ndmu_1stbump->SetLineWidth(2);
  hpt_2ndmu_1stbump->SetLineStyle(2);
  hpt_2ndmu_1stbump->Draw("same");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{R}, (V-B-q)_{V}","L");

  TFile* file = new TFile("vn30GeV_A_V_histos.root");
  hpt_2ndmu_1stbump->SetLineWidth(2);
  hpt_2ndmu_1stbump->SetLineStyle(3);
  hpt_2ndmu_1stbump->Draw("same");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{A}, (V-B-q)_{V}","L");

  TFile* file = new TFile("vn30GeV_V_V_histos.root");
  hpt_2ndmu_1stbump->SetLineWidth(1);
  hpt_2ndmu_1stbump->SetLineStyle(1);
  hpt_2ndmu_1stbump->Draw("same");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{V}, (V-B-q)_{V}","L");

  leg->Draw();
  TLatex *t = new TLatex();
  t->SetTextSize(0.040);
  t->DrawLatex(20,47.0,"1st di-muon excess, VLQ models");
  //
  TLatex *tex = new TLatex(0.15,0.96,"CMS di-muon excess, L=19.8 fb^{-1} (8 TeV)");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.8,0.96,"Simulation");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  c2->SaveAs("pt_2ndmuon_VLQ_2.pdf"); 

  TCanvas* c3 = new TCanvas("X","Y",1);

  setTDRStyle(0,0,0);

  TLegend *leg = new TLegend(0.50,0.65,0.90,0.85,NULL,"brNDC");
  leg->SetFillColor(10);

  TFile* file = new TFile("vn30GeV_L_L_histos.root");
  hpt_2ndmu_1stbump->GetXaxis()->SetTitle("sub-leading muon p_{T}, GeV");
  hpt_2ndmu_1stbump->GetXaxis()->SetTitleSize(0.05);
  hpt_2ndmu_1stbump->GetXaxis()->SetLabelSize(0.05);
  hpt_2ndmu_1stbump->GetYaxis()->SetTitle("N events");
  hpt_2ndmu_1stbump->GetYaxis()->SetTitleSize(0.05); 
  hpt_2ndmu_1stbump->GetYaxis()->SetLabelSize(0.05);
  hpt_2ndmu_1stbump->SetMaximum(50.);
  hpt_2ndmu_1stbump->SetLineWidth(4);
  hpt_2ndmu_1stbump->SetLineStyle(1);
  hpt_2ndmu_1stbump->Draw("hist");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{L}, (V-B-q)_{L}","L");

  TFile* file = new TFile("vn30GeV_R_L_histos.root");
  hpt_2ndmu_1stbump->SetLineWidth(2);
  hpt_2ndmu_1stbump->SetLineStyle(2);
  hpt_2ndmu_1stbump->Draw("same");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{R}, (V-B-q)_{L}","L");

  TFile* file = new TFile("vn30GeV_A_L_histos.root");
  hpt_2ndmu_1stbump->SetLineWidth(2);
  hpt_2ndmu_1stbump->SetLineStyle(3);
  hpt_2ndmu_1stbump->Draw("same");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{A}, (V-B-q)_{L}","L");

  TFile* file = new TFile("vn30GeV_V_L_histos.root");
  hpt_2ndmu_1stbump->SetLineWidth(1);
  hpt_2ndmu_1stbump->SetLineStyle(1);
  hpt_2ndmu_1stbump->Draw("same");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{V}, (V-B-q)_{L}","L");

  leg->Draw();
  TLatex *t = new TLatex();
  t->SetTextSize(0.040);
  t->DrawLatex(20,47.0,"1st di-muon excess, VLQ models");
  //
  TLatex *tex = new TLatex(0.15,0.96,"CMS di-muon excess, L=19.8 fb^{-1} (8 TeV)");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.8,0.96,"Simulation");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  c3->SaveAs("pt_2ndmuon_VLQ_3.pdf"); 

  TCanvas* c4 = new TCanvas("X","Y",1);

  setTDRStyle(0,0,0);

  TLegend *leg = new TLegend(0.50,0.65,0.90,0.85,NULL,"brNDC");
  leg->SetFillColor(10);

  TFile* file = new TFile("vn30GeV_L_R_histos.root");
  hpt_2ndmu_1stbump->GetXaxis()->SetTitle("sub-leading muon p_{T}, GeV");
  hpt_2ndmu_1stbump->GetXaxis()->SetTitleSize(0.05);
  hpt_2ndmu_1stbump->GetXaxis()->SetLabelSize(0.05);
  hpt_2ndmu_1stbump->GetYaxis()->SetTitle("N events");
  hpt_2ndmu_1stbump->GetYaxis()->SetTitleSize(0.05); 
  hpt_2ndmu_1stbump->GetYaxis()->SetLabelSize(0.05);
  hpt_2ndmu_1stbump->SetMaximum(50.);
  hpt_2ndmu_1stbump->SetLineWidth(4);
  hpt_2ndmu_1stbump->SetLineStyle(1);
  hpt_2ndmu_1stbump->Draw("hist");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{L}, (V-B-q)_{R}","L");

  TFile* file = new TFile("vn30GeV_R_R_histos.root");
  hpt_2ndmu_1stbump->SetLineWidth(2);
  hpt_2ndmu_1stbump->SetLineStyle(2);
  hpt_2ndmu_1stbump->Draw("same");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{R}, (V-B-q)_{R}","L");

  TFile* file = new TFile("vn30GeV_A_R_histos.root");
  hpt_2ndmu_1stbump->SetLineWidth(2);
  hpt_2ndmu_1stbump->SetLineStyle(3);
  hpt_2ndmu_1stbump->Draw("same");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{A}, (V-B-q)_{R}","L");

  TFile* file = new TFile("vn30GeV_V_R_histos.root");
  hpt_2ndmu_1stbump->SetLineWidth(1);
  hpt_2ndmu_1stbump->SetLineStyle(1);
  hpt_2ndmu_1stbump->Draw("same");
  leg->AddEntry(hpt_2ndmu_1stbump,"(V-f)_{V}, (V-B-q)_{R}","L");

  leg->Draw();
  TLatex *t = new TLatex();
  t->SetTextSize(0.040);
  t->DrawLatex(20,47.0,"1st di-muon excess, VLQ models");
  //
  TLatex *tex = new TLatex(0.15,0.96,"CMS di-muon excess, L=19.8 fb^{-1} (8 TeV)");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.8,0.96,"Simulation");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  c4->SaveAs("pt_2ndmuon_VLQ_4.pdf"); 
}
