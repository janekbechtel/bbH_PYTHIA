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


  TFile* file = new TFile("CPtestMG5_H_hist.root");
  TH1F* hpth_H = (TH1F*)hpth->Clone();
  TH1F* hmjj_H = (TH1F*)hmjj->Clone();
  TH1F* hdetajj_H = (TH1F*)hdetajj->Clone();
  // after mjj > 300 GeV cut
  TH1F* hdetajj_mjjcut_H = (TH1F*)hdetajj_mjjcut->Clone();
  Double_t scale = 1./hdetajj_mjjcut_H->Integral();
  hdetajj_mjjcut_H->Scale(scale);
  // after deta_jj > 3.0 cut
  TH1F* hpth_detacut_H = (TH1F*)hpth_detacut->Clone();
  TH1F* hdphijj_detacut_H = (TH1F*)hdphijj->Clone();
  TH1F* hmjj_detacut_H = (TH1F*)hmjj_detacut->Clone();
  Double_t scaleH = 1./ hdetajj_H->Integral();
  hpth_H->Scale(scaleH);
  hdetajj_H->Scale(scaleH);
  hmjj_H->Scale(scaleH);
  Double_t scaleH_detacut = 1./ hdphijj_detacut_H->Integral(); 
  hpth_detacut_H->Scale(scaleH_detacut);
  hmjj_detacut_H->Scale(scaleH_detacut);
  hdphijj_detacut_H->Scale(scaleH_detacut);

  TFile* file = new TFile("CPtestMG5_A_hist.root");
  TH1F* hpth_A = (TH1F*)hpth->Clone();
  TH1F* hdetajj_A = (TH1F*)hdetajj->Clone();
  TH1F* hmjj_A = (TH1F*)hmjj->Clone();
  // after mjj > 300 GeV cut
  TH1F* hdetajj_mjjcut_A = (TH1F*)hdetajj_mjjcut->Clone();
  scale = 1./hdetajj_mjjcut_A->Integral();
  hdetajj_mjjcut_A->Scale(scale);
  // after deta_jj > 3.0 cut
  TH1F* hpth_detacut_A = (TH1F*)hpth_detacut->Clone();
  TH1F* hdphijj_detacut_A = (TH1F*)hdphijj->Clone();
  TH1F* hmjj_detacut_A = (TH1F*)hmjj_detacut->Clone();
  Double_t scaleA = 1./ hdetajj_A->Integral();
  hpth_A->Scale(scaleA);
  hdetajj_A->Scale(scaleA);
  hmjj_A->Scale(scaleA);

  Double_t scaleA_detacut = 1./ hdphijj_detacut_A->Integral(); 
  hpth_detacut_A->Scale(scaleA_detacut);
  hmjj_detacut_A->Scale(scaleA_detacut);
  hdphijj_detacut_A->Scale(scaleA_detacut);

  setTDRStyle(0,0,0);

  TCanvas* c1 = new TCanvas("X","Y",1);

  hpth_H->GetXaxis()->SetTitle("Higgs boson p_{T} (GeV)");
  hpth_H->GetYaxis()->SetTitle("probability");
  hpth_H->GetXaxis()->SetTitleSize(0.05);
  hpth_H->GetXaxis()->SetLabelSize(0.05);
  hpth_H->GetYaxis()->SetTitleSize(0.05); 
  hpth_H->GetYaxis()->SetLabelSize(0.05);
  hpth_H->GetYaxis()->SetTitleOffset(1.35);
  hpth_H->SetMaximum(0.06);
  hpth_H->SetLineWidth(2);
  hpth_H->SetLineStyle(1);
  hpth_H->Draw("hist");

  hpth_A->SetLineWidth(2);
  hpth_A->SetLineStyle(2);
  hpth_A->Draw("samehist");

  TLegend *leg = new TLegend(0.60,0.55,0.80,0.65,NULL,"brNDC");
  leg->SetFillColor(10);
  leg->AddEntry(hpth_H,"CP-even h","L");
  leg->AddEntry(hpth_A,"CP-odd h","L");
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetTextSize(0.040);
  t->DrawLatex(30,0.054,"HC MG5_aMC@NLO, gg#rightarrowh(125)+2jets");
  t->DrawLatex(130,0.046,"p_{T}^{jet}>25 GeV, |#eta^{jet}|<4.7");
  //
  TLatex *tex = new TLatex(0.15,0.96,"#bf{CMS}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.65,0.96,"Simulation (13 TeV)");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  c1->SaveAs("CPpth.pdf"); 

  TCanvas* c2 = new TCanvas("X","Y",1);
  hdetajj_H->GetXaxis()->SetTitle("#Delta#eta(j1,j2)");
  hdetajj_H->GetXaxis()->SetTitleSize(0.05);
  hdetajj_H->GetXaxis()->SetLabelSize(0.05);
  hdetajj_H->GetYaxis()->SetTitle("probability");
  hdetajj_H->GetYaxis()->SetTitleOffset(1.35);
  hdetajj_H->GetYaxis()->SetTitleSize(0.05); 
  hdetajj_H->GetYaxis()->SetLabelSize(0.05);
  hdetajj_H->SetMaximum(0.1);
  //  hdetajj_H->SetMinimum(1000.);
  hdetajj_H->SetLineWidth(2);
  hdetajj_H->SetLineStyle(1);
  hdetajj_H->Draw("hist");

  hdetajj_A->SetLineWidth(2);
  hdetajj_A->SetLineStyle(2);
  hdetajj_A->Draw("samehist");

  TLegend *leg = new TLegend(0.60,0.55,0.80,0.65,NULL,"brNDC");
  leg->SetFillColor(10);
  leg->AddEntry(hdetajj_H,"CP-even h","L");
  leg->AddEntry(hdetajj_A,"CP-odd h","L");
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetTextSize(0.040);
  t->DrawLatex(1.0,0.09,"HC MG5_aMC@NLO, gg#rightarrowh(125)+2jets");
  t->DrawLatex(3.0,0.08,"p_{T}^{jet}>25 GeV, |#eta^{jet}|<4.7");
  //
  TLatex *tex = new TLatex(0.15,0.96,"#bf{CMS}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.65,0.96,"Simulation (13 TeV)");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  c2->SaveAs("CPdetajj.pdf"); 

  TCanvas* c3 = new TCanvas("X","Y",1);
  hmjj_H->GetXaxis()->SetTitle("m(j1,j2) (GeV)");
  hmjj_H->GetXaxis()->SetTitleSize(0.05);
  hmjj_H->GetXaxis()->SetLabelSize(0.05);
  hmjj_H->GetYaxis()->SetTitle("probability");
  hmjj_H->GetYaxis()->SetTitleOffset(1.35);
  hmjj_H->GetYaxis()->SetTitleSize(0.05); 
  hmjj_H->GetYaxis()->SetLabelSize(0.05);
  hmjj_H->SetMaximum(0.07);
  //  hmjj_H->SetMinimum(1000.);
  hmjj_H->SetLineWidth(2);
  hmjj_H->SetLineStyle(1);
  hmjj_H->Draw("hist");

  hmjj_A->SetLineWidth(2);
  hmjj_A->SetLineStyle(2);
  hmjj_A->Draw("samehist");

  TLegend *leg = new TLegend(0.60,0.55,0.80,0.65,NULL,"brNDC");
  leg->SetFillColor(10);
  leg->AddEntry(hmjj_H,"CP-even h","L");
  leg->AddEntry(hmjj_A,"CP-odd h","L");
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetTextSize(0.040);
  t->DrawLatex(140.0,0.064,"HC MG5_aMC@NLO, gg#rightarrowh(125)+2jets");
  t->DrawLatex(450.0,0.058,"p_{T}^{jet}>25 GeV, |#eta^{jet}|<4.7");
  //
  TLatex *tex = new TLatex(0.15,0.96,"#bf{CMS}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.65,0.96,"Simulation (13 TeV)");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  c3->SaveAs("CPmjj.pdf"); 


  TCanvas* c4 = new TCanvas("X","Y",1);
  hmjj_detacut_H->GetXaxis()->SetTitle("m(j1,j2) (GeV)");
  hmjj_detacut_H->GetXaxis()->SetTitleSize(0.05);
  hmjj_detacut_H->GetXaxis()->SetLabelSize(0.05);
  hmjj_detacut_H->GetYaxis()->SetTitle("probability");
  hmjj_detacut_H->GetYaxis()->SetTitleOffset(1.35);
  hmjj_detacut_H->GetYaxis()->SetTitleSize(0.05); 
  hmjj_detacut_H->GetYaxis()->SetLabelSize(0.05);
  hmjj_detacut_H->SetMaximum(0.20);
  //  hmjj_detacut_H->SetMinimum(1000.);
  hmjj_detacut_H->SetLineWidth(2);
  hmjj_detacut_H->SetLineStyle(1);
  hmjj_detacut_H->Draw("hist");

  hmjj_detacut_A->SetLineWidth(2);
  hmjj_detacut_A->SetLineStyle(2);
  hmjj_detacut_A->Draw("samehist");

  TLegend *leg = new TLegend(0.60,0.55,0.80,0.65,NULL,"brNDC");
  leg->SetFillColor(10);
  leg->AddEntry(hmjj_detacut_H,"CP-even h","L");
  leg->AddEntry(hmjj_detacut_A,"CP-odd h","L");
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetTextSize(0.040);
  t->DrawLatex(140.0,0.18,"HC MG5_aMC@NLO, gg#rightarrowh(125)+2jets");
  t->DrawLatex(450.0,0.16,"p_{T}^{jet}>25 GeV, |#eta^{jet}|<4.7");
  t->DrawLatex(450.0,0.14,"|#Delta#eta_{j1j2}|>3.0");
  //
  TLatex *tex = new TLatex(0.15,0.96,"#bf{CMS}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.65,0.96,"Simulation (13 TeV)");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  c4->SaveAs("CPmjj_detacut.pdf"); 

  TCanvas* c5 = new TCanvas("X","Y",1);
  hpth_detacut_H->GetXaxis()->SetTitle("Higgs boson p_{T} (GeV)");
  hpth_detacut_H->GetXaxis()->SetTitleSize(0.05);
  hpth_detacut_H->GetXaxis()->SetLabelSize(0.05);
  hpth_detacut_H->GetYaxis()->SetTitle("probability");
  hpth_detacut_H->GetYaxis()->SetTitleOffset(1.35);
  hpth_detacut_H->GetYaxis()->SetTitleSize(0.05); 
  hpth_detacut_H->GetYaxis()->SetLabelSize(0.05);
  hpth_detacut_H->SetMaximum(0.20);
  //  hpth_detacut_H->SetMinimum(1000.);
  hpth_detacut_H->SetLineWidth(2);
  hpth_detacut_H->SetLineStyle(1);
  hpth_detacut_H->Draw("hist");

  hpth_detacut_A->SetLineWidth(2);
  hpth_detacut_A->SetLineStyle(2);
  hpth_detacut_A->Draw("samehist");

  TLegend *leg = new TLegend(0.60,0.50,0.80,0.60,NULL,"brNDC");
  leg->SetFillColor(10);
  leg->AddEntry(hpth_detacut_H,"CP-even h","L");
  leg->AddEntry(hpth_detacut_A,"CP-odd h","L");
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetTextSize(0.040);
  t->DrawLatex(30.0,0.18,"HC MG5_aMC@NLO, gg#rightarrowh(125)+2jets");
  t->DrawLatex(100.0,0.16,"p_{T}^{jet}>25 GeV, |#eta^{jet}|<4.7");
  t->DrawLatex(100.0,0.14,"|#Delta#eta_{j1j2}|>3.0");
  //
  TLatex *tex = new TLatex(0.15,0.96,"#bf{CMS}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.65,0.96,"Simulation (13 TeV)");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  c5->SaveAs("CPpth_detacut.pdf"); 

  TCanvas* c6 = new TCanvas("X","Y",1);
  hdetajj_mjjcut_H->GetXaxis()->SetTitle("#Delta#eta(j1,j2)");
  hdetajj_mjjcut_H->GetXaxis()->SetTitleSize(0.05);
  hdetajj_mjjcut_H->GetXaxis()->SetLabelSize(0.05);
  hdetajj_mjjcut_H->GetYaxis()->SetTitle("probability");
  hdetajj_mjjcut_H->GetYaxis()->SetTitleOffset(1.35);
  hdetajj_mjjcut_H->GetYaxis()->SetTitleSize(0.05); 
  hdetajj_mjjcut_H->GetYaxis()->SetLabelSize(0.05);
  hdetajj_mjjcut_H->SetMaximum(0.1);
  //  hdetajj_mjjcut_H->SetMinimum(1000.);
  hdetajj_mjjcut_H->SetLineWidth(2);
  hdetajj_mjjcut_H->SetLineStyle(1);
  hdetajj_mjjcut_H->Draw("hist");

  hdetajj_mjjcut_A->SetLineWidth(2);
  hdetajj_mjjcut_A->SetLineStyle(2);
  hdetajj_mjjcut_A->Draw("samehist");

  TLegend *leg = new TLegend(0.70,0.55,0.90,0.65,NULL,"brNDC");
  leg->SetFillColor(10);
  leg->AddEntry(hdetajj_mjjcut_H,"CP-even h","L");
  leg->AddEntry(hdetajj_mjjcut_A,"CP-odd h","L");
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetTextSize(0.040);
  t->DrawLatex(1.0,0.09,"HC MG5_aMC@NLO, gg#rightarrowh(125)+2jets");
  t->DrawLatex(3.0,0.08,"p_{T}^{jet}>25 GeV, |#eta^{jet}|<4.7");
  t->DrawLatex(3.0,0.07,"m_{j1j2}>300 GeV");
  //
  TLatex *tex = new TLatex(0.15,0.96,"#bf{CMS}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.65,0.96,"Simulation (13 TeV)");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  c6->SaveAs("CPdetajj_mjjcut.pdf"); 

  TCanvas* c7 = new TCanvas("X","Y",1);
  hdphijj_detacut_H->GetXaxis()->SetTitle("#Delta#phi(j1,j2) (degree)");
  hdphijj_detacut_H->GetXaxis()->SetTitleSize(0.05);
  hdphijj_detacut_H->GetXaxis()->SetLabelSize(0.05);
  hdphijj_detacut_H->GetYaxis()->SetTitle("probability");
  hdphijj_detacut_H->GetYaxis()->SetTitleOffset(1.35);
  hdphijj_detacut_H->GetYaxis()->SetTitleSize(0.05); 
  hdphijj_detacut_H->GetYaxis()->SetLabelSize(0.05);
  hdphijj_detacut_H->SetMaximum(0.14);
  hdphijj_detacut_H->SetMinimum(0.00);
  hdphijj_detacut_H->SetLineWidth(2);
  hdphijj_detacut_H->SetLineStyle(1);
  hdphijj_detacut_H->Draw("hist");

  hdphijj_detacut_A->SetLineWidth(2);
  hdphijj_detacut_A->SetLineStyle(2);
  hdphijj_detacut_A->Draw("samehist");

  TLegend *leg = new TLegend(0.70,0.65,0.90,0.75,NULL,"brNDC");
  leg->SetFillColor(10);
  leg->AddEntry(hdphijj_detacut_H,"CP-even h","L");
  leg->AddEntry(hdphijj_detacut_A,"CP-odd h","L");
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetTextSize(0.040);
  t->DrawLatex(-100.,0.13,"HC MG5_aMC@NLO, gg#rightarrowh(125)+2jets");
  t->DrawLatex(-50.,0.12,"p_{T}^{jet}>25 GeV, |#eta^{jet}|<4.7");
  t->DrawLatex(-50.,0.11,"|#Delta#eta_{j1j2}|>3.0");
  //
  TLatex *tex = new TLatex(0.15,0.96,"#bf{CMS}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.65,0.96,"Simulation (13 TeV)");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  c7->SaveAs("CPdphijj_detacut.pdf"); 
}
