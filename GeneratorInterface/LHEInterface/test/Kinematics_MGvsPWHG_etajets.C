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

  setTDRStyle(0,0,0);

  TCanvas* c1 = new TCanvas("X","Y",1);
  TFile* file = new TFile("bbh_tautau_700gev_powheg_hist.root");
  Double_t scale = 1./hetab1->Integral();
  hetab1->Scale(scale);
  hetab1->GetXaxis()->SetTitle("leading b-jet #eta");
  hetab1->GetXaxis()->SetTitleSize(0.05);
  hetab1->GetXaxis()->SetLabelSize(0.05);
  hetab1->GetYaxis()->SetTitle("probability"); 
  hetab1->GetYaxis()->SetTitleSize(0.05); 
  hetab1->GetYaxis()->SetLabelSize(0.05);
  hetab1->GetXaxis()->SetTitleOffset(1.2);
  hetab1->GetYaxis()->SetTitleOffset(1.2);
  hetab1->SetMaximum(0.15);
  hetab1->SetMinimum(0.0);
  hetab1->SetLineStyle(1);
  hetab1->SetLineWidth(3);
  hetab1->Draw("hist");
  TLegend *leg = new TLegend(0.50,0.70,0.90,0.80,NULL,"brNDC");
  leg->SetFillColor(10);
  leg->AddEntry(hetab1,"POWHEG+PY8","L");
  TFile* file = new TFile("bbh_tautau_700gev_mg5_hist.root");
  scale = 1./hetab1->Integral();
  hetab1->Scale(scale);
  hetab1->SetLineStyle(2);
  hetab1->SetLineWidth(3);
  hetab1->Draw("samehist");
  leg->AddEntry(hetab1,"MG5_aMC@NLO+PY8","L");
  leg->Draw();
  //  Latex
  TLatex *tex = new TLatex(0.15,0.96,"#bf{CMS}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  TLatex *tex = new TLatex(0.23,0.96,"#it{Simulation}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(18);
  tex->SetLineWidth(2);
  tex->Draw();
  TLatex *t = new TLatex();
  t->SetTextSize(0.05);
  t->DrawLatex(-3.0,0.135,"gg#rightarrowbbA, m_{A}=700 GeV");
  c1->SaveAs("etab1_mA700_MGvsPWHG.pdf"); 

  TCanvas* c2 = new TCanvas("X","Y",1);
  TFile* file = new TFile("bbh_tautau_700gev_powheg_hist.root");
  Double_t scale = 1./hetab2->Integral();
  hetab2->Scale(scale);
  hetab2->GetXaxis()->SetTitle("sub-leading b-jet #eta");
  hetab2->GetXaxis()->SetTitleSize(0.05);
  hetab2->GetXaxis()->SetLabelSize(0.05);
  hetab2->GetYaxis()->SetTitle("probability"); 
  hetab2->GetYaxis()->SetTitleSize(0.05); 
  hetab2->GetYaxis()->SetLabelSize(0.05);
  hetab2->GetXaxis()->SetTitleOffset(1.2);
  hetab2->GetYaxis()->SetTitleOffset(1.2);
  hetab2->SetMaximum(0.17);
  hetab2->SetMinimum(0.0);
  hetab2->SetLineStyle(1);
  hetab2->SetLineWidth(3);
  hetab2->Draw("hist");
  TLegend *leg = new TLegend(0.50,0.70,0.90,0.80,NULL,"brNDC");
  leg->SetFillColor(10);
  leg->AddEntry(hetab2,"POWHEG+PY8","L");
  TFile* file = new TFile("bbh_tautau_700gev_mg5_hist.root");
  scale = 1./hetab2->Integral();
  hetab2->Scale(scale);
  hetab2->SetLineStyle(2);
  hetab2->SetLineWidth(3);
  hetab2->Draw("samehist");
  leg->AddEntry(hetab2,"MG5_aMC@NLO+PY8","L");
  leg->Draw();
  //  Latex
  TLatex *tex = new TLatex(0.15,0.96,"#bf{CMS}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  TLatex *tex = new TLatex(0.23,0.96,"#it{Simulation}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(18);
  tex->SetLineWidth(2);
  tex->Draw();
  TLatex *t = new TLatex();
  t->SetTextSize(0.05);
  t->DrawLatex(-3.0,0.15,"gg#rightarrowbbA, m_{A}=700 GeV");
  c2->SaveAs("etab2_mA700_MGvsPWHG.pdf"); 

  TCanvas* c3 = new TCanvas("X","Y",1);
  TFile* file = new TFile("bbh_tautau_700gev_powheg_hist.root");
  Double_t scale = 1./hetaj->Integral();
  hetaj->Scale(scale);
  hetaj->GetXaxis()->SetTitle("leading non-b-jet #eta");
  hetaj->GetXaxis()->SetTitleSize(0.05);
  hetaj->GetXaxis()->SetLabelSize(0.05);
  hetaj->GetYaxis()->SetTitle("probability"); 
  hetaj->GetYaxis()->SetTitleSize(0.05); 
  hetaj->GetYaxis()->SetLabelSize(0.05);
  hetaj->GetXaxis()->SetTitleOffset(1.2);
  hetaj->GetYaxis()->SetTitleOffset(1.2);
  hetaj->SetMaximum(0.14);
  hetaj->SetMinimum(0.0);
  hetaj->SetLineStyle(1);
  hetaj->SetLineWidth(3);
  hetaj->Draw("hist");
  TLegend *leg = new TLegend(0.50,0.70,0.90,0.80,NULL,"brNDC");
  leg->SetFillColor(10);
  leg->AddEntry(hetaj,"POWHEG+PY8","L");
  TFile* file = new TFile("bbh_tautau_700gev_mg5_hist.root");
  scale = 1./hetaj->Integral();
  hetaj->Scale(scale);
  hetaj->SetLineStyle(2);
  hetaj->SetLineWidth(3);
  hetaj->Draw("samehist");
  leg->AddEntry(hetaj,"MG5_aMC@NLO+PY8","L");
  leg->Draw();
  //  Latex
  TLatex *tex = new TLatex(0.15,0.96,"#bf{CMS}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  TLatex *tex = new TLatex(0.23,0.96,"#it{Simulation}");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(18);
  tex->SetLineWidth(2);
  tex->Draw();
  TLatex *t = new TLatex();
  t->SetTextSize(0.05);
  t->DrawLatex(-4.0,0.12,"gg#rightarrowbbA, m_{A}=700 GeV; #geq 1b selection");
  c3->SaveAs("etaj_mA700_MGvsPWHG.pdf"); 
}
