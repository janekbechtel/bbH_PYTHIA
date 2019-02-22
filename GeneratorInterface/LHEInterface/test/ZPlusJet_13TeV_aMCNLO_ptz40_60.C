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


  // di-muons
  TCanvas* c1 = new TCanvas("X","Y",1);
  TFile* file = new TFile("mumuj_mg5_amcnlo_histos_ptz40_60.root");
  Double_t scale_bar = hptj_div_ptz_barrel->Integral();
  cout <<"  Nev barrel = " << scale_bar << endl;
  hptj_div_ptz_barrel->Scale(1./scale_bar);
  //  hptj_div_ptz_barrel->Rebin();
  hptj_div_ptz_barrel->GetXaxis()->SetTitle("p_{T}^{j1}/p_{T}^{Z}");
  hptj_div_ptz_barrel->GetXaxis()->SetTitleSize(0.05);
  hptj_div_ptz_barrel->GetXaxis()->SetLabelSize(0.05);
  hptj_div_ptz_barrel->GetYaxis()->SetTitle("probability"); 
  hptj_div_ptz_barrel->GetYaxis()->SetTitleSize(0.05); 
  hptj_div_ptz_barrel->GetYaxis()->SetLabelSize(0.05);
  hptj_div_ptz_barrel->GetXaxis()->SetTitleOffset(1.2);
  hptj_div_ptz_barrel->GetYaxis()->SetTitleOffset(1.3);
  hptj_div_ptz_barrel->SetMaximum(0.3);
  hptj_div_ptz_barrel->SetMinimum(0.0);
  hptj_div_ptz_barrel->SetLineStyle(1);
  hptj_div_ptz_barrel->SetLineWidth(3);
  hptj_div_ptz_barrel->Draw("hist");
  TLegend *leg = new TLegend(0.70,0.4,0.90,0.5,NULL,"brNDC");
  leg->SetFillColor(10);
  leg->AddEntry(hptj_div_ptz_barrel,"|#eta_{j1}|<2.5","L");
  
  Double_t scale_hf = hptj_div_ptz_forward->Integral();
  cout <<"  Nev forward = " << scale_hf << endl; 
  hptj_div_ptz_forward->Scale(1./scale_hf);
  //  hptj_div_ptz_forward->Rebin();
  hptj_div_ptz_forward->SetLineStyle(2);
  hptj_div_ptz_forward->SetLineWidth(3);
  hptj_div_ptz_forward->Draw("samehist");
  leg->AddEntry(hptj_div_ptz_forward,"|#eta_{j1}|>3.2","L");
  leg->Draw();

  /*
  hptj_div_ptz_forwardp->Scale(1./scale_hf);
  hptj_div_ptz_forwardm->Scale(1./scale_hf);

  hptj_div_ptz_forwardp->SetLineStyle(1);
  hptj_div_ptz_forwardp->SetLineWidth(1);
  hptj_div_ptz_forwardp->Draw("samehist");
  hptj_div_ptz_forwardm->SetLineStyle(3);
  hptj_div_ptz_forwardm->SetLineWidth(1);
  hptj_div_ptz_forwardm->Draw("samehist");
  TLegend *leg1 = new TLegend(0.15,0.20,0.35,0.30,NULL,"brNDC");
  leg1->SetFillColor(10);
  leg1->AddEntry(hptj_div_ptz_forwardp,"#eta_{J1}>3.2","L");
  leg1->AddEntry(hptj_div_ptz_forwardm,"#eta_{J1}<-3.2","L");
  leg1->Draw();
  */

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
  //  TLatex *tex = new TLatex(0.63,0.96,"36 fb^{-1} re-reco (13 TeV)");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  tex->SetLineWidth(2);
  tex->Draw();
  // end Latex

  TLatex *t = new TLatex();
  t->SetTextSize(0.05);
  t->DrawLatex(0.3,0.28,"Z(#mu#mu)+1 jet NLO, MG5_aMC@NLO+PY8");
  t->SetTextSize(0.04);
  t->DrawLatex(0.15,0.25,"#bullet p_{T}^{#mu1,2}>25, 20 GeV, |#eta^{#mu}|<2.1, 80<m_{#mu#mu}<100 GeV");
  t->DrawLatex(0.15,0.22,"#bullet 40<p_{T}^{Z}<60 GeV, #Delta#phi(Z-j1)>2.6");
  t->DrawLatex(0.15,0.19,"#bullet p_{T}^{j}>20 GeV, |#eta|<4.7");
  t->SetTextSize(0.04);
  t->DrawLatex(1.35,0.05,"j1 - leading p_{T} jet");
  t->SetTextSize(0.05);
  t->DrawLatex(1.0,0.17,"#it{Particle level analysis }");
  c1->SaveAs("ZplusJet_8TeV_aMCNLOptz40_60.pdf"); 
}
