#define single_top_cxx
#include "single_top.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

void single_top::Weights()
{
  std::cout <<"  I am in method Weights " << std::endl;

  TFile* file = new TFile("costhetaw_for_weights.root");
  TH1F* hcosthetaw = (TH1F*)file->Get("hcosthetaw");
  TH1F* hcosthetaw_for_weights = (TH1F*)hcosthetaw->Clone();
  nbins = hcosthetaw_for_weights->GetNbinsX();
  Double_t hcosthetaw_for_weights_nev = hcosthetaw_for_weights->Integral();
  Double_t nev_per_bin = hcosthetaw_for_weights_nev/nbins;

  // histogram with flat distribution
  TH1F* hcosthetas_for_weights  = new TH1F("hcosthetas_for_weights", "costhetas_for_weights", nbins, -1.0, 1.0);
  // fill it uniformly
  for (Int_t i = 1; i <= nbins; i++) {
    hcosthetas_for_weights->SetBinContent(i,nev_per_bin);
  }

  // ratio of histograms "flat"/"w"
  TH1F* hweighthistos  = new TH1F("hweighthistos", "weighthistos", nbins, -1.0, 1.0);
  hweighthistos->Divide(hcosthetas_for_weights,hcosthetaw_for_weights,1.,1.);
  // evaluate weights
  std::cout <<"  Initial histos for weight calculation " << std::endl;
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t n_hcosthetaw_for_weights = hcosthetaw_for_weights->GetBinContent(i);
    Double_t n_hcosthetas_for_weights = hcosthetas_for_weights->GetBinContent(i);
    Double_t r_hweighthistos = hweighthistos->GetBinContent(i);
    std::cout <<"  bin i = " << i 
	 <<" " << hweighthistos->GetBinLowEdge(i) 
	 <<" < costheta < " << hweighthistos->GetBinLowEdge(i)+hweighthistos->GetBinWidth(i)
	 <<" nw= " << n_hcosthetaw_for_weights 
	 <<" ns= " << n_hcosthetas_for_weights 
	 <<" r= " << r_hweighthistos << std::endl; 
    lowedge[i-1] = hweighthistos->GetBinLowEdge(i);
    width[i-1] = hweighthistos->GetBinWidth(i);
    weights[i-1] = r_hweighthistos; 
  }

  std::cout <<"  My final weights " << std::endl;
  for (Int_t i = 0; i < nbins; i++) {
    std::cout <<"  -> i = " << i 
	 <<" low edge= " << lowedge[i]
	 <<" width= " << width[i]
	 <<" weight= " << weights[i] << std::endl;
  }
}

void single_top::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L single_top.C
//      Root > single_top
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Weights();

   // parton level plots
   TH1F * hpt_b = new TH1F("hpt_b", "hpt_b", 100, 0., 200.);
   TH1F * heta_b = new TH1F("heta_b", "heta_b", 100, -5., 5.);
   TH1F * hpt_b_1stbump = new TH1F("hpt_b_1stbump", "hpt_b_1stbump", 100, 0., 200.);
   TH1F * heta_b_1stbump = new TH1F("heta_b_1stbump", "heta_b_1stbump", 100, -5., 5.);
   TH1F * hpt_b_2ndbump = new TH1F("hpt_b_2ndbump", "hpt_b_2ndbump", 100, 0., 200.);
   TH1F * heta_b_2ndbump = new TH1F("heta_b_2ndbump", "heta_b_2ndbump", 100, -5., 5.);


   TH1F * hpt_q = new TH1F("hpt_q", "hpt_q", 100, 0., 200.);
   TH1F * heta_q = new TH1F("heta_q", "heta_q", 100, -5., 5.);
   TH1F * hpt_q_bcut = new TH1F("hpt_q_bcut", "hpt_q_bcut", 100, 0., 200.);
   TH1F * heta_q_bcut = new TH1F("heta_q_bcut", "heta_q_bcut", 100, -5., 5.);

   TH1F * hpt_q_bcut_1stbump = new TH1F("hpt_q_bcut_1stbump", "hpt_q_bcut_1stbump", 100, 0., 200.);
   TH1F * hpt_q_bcut_2ndbump = new TH1F("hpt_q_bcut_2ndbump", "hpt_q_bcut_2ndbump", 100, 0., 200.);

   TH1F * heta_q_bcut_qptcut = new TH1F("heta_q_bcut_qptcut", "heta_q_bcut_qptcut", 100, -5., 5.);
   
   // cos theta muon in SCM of di-muon system
   TH1F * hptw = new TH1F( "hptw", "ptw", 40, 0., 200.);
   TH1F * hcosthetaw = new TH1F("hcosthetaw", "costhetaw", 50, -1.0, 1.0);
   TH1F * hcostheta_scalar = new TH1F("hcostheta_scalar", "costheta_scalar", 50, -1.0, 1.0);
   TH1F * hcostheta_longve = new TH1F("hcostheta_longve", "costheta_longve", 50, -1.0, 1.0);
   TH1F * hcostheta_trnpve = new TH1F("hcostheta_trnpve", "costheta_trnpve", 50, -1.0, 1.0);
   TH1F * hcostheta_trnmve = new TH1F("hcostheta_trnmve", "costheta_trnmve", 50, -1.0, 1.0);
   TH1F * hcostheta_s2m1 = new TH1F("hcostheta_s2m1", "costheta_s2m1", 50, -1.0, 1.0);

   TH1F * hdrmub = new TH1F("hdrmub", "drmub", 60, 0, 6);

   // histos for 1st bump
   TH1F * hptw_1stbump = new TH1F( "hptw_1stbump", "ptw_1stbump", 40, 0., 200.);
   TH1F * hpt_1stmu_1stbump = new TH1F("hpt_1stmu_1stbump", "pt_1stmu_1stbump", 150, 0., 150.);

   TH1F * hpt_2ndmu_1stbump = new TH1F("hpt_2ndmu_1stbump", "pt_2ndmu_1stbump", 80, 0., 80.);
   TH1F * hpt_2ndmu_1stbump_scalar = new TH1F("hpt_2ndmu_1stbump_scalar", "pt_2ndmu_1stbump_scalar", 80, 0., 80.);
   TH1F * hpt_2ndmu_1stbump_trnpve = new TH1F("hpt_2ndmu_1stbump_trnpve", "pt_2ndmu_1stbump_trnpve", 80, 0., 80.);
   TH1F * hpt_2ndmu_1stbump_trnmve = new TH1F("hpt_2ndmu_1stbump_trnmve", "pt_2ndmu_1stbump_trnmve", 80, 0., 80.);
   TH1F * hpt_2ndmu_1stbump_s2m1 = new TH1F("hpt_2ndmu_1stbump_s2m1", "pt_2ndmu_1stbump_s2m1", 80, 0., 80.);

   TH1F * hcostheta_1stbump = new TH1F("hcostheta_1stbump", "costheta_1stbump", 10, -1.0, 1.0);
   TH1F * hpt_dimuon_1stbump = new TH1F("hpt_dimuon_1stbump", "pt_dimuon_1stbump", 20, 0., 200.);
   TH1F * hpt_dimuon_1stbump_s2m1 = new TH1F("hpt_dimuon_1stbump_s2m1", "pt_dimuon_1stbump_s2m1", 20, 0., 200.);

   TH1F * hnbj_endcap_1stbump = new TH1F("hnbj_endcap_1stbump", "nbj_endcap_1stbump", 5, 0.0, 5.0);
   TH1F * hnlj_endcap_1stbump = new TH1F("hnlj_endcap_1stbump", "nlj_endcap_1stbump", 5, 0.0, 5.0);
   TH1F * hnj_endcap_1stbump  = new TH1F("hnj_endcap_1stbump", "nj_endcap_1stbump", 5, 0.0, 5.0);

   TH1F * hmass_dijet_1stbump = new TH1F("hmass_dijet_1stbump", "mass_dijet_1stbump", 25, 0., 1000.);
   TH1F * hm_mumuj_1stbump = new TH1F("hm_mumuj_1stbump", "m_mumuj_1stbump", 25, 0., 1000.);
   TH1F * hmass_dimuondijet_1stbump  
     = new TH1F("hmass_dimuondijet_1stbump", "mass_dimuondijet_1stbump", 25, 0., 1000.);
   TH1F * heta_dimuondijet_1stbump = new TH1F("heta_dimuondijet_1stbump", "eta_dimuondijet_1stdbump", 20, -5., 5.);

   TH1F * hpt_bjet_1stbump  = new TH1F("hpt_bjet_1stbump", "pt_bjet_1stbump", 20, 0., 200.);
   TH1F * heta_bjet_1stbump  = new TH1F("heta_bjet_1stbump", "eta_bjet_1stbump", 20, -5.0, 5.0);
   TH1F * hpt_fwdj_1stbump  = new TH1F("hpt_fwdj_1stbump", "pt_fwdj_1stbump", 20, 0., 200.);
   TH1F * heta_fwdj_1stbump  = new TH1F("heta_fwdj_1stbump", "eta_fwdj_1stbump", 20, -5.0, 5.0);

   // histos for 2nd bump
   TH1F * hptw_2ndbump = new TH1F( "hptw_2ndbump", "ptw_2ndbump", 40, 0., 200.);
   TH1F * hpt_1stmu_2ndbump = new TH1F( "hpt_1stmu_2ndbump", "pt_1stmu_2ndbump", 150, 0., 150.);

   TH1F * hpt_2ndmu_2ndbump = new TH1F( "hpt_2ndmu_2ndbump", "pt_2ndmu_2ndbump", 80, 0., 80.);
   TH1F * hpt_2ndmu_2ndbump_scalar = new TH1F("hpt_2ndmu_2ndbump_scalar", "pt_2ndmu_2ndbump_scalar", 80, 0., 80.);
   TH1F * hpt_2ndmu_2ndbump_trnpve = new TH1F("hpt_2ndmu_2ndbump_trnpve", "pt_2ndmu_2ndbump_trnpve", 80, 0., 80.);
   TH1F * hpt_2ndmu_2ndbump_trnmve = new TH1F("hpt_2ndmu_2ndbump_trnmve", "pt_2ndmu_2ndbump_trnmve", 80, 0., 80.);


   TH1F * hcostheta_2ndbump = new TH1F("hcostheta_2ndbump", "costheta_2ndbump", 10, -1.0, 1.0);
   TH1F * hpt_dimuon_2ndbump  = new TH1F("hpt_dimuon_2ndbump", "pt_dimuon_2ndbump", 20, 0., 200.);

   TH1F * hnbj_barrel_2ndbump = new TH1F("hnbj_barrel_2ndbump", "nbj_barrel_2ndbump", 5, 0.0, 5.0);
   TH1F * hnlj_barrel_2ndbump = new TH1F("hnlj_barrel_2ndbump", "nlj_barrel_2ndbump", 5, 0.0, 5.0);
   TH1F * hnj_barrel_2ndbump  = new TH1F("hnj_barrel_2ndbump", "nj_barrel_2ndbump", 5, 0.0, 5.0);

   TH1F * hmass_dijet_2ndbump = new TH1F("hmass_dijet_2ndbump", "mass_dijet_2ndbump", 25, 0., 1000.);
   TH1F * hm_mumuj_2ndbump = new TH1F("hm_mumuj_2ndbump", "m_mumuj_2ndbump", 25, 0., 1000.);
   TH1F * hmass_dimuondijet_2ndbump  
     = new TH1F("hmass_dimuondijet_2ndbump", "mass_dimuondijet_2ndbump", 25, 0., 1000.);
   TH1F * heta_dimuondijet_2ndbump = new TH1F("heta_dimuondijet_2ndbump", "eta_dimuondijet_2ndbump", 20, -5., 5.);

   TH1F * hpt_bjet_2ndbump  = new TH1F("hpt_bjet_2ndbump", "pt_bjet_2ndbump", 20, 0., 200.);
   TH1F * heta_bjet_2ndbump  = new TH1F("heta_bjet_2ndbump", "eta_bjet_2ndbump", 20, -5.0, 5.0);
   TH1F * hpt_barj_2ndbump  = new TH1F("hpt_barj_2ndbump", "pt_barj_2ndbump", 20, 0., 200.);
   TH1F * heta_barj_2ndbump  = new TH1F("heta_barj_2ndbump", "eta_barj_2ndbump", 20, -5.0, 5.0);

   TH1F * hdphi_dimuondijet_2ndbump = new TH1F("hdphi_dimuondijet_2ndbump","dphi_dimuondijet_2ndbump", 17, 0., 3.4);

   //
   Double_t ptmucutmax = 15.;
   Double_t ptmucutmin = 10.;
   Double_t etamucut1  = 2.1;
   Double_t etamucut2  = 2.1;
   Double_t ptj_cut = 30.;
   Double_t etaj_cut = 2.4;
   Double_t dphi_mumu_dijet_2ndbump_cut = 2.5;
   //   Double_t costheta_cut = 0.4;
   Double_t costheta_cut = 1.0;

   Double_t n_total = 0.;
   Double_t w_total = 0.;

   Double_t n_muons = 0.;
   Double_t w_muons = 0.;

   Double_t n_1stbump = 0;
   Double_t w_1stbump = 0;

   Double_t n_2ndbump = 0;
   Double_t w_2ndbump = 0;


   Double_t n_1stbump_q = 0;
   Double_t w_1stbump_q = 0;

   Double_t n_2ndbump_q = 0;
   Double_t w_2ndbump_q = 0;


   Double_t Rmujetmatch = 0.5;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      Double_t ptmu1 = ptmu;
      Double_t etamu1 = etamu;
      Double_t phimu1 = phimu;
      Double_t ptmu2 = ptnu;
      Double_t etamu2 = etanu;
      Double_t phimu2 = phinu;

      Double_t dr_mu_b = deltaR(etamu,phimu,etab,phib);
      hdrmub->Fill(dr_mu_b);
      // prepare for cos calculation in the dimuon center mass
      // 1st muon is mu- (13), second muon is mu+ (-13)

      TLorentzVector mu1st, mu2nd, dimuon;

      mu1st.SetPx(ptmu1*cos(phimu1));
      mu1st.SetPy(ptmu1*sin(phimu1));
      Double_t theta1 = 2. * atan(exp(-etamu1));
      mu1st.SetPz(ptmu1/tan(theta1));
      mu1st.SetE(ptmu1/sin(theta1));

      mu2nd.SetPx(ptmu2*cos(phimu2));
      mu2nd.SetPy(ptmu2*sin(phimu2));
      Double_t theta2 = 2. * atan(exp(-etamu2));
      mu2nd.SetPz(ptmu2/tan(theta2));
      mu2nd.SetE(ptmu2/sin(theta2));
      
      dimuon = mu1st + mu2nd;

      // go to di-muon system frame
      Double_t bx = -dimuon.Px()/dimuon.E();
      Double_t by = -dimuon.Py()/dimuon.E();
      Double_t bz = -dimuon.Pz()/dimuon.E();

      // move mu- to di-muon center of mass
      mu1st.Boost(bx,by,bz); 
      // calculate costheta
      Double_t costheta = (dimuon.Px()* mu1st.Px()+
			   dimuon.Py()* mu1st.Py()+
			   dimuon.Pz()* mu1st.Pz()) / (dimuon.P()*mu1st.P());
      // new weight to change cons theta
      // spin 1
      Double_t lpolar   = (1-costheta*costheta);
      Double_t trpolarp = (1-costheta)*(1-costheta);
      Double_t trpolarm = (1+costheta)*(1+costheta);
      // spin 2 projection 1
      //      Double_t s2m1 = (1-3*costheta*costheta+4*costheta*costheta*costheta*costheta);
      // spin2 projection 2
      //      Double_t s2m1 = (1-costheta*costheta*costheta*costheta);
      Double_t s2m1 = (1-costheta*costheta)*(1-costheta*costheta);

      // weights
      Double_t weight_scalar = 0.;
      for (Int_t i = 0; i < 50; i++) {
	if(costheta >= lowedge[i] && costheta < lowedge[i]+width[i]) {
	  weight_scalar = weights[i];
	} 
      }

      Double_t weight_longve = weight_scalar*lpolar;
      Double_t weight_trnpve = weight_scalar*trpolarp;
      Double_t weight_trnmve = weight_scalar*trpolarm;
      Double_t weight_s2m1   = weight_scalar*s2m1;

      // choose the model
      if(jentry==0) {
	std::cout <<" " << std::endl;
	std::cout <<" Longitudinal Vector Model" << std::endl;
	//	std::cout <<" Scalar Model" << std::endl;
	//	std::cout <<" Transverse Vactor Model. Polarization+" << std::endl;
	//       	std::cout <<" Transverse Vactor Model. Polarization-" << std::endl;
      }
      Double_t event_weight = weight_longve;
      //      Double_t event_weight = weight_scalar;
      //      Double_t event_weight = weight_trnpve 
      //      Double_t event_weight = weight_trnmve; 

      hcostheta_scalar->Fill(costheta,weight_scalar);
      hcostheta_longve->Fill(costheta,weight_longve);
      hcostheta_trnpve->Fill(costheta,weight_trnpve);
      hcostheta_trnmve->Fill(costheta,weight_trnmve);
      hcostheta_s2m1->Fill(costheta,weight_s2m1);


      Double_t ptmumax = ptmu1;
      Double_t etamumax = etamu1;
      Double_t ptmumin = ptmu2;
      Double_t etamumin = etamu2;

      if(ptmu1 < ptmu2 ){
	ptmumax = ptmu2;
	etamumax = etamu2;
	ptmumin = ptmu1;
	etamumin = etamu1;
      }

      n_total += 1;
      w_total += event_weight;
      // selections
      if(ptmumax<ptmucutmax || ptmumin<ptmucutmin || fabs(etamumax)>etamucut1 || fabs(etamumin)>etamucut2) 
	{continue;}

      n_muons += 1;
      w_muons += event_weight;

      hcosthetaw->Fill(costheta,event_weight);
      hptw->Fill(ptw,event_weight);

      // 1st bump parton selections
      if( ptb > ptj_cut && fabs(etab) < etaj_cut && ptq > ptj_cut) { 
	if(fabs(etaq) > etaj_cut) {
	  n_1stbump_q += 1;
	  w_1stbump_q += event_weight;
	} else {
	  n_2ndbump_q += 1;
	  w_2ndbump_q += event_weight;
	}
      }

      // partons
      hpt_b->Fill(ptb);
      heta_b->Fill(etab);
      hpt_q->Fill(ptq);
      heta_q->Fill(etaq);
      
      if(ptq > ptj_cut) {
	if(fabs(etaq) > etaj_cut) {
	  hpt_b_1stbump->Fill(ptb);
	  heta_b_1stbump->Fill(etab);
	} else {
	  hpt_b_2ndbump->Fill(ptb);
	  heta_b_2ndbump->Fill(etab);
	}
      }

      if(ptb > ptj_cut && fabs(etab) < etaj_cut) {
	hpt_q_bcut->Fill(ptq);
	heta_q_bcut->Fill(etaq);
	if(fabs(etaq) > etaj_cut) {
	  hpt_q_bcut_1stbump->Fill(ptq);
	} else {
	  hpt_q_bcut_2ndbump->Fill(ptq);
	}
	if(ptq > ptj_cut) {
	  heta_q_bcut_qptcut->Fill(etaq);
	}
      }

      // parton jet selections
      Int_t njets = pTJ->size();
      //      std::cout <<" ==> number of parton jets = " << njets << std::endl;

      Int_t nbj_barrel = 0;
      Int_t nj_barrel  = 0;
      Int_t nj_endcap  = 0;

      Int_t nbj_endcap = 0;
      Int_t nlj_endcap = 0;
      Int_t nlj_barrel = 0;

      // jet classification
      // b-jet in barrel with maximal pT
      Double_t pTmaxBJbarrel = -10.;
      Double_t etaBJbarrel(0.), phiBJbarrel(0.);
      // no tagged jet in barrel with maximal pT
      Double_t pTmaxLJbarrel = -10.;
      Double_t etaLJbarrel(0.), phiLJbarrel(0.);
      // jet in forward with maximal pT
      Double_t pTmaxJforward = -10.;
      Double_t etaJforward(0.), phiJforward(0.);

      TLorentzVector b_jet, ljet, fjet;
      TLorentzVector barjet;
      TLorentzVector b_jet1, b_jet2;  // needed in case of two b-jets in the 2nd bump
      TLorentzVector dijet_1stbump, dimuonj_1stbump, dimuondijet_1stbump;
      TLorentzVector dijet_2ndbump, dimuonj_2ndbump, dimuondijet_2ndbump;

      Double_t theta;

      for( Int_t i = 0; i < njets; i++) {

	if( (*pTJ)[i] < ptj_cut) {continue;}

	Double_t DRjmu1 = deltaR((*EtaJ)[i], (*PhiJ)[i], etamu1, phimu1);
	Double_t DRjmu2 = deltaR((*EtaJ)[i], (*PhiJ)[i], etamu2, phimu2);
	
	/*
	  std::cout <<"    ij = " << i+1
	  <<" ptj = " <<  (*pTJ)[i]
	  <<" ib = " << (*bjet)[i]
	  <<" DRjmu1 = " << DRjmu1
	  <<" DRjmu2 = " << DRjmu2 << std::endl;
	*/

	if (DRjmu1 <= Rmujetmatch || DRjmu2 <= Rmujetmatch) {continue;}

	if(fabs((*EtaJ)[i]) < etaj_cut && (*bjet)[i] == 1) 
	  {
	    nbj_barrel += 1;
	    if(nbj_barrel == 1) {
	      b_jet1.SetPx((*pTJ)[i]*cos((*PhiJ)[i]));
	      b_jet1.SetPy((*pTJ)[i]*sin((*PhiJ)[i]));
	      theta = 2. * atan(exp(-(*EtaJ)[i]));
	      b_jet1.SetPz((*pTJ)[i]/tan(theta));
	      b_jet1.SetE((*pTJ)[i]/sin(theta));
	    }
	    if(nbj_barrel == 2) {
	      b_jet2.SetPx((*pTJ)[i]*cos((*PhiJ)[i]));
	      b_jet2.SetPy((*pTJ)[i]*sin((*PhiJ)[i]));
	      theta = 2. * atan(exp(-(*EtaJ)[i]));
	      b_jet2.SetPz((*pTJ)[i]/tan(theta));
	      b_jet2.SetE((*pTJ)[i]/sin(theta));
	    }
	    if((*pTJ)[i] > pTmaxBJbarrel) {
	      pTmaxBJbarrel = (*pTJ)[i];
	      etaBJbarrel   = (*EtaJ)[i];
	      phiBJbarrel   = (*PhiJ)[i];
	    }
	  }

	if(fabs((*EtaJ)[i]) < etaj_cut) {nj_barrel += 1;}
	
	if(fabs((*EtaJ)[i]) > etaj_cut && fabs((*EtaJ)[i]) < 4.7) 
	  {
	    nj_endcap += 1;
	    if((*pTJ)[i] > pTmaxJforward) {
	      pTmaxJforward = (*pTJ)[i];
	      etaJforward   = (*EtaJ)[i];
	      phiJforward   = (*PhiJ)[i];
	    }
	  }

	// count numbers of b-jets and light jets in barrel and endcap
	if(fabs((*EtaJ)[i]) > etaj_cut && fabs((*EtaJ)[i]) < 4.7 && (*bjet)[i] == 1) {nbj_endcap += 1;}
	if(fabs((*EtaJ)[i]) > etaj_cut && fabs((*EtaJ)[i]) < 4.7 && (*bjet)[i] == 0) {nlj_endcap += 1;}
	if(fabs((*EtaJ)[i]) < etaj_cut && (*bjet)[i] == 0) 
	  {
	    nlj_barrel += 1;
	    if((*pTJ)[i] > pTmaxLJbarrel) {
	      pTmaxLJbarrel = (*pTJ)[i];
	      etaLJbarrel   = (*EtaJ)[i];
	      phiLJbarrel   = (*PhiJ)[i];
	    }
	  }
      }
      
      if(nbj_barrel >= 1) {
	b_jet.SetPx(pTmaxBJbarrel*cos(phiBJbarrel));
	b_jet.SetPy(pTmaxBJbarrel*sin(phiBJbarrel));
	theta = 2. * atan(exp(-etaBJbarrel));
	b_jet.SetPz(pTmaxBJbarrel/tan(theta));
	b_jet.SetE(pTmaxBJbarrel/sin(theta));
      }

      if(nlj_barrel >= 1) {
	ljet.SetPx(pTmaxLJbarrel*cos(phiLJbarrel));
	ljet.SetPy(pTmaxLJbarrel*sin(phiLJbarrel));
	theta = 2. * atan(exp(-etaLJbarrel));
	ljet.SetPz(pTmaxLJbarrel/tan(theta));
	ljet.SetE(pTmaxLJbarrel/sin(theta));
      }

      if(nj_endcap >= 1) {
	fjet.SetPx(pTmaxJforward*cos(phiJforward));
	fjet.SetPy(pTmaxJforward*sin(phiJforward));
	theta = 2. * atan(exp(-etaJforward));
	fjet.SetPz(pTmaxJforward/tan(theta));
	fjet.SetE(pTmaxJforward/sin(theta));
      }

      // 1st bump selections
      if( nbj_barrel == 1 && nj_barrel == 1 && nj_endcap >= 1) {
	n_1stbump += 1;
	w_1stbump += event_weight;

	dijet_1stbump = b_jet + fjet;
	dimuonj_1stbump = dimuon+fjet;
	dimuondijet_1stbump = dimuon+dijet_1stbump;

	if(fabs(costheta) < costheta_cut) {
	  hptw_1stbump->Fill(ptw);
	  hpt_1stmu_1stbump->Fill(ptmumax,event_weight);

	  hpt_2ndmu_1stbump->Fill(ptmumin,event_weight);
	  hpt_2ndmu_1stbump_scalar->Fill(ptmumin,weight_scalar);
	  hpt_2ndmu_1stbump_trnpve->Fill(ptmumin,weight_trnpve);
	  hpt_2ndmu_1stbump_trnmve->Fill(ptmumin,weight_trnmve);
	  hpt_2ndmu_1stbump_s2m1->Fill(ptmumin,weight_s2m1);

	  hcostheta_1stbump->Fill(costheta,event_weight);
	  hpt_dimuon_1stbump->Fill(dimuon.Pt(),event_weight);
	  hpt_dimuon_1stbump_s2m1->Fill(dimuon.Pt(),weight_s2m1);
	  
	  hnbj_endcap_1stbump->Fill(nbj_endcap*1.,event_weight);
	  hnlj_endcap_1stbump->Fill(nlj_endcap*1.,event_weight);
	  hnj_endcap_1stbump->Fill(nj_endcap*1.,event_weight);
	  
	  hmass_dijet_1stbump->Fill(dijet_1stbump.M(),event_weight);
	  hm_mumuj_1stbump->Fill(dimuonj_1stbump.M(),event_weight);
	  hmass_dimuondijet_1stbump->Fill(dimuondijet_1stbump.M(),event_weight);  
	  heta_dimuondijet_1stbump->Fill(dimuondijet_1stbump.Eta(),event_weight);  
	  
	  hpt_bjet_1stbump->Fill(b_jet.Pt(),event_weight);
	  heta_bjet_1stbump->Fill(b_jet.Eta(),event_weight);
	  hpt_fwdj_1stbump->Fill(fjet.Pt(),event_weight);
	  heta_fwdj_1stbump->Fill(fjet.Eta(),event_weight);
	}
      }

      // 2nd bump selections
      if( nbj_barrel >= 1 && nj_barrel == 2 && nj_endcap == 0) {

	if(nbj_barrel == 1) {
	  barjet = ljet;
	  //	  dphi_jj_2ndbump = deltaPhi(b_jet.Phi(),ljet.Phi());
	}
	if(nbj_barrel == 2) {
	  barjet = b_jet2;
	  //	  dphi_jj_2ndbump = deltaPhi(b_jet1.Phi(),b_jet2.Phi());
	}	
	dijet_2ndbump = b_jet + barjet;
	dimuonj_2ndbump = dimuon+barjet;
	dimuondijet_2ndbump = dimuon+dijet_2ndbump;

	Double_t Dphi_dimuon_dijet_2ndbump = deltaPhi(dimuon.Phi(),dijet_2ndbump.Phi());
	hdphi_dimuondijet_2ndbump->Fill(Dphi_dimuon_dijet_2ndbump,event_weight);
	if(Dphi_dimuon_dijet_2ndbump > dphi_mumu_dijet_2ndbump_cut) {

	  n_2ndbump += 1;
	  w_2ndbump += event_weight;
	
	  if(fabs(costheta) < costheta_cut) {
	    hptw_2ndbump->Fill(ptw);
	    hpt_1stmu_2ndbump->Fill(ptmumax,event_weight);

	    hpt_2ndmu_2ndbump->Fill(ptmumin,event_weight);
	    hpt_2ndmu_2ndbump_scalar->Fill(ptmumin,weight_scalar);
	    hpt_2ndmu_2ndbump_trnpve->Fill(ptmumin,weight_trnpve);
	    hpt_2ndmu_2ndbump_trnmve->Fill(ptmumin,weight_trnmve);

	    hcostheta_2ndbump->Fill(costheta,event_weight);
	    hpt_dimuon_2ndbump->Fill(dimuon.Pt(),event_weight);
	    
	    hnbj_barrel_2ndbump->Fill(nbj_barrel*1.,event_weight);
	    hnlj_barrel_2ndbump->Fill(nlj_barrel*1.,event_weight);
	    hnj_barrel_2ndbump->Fill(nj_barrel*1.,event_weight);
	    
	    hmass_dijet_2ndbump->Fill(dijet_2ndbump.M(),event_weight);
	    hm_mumuj_2ndbump->Fill(dimuonj_2ndbump.M(),event_weight);
	    hmass_dimuondijet_2ndbump->Fill(dimuondijet_2ndbump.M(),event_weight);
	    heta_dimuondijet_2ndbump->Fill(dimuondijet_2ndbump.Eta(),event_weight);
	    
	    hpt_bjet_2ndbump->Fill(b_jet.Pt(),event_weight);
	    heta_bjet_2ndbump->Fill(b_jet.Eta(),event_weight);
	    hpt_barj_2ndbump->Fill(barjet.Pt(),event_weight);
	    heta_barj_2ndbump->Fill(barjet.Eta(),event_weight);
	  }
	}
      }
   }

   std::cout <<"  Total number of events = " << n_total << std::endl;
   std::cout <<" " << std::endl;
   std::cout <<"  Cuts on muons: " << ptmucutmax <<", "<< ptmucutmin <<"  GeV" 
	     <<" eta1 < "<< etamucut1 
	     <<" eta2 < "<< etamucut2 
	     <<" costheta_cut = " << costheta_cut << std::endl;
   std::cout <<" " << std::endl;
   std::cout <<"  Cuts on jets, pt > " << ptj_cut <<" GeV, eta cut= "<< etaj_cut << std::endl; 
   std::cout <<" " << std::endl;
   std::cout <<"  Dphi dimuon-dijet > " << dphi_mumu_dijet_2ndbump_cut << std::endl;

   Double_t eff_muons = w_muons/w_total;

   Double_t eff_1stbump = w_1stbump/w_total;
   Double_t eff_2ndbump = w_2ndbump/w_total;

   Double_t eff_1stbump_q = w_1stbump_q/w_total;
   Double_t eff_2ndbump_q = w_2ndbump_q/w_total;

   std::cout <<"  " << std::endl;
   std::cout <<"  Efficiency of di-muon selections " 
	<<" = " << eff_muons 
	<<" ("<< n_muons<<"/"<< n_total <<")"<<  std::endl;
   std::cout <<" " << std::endl; 
    std::cout <<"  Efficiency of 1st bump selections with gen jets " 
	<<" = " << eff_1stbump 
	<<" ("<< n_1stbump<<"/"<< n_total <<")"<<  std::endl;
   std::cout <<" " << std::endl;
   std::cout <<"  Efficiency of 2nd bump selections with gen jets" 
	<<" = " << eff_2ndbump 
	<<" ("<< n_2ndbump<<"/"<< n_total <<")"<<  std::endl;

   std::cout << " " << std::endl;
   std::cout <<"  Efficiency of 1st bump selections with partons" 
	<<" = " << eff_1stbump_q 
	<<" ("<< n_1stbump_q<<"/"<< n_total <<")"<<  std::endl;
   std::cout <<" " << std::endl;
   std::cout <<"  Efficiency of 2nd bump selections with gen partons (no dphi cut)" 
	<<" = " << eff_2ndbump_q 
	<<" ("<< n_2ndbump_q<<"/"<< n_total <<")"<<  std::endl;

   TFile efile("test.root","recreate");

   hcosthetaw->Write();
   hptw->Write();

   hdrmub->Write();

   hcostheta_scalar->Write();
   hcostheta_longve->Write();
   hcostheta_trnpve->Write();
   hcostheta_trnmve->Write();
   hcostheta_s2m1->Write();

   // 1st bump histos
   hptw_1stbump->Write();
   hpt_1stmu_1stbump->Write();

   hpt_2ndmu_1stbump->Write();
   hpt_2ndmu_1stbump_scalar->Write();
   hpt_2ndmu_1stbump_trnpve->Write();
   hpt_2ndmu_1stbump_trnmve->Write();
   hpt_2ndmu_1stbump_s2m1->Write();

   hcostheta_1stbump->Write();
   hpt_dimuon_1stbump->Write();
   hpt_dimuon_1stbump_s2m1->Write();
   hnbj_endcap_1stbump->Write();
   hnlj_endcap_1stbump->Write();
   hnj_endcap_1stbump->Write();
   hmass_dijet_1stbump->Write();
   hm_mumuj_1stbump->Write();
   hmass_dimuondijet_1stbump->Write();
   heta_dimuondijet_1stbump->Write();
   hpt_bjet_1stbump->Write();
   heta_bjet_1stbump->Write();
   hpt_fwdj_1stbump->Write();
   heta_fwdj_1stbump->Write();
   //
   // 2nd bump histos
   hptw_2ndbump->Write();
   hpt_1stmu_2ndbump->Write();
  
   hpt_2ndmu_2ndbump->Write();
   hpt_2ndmu_2ndbump_scalar->Write();
   hpt_2ndmu_2ndbump_trnpve->Write();
   hpt_2ndmu_2ndbump_trnmve->Write();

   hcostheta_2ndbump->Write();
   hpt_dimuon_2ndbump->Write();
   hnbj_barrel_2ndbump->Write();
   hnlj_barrel_2ndbump->Write();
   hnj_barrel_2ndbump->Write();
   hmass_dijet_2ndbump->Write();
   hm_mumuj_2ndbump->Write();
   hmass_dimuondijet_2ndbump->Write();
   heta_dimuondijet_2ndbump->Write();
   hpt_bjet_2ndbump->Write();
   heta_bjet_2ndbump->Write();
   hpt_barj_2ndbump->Write();
   heta_barj_2ndbump->Write();
   hdphi_dimuondijet_2ndbump->Write();
   //
   // partons
   hpt_b->Write();
   heta_b->Write();
   hpt_q->Write();
   heta_q->Write();
   hpt_q_bcut->Write();
   heta_q_bcut->Write();
   heta_q_bcut_qptcut->Write();
   
   hpt_b_1stbump->Write();
   heta_b_1stbump->Write();
   hpt_b_2ndbump->Write();
   heta_b_2ndbump->Write();

   hpt_q_bcut_1stbump->Write();
   hpt_q_bcut_2ndbump->Write();
   
   efile.Close();

}
