#define bbh_tautau_cxx
#include "bbh_tautau.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

Double_t bbh_tautau::deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi; 
  return dphi;
}

Double_t bbh_tautau::deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t bbh_tautau::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t deta = deltaEta(eta1, eta2);
  Double_t dphi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

void bbh_tautau::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L bbh_tautau.C
//      Root > bbh_tautau t
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

   TH1F * hpth62 = new TH1F( "hpth62", "pth62", 20, 0., 400.);
   TH1F * hpth62pw = new TH1F( "hpth62pw", "pth62pw", 20, 0., 400.);
   TH1F * hpth62nw = new TH1F( "hpth62nw", "pth62nw", 20, 0., 400.);

   TH1F * hptb1  = new TH1F( "hptb1", "pthb1", 30, 0., 150.);
   TH1F * hetab1 = new TH1F( "hetab1", "petab1", 25, -5.0, 5.0);

   TH1F * hptb2  = new TH1F( "hptb2", "pthb2", 30, 0., 150.);
   TH1F * hetab2 = new TH1F( "hetab2", "petab2", 25, -5.0, 5.0);

   TH1F * hptj  = new TH1F( "hptj", "pthj", 30, 0., 150.);
   TH1F * hetaj = new TH1F( "hetaj", "petaj", 25, -5.0, 5.0);

   TH1F * hnlj = new TH1F( "hnlj", "nlj", 5, 0., 5.0);

   Double_t ptmucutmax = 0.;
   Double_t ptmucutmin = 0.;
   Double_t etamucut1  = 12.1;
   Double_t etamucut2  = 12.1;
   Double_t ptbcut = 25.0;
   Double_t ptjcut = 25.0;
   Double_t etabcut = 2.5;
   Double_t etajcut = 4.7;

   Double_t Rmujetmatch = 0.4;

   Double_t sum_weights_0[40];
   Double_t sum_weights_s[40];
   Double_t eff_weights[40];
   // total number of events processed
   Int_t ntot = 0.;
   // number of selected events
   Int_t nsel = 0.;
   // countings event weights (for aMC@NLO != 1.)
   Double_t w_ntot = 0.;
   Double_t w_nsel_muons  = 0.;
   Double_t w_nsel_bjet   = 0.;
   Double_t w_nselbr_jets = 0.;
   Double_t w_nselfw_jets = 0.;

   for( Int_t i = 0; i <= 39; i++) {
     sum_weights_0[i] = 0.;
     sum_weights_s[i] = 0.;
   }

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
   
   // if (Cut(ientry) < 0) continue;

      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // weight to fill histograms
      // PY8 and POWHEG
      Double_t weight_evt = 1.;
      //      cout <<"  weight size = " << evtweight->size() << endl; 
      // MG5_aMC@NLO
  //     if( evtweight->size() != 0) {
	// weight_evt = (*evtweight)[0];
	// for( Int_t i = 0; i < evtweight->size(); i++) {
	//   sum_weights_0[i] += (*evtweight)[i];
	// }
  //     }
      
      ntot += 1.;

      w_ntot += weight_evt;

      hpth62->Fill(pth62,weight_evt);

      if(weight_evt > 0.) {hpth62pw->Fill(pth62,weight_evt);}
      if(weight_evt < 0.) {hpth62nw->Fill(pth62,weight_evt);}

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

      // selections
      if(ptmumax<ptmucutmax || ptmumin<ptmucutmin || fabs(etamu1)>etamucut1 || fabs(etamu2)>etamucut2) {continue;}
      w_nsel_muons += (*evtweight)[0];

      // jets
      Int_t nalljets = pTJ->size();
      Int_t nbjets = 0;
      Int_t njets_br = 0;
      Int_t njets_fw = 0;

      Double_t ptj = 0.;
      Double_t etaj = 0.;

      for( Int_t i = 0; i < nalljets; i++) {

	Double_t DRjmu1 = deltaR((*EtaJ)[i], (*PhiJ)[i], etamu1, phimu1);
	Double_t DRjmu2 = deltaR((*EtaJ)[i], (*PhiJ)[i], etamu2, phimu2);
	// take jets not coinciding with tau_h
	if (DRjmu1 >= Rmujetmatch && DRjmu2 >= Rmujetmatch) {
	  // count b-jets
	  if( (*pTJ)[i] > ptbcut && fabs((*EtaJ)[i]) < etabcut && (*bjet)[i] == 1 ) {
	    nbjets += 1;
	    if(nbjets == 1) {
	      hptb1->Fill((*pTJ)[i],weight_evt);
	      hetab1->Fill((*EtaJ)[i],weight_evt);
	    }
	    if(nbjets == 2) {
	      hptb2->Fill((*pTJ)[i],weight_evt);
	      hetab2->Fill((*EtaJ)[i],weight_evt);
	    }
	  }
	  // count jets in barrel
	  if( (*pTJ)[i] > ptjcut && fabs((*EtaJ)[i]) < etabcut) {
	    njets_br += 1;
	    if(njets_br == 1) {
	      ptj = (*pTJ)[i];
	      etaj = (*EtaJ)[i];
	    }
	  }
	  // count jets in forward
	  if( (*pTJ)[i] > ptjcut && fabs((*EtaJ)[i]) > etabcut && fabs((*EtaJ)[i]) < etajcut) {
	    njets_fw += 1;
	    if(njets_fw == 1) {
	      ptj = (*pTJ)[i];
	      etaj = (*EtaJ)[i];
	    }
	  }
	}
      }

      if(nbjets != 2) {continue;}
      w_nsel_bjet += (*evtweight)[0];
      
      if(njets_br < 0) {continue;}
      w_nselbr_jets += (*evtweight)[0];

      if(njets_fw < 0) {continue;}
      w_nselfw_jets += (*evtweight)[0];

      nsel += 1.;

      // weights for selected events
  //     if( evtweight->size() != 0) {
	// for( Int_t i = 0; i < evtweight->size(); i++) {
	//   sum_weights_s[i] += (*evtweight)[i];
	// }
  //     }
   }

   // Event selections
   cout <<"  Event selections: " << endl;
   cout <<"         max pT muon > " << ptmucutmax <<"  |eta| < " << etamucut1 << endl;
   cout <<"         min pT muon > " << ptmucutmin <<"  |eta| < " << etamucut2 << endl;
   cout <<"         at least one b-jet pT > " << ptbcut <<", |eta| < " << etabcut << endl; 
   //   cout <<"         only one jet pT > " << ptjcut <<", |eta| < " << etabcut << endl; 
   //   cout <<"         at least one jet pT > " << ptjcut 
   //	<<", |eta| > " << etabcut 
   //	<<" and < " << etajcut << endl; 
   cout <<" " << endl;
   //
   cout <<" Total number of events processed: " << ntot << endl;
   cout <<" Total number of events selected:  " << nsel << endl;
   cout <<" " << endl;
  
   cout <<"      Sum of event weights at each selection : " << endl;
   cout <<" before selections: " << w_ntot << endl;
   cout <<" muons selections: " << w_nsel_muons << endl;   
   cout <<" >= 1b selection: " << w_nsel_bjet << endl;   
   cout <<" 1 barrel jet   : " << w_nselbr_jets << endl;   
   cout <<" > 0 forward jet: " << w_nselfw_jets << endl;   
   cout <<" " << endl;

   cout <<"      Efficiency of selection relative to previous selection: " << endl;
   cout <<" muons selections: " << w_nsel_muons/w_ntot << endl;   
   cout <<" >= 1b selection: " << w_nsel_bjet/w_nsel_muons << endl;   
   cout <<" 1 barrel jet   : " << w_nselbr_jets/w_nsel_bjet << endl;   
   cout <<" > 0 forward jet: " << w_nselfw_jets/w_nselbr_jets << endl;   
   cout <<" " << endl;
   cout <<"      Total efficiency: " << w_nselfw_jets/w_ntot << endl;

   Double_t eff0   = w_nselfw_jets/w_ntot;  
   Double_t sigmat = w_ntot/ntot;
   Double_t sigmas = w_nselfw_jets/ntot;
   cout <<"  " << endl;
   cout <<"  ===> sigma before selection (pb) = " << sigmat << endl;
   cout <<"  ===> sigma  after selection (pb) = " << sigmas << endl;
   cout <<"  ===> eff = " << eff0 << endl;

   // pdf uncertainty
   Double_t sigmatpdf_min = 100000.;
   Double_t sigmatpdf_max = 0.;

   Double_t sigmaspdf_min = 100000.;
   Double_t sigmaspdf_max = 0.;

   Double_t effpdf_min = 100000.;
   Double_t effpdf_max = 0.;

   cout <<"   PDFs " << endl;
   for( Int_t i = 9; i <= 39; i++) {
     Double_t effi = sum_weights_s[i]/sum_weights_0[i];
     Double_t sigmati = sum_weights_0[i]/ntot;
     Double_t sigmasi = sum_weights_s[i]/ntot;
     //     cout <<" i = " << i <<" effi = " << effi <<" sigmasi (pb) = " << sigmasi << endl; 

     if(effi < effpdf_min) {
       effpdf_min = effi;
     }
     if(effi > effpdf_max) {
       effpdf_max = effi;
     }

     if(sigmasi < sigmaspdf_min) {
       sigmaspdf_min = sigmasi;
     }
     if(sigmasi > sigmaspdf_max) {
       sigmaspdf_max = sigmasi;
     }

     if(sigmati < sigmatpdf_min) {
       sigmatpdf_min = sigmati;
     }
     if(sigmati > sigmatpdf_max) {
       sigmatpdf_max = sigmati;
     }

   }
   
   cout <<" SCALE " << endl;
   // QCD scal uncertainty
   Double_t effqcd_min = 100000.;
   Double_t effqcd_max = 0.;

   Double_t sigmasqcd_min = 100000.;
   Double_t sigmasqcd_max = 0.;

   Double_t sigmatqcd_min = 100000.;
   Double_t sigmatqcd_max = 0.;

   for( Int_t i = 0; i <= 8; i++) {
     if(i == 5 || i == 7) {continue;}
     Double_t effi = sum_weights_s[i]/sum_weights_0[i];
     Double_t sigmati = sum_weights_0[i]/ntot;
     Double_t sigmasi = sum_weights_s[i]/ntot;
     //     cout <<" i = " << i <<" effi = " << effi <<" sigmasi (pb) = " << sigmasi << endl; 
     if(effi < effqcd_min) {
       effqcd_min = effi;
     }
     if(effi > effqcd_max) {
       effqcd_max = effi;
     }
     if(sigmasi < sigmasqcd_min) {
       sigmasqcd_min = sigmasi;
     }
     if(sigmasi > sigmasqcd_max) {
       sigmasqcd_max = sigmasi;
     }

     if(sigmati < sigmatqcd_min) {
       sigmatqcd_min = sigmati;
     }
     if(sigmati > sigmatqcd_max) {
       sigmatqcd_max = sigmati;
     }

   }
   
   cout <<"           QCD uncertainty" << endl;

   cout <<"      Cross-section before selections " << endl;
   cout <<"   sigmatqcd max = " << sigmatqcd_max <<" sigmatqcd min = " << sigmatqcd_min << endl;
   cout <<" in %:  + " << 100.*(sigmatqcd_max-sigmat)/sigmat 
	       <<" - " << 100.*(sigmat-sigmatqcd_min)/sigmat << endl;

   cout <<"      Cross-section after selections " << endl;
   cout <<"   sigmasqcd max = " << sigmasqcd_max <<" sigmasqcd min = " << sigmasqcd_min << endl;
   cout <<" in %:  + " << 100.*(sigmasqcd_max-sigmas)/sigmas 
	       <<" - " << 100.*(sigmas-sigmasqcd_min)/sigmas << endl;

   cout <<"      Efficiency " << endl;
   cout <<"   effqcd max = " << effqcd_max <<" effqcd min = " << effqcd_min << endl;
   cout <<" in %:  + " << 100.*(effqcd_max-eff0)/eff0 <<" - " << 100.*(eff0-effqcd_min)/eff0 << endl;
   cout <<" " << endl;

   cout <<"           PDF uncertainty" << endl;

   cout <<"      Cross-section before selections " << endl;
   cout <<"   sigmatpdf max = " << sigmatpdf_max <<" sigmatpdf min = " << sigmatpdf_min << endl;
   cout <<" in %:  + " << 100.*(sigmatpdf_max-sigmat)/sigmat 
	       <<" - " << 100.*(sigmat-sigmatpdf_min)/sigmat << endl;

   cout <<"      Cross-section after selections " << endl;
   cout <<"   sigmaspdf max = " << sigmaspdf_max <<" sigmaspdf min = " << sigmaspdf_min << endl;
   cout <<" in %:  + " << 100.*(sigmaspdf_max-sigmas)/sigmas 
	       <<" - " << 100.*(sigmas-sigmaspdf_min)/sigmas << endl;

   cout <<"      Efficiency " << endl;
   cout <<"   effpdf max = " << effpdf_max <<" effpdf min = " << effpdf_min << endl;
   cout <<" in %:  + " << 100.*(effpdf_max-eff0)/eff0 <<" - " << 100.*(eff0-effpdf_min)/eff0 << endl;
   cout <<" " << endl;
   

   hpth62->Draw("hist");
   //   TFile efile("bbh_tautau_700gev_mg5_yb2_alphadef025_pdfrwt_hist.root","recreate");
   //   TFile efile("bbh_tautau_700gev_mg5_yb2_alphadef025_div_sqrt2_hist.root","recreate");
   TFile efile("test_hist.root","recreate");

   //   TFile efile("bbh_tautau_700gev_py8_hist.root","recreate");
   //   TFile efile("bbh_tautau_700gev_powheg_hist.root","recreate");
   hpth62->Write();  
   hpth62pw->Write();
   hpth62nw->Write();
   hptb1->Write();
   hetab1->Write();
   hptb2->Write();
   hetab2->Write();
   hptj->Write();
   hetaj->Write();
   hnlj->Write();
   efile.Close();
}
