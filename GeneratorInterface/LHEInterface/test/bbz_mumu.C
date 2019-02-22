#define bbz_mumu_cxx
#include "bbz_mumu.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

Double_t bbz_mumu::deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi; 
  return dphi;
}

Double_t bbz_mumu::deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t bbz_mumu::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t deta = deltaEta(eta1, eta2);
  Double_t dphi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

void bbz_mumu::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L bbz_mumu.C
//      Root > bbz_mumu t
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

   TH1F * hdr_bb = new TH1F("hdr_bb", "dr_bb", 25, 0., 5.0);

   TH1F * hpth62  = new TH1F( "hpth62", "pth62", 20, 0., 200.);
   TH1F * hpth62s = new TH1F( "hpth62s", "pth62s", 20, 0., 200.);
   TH1F * hpth62pw = new TH1F( "hpth62pw", "pth62pw", 40, 0., 200.);
   TH1F * hpth62nw = new TH1F( "hpth62nw", "pth62nw", 40, 0., 200.);

   TH1F * hptb1  = new TH1F( "hptb1", "pthb1", 30, 0., 150.);
   TH1F * hetab1 = new TH1F( "hetab1", "petab1", 25, -5.0, 5.0);

   TH1F * hptb2  = new TH1F( "hptb2", "pthb2", 30, 0., 150.);
   TH1F * hetab2 = new TH1F( "hetab2", "petab2", 25, -5.0, 5.0);

   TH1F * hptj  = new TH1F( "hptj", "pthj", 30, 0., 150.);
   TH1F * hetaj = new TH1F( "hetaj", "petaj", 25, -5.0, 5.0);

   TH1F * hnlj = new TH1F( "hnlj", "nlj", 5, 0., 5.0);

   Double_t ptmucutmax = 20.;
   Double_t ptmucutmin = 20.;
   Double_t etamucut1  = 2.4;
   Double_t etamucut2  = 2.4;
   Double_t ptbcut = 30.0;
   Double_t ptjcut = 30.0;
   Double_t etabcut = 2.4;
   Double_t etajcut = 4.7;

   Double_t Rmujetmatch = 0.5;

   Double_t sum_weights_0[39];
   Double_t sum_weights_s[39];
   Double_t eff_weights[39];

   // countings for PY8
   Int_t ntot = 0;
   Int_t nsel_taus = 0;
   Int_t nsel_bjet = 0;
   Int_t nsel_jets = 0;

   // countings for aMC@NLO
   Double_t w_ntot = 0.;
   Double_t w_nsel_taus = 0.;
   Double_t w_nsel_bjet = 0.;
   Double_t w_nsel_jets = 0.;

   for( Int_t i = 0; i <= 38; i++) {
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
      // MG5
      //      Double_t weight_evt = (*evtweight)[0];
      //      cout <<"   weight = " <<  weight_evt << endl;
      // PY8 and POWHEG
      Double_t weight_evt = 1.;
      hpth62->Fill(pth62,weight_evt);

      // sum of weights before selections
      if( evtweight->size() != 0) {
	for( Int_t i = 0; i <= 38; i++) {
	  sum_weights_0[i] += (*evtweight)[i];
	}
      }

      ntot += 1.;
      if( evtweight->size() != 0) {
	w_ntot += (*evtweight)[0];
      }

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
      nsel_taus += 1;

      if(fabs(etab1) < etabcut && fabs(etab2) < etabcut && ptb1 > ptbcut && ptb2 > ptbcut) {
	double drbb = deltaR(etab1,phib1,etab2,phib2);
        hdr_bb->Fill(drbb,weight_evt);
      }


      if(evtweight->size() != 0) {
	w_nsel_taus += (*evtweight)[0];
      }

      // jets
      Int_t nalljets = pTJ->size();
      Int_t nbjets = 0;
      Int_t njets = 0;
      Int_t nljets = 0;

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
	  // count light jets
	  if( (*pTJ)[i] > ptjcut && fabs((*EtaJ)[i]) < etajcut && (*bjet)[i] == 0 ) {
	    nljets += 1;
	    if(nljets == 1) {
	      ptj = (*pTJ)[i];
	      etaj = (*EtaJ)[i];
	    }
	  }
	  // count light and b jets
	  if( (*pTJ)[i] > ptjcut && fabs((*EtaJ)[i]) < etajcut ) {njets += 1;}
	}
      }

      if(nbjets < 1) {continue;}
      nsel_bjet += 1;
      hpth62s->Fill(pth62,weight_evt);
      if(evtweight->size() != 0) {
	w_nsel_bjet += (*evtweight)[0];
      }
      
      hnlj->Fill(1.*nljets,weight_evt);
   
      if(nljets >= 1) {
	hptj->Fill(ptj,weight_evt);
	hetaj->Fill(etaj,weight_evt);
      }

      if(njets >= 2) {continue;}
      nsel_jets += 1;
      if(evtweight->size() != 0) {
	w_nsel_jets += (*evtweight)[0];
      }
      // weights for selected events
      if( evtweight->size() != 0) {
	for( Int_t i = 0; i <= 38; i++) {
	  sum_weights_s[i] += (*evtweight)[i];
	}
      }
   }

   // Event selections
   cout <<"  Event selections: " << endl;
   cout <<"         max pT tauh > " << ptmucutmax <<"  eta < " << etamucut1 << endl;
   cout <<"         min pT tauh > " << ptmucutmin <<"  eta < " << etamucut2 << endl;
   cout <<"         at least one b-jet pT > " << ptbcut <<" eta < " << etabcut << endl; 
   cout <<"         no more than one jet pT > " << ptjcut <<" eta < " << etajcut << endl; 
   cout <<" " << endl;
   //
   cout <<" Total number of events processed: " << ntot << endl;
   cout <<" Passed tauh selections: " << nsel_taus << endl;   
   cout <<" Passed >= 1b selection: " << nsel_bjet << endl;   
   cout <<" Passed <= 1j selection: " << nsel_jets << endl;   
   cout <<" " << endl;
  
   if( evtweight->size() == 0) {
     cout <<"      Efficiency of every selection for PYTHIA8/POWHEG: " << endl;
     cout <<" tauh selections: " << 1.*nsel_taus/ntot << endl;   
     cout <<" >= 1b selection: " << 1.*nsel_bjet/nsel_taus << endl;   
     cout <<" <= 1j selection: " << 1.*nsel_jets/nsel_bjet << endl;   
     cout <<" Total efficiency " << 1.*nsel_jets/ntot << endl;   
     cout <<" " << endl;
   }

   if( evtweight->size() != 0) {

     cout <<"  ============ THIS IS CALCULATIONS FOR aMC@NLO =============" << endl;
     // Efficiencies
     Double_t eff0 = sum_weights_s[0]/sum_weights_0[0];
     cout <<"           Efficiency = " << eff0 << endl;
     cout <<" " << endl;

     // pdf uncertainty
     Double_t effpdf_min = 10000.;
     Double_t effpdf_max = 0.;
     for( Int_t i = 9; i <= 38; i++) {
       Double_t effi = sum_weights_s[i]/sum_weights_0[i];
       //     cout <<" i = " << i <<" effi = " << effi << endl; 
       if(effi < effpdf_min) {
	 effpdf_min = effi;
       }
       if(effi > effpdf_max) {
	 effpdf_max = effi;
       }
     }

     // QCD scal uncertainty
     Double_t effqcd_min = 10000.;
     Double_t effqcd_max = 0.;
     for( Int_t i = 0; i <= 8; i++) {
       if(i == 5 || i == 7) {continue;}
       Double_t effi = sum_weights_s[i]/sum_weights_0[i];
       //     cout <<" i = " << i <<" effi = " << effi << endl; 
       if(effi < effqcd_min) {
	 effqcd_min = effi;
       }
       if(effi > effqcd_max) {
	 effqcd_max = effi;
       }
     }

     cout <<"           QCD uncertainty" << endl;
     cout <<"   effqcd max = " << effqcd_max <<" effqcd min = " << effqcd_min << endl;
     cout <<" in %:  + " << 100.*(effqcd_max-eff0)/eff0 <<" - " << 100.*(eff0-effqcd_min)/eff0 << endl;
     
     cout <<" " << endl;
     
     cout <<"           PDF uncertainty" << endl;
     cout <<"   effpdf max = " << effpdf_max <<" effpdf min = " << effpdf_min << endl;
     cout <<" in %:  + " << 100.*(effpdf_max-eff0)/eff0 <<" - " << 100.*(eff0-effpdf_min)/eff0 << endl;

     cout <<" " << endl;
     cout <<"      Efficiency of every selection for aMC@NLO: " << endl;
     cout <<" tauh selections: " << w_nsel_taus/w_ntot << endl;   
     cout <<" >= 1b selection: " << w_nsel_bjet/w_nsel_taus << endl;   
     cout <<" <= 1j selection: " << w_nsel_jets/w_nsel_bjet << endl;   
     cout <<" " << endl;
   }

   hdr_bb->Draw("hist");
   //   TFile efile("bbz_mumu_alphadef025_pdfrwt_histos.root","recreate");
   TFile efile("bbz_mumu_alphadef_x_2_pdfrwt_histos.root","recreate");
   //   TFile efile("bbz_mumu_alphadef_div_sqrt2_pdfrwt_histos.root","recreate");

   //   TFile efile("bbz_mumu_via_gbb_hist.root","recreate");
   //   TFile efile("bbz_mumu_via_gg_hist.root","recreate");
   hpth62->Write();  
   hpth62s->Write();  
   hpth62pw->Write();
   hpth62nw->Write();
   hdr_bb->Write();
   hptb1->Write();
   hetab1->Write();
   hptb2->Write();
   hetab2->Write();
   hptj->Write();
   hetaj->Write();
   hnlj->Write();
   efile.Close();
}
