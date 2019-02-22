#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

Double_t test::deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi; 
  return dphi;
}

Double_t test::deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t test::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t deta = deltaEta(eta1, eta2);
  Double_t dphi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

void test::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L test.C
//      Root > test t
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

   MyVariables();

   TH1F * hptj1   = new TH1F( "hptj1", "ptj1", 30, 0.0, 150.0);
   TH1F * hptj2   = new TH1F( "hptj2", "ptj2", 30, 0.0, 150.0);
   TH1F * hdetajj = new TH1F( "hdetajj", "detajj", 35, 0.0, 7.0);

   TH1F * hdphijj = new TH1F( "hdphijj", "dphijj", 18, -180., 180.0);
   TH1F * hnjets  = new TH1F( "hnjets", "njets", 7, 0., 7.0);

   TH1F * hdphipp = new TH1F( "hdphipp", "dphipp", 18, -180., 180.0);

   Double_t ptjcut = 25.0;
   Double_t etajcut = 4.7;
   Double_t Rmujetmatch = 0.5;
   Double_t Detajj_cut = 3.0;

   // countings for aMC@NLO
   Double_t sum_weights_0[9];
   Double_t sum_weights_s[9];
   Double_t eff_weights[9];

   // countings for PY8
   Double_t ntot = 0;
   Double_t ntotpw = 0;
   Double_t ntotnw = 0;
   Double_t nsel = 0;

   for( Int_t i = 0; i <= 8; i++) {
     sum_weights_0[i] = 0.;
     sum_weights_s[i] = 0.;
   }

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      // weight to fill histograms
      // MG5
      //      cout <<"  evt size = " <<  evtweight->size() << endl;

      //      Double_t weight_evt = (*evtweight)[0];
      Double_t weight_evt = 1.;

      //      cout <<"  N partons = " <<  idp->size() << endl;
      
      if((*ptp)[0] > ptjcut && (*ptp)[1] > ptjcut) {

	Double_t pTp1  = (*ptp)[0];
	Double_t etap1 = (*etap)[0];
	Double_t phip1 = (*phip)[0];
	
	Double_t pTp2  = (*ptp)[1];
	Double_t etap2 = (*etap)[1];
	Double_t phip2 = (*phip)[1];
	
	if(etap1 > etap2) {
	  pTp1  = (*ptp)[1];
	  etap1 = (*etap)[1];
	  phip1 = (*phip)[1];
	  
	  pTp2  = (*ptp)[0];
	  etap2 = (*etap)[0];
	  phip2 = (*phip)[0];
	} 
	
	Double_t Detap1p2 = etap1 - etap2;
	Double_t pi = 3.1415927;
	Double_t Dphip1p2 = (phip1 - phip2)*180.0/pi;

	if(fabs(Detap1p2) < Detajj_cut) {continue;}
	//	hdphipp->Fill(Dphip1p2,weight_evt);
      }  

      // sum of weights before selections
      if( evtweight->size() != 0) {
	for( Int_t i = 0; i <= 8; i++) {
	  sum_weights_0[i] += (*evtweight)[i];
	}
      }

      ntot += 1.;
      //      if((*evtweight)[0] > 0.) {ntotpw += 1.;}
      //      if((*evtweight)[0] < 0.) {ntotnw += 1.;}

      Int_t nalljets = pTJ->size();
      Int_t njets = 0;

      EtaJet->clear();
      PhiJet->clear();
      pTJet->clear();

      //      cout <<"  ===> Event = " << endl;

      for( Int_t i = 0; i < nalljets; i++) {
	Double_t DRjmu1 = deltaR((*EtaJ)[i], (*PhiJ)[i], etamu1, phimu1);
	Double_t DRjmu2 = deltaR((*EtaJ)[i], (*PhiJ)[i], etamu2, phimu2);
	// take jets not coinciding with muons
	if (DRjmu1 >= Rmujetmatch && DRjmu2 >= Rmujetmatch) {
	  // count light jets
	  if( (*pTJ)[i] > ptjcut && fabs((*EtaJ)[i]) < etajcut) {
	    njets += 1;
	    //	    cout <<"  jet i = " << i <<" pT = " << (*pTJ)[i] <<" eta = " << (*EtaJ)[i] << endl;
	    EtaJet->push_back((*EtaJ)[i]);
	    PhiJet->push_back((*PhiJ)[i]);
	    pTJet->push_back((*pTJ)[i]);
	  }
	}
      }

      Int_t nsjets = pTJet->size();
      for( Int_t i = 0; i < nsjets; i++) {
	//	cout <<"  jet j = " << i <<" pT = " << (*pTJet)[i] <<" eta = " << (*EtaJet)[i] << endl;
      }

      if(njets < 2) {continue;}
      
      Double_t pTj1  = (*pTJet)[0];
      Double_t etaj1 = (*EtaJet)[0];
      Double_t phij1 = (*PhiJet)[0];

      Double_t pTj2  = (*pTJet)[1];
      Double_t etaj2 = (*EtaJet)[1];
      Double_t phij2 = (*PhiJet)[1];

      if(etaj1 > etaj2) {
	pTj1  = (*pTJet)[1];
	etaj1 = (*EtaJet)[1];
	phij1 = (*PhiJet)[1];

	pTj2  = (*pTJet)[0];
	etaj2 = (*EtaJet)[0];
	phij2 = (*PhiJet)[0];
      } 

      Double_t Detaj1j2 = etaj1 - etaj2;
      Double_t pi = 3.1415927;
      Double_t Dphij1j2 = (phij1 - phij2)*180.0/pi;

      hptj1->Fill((*pTJet)[0],weight_evt);
      hptj2->Fill((*pTJet)[1],weight_evt);
      hdetajj->Fill(fabs(Detaj1j2),weight_evt);
 
      if(fabs(Detaj1j2) < Detajj_cut) {continue;}

      nsel += 1;

      hdphijj->Fill(Dphij1j2,weight_evt);
      hdphipp->Fill(Dphip1p2,weight_evt);
      hnjets->Fill(njets,weight_evt);

      // weights for selected events
      if( evtweight->size() != 0) {
	for( Int_t i = 0; i <= 8; i++) {
	  sum_weights_s[i] += (*evtweight)[i];
	}
      }
   }

   hdetajj->Draw("hist");
   hdphijj->SetMinimum(0.);
   hdphijj->Draw("hist");
   //   hdphipp->Draw("hist");
   hnjets->SetMinimum(0.);
   hnjets->Draw("hist");

   TFile efile("CPtestMG5_H_LO_hist.root","recreate");
   //   TFile efile("CPtestMG5_H_hist.root","recreate");
   //   TFile efile("CPtestMG5_A_hist.root","recreate");
   //   TFile efile("CPtestMG5_CPmix_hist.root","recreate");
   hptj1->Write();
   hptj2->Write();
   hdetajj->Write();
   hdphijj->Write();
   hdphipp->Write();
   hnjets->Write();
   efile.Close();

   // Event selections
   cout <<"  Event selections: " << endl;
   cout <<"         more than one jet pT > " << ptjcut <<" eta < " << etajcut << endl; 
   cout <<"         Detajj > " << Detajj_cut << endl; 
   cout <<" " << endl;
   //
   cout <<" Total number of events processed: " << ntot << endl;
   cout <<"      with positive weights: " << ntotpw << endl;
   cout <<"      with negative weights: " << ntotnw << endl;
   cout <<" Passed selection: " << nsel << endl;   
   cout <<" " << endl;
  
   if( evtweight->size() == 0) {
     cout <<"      Efficiency of every selection for PYTHIA8/POWHEG: " << endl;
     cout <<" Total efficiency " << 1.*nsel/ntot << endl;   
     cout <<" " << endl;
   }

   if( evtweight->size() != 0) {

     cout <<"  ============ THIS IS CALCULATIONS FOR aMC@NLO =============" << endl;
     // Efficiencies
     Double_t eff0 = sum_weights_s[0]/sum_weights_0[0];
     Double_t sigma0 = sum_weights_0[0]/ntot;
     cout <<" sigma0 = " <<  sigma0
	  <<" sigmas = " << sum_weights_s[0]/ntot << endl;  
     cout <<"           Efficiency = " << eff0 << endl;
     cout <<" " << endl;

     // QCD scal uncertainty
     Double_t effqcd_min = 10000.;
     Double_t effqcd_max = 0.;

     Double_t sigma0_min = 10000.;
     Double_t sigma0_max = 0.;

     for( Int_t i = 0; i <= 8; i++) {
       if(i == 5 || i == 7) {continue;}
       Double_t effi = sum_weights_s[i]/sum_weights_0[i];
       Double_t sigma = sum_weights_0[i]/ntot;
       cout <<" i = " << i <<" effi = " << effi <<" sigma = " << sigma << endl;
       if(effi < effqcd_min) {
	 effqcd_min = effi;
       }
       if(effi > effqcd_max) {
	 effqcd_max = effi;
       }

       if(sigma < sigma0_min) {
	 sigma0_min = sigma;
       }
       if(sigma > sigma0_max) {
	 sigma0_max = sigma;
       }
     }

     cout <<"           total cross-section uncertainty due to scale" << endl;
     cout <<"   sigma0 max = " << sigma0_max <<" sigma0 min = " << sigma0_min << endl;
     cout <<" in %:  + " << 100.*(sigma0_max-sigma0)/sigma0 
	  <<" - " << 100.*(sigma0-sigma0_min)/sigma0 << endl;

     cout <<"           selection uncertainty due to scale" << endl;
     cout <<"   effqcd max = " << effqcd_max <<" effqcd min = " << effqcd_min << endl;
     cout <<" in %:  + " << 100.*(effqcd_max-eff0)/eff0 <<" - " << 100.*(eff0-effqcd_min)/eff0 << endl;
     
     cout <<" " << endl;
     
     /*
     cout <<"           PDF uncertainty" << endl;
     cout <<"   effpdf max = " << effpdf_max <<" effpdf min = " << effpdf_min << endl;
     cout <<" in %:  + " << 100.*(effpdf_max-eff0)/eff0 <<" - " << 100.*(eff0-effpdf_min)/eff0 << endl;

     cout <<" " << endl;
     cout <<"      Efficiency of every selection for aMC@NLO: " << endl;
     cout <<" tauh selections: " << w_nsel_taus/w_ntot << endl;   
     cout <<" >= 1b selection: " << w_nsel_bjet/w_nsel_taus << endl;   
     cout <<" <= 1j selection: " << w_nsel_jets/w_nsel_bjet << endl;   
     cout <<" " << endl;
     */
   }
}
