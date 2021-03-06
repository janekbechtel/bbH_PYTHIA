#define vlq_cxx
#include "vlq.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>

void vlq::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L vlq.C
//      Root > vlq t
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

   vector<double> *EtaJ         = new std::vector<double>();
   vector<double> *PhiJ         = new std::vector<double>();
   vector<double> *pTJ          = new std::vector<double>();
   vector<double> *bjet         = new std::vector<double>();

   // eta of not b-tagged jet after muon selections
   TH1F * heta_ljet  = new TH1F("heta_ljet", "eta_ljet", 20, -5.0, 5.0);

   TH1F * hdr_dimuon_1stbump = new TH1F("hdr_dimuon_1stbump", "dr_dimuon_1stbump", 25, 0., 2.5);
   TH1F * hdr_dimuon_2ndbump = new TH1F("hdr_dimuon_2ndbump", "dr_dimuon_2ndbump", 25, 0., 2.5);

   TH1F * hisomu1  = new TH1F("hisomu1", "isomu1", 20, 0.0, 0.5);
   TH1F * hisomu2  = new TH1F("hisomu2", "isomu2", 20, 0.0, 0.5);

   TH1F * hmass_mumubb = new TH1F("hmass_mumubb", "mass_mumubb", 40, 0., 200.);

   // DR mumu, bb, min DRmub
   TH1F * hdr_mumu = new TH1F("hdr_mumu", "dr_mumu", 25, 0., 2.5);
   TH1F * hdr_bb = new TH1F("hdr_bb", "dr_bb", 25, 0., 2.5);
   TH1F * hdr_mumu_bb = new TH1F("hdr_mumu_bb", "dr_mumu_bb", 25, 0., 2.5);
   TH1F * hpt_dimuon = new TH1F("hpt_dimuon", "pt_dimuon", 20, 0., 200.);
   TH1F * hpt_dibquark = new TH1F("hpt_dibquark", "pt_dibquark", 20, 0., 200.);
   TH1F * hpt_2ndmu = new TH1F("hpt_2ndmu", "pt_2ndmu", 50, 0., 50.);
   TH1F * hcostheta = new TH1F("hcostheta", "costheta", 50, -1.0, 1.0);

   // histos for 1st bump
   TH1F * hpt_1stmu_1stbump = new TH1F("hpt_1stmu_1stbump", "pt_1stmu_1stbump", 150, 0., 150.);
   TH1F * hpt_2ndmu_1stbump = new TH1F("hpt_2ndmu_1stbump", "pt_2ndmu_1stbump", 80, 0., 80.);
   TH1F * hpt_dimuon_1stbump = new TH1F("hpt_dimuon_1stbump", "pt_dimuon_1stbump", 20, 0., 200.);
   TH1F * hcostheta_1stbump = new TH1F("hcostheta_1stbump", "costheta_1stbump", 10, -1.0, 1.0);
   TH1F * hmass_mumu_1stbump = new TH1F("hmass_mumu_1stbump", "mass_mumu_1stbump", 100, 0., 100.);
   TH1F * hmass_bb_1stbump = new TH1F("hmass_bb_1stbump", "mass_bb_1stbump", 100, 0., 100.);
   TH1F * hmass_B_1stbump = new TH1F("hmass_B_1stbump", "mass_B_1stbump", 60, 0., 300.);
   TH1F * hmass_dijet_1stbump = new TH1F("hmass_dijet_1stbump", "mass_dijet_1stbump", 25, 0., 1000.);
   TH1F * hmass_mumuj_1stbump = new TH1F("hmass_mumuj_1stbump", "mass_mumuj_1stbump", 25, 0., 1000.);
   TH1F * hmass_dimuondijet_1stbump  
     = new TH1F("hmass_dimuondijet_1stbump", "mass_dimuondijet_1stbump", 25, 0., 1000.);
   TH1F * hpt_bjet_1stbump  = new TH1F("hpt_bjet_1stbump", "pt_bjet_1stbump", 20, 0., 200.);
   TH1F * heta_bjet_1stbump  = new TH1F("heta_bjet_1stbump", "eta_bjet_1stbump", 20, -5.0, 5.0);
   TH1F * hpt_fwdj_1stbump  = new TH1F("hpt_fwdj_1stbump", "pt_fwdj_1stbump", 20, 0., 200.);
   TH1F * heta_fwdj_1stbump  = new TH1F("heta_fwdj_1stbump", "eta_fwdj_1stbump", 20, -5.0, 5.0);
   TH1F * hptmumu_div_mmumu_signal  = new TH1F("hptmumu_div_mmumu_signal", "ptmumu_div_mmumu_signal", 50, 0., 10.);

   // histos for 2nd bump
   TH1F * hpt_1stmu_2ndbump = new TH1F( "hpt_1stmu_2ndbump", "pt_1stmu_2ndbump", 150, 0., 150.);
   TH1F * hpt_2ndmu_2ndbump = new TH1F( "hpt_2ndmu_2ndbump", "pt_2ndmu_2ndbump", 80, 0., 80.);
   TH1F * hpt_dimuon_2ndbump  = new TH1F("hpt_dimuon_2ndbump", "pt_dimuon_2ndbump", 20, 0., 200.);
   TH1F * hcostheta_2ndbump = new TH1F("hcostheta_2ndbump", "costheta_2ndbump", 10, -1.0, 1.0);
   TH1F * hmass_mumu_2ndbump = new TH1F("hmass_mumu_2ndbump", "mass_mumu_2ndbump", 100, 0., 100.);
   TH1F * hmass_bb_2ndbump = new TH1F("hmass_bb_2ndbump", "mass_bb_2ndbump", 100, 0., 100.);
   TH1F * hmass_B_2ndbump = new TH1F("hmass_B_2ndbump", "mass_B_2ndbump", 60, 0., 300.);
   TH1F * hmass_dijet_2ndbump = new TH1F("hmass_dijet_2ndbump", "mass_dijet_2ndbump", 25, 0., 1000.);
   TH1F * hmass_mumuj_2ndbump = new TH1F("hmass_mumuj_2ndbump", "mass_mumuj_2ndbump", 25, 0., 1000.);
   TH1F * hmass_dimuondijet_2ndbump  
     = new TH1F("hmass_dimuondijet_2ndbump", "mass_dimuondijet_2ndbump", 25, 0., 1000.);
   TH1F * hpt_bjet_2ndbump  = new TH1F("hpt_bjet_2ndbump", "pt_bjet_2ndbump", 20, 0., 200.);
   TH1F * heta_bjet_2ndbump  = new TH1F("heta_bjet_2ndbump", "eta_bjet_2ndbump", 20, -5.0, 5.0);
   TH1F * hpt_barj_2ndbump  = new TH1F("hpt_barj_2ndbump", "pt_barj_2ndbump", 20, 0., 200.);
   TH1F * heta_barj_2ndbump  = new TH1F("heta_barj_2ndbump", "eta_barj_2ndbump", 20, -5.0, 5.0);
   TH1F * hdphi_dimuondijet_2ndbump = new TH1F("hdphi_dimuondijet_2ndbump","dphi_dimuondijet_2ndbump", 17, 0., 3.4);

   //
   Double_t ptmucutmax = 25.;
   Double_t ptmucutmin = 25.;
   Double_t etamucut1  = 2.1;
   Double_t etamucut2  = 2.1;
   Double_t ptj_cut = 30.;
   Double_t etaj_cut = 2.4;
   Double_t dphi_mumu_dijet_2ndbump_cut = 2.5;
   //   Double_t costheta_cut = 0.4;
   Double_t costheta_cut = 1.0;

   Double_t n_total = 0.;
   Double_t n_muons = 0.;
   Double_t n_1stbump = 0;
   Double_t n_2ndbump = 0;

   Double_t w_total = 0.;
   Double_t w_muons = 0.;
   Double_t w_1stbump = 0;
   Double_t w_2ndbump = 0;

   Double_t Rmujetmatch = 0.5;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      // prepare for cos calculation in the dimuon center mass
      // 1st muon is mu- (13), second muon is mu+ (-13)

      EtaJ->clear();
      PhiJ->clear();
      pTJ->clear();
      bjet->clear();

      Int_t njets = pTJ04->size();
      for( Int_t i = 0; i < njets; i++) {
	EtaJ->push_back((*EtaJ04)[i]);		
	PhiJ->push_back((*PhiJ04)[i]);		
	pTJ->push_back((*pTJ04)[i]);		
	bjet->push_back((*bjet04)[i]);		
      }

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


      TLorentzVector b1st, b2nd, dibquark;

      b1st.SetPx(ptb1*cos(phib1));
      b1st.SetPy(ptb1*sin(phib1));
      theta1 = 2. * atan(exp(-etab1));
      b1st.SetPz(ptb1/tan(theta1));
      b1st.SetE(ptb1/sin(theta1));

      b2nd.SetPx(ptb2*cos(phib2));
      b2nd.SetPy(ptb2*sin(phib2));
      theta2 = 2. * atan(exp(-etab2));
      b2nd.SetPz(ptb2/tan(theta2));
      b2nd.SetE(ptb2/sin(theta2));
      
      dibquark = b1st + b2nd;

      TLorentzVector mumubb = dimuon+dibquark;

      hmass_mumubb->Fill(mumubb.M());
      
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

      Double_t event_weight = 1.;

      hcostheta->Fill(costheta);

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
      
      double drmumu = deltaR(etamu1, phimu1, etamu2, phimu2);
      double drbb = deltaR(b1st.Eta(),b1st.Phi(),b2nd.Eta(),b2nd.Phi());
      double drmumubb = deltaR(dimuon.Eta(),dimuon.Phi(),dibquark.Eta(),dibquark.Phi());

      hdr_mumu->Fill(drmumu);
      hdr_bb->Fill(drbb);
      hdr_mumu_bb->Fill(drmumubb);
      hpt_dimuon->Fill(dimuon.Pt());
      hpt_dibquark->Fill(dibquark.Pt());
      hpt_2ndmu->Fill(ptmumin);

      // isolation

      Double_t isomuon1 = sumpt_trk_mu1/ptmu1;
      Double_t isomuon2 = sumpt_trk_mu2/ptmu2;

      n_muons += 1;
      w_muons += event_weight;

      // parton jet selections
      njets = pTJ->size();

      Int_t nbj_barrel = 0;
      Int_t nbj_endcap = 0;

      Int_t nlj_barrel = 0;
      Int_t nlj_endcap = 0;

      Int_t nj_barrel  = 0;
      Int_t nj_endcap  = 0;

      // jet classification
      // b-jet in barrel with maximal pT
      Double_t pTmaxBJbarrel = -10.;
      Double_t etaBJbarrel(0.), phiBJbarrel(0.);
      // no tagged jet in barrel with maximal pT
      Double_t pTmaxLJbarrel = -10.;
      Double_t etaLJbarrel(0.), phiLJbarrel(0.);
      // no tagged jet in the whole acceptance
      Double_t pTmaxLJ = -10.;
      Double_t etaLJ(0.), phiLJ(0.);
      // jet in forward with maximal pT
      Double_t pTmaxJforward = -10.;
      Double_t etaJforward(0.), phiJforward(0.);

      TLorentzVector b_jet, ljet, fjet;
      TLorentzVector barjet;
      TLorentzVector b_jet1, b_jet2;  // needed in case of two b-jets in the 2nd bump
      TLorentzVector dijet_1stbump, dimuonj_1stbump, dimuonb_1stbump, dimuondijet_1stbump;
      TLorentzVector dijet_2ndbump, dimuonj_2ndbump, dimuonb_2ndbump, dimuondijet_2ndbump;

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

	//	if (DRjmu1 <= Rmujetmatch || DRjmu2 <= Rmujetmatch) {continue;}

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

	// light jet in the whole acceptance
	if(fabs((*EtaJ)[i]) < 4.7 && (*bjet)[i] == 0) 
	  {
	    if((*pTJ)[i] > pTmaxLJ) {
	      pTmaxLJ = (*pTJ)[i];
	      etaLJ   = (*EtaJ)[i];
	      phiLJ   = (*PhiJ)[i];
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


      if(nbj_barrel == 1) {
	TLorentzVector dimuonb = dimuon+b_jet;
	if(dimuonb.M() > 160. && dimuonb.M() < 180.) {
	  if(nlj_barrel+nlj_endcap >= 1) {heta_ljet->Fill(etaLJ,event_weight);}      
	}
      }
      // 1st bump selections
      if(nbj_barrel == 1  && nj_barrel == 1 && nj_endcap >= 1) {
	n_1stbump += 1;
	w_1stbump += event_weight;

	dimuonb_1stbump = dimuon+b_jet;
	dijet_1stbump = b_jet + fjet;
	dimuonj_1stbump = dimuon+fjet;
	dimuondijet_1stbump = dimuon+dijet_1stbump;

	hisomu1->Fill(isomuon1);
	hisomu2->Fill(isomuon2);

	hptmumu_div_mmumu_signal->Fill(dimuon.Pt()/dimuon.M());
	
	hdr_dimuon_1stbump->Fill(drmumu);

	if(fabs(costheta) < costheta_cut) {
	  hpt_1stmu_1stbump->Fill(ptmumax,event_weight);
	  hpt_2ndmu_1stbump->Fill(ptmumin,event_weight);
	  hpt_dimuon_1stbump->Fill(dimuon.Pt(),event_weight);
	  hcostheta_1stbump->Fill(costheta,event_weight);
	  hmass_mumu_1stbump->Fill(dimuon.M(),event_weight);
	  hmass_bb_1stbump->Fill(dibquark.M(),event_weight);
	  hmass_B_1stbump->Fill(dimuonb_1stbump.M(),event_weight);
	  hmass_dijet_1stbump->Fill(dijet_1stbump.M(),event_weight);
	  hmass_mumuj_1stbump->Fill(dimuonj_1stbump.M(),event_weight);
	  hmass_dimuondijet_1stbump->Fill(dimuondijet_1stbump.M(),event_weight);  
	  hpt_bjet_1stbump->Fill(b_jet.Pt(),event_weight);
	  heta_bjet_1stbump->Fill(b_jet.Eta(),event_weight);
	  hpt_fwdj_1stbump->Fill(fjet.Pt(),event_weight);
	  heta_fwdj_1stbump->Fill(fjet.Eta(),event_weight);
	}
      }

      // 2nd bump selections
      if( nbj_barrel >= 1 && nj_barrel == 2 && nj_endcap == 0) {
      //      if( nbj_barrel == 2 && nj_barrel == 2 && nj_endcap == 0) {

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
	dimuonb_2ndbump = dimuon+b_jet;
	dimuondijet_2ndbump = dimuon+dijet_2ndbump;

	Double_t Dphi_dimuon_dijet_2ndbump = deltaPhi(dimuon.Phi(),dijet_2ndbump.Phi());
	hdphi_dimuondijet_2ndbump->Fill(Dphi_dimuon_dijet_2ndbump,event_weight);
	if(Dphi_dimuon_dijet_2ndbump > dphi_mumu_dijet_2ndbump_cut) {

	  n_2ndbump += 1;
	  w_2ndbump += event_weight;
	
	  hdr_dimuon_2ndbump->Fill(drmumu);

	  hisomu1->Fill(isomuon1);
	  hisomu2->Fill(isomuon2);

	  if(fabs(costheta) < costheta_cut) {
	    hpt_1stmu_2ndbump->Fill(ptmumax,event_weight);
	    hpt_2ndmu_2ndbump->Fill(ptmumin,event_weight);
	    hpt_dimuon_2ndbump->Fill(dimuon.Pt(),event_weight);
	    hmass_mumu_2ndbump->Fill(dimuon.M(),event_weight);
	    hmass_bb_2ndbump->Fill(dibquark.M(),event_weight);
	    hmass_B_2ndbump->Fill(dimuonb_2ndbump.M(),event_weight);
	    hcostheta_2ndbump->Fill(costheta,event_weight);
	    hmass_dijet_2ndbump->Fill(dijet_2ndbump.M(),event_weight);
	    hmass_mumuj_2ndbump->Fill(dimuonj_2ndbump.M(),event_weight);
	    hmass_dimuondijet_2ndbump->Fill(dimuondijet_2ndbump.M(),event_weight);
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

   TFile efile("test.root","recreate");
   heta_ljet->Write();

   hisomu1->Write();
   hisomu2->Write();

   hdr_dimuon_1stbump->Write();
   hdr_dimuon_2ndbump->Write();

   hmass_mumubb->Write();

   hdr_mumu->Write();
   hdr_bb->Write();
   hdr_mumu_bb->Write();
   hpt_dimuon->Write();
   hpt_dibquark->Write();
   hpt_2ndmu->Write();
   hcostheta->Write();
  
   // 1st bump histos
   hpt_1stmu_1stbump->Write();
   hpt_2ndmu_1stbump->Write();
   hpt_dimuon_1stbump->Write();
   hcostheta_1stbump->Write();
   hmass_mumu_1stbump->Write();
   hmass_bb_1stbump->Write();
   hmass_B_1stbump->Write();
   hmass_dijet_1stbump->Write();
   hmass_mumuj_1stbump->Write();
   hmass_dimuondijet_1stbump->Write();
   hpt_bjet_1stbump->Write();
   heta_bjet_1stbump->Write();
   hpt_fwdj_1stbump->Write();
   heta_fwdj_1stbump->Write();
   hptmumu_div_mmumu_signal->Write();

   // 2nd bump histos
   hpt_1stmu_2ndbump->Write();
   hpt_2ndmu_2ndbump->Write();
   hpt_dimuon_2ndbump->Write();
   hcostheta_2ndbump->Write();
   hmass_mumu_2ndbump->Write();
   hmass_bb_2ndbump->Write();
   hmass_B_2ndbump->Write();
   hmass_dijet_2ndbump->Write();
   hmass_mumuj_2ndbump->Write();
   hmass_dimuondijet_2ndbump->Write();
   hpt_bjet_2ndbump->Write();
   heta_bjet_2ndbump->Write();
   hpt_barj_2ndbump->Write();
   heta_barj_2ndbump->Write();
   hdphi_dimuondijet_2ndbump->Write();
   //
   efile.Close();
}
