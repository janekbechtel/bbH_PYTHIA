#define cone04vs05_cxx
#include "cone04vs05.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void cone04vs05::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L cone04vs05.C
//      Root > cone04vs05 t
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

   TH1F * hjet04VS05 = new TH1F("hjet04VS05", "jet04VS05", 40, 0.8, 1.2);

   //
   Double_t ptmucutmax = 25.;
   Double_t ptmucutmin = 25.;
   Double_t etamucut1  = 2.1;
   Double_t etamucut2  = 2.1;
   Double_t ptj_cut = 30.;
   Double_t etaj_cut = 2.4;
   Double_t dphi_mumu_dijet_2ndbump_cut = 2.5;

   Double_t n_total = 0.;
   Double_t n_muons = 0.;
   Double_t n_1stbump = 0;
   Double_t n_2ndbump = 0;

   Double_t Rmujetmatch = 0.5;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      // parton jet selections

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

      TLorentzVector b_jet, ljet, fjet;
      TLorentzVector barjet;
      TLorentzVector b_jet1, b_jet2;  // needed in case of two b-jets in the 2nd bump
      TLorentzVector dijet_1stbump, dimuonj_1stbump, dimuonb_1stbump, dimuondijet_1stbump;
      TLorentzVector dijet_2ndbump, dimuonj_2ndbump, dimuonb_2ndbump, dimuondijet_2ndbump;

      // parton jet selections
      Int_t njets = pTJ04->size();

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
      // no tagged jet in the whole acceptance
      Double_t pTmaxLJ = -10.;
      Double_t etaLJ(0.), phiLJ(0.);
      // jet in forward with maximal pT
      Double_t pTmaxJforward = -10.;
      Double_t etaJforward(0.), phiJforward(0.);

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
      // selections
      if(ptmumax<ptmucutmax || ptmumin<ptmucutmin || fabs(etamumax)>etamucut1 || fabs(etamumin)>etamucut2) 
	{continue;}

      n_muons += 1;

      for( Int_t i = 0; i < njets; i++) {

	if( (*pTJ04)[i] < ptj_cut) {continue;}

	Double_t DRjmu1 = deltaR((*EtaJ04)[i], (*PhiJ04)[i], etamu1, phimu1);
	Double_t DRjmu2 = deltaR((*EtaJ04)[i], (*PhiJ04)[i], etamu2, phimu2);
	Double_t theta;

	if (DRjmu1 <= Rmujetmatch || DRjmu2 <= Rmujetmatch) {continue;}

	if(fabs((*EtaJ04)[i]) < etaj_cut && (*bjet04)[i] == 1) 
	  {
	    nbj_barrel += 1;
	    if(nbj_barrel == 1) {
	      b_jet1.SetPx((*pTJ04)[i]*cos((*PhiJ04)[i]));
	      b_jet1.SetPy((*pTJ04)[i]*sin((*PhiJ04)[i]));
	      theta = 2. * atan(exp(-(*EtaJ04)[i]));
	      b_jet1.SetPz((*pTJ04)[i]/tan(theta));
	      b_jet1.SetE((*pTJ04)[i]/sin(theta));
	    }
	    if(nbj_barrel == 2) {
	      b_jet2.SetPx((*pTJ04)[i]*cos((*PhiJ04)[i]));
	      b_jet2.SetPy((*pTJ04)[i]*sin((*PhiJ04)[i]));
	      theta = 2. * atan(exp(-(*EtaJ04)[i]));
	      b_jet2.SetPz((*pTJ04)[i]/tan(theta));
	      b_jet2.SetE((*pTJ04)[i]/sin(theta));
	    }
	    if((*pTJ04)[i] > pTmaxBJbarrel) {
	      pTmaxBJbarrel = (*pTJ04)[i];
	      etaBJbarrel   = (*EtaJ04)[i];
	      phiBJbarrel   = (*PhiJ04)[i];
	    }
	  }

	if(fabs((*EtaJ04)[i]) < etaj_cut) {nj_barrel += 1;}
	
	if(fabs((*EtaJ04)[i]) > etaj_cut && fabs((*EtaJ04)[i]) < 4.7) 
	  {
	    nj_endcap += 1;
	    if((*pTJ04)[i] > pTmaxJforward) {
	      pTmaxJforward = (*pTJ04)[i];
	      etaJforward   = (*EtaJ04)[i];
	      phiJforward   = (*PhiJ04)[i];
	    }
	  }

	// count numbers of b-jets and light jets in barrel and endcap
	if(fabs((*EtaJ04)[i]) > etaj_cut && fabs((*EtaJ04)[i]) < 4.7 && (*bjet04)[i] == 1) {nbj_endcap += 1;}
	if(fabs((*EtaJ04)[i]) > etaj_cut && fabs((*EtaJ04)[i]) < 4.7 && (*bjet04)[i] == 0) {nlj_endcap += 1;}

	if(fabs((*EtaJ04)[i]) < etaj_cut && (*bjet04)[i] == 0) 
	  {
	    nlj_barrel += 1;
	    if((*pTJ04)[i] > pTmaxLJbarrel) {
	      pTmaxLJbarrel = (*pTJ04)[i];
	      etaLJbarrel   = (*EtaJ04)[i];
	      phiLJbarrel   = (*PhiJ04)[i];
	    }
	  }

	// light jet in the whole acceptance
	if(fabs((*EtaJ04)[i]) < 4.7 && (*bjet04)[i] == 0) 
	  {
	    if((*pTJ04)[i] > pTmaxLJ) {
	      pTmaxLJ = (*pTJ04)[i];
	      etaLJ   = (*EtaJ04)[i];
	      phiLJ   = (*PhiJ04)[i];
	    }
	  }

      }

      // 1st bump selections
      if(nbj_barrel == 1  && nj_barrel == 1 && nj_endcap >= 1) {
	n_1stbump += 1;
      }

      // 2nd bump selections
      if( nbj_barrel >= 1 && nj_barrel == 2 && nj_endcap == 0) {

	if(nbj_barrel == 1) {
	  barjet = ljet;
	}
	if(nbj_barrel == 2) {
	  barjet = b_jet2;
	}	
	dijet_2ndbump = b_jet + barjet;
	Double_t Dphi_dimuon_dijet_2ndbump = deltaPhi(dimuon.Phi(),dijet_2ndbump.Phi());
	if(Dphi_dimuon_dijet_2ndbump > dphi_mumu_dijet_2ndbump_cut) {
	  n_2ndbump += 1;
	}
      }
      //


      TLorentzVector qjet;

      qjet.SetPx(q_px);
      qjet.SetPy(q_py);
      qjet.SetPz(q_pz);
      qjet.SetE(q_en);
      
      /*
      cout <<" qjet px = " << qjet.Px() 
	   <<" py = " << qjet.Py() 
	   <<" pz = " << qjet.Pz() 
	   <<" pt = " << qjet.Pt() 
	   <<" eta = " << qjet.Eta() 
	   <<" phi = " << qjet.Phi() 
	   <<" id = " << q_id << endl;
      */
      // matching
      Double_t ptj04  = 0.;
      Double_t etaj04 = 0.;
      Double_t phij04 = 0.;
      Int_t njets04 = pTJ04->size();
      for( Int_t i = 0; i < njets04; i++) {
	Double_t DRjq = deltaR((*EtaJ04)[i], (*PhiJ04)[i], qjet.Eta(), qjet.Phi());
	/*
	cout <<"  jet04 " << i 
	     <<" pT = " << (*pTJ04)[i] 
	     <<" eta = " << (*EtaJ04)[i]
	     <<" phi = " << (*PhiJ04)[i] 
	     <<" DRjq = " << DRjq << endl;
	*/
	if(DRjq < 0.5) {
	  ptj04 = (*pTJ04)[i];
	  etaj04 = (*EtaJ04)[i];
	  phij04 = (*PhiJ04)[i];
	  break;
	  //	  cout <<"  this jet04 matches quark ! " << endl;
	}
      }

      Double_t ptj05  = 0.;
      Double_t etaj05 = 0.;
      Double_t phij05 = 0.;
      Int_t njets05 = pTJ05->size();
      for( Int_t i = 0; i < njets05; i++) {
	Double_t DRjq = deltaR((*EtaJ05)[i], (*PhiJ05)[i], qjet.Eta(), qjet.Phi());
	/*
	cout <<"  jet05 " << i 
	     <<" pT = " << (*pTJ05)[i] 
	     <<" eta = " << (*EtaJ05)[i]
	     <<" phi = " << (*PhiJ05)[i] 
	     <<" DRjq = " << DRjq << endl;
	*/
	if(DRjq < 0.5) {
	  ptj05  = (*pTJ05)[i];
	  etaj05 = (*EtaJ05)[i];
	  phij05 = (*PhiJ05)[i];
	  break;
	  //	  cout <<"  this jet05 matches quark ! " << endl;
	}
      }
      //      
      Double_t qjet_pt = qjet.Pt();
      if(qjet_pt > 40. && qjet_pt < 50000. && fabs(qjet.Eta()) > 2.4 && fabs(qjet.Eta()) < 4.7) { 
	if(ptj04 > 0. &&  ptj05 > 0.) {
	  hjet04VS05->Fill(ptj04/ptj05);
	}
      }
   }

   std::cout <<"  Total number of events = " << n_total << std::endl;
   std::cout <<" " << std::endl;
   std::cout <<"  Cuts on muons: " << ptmucutmax <<", "<< ptmucutmin <<"  GeV" 
	     <<" eta1 < "<< etamucut1 
	     <<" eta2 < "<< etamucut2 << std::endl;

   std::cout <<" " << std::endl;
   std::cout <<"  Cuts on jets, pt > " << ptj_cut <<" GeV, eta cut= "<< etaj_cut << std::endl; 
   std::cout <<" " << std::endl;
   std::cout <<"  Dphi dimuon-dijet > " << dphi_mumu_dijet_2ndbump_cut << std::endl;

   Double_t eff_muons = n_muons/n_total;
   Double_t eff_1stbump = n_1stbump/n_total;
   Double_t eff_2ndbump = n_2ndbump/n_total;

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

   TFile efile("test_histos.root","recreate");
   hjet04VS05->Write();
   efile.Close();
}
