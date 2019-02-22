#define mumuj_cxx
#include "mumuj.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

Double_t mumuj::deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi; 
  return dphi;
}

Double_t mumuj::deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t mumuj::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t deta = deltaEta(eta1, eta2);
  Double_t dphi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

void mumuj::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L mumuj.C
//      Root > mumuj t
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

   TH1F * hpth62  = new TH1F( "hpth62", "pth62", 20, 0., 200.);
   TH1F * hnljets = new TH1F("hnljets", "nljets", 5, 0., 5.);

   TH1F * hmass_mumu = new TH1F("hmass_mumu", "mass_mumu", 200, 0., 200.);
   TH2F * hnj_bf = new TH2F("hnj_bf", "nj_bf", 5, 0., 5., 5, 0., 5.);
   TH1F * hptj_div_ptz_barrel = new TH1F("hptj_div_ptz_barrel", "ptj_div_ptz_barrel", 20, 0., 2.);
   TH1F * hdphi_jz_barrel     = new TH1F("hdphi_jz_barrel", "hdphi_jz_barrel", 64, 0., 3.2);
   TH1F * hptj_div_ptz_forward = new TH1F("hptj_div_ptz_forward", "ptj_div_ptz_forward", 20, 0., 2.);
   TH1F * hdphi_jz_forward     = new TH1F("hdphi_jz_forward", "hdphi_jz_forward", 64, 0., 3.2);

   Double_t ptmucutmax = 25.;
   Double_t ptmucutmin = 20.;
   Double_t etamucut1  = 2.1;
   Double_t etamucut2  = 2.1;
   Double_t ptjcut = 20.0;
   Double_t etajcut = 4.7;
   Double_t etajbarrel = 2.5;
   Double_t etajforward = 3.2;
   Double_t dimuonmasslowcut  = 80.;
   Double_t dimuonmasshighcut = 100.;
   Double_t Dphi_dimuon_jet_cut = 2.6;
   Double_t Rmujetmatch = 0.5;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
   
   // if (Cut(ientry) < 0) continue;

      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // weight to fill histograms
      // MG5
      Double_t weight_evt = (*evtweight)[0];

      Double_t ptmumax = ptmu1;
      Double_t etamumax = etamu1;
      Double_t phimumax = phimu1;
      Double_t ptmumin = ptmu2;
      Double_t etamumin = etamu2;
      Double_t phimumin = phimu2;

      if(ptmu1 < ptmu2 ){

	ptmumax = ptmu2;
	etamumax = etamu2;
	phimumax = phimu2;

	ptmumin = ptmu1;
	etamumin = etamu1;
	phimumin = phimu1;

      }

      TLorentzVector mu1st, mu2nd, dimuon;
      mu1st.SetPx(ptmumax*cos(phimumax));
      mu1st.SetPy(ptmumax*sin(phimumax));
      Double_t theta = 2. * atan(exp(-etamumax));
      mu1st.SetPz(ptmumax/tan(theta));
      mu1st.SetE(ptmumax/sin(theta));

      mu2nd.SetPx(ptmumin*cos(phimumin));
      mu2nd.SetPy(ptmumin*sin(phimumin));
      theta = 2. * atan(exp(-etamumin));
      mu2nd.SetPz(ptmumin/tan(theta));
      mu2nd.SetE(ptmumin/sin(theta));

      dimuon = mu1st + mu2nd;
      Double_t dimuonmass = dimuon.M();

      // selections
      if(ptmumax<ptmucutmax || ptmumin<ptmucutmin || fabs(etamu1)>etamucut1 || fabs(etamu2)>etamucut2) {continue;}
      if(dimuonmass < dimuonmasslowcut || dimuonmass > dimuonmasshighcut) {continue;}

      // jets
      Int_t njet = pTJ->size();
      Double_t pTmaxLJ = -10.;
      Double_t etaLJ(0.), phiLJ(0.);
      // count barrel and forward jets
      Double_t nj_barrel = 0.;
      Double_t nj_forward = 0.;
      Int_t nljets = 0;

      for( Int_t i = 0; i < njet; i++) {

	if((*pTJ)[i] < ptjcut) {continue;}
	if(fabs((*EtaJ)[i]) > etajcut) {continue;}
 
	Double_t DRjmu1 = deltaR((*EtaJ)[i], (*PhiJ)[i], etamu1, phimu1);
	Double_t DRjmu2 = deltaR((*EtaJ)[i], (*PhiJ)[i], etamu2, phimu2);
	if (DRjmu1 < Rmujetmatch || DRjmu2 < Rmujetmatch) {continue;}
	nljets += 1;
	if((*pTJ)[i] > pTmaxLJ) {
	  pTmaxLJ = (*pTJ)[i];
	  etaLJ   = (*EtaJ)[i];
	  phiLJ   = (*PhiJ)[i];
	}
	Double_t Dphi_dimuon_jetc = deltaPhi(dimuon.Phi(),(*PhiJ)[i]);
	if(Dphi_dimuon_jetc > Dphi_dimuon_jet_cut) {
	  if(fabs((*EtaJ)[i]) < etajbarrel) {
	    nj_barrel += 1;
	  }
	  if(fabs((*EtaJ)[i]) > etajforward) {
	    nj_forward += 1;
	  }
	}
      }

      hnljets->Fill(nljets,weight_evt);

      TLorentzVector ljet;
      if(nljets < 1) {continue;}
      ljet.SetPx(pTmaxLJ*cos(phiLJ));
      ljet.SetPy(pTmaxLJ*sin(phiLJ));
      theta = 2. * atan(exp(-etaLJ));
      ljet.SetPz(pTmaxLJ/tan(theta));
      ljet.SetE(pTmaxLJ/sin(theta));
    
      hpth62->Fill(pth62,weight_evt);

      Double_t Dphi_dimuon_jet = deltaPhi(dimuon.Phi(),ljet.Phi());
      if(Dphi_dimuon_jet < Dphi_dimuon_jet_cut) {continue;}

    	hmass_mumu->Fill(dimuon.M(),weight_evt);

      //      if(dimuon.Pt() > 100. && dimuon.Pt() < 150. && nljets > 0) {
      //      if(dimuon.Pt() > 60. && dimuon.Pt() < 100. && nljets > 0) {
      if(dimuon.Pt() > 40. && dimuon.Pt() < 60.) {

	hnj_bf->Fill(nj_barrel,nj_forward,weight_evt);

	Double_t ptj_div_ptz = ljet.Pt()/ dimuon.Pt();
	if(fabs(ljet.Eta()) < etajbarrel) {
	  hptj_div_ptz_barrel->Fill(ptj_div_ptz,weight_evt);
	}
	if(fabs(ljet.Eta()) > etajforward) {
	  hptj_div_ptz_forward->Fill(ptj_div_ptz,weight_evt);
	}
      }

   }

   TFile efile("mumuj_mg5_amcnlo_histos.root","recreate");

   hpth62->Write();  
   hnljets->Write();

   hmass_mumu->Write();
   hnj_bf->Write();
   hptj_div_ptz_barrel->Write();
   hdphi_jz_barrel->Write();
   hptj_div_ptz_forward->Write();
   hdphi_jz_forward->Write();

   efile.Close();


   Double_t Rmujetmatch = 0.5;

   std::cout <<"             8 TeV Z+1j NLO " << std::endl;
   std::cout <<"        - pT 1st muon > " << ptmucutmax <<" eta < " << etamucut1 << std::endl; 
   std::cout <<"        - pT 2nd muon > " << ptmucutmin <<" eta < " << etamucut2 << std::endl; 
   std::cout <<"        - dimuon mass ["<< dimuonmasslowcut <<"-" << dimuonmasshighcut<<"] GeV" << std::endl;
   std::cout <<"        - DR(mu-jet) > " << Rmujetmatch << std::endl;
   std::cout <<"        - at least one jet pT >" << ptjcut <<" eta < " << etajcut <<  std::endl; 
   std::cout <<"        - barrel jets, eta <" << etajbarrel <<  std::endl; 
   std::cout <<"        - forward jets: " << etajforward <<" < eta < " << etajcut << std::endl; 
   std::cout <<"        - Dphi_dimuon_jet >" << Dphi_dimuon_jet_cut<<  std::endl; 
}
