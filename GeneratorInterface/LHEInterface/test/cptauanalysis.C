#define cptauanalysis_cxx
#include "cptauanalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

Double_t cptauanalysis::deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi; 
  return dphi;
}

Double_t cptauanalysis::deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t cptauanalysis::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t deta = deltaEta(eta1, eta2);
  Double_t dphi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

void cptauanalysis::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L cptauanalysis.C
//      Root > cptauanalysis t
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

   // CP sensitive variables
   TH1F * hphip  = new TH1F( "hphip", "Dphi for y+y- > 0", 10, 0., 3.14);
   TH1F * hphim  = new TH1F( "hphim", "Dphi for y+y- < 0", 10, 0., 3.14);
   // pT pi- and photons
   TH1F * hptpipm  = new TH1F( "hptpipm", "pT pi+/- after selections", 50, 0., 100.);
   TH1F * hptgam  = new TH1F( "hptgam", "pT photons after selections", 50, 0., 100.);
   
   // plots
   TH1F * hdrpimpi01  = new TH1F( "hdrpimpi01", "DR(pi-,pi0), after selections", 25, 0., 0.435);
   TH1F * hdrgam1gam2_pi01  = new TH1F( "hdrgam1gam2_pi01", "DR(gam1-gam2) after selections", 25, 0., 0.435);
   TH1F * hdrmin_pim_gam  = new TH1F( "hdrmin_pim_gam", "DR min (pi-,gam), after selections", 25, 0., 0.435);

   // accceptance and off line cuts
   // pt of tau_h in off-line
   Double_t ptrocut = 40.;
   // eta acceptance for pi+/- and photons
   Double_t etacut = 1.479;
   // upper cut on Higgs pT
   Double_t pthuppercut = 1000.;
   // geometry and magnetic field
   Double_t RECAL = 1.29;
   Double_t B = 3.8;
   // pi 
   Double_t pi = 3.1415927;
   //
   //
   // counts
   Double_t ntot = 0.;
   Double_t nsel = 0.;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      ntot += 1.;

      TLorentzVector pip, pim, pi01, pi02, gamma1_pi01, gamma2_pi01, rop, rom, tworo;
      // pi+
      pip.SetPx(ptpip*cos(phipip));
      pip.SetPy(ptpip*sin(phipip));
      Double_t theta = 2.*atan(exp(-etapip));
      pip.SetPz(ptpip/tan(theta));
      pip.SetE(ptpip/sin(theta));
      // pi-
      pim.SetPx(ptpim*cos(phipim));
      pim.SetPy(ptpim*sin(phipim));
      theta = 2.*atan(exp(-etapim));
      pim.SetPz(ptpim/tan(theta));
      pim.SetE(ptpim/sin(theta));
      // pi01
      pi01.SetPx(ptpi01*cos(phipi01));
      pi01.SetPy(ptpi01*sin(phipi01));
      theta = 2.*atan(exp(-etapi01));
      pi01.SetPz(ptpi01/tan(theta));
      pi01.SetE(ptpi01/sin(theta));
      // pi02
      pi02.SetPx(ptpi02*cos(phipi02));
      pi02.SetPy(ptpi02*sin(phipi02));
      theta = 2.*atan(exp(-etapi02));
      pi02.SetPz(ptpi02/tan(theta));
      pi02.SetE(ptpi02/sin(theta));
      // gamma1 from pi01
      gamma1_pi01.SetPx(ptgam1pi01*cos(phigam1pi01));
      gamma1_pi01.SetPy(ptgam1pi01*sin(phigam1pi01));
      theta = 2.*atan(exp(-etagam1pi01));
      gamma1_pi01.SetPz(ptgam1pi01/tan(theta));
      gamma1_pi01.SetE(ptgam1pi01/sin(theta));
      // gamma2 from pi02
      gamma2_pi01.SetPx(ptgam2pi01*cos(phigam2pi01));
      gamma2_pi01.SetPy(ptgam2pi01*sin(phigam2pi01));
      theta = 2.*atan(exp(-etagam2pi01));
      gamma2_pi01.SetPz(ptgam2pi01/tan(theta));
      gamma2_pi01.SetE(ptgam2pi01/sin(theta));
      // variables to plot
      // DR pi- pi0
      Double_t DRpimpi01 = pim.DeltaR(pi01);
      // DR gam1-gam2 from pi01
      Double_t DRgam1gam2_pi01 =gamma1_pi01.DeltaR(gamma2_pi01);
      // Dphi shift of pi- in magnetic field when it reaches ECAL surface
      // and calculate a minimal distance betweej gammas from pi01 and pi-
      Double_t sinphi = 0.5*RECAL*0.3*B/ptpim;
      Double_t dphi_pim = asin(sinphi);
      Double_t new_phipim = phipim+dphi_pim;
      Double_t phipim_at_ECAL = new_phipim;
      if(new_phipim > pi) {phipim_at_ECAL = new_phipim-2.0*pi;}

      Double_t DR_pim_gamma1 = deltaR(etapim, phipim_at_ECAL, etagam1pi01, phigam1pi01);
      Double_t DR_pim_gamma2 = deltaR(etapim, phipim_at_ECAL, etagam2pi01, phigam2pi01);

      Double_t DRMin_pim_gamma = DR_pim_gamma1; 
      if(DR_pim_gamma1 > DR_pim_gamma2) {DRMin_pim_gamma = DR_pim_gamma2;}

      rom = pim + pi01;
      rop = pip + pi02;
      tworo = rom + rop;
      Double_t ym = (pim.E()-pi01.E())/(pim.E()+pi01.E());
      Double_t yp = (pip.E()-pi02.E())/(pip.E()+pi02.E());
      // two ro center frame
      Double_t bx = -tworo.Px()/tworo.E();
      Double_t by = -tworo.Py()/tworo.E();
      Double_t bz = -tworo.Pz()/tworo.E();
      // move pi+/- and pi0s to this system
      pim.Boost(bx,by,bz);
      pi01.Boost(bx,by,bz);
      pip.Boost(bx,by,bz);
      pi02.Boost(bx,by,bz);
      rom.Boost(bx,by,bz);
      rop.Boost(bx,by,bz);
      //
      Double_t cosromrop = (rom.Px()*rop.Px()+
			    rom.Py()*rop.Py()+
			    rom.Pz()*rop.Pz())/(rom.P()*rop.P());
      //      cout <<" cosromeom = " << cosromrop 
      //	   <<" p rom = " << rom.P()
      //	   <<" p rop = " << rom.P() << endl;
      // calculate vector product
      // prod  x     y     z
      // pim   x     y     z
      // pi0   x     y     z
      TLorentzVector romplane, ropplane;
      // for ro- plane
      romplane.SetPx(pim.Py()*pi01.Pz()-pim.Pz()*pi01.Py()); 
      romplane.SetPy(-1.0*(pim.Px()*pi01.Pz()-pim.Pz()*pi01.Px())); 
      romplane.SetPz(pim.Px()*pi01.Py()-pim.Py()*pi01.Px()); 
      // for ro+ plane
      ropplane.SetPx(pip.Py()*pi02.Pz()-pip.Pz()*pi02.Py()); 
      ropplane.SetPy(-1.0*(pip.Px()*pi02.Pz()-pip.Pz()*pi02.Px())); 
      ropplane.SetPz(pip.Px()*pi02.Py()-pip.Py()*pi02.Px()); 
      // define cos between two planes
      Double_t cosphi = (romplane.Px()*ropplane.Px()+
      			 romplane.Py()*ropplane.Py()+
      			 romplane.Pz()*ropplane.Pz())/(romplane.P()*ropplane.P());
      Double_t phi = acos(cosphi);
      //
      // selections
      if(pth62 > pthuppercut) {continue;}
      if(rom.Pt() < ptrocut || rop.Pt() < ptrocut) {continue;}
      if(fabs(etapip) > etacut || fabs(etapim) > etacut) {continue;}
      if(fabs(etagam1pi01) > etacut || fabs(etagam2pi01) > etacut) {continue;}    
      if(fabs(etagam1pi02) > etacut || fabs(etagam2pi02) > etacut) {continue;}    

      hptpipm->Fill(ptpim);
      hptpipm->Fill(ptpip);
      hptgam->Fill(ptgam1pi01);
      hptgam->Fill(ptgam2pi01);
      hptgam->Fill(ptgam1pi02);
      hptgam->Fill(ptgam2pi02);
      // fill histos for plots
      if(ym*yp > 0) {
	hphip->Fill(phi);
      }
      if(ym*yp < 0) {
	hphim->Fill(phi);
      }
      hdrpimpi01->Fill(DRpimpi01);
      hdrgam1gam2_pi01->Fill(DRgam1gam2_pi01);
      hdrmin_pim_gam->Fill(DRMin_pim_gamma); 

      nsel += 1.;
   }

   TCanvas* c1 = new TCanvas("X","Y",1);
   hphip->SetMaximum(1000.);
   hphip->SetMinimum(0.);
   hphip->Draw("hist");
   c1->SaveAs("hphip.pdf");

   TCanvas* c2 = new TCanvas("X","Y",1);
   hphim->SetMaximum(1000.);
   hphim->SetMinimum(0.);
   hphim->Draw("hist");
   c2->SaveAs("hphim.pdf");

   TCanvas* c3 = new TCanvas("X","Y",1);
   hdrgam1gam2_pi01->Draw("hist");
   c3->SaveAs("hdrgam1gam2.pdf");

   TCanvas* c4 = new TCanvas("X","Y",1);
   hdrpimpi01->Draw("hist");
   c4->SaveAs("hdrpimpi01.pdf");

   TCanvas* c5 = new TCanvas("X","Y",1);
   hdrmin_pim_gam->Draw("hist");
   c5->SaveAs("hdrmin_pim_gam.pdf");
   
   cout <<"  N events processed " << ntot << endl;
   cout << " " << endl;
   cout <<"  ========= Selections ===========" << endl;
   cout <<"    pT Higgs < " << pthuppercut << endl;
   cout <<"     pT ro1, ro2 > " << ptrocut << endl;
   cout <<"     eta of pi+/- and photons < " << etacut << endl;
   cout << " " << endl;
   cout <<"  N selected events = " << nsel << endl;
}
