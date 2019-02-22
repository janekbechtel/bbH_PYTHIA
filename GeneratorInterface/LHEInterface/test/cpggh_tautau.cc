// system include files
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
//
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/InputTag.h"
// gen particles
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/WeightsInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// parton jets
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
// MC info
#include "DataFormats/Math/interface/LorentzVector.h"
//#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "CLHEP/HepPDT/DefaultConfig.hh"
//
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "DataFormats/Math/interface/deltaR.h"
//double dR = deltaR( c1, c2 );
//
// root
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
//
using namespace std;
using namespace reco;

//
// class declaration
//

class cpggh_tautau : public edm::EDAnalyzer {
public:
  explicit cpggh_tautau(const edm::ParameterSet&);
  ~cpggh_tautau();


private:
  //      virtual void beginJob(const edm::EventSetup&) ;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  // output root file
  string fOutputFileName ;

  // variables to store in ntpl
  double pth22, pth44, pth62, Hmass;
  // pi+/-
  double ptpip, etapip, phipip;
  double ptpim, etapim, phipim;
  // pi0s and gammas
  double ptpi01, etapi01, phipi01;
  double ptgam1pi01, etagam1pi01, phigam1pi01;
  double ptgam2pi01, etagam2pi01, phigam2pi01;
  // 
  double ptpi02, etapi02, phipi02;
  double ptgam1pi02, etagam1pi02, phigam1pi02;
  double ptgam2pi02, etagam2pi02, phigam2pi02;
  // event weights
  std::vector<double> *evtweight;

  //
  TFile*      hOutputFile ;
  TTree*      t1;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

// ------------ method called once each job just before starting event loop  ------------
void 
//cpggh_tautau::beginJob(const edm::EventSetup&)
cpggh_tautau::beginJob()
{
  using namespace edm;

  // creating a simple tree

  hOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;

  evtweight    = new std::vector<double>();

  t1 = new TTree("t1","analysis tree");

  t1->Branch("pth22",&pth22,"pth22/D");
  t1->Branch("pth44",&pth44,"pth44/D");
  t1->Branch("pth62",&pth62,"pth62/D");
  t1->Branch("Hmass",&Hmass,"Hmass/D");

  // pi+
  t1->Branch("ptpip",&ptpip,"ptpip/D");
  t1->Branch("etapip",&etapip,"etapip/D");
  t1->Branch("phipip",&phipip,"phipip/D");
  // pi-
  t1->Branch("ptpim",&ptpim,"ptpim/D");
  t1->Branch("etapim",&etapim,"etapim/D");
  t1->Branch("phipim",&phipim,"phipim/D");
  // pi01 from ro+
  t1->Branch("ptpi01",&ptpi01,"ptpi01/D");
  t1->Branch("etapi01",&etapi01,"etapi01/D");
  t1->Branch("phipi01",&phipi01,"phipi01/D");
  // pi02 from ro-
  t1->Branch("ptpi02",&ptpi02,"ptpi02/D");
  t1->Branch("etapi02",&etapi02,"etapi02/D");
  t1->Branch("phipi02",&phipi02,"phipi02/D");
  // gamma1 from pi01 from ro+
  t1->Branch("ptgam1pi01",&ptgam1pi01,"ptgam1pi01/D");
  t1->Branch("etagam1pi01",&etagam1pi01,"etagam1pi01/D");
  t1->Branch("phigam1pi01",&phigam1pi01,"phigam1pi01/D");
  // gamma2 from pi01 from ro+
  t1->Branch("ptgam2pi01",&ptgam2pi01,"ptgam2pi01/D");
  t1->Branch("etagam2pi01",&etagam2pi01,"etagam2pi01/D");
  t1->Branch("phigam2pi01",&phigam2pi01,"phigam2pi01/D");
  // gamma1 from pi02 from ro-
  t1->Branch("ptgam1pi02",&ptgam1pi02,"ptgam1pi02/D");
  t1->Branch("etagam1pi02",&etagam1pi02,"etagam1pi02/D");
  t1->Branch("phigam1pi02",&phigam1pi02,"phigam1pi02/D");
  // gamma2 from pi02 from ro-
  t1->Branch("ptgam2pi02",&ptgam2pi02,"ptgam2pi02/D");
  t1->Branch("etagam2pi02",&etagam2pi02,"etagam2pi02/D");
  t1->Branch("phigam2pi02",&phigam2pi02,"phigam2pi02/D");
  //
  t1->Branch("evtweight","vector<double>",&evtweight);

  return ;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
cpggh_tautau::endJob() {

  hOutputFile->Write() ;
  hOutputFile->Close() ;

  return ;
}

//
// constructors and destructor
//
cpggh_tautau::cpggh_tautau(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  using namespace edm;
  // 
  // get name of output file with histogramms
  fOutputFileName = iConfig.getUntrackedParameter<string>("HistOutFile");
}


cpggh_tautau::~cpggh_tautau()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
cpggh_tautau::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  //
  pth22    = 0.; 
  pth44    = 0.; 
  pth62    = 0.; 
  Hmass    = 0.;
  // pi+
  ptpip = 0.; 
  etapip = 0.;
  phipip = 0.;
  // pi-
  ptpim = 0.;
  etapim = 0.;
  phipim = 0.;
  // pi0s and gammas
  ptpi01 = 0.; 
  etapi01 = 0.;
  phipi01 = 0.;
  ptgam1pi01 = 0.;
  etagam1pi01 = 0.;
  phigam1pi01 = 0.;
  ptgam2pi01 = 0.;
  etagam2pi01 = 0.;
  phigam2pi01 = 0.;
  // 
  ptpi02 = 0.;
  etapi02 = 0.;
  phipi02 = 0.;
  ptgam1pi02 = 0.;
  etagam1pi02 = 0.;
  phigam1pi02 = 0.;
  ptgam2pi02 = 0.;
  etagam2pi02 = 0.;
  phigam2pi02 = 0.;
  //
  evtweight->clear();

  edm::Handle<GenEventInfoProduct> genEvt;
  iEvent.getByLabel("generator",genEvt);

  // event weight
  //  double weightevt=genEvt->weight();
  //  cout <<"  weightevt = " << weightevt << endl;

  //
  //  for( unsigned int i = 0; i < 50; i++) {
  //    evtweight->push_back(1.);
  //  }

  edm::Handle<LHEEventProduct > lheEvt;
  iEvent.getByLabel("source",lheEvt);
  std::vector<gen::WeightsInfo> weights = lheEvt->weights();
  if(weights.size() != 0) {evtweight->clear();}
  for( unsigned int i = 0; i < weights.size(); i++) {
    gen::WeightsInfo lhew = weights[i];
    evtweight->push_back(lhew.wgt);
    //    cout <<"  weights -> i = " << lhew.id <<"  weight = " << lhew.wgt << endl;
  }

  /*
  iEvent.getByLabel("source",lheEvt);
  if(lheEvt.isValid()) {
  int nup = lheEvt->hepeup().NUP;
  for( int i = 0; i < nup; i++) {
  int lhep_id = lheEvt->hepeup().IDUP[i];
  lhef::HEPEUP::FiveVector p = lheEvt->hepeup().PUP[i];
  int motherl = lheEvt->hepeup().MOTHUP[i].first;
  int mother = 0.;
  if(motherl !=  0) {mother  = lheEvt->hepeup().IDUP[motherl-1];}
  cout <<" lhep_id " << lhep_id
  <<" motherl = " << motherl
  <<" motherID = " << mother << endl;
  double ppx = p[0];
  double ppy = p[1];
  double ppz = p[2];
  double ep  = p[3];
  if(i > 2) {	
  math::XYZTLorentzVector parton(ppx, ppy, ppz, ep);
  etap->push_back(parton.Eta());
  phip->push_back(parton.Phi());
  ptp->push_back(parton.Pt());
  idp->push_back(lhep_id);
  cout <<" id = " << lhep_id
  <<" px = " << ppx
  <<" py = " << ppy
  <<" pz = " << ppz
  <<" ep = " << ep << endl;
  }
  }
  }
  */

  edm::Handle<reco::GenParticleCollection> genparticles;
  iEvent.getByLabel("genParticles", genparticles);

  if (!genparticles.isValid()) {return;}

  for( size_t i = 0; i < genparticles->size(); i++)
    {
      const reco::GenParticle & p = (*genparticles)[i];
      
      int motherID = 0;
      //      int grmotherID = 0;
      //      int motherSt = 0;
      
      const Candidate * moth = p.mother();
      if(moth) {
	motherID = moth->pdgId();
	/*
	const Candidate * grmoth = moth->mother();
	if(grmoth) {
          grmotherID = grmoth->pdgId();
	}
	*/
      }
      
      /*
      cout <<" i = " << i <<" id = " << p.pdgId() 
	   <<" status = " << p.status() 
	   <<" mother ID = " << motherID
	   <<" grmotherID = " << grmotherID 
	   <<" ndouthers = " << p.numberOfDaughters() 
	   <<" pt/eta/phi = " << p.pt() <<" " << p.eta() <<" "<< p.phi() << endl;
      */
      //
      // ro-
      if(p.pdgId() == -213 && (abs(motherID) == 15)) {
	//	cout <<" i = " << i <<" id = " << p.pdgId() 
	//	     <<" status = " << p.status() 
	//	     <<" mother ID = " << motherID 
	//	     <<" ndouthers = " << p.numberOfDaughters() << endl;
	int nd = p.numberOfDaughters();
	if(nd != 0) {
	  for(int i = 0; i < nd; i++) {
	    const Candidate * daug = p.daughter(i);
	    int id =daug->pdgId();
	    //	    cout <<"  douther id = " << id << endl;
	    // pi-
	    if(id == -211) {
	      ptpim = daug->pt();
	      etapim = daug->eta();
	      phipim = daug->phi();
	      //	      cout <<" pi-, pt/eta/phi = " << ptpim <<" " << etapim <<" "<< phipim << endl;
	    }
	    if(id == 111) {
	      ptpi01 = daug->pt(); 
	      etapi01 = daug->eta();
	      phipi01 = daug->phi();
	      int nd = daug->numberOfDaughters();
	      if(nd != 0) {
		for(int i = 0; i < nd; i++) {
		  const Candidate * daug111 = daug->daughter(i);
		  //		  int id =daug111->pdgId();
		  //		  cout <<"      douther id = " << id << endl;
		  if(i == 0) {
		    ptgam1pi01 = daug111->pt();
		    etagam1pi01 = daug111->eta();
		    phigam1pi01 = daug111->phi();
		    //		    cout <<" gam1 from po-, pt/eta/phi = " << ptgam1pi01 
		    //			 <<" " << etagam1pi01
		    //			 <<" "<< phigam1pi01 << endl;
		  } else {
		    ptgam2pi01 = daug111->pt();
		    etagam2pi01 = daug111->eta();
		    phigam2pi01 = daug111->phi();
		    //		    cout <<" gam2 from po-, pt/eta/phi = " << ptgam2pi01 
		    //			 <<" " << etagam2pi01
		    //			 <<" "<< phigam2pi01 << endl;
		  }
		}
	      }
	    }
	  }
	}
      }

      // ro+
      if(p.pdgId() == 213 && (abs(motherID) == 15)) {
	//	cout <<" i = " << i <<" id = " << p.pdgId() 
	//	     <<" status = " << p.status() 
	//	     <<" mother ID = " << motherID 
	//	     <<" ndouthers = " << p.numberOfDaughters() << endl;
	int nd = p.numberOfDaughters();
	if(nd != 0) {
	  for(int i = 0; i < nd; i++) {
	    const Candidate * daug = p.daughter(i);
	    int id =daug->pdgId();
	    //	    cout <<"  douther id = " << id << endl;
	    // pi+
	    if(id == 211) {
	      ptpip = daug->pt();
	      etapip = daug->eta();
	      phipip = daug->phi();
	      //	      cout <<" pi+, pt/eta/phi = " << ptpip <<" " << etapip <<" "<< phipip << endl;
	    }
	    if(id == 111) {
	      ptpi02 = daug->pt(); 
	      etapi02 = daug->eta();
	      phipi02 = daug->phi();
	      int nd = daug->numberOfDaughters();
	      if(nd != 0) {
		for(int i = 0; i < nd; i++) {
		  const Candidate * daug111 = daug->daughter(i);
		  //		  int id =daug111->pdgId();
		  //		  cout <<"      douther id = " << id << endl;
		  if(i == 0) {
		    ptgam1pi02 = daug111->pt();
		    etagam1pi02 = daug111->eta();
		    phigam1pi02 = daug111->phi();
		    //		    cout <<" gam1 from po+, pt/eta/phi = " << ptgam1pi02 
		    //			 <<" " << etagam1pi02
		    //			 <<" "<< phigam1pi02 << endl;
		  } else {
		    ptgam2pi02 = daug111->pt();
		    etagam2pi02 = daug111->eta();
		    phigam2pi02 = daug111->phi();
		    //		    cout <<" gam2 from po+, pt/eta/phi = " << ptgam2pi02 
		    //			 <<" " << etagam2pi02
		    //			 <<" "<< phigam2pi02 << endl;
		  }
		}
	      }
	    }
	  }
	}
      }
      if( (p.pdgId() == 25 || p.pdgId() == 35) || p.pdgId() == 36 )
	{
	  //	  cout <<" Higgs id = " << p.pdgId() <<" status = " << fabs(p.status()) << endl;
	  if( fabs(p.status()) == 22) {pth22 = p.pt();}
	  if( fabs(p.status()) == 44) {pth44 = p.pt();}
	  if( fabs(p.status()) == 62) {
	    pth62 = p.pt();
	    Hmass = p.mass();
	  }
	}
    }
  t1->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(cpggh_tautau);
