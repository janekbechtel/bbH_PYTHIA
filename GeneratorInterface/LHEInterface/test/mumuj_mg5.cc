// system include files
// http://cmslxr.fnal.gov/lxr/source/RecoJets/JetAlgorithms/src/JetIDHelper.cc?v=CMSSW_3_3_6
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
//#include "CLHEP/Vector/LorentzVector.h"
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

class mumuj_mg5 : public edm::EDAnalyzer {
public:
  explicit mumuj_mg5(const edm::ParameterSet&);
  ~mumuj_mg5();


private:
  //      virtual void beginJob(const edm::EventSetup&) ;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  // output root file
  string fOutputFileName ;
  edm::InputTag partonjetsSrc; 
  // names of modules, producing object collections
  //  edm::InputTag partonjetsSrc; 
  // variables to store in ntpl
  double pth22, pth44, pth62, Hmass;
  double ptb1, etab1, phib1;
  double ptb2, etab2, phib2;
  double ptmu1, etamu1, phimu1;
  double ptmu2, etamu2, phimu2;  
  //  double evtweight;

  std::vector<double> *EtaJ;
  std::vector<double> *PhiJ;
  std::vector<double> *pTJ;
  std::vector<double> *bjet;
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
//mumuj_mg5::beginJob(const edm::EventSetup&)
mumuj_mg5::beginJob()
{
  using namespace edm;

  // creating a simple tree

  hOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;

  EtaJ         = new std::vector<double>();
  PhiJ         = new std::vector<double>();
  pTJ          = new std::vector<double>();
  bjet         = new std::vector<double>();
  bjet         = new std::vector<double>();
  evtweight    = new std::vector<double>();
  
  t1 = new TTree("t1","analysis tree");

  t1->Branch("pth22",&pth22,"pth22/D");
  t1->Branch("pth44",&pth44,"pth44/D");
  t1->Branch("pth62",&pth62,"pth62/D");
  t1->Branch("Hmass",&Hmass,"Hmass/D");

  t1->Branch("ptb1",&ptb1,"ptb1/D");
  t1->Branch("etab1",&etab1,"etab1/D");
  t1->Branch("phib1",&phib1,"phib1/D");
  //
  t1->Branch("ptb2",&ptb2,"ptb2/D");
  t1->Branch("etab2",&etab2,"etab2/D");
  t1->Branch("phib2",&phib2,"phib2/D");
  //
  t1->Branch("ptmu1",&ptmu1,"ptmu1/D");
  t1->Branch("etamu1",&etamu1,"etamu1/D");
  t1->Branch("phimu1",&phimu1,"phimu1/D");
  //
  t1->Branch("ptmu2",&ptmu2,"ptmu2/D");
  t1->Branch("etamu2",&etamu2,"etamu2/D");
  t1->Branch("phimu2",&phimu2,"phimu2/D");

  t1->Branch("EtaJ","vector<double>",&EtaJ);
  t1->Branch("PhiJ","vector<double>",&PhiJ);
  t1->Branch("pTJ" ,"vector<double>",&pTJ);
  t1->Branch("bjet" ,"vector<double>",&bjet);
  t1->Branch("evtweight","vector<double>",&evtweight);


  return ;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
mumuj_mg5::endJob() {

  hOutputFile->Write() ;
  hOutputFile->Close() ;

  return ;
}

//
// constructors and destructor
//
mumuj_mg5::mumuj_mg5(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  using namespace edm;
  // 
  // get name of output file with histogramms
  fOutputFileName = iConfig.getUntrackedParameter<string>("HistOutFile");
  //
  // get parton jets
  partonjetsSrc      = iConfig.getParameter<edm::InputTag>("parton_jets");
}


mumuj_mg5::~mumuj_mg5()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
mumuj_mg5::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  //
  pth22    = 0.; 
  pth44    = 0.; 
  pth62    = 0.; 
  Hmass    = 0.;

  ptb1    = 0.; 
  etab1   = 0.; 
  phib1   = 0.;
  ptb2    = 0.; 
  etab2   = 0.; 
  phib2   = 0.;
  ptmu1   = 0.; 
  etamu1  = 0.; 
  phimu1  = 0.;
  ptmu2   = 0.; 
  etamu2  = 0.; 
  phimu2  = 0.;

  EtaJ->clear();
  PhiJ->clear();
  pTJ->clear();
  bjet->clear();
  evtweight->clear();

  double Zmass = 0.;

  edm::Handle<GenEventInfoProduct> genEvt;
  iEvent.getByLabel("generator",genEvt);

  // event weight
  //  evtweight=genEvt->weight();
  //  cout <<" ==> event weight = " << evtweight << endl;

  edm::Handle<LHEEventProduct > lheEvt;
  iEvent.getByLabel("source",lheEvt);
  std::vector<gen::WeightsInfo> weights = lheEvt->weights();  
  for( unsigned int i = 0; i < weights.size(); i++) {
    gen::WeightsInfo lhew = weights[i];
    evtweight->push_back(lhew.wgt);
    //    cout <<"  weights -> i = " << lhew.id <<"  weight = " << lhew.wgt << endl;
  }
 // parton jets
  edm::Handle<GenJetCollection> partonjets;
  iEvent.getByLabel(partonjetsSrc, partonjets);

  if( partonjets->size() != 0) {

    for(GenJetCollection::const_iterator partonjet = partonjets->begin();  partonjet != partonjets->end(); ++partonjet ) { 
      if( ( partonjet->pt() > 30. ) && ( fabs(partonjet->eta()) ) < 4.7 ) {
	double ib = 0;
	std::vector<const reco::Candidate*> partons = partonjet->getJetConstituentsQuick();
	for (unsigned int i = 0; i < partons.size(); ++i) {
	  const reco::Candidate* parton = partons[i];
	  //	      cout <<"   jet constituent i = " << i <<" ID " << parton->pdgId() << endl;
	  if( fabs(parton->pdgId()) == 5 ) {
	    ib = 1.;
	  }
	}
	EtaJ->push_back(partonjet->eta());
	PhiJ->push_back(partonjet->phi());
	pTJ->push_back(partonjet->pt());
	bjet->push_back(ib);
	// cout <<"      jet in acceptance jet pt = " <<  partonjet->pt() <<" eta = " << partonjet->eta() <<" b " << ib << endl; 
      }
    }
  }

  edm::Handle<reco::GenParticleCollection> genparticles;
  iEvent.getByLabel("genParticles", genparticles);

  //  math::XYZTLorentzVector  tauh(0.,0.,0.,0.);

  //  cout <<"  ==> genparticles size = " << genparticles->size() << endl;

   math::XYZTLorentzVector  muon1(0.,0.,0.,0.);
   math::XYZTLorentzVector  muon2(0.,0.,0.,0.);

  for( size_t i = 0; i < genparticles->size(); i++)
    {
      const reco::GenParticle & p = (*genparticles)[i];

      int motherID = 0;
      int motherSt = 0;

      const Candidate * moth = p.mother();
      if(moth) {
	motherID = moth->pdgId();
	motherSt = moth->status();	
      }

      /*      
      if(fabs(p.pdgId()) == 13) {cout <<" i = " << i <<" id = " << p.pdgId() 
				      <<" status = " << p.status()
				      <<" mother = " << motherID 
				      <<" mother status = " << motherSt << endl;}
      */
      
      int imuplus = 0;

      if(((p.pdgId() == 13 || p.pdgId() == 14) && motherID == 23 && motherSt == 62) ||
	  (p.pdgId() == 13 && p.status() == 23) ) {
	imuplus += 1;
	math::XYZTLorentzVector muonc(p.px(), p.py(), p.pz(), p.p());
	if(imuplus == 1) {
	  ptmu1  = p.pt();
	  etamu1 = p.eta();
	  phimu1 = p.phi();
	  muon1 = muonc;
	  //	  cout <<" -> mu1 px= " << p.px() <<" status = " << p.status() << endl; 
	}
      }

      //      if(p.pdgId() == -13 && motherID == 23 && motherSt == 62) {
      int imuminus = 0;
      if(((p.pdgId() == -13 || p.pdgId() == -14) && motherID == 23 && motherSt == 62) ||
	  (p.pdgId() == -13 && p.status() == 23) ) {
	imuminus += 1;
	math::XYZTLorentzVector muonc(p.px(), p.py(), p.pz(), p.p());
	if(imuminus == 1) {
	  ptmu2  = p.pt();
	  etamu2 = p.eta();
	  phimu2 = p.phi();
	  muon2 = muonc;
	  //	  cout <<" -> mu2 px= " << p.px() <<" status = " << p.status() << endl; 
	}
      }

      // select b quarks
      int ibquark = 0;
      if(p.pdgId() == 5 && p.status() == 23) {
	ibquark += 1;
	if(ibquark == 1) {
	  ptb1  = p.pt();
	  etab1 = p.eta();
	  phib1 = p.phi();
	}
	//       	cout <<" -> b quark, pt = " << p.pt() <<" eta = " << p.eta() << endl; 
      }
      //
      int ib_barquark = 0;
      if(p.pdgId() == -5 && p.status() == 23) {
	ib_barquark += 1;
	if(ib_barquark == 1) {
	  ptb2  = p.pt();
	  etab2 = p.eta();
	  phib2 = p.phi();
	}
	//       	cout <<" -> b~quark, pt = " << p.pt() <<" eta = " << p.eta() << endl; 
	//	cout <<" ->b2 px= " << p.px() << endl; 
      }
      if(p.pdgId() == 23 && p.status() == 62) {
	Zmass = p.mass();
      }
    }
  math::XYZTLorentzVector twomuons = muon1 + muon2;
  pth62 = twomuons.Pt();
  Hmass = twomuons.M();
  if(Hmass < 10.) {
    cout <<"   ===> di muon mass = " << Hmass <<" Z mass " << Zmass << endl;
  }
  t1->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(mumuj_mg5);
