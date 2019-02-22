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

class cpggh : public edm::EDAnalyzer {
public:
  explicit cpggh(const edm::ParameterSet&);
  ~cpggh();


private:
  //      virtual void beginJob(const edm::EventSetup&) ;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  // output root file
  string fOutputFileName ;
  edm::InputTag partonjetsSrc; 

  // variables to store in ntpl
  double pth22, pth44, pth62, Hmass;
  double ptmu1, etamu1, phimu1;
  double ptmu2, etamu2, phimu2;  

  std::vector<double> *EtaJ;
  std::vector<double> *PhiJ;
  std::vector<double> *pTJ;

  std::vector<double> *etap;
  std::vector<double> *phip;
  std::vector<double> *ptp;
  std::vector<double> *idp;

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
//cpggh::beginJob(const edm::EventSetup&)
cpggh::beginJob()
{
  using namespace edm;

  // creating a simple tree

  hOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;

  EtaJ         = new std::vector<double>();
  PhiJ         = new std::vector<double>();
  pTJ          = new std::vector<double>();

  etap         = new std::vector<double>();
  phip         = new std::vector<double>();
  ptp          = new std::vector<double>();
  idp          = new std::vector<double>();

  evtweight    = new std::vector<double>();

  t1 = new TTree("t1","analysis tree");

  t1->Branch("pth22",&pth22,"pth22/D");
  t1->Branch("pth44",&pth44,"pth44/D");
  t1->Branch("pth62",&pth62,"pth62/D");
  t1->Branch("Hmass",&Hmass,"Hmass/D");

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

  t1->Branch("etap","vector<double>",&etap);
  t1->Branch("phip","vector<double>",&phip);
  t1->Branch("ptp" ,"vector<double>",&ptp);
  t1->Branch("idp" ,"vector<double>",&idp);

  t1->Branch("evtweight","vector<double>",&evtweight);

  return ;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
cpggh::endJob() {

  hOutputFile->Write() ;
  hOutputFile->Close() ;

  return ;
}

//
// constructors and destructor
//
cpggh::cpggh(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  using namespace edm;
  // 
  // get name of output file with histogramms
  fOutputFileName = iConfig.getUntrackedParameter<string>("HistOutFile");
  //
  // get parton jets
  partonjetsSrc      = iConfig.getParameter<edm::InputTag>("gen_jets_ak4");
}


cpggh::~cpggh()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
cpggh::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  //
  pth22    = 0.; 
  pth44    = 0.; 
  pth62    = 0.; 
  Hmass    = 0.;
  //
  ptmu1   = 0.; 
  etamu1  = 0.; 
  phimu1  = 0.;
  ptmu2   = 0.; 
  etamu2  = 0.; 
  phimu2  = 0.;
  //
  EtaJ->clear();
  PhiJ->clear();
  pTJ->clear();
  //
  etap->clear();
  phip->clear();
  ptp->clear();
  idp->clear();
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

  iEvent.getByLabel("source",lheEvt);
  if(lheEvt.isValid()) {
    int nup = lheEvt->hepeup().NUP;
    for( int i = 0; i < nup; i++) {
      int lhep_id = lheEvt->hepeup().IDUP[i];
      lhef::HEPEUP::FiveVector p = lheEvt->hepeup().PUP[i];
      //  int motherl = lheEvt->hepeup().MOTHUP[i].first;
      //      int mother = 0.;
      //      if(motherl !=  0) {mother  = lheEvt->hepeup().IDUP[motherl-1];}
      //      cout <<" lhep_id " << lhep_id
      //	   <<" motherl = " << motherl
      //	   <<" motherID = " << mother << endl;
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
	/*
	cout <<" id = " << lhep_id
	     <<" px = " << ppx
	     <<" py = " << ppy
	     <<" pz = " << ppz
	     <<" ep = " << ep << endl;
	*/
      }
    }
  }

  // parton jets
  edm::Handle<GenJetCollection> partonjets;
  iEvent.getByLabel(partonjetsSrc, partonjets);
  
  if (!partonjets.isValid()) {return;}
  
  if( partonjets->size() != 0) {

    for(GenJetCollection::const_iterator partonjet = partonjets->begin();  
	partonjet != partonjets->end(); ++partonjet ) { 

      if( ( partonjet->pt() > 20. ) && ( fabs(partonjet->eta()) ) < 4.7 ) {
	
	/*
	  double ip = 0;
	  std::vector<const reco::Candidate*> partons = partonjet->getJetConstituentsQuick();
	  for (unsigned int i = 0; i < partons.size(); ++i) {
	  const reco::Candidate* parton = partons[i];
	  cout <<"   jet constituent i = " << i 
	  <<" ID " << parton->pdgId() 
	  <<" px=" << parton->px()
	  <<" py=" << parton->py()
	  <<" pz=" << parton->pz() << endl;
	  
	  if( fabs(parton->pdgId()) != 5 ) {
	  ip = 1.;
	  }
	  }
	*/
	EtaJ->push_back(partonjet->eta());
	PhiJ->push_back(partonjet->phi());
	pTJ->push_back(partonjet->pt());
      }
    }
  }

  edm::Handle<reco::GenParticleCollection> genparticles;
  iEvent.getByLabel("genParticles", genparticles);

  if (!genparticles.isValid()) {return;}

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
      
      cout <<" i = " << i <<" id = " << p.pdgId() 
	   <<" status = " << p.status() 
	   <<" mother ID = " << motherID 
	   <<" ndouthers = " << p.numberOfDaughters() << endl;
	      

      if(p.pdgId() == 13 && (motherID == 25 || motherID == 36) && motherSt == 62) {
	ptmu1  = p.pt();
	etamu1 = p.eta();
	phimu1 = p.phi();
	//	cout <<" ->mu1 px= " << p.px() << endl; 
      }

      if(p.pdgId() == -13 && (motherID == 25 || motherID == 36) && motherSt == 62) {
	ptmu2  = p.pt();
	etamu2 = p.eta();
	phimu2 = p.phi();
	//	cout <<" ->mu2 px= " << p.px() << endl; 
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
DEFINE_FWK_MODULE(cpggh);
