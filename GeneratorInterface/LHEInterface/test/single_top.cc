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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// parton jets
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
// MC info
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
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

class single_top : public edm::EDAnalyzer {
public:
  explicit single_top(const edm::ParameterSet&);
  ~single_top();


private:
  //      virtual void beginJob(const edm::EventSetup&) ;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  // output root file
  string fOutputFileName ;
  // names of modules, producing object collections
  edm::InputTag parton_or_get_jets_Scr; 
  // variables to store in ntpl
  double ptw, etaw, phiw, chw;
  double ptb, etab, phib;
  double ptq, etaq, phiq;
  double ptmu, etamu, phimu;
  double ptnu, etanu, phinu;
  std::vector<double> *EtaJ;
  std::vector<double> *PhiJ;
  std::vector<double> *pTJ;
  std::vector<double> *bjet;

  //  Int_t BIDS[12]={511,521,531,541, 5112,5122,5132,5142, 5212,5222,5232,5242}
  
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
//single_top::beginJob(const edm::EventSetup&)
single_top::beginJob()
{
  using namespace edm;

  //  EtaRaw       = new std::vector<double>();

  // creating a simple tree

  hOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;

  EtaJ         = new std::vector<double>();
  PhiJ         = new std::vector<double>();
  pTJ          = new std::vector<double>();
  bjet         = new std::vector<double>();

  t1 = new TTree("t1","analysis tree");

  t1->Branch("ptw",&ptw,"ptw/D");
  t1->Branch("etaw",&etaw,"etaw/D");
  t1->Branch("phiw",&phiw,"phiw/D");
  t1->Branch("chw",&chw,"chw/D");
  //
  t1->Branch("ptb",&ptb,"ptb/D");
  t1->Branch("etab",&etab,"etab/D");
  t1->Branch("phib",&phib,"phib/D");
  //
  t1->Branch("ptq",&ptq,"ptq/D");
  t1->Branch("etaq",&etaq,"etaq/D");
  t1->Branch("phiq",&phiq,"phiq/D");
  //
  t1->Branch("ptmu",&ptmu,"ptmu/D");
  t1->Branch("etamu",&etamu,"etamu/D");
  t1->Branch("phimu",&phimu,"phimu/D");
  //
  t1->Branch("ptnu",&ptnu,"ptnu/D");
  t1->Branch("etanu",&etanu,"etanu/D");
  t1->Branch("phinu",&phinu,"phinu/D");
  //
  t1->Branch("EtaJ","vector<double>",&EtaJ);
  t1->Branch("PhiJ","vector<double>",&PhiJ);
  t1->Branch("pTJ" ,"vector<double>",&pTJ);
  t1->Branch("bjet" ,"vector<double>",&bjet);

  return ;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
single_top::endJob() {

  hOutputFile->Write() ;
  hOutputFile->Close() ;
  
  return ;
}

//
// constructors and destructor
//
single_top::single_top(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  using namespace edm;
  // 
  // get name of output file with histogramms
  fOutputFileName = iConfig.getUntrackedParameter<string>("HistOutFile");
  //
  // get parton jets
  parton_or_get_jets_Scr   = iConfig.getParameter<edm::InputTag>("parton_or_gen_jets");
}


single_top::~single_top()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
single_top::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<GenEventInfoProduct> genEvt;
  iEvent.getByLabel("generator",genEvt);

  ptw = 0.;
  etaw = 0.;
  phiw = 0.;
  chw  = 0.;

  ptb = 0.;
  etab = 0.;
  phib = 0.;

  ptq = 0.;
  etaq = 0.;
  phiq = 0.;

  ptmu   = 0.; 
  etamu  = 0.; 
  phimu  = 0.;
  ptnu   = 0.; 
  etanu  = 0.; 
  phinu  = 0.;

  EtaJ->clear();
  PhiJ->clear();
  pTJ->clear();
  bjet->clear();


  // parton jets

  edm::Handle<GenJetCollection> parton_or_gen_jets;
  iEvent.getByLabel(parton_or_get_jets_Scr, parton_or_gen_jets);

  if (!parton_or_gen_jets.isValid()) {return;}

  if( parton_or_gen_jets->size() != 0) {

    for(GenJetCollection::const_iterator parton_or_gen_jet = parton_or_gen_jets->begin();  
	parton_or_gen_jet != parton_or_gen_jets->end(); ++parton_or_gen_jet ) { 
      if( ( parton_or_gen_jet->pt() > 20. ) && ( fabs(parton_or_gen_jet->eta()) ) < 4.7 ) {

	/*
	double ib = 0;
	std::vector<const reco::Candidate*> parton_or_genparticles = parton_or_gen_jet->getJetConstituentsQuick();
	for (unsigned int i = 0; i < parton_or_genparticles.size(); ++i) {
	  const reco::Candidate* parton_or_genparticle = parton_or_genparticles[i];

//	  cout <<"   jet constituent i = " << i 
//	       <<" ID " << parton_or_genparticle->pdgId() 
//	       <<" px=" << parton_or_genparticle->px()
//	       <<" py=" << parton_or_genparticle->py()
//	       <<" pz=" << parton_or_genparticle->pz() << endl;

	  if( fabs(parton_or_genparticle->pdgId()) == 5 ) {
	    ib = 1.;
	  }
	}
	*/

	EtaJ->push_back(parton_or_gen_jet->eta());
	PhiJ->push_back(parton_or_gen_jet->phi());
	pTJ->push_back(parton_or_gen_jet->pt());
      }
    }
  }

  // define the jet flavour
  edm::Handle<reco::GenParticleCollection> genparticles;
  iEvent.getByLabel("genParticles", genparticles);

  for (unsigned int j = 0; j < pTJ->size(); ++j) {

    double ib = 0;

    //    cout <<"  take jet " << j << endl;

    for( size_t i = 0; i < genparticles->size(); i++) {
      const reco::GenParticle & p = (*genparticles)[i];
      if( abs(p.pdgId()) == 511 ||
	  abs(p.pdgId()) == 521 ||
	  abs(p.pdgId()) == 531 ||
	  abs(p.pdgId()) == 541 ||
	  abs(p.pdgId()) == 5112 ||	    
	  abs(p.pdgId()) == 5122 ||	    
	  abs(p.pdgId()) == 5132 ||	    
	  abs(p.pdgId()) == 5142 ||	    
	  abs(p.pdgId()) == 5212 ||	    
	  abs(p.pdgId()) == 5222 ||	    
	  abs(p.pdgId()) == 5232 ||	    
	  abs(p.pdgId()) == 5242 ) {
	double djb = deltaR((*EtaJ)[j],(*PhiJ)[j],p.eta(),p.phi());
	//	cout <<"  B hadron found for this jet " << i <<" ID = " << p.pdgId()
	//	     <<" eta = " << p.eta()
	//	     <<" phi = " << p.phi() << endl;
	if(djb < 0.5) {
	  ib = 1;
	  //	  cout <<"   -> jet below is matched with B hadron ID = " << p.pdgId() << endl;
	  continue;
	}
      }	    
    }
    bjet->push_back(ib);
    /*
    cout <<" jet " << j+1
	 <<" pT = " << (*pTJ)[j]
	 <<" eta = " << (*EtaJ)[j]
	 <<" phi = " << (*PhiJ)[j] 
	 <<" ib = " << (*bjet)[j] << endl;
    */
  }

  Double_t ptq_i=0.;
  Double_t etaq_i=0.;
  Double_t phiq_i=0.;
  //  Int_t qid_i=0;

  for( size_t i = 0; i < genparticles->size(); i++)
    {
      const reco::GenParticle & p = (*genparticles)[i];

      int mothid = 0;
      //      int mothst = 0;

      const Candidate * moth = p.mother();
      if(moth) {
	mothid = moth->pdgId();
	//	mothst = moth->status();
      }
      /*
      cout <<" i = " << i+1 
	   <<" ID = " << p.pdgId() 
	   <<" status = " << p.status() 
	   <<" eta = " << p.eta()
	   <<" mother id= " << mothid << endl; 
      */
      if( abs(p.pdgId()) < 5 && p.status() == 3) 
	{
	  ptq_i=p.pt();
	  etaq_i=p.eta();
	  phiq_i=p.phi();
	  //	  qid_i=p.pdgId();
	}
	//	   <<" mother stat= " << mothst << endl;

      // W from top decay and forward wuark
      if( abs(p.pdgId()) == 24 &&  abs(mothid) == 6 && p.status() == 3) {
	ptw = p.pt();
	etaw = p.eta();
	phiw = p.phi();
	chw  = p.charge();

	ptq = ptq_i;
	etaq = etaq_i;
	phiq = phiq_i;

	//	cout <<" ==> W ID = " << p.pdgId() <<" pTW = " << ptw <<" etaw = " << etaw << endl;
	//	cout <<" ==> q ID = " << qid_i <<" pTq = " << ptq <<" etaq = " << etaq << endl;
      }
      //  b from top decay
      if(abs(p.pdgId()) == 5 && abs(mothid) == 6 && p.status() == 3) {
	ptb = p.pt();
	etab = p.eta();
	phib = p.phi();
	//	cout <<" ==> b- ID = " << p.pdgId() <<" pTb = " << ptb <<" etab = " << etab << endl;
      }
      //  mu from W decay
      if(abs(p.pdgId()) == 13 && abs(mothid) == 24 && p.status() == 3) {
	ptmu = p.pt();
	etamu = p.eta();
	phimu = p.phi();
	//	cout <<" ==> mu ID = " << p.pdgId() <<" pTmu=" << ptmu <<" etamu=" << etamu <<" phimu=" << phimu << endl;
      }
      //  nu from W decay
      if(abs(p.pdgId()) == 14 && abs(mothid) == 24 && p.status() == 3) {
	ptnu = p.pt();
	etanu = p.eta();
	phinu = p.phi();
	//	cout <<" ==> nu ID = " << p.pdgId() <<" pTn =" << ptnu <<" etanu= " << etanu <<" phinu=" << phinu << endl;
      }
    }
  t1->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(single_top);
