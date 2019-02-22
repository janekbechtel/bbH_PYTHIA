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

class vlq : public edm::EDAnalyzer {
public:
  explicit vlq(const edm::ParameterSet&);
  ~vlq();


private:
  //      virtual void beginJob(const edm::EventSetup&) ;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  // output root file
  string fOutputFileName ;
  edm::InputTag parton_or_gen_jets_Src; 
  // names of modules, producing object collections
  //  edm::InputTag partonjetsSrc; 
  // variables to store in ntpl
  double ptmu1, etamu1, phimu1;
  double ptmu2, etamu2, phimu2;
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
//vlq::beginJob(const edm::EventSetup&)
vlq::beginJob()
{
  using namespace edm;

  // creating a simple tree

  hOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;

  EtaJ         = new std::vector<double>();
  PhiJ         = new std::vector<double>();
  pTJ          = new std::vector<double>();
  bjet         = new std::vector<double>();
  evtweight    = new std::vector<double>();

  t1 = new TTree("t1","analysis tree");

  //
  t1->Branch("ptmu1",&ptmu1,"ptmu1/D");
  t1->Branch("etamu1",&etamu1,"etamu1/D");
  t1->Branch("phimu1",&phimu1,"phimu1/D");
  //
  t1->Branch("ptmu2",&ptmu2,"ptmu2/D");
  t1->Branch("etamu2",&etamu2,"etamu2/D");
  t1->Branch("phimu2",&phimu2,"phimu2/D");
  //
  t1->Branch("EtaJ","vector<double>",&EtaJ);
  t1->Branch("PhiJ","vector<double>",&PhiJ);
  t1->Branch("pTJ" ,"vector<double>",&pTJ);
  t1->Branch("bjet" ,"vector<double>",&bjet);

  return ;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
vlq::endJob() {

  hOutputFile->Write() ;
  hOutputFile->Close() ;

  return ;
}

//
// constructors and destructor
//
vlq::vlq(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  using namespace edm;
  // 
  // get name of output file with histogramms
  fOutputFileName = iConfig.getUntrackedParameter<string>("HistOutFile");
  //
  // get parton jets
  parton_or_gen_jets_Src   = iConfig.getParameter<edm::InputTag>("parton_or_gen_jets");
}


vlq::~vlq()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
vlq::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

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

  //  cout <<" ==> event weight = " << evtweight << endl;

 // parton jets
  edm::Handle<GenJetCollection> parton_or_gen_jets;
  iEvent.getByLabel(parton_or_gen_jets_Src, parton_or_gen_jets);

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


  //  edm::Handle<reco::GenParticleCollection> genparticles;
  //  iEvent.getByLabel("genParticles", genparticles);

  //  math::XYZTLorentzVector  tauh(0.,0.,0.,0.);

  //  cout <<"  ==> genparticles size = " << genparticles->size() << endl;

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
      cout <<" i = " << i <<" id = " << p.pdgId() <<" status = " << p.status() 
	   <<" motherID = " << motherID <<" motherS = " << motherSt << endl;
      */

      if(p.pdgId() == 13 && (motherID == 7000001 || motherID == 7000021) && motherSt >= 50) {
	ptmu1  = p.pt();
	etamu1 = p.eta();
	phimu1 = p.phi();
	//	cout <<" ->mu1 px= " << p.px() << endl; 
      }

      if(p.pdgId() == -13 && (motherID == 7000001 || motherID == 7000021) && motherSt >= 50) {
	ptmu2  = p.pt();
	etamu2 = p.eta();
	phimu2 = p.phi();
	//	cout <<" ->mu2 px= " << p.px() << endl; 
      }
    }
  t1->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(vlq);
