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

class cone04vs05 : public edm::EDAnalyzer {
public:
  explicit cone04vs05(const edm::ParameterSet&);
  ~cone04vs05();


private:
  //      virtual void beginJob(const edm::EventSetup&) ;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  // output root file
  string fOutputFileName ;
  edm::InputTag genjets04_Src; 
  edm::InputTag genjets05_Src; 
  // names of modules, producing object collections
  //  edm::InputTag partonjetsSrc; 
  // variables to store in ntpl
  int q_id;
  double q_px, q_py, q_pz, q_en;
  double b_px, b_py, b_pz, b_en;
  double ptmu1, etamu1, phimu1;
  double ptmu2, etamu2, phimu2;
  double ptb1, etab1, phib1;
  double ptb2, etab2, phib2;

  double sumpt_trk_mu1, sumpt_trk_mu2;

  std::vector<double> *EtaJ04;
  std::vector<double> *PhiJ04;
  std::vector<double> *pTJ04;
  std::vector<double> *bjet04;

  std::vector<double> *EtaJ05;
  std::vector<double> *PhiJ05;
  std::vector<double> *pTJ05;
  std::vector<double> *bjet05;

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
//cone04vs05::beginJob(const edm::EventSetup&)
cone04vs05::beginJob()
{
  using namespace edm;

  // creating a simple tree

  hOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;

  EtaJ04         = new std::vector<double>();
  PhiJ04         = new std::vector<double>();
  pTJ04          = new std::vector<double>();
  bjet04         = new std::vector<double>();

  EtaJ05         = new std::vector<double>();
  PhiJ05         = new std::vector<double>();
  pTJ05          = new std::vector<double>();
  bjet05         = new std::vector<double>();

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
  t1->Branch("ptb1",&ptb1,"ptb1/D");
  t1->Branch("etab1",&etab1,"etab1/D");
  t1->Branch("phib1",&phib1,"phib1/D");
  //
  t1->Branch("ptb2",&ptb2,"ptb2/D");
  t1->Branch("etab2",&etab2,"etab2/D");
  t1->Branch("phib2",&phib2,"phib2/D");
  //
  t1->Branch("sumpt_trk_mu1",&sumpt_trk_mu1,"sumpt_trk_mu1/D");
  t1->Branch("sumpt_trk_mu2",&sumpt_trk_mu2,"sumpt_trk_mu2/D");
  //
  t1->Branch("q_id",&q_id,"q_id/I");
  t1->Branch("q_px",&q_px,"q_px/D");
  t1->Branch("q_py",&q_py,"q_py/D");
  t1->Branch("q_pz",&q_pz,"q_pz/D");
  t1->Branch("q_en",&q_en,"q_en/D");
  //
  t1->Branch("b_px",&b_px,"b_px/D");
  t1->Branch("b_py",&b_py,"b_py/D");
  t1->Branch("b_pz",&b_pz,"b_pz/D");
  t1->Branch("b_en",&b_en,"b_en/D");
  //
  t1->Branch("EtaJ04","vector<double>",&EtaJ04);
  t1->Branch("PhiJ04","vector<double>",&PhiJ04);
  t1->Branch("pTJ04" ,"vector<double>",&pTJ04);
  t1->Branch("bjet04" ,"vector<double>",&bjet04);
  //
  t1->Branch("EtaJ05","vector<double>",&EtaJ05);
  t1->Branch("PhiJ05","vector<double>",&PhiJ05);
  t1->Branch("pTJ05" ,"vector<double>",&pTJ05);
  t1->Branch("bjet05" ,"vector<double>",&bjet05);

  return ;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
cone04vs05::endJob() {

  hOutputFile->Write() ;
  hOutputFile->Close() ;

  return ;
}

//
// constructors and destructor
//
cone04vs05::cone04vs05(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  using namespace edm;
  // 
  // get name of output file with histogramms
  fOutputFileName = iConfig.getUntrackedParameter<string>("HistOutFile");
  //
  // get parton jets
  genjets04_Src   = iConfig.getParameter<edm::InputTag>("gen_jets_ak4");
  genjets05_Src   = iConfig.getParameter<edm::InputTag>("gen_jets_ak5");
}


cone04vs05::~cone04vs05()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
cone04vs05::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  ptmu1   = 0.; 
  etamu1  = 0.; 
  phimu1  = 0.;

  ptmu2   = 0.; 
  etamu2  = 0.; 
  phimu2  = 0.;

  ptb1   = 0.; 
  etab1  = 0.; 
  phib1  = 0.;

  ptb2   = 0.; 
  etab2  = 0.; 
  phib2  = 0.;

  sumpt_trk_mu1 = 0.;
  sumpt_trk_mu2 = 0.;

  q_id    = 0.; 
  q_px    = 0.; 
  q_py    = 0.;
  q_pz    = 0.;
  q_en    = 0.;

  b_px    = 0.; 
  b_py    = 0.;
  b_pz    = 0.;
  b_en    = 0.;

  EtaJ04->clear();
  PhiJ04->clear();
  pTJ04->clear();
  bjet04->clear();

  EtaJ05->clear();
  PhiJ05->clear();
  pTJ05->clear();
  bjet05->clear();

  evtweight->clear();

  edm::Handle<GenEventInfoProduct> genEvt;
  iEvent.getByLabel("generator",genEvt);

  // event weight
  //  double weightevt=genEvt->weight();

  edm::Handle<LHEEventProduct > lheEvt;
  iEvent.getByLabel("source",lheEvt);
  if(lheEvt.isValid()) {
    int nup = lheEvt->hepeup().NUP;
    for( int i = 0; i < nup; i++) {
      int lhep_id = lheEvt->hepeup().IDUP[i];
      lhef::HEPEUP::FiveVector p = lheEvt->hepeup().PUP[i];
      int motherl = lheEvt->hepeup().MOTHUP[i].first;
      int mother = 0.;
      if(motherl !=  0) {mother  = lheEvt->hepeup().IDUP[motherl-1];}
      //      cout <<" lhep_id " << lhep_id
      //	   <<" motherl = " << motherl
      //	   <<" motherID = " << mother << endl;
      double ppx = p[0];
      double ppy = p[1];
      double ppz = p[2];
      double ep  = p[3];
      if(abs(lhep_id) == 5 && abs(mother) == 6000007) {
	b_px = ppx;
	b_py = ppy;
	b_pz = ppz;
	b_en = ep;
	//	cout <<"  B fould " << endl;
      }
      if(i == (nup-1)) {
	q_px = ppx;
	q_py = ppy;
	q_pz = ppz;
	q_en = ep;
	q_id = lhep_id;
	//	cout <<"  q_if = " << q_id << endl;
      }
    }
    //
    /*
    if(noX == 0) {
      if(lheEvt.isValid()) {
	int nup = lheEvt->hepeup().NUP;
	for( int i = 0; i < nup; i++) {
	  int lhep_id = lheEvt->hepeup().IDUP[i];
	  cout <<" i = " << i <<" lhep_id = " << lhep_id << endl;
	}
      }
    }
    */
    //
    std::vector<gen::WeightsInfo> weights = lheEvt->weights();
    if(weights.size() != 0) {evtweight->clear();}
    for( unsigned int i = 0; i < weights.size(); i++) {
      gen::WeightsInfo lhew = weights[i];
      evtweight->push_back(lhew.wgt);
    }
  }

  //  cout <<" ==> event weight = " << evtweight << endl;

 // gen jets 04
  edm::Handle<GenJetCollection> genjets04;
  iEvent.getByLabel(genjets04_Src, genjets04);

  if (!genjets04.isValid()) {return;}

  if( genjets04->size() != 0) {

    for(GenJetCollection::const_iterator parton_or_gen_jet = genjets04->begin();  
	parton_or_gen_jet != genjets04->end(); ++parton_or_gen_jet ) { 
      if( ( parton_or_gen_jet->pt() > 20. ) && ( fabs(parton_or_gen_jet->eta()) ) < 4.7 ) {
	EtaJ04->push_back(parton_or_gen_jet->eta());
	PhiJ04->push_back(parton_or_gen_jet->phi());
	pTJ04->push_back(parton_or_gen_jet->pt());
      }
    }
  }

 // gen jets 05
  edm::Handle<GenJetCollection> genjets05;
  iEvent.getByLabel(genjets05_Src, genjets05);

  if (!genjets05.isValid()) {return;}

  if( genjets05->size() != 0) {

    for(GenJetCollection::const_iterator parton_or_gen_jet = genjets05->begin();  
	parton_or_gen_jet != genjets05->end(); ++parton_or_gen_jet ) { 
      if( ( parton_or_gen_jet->pt() > 20. ) && ( fabs(parton_or_gen_jet->eta()) ) < 4.7 ) {
	EtaJ05->push_back(parton_or_gen_jet->eta());
	PhiJ05->push_back(parton_or_gen_jet->phi());
	pTJ05->push_back(parton_or_gen_jet->pt());
      }
    }
  }


  edm::Handle<reco::GenParticleCollection> genparticles;
  iEvent.getByLabel("genParticles", genparticles);

  // define the jet 04 flavour
  for (unsigned int j = 0; j < pTJ04->size(); ++j) {
    double ib = 0;
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
	double djb = deltaR((*EtaJ04)[j],(*PhiJ04)[j],p.eta(),p.phi());
	if(djb < 0.5) {
	  ib = 1;
	  continue;
	}
      }	    
    }
    bjet04->push_back(ib);
  }

  // define the jet 05 flavour
  for (unsigned int j = 0; j < pTJ05->size(); ++j) {
    double ib = 0;
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
	double djb = deltaR((*EtaJ05)[j],(*PhiJ05)[j],p.eta(),p.phi());
	if(djb < 0.5) {
	  ib = 1;
	  continue;
	}
      }	    
    }
    bjet05->push_back(ib);
  }

  // select muons
  for( size_t i = 0; i < genparticles->size(); i++)
    {
      const reco::GenParticle & p = (*genparticles)[i];

      int motherID = 0;
      //      int motherSt = 0;

      const Candidate * moth = p.mother();
      if(moth) {
	motherID = moth->pdgId();
	//	motherSt = moth->status();	
      }

      
      if(p.pdgId() == 5  && (motherID == 23 || motherID == 35 || motherID == 36)) {
	ptb1  = p.pt();
	etab1 = p.eta();
	phib1 = p.phi();
	//	cout <<" ->b1 px= " << p.px() << endl; 
      }

      if(p.pdgId() == -5 && (motherID == 23 || motherID == 35 || motherID == 36)) {
	ptb2  = p.pt();
	etab2 = p.eta();
	phib2 = p.phi();
	//	cout <<" ->b2 px= " << p.px() << endl; 
      }

      if(p.pdgId() == 13  && (motherID == 23 || abs(motherID) == 7000021 || motherID == 35 || motherID == 36)) {
	ptmu1  = p.pt();
	etamu1 = p.eta();
	phimu1 = p.phi();
      }

      if(p.pdgId() == -13 && (motherID == 23 || abs(motherID) == 7000021 || motherID == 35 || motherID == 36)) {
	ptmu2  = p.pt();
	etamu2 = p.eta();
	phimu2 = p.phi();
      }
    }

  // calculate sum pT tracks, pT > 0.5 GeV around muon
  for( size_t i = 0; i < genparticles->size(); i++)
    {
      const reco::GenParticle & p = (*genparticles)[i];

      if(p.status() == 1 && p.charge() != 0 && abs(p.pdgId()) != 13 && p.pt() > 0.5) 
	{
	  double dr_trk_mu1 = deltaR(p.eta(),p.phi(),etamu1,phimu1);
	  double dr_trk_mu2 = deltaR(p.eta(),p.phi(),etamu2,phimu2);
	  if(dr_trk_mu1 < 0.3) {sumpt_trk_mu1 += p.pt();}
	  if(dr_trk_mu2 < 0.3) {sumpt_trk_mu2 += p.pt();}
	  /*
	  cout <<" i = " << i <<" status = " << p.status()
	       <<" charge = " <<  p.charge()
	       <<" ID = " << abs(p.pdgId())
	       <<" pT = " << p.pt()
	       <<" drmu1 = " << dr_trk_mu1
	       <<" drmu2 = " << dr_trk_mu2
	       <<" sumpt_trk_mu1 = " << sumpt_trk_mu1
	       <<" sumpt_trk_mu2 = " << sumpt_trk_mu2 << endl;
	  */
	}
    }

  t1->Fill();
}
//define this as a plug-in
DEFINE_FWK_MODULE(cone04vs05);
