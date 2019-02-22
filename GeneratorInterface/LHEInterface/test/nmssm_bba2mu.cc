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

class nmssm_bba2mu : public edm::EDAnalyzer {
public:
  explicit nmssm_bba2mu(const edm::ParameterSet&);
  ~nmssm_bba2mu();


private:
  //      virtual void beginJob(const edm::EventSetup&) ;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  // output root file
  string fOutputFileName ;
  // names of modules, producing object collections
  edm::InputTag partonjetsSrc; 
  // variables to store in ntpl
  double pth22, pth44, pth62, Hmass;
  double ptb1, etab1, phib1;
  double ptb2, etab2, phib2;
  double ptmu1, etamu1, phimu1;
  double ptmu2, etamu2, phimu2;
  double evtweight;

  std::vector<double> *EtaJ;
  std::vector<double> *PhiJ;
  std::vector<double> *pTJ;
  std::vector<double> *bjet;

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
//nmssm_bba2mu::beginJob(const edm::EventSetup&)
nmssm_bba2mu::beginJob()
{
  using namespace edm;

  // creating a simple tree

  hOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;

  EtaJ         = new std::vector<double>();
  PhiJ         = new std::vector<double>();
  pTJ          = new std::vector<double>();
  bjet         = new std::vector<double>();

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
  //  t1->Branch("EtaRaw","vector<double>",&EtaRaw);

  t1->Branch("evtweight",&evtweight,"evtweight/D");

  t1->Branch("EtaJ","vector<double>",&EtaJ);
  t1->Branch("PhiJ","vector<double>",&PhiJ);
  t1->Branch("pTJ" ,"vector<double>",&pTJ);
  t1->Branch("bjet" ,"vector<double>",&bjet);

  return ;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
nmssm_bba2mu::endJob() {

  hOutputFile->Write() ;
  hOutputFile->Close() ;

  return ;
}

//
// constructors and destructor
//
nmssm_bba2mu::nmssm_bba2mu(const edm::ParameterSet& iConfig)

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


nmssm_bba2mu::~nmssm_bba2mu()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
nmssm_bba2mu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

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

  evtweight = 1.;

  EtaJ->clear();
  PhiJ->clear();
  pTJ->clear();
  bjet->clear();

  edm::Handle<GenEventInfoProduct> genEvt;;
  iEvent.getByLabel("generator",genEvt);
 // parton jets
  edm::Handle<GenJetCollection> partonjets;
  iEvent.getByLabel(partonjetsSrc, partonjets);

  if (!partonjets.isValid()) {return;}

  if( partonjets->size() != 0) {

    for(GenJetCollection::const_iterator partonjet = partonjets->begin(); partonjet != partonjets->end(); ++partonjet ) { 
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
      }
    }
  }


  //  int procid = (int) genEvt->signalProcessID();

  //  edm::Handle<std::vector<PileupSummaryInfo> > infoPU;
  //  iEvent.getByLabel("addPileupInfo",infoPU);
  //  for(std::vector<PileupSummaryInfo>::const_iterator it = infoPU->begin(); it != infoPU->end(); it++)
  //    {
  //      genPU = it->getPU_NumInteractions();
  //      if (verbosity > 0) {   std::cout<<" Pileup "<<genPU<<std::endl;}
  //    }

  //  edm::Handle<reco::CandidateView> mctruth;
  //  for (unsigned int i = 0; i < mctruth->size(); ++ i) {
  //    if(ii > 29) break;
  //    const reco::Candidate& p = (*mctruth)[i];
  //  }
  


  //  cout <<"  Event particles " << endl;

  edm::Handle<reco::GenParticleCollection> genparticles;
  iEvent.getByLabel("genParticles", genparticles);

  for( size_t i = 0; i < genparticles->size(); i++)
    {
      const reco::GenParticle & p = (*genparticles)[i];

      int motherID = 0;
      int motherSt = 0;

      /*
      int grmotherID = 0;
      int grmotherSt = 0;

      int grgrmotherID = 0;
      int grgrmotherSt = 0;
      double grgrmotherpx = 0.;
      */

      const Candidate * moth = p.mother();
      if(moth) {
	motherID = moth->pdgId();
	motherSt = moth->status();	
	/*
	const Candidate * grmoth = moth->mother();
	if(grmoth) {
	  grmotherID = grmoth->pdgId();
	  grmotherSt = grmoth->status();
	  const Candidate * grgrmoth = grmoth->mother();
	  if(grgrmoth) {
	    grgrmotherID = grgrmoth->pdgId();
	    grgrmotherSt = grgrmoth->status();
	    grgrmotherpx = grgrmoth->px();
	  }	  
	}
	*/
      }

      /*
      cout <<" i " << i 
	   <<" ID " << p.pdgId() 
	   <<" status " << p.status() 
	   <<" motherID " << motherID 
	   <<" motherST " << motherSt
	   <<" grmotherID " << grmotherID 
	   <<" grmotherST " << grmotherSt 
	   <<" grgrmotherID " << grgrmotherID 
	   <<" grgrmotherST " << grgrmotherSt
	   <<" grgrmotherpx " << grgrmotherpx << endl;
      */

      // select Higgs
      if(p.pdgId() == 36) {
	if(p.status() == 3) {pth44 = p.pt();}
	if(p.status() == 2) {pth62 = p.pt();}
	Hmass = p.mass();
      }

      // select muon1
      if(p.pdgId() == 13 && p.status() == 3 &&
	 abs(motherID) == 36 && motherSt == 3) {
	ptmu1  = p.pt();
	etamu1 = p.eta();
	phimu1 = p.phi();
	//	cout <<"   --> muon1 " << i << endl;
      }
      // select muon2
      if(p.pdgId() == -13 && p.status() == 3 &&
	 abs(motherID) == 36 && motherSt == 3) {
	ptmu2  = p.pt();
	etamu2 = p.eta();
	phimu2 = p.phi();
	//	cout <<"   --> muon2 " << i << endl;
      }


      // select b quarks
      if(p.pdgId() == 5 && p.status() == 3) {
	ptb1  = p.pt();
	etab1 = p.eta();
	phib1 = p.phi();
	//	cout <<"   --> b-quark " << i <<" ptx = " << p.px() << endl;
      }

      // select b_bar quarks
      if(p.pdgId() == -5 && p.status() == 3) {
	ptb2  = p.pt();
	etab2 = p.eta();
	phib2 = p.phi();
	//	cout <<"   --> unti-b-quark " << i <<" ptx = " << p.px() << endl;
      }
    }

  t1->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(nmssm_bba2mu);
