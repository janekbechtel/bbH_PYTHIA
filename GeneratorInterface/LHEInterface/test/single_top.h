//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 25 14:59:04 2016 by ROOT version 5.32/00
// from TTree t1/analysis tree
// found on file: single_top.root
//////////////////////////////////////////////////////////

#ifndef single_top_h
#define single_top_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class single_top {
public :

  // weights stuff
  Int_t nbins;
  Double_t weights[50];
  Double_t lowedge[50];
  Double_t width[50];

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        ptw;
   Double_t        etaw;
   Double_t        phiw;
   Double_t        chw;
   Double_t        ptb;
   Double_t        etab;
   Double_t        phib;
   Double_t        ptq;
   Double_t        etaq;
   Double_t        phiq;
   Double_t        ptmu;
   Double_t        etamu;
   Double_t        phimu;
   Double_t        ptnu;
   Double_t        etanu;
   Double_t        phinu;
   vector<double>  *EtaJ;
   vector<double>  *PhiJ;
   vector<double>  *pTJ;
   vector<double>  *bjet;

   // List of branches
   TBranch        *b_ptw;   //!
   TBranch        *b_etaw;   //!
   TBranch        *b_phiw;   //!
   TBranch        *b_chw;   //!
   TBranch        *b_ptb;   //!
   TBranch        *b_etab;   //!
   TBranch        *b_phib;   //!
   TBranch        *b_ptq;   //!
   TBranch        *b_etaq;   //!
   TBranch        *b_phiq;   //!
   TBranch        *b_ptmu;   //!
   TBranch        *b_etamu;   //!
   TBranch        *b_phimu;   //!
   TBranch        *b_ptnu;   //!
   TBranch        *b_etanu;   //!
   TBranch        *b_phinu;   //!
   TBranch        *b_EtaJ;   //!
   TBranch        *b_PhiJ;   //!
   TBranch        *b_pTJ;   //!
   TBranch        *b_bjet;   //!

   single_top(TTree *tree=0);
   virtual ~single_top();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

  virtual Double_t deltaEta(Double_t phi1, Double_t phi2);
  virtual Double_t deltaPhi(Double_t phi1, Double_t phi2);
  virtual Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
  virtual void Weights();
};

#endif

#ifdef single_top_cxx
single_top::single_top(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     //     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("single_top_genjets_mt150.root");
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("single_top_genjets.root");
      if (!f || !f->IsOpen()) {
	//	f = new TFile("single_top_genjets_mt150.root");
	f = new TFile("single_top_genjets.root");
      }
      f->GetObject("t1",tree);

   }
   Init(tree);
}

single_top::~single_top()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t single_top::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t single_top::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void single_top::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   EtaJ = 0;
   PhiJ = 0;
   pTJ = 0;
   bjet = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ptw", &ptw, &b_ptw);
   fChain->SetBranchAddress("etaw", &etaw, &b_etaw);
   fChain->SetBranchAddress("phiw", &phiw, &b_phiw);
   fChain->SetBranchAddress("chw", &chw, &b_chw);
   fChain->SetBranchAddress("ptb", &ptb, &b_ptb);
   fChain->SetBranchAddress("etab", &etab, &b_etab);
   fChain->SetBranchAddress("phib", &phib, &b_phib);
   fChain->SetBranchAddress("ptq", &ptq, &b_ptq);
   fChain->SetBranchAddress("etaq", &etaq, &b_etaq);
   fChain->SetBranchAddress("phiq", &phiq, &b_phiq);
   fChain->SetBranchAddress("ptmu", &ptmu, &b_ptmu);
   fChain->SetBranchAddress("etamu", &etamu, &b_etamu);
   fChain->SetBranchAddress("phimu", &phimu, &b_phimu);
   fChain->SetBranchAddress("ptnu", &ptnu, &b_ptnu);
   fChain->SetBranchAddress("etanu", &etanu, &b_etanu);
   fChain->SetBranchAddress("phinu", &phinu, &b_phinu);
   fChain->SetBranchAddress("EtaJ", &EtaJ, &b_EtaJ);
   fChain->SetBranchAddress("PhiJ", &PhiJ, &b_PhiJ);
   fChain->SetBranchAddress("pTJ", &pTJ, &b_pTJ);
   fChain->SetBranchAddress("bjet", &bjet, &b_bjet);
   Notify();
}

Bool_t single_top::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void single_top::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t single_top::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Double_t single_top::deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi; 
  return dphi;
}

Double_t single_top::deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t single_top::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t deta = deltaEta(eta1, eta2);
  Double_t dphi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

#endif // #ifdef single_top_cxx
