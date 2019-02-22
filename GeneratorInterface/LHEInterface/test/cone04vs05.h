//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov  4 16:19:01 2016 by ROOT version 5.34/18
// from TTree t1/analysis tree
// found on file: cone04vs05.root
//////////////////////////////////////////////////////////

#ifndef cone04vs05_h
#define cone04vs05_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class cone04vs05 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        ptmu1;
   Double_t        etamu1;
   Double_t        phimu1;
   Double_t        ptmu2;
   Double_t        etamu2;
   Double_t        phimu2;
   Int_t           q_id;
   Double_t        q_px;
   Double_t        q_py;
   Double_t        q_pz;
   Double_t        q_en;
   Double_t        b_px;
   Double_t        b_py;
   Double_t        b_pz;
   Double_t        b_en;
   vector<double>  *EtaJ04;
   vector<double>  *PhiJ04;
   vector<double>  *pTJ04;
   vector<double>  *bjet04;
   vector<double>  *EtaJ05;
   vector<double>  *PhiJ05;
   vector<double>  *pTJ05;
   vector<double>  *bjet05;

   // List of branches
   TBranch        *b_ptmu1;   //!
   TBranch        *b_etamu1;   //!
   TBranch        *b_phimu1;   //!
   TBranch        *b_ptmu2;   //!
   TBranch        *b_etamu2;   //!
   TBranch        *b_phimu2;   //!
   TBranch        *b_q_id;   //!
   TBranch        *b_q_px;   //!
   TBranch        *b_q_py;   //!
   TBranch        *b_q_pz;   //!
   TBranch        *b_q_en;   //!
   TBranch        *b_b_px;   //!
   TBranch        *b_b_py;   //!
   TBranch        *b_b_pz;   //!
   TBranch        *b_b_en;   //!
   TBranch        *b_EtaJ04;   //!
   TBranch        *b_PhiJ04;   //!
   TBranch        *b_pTJ04;   //!
   TBranch        *b_bjet04;   //!
   TBranch        *b_EtaJ05;   //!
   TBranch        *b_PhiJ05;   //!
   TBranch        *b_pTJ05;   //!
   TBranch        *b_bjet05;   //!

   cone04vs05(TTree *tree=0);
   virtual ~cone04vs05();
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
};

#endif

#ifdef cone04vs05_cxx
cone04vs05::cone04vs05(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("oldsample.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("oldsample.root");
      }
      f->GetObject("t1",tree);

   }
   Init(tree);
}

cone04vs05::~cone04vs05()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t cone04vs05::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t cone04vs05::LoadTree(Long64_t entry)
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

void cone04vs05::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   EtaJ04 = 0;
   PhiJ04 = 0;
   pTJ04 = 0;
   bjet04 = 0;
   EtaJ05 = 0;
   PhiJ05 = 0;
   pTJ05 = 0;
   bjet05 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ptmu1", &ptmu1, &b_ptmu1);
   fChain->SetBranchAddress("etamu1", &etamu1, &b_etamu1);
   fChain->SetBranchAddress("phimu1", &phimu1, &b_phimu1);
   fChain->SetBranchAddress("ptmu2", &ptmu2, &b_ptmu2);
   fChain->SetBranchAddress("etamu2", &etamu2, &b_etamu2);
   fChain->SetBranchAddress("phimu2", &phimu2, &b_phimu2);
   fChain->SetBranchAddress("q_id", &q_id, &b_q_id);
   fChain->SetBranchAddress("q_px", &q_px, &b_q_px);
   fChain->SetBranchAddress("q_py", &q_py, &b_q_py);
   fChain->SetBranchAddress("q_pz", &q_pz, &b_q_pz);
   fChain->SetBranchAddress("q_en", &q_en, &b_q_en);
   fChain->SetBranchAddress("b_px", &b_px, &b_b_px);
   fChain->SetBranchAddress("b_py", &b_py, &b_b_py);
   fChain->SetBranchAddress("b_pz", &b_pz, &b_b_pz);
   fChain->SetBranchAddress("b_en", &b_en, &b_b_en);
   fChain->SetBranchAddress("EtaJ04", &EtaJ04, &b_EtaJ04);
   fChain->SetBranchAddress("PhiJ04", &PhiJ04, &b_PhiJ04);
   fChain->SetBranchAddress("pTJ04", &pTJ04, &b_pTJ04);
   fChain->SetBranchAddress("bjet04", &bjet04, &b_bjet04);
   fChain->SetBranchAddress("EtaJ05", &EtaJ05, &b_EtaJ05);
   fChain->SetBranchAddress("PhiJ05", &PhiJ05, &b_PhiJ05);
   fChain->SetBranchAddress("pTJ05", &pTJ05, &b_pTJ05);
   fChain->SetBranchAddress("bjet05", &bjet05, &b_bjet05);
   Notify();
}

Bool_t cone04vs05::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void cone04vs05::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t cone04vs05::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Double_t cone04vs05::deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi; 
  return dphi;
}

Double_t cone04vs05::deltaEta(Double_t eta1, Double_t eta2)
{
  Double_t deta = fabs(eta1 - eta2);
  return deta;
}

Double_t cone04vs05::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t deta = deltaEta(eta1, eta2);
  Double_t dphi = deltaPhi(phi1, phi2);
  Double_t dr = sqrt(deta*deta + dphi*dphi);
  return dr;
}

#endif // #ifdef cone04vs05_cxx
