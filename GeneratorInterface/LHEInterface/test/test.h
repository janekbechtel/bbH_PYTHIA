//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 23 13:02:47 2018 by ROOT version 5.34/18
// from TTree t1/analysis tree
// found on file: CPtestMG5_H_LO.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class test {
public :
  
  // jets pT > 30 GeV, |eta| < 4.7
  vector<double> *EtaJet;
  vector<double> *PhiJet;
  vector<double> *pTJet;

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        pth22;
   Double_t        pth44;
   Double_t        pth62;
   Double_t        Hmass;
   Double_t        ptmu1;
   Double_t        etamu1;
   Double_t        phimu1;
   Double_t        ptmu2;
   Double_t        etamu2;
   Double_t        phimu2;
   vector<double>  *EtaJ;
   vector<double>  *PhiJ;
   vector<double>  *pTJ;
   vector<double>  *etap;
   vector<double>  *phip;
   vector<double>  *ptp;
   vector<double>  *idp;
   vector<double>  *evtweight;

   // List of branches
   TBranch        *b_pth22;   //!
   TBranch        *b_pth44;   //!
   TBranch        *b_pth62;   //!
   TBranch        *b_Hmass;   //!
   TBranch        *b_ptmu1;   //!
   TBranch        *b_etamu1;   //!
   TBranch        *b_phimu1;   //!
   TBranch        *b_ptmu2;   //!
   TBranch        *b_etamu2;   //!
   TBranch        *b_phimu2;   //!
   TBranch        *b_EtaJ;   //!
   TBranch        *b_PhiJ;   //!
   TBranch        *b_pTJ;   //!
   TBranch        *b_etap;   //!
   TBranch        *b_phip;   //!
   TBranch        *b_ptp;   //!
   TBranch        *b_idp;   //!
   TBranch        *b_evtweight;   //!

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void MyVariables();
};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("CPtestMG5_H_LO.root");
     //     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("CPtestMG5_H.root");
      if (!f || !f->IsOpen()) {
	f = new TFile("CPtestMG5_H_LO.root");
	//	f = new TFile("CPtestMG5_H.root");
      }
      f->GetObject("t1",tree);

   }
   Init(tree);
}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
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

void test::Init(TTree *tree)
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
   etap = 0;
   phip = 0;
   ptp = 0;
   idp = 0;
   evtweight = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("pth22", &pth22, &b_pth22);
   fChain->SetBranchAddress("pth44", &pth44, &b_pth44);
   fChain->SetBranchAddress("pth62", &pth62, &b_pth62);
   fChain->SetBranchAddress("Hmass", &Hmass, &b_Hmass);
   fChain->SetBranchAddress("ptmu1", &ptmu1, &b_ptmu1);
   fChain->SetBranchAddress("etamu1", &etamu1, &b_etamu1);
   fChain->SetBranchAddress("phimu1", &phimu1, &b_phimu1);
   fChain->SetBranchAddress("ptmu2", &ptmu2, &b_ptmu2);
   fChain->SetBranchAddress("etamu2", &etamu2, &b_etamu2);
   fChain->SetBranchAddress("phimu2", &phimu2, &b_phimu2);
   fChain->SetBranchAddress("EtaJ", &EtaJ, &b_EtaJ);
   fChain->SetBranchAddress("PhiJ", &PhiJ, &b_PhiJ);
   fChain->SetBranchAddress("pTJ", &pTJ, &b_pTJ);
   fChain->SetBranchAddress("etap", &etap, &b_etap);
   fChain->SetBranchAddress("phip", &phip, &b_phip);
   fChain->SetBranchAddress("ptp", &ptp, &b_ptp);
   fChain->SetBranchAddress("idp", &idp, &b_idp);
   fChain->SetBranchAddress("evtweight", &evtweight, &b_evtweight);
   Notify();
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void test::MyVariables()
{
  EtaJet         = new std::vector<double>();
  PhiJet         = new std::vector<double>();
  pTJet          = new std::vector<double>();
}
#endif // #ifdef test_cxx
