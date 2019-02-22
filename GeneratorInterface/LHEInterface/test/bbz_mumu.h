//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 27 10:01:27 2014 by ROOT version 5.34/07
// from TTree t1/analysis tree
// found on file: bbz_mumu_py6_30gev_yb2.root
//////////////////////////////////////////////////////////

#ifndef bbz_mumu_h
#define bbz_mumu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class bbz_mumu {
public :
  // weights stuff
  Double_t pthweights[31];
  Double_t pthlowedge[31];
  Double_t binwidth[31];
  Int_t nbins;

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        pth22;
   Double_t        pth44;
   Double_t        pth62;
   Double_t        Hmass;
   Double_t        ptb1;
   Double_t        etab1;
   Double_t        phib1;
   Double_t        ptb2;
   Double_t        etab2;
   Double_t        phib2;
   Double_t        ptmu1;
   Double_t        etamu1;
   Double_t        phimu1;
   Double_t        ptmu2;
   Double_t        etamu2;
   Double_t        phimu2;
   vector<double>  *EtaJ;
   vector<double>  *PhiJ;
   vector<double>  *pTJ;
   vector<double>  *bjet;
   vector<double>   *evtweight;

   // List of branches
   TBranch        *b_pth22;   //!
   TBranch        *b_pth44;   //!
   TBranch        *b_pth62;   //!
   TBranch        *b_Hmass;   //!
   TBranch        *b_ptb1;   //!
   TBranch        *b_etab1;   //!
   TBranch        *b_phib1;   //!
   TBranch        *b_ptb2;   //!
   TBranch        *b_etab2;   //!
   TBranch        *b_phib2;   //!
   TBranch        *b_ptmu1;   //!
   TBranch        *b_etamu1;   //!
   TBranch        *b_phimu1;   //!
   TBranch        *b_ptmu2;   //!
   TBranch        *b_etamu2;   //!
   TBranch        *b_phimu2;   //!
   TBranch        *b_EtaJ;   //!
   TBranch        *b_PhiJ;   //!
   TBranch        *b_pTJ;   //!
   TBranch        *b_bjet;   //!
   TBranch        *b_evtweight;   //!

   bbz_mumu(TTree *tree=0);
   virtual ~bbz_mumu();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
  virtual void Weights();

};

#endif

#ifdef bbz_mumu_cxx
bbz_mumu::bbz_mumu(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

     //     TChain * chain = new TChain("t1","bbz_mumu_via_gbb.root");
     //     chain->Add("bbz_mumu_via_gbb.root");

     //     TChain * chain = new TChain("t1","bbz_mumu_via_gg.root");
     //     chain->Add("bbz_mumu_via_gg.root");

     //     TChain * chain = new TChain("t1","bbz_mumu_alphadef025_pdfrwt.root");
     //     chain->Add("bbz_mumu_alphadef025_pdfrwt.root");

     TChain * chain = new TChain("t1","bbz_mumu_alphadef_x_2_pdfrwt.root");
     chain->Add("bbz_mumu_alphadef_x_2_pdfrwt.root");

     tree = chain;
   }
   Init(tree);
}

bbz_mumu::~bbz_mumu()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t bbz_mumu::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t bbz_mumu::LoadTree(Long64_t entry)
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

void bbz_mumu::Init(TTree *tree)
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

   fChain->SetBranchAddress("pth22", &pth22, &b_pth22);
   fChain->SetBranchAddress("pth44", &pth44, &b_pth44);
   fChain->SetBranchAddress("pth62", &pth62, &b_pth62);
   fChain->SetBranchAddress("Hmass", &Hmass, &b_Hmass);
   fChain->SetBranchAddress("ptb1", &ptb1, &b_ptb1);
   fChain->SetBranchAddress("etab1", &etab1, &b_etab1);
   fChain->SetBranchAddress("phib1", &phib1, &b_phib1);
   fChain->SetBranchAddress("ptb2", &ptb2, &b_ptb2);
   fChain->SetBranchAddress("etab2", &etab2, &b_etab2);
   fChain->SetBranchAddress("phib2", &phib2, &b_phib2);
   fChain->SetBranchAddress("ptmu1", &ptmu1, &b_ptmu1);
   fChain->SetBranchAddress("etamu1", &etamu1, &b_etamu1);
   fChain->SetBranchAddress("phimu1", &phimu1, &b_phimu1);
   fChain->SetBranchAddress("ptmu2", &ptmu2, &b_ptmu2);
   fChain->SetBranchAddress("etamu2", &etamu2, &b_etamu2);
   fChain->SetBranchAddress("phimu2", &phimu2, &b_phimu2);
   fChain->SetBranchAddress("EtaJ", &EtaJ, &b_EtaJ);
   fChain->SetBranchAddress("PhiJ", &PhiJ, &b_PhiJ);
   fChain->SetBranchAddress("pTJ", &pTJ, &b_pTJ);
   fChain->SetBranchAddress("bjet", &bjet, &b_bjet);
   fChain->SetBranchAddress("evtweight", &evtweight, &b_evtweight);
   Notify();
}

Bool_t bbz_mumu::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void bbz_mumu::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t bbz_mumu::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


#endif // #ifdef bbz_mumu_cxx
