//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 23 12:29:21 2018 by ROOT version 5.34/36
// from TTree t1/analysis tree
// found on file: cpgghntpl.root
//////////////////////////////////////////////////////////

#ifndef cptauanalysis_h
#define cptauanalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class cptauanalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        pth22;
   Double_t        pth44;
   Double_t        pth62;
   Double_t        Hmass;
   Double_t        ptpip;
   Double_t        etapip;
   Double_t        phipip;
   Double_t        ptpim;
   Double_t        etapim;
   Double_t        phipim;
   Double_t        ptpi01;
   Double_t        etapi01;
   Double_t        phipi01;
   Double_t        ptpi02;
   Double_t        etapi02;
   Double_t        phipi02;
   Double_t        ptgam1pi01;
   Double_t        etagam1pi01;
   Double_t        phigam1pi01;
   Double_t        ptgam2pi01;
   Double_t        etagam2pi01;
   Double_t        phigam2pi01;
   Double_t        ptgam1pi02;
   Double_t        etagam1pi02;
   Double_t        phigam1pi02;
   Double_t        ptgam2pi02;
   Double_t        etagam2pi02;
   Double_t        phigam2pi02;
   vector<double>  *evtweight;

   // List of branches
   TBranch        *b_pth22;   //!
   TBranch        *b_pth44;   //!
   TBranch        *b_pth62;   //!
   TBranch        *b_Hmass;   //!
   TBranch        *b_ptpip;   //!
   TBranch        *b_etapip;   //!
   TBranch        *b_phipip;   //!
   TBranch        *b_ptpim;   //!
   TBranch        *b_etapim;   //!
   TBranch        *b_phipim;   //!
   TBranch        *b_ptpi01;   //!
   TBranch        *b_etapi01;   //!
   TBranch        *b_phipi01;   //!
   TBranch        *b_ptpi02;   //!
   TBranch        *b_etapi02;   //!
   TBranch        *b_phipi02;   //!
   TBranch        *b_ptgam1pi01;   //!
   TBranch        *b_etagam1pi01;   //!
   TBranch        *b_phigam1pi01;   //!
   TBranch        *b_ptgam2pi01;   //!
   TBranch        *b_etagam2pi01;   //!
   TBranch        *b_phigam2pi01;   //!
   TBranch        *b_ptgam1pi02;   //!
   TBranch        *b_etagam1pi02;   //!
   TBranch        *b_phigam1pi02;   //!
   TBranch        *b_ptgam2pi02;   //!
   TBranch        *b_etagam2pi02;   //!
   TBranch        *b_phigam2pi02;   //!
   TBranch        *b_evtweight;   //!

   cptauanalysis(TTree *tree=0);
   virtual ~cptauanalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef cptauanalysis_cxx
cptauanalysis::cptauanalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("cpevengghntpl.root");
     //     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("cpoddgghntpl.root");
      if (!f || !f->IsOpen()) {
	f = new TFile("cpevengghntpl.root");
	//	f = new TFile("cpoddgghntpl.root");
      }
      f->GetObject("t1",tree);

   }
   Init(tree);
}

cptauanalysis::~cptauanalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t cptauanalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t cptauanalysis::LoadTree(Long64_t entry)
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

void cptauanalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
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
   fChain->SetBranchAddress("ptpip", &ptpip, &b_ptpip);
   fChain->SetBranchAddress("etapip", &etapip, &b_etapip);
   fChain->SetBranchAddress("phipip", &phipip, &b_phipip);
   fChain->SetBranchAddress("ptpim", &ptpim, &b_ptpim);
   fChain->SetBranchAddress("etapim", &etapim, &b_etapim);
   fChain->SetBranchAddress("phipim", &phipim, &b_phipim);
   fChain->SetBranchAddress("ptpi01", &ptpi01, &b_ptpi01);
   fChain->SetBranchAddress("etapi01", &etapi01, &b_etapi01);
   fChain->SetBranchAddress("phipi01", &phipi01, &b_phipi01);
   fChain->SetBranchAddress("ptpi02", &ptpi02, &b_ptpi02);
   fChain->SetBranchAddress("etapi02", &etapi02, &b_etapi02);
   fChain->SetBranchAddress("phipi02", &phipi02, &b_phipi02);
   fChain->SetBranchAddress("ptgam1pi01", &ptgam1pi01, &b_ptgam1pi01);
   fChain->SetBranchAddress("etagam1pi01", &etagam1pi01, &b_etagam1pi01);
   fChain->SetBranchAddress("phigam1pi01", &phigam1pi01, &b_phigam1pi01);
   fChain->SetBranchAddress("ptgam2pi01", &ptgam2pi01, &b_ptgam2pi01);
   fChain->SetBranchAddress("etagam2pi01", &etagam2pi01, &b_etagam2pi01);
   fChain->SetBranchAddress("phigam2pi01", &phigam2pi01, &b_phigam2pi01);
   fChain->SetBranchAddress("ptgam1pi02", &ptgam1pi02, &b_ptgam1pi02);
   fChain->SetBranchAddress("etagam1pi02", &etagam1pi02, &b_etagam1pi02);
   fChain->SetBranchAddress("phigam1pi02", &phigam1pi02, &b_phigam1pi02);
   fChain->SetBranchAddress("ptgam2pi02", &ptgam2pi02, &b_ptgam2pi02);
   fChain->SetBranchAddress("etagam2pi02", &etagam2pi02, &b_etagam2pi02);
   fChain->SetBranchAddress("phigam2pi02", &phigam2pi02, &b_phigam2pi02);
   fChain->SetBranchAddress("evtweight", &evtweight, &b_evtweight);
   Notify();
}

Bool_t cptauanalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void cptauanalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t cptauanalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef cptauanalysis_cxx
