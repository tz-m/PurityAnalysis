//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec  8 13:07:00 2016 by ROOT version 6.08/00
// from TTree mcanabins/mcanabins
// found on file: reco_hist.root
//////////////////////////////////////////////////////////

#ifndef mcanabins_h
#define mcanabins_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class mcanabins {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           bin;
   Double_t        bincenter;
   Double_t        purity;
   Double_t        efficiency;
   Int_t           tp;
   Int_t           fp;
   Int_t           fn;

   // List of branches
   TBranch        *b_bin;   //!
   TBranch        *b_bincenter;   //!
   TBranch        *b_purity;   //!
   TBranch        *b_efficiency;   //!
   TBranch        *b_tp;   //!
   TBranch        *b_fp;   //!
   TBranch        *b_fn;   //!

   mcanabins(TTree *tree=0);
   virtual ~mcanabins();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mcanabins_cxx
mcanabins::mcanabins(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("robustmcana_3ms_2.0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("robustmcana_3ms_2.0.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("robustmcana_3ms_2.0.root:/robustmcana");
      dir->GetObject("mcanabins",tree);

   }
   Init(tree);
}

mcanabins::~mcanabins()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mcanabins::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mcanabins::LoadTree(Long64_t entry)
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

void mcanabins::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("bin", &bin, &b_bin);
   fChain->SetBranchAddress("bincenter", &bincenter, &b_bincenter);
   fChain->SetBranchAddress("purity", &purity, &b_purity);
   fChain->SetBranchAddress("efficiency", &efficiency, &b_efficiency);
   fChain->SetBranchAddress("tp", &tp, &b_tp);
   fChain->SetBranchAddress("fp", &fp, &b_fp);
   fChain->SetBranchAddress("fn", &fn, &b_fn);
   Notify();
}

Bool_t mcanabins::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mcanabins::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mcanabins::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mcanabins_cxx
