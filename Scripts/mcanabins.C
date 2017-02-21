#define mcanabins_cxx
#include "mcanabins.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void mcanabins::Loop()
{
//   In a ROOT session, you can do:
//      root> .L mcanabins.C
//      root> mcanabins t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   std::map<int,std::vector<double> > puritymap;
   std::map<int,std::vector<double> > efficiencymap;
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (purity>-0.5) puritymap[bin].push_back(purity);
      if (efficiency > -0.5 && efficiency < 1.5) efficiencymap[bin].push_back(efficiency);
   }

   TCanvas * canv = new TCanvas("c","c",2000,1600);
   canv->Divide(2);
   TGraphErrors * pur = new TGraphErrors();
   TGraphErrors * eff = new TGraphErrors();
   for (int i = 0; i < 22; i++)
     {
       pur->SetPoint(i,i,TMath::Mean(puritymap[i].size(),puritymap[i].data()));
       pur->SetPointError(i,0,TMath::RMS(puritymap[i].size(),puritymap[i].data())/sqrt(puritymap[i].size()));
       eff->SetPoint(i,i,TMath::Mean(efficiencymap[i].size(),efficiencymap[i].data()));
       eff->SetPointError(i,0,TMath::RMS(efficiencymap[i].size(),efficiencymap[i].data())/sqrt(efficiencymap[i].size()));
     }
   canv->cd(1);
   pur->SetMarkerStyle(20);
   pur->Draw("alpe");
   canv->cd(2);
   eff->SetMarkerStyle(20);
   eff->Draw("alpe");
   canv->WaitPrimitive();
}
