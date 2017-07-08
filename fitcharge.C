#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"

void fitcharge(int iflag = 5) {


  gROOT->Reset(); 
  gROOT->ProcessLine(".L LGfitter.C");
  gStyle->SetOptFit(0111);
   TFile *MyInFile = new TFile("PurityAnalysis_data_prepostoutage_oldgain.root");
  if ( MyInFile->IsOpen() ) cout << "File opened successfully" << endl; 
  TString outfile = "fitresults.txt";
  ofstream fout;
  fout.open(outfile);
   

    TString stemp ="roodata_clone_h10__charge";
    //  TString stemp ="charge11";
    TH1 *bg; MyInFile->GetObject(stemp,bg);
    //bg->Rebin();
    // starting parameters for the fit - TUNE ME -
    Double_t bfp[4];
    Double_t bfpe[4];
    Double_t bcov[4];    
    Int_t ii = bg->GetEntries(); 
    Double_t mm = bg->GetBinCenter(bg->GetMaximumBin());
    Double_t amp = ii*mm/100.0;
    Double_t rms = bg->GetRMS();
    cout << " rms = " << rms << endl;
    Double_t mw = mm/100.0;
    // starting parameter values
    bfp[0]=mw*10.0; bfp[1]=mm; bfp[2]=amp; bfp[3]=mw*5.0;
    // step size for fitter
    bfpe[0]=0.1*mw; bfpe[1]=0.01*mm; bfpe[2]=0.05*amp; bfpe[3]=0.1*mw;    
     
    float lb,ub;
    int ibpeak = bg->GetMaximumBin();
    int iph = bg->GetBinContent(ibpeak);
    float low = 0.3*iph;
    float high = 0.5*iph;
    for (int ib=ibpeak;ib>0;ib--) {
      int itest = bg->GetBinContent(ib);
      if (itest<low) {lb=bg->GetBinCenter(ib); break;}
    }
    // cout << stemp << endl;
    cout << "lower bound " << lb << endl;
    for (int ib=ibpeak;ib<2*ibpeak;ib++) {
      int itest = bg->GetBinContent(ib);
      if (itest<high) {ub=bg->GetBinCenter(ib); break;}
    }
     cout << "upper bound " << ub << endl;
    
    // comment for now, narrow peak because no noise, no diff
    // if ((ub-lb)<200) {  bfp[3]=15.0; bfpe[3]=0.3; }    
    // else if ((ub-lb)<300) {  bfp[3]=20.0; bfpe[3]=0.5; }    
    // else if ((ub-lb)<400) {  bfp[3]=30.0; bfpe[3]=1.0; }    
    //    cout << "Histogram maximum " << mm << endl;
    //  Fit to Landau-Gaussian convolution
    TF1 *fitres = LGfitter(bg,lb,ub,bfp,bfpe,bcov,1);
    // Need to check for good fit results here    
    double totwid = bfp[0]+bfp[3];
    cout << "cov piece  " << bcov[3] << " " << bfpe[0] << " " << bfpe[3] << endl;
    
    double errtemp = bfpe[0]*bfpe[0]+bfpe[3]*bfpe[3]+2.0*bcov[3];
    double errtw;
    if (errtemp>0) {errtw=sqrt(errtemp);} 
    else {errtw=0.0;}
    cout << errtw << endl;
    // Write fit results to output file 
    TString myres=Form("%8.2f %6.2f %6.2f %8.2f %8.2f %6.2f %6.2f %8.2f %8.2f %8.2f %8.2f", 
			 1.0,bfp[1],bfpe[1],bfp[0],bfpe[0],
		       1e-4*bfp[2],1e-4*bfpe[2],bfp[3],bfpe[3],totwid,errtw);
    fout << myres << endl;
    
    // Draw it  
    //  bg->GetXaxis()->SetRange(0,75);
    TCanvas *c1 = new TCanvas("c1","bkgd",1600,1600);
    bg->Draw();
    c1->Print(stemp+".eps");  
  
  
  //  MyInFile->Close();  
    fout.close();

}
