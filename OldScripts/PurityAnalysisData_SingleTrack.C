#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TImage.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooGlobalFunc.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooArgSet.h"

#include <stdio.h>
#include <string>
#include <algorithm>
#include <map>

#ifdef __MAKECINT__
#pragma link C++ class vector<double>+;
#pragma link C++ class vector<int>+;
#endif

struct ChargeHistogram
{
  TH1D * hist;
  Double_t binlow;
  Double_t binhigh;
};

struct FitResult
{
  Double_t meanlandau;
  Double_t widthlandau;
  Double_t meangauss;
  Double_t widthgauss;
  Int_t status;
};

double hitGeomDist(TVector3 hitloc, TVector3 trigloc1, TVector3 trigloc2)
{
  return (((hitloc-trigloc1).Cross(hitloc-trigloc2)).Mag()/(trigloc2-trigloc1).Mag());
}

int CombineHits(std::vector<std::pair<double,double> > chanhits, std::vector<std::pair<double,double> > & consolidated, double threshold)
{
  Int_t numconsolidated = 0;
  size_t nhits = chanhits.size();
  for (size_t i = 0; i < nhits; )
    {
      double hittime = chanhits[i].first;
      double hitsize = chanhits[i].second;
      int n = 1;
      for (size_t j = i+1; fabs(chanhits[j-1].first-chanhits[j].first)<threshold && j < nhits; j++)
	{
	  n++;
	  hittime += chanhits[j].first;
	  hitsize += chanhits[j].second;
	}
      numconsolidated += n-1;
      hittime /= n;
      consolidated.push_back(std::make_pair(hittime,hitsize));
      i += n;
    }
  return numconsolidated;
}

void PurityAnalysisData_SingleTrack()
{
  TFile * file = TFile::Open("15xxx_reco_histfix.root","READ");
  TTreeReader reader("counterhits/CounterHitSelection",file);

  TTreeReaderValue<Int_t> nhits(reader,"nhits");
  TTreeReaderArray<Int_t> tpc(reader,"tpc");
  TTreeReaderArray<Double_t> hitx(reader,"hitx");
  TTreeReaderArray<Double_t> hity(reader,"hity");
  TTreeReaderArray<Double_t> hitz(reader,"hitz");
  TTreeReaderValue<Int_t> c1(reader,"c1");
  TTreeReaderValue<Int_t> c2(reader,"c2");
  TTreeReaderValue<Int_t> trignum(reader,"trignum");
  TTreeReaderValue<Double_t> c1x(reader,"c1x");
  TTreeReaderValue<Double_t> c1y(reader,"c1y");
  TTreeReaderValue<Double_t> c1z(reader,"c1z");
  TTreeReaderValue<Double_t> c2x(reader,"c2x");
  TTreeReaderValue<Double_t> c2y(reader,"c2y");
  TTreeReaderValue<Double_t> c2z(reader,"c2z");
  TTreeReaderValue<Int_t> run(reader,"run");
  TTreeReaderValue<Int_t> event(reader,"event");
  TTreeReaderArray<Double_t> perpdist(reader,"perpdist");
  TTreeReaderArray<Int_t> countercut(reader,"countercut");
  TTreeReaderArray<Double_t> peaktick(reader,"peaktick");
  TTreeReaderArray<Double_t> peaktime(reader,"peaktime");
  TTreeReaderArray<Double_t> hitt(reader,"hitt");
  TTreeReaderArray<Int_t> hitchan(reader,"hitchan");
  TTreeReaderArray<Int_t> hitsonchan(reader,"hitsonchan");
  TTreeReaderArray<Double_t> hitintegral(reader,"hitintegral");
  TTreeReaderArray<Double_t> hitsumadc(reader,"hitsumadc");
  TTreeReaderArray<Double_t> hitamplitude(reader,"hitamplitude");
  TTreeReaderArray<Double_t> driftdist(reader,"driftdist");
  TTreeReaderArray<Double_t> hitsigmaintegral(reader,"hitsigmaintegral");
  TTreeReaderArray<Double_t> hitsigmaamplitude(reader,"hitsigmaamplitude");
  TTreeReaderValue<Double_t> ransac_constant(reader,"ransac_constant");
  TTreeReaderValue<Double_t> ransac_constanterr(reader,"ransac_constanterr");
  TTreeReaderValue<Double_t> ransac_linear(reader,"ransac_linear");
  TTreeReaderValue<Double_t> ransac_linearerr(reader,"ransac_linearerr");
  TTreeReaderValue<Double_t> ransac_quadratic(reader,"ransac_quadratic");
  TTreeReaderValue<Double_t> ransac_quadraticerr(reader,"ransac_quadraticerr");
  TTreeReaderValue<Double_t> ransac_chi2(reader,"ransac_chi2");
  TTreeReaderValue<Double_t> ransac_ndf(reader,"ransac_ndf");
  TTreeReaderValue<Int_t> ransac_success(reader,"ransac_success");
  TTreeReaderValue<Double_t> ransac_sumsqrresidual(reader,"ransac_sumsqrresidual");
  TTreeReaderArray<Int_t> ransac_realhit(reader,"ransac_realhit");

  Int_t nentries = reader.GetEntries(true);

  UInt_t Nbins = 22;
  Double_t ADCcutoff = 1000;
  Double_t drifttimemax = 2011;//ticks

  // hitmap[run][event][chan](vector<hittime,hitsize>)
  std::map<Int_t, std::map<Int_t, std::map<Int_t, std::vector<std::pair<double,double> > > > > hitmap, consolidatedhitmap;

  Double_t pitch = 0.449;
  while (reader.Next())
    {
      Double_t sinthetacounter=1;
      if (*trignum == 111) sinthetacounter = TMath::Sin(TMath::ATan2(*c1z-*c2z,*c1x-*c2x));
      else if (*trignum == 112 || *trignum == 113) sinthetacounter = TMath::Cos(TMath::ATan2(*c1z-*c2z,*c1x-*c2x));
      for (Int_t i_hit = 0; i_hit < *nhits; i_hit++)
	{
	  if (hitx[i_hit] < -400 || hitx[i_hit] > 400) continue;
	  if (hitsigmaintegral[i_hit] > 10 || hitsumadc[i_hit] > 3000) continue;
	  if (!countercut[i_hit]) continue;
	  if (*ransac_sumsqrresidual/(*ransac_ndf) > 0.25) continue;
	  if (*ransac_chi2/(*ransac_ndf) > 0.5) continue;
	  if (tpc[i_hit] % 2 == 0) continue;
	  if (!ransac_realhit[i_hit]) continue;
	  if (*trignum == 111) continue;
	  if (fabs(pitch/sinthetacounter) < 1e-4) continue;
	  hitmap[*run][*event][hitchan[i_hit]].push_back(std::make_pair(hitt[i_hit],hitintegral[i_hit]/(pitch/sinthetacounter)));
	  //hitmap[*run][*event][hitchan[i_hit]].push_back(std::make_pair(hitt[i_hit],hitintegral[i_hit]));
	}
    }

  for (auto & i_run : hitmap)
    {
      for (auto & i_evt : i_run.second)
	{
	  for (auto & i_chan : i_evt.second)
	    {
	      std::sort(i_chan.second.begin(),i_chan.second.end());
	      CombineHits(i_chan.second,consolidatedhitmap[i_run.first][i_evt.first][i_chan.first],10);
	    }
	}
    }

  Double_t avgslope = 0;
  Int_t i_slope = 0;
  for (auto & i_run : consolidatedhitmap)
    {
      for (auto & i_evt : i_run.second)
	{
	  Int_t i_pt = 0;
	  TGraphErrors * trackgraph = new TGraphErrors();
	  for (auto & i_chan : i_evt.second)
	    {
	      for (auto & i_hit : i_chan.second)
		{
		  trackgraph->SetPoint(i_pt,i_hit.first,i_hit.second);
		  trackgraph->SetPointError(i_pt,10,sqrt(i_hit.second));
		  i_pt++;
		}
	    }
	  if (i_pt == 0) continue;
	  //TCanvas * canv = new TCanvas("c","canv",2400,1600);
	  //trackgraph->SetMarkerStyle(20);
	  //trackgraph->SetMarkerSize(2);
	  //trackgraph->Draw("AP");
	  trackgraph->Fit("expo","Q");
	  TF1 * ex = (TF1*)trackgraph->GetFunction("expo");
	  if (ex->GetParameter(1) > 0) continue;
	  //if (ex->GetChisquare()/ex->GetNDF() > 10) continue;
	  avgslope+=ex->GetParameter(1);
	  i_slope++;
	  //canv->Update();
	  //canv->WaitPrimitive();
	  //delete trackgraph;
	  //delete ex;
	  //delete canv;
	}
    }

  
  if (i_slope != 0) std::cout << "avg slope = " << avgslope/i_slope << "  i_slope=" << i_slope << std::endl;
  else std::cout << "no data" << std::endl;
}
