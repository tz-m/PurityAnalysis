#ifndef PURITYANALYSIS_H
#define PURITYANALYSIS_H

#include "Analysis.h"

#include "THStack.h"
#include "TFitResult.h"
#include "TSystem.h"
#include "TImage.h"
#include "TPaveLabel.h"
#include "TVirtualFitter.h"
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
#include "RooAddition.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooMinuit.h"

struct ChargeHistogram
{
  TH1F * hist;
  Float_t binlow;
  Float_t bincenter;
  Float_t binhigh;
  Float_t binwidth;
};

struct FitResult
{
  Float_t meanlandau;
  Float_t widthlandau;
  Float_t meangauss;
  Float_t widthgauss;
  Float_t meangauss2;
  Float_t widthgauss2;
  Int_t status;
  Float_t chi2ndf;
};

class PurityAnalysis : public Analysis {
public:
  void run();

  //void reconfigureCuts();

private:


};

void PurityAnalysis::run()
{
  //reconfigureCuts();

  UInt_t Nbins = 22;
  Float_t drifttimemax = 2012;
  Float_t ADCcutofflow = 400; //-700;
  Float_t ADCcutoffhigh = 5000; //23000;

  TH1F * found = new TH1F("found","Hits Found",200,ADCcutofflow,ADCcutoffhigh);
  TH1F * assumed = new TH1F("assumed","Assumed Hits",200,ADCcutofflow,ADCcutoffhigh);

  RooRealVar charge("charge","Hit Integral (ADC*tick)",ADCcutofflow,ADCcutoffhigh);
  RooRealVar tdrift("tdrift","tdrift",0,drifttimemax);
  RooDataSet * chgdata = new RooDataSet("chgdata","chgdata",RooArgSet(tdrift,charge));

  std::vector<ChargeHistogram> timeBinHist(Nbins);
  for (UInt_t i = 0; i < Nbins; i++)
    {
      timeBinHist[i].binlow = i*(drifttimemax/Nbins);
      timeBinHist[i].binhigh = (i+1)*(drifttimemax/Nbins);
      timeBinHist[i].bincenter = (timeBinHist[i].binlow + timeBinHist[i].binhigh)/2.0;
      timeBinHist[i].binwidth = timeBinHist[i].binhigh - timeBinHist[i].binlow;
      timeBinHist[i].hist = new TH1F(TString::Format("h%i",i),TString::Format("RobustHitFinder, Bin %i, %.2f <= t < %.2f (us)",i,timeBinHist[i].binlow,timeBinHist[i].binhigh),100,0,drifttimemax);
    }

  TH1F * hlong = new TH1F("hlong","Hit Integral (ADC)",1000,ADCcutofflow,ADCcutoffhigh);

  std::map<Int_t, std::vector<std::pair<Float_t,Float_t> > > hitmap;

  for (auto const & hititr : *(file->GetHitMap()))
    {
      const types::HitInfo * hit = &(hititr.second);
      if (hit->run >= 15483 && hit->run <= 15502) continue;
      if (hit->run == 15593) continue;
      if (hit->run >= 15634 && hit->run <= 15643) continue;
      if (hit->run >= 15664 && hit->run <= 15674) continue;
      if (hit->run == 15823 || hit->run == 15914 || hit->run == 15940 || hit->run == 15980 || hit->run == 16025 || hit->run == 16083 || hit->run == 16147 || hit->run == 16586) continue;
      if (hit->run >= 16593 && hit->run <= 16596) continue;
      if (hit->channel == 566 || hit->channel == 885 || hit->channel == 1547) continue;
      //if (hit->pedmean < 200 || hit->pedmean > 1500) continue;
      //if (hit->pedrms < 7 || hit->pedrms > 30) continue;
      //if (!(*countercut)) continue;
      //if ((*tpc) % 2 == 0) continue;
      //if (!(*fitrealhit)) continue;
      if (cuts.ChannelPass(hit) &&
          //cuts.CounterPass(hit) &&
          cuts.HitPass(hit))
        {
          if (!(hit->assumedhit)) found->Fill(hit->integral/hit->segmentlength);
          if (hit->assumedhit) assumed->Fill(hit->integral/hit->segmentlength);
          //std::cout << "making hit at tick " << hit->peaktick << " with dqdx " << hit->integral / hit->segmentlength << std::endl;
          if (!(hit->assumedhit)) hitmap[hit->channel].push_back(std::make_pair(hit->peaktime,(hit->integral)/(hit->segmentlength)));
        }
    }

  TCanvas * ctest = new TCanvas("ctest","ctest",2000,1600);
  ctest->cd();
  std::cout << "Found Number = " << found->GetEntries() << "  Assumed Number = " << assumed->GetEntries() << std::endl;
  THStack * hs = new THStack("hs","Found and Assumed Hits");
  found->SetFillColor(kRed);
  assumed->SetFillColor(kBlue);
  hs->Add(found);
  hs->Add(assumed);
  hs->Draw();
  //ctest->WaitPrimitive();

  for (auto & i_chan : hitmap)
    {
      for (auto & i_hit : i_chan.second)
        {
          if (i_hit.second < ADCcutoffhigh && i_hit.second > ADCcutofflow)
            {
              hlong->Fill(i_hit.second);
            }
          charge = (Float_t)(i_hit.second);
          tdrift = (Float_t)(i_hit.first);
          chgdata->add(RooArgSet(tdrift,charge));
          for (UInt_t i_hist = 0; i_hist < timeBinHist.size(); i_hist++)
            {
              if (i_hit.first >= timeBinHist[i_hist].binlow && i_hit.first < timeBinHist[i_hist].binhigh)
                {
                  timeBinHist[i_hist].hist->Fill(i_hit.second);
                  break;
                }
            }
        }
    }



  std::vector<FitResult> results(Nbins);

//RooWorkspace *workspace = new RooWorkspace("workspace","workspace");

  RooRealVar ml("LandMPV","mean landau",2800,1,10000);
  RooRealVar sl("LandWidth","sigma landau",1000,1,10000);
  RooLandau landau("lx","lx",charge,ml,sl);
  RooRealVar mg("mg","mean gauss",0);
  RooRealVar sg("GaussWidth","sigma gauss",50,1,1000);
  RooGaussian gauss("gauss","gauss",charge,mg,sg);
  RooFFTConvPdf lxg("lxg","landau (x) gauss",charge,landau,gauss);

/*
   workspace->import(lxg);
   workspace->import(landau);
   workspace->import(gauss);
   workspace->import(*chgdata);
 */
  TGraphAsymmErrors * hbin = new TGraphAsymmErrors(Nbins);
  hbin->SetTitle("Most Probable dQ/dx");
  //TH1F * hbin = new TH1F("hbin","Most Probable dQ/dx",Nbins,0,drifttimemax);
  Double_t minMPV = 99999;
  Double_t maxMPV = -99999;

  TGraphAsymmErrors * gw = new TGraphAsymmErrors(Nbins);
  gw->SetTitle("Width of Gaussian in Convolution");
  TGraphAsymmErrors * lw = new TGraphAsymmErrors(Nbins);
  lw->SetTitle("Width of Landau in Convolution");

  for (UInt_t i = 0; i < Nbins; i++)
    {
      if (timeBinHist[i].hist->GetEntries() > 0)
        {
          char shortname[100];
          sprintf(shortname,"dQdx_bin%d",i);

          Float_t maxbin = timeBinHist[i].hist->GetBinCenter(timeBinHist[i].hist->GetMaximumBin());
          //Float_t maxbin = 2800;

          // Setup component pdfs
          // --------------------

          //setup observable
          charge.setRange(shortname,maxbin-1000,maxbin+1500);
          //charge.setRange(shortname,ADCcutofflow,ADCcutoffhigh);

          //setup landau(t,ml,sl)
          ml.setVal(maxbin);
          ml.setMin(maxbin-1000); //0);
          ml.setMax(maxbin+1000); //0);
          sl.setVal(70); //850);
          sl.setMin(50); //60);
          sl.setMax(200); //2000);

          //setup gauss(t,mg,sg)
          //mg.setVal(1);
          //mg.setMin(-100);
          //mg.setMax(100);
          sg.setVal(300); //550);
          sg.setMin(100); //10);
          sg.setMax(500); //2000);

          // construct convolution pdf
          // -------------------------

          // set num bins to be used for FFT sampling
          charge.setBins(10000,"fft");

          // fit convoluted pdf to binHist
          // -----------------------------

          char cut[100];
          sprintf(cut,"tdrift > %f && tdrift < %f",timeBinHist[i].binlow,timeBinHist[i].binhigh);
          RooDataSet* roodata = (RooDataSet*) chgdata->reduce(RooFit::SelectVars(charge),RooFit::Cut(cut),RooFit::Name(shortname),RooFit::Title(shortname));
          //RooDataHist* roodata = new RooDataHist("roodata","roodata histogram",RooArgSet(charge),*roodataset);


          //RooAbsReal* nll = lxg.createNLL(*roodata,RooFit::NumCPU(4));
          //RooMinuit m(*nll);
          //m.migrad();
          //m.hesse();
          //RooFitResult * rooresult = m.save();
          //RooFitResult * rooresult =
          lxg.fitTo(*roodata,RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range(shortname));
          RooPlot * frame = charge.frame(RooFit::Title(timeBinHist[i].hist->GetTitle()));
          roodata->plotOn(frame,RooFit::Binning(200));
          lxg.plotOn(frame,RooFit::LineColor(kRed),RooFit::ShiftToZero());

          TPaveLabel * tl = new TPaveLabel(0.6,0.82,0.99,0.9,Form(" #chi^{2}/NDF = %f",frame->chiSquare()),"NDC");
          tl->SetFillColor(0);
          tl->SetBorderSize(1);
          tl->SetTextAlign(12);
          frame->addObject(tl);
          frame->getAttText()->SetTextSize(0.32);
          roodata->statOn(frame,RooFit::Layout(0.6,0.99,0.8));
          frame->getAttText()->SetTextSize(0.025);
          lxg.paramOn(frame,RooFit::Layout(0.6,0.99,0.6));
          frame->getAttText()->SetTextSize(0.025);

          results[i].meanlandau = ml.getVal();
          results[i].widthlandau = sl.getVal();
          results[i].meangauss = mg.getVal();
          results[i].widthgauss = sg.getVal();
          results[i].chi2ndf = frame->chiSquare();

          if (minMPV > results[i].meanlandau) minMPV = results[i].meanlandau;
          if (maxMPV < results[i].meanlandau) maxMPV = results[i].meanlandau;

          char buff[100];
          sprintf(buff,"bin%d",i);
          TCanvas* c = new TCanvas(buff,buff,1600,1600);
          gPad->SetLeftMargin(0.15);
          frame->GetYaxis()->SetTitleOffset(1.4);
          frame->Draw();
          gSystem->ProcessEvents();
          TImage *img = TImage::Create();
          img->FromPad(c);
          char newbuff[100];
          sprintf(newbuff,"%s.png",buff);
          img->WriteImage(newbuff);

          char buffnew[100];
          sprintf(buffnew,"bin%d_resid",i);
          TCanvas* cresid = new TCanvas(buffnew,buffnew,1600,1600);
          gPad->SetLeftMargin(0.15);
          RooPlot* frame2 = charge.frame(RooFit::Title("residual"));
          RooHist* hresid = frame->residHist();
          frame2->addPlotable(hresid,"P");
          frame2->GetYaxis()->SetTitleOffset(1.4);
          frame2->Draw();
          gSystem->ProcessEvents();
          TImage *img2 = TImage::Create();
          img2->FromPad(cresid);
          char buffnew2[100];
          sprintf(buffnew2,"%s.png",buffnew);
          img2->WriteImage(buffnew2);

          char buffpull[100];
          sprintf(buffpull,"bin%d_pull",i);
          TCanvas* cpull = new TCanvas(buffpull,buffpull,1600,1600);
          gPad->SetLeftMargin(0.15);
          RooPlot* frame3 = charge.frame(RooFit::Title("pull"));
          RooHist* hpull2 = frame->pullHist();
          TH1F* hpull = new TH1F("hpull","pull",100,0,0);
          for (Int_t ibin=0; ibin<hpull2->GetN(); ibin++)
            {
              Double_t x,y;
              hpull2->GetPoint(ibin,x,y);
              hpull->Fill(y);
            }
          hpull->Draw();
          gSystem->ProcessEvents();
          TImage *img3 = TImage::Create();
          img3->FromPad(cpull);
          char buffpull2[100];
          sprintf(buffpull2,"%s.png",buffpull);
          img3->WriteImage(buffpull2);

          //c->WaitPrimitive();
          delete c;
          delete cresid;
          delete cpull;
          delete frame;
          delete frame2;
          delete frame3;
          delete img;
          delete img2;
          delete img3;

          hbin->SetPoint(i,timeBinHist[i].bincenter,ml.getVal());
          hbin->SetPointError(i,(timeBinHist[i].binwidth)/2.0, (timeBinHist[i].binwidth)/2.0,ml.getError(),ml.getError());
          gw->SetPoint(i,timeBinHist[i].bincenter,sg.getVal());
          gw->SetPointError(i,(timeBinHist[i].binwidth)/2.0,(timeBinHist[i].binwidth)/2.0,sg.getError(),sg.getError());
          lw->SetPoint(i,timeBinHist[i].bincenter,sl.getVal());
          lw->SetPointError(i,(timeBinHist[i].binwidth)/2.0,(timeBinHist[i].binwidth)/2.0,sl.getError(),sl.getError());
        }
    }

  for (UInt_t i = 0; i < Nbins; i++)
    {
      std::cout << "Bin " << i << ": ml=" << results[i].meanlandau << "  sl=" << results[i].widthlandau << "  mg=" << results[i].meangauss << "  sg=" << results[i].widthgauss << "  chi2ndf=" << results[i].chi2ndf << std::endl;
    }

  TF1 * expo = new TF1("expo","[0]*exp(-x/[1])",0,drifttimemax);
  expo->SetParNames("dQdx0","eLifetime");
  expo->SetParameters(3000,3000);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  TCanvas * canv2 = new TCanvas("canv2","canv2",1600,800);
  canv2->cd();
  hbin->Fit("expo","R");
  hbin->Draw("ap");
  hbin->SetMarkerStyle(kFullDotLarge);
  hbin->SetMarkerSize(2);

  TH1F * hint = new TH1F("hint","95#% confidence band",100,0,drifttimemax);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint,0.68);
  hint->SetStats(kFALSE);
  hint->SetFillColorAlpha(kRed,0.35);
  hint->Draw("e3 same");

  hbin->GetYaxis()->SetRangeUser(minMPV-100,maxMPV+100);
  hbin->GetXaxis()->SetTitle("Drift Time (#mu s)");
  hbin->GetYaxis()->SetTitle("dQ/dx (ADC/cm)");

  TCanvas * canv3 = new TCanvas("canv3","canv3",1600,800);
  canv3->cd();
  gw->Draw("ap");
  gw->SetMarkerStyle(kFullDotLarge);
  gw->SetMarkerSize(2);
  gw->GetXaxis()->SetTitle("Drift Time (#mu s)");

  TCanvas * canv4 = new TCanvas("canv4","canv4",1600,800);
  canv4->cd();
  lw->Draw("ap");
  lw->SetMarkerStyle(kFullDotLarge);
  lw->SetMarkerSize(2);
  lw->GetXaxis()->SetTitle("Drift Time (#mu s)");

/*
   TCanvas * canv3 = new TCanvas("canv3","canv3",1600,800);
   canv3->cd();
   hlong->Draw();
 */
}

#endif
