#ifndef PURITYANALYSIS_H
#define PURITYANALYSIS_H

#include "Analysis.h"

#include "TFitResult.h"
#include "TSystem.h"
#include "TImage.h"
#include "TPaveLabel.h"
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
  Float_t binhigh;
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
  Float_t g2amp;
  Float_t g3amp;
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
  Float_t ADCcutoff = 10000;

  RooRealVar charge("charge","Hit Integral (ADC*tick)",0,ADCcutoff);
  RooRealVar tdrift("tdrift","tdrift",0,drifttimemax);
  RooDataSet * chgdata = new RooDataSet("chgdata","chgdata",RooArgSet(tdrift,charge));

  std::vector<ChargeHistogram> timeBinHist(Nbins);
  for (UInt_t i = 0; i < Nbins; i++)
    {
      timeBinHist[i].binlow = i*(drifttimemax/Nbins);
      timeBinHist[i].binhigh = (i+1)*(drifttimemax/Nbins);
      timeBinHist[i].hist = new TH1F(TString::Format("h%i",i),TString::Format("RobustHitFinder, Bin %i, %.2f <= t < %.2f (us)",i,timeBinHist[i].binlow,timeBinHist[i].binhigh),100,0,ADCcutoff);
    }

  TH1F * hlong = new TH1F("hlong","Hit Integral (ADC)",1000,0,0);

  std::map<Int_t, std::vector<std::pair<Float_t,Float_t> > > hitmap;

  for (auto const & hititr : *(file.GetHitMap()))
    {
      const types::HitInfo * hit = &(hititr.second);
      //if (!(*countercut)) continue;
      //if ((*tpc) % 2 == 0) continue;
      //if (!(*fitrealhit)) continue;
      if (cuts.ChannelPass(hit) &&
          cuts.CounterPass(hit) &&
          cuts.HitPass(hit))
        {
          hitmap[hit->channel].push_back(std::make_pair(hit->hitt,(hit->integral)/(hit->segmentlength)));
        }
    }

  for (auto & i_chan : hitmap)
    {
      for (auto & i_hit : i_chan.second)
        {
          if (i_hit.second < ADCcutoff)
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
  TH1F * hbin = new TH1F("hbin","LandMPV",Nbins,0,2011);

  for (UInt_t i = 0; i < Nbins; i++)
    {
      if (timeBinHist[i].hist->GetEntries() > 0)
        {
          char shortname[100];
          sprintf(shortname,"dQdx_bin%d",i);

          //Float_t maxbin = timeBinHist[i].hist->GetBinCenter(timeBinHist[i].hist->GetMaximumBin());
          Float_t maxbin = 2800;

          // Setup component pdfs
          // --------------------

          //setup observable
          //charge.setRange(shortname,0.1*maxbin,3*maxbin);
          charge.setRange(shortname,0,ADCcutoff);

          //setup landau(t,ml,sl)
          ml.setVal(maxbin);
          ml.setMin(maxbin-1000);
          ml.setMax(maxbin+1000);
          sl.setVal(850);
          sl.setMin(60);
          sl.setMax(1000);

          //setup gauss(t,mg,sg)
          //mg.setVal(1);
          //mg.setMin(-100);
          //mg.setMax(100);
          sg.setVal(350);
          sg.setMin(100);
          sg.setMax(800);

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
          lxg.fitTo(*roodata,RooFit::Extended(true),RooFit::Minimizer("Minuit2","migrad"),RooFit::Save(),RooFit::PrintLevel(-1),RooFit::Range(shortname));
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

          hbin->SetBinContent(i,ml.getVal());
          hbin->SetBinError(i,ml.getError());

        }
    }

  for (UInt_t i = 0; i < Nbins; i++)
    {
      std::cout << "Bin " << i << ": ml=" << results[i].meanlandau << "  sl=" << results[i].widthlandau << "  mg=" << results[i].meangauss << "  sg=" << results[i].widthgauss << "  chi2ndf=" << results[i].chi2ndf << std::endl;
    }

  TCanvas * canv2 = new TCanvas("canv2","canv2",1600,800);
  canv2->cd();
  hbin->Draw();

  TCanvas * canv3 = new TCanvas("canv3","canv3",1600,800);
  canv3->cd();
  hlong->Draw();

}

#endif
