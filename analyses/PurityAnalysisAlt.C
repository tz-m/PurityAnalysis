#include "THStack.h"
#include "TFitResult.h"
#include "TSystem.h"
#include "TImage.h"
#include "TPaveLabel.h"
#include "TVirtualFitter.h"
#include "TColor.h"
#include "TFrame.h"
#include "TProfile.h"
#include "TUnfoldDensity.h"
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
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooMinuit.h"
#include "RooAbsPdf.h"
#include "RooPolynomial.h"
#include "RooProduct.h"
#include "Math/VavilovAccurate.h"
#include "RooFunctorBinding.h"
#include "Math/Functor.h"
#include "TLatex.h"
#include "RooFormulaVar.h"

#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"

#include <string>
#include <sstream>
#include <map>
#include <algorithm>

#ifdef __MAKECINT__
#pragma link C++ class vector<double>+;
#pragma link C++ class vector<int>+;
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<bool>+;
#endif

struct ChargeHistogram
{
  TH1F * hist;
  TH1F * histMCT;
  Float_t binlow;
  Float_t bincenter;
  Float_t binhigh;
  Float_t binwidth;
};

struct FitResult
{
  Int_t binnum;
  Float_t landmpv;
  Float_t landwidth;
  Float_t gaussmean;
  Float_t gausswidth;
  Float_t gaussmean_l;
  Float_t gausswidth_l;
  Float_t landmpverr;
  Float_t landwidtherr;
  Float_t gaussmeanerr;
  Float_t gausswidtherr;
  Float_t gaussmean_lerr;
  Float_t gausswidth_lerr;
  Float_t bincenter;
  Float_t binwidtherr;
  Float_t lxg_chi2ndf;
  Float_t pg_chi2ndf;
  Float_t lxglow;
  Float_t lxghigh;
  Bool_t fitsuccess;
  Bool_t MCT;
};


struct HitSave {
  Float_t hitt;
  Float_t dqdx;
  Bool_t asshit;
};

struct HitInfo {
  Int_t run;
  Int_t subrun;
  Int_t event;
  Double_t t0;
  UInt_t c1;
  UInt_t c2;
  UInt_t trignum;
  Float_t c1x;
  Float_t c1y;
  Float_t c1z;
  Float_t c2x;
  Float_t c2y;
  Float_t c2z;
  Float_t distancecut;
  std::vector<Int_t> channel;
  std::vector<Int_t> wire;
  std::vector<Int_t> tpc;
  std::vector<Int_t> signalsize;
  //std::vector<Float_t> signal;
  //std::vector<Float_t> signalFilter;
  std::vector<Float_t> baseline;
  std::vector<Float_t> rms;
  std::vector<Float_t> baselineFilter;
  std::vector<Float_t> rmsFilter;
  std::vector<Float_t> pedmean;
  std::vector<Float_t> pedrms;
  std::vector<Float_t> integral;
  std::vector<Float_t> integralFilter;
  std::vector<Float_t> sumADC;
  std::vector<Float_t> sigmaintegral;
  std::vector<Float_t> sigmaintegralFilter;
  std::vector<Float_t> amplitude;
  std::vector<Float_t> amplitudeFilter;
  std::vector<Float_t> peaktick;
  std::vector<Float_t> peaktickFilter;
  std::vector<Float_t> peaktime;
  std::vector<Float_t> peaktimeFilter;
  std::vector<Int_t> begintick;
  std::vector<Int_t> endtick;
  std::vector<Int_t> width;
  std::vector<Float_t> hitx;
  std::vector<Float_t> hity;
  std::vector<Float_t> hitz;
  std::vector<Float_t> hiterrxlo;
  std::vector<Float_t> hiterrxhi;
  std::vector<Float_t> hiterrylo;
  std::vector<Float_t> hiterryhi;
  std::vector<Float_t> hiterrzlo;
  std::vector<Float_t> hiterrzhi;
  std::vector<Float_t> perpdist;
  std::vector<Float_t> hitt;
  std::vector<Float_t> driftdist;
  std::vector<Bool_t> countercut;
  std::vector<Bool_t> ismctruth;
  Float_t fitconstant;
  Float_t fitconstanterr;
  Float_t fitlinear;
  Float_t fitlinearerr;
  Float_t fitquadratic;
  Float_t fitquadraticerr;
  Float_t fitchi2;
  Float_t fitsumsqrresidual;
  Float_t fitndf;
  Float_t fitmle;
  Bool_t fitsuccess;
  std::vector<Bool_t> fitrealhit;
  std::vector<Float_t> segmentlength;
  std::vector<Bool_t> assumedhit;
  std::vector<Int_t> numGoodHitsChan;
  Int_t nwiresTPC0;
  Int_t nwiresTPC1;
  Int_t nwiresTPC2;
  Int_t nwiresTPC3;
  Int_t nwiresTPC4;
  Int_t nwiresTPC5;
  Int_t nwiresTPC6;
  Int_t nwiresTPC7;
};



class PurityAnalysisAlt {
public:
  PurityAnalysisAlt();
  PurityAnalysisAlt(Float_t mcs, Int_t life, Bool_t doMCT);
  void initialize(Bool_t data, Float_t mcs, Int_t life, Bool_t doMCT);
  void run();

private:

  std::vector<HitInfo> hitvec;

  int getBinNumber( Float_t t, const std::vector<ChargeHistogram> & tbh );

  bool analyse_MCTruth;
  TString infname;
  TString outfname;

  Float_t mcscale;
  Int_t lifetimesim;

};

PurityAnalysisAlt::PurityAnalysisAlt()
{
  // Default constructor -- DATA
  // bogus parameters for mcs life and doMCT
  this->initialize(true,1.0,4000,false);
  this->run();
}

PurityAnalysisAlt::PurityAnalysisAlt(Float_t mcs, Int_t life, Bool_t doMCT)
{
  // Alt constructor -- SIMULATION
  // By default, do MC truth fits as well
  this->initialize(false,mcs,life,doMCT);
  this->run();
}

void PurityAnalysisAlt::initialize(Bool_t data, Float_t mcs, Int_t life, Bool_t doMCT)
{
  analyse_MCTruth = doMCT;
  mcscale = mcs;
  lifetimesim = life;
  if (data)
    {
      infname = "/media/mthiesse/Dell Portable Hard Drive/PurityData/robustreco_data_hist.root";
      outfname = "PurityAnalysis_data.root";
    }
  else
    {
      infname = TString::Format("/media/mthiesse/Dell Portable Hard Drive/PurityData/robust_oldgain_%ius_mcscale%.1f_hist2.root",life,mcscale);
      outfname = TString::Format("/home/mthiesse/PurityAnalysis/PurityAnalysis_%ius_mcscale%.1f.root",life,mcscale);
    }
  std::cout << "infname = " << infname << std::endl;
  std::cout << "outfname = " << outfname << std::endl;

  TFile * file = TFile::Open(infname,"READ");
  if (!file || file->IsZombie())
    {
      std::stringstream ss;
      ss << "ReadHistFile::ReadFile() -- Input file is not read";
      throw std::runtime_error(ss.str());
    }
  TTreeReader reader("robusthit/RobustHitFinder",file);
  TTreeReaderValue<Int_t> run(reader,"run");
  TTreeReaderValue<Int_t> subrun(reader,"subrun");
  TTreeReaderValue<Int_t> event(reader,"event");
  TTreeReaderValue<Double_t> t0(reader,"t0");
  TTreeReaderValue<UInt_t> c1(reader,"c1");
  TTreeReaderValue<UInt_t> c2(reader,"c2");
  TTreeReaderValue<UInt_t> trignum(reader,"trignum");
  TTreeReaderValue<Float_t> c1x(reader,"c1x");
  TTreeReaderValue<Float_t> c1z(reader,"c1z");
  TTreeReaderValue<Float_t> c2x(reader,"c2x");
  TTreeReaderValue<Float_t> c2z(reader,"c2z");
  TTreeReaderValue<Float_t> distancecut(reader,"distancecut");
  TTreeReaderArray<Int_t> channel(reader,"channel");
  TTreeReaderArray<Int_t> wire(reader,"wire");
  TTreeReaderArray<Int_t> tpc(reader,"tpc");
  TTreeReaderArray<Float_t> rms(reader,"rms");
  TTreeReaderArray<Float_t> baseline(reader,"baseline");
  TTreeReaderArray<Float_t> pedmean(reader,"pedmean");
  TTreeReaderArray<Float_t> pedrms(reader,"pedrms");
  TTreeReaderArray<Float_t> integral(reader,"integral");
  TTreeReaderArray<Float_t> sumADC(reader,"sumADC");
  TTreeReaderArray<Float_t> sigmaintegral(reader,"sigmaintegral");
  TTreeReaderArray<Float_t> amplitude(reader,"amplitude");
  TTreeReaderArray<Float_t> peaktick(reader,"peaktick");
  TTreeReaderArray<Float_t> peaktime(reader,"peaktime");
  TTreeReaderArray<Int_t> begintick(reader,"begintick");
  TTreeReaderArray<Int_t> endtick(reader,"endtick");
  TTreeReaderArray<Int_t> width(reader,"width");
  TTreeReaderArray<Float_t> hitx(reader,"hitx");
  TTreeReaderArray<Float_t> hitz(reader,"hitz");
  TTreeReaderArray<Float_t> hiterrxlo(reader,"hiterrxlo");
  TTreeReaderArray<Float_t> hiterrxhi(reader,"hiterrxhi");
  TTreeReaderArray<Float_t> hiterrzlo(reader,"hiterrzlo");
  TTreeReaderArray<Float_t> hiterrzhi(reader,"hiterrzhi");
  TTreeReaderArray<Float_t> perpdist(reader,"perpdist");
  TTreeReaderArray<Float_t> hitt(reader,"hitt");
  TTreeReaderArray<Float_t> driftdist(reader,"driftdist");
  TTreeReaderValue<std::vector<Bool_t> > countercut(reader,"countercut");
  TTreeReaderValue<Float_t> fitconstant(reader,"fitconstant");
  TTreeReaderValue<Float_t> fitconstanterr(reader,"fitconstanterr");
  TTreeReaderValue<Float_t> fitlinear(reader,"fitlinear");
  TTreeReaderValue<Float_t> fitlinearerr(reader,"fitlinearerr");
  TTreeReaderValue<Float_t> fitquadratic(reader,"fitquadratic");
  TTreeReaderValue<Float_t> fitquadraticerr(reader,"fitquadraticerr");
  TTreeReaderValue<Float_t> fitchi2(reader,"fitchi2");
  TTreeReaderValue<Float_t> fitsumsqrresidual(reader,"fitsumsqrresidual");
  TTreeReaderValue<Float_t> fitndf(reader,"fitndf");
  TTreeReaderValue<Float_t> fitmle(reader,"fitmle");
  TTreeReaderValue<Bool_t> fitsuccess(reader,"fitsuccess");
  TTreeReaderValue<std::vector<Bool_t> > fitrealhit(reader,"fitrealhit");
  TTreeReaderValue<std::vector<Bool_t> > assumedhit(reader,"assumedhit");
  TTreeReaderArray<Float_t> segmentlength(reader,"segmentlength");
  TTreeReaderArray<Int_t> numGoodHitsChan(reader,"numGoodHitsChan");
  TTreeReaderValue<std::vector<Bool_t> > ismctruth(reader,"ismctruth");

  while (reader.Next())
    {
      HitInfo hf;
      hf.run = *run;
      //hf.subrun = *subrun;
      hf.event = *event;
      hf.t0 = *t0;
      hf.c1 = *c1;
      hf.c2 = *c2;
      hf.trignum = *trignum;
      //hf.c1x = *c1x;
      //hf.c1z = *c1z;
      //hf.c2x = *c2x;
      //hf.c2z = *c2z;
      //hf.distancecut = *distancecut;

      int nchan = 0;
      for (auto j : channel) ++nchan;
      for (int j = 0; j < nchan; ++j) hf.channel.push_back(channel.At(j));

      //int nwire = 0;
      //for (auto j : wire) ++nwire;
      //for (int j = 0; j < nwire; ++j) hf.wire.push_back(wire.At(j));

      int ntpc = 0;
      for (auto j : tpc) ++ntpc;
      for (int j = 0; j < ntpc; ++j) hf.tpc.push_back(tpc.At(j));

      int nrms = 0;
      for (auto j : rms) ++nrms;
      for (int j = 0; j < nrms; ++j) hf.rms.push_back(rms.At(j));

      int nbase = 0;
      for (auto j : baseline) ++nbase;
      for (int j = 0; j < nbase; ++j) hf.baseline.push_back(baseline.At(j));

      //int npm = 0;
      //for (auto j : pedmean) ++npm;
      //for (int j = 0; j < npm; ++j) hf.pedmean.push_back(pedmean.At(j));

      //int npr = 0;
      //for (auto j : pedrms) ++npr;
      //for (int j = 0; j < npr; ++j) hf.pedrms.push_back(pedrms.At(j));

      int nint = 0;
      for (auto j : integral) ++nint;
      for (int j = 0; j < nint; ++j) hf.integral.push_back(integral.At(j));

      //for (auto j : sumADC) hf.sumADC.push_back(*j);
      //for (auto j : sigmaintegral) hf.sigmaintegral.push_back(*j);
      //for (auto j : amplitude) hf.amplitude.push_back(*j);
      //for (auto j : peaktick) hf.peaktick.push_back(*j);
      //for (auto j : peaktime) hf.peaktime.push_back(*j);
      //for (auto j : begintick) hf.begintick.push_back(*j);
      //for (auto j : endtick) hf.endtick.push_back(*j);

      int nwidth = 0;
      for (auto j : width) ++nwidth;
      for (int j = 0; j < nwidth; ++j) hf.width.push_back(width.At(j));

      //for (auto j : hitx) hf.hitx.push_back(*j);
      //for (auto j : hitz) hf.hitz.push_back(*j);
      //for (auto j : hiterrxlo) hf.hiterrxlo.push_back(*j);
      //for (auto j : hiterrxhi) hf.hiterrxhi.push_back(*j);
      //for (auto j : hiterrzlo) hf.hiterrzlo.push_back(*j);
      //for (auto j : hiterrzhi) hf.hiterrzhi.push_back(*j);
      //for (auto j : perpdist) hf.perpdist.push_back(*j);

      int nhitt = 0;
      for (auto j : hitt) ++nhitt;
      for (int j = 0; j < nhitt; ++j) hf.hitt.push_back(hitt.At(j));

      //int ndist = 0;
      //for (auto j : driftdist) ++ndist;
      //for (int j = 0; j < ndist; ++j) hf.driftdist.push_back(driftdist.At(j));

      //for (auto j : countercut) hf.countercut.push_back(*j);

      hf.fitconstant = *fitconstant;
      //hf.fitconstanterr = *fitconstanterr;
      //hf.fitlinear = *fitlinear;
      //hf.fitlinearerr = *fitlinearerr;
      hf.fitquadratic = *fitquadratic;
      //hf.fitquadraticerr = *fitquadraticerr;
      //hf.fitchi2 = *fitchi2;
      hf.fitsumsqrresidual = *fitsumsqrresidual;
      hf.fitndf = *fitndf;
      //hf.fitmle = *fitmle;
      hf.fitsuccess = *fitsuccess;

      int nreal = 0;
      for (auto j : *fitrealhit) ++nreal;
      for (int j = 0; j < nreal; ++j) hf.fitrealhit.push_back((*fitrealhit)[j]);

      int nass = 0;
      for (auto j : *assumedhit) ++nass;
      for (int j = 0; j < nass; ++j) hf.assumedhit.push_back((*assumedhit)[j]);

      int nseg = 0;
      for (auto j : segmentlength) ++nseg;
      for (int j = 0; j < nseg; ++j) hf.segmentlength.push_back(segmentlength.At(j));

      int nismc = 0;
      for (auto j : *ismctruth) ++nismc;
      for (int j = 0; j < nismc; ++j) hf.ismctruth.push_back((*ismctruth)[j]);

      //for (auto j : numGoodHitsChan) hf.numGoodHitsChan.push_back(*j);

      if (nchan != ntpc || ntpc != nrms || nrms != nbase || nbase != nint || nint != nwidth || nwidth != nhitt || nhitt != nreal || nreal != nass || nass != nseg)
        {
          std::cout << "Shit's fuk'd yo!" << std::endl;
          return;
        }

      hitvec.push_back(hf);
    }

  std::cout << hitvec.size() << " events read" << std::endl;
  delete file;
}

int PurityAnalysisAlt::getBinNumber( Float_t t, const std::vector<ChargeHistogram> & tbh )
{
  for ( UInt_t i_hist = 0; i_hist < tbh.size(); i_hist++ )
    {
      if ( t >= tbh[i_hist].binlow && t < tbh[i_hist].binhigh )
        {
          return ( int )i_hist;
        }
    }
  return -1;
}

void PurityAnalysisAlt::run()
{
  bool DoPull = false;
  bool DoResidual = false;
  bool DoLxG = true;
  bool DoPG = false;
  float Nbinsfloat = 22;

  UInt_t Nbins = ( UInt_t )Nbinsfloat;
  Float_t drifttimemax = 2012;
  Float_t ADCcutofflow = 0;
  Float_t ADCcutoffhigh = 6000;
  Float_t ADCcutofflowMCT = 800;
  Float_t ADCcutoffhighMCT = 3000;

  TFile * fileout = TFile::Open( outfname, "recreate" );
  Int_t runnum;
  Int_t event;
  Float_t hitt;
  Int_t binnum;
  Bool_t foundrealhit;
  Bool_t assumedhit;
  Float_t dqdx;
  Int_t channel;
  Int_t tpc;
  Bool_t ismctruth;
  TTree * hitstree = new TTree( "hits", "Hit Information" );
  hitstree->Branch( "runnum", &runnum, "runnum/I" );
  hitstree->Branch( "event", &event, "event/I" );
  hitstree->Branch( "hitt", &hitt, "hitt/F" );
  hitstree->Branch( "binnum", &binnum, "binnum/I" );
  hitstree->Branch( "foundrealhit", &foundrealhit, "foundrealhit/O" );
  hitstree->Branch( "assumedhit", &assumedhit, "assumedhit/O" );
  hitstree->Branch( "ismctruth", &ismctruth, "ismctruth/O");
  hitstree->Branch( "dqdx", &dqdx, "dqdx/F" );
  hitstree->Branch( "channel", &channel, "channel/I" );
  hitstree->Branch( "tpc", &tpc, "tpc/I" );

  TH1F * found = new TH1F( "found", "Hits Found", 200, ADCcutofflow, ADCcutoffhigh );
  TH1F * assumed = new TH1F( "assumed", "Assumed Hits", 200, ADCcutofflow, ADCcutoffhigh );

  RooRealVar charge( "charge", "dQ/dx [ADC/cm]", ADCcutofflow, ADCcutoffhigh );
  RooRealVar chargeMCT( "chargeMCT","dQ/dx [ADC/cm]", ADCcutofflowMCT,ADCcutoffhighMCT );
  RooRealVar tdrift( "tdrift", "tdrift", 0, drifttimemax );
  RooRealVar tdriftMCT( "tdriftMCT","tdriftMCT",0,drifttimemax);
  RooDataSet * chgdata = new RooDataSet( "chgdata", "chgdata", RooArgSet( tdrift, charge ), RooFit::StoreError( RooArgSet( charge ) ));
  chgdata->Print();
  RooDataSet * chgdataMCT = new RooDataSet( "chgdataMCT", "chgdataMCT", RooArgSet(tdriftMCT,chargeMCT),RooFit::StoreError(RooArgSet(chargeMCT)));
  chgdataMCT->Print();

  std::vector<ChargeHistogram> timeBinHist( Nbins );
  for ( UInt_t i = 0; i < Nbins; i++ )
    {
      timeBinHist[i].binlow = i*( drifttimemax/Nbins );
      timeBinHist[i].binhigh = ( i+1 )*( drifttimemax/Nbins );
      timeBinHist[i].bincenter = ( timeBinHist[i].binlow + timeBinHist[i].binhigh )/2.0;
      timeBinHist[i].binwidth = timeBinHist[i].binhigh - timeBinHist[i].binlow;
      timeBinHist[i].hist = new TH1F( TString::Format( "h%i", i ), TString::Format( "RobustHitFinder, Bin %i, %.2f <= t < %.2f ( us )", i, timeBinHist[i].binlow, timeBinHist[i].binhigh ), 200, ADCcutofflow, ADCcutoffhigh );
      timeBinHist[i].hist->Sumw2();
      timeBinHist[i].histMCT = new TH1F( TString::Format( "h%iMCT", i ), TString::Format( "RobustHitFinder, Bin %i, %.2f <= t < %.2f ( us ), MCT", i, timeBinHist[i].binlow, timeBinHist[i].binhigh ), 200, ADCcutofflowMCT, ADCcutoffhighMCT );
      timeBinHist[i].histMCT->Sumw2();
    }

  std::map<Int_t, std::vector<HitSave> > hitmap;
  std::map<Int_t, std::vector<HitSave> > hitmapMCT;

  std::map<UInt_t, std::pair<Int_t, Int_t> > fabinmap;
  for ( UInt_t i = 0; i < Nbins; i++ )
    {
      fabinmap[i] = std::make_pair( 0, 0 );
    }

  std::map<int,std::map<int,int> > runseventshits;
  std::set<int> runset;
  for ( auto hit : hitvec)
    {
      if (hit.run > 17512) continue;      // pump broke immediately following this run
      if (hit.run > 15000) continue;  // not sufficient data in later runs to justify using them
      if (hit.trignum != 111) continue;
      if ( (hit.c1 == 28 && (hit.c2 == 6 || hit.c2 == 7 || hit.c2 == 8 || hit.c2 == 9)) ||
           (hit.c1 == 29 && (hit.c2 == 6 || hit.c2 == 7 || hit.c2 == 8)) ||
           (hit.c1 == 30 && (hit.c2 == 6 || hit.c2 == 7)) ||
           (hit.c1 == 31 && (hit.c2 == 6)) ||
           (hit.c1 == 36 && (hit.c2 == 15)) ||
           (hit.c1 == 37 && (hit.c2 == 15 || hit.c2 == 14))) continue;
      runset.insert(hit.run);

      if ( !(hit.fitsuccess) ) continue;
      if ( (hit.fitsumsqrresidual/hit.fitndf) > 3 ) continue;
      //if ( hit.fitconstant<-49 ) continue;
      //if ( hit.fitconstant>249 ) continue;
      //if ( hit.fitndf<50 ) continue;     // basically cut on 50 hits per track
      //if ( fabs( hit.fitquadratic )>0.000199 ) continue;     // skip if fit parameter is maximal

      int nhits = hit.channel.size();
      for (int i = 0; i < nhits; ++i)
        {
          if ( hit.ismctruth[i] ) continue;
          //if ( hit.channel[i] == 566 || hit.channel[i] == 885 || hit.channel[i] == 1547 ) continue;
          //if ( hit.tpc != 1 ) continue;//if you uncomment this, then must re-consider the fitndf cut below

          //if ( (hit.tpc[i])%2 == 0 ) continue;
          //if ( hit.width[i] > 400 ) continue;
          //if ( hit.integral[i] < -500e4 ) continue;
          //if ( hit.rms[i] > 40 ) continue;
          //if ( hit.rms[i] < 10 ) continue;
          //if ( fabs( hit.baseline[i] ) > 20 ) continue;


          if (hit.channel[i] == 288 || hit.channel[i] == 399 ||
              hit.channel[i] == 400 || hit.channel[i] == 511 ||
              hit.channel[i] == 800 || hit.channel[i] == 911 ||
              hit.channel[i] == 912 || hit.channel[i] == 1023 ||
              hit.channel[i] == 1312 || hit.channel[i] == 1423 ||
              hit.channel[i] == 1424 || hit.channel[i] == 1535 ||
              hit.channel[i] == 1824 || hit.channel[i] == 1935 ||
              hit.channel[i] == 1936 || hit.channel[i] == 2047) continue;



          runnum = hit.run;
          event = hit.event;
          hitt = hit.hitt[i];
          binnum = getBinNumber( hitt, timeBinHist );
          foundrealhit = hit.fitrealhit[i];
          assumedhit = hit.assumedhit[i];
          Double_t hitintegral = hit.integral[i];     // / mcscale;
          dqdx =  hitintegral / hit.segmentlength[i];
          channel = hit.channel[i];
          tpc = hit.tpc[i];
          ismctruth = false;

          if ( foundrealhit && !assumedhit )
            {
              found->Fill( dqdx );
              fabinmap[binnum].first = 1 + fabinmap[binnum].first;
            }
          if ( assumedhit )
            {
              assumed->Fill( dqdx );
              fabinmap[binnum].second = 1 + fabinmap[binnum].second;
            }
          if ( assumedhit || foundrealhit )
            {
              HitSave hs;
              hs.hitt = hitt;
              hs.dqdx = dqdx;
              hs.asshit = assumedhit;
              hitmap[channel].push_back( hs );

              hitstree->Fill();

              if (runseventshits.find(runnum) == runseventshits.end())
                {
                  runseventshits[runnum][event] = 1;
                }
              else
                {
                  if (runseventshits[runnum].find(event) == runseventshits[runnum].end())
                    {
                      runseventshits[runnum][event] = 1;
                    }
                  else
                    {
                      runseventshits[runnum][event] += 1;
                    }
                }

            }
        }
      if (analyse_MCTruth)
        {
          int nhits = hit.channel.size();
          for (int i = 0; i < nhits; ++i)
            {
              if (!(hit.ismctruth[i])) continue;
              if (hit.channel[i] == 288 || hit.channel[i] == 399 ||
                  hit.channel[i] == 400 || hit.channel[i] == 511 ||
                  hit.channel[i] == 800 || hit.channel[i] == 911 ||
                  hit.channel[i] == 912 || hit.channel[i] == 1023 ||
                  hit.channel[i] == 1312 || hit.channel[i] == 1423 ||
                  hit.channel[i] == 1424 || hit.channel[i] == 1535 ||
                  hit.channel[i] == 1824 || hit.channel[i] == 1935 ||
                  hit.channel[i] == 1936 || hit.channel[i] == 2047) continue;
              HitSave hs;
              hs.hitt = (hit.hitt[i]-64)*0.5 - hit.t0/1000;
              hs.dqdx = (hit.integral[i])/hit.segmentlength[i];
              hs.asshit = hit.assumedhit[i];
              hitmapMCT[hit.channel[i]].push_back( hs );
              runnum = hit.run;
              event = hit.event;
              hitt = hs.hitt;
              binnum = getBinNumber( hitt, timeBinHist );
              foundrealhit = true;
              assumedhit = false;
              ismctruth = true;
              dqdx = hs.dqdx;
              channel = hit.channel[i];
              tpc = hit.tpc[i];
              hitstree->Fill();
            }
        }
    }
  hitstree->Write();

  std::ofstream runsfile;
  runsfile.open("PurityRuns.txt");
  for (auto it = runset.begin(); it != runset.end(); ++it)
    {
      runsfile << *it << "\n";
    }
  runsfile.close();

  TCanvas * canvFAratio = new TCanvas( "canvFAratio", "FAratio vs drift", 2000, 1600 );
  TPad *pad = new TPad( "pad", "", 0, 0, 1, 1 );
  pad->SetFillColor( 0 );
  pad->Draw();
  pad->cd();
  TH1F *hr = pad->DrawFrame( 0, 0, drifttimemax, 1 );
  hr->SetXTitle( "Drift Time [ #mus]" );
  hr->SetYTitle( "#Assumed / ( #Assumed + #Found )" );
  pad->GetFrame()->SetFillColor( 0 );
  pad->GetFrame()->SetBorderSize( 12 );
  TGraph * FAratio = new TGraph();
  TGraph * sumFA = new TGraph();
  Double_t maxsum = 0, minsum = 99999;
  for ( UInt_t b = 0; b < Nbins; ++b )
    {
      std::cout << "fabinmap[" << b << "].first=" << fabinmap[b].first << "  second=" << fabinmap[b].second << std::endl;
      if ( fabinmap[b].first != 0 )
        {
          FAratio->SetPoint( b, timeBinHist[b].bincenter, ( double )fabinmap[b].second / ( ( double )fabinmap[b].second+( double )fabinmap[b].first ) );
          Double_t sum = ( double )fabinmap[b].second+( double )fabinmap[b].first;
          sumFA->SetPoint( b, timeBinHist[b].bincenter, sum );
          if ( sum > maxsum ) maxsum = sum;
          if ( sum < minsum ) minsum = sum;
        }
    }
  FAratio->Draw( "lp" );
  //FAratio->GetXaxis()->SetTitle( "Drift Time ( us )" );
  //FAratio->GetYaxis()->SetTitle( "# Assumed / ( # Assumed + # Found )" );
  FAratio->SetMarkerStyle( 21 );
  FAratio->SetMarkerSize( 2 );
  //FAratio->SetMarkerColor( 2 );
  FAratio->SetLineWidth( 2 );
  canvFAratio->cd();
  TPad *overlay = new TPad( "overlay", "", 0, 0, 1, 1 );
  overlay->SetFillStyle( 4000 );
  overlay->SetFillColor( 0 );
  overlay->SetFrameFillStyle( 4000 );
  overlay->Draw();
  overlay->cd();
  Double_t xmin1 = 0;
  Double_t ymin1 = minsum*0.9;
  Double_t xmax1 = drifttimemax;
  Double_t ymax1 = maxsum*1.1;
  TH1F *hframe = overlay->DrawFrame( xmin1, ymin1, xmax1, ymax1 );
  hframe->GetXaxis()->SetLabelOffset( 99 );
  hframe->GetYaxis()->SetLabelOffset( 99 );
  sumFA->Draw( "lp" );
  sumFA->SetMarkerStyle( 21 );
  sumFA->SetMarkerSize( 2 );
  sumFA->SetMarkerColor( 4 );
  sumFA->SetLineWidth( 2 );
  TGaxis *axis1 = new TGaxis( xmax1, ymin1, xmax1, ymax1, ymin1, ymax1, 510, "+L" );
  axis1->SetTitle( "Total # of hits" );
  axis1->SetTitleOffset( 1.1 );
  axis1->SetLineColor( 4 );
  axis1->SetLabelColor( 4 );
  axis1->Draw();
  canvFAratio->Update();
  canvFAratio->Write();
  canvFAratio->SaveAs( "canvFAratio.png" );
  delete canvFAratio;

  gStyle->SetOptStat( 11 );
  TCanvas * canvtest = new TCanvas( "fahist", "Found and Assumed Hits", 2000, 1600 );
  TPad * pad1 = new TPad( "pad1", "", 0, 0, 1, 1 );
  TPad * pad2 = new TPad( "pad2", "", 0, 0, 1, 1 );
  pad2->SetFillStyle( 4000 );
  pad1->Draw();
  pad1->cd();
  Double_t maxval = std::max( found->GetBinContent( found->GetMaximumBin() ), assumed->GetBinContent( assumed->GetMaximumBin() ) );
  found->SetLineColor( kRed );
  found->SetLineWidth( 2 );
  found->SetAxisRange( 0, maxval*1.1, "Y" );
  found->Draw();
  found->GetXaxis()->SetTitle( "Charge" );
  pad1->Update();
  TPaveStats * ps1 = ( TPaveStats* )found->GetListOfFunctions()->FindObject( "stats" );
  if (ps1)
    {
      ps1->SetX1NDC( 0.4 );
      ps1->SetX2NDC( 0.6 );
      ps1->SetTextColor( kRed );
    }
  pad1->Modified();
  canvtest->cd();
  Double_t ymin = 0;
  Double_t ymax = maxval*1.1;
  Double_t dy = ( ymax-ymin )/0.8;
  Double_t xmin = ADCcutofflow;
  Double_t xmax = ADCcutoffhigh;
  Double_t dx = ( xmax-xmin )/0.8;
  pad2->Range( xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy );
  pad2->Draw();
  pad2->cd();
  assumed->SetLineColor( kBlue );
  assumed->SetLineWidth( 2 );
  assumed->Draw( "][sames" );
  assumed->GetXaxis()->SetTitle( "Charge" );
  pad2->Update();
  TPaveStats * ps2 = ( TPaveStats* )assumed->GetListOfFunctions()->FindObject( "stats" );
  if (ps2)
    {
      ps2->SetX1NDC( 0.65 );
      ps2->SetX2NDC( 0.85 );
      ps2->SetTextColor( kBlue );
    }
  TGaxis * axis = new TGaxis( xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L" );
  axis->SetLabelFont( 42 );
  axis->SetLabelSize( 0.035 );
  axis->Draw();
  TLegend * leg = new TLegend( 0.55, 0.6, 0.85, 0.75 );
  leg->AddEntry( found, "Found Hits", "l" );
  leg->AddEntry( assumed, "Assumed Hits", "l" );
  leg->Draw();
  canvtest->Update();
  canvtest->Write();
  canvtest->SaveAs( "fahist.png" );
  delete canvtest;


  for ( auto & i_chan : hitmap )
    {
      for ( auto & i_hit : i_chan.second )
        {
          if ( i_hit.dqdx < ADCcutoffhigh && i_hit.dqdx > ADCcutofflow )
            {
              UInt_t ibin = 0;
              for (UInt_t i_hist = 0; i_hist < timeBinHist.size(); i_hist++ )
                {
                  if (i_hit.hitt >= timeBinHist[i_hist].binlow && i_hit.hitt < timeBinHist[i_hist].binhigh) {
                      ibin = i_hist;
                      break;
                    }
                }
              Float_t dqdx = i_hit.dqdx;
              charge = dqdx;
              tdrift = ( Float_t )( i_hit.hitt );
              chgdata->add( RooArgSet( tdrift, charge ) );
              timeBinHist[ibin].hist->Fill(dqdx);
            }
        }
    }

  for ( auto & i_chan : hitmapMCT )
    {
      for ( auto & i_hit : i_chan.second )
        {
          if ( i_hit.dqdx < ADCcutoffhighMCT && i_hit.dqdx > ADCcutofflowMCT )
            {
              UInt_t ibin = 0;
              for (UInt_t i_hist = 0; i_hist < timeBinHist.size(); i_hist++ )
                {
                  if (i_hit.hitt >= timeBinHist[i_hist].binlow && i_hit.hitt < timeBinHist[i_hist].binhigh) {
                      ibin = i_hist;
                      break;
                    }
                }
              Float_t dqdx = i_hit.dqdx;
              chargeMCT = dqdx;
              tdriftMCT = ( Float_t )( i_hit.hitt );
              chgdataMCT->add( RooArgSet( tdriftMCT, chargeMCT ) );
              timeBinHist[ibin].histMCT->Fill(dqdx);
            }
        }
    }

  chgdata->Print();
  chgdataMCT->Print();

  std::vector<FitResult> results( Nbins );
  std::vector<FitResult> resultsMCT( Nbins );

  Double_t mlscaleVal = 1000;
  Double_t slscaleVal = 100;
  Double_t sgscaleVal = 100;
  RooRealVar mpvshift("MPVShift","mpv shift landau",0.-22278298);
  RooRealVar one("one","one",1.0);
  RooRealVar mlvar( "LandMean/1000", "mean landau var", 2000/mlscaleVal, ADCcutofflow/mlscaleVal, ADCcutoffhigh/mlscaleVal );
  RooRealVar mlscale("LandMPVScale","scale mean landau", mlscaleVal);
  RooProduct ml("LandMean","mlvar*mlscale",RooArgSet(mlvar,mlscale));
  RooRealVar slvar( "LandWidth/100", "sigma landau var", 100/slscaleVal, 10/slscaleVal, 2000/slscaleVal );
  RooRealVar slscale("LandWidthScale","scale sigma landau", slscaleVal);
  RooProduct sl("LandWidth","slvar*slscale",RooArgSet(slvar,slscale));
  ////RooRealVar sl("LandWidth","sigma landau",608);
  //RooRealVar smlratio("LandWidthMPV","sigma/mpv",0.2059);
  //RooProduct sl("LandWidth","smlratio*ml",RooArgSet(smlratio,ml));
  ///////////////RooFormulaVar mlm("LandMPV","@0+@1*@2",RooArgSet(ml,mpvshift,sl));
  RooLandau landau( "lx", "lx", charge, ml, sl );
  RooRealVar mg( "mg", "mean gauss", 0 );
  RooRealVar sgvar( "GaussWidth/100", "sigma gauss var", 100/sgscaleVal, 10/sgscaleVal, 2000/sgscaleVal );
  RooRealVar sgscale( "GaussWidthScale","scale sigma gauss",sgscaleVal);
  RooProduct sg("GaussWidth","sgvar*sgscale",RooArgSet(sgvar,sgscale));
  RooGaussian gauss( "gauss", "gauss", charge, mg, sg );
  RooFFTConvPdf lxg( "lxg", "landau ( x ) gauss", charge, landau, gauss );
  RooPolynomial pol0_2( "pol0_2", "pol0", charge, RooArgList() );
  RooRealVar lxg_yield( "BkgCoef", "BkgCoef", 0 );
  RooAddPdf clxg( "clxg", "const + ( landau ( x ) gauss )", RooArgList( pol0_2, lxg ), RooArgList( lxg_yield ) );
  RooFormulaVar ratio( "ratio","@0 / @1",RooArgList(sl,ml));
  RooGaussian ratioconstraint("ratioconstraint","ratioconstraint",ratio,RooFit::RooConst(0.22),RooFit::RooConst(0.5));

  RooRealVar mlvarMCT( "LandMeanMCT/1000", "mean landau var MCT", 2000/mlscaleVal, ADCcutofflowMCT/mlscaleVal, ADCcutoffhighMCT/mlscaleVal );
  RooProduct mlMCT("LandMeanMCT","mlvarMCT*mlscale",RooArgSet(mlvarMCT,mlscale));
  RooRealVar slvarMCT( "LandWidthMCT/100", "sigma landau var MCT", 100/slscaleVal, 10/slscaleVal, 2000/slscaleVal );
  RooProduct slMCT("LandWidthMCT","slvarMCT*slscale",RooArgSet(slvarMCT,slscale));
  RooLandau landauMCT( "lxMCT", "lxMCT", chargeMCT, mlMCT, slMCT );
  RooRealVar sgvarMCT( "GaussWidthMCT/100", "sigma gauss var MCT", 100/sgscaleVal, 10/sgscaleVal, 2000/sgscaleVal );
  RooProduct sgMCT("GaussWidthMCT","sgvarMCT*sgscale",RooArgSet(sgvarMCT,sgscale));
  RooGaussian gaussMCT( "gaussMCT", "gaussMCT", chargeMCT, mg, sgMCT );
  RooFFTConvPdf lxgMCT( "lxgMCT", "landau ( x ) gauss MCT", chargeMCT, landauMCT, gaussMCT );


  RooRealVar pgmean( "PeakGaussMean", "peak gauss mean", 2800, -1000, 200000 );
  RooRealVar pgwidth( "PeakGaussWidth", "peak gauss width", 1000, 1, 200000 );
  RooGaussian pg( "PeakGaussian", "peak gaussian", charge, pgmean, pgwidth );
  RooPolynomial pol0_3( "pol0_3", "pol0", charge, RooArgList() );
  RooRealVar pg_yield( "BkgCoef", "BkgCoef", 0.5, 0, 1 );
  RooAddPdf cpg( "cpg", "const + peak gaussian", RooArgList( pol0_3, pg ), RooArgList( pg_yield ) );

  TGraphErrors * hbin_lxg = new TGraphErrors( Nbins );
  hbin_lxg->SetTitle( "Most Probable dQ/dx in L( x )g Fit" );
  Double_t minMPVlxg = 99999;
  Double_t maxMPVlxg = -99999;

  TGraphErrors * hbin_lxgMCT = new TGraphErrors( Nbins );
  hbin_lxgMCT->SetTitle( "Most Probable dQ/dx in L( x )g Fit, MCT" );
  Double_t minMPVlxgMCT = 99999;
  Double_t maxMPVlxgMCT = -99999;

  TGraphErrors * hbin_pg = new TGraphErrors( Nbins );
  hbin_pg->SetTitle( "Peak dQ/dx in PeakGauss Fit" );
  Double_t minMPVpg = 99999;
  Double_t maxMPVpg = -99999;

  TGraphErrors * gw_lxg = new TGraphErrors( Nbins );
  gw_lxg->SetTitle( "Width of Gaussian in L( x )g Convolution" );
  TGraphErrors * lw = new TGraphErrors( Nbins );
  lw->SetTitle( "Width of Landau in L( x )g Convolution" );
  TGraphErrors * gw_pg = new TGraphErrors( Nbins );
  gw_pg->SetTitle( "Width of Gaussian Fit to Peak" );

  TGraphErrors * landwidmpv = new TGraphErrors(Nbins);
  landwidmpv->SetTitle("Landau Width / Landau MPV");

  TGraph * lxgchi2 = new TGraph( Nbins );
  lxgchi2->SetTitle( "Chi2/NDF of L( x )g Best Fit" );

  TGraphErrors * gw_lxgMCT = new TGraphErrors( Nbins );
  gw_lxgMCT->SetTitle( "Width of Gaussian in L( x )g Convolution, MCT" );
  TGraphErrors * lwMCT = new TGraphErrors( Nbins );
  lwMCT->SetTitle( "Width of Landau in L( x )g Convolution, MCT" );

  TGraphErrors * landwidmpvMCT = new TGraphErrors(Nbins);
  landwidmpvMCT->SetTitle("Landau Width / Landau MPV, MCT");

  TGraph * lxgchi2MCT = new TGraph( Nbins );
  lxgchi2MCT->SetTitle( "Chi2/NDF of L( x )g Best Fit, MCT" );

  for ( UInt_t i = 0; i < Nbins; i++ )
    {
      if ( timeBinHist[i].hist->GetEntries() > 0 )
        {
          char shortname[100];
          sprintf( shortname, "dQdx_bin%d", i );

          char cut[100];
          sprintf( cut, "tdrift > %f && tdrift < %f", timeBinHist[i].binlow, timeBinHist[i].binhigh );
          RooDataSet* roodata = ( RooDataSet* ) chgdata->reduce( RooFit::SelectVars( RooArgSet(charge) ), RooFit::Cut( cut ), RooFit::Name( shortname ), RooFit::Title( shortname ) );

          TH1D * clonehist = (TH1D*)roodata->createHistogram(TString::Format("roodata_clone_%s",timeBinHist[i].hist->GetName()),charge,RooFit::Binning(200,ADCcutofflow,ADCcutoffhigh));

          clonehist->Write();

          Float_t maxbin = clonehist->GetBinCenter( clonehist->GetMaximumBin() );
          //Float_t maxbin = 2000;

          Float_t maxvalue = clonehist->GetBinContent( clonehist->GetMaximumBin() );
          Float_t cutoffheightlow = maxvalue*0.66;
          Float_t cutoffheighthigh = maxvalue*0.66;
          Float_t lowval = clonehist->GetBinCenter( clonehist->FindFirstBinAbove( cutoffheightlow ) );
          Float_t highval = clonehist->GetBinCenter( clonehist->FindLastBinAbove( cutoffheighthigh ) );

          std::cout << "maxbin=" << maxbin << "  maxvalue=" << maxvalue << "  cutoffheight=( " << cutoffheightlow << ", " << cutoffheighthigh << " )  lowval=" << lowval << "  highval=" << highval << std::endl;

          pgmean.setVal( maxbin );
          pgmean.setMin( lowval );
          pgmean.setMax( highval );
          pgwidth.setVal( 1000 );
          pgwidth.setMin( 1 );
          pgwidth.setMax( 10000 );


          Float_t lxgfitlow = maxvalue*0.4;
          Float_t lxgfithigh = maxvalue*0.4;
          Float_t lxglowval = clonehist->GetBinCenter( clonehist->FindFirstBinAbove( lxgfitlow ) );
          Float_t lxghighval = clonehist->GetBinCenter( clonehist->FindLastBinAbove( lxgfithigh ) );

          // Setup component pdfs
          // --------------------

          //setup observable
          //charge.setRange( shortname, maxbin-2000, maxbin+3000 );
          charge.setRange( shortname, ADCcutofflow, ADCcutoffhigh );

          //setup landau( t, ml, sl )
          mlvar.setVal( maxbin/mlscaleVal );
          //mlvar.setMin( (maxbin-400)/mlscaleVal ); //0 );
          //mlvar.setMax( (maxbin+400)/mlscaleVal ); //0 );
          //ml.setVal(maxbin);
          //ml.setMin(maxbin-2000);
          //ml.setMax(maxbin+2000);
          slvar.setVal( ((i>0) ? results[i-1].landwidth : 600)/slscaleVal );
          //sl.setMin( 10 );
          //sl.setMax( 2000 );

          //setup gauss( t, mg, sg )
          //mg.setVal( 1 );
          //mg.setMin( -100 );
          //mg.setMax( 100 );
          sgvar.setVal( ((i>0) ? results[i-1].gausswidth_l : 500)/slscaleVal );
          //sg.setMin( 10 );
          //sg.setMax( 2000 );

          // construct convolution pdf
          // -------------------------

          // set num bins to be used for FFT sampling
          charge.setBins( 20000, "cache" );

          // fit convoluted pdf to binHist
          // -----------------------------



          if ( DoLxG )
            {
              RooFitResult * fr = lxg.fitTo( *roodata,/* RooFit::ExternalConstraints(ratioconstraint),*/ RooFit::Save(), RooFit::PrintLevel( 1 ), RooFit::Range( lxglowval, lxghighval ), RooFit::Minimizer("Minuit2"), RooFit::Minos(kTRUE), RooFit::Offset(kTRUE)); //, RooFit::SumW2Error(kTRUE));
              RooPlot * lxgframe = charge.frame( RooFit::Title("" /*TString::Format( "%s L( x )g", timeBinHist[i].hist->GetTitle()) */) );
              roodata->plotOn( lxgframe, RooFit::Binning( 200 ), RooFit::DataError(RooAbsData::SumW2), RooFit::LineColor(kBlack), RooFit::MarkerColor(kBlack) );
              lxg.plotOn( lxgframe, RooFit::Name( "lxg" ), RooFit::LineColor( kRed ), RooFit::Range( lxglowval, lxghighval ) );
              results[i].lxg_chi2ndf = lxgframe->chiSquare();

              std::cout << "status = " << fr->status() << "  covQual = " << fr->covQual() << "  numInvalidNLL = " << fr->numInvalidNLL() << "  EDM = " << fr->edm() << "  minNll = " << fr->minNll() << std::endl;

              std::cout << "Correlation Matrix:" << std::endl;
              fr->correlationMatrix().Print();
              std::cout << "Covariance Matrix:" << std::endl;
              fr->covarianceMatrix().Print();


              TPaveLabel * tlxg = new TPaveLabel( 0.7, 0.83, 0.99, 0.9, Form( " #chi^{2}/NDF = %f", lxgframe->chiSquare() ), "NDC" );
              tlxg->SetFillColor( 0 );
              tlxg->SetBorderSize( 1 );
              tlxg->SetTextAlign( 12 );
              lxgframe->addObject( tlxg );
              lxgframe->getAttText()->SetTextSize( 0.30 );
              roodata->statOn( lxgframe, RooFit::Layout( 0.7, 0.99, 0.8 ) );
              lxgframe->getAttText()->SetTextSize( 0.019 );
              lxg.paramOn( lxgframe, RooFit::Layout( 0.7, 0.99, 0.6 ) );
              lxgframe->getAttText()->SetTextSize( 0.019 );
              TPaveLabel * tlxgmpv = new TPaveLabel(0.7, 0.3,0.99,0.4,Form("LandMPV=%.1f +/- %.1f",ml.getVal()-0.22278298*sl.getVal(),sqrt(pow(mlvar.getError()*mlscaleVal,2)+pow(0.22278298,2)*pow(slvar.getError()*slscaleVal,2)+2*(-0.22278298)*fr->covarianceMatrix()[1][2])),"NDC");
              tlxgmpv->SetFillColor(0);
              tlxgmpv->SetBorderSize(1);
              tlxgmpv->SetTextAlign(12);
              lxgframe->addObject(tlxgmpv);
              //lxgframe->getAttText()->SetTextSize(0.2);

              char buff[100];
              sprintf( buff, "lxg_bin%d", i );
              TCanvas* canvlxg2 = new TCanvas( buff, buff, 1600, 1600 );
              gPad->SetLeftMargin( 0.15 );
              lxgframe->GetYaxis()->SetTitleOffset( 1.6 );
              lxgframe->SetTitle("");
              lxgframe->Draw();
              canvlxg2->Write();
              canvlxg2->SaveAs( TString::Format( "%s.png", buff ) );

              if ( DoResidual )
                {
                  char buffresid[100];
                  sprintf( buffresid, "bin%d_resid_lxg", i );
                  TCanvas* cresidlxg = new TCanvas( buffresid, buffresid, 1600, 1600 );
                  gPad->SetLeftMargin( 0.15 );
                  RooPlot* frameresidlxg = charge.frame( RooFit::Title( "lxg_residual" ) );
                  RooHist* hresidlxg = lxgframe->residHist();
                  frameresidlxg->addPlotable( hresidlxg, "P" );
                  frameresidlxg->GetYaxis()->SetTitleOffset( 1.4 );
                  frameresidlxg->Draw();
                  cresidlxg->Write();
                  gSystem->ProcessEvents();

                  delete frameresidlxg;
                  delete cresidlxg;
                }

              if ( DoPull )
                {
                  char buffpull[100];
                  sprintf( buffpull, "bin%d_pull_lxg", i );
                  TCanvas* cpulllxg = new TCanvas( buffpull, buffpull, 1600, 1600 );
                  gPad->SetLeftMargin( 0.15 );
                  RooPlot* framepulllxg = charge.frame( RooFit::Title( "lxg_pull" ) );
                  RooHist* hpull2lxg = lxgframe->pullHist();
                  TH1F* hpulllxg = new TH1F( "hpull_lxg", "pull_lxg", 100, 0, 0 );
                  for ( Int_t ibin=0; ibin<hpull2lxg->GetN(); ibin++ )
                    {
                      Double_t x, y;
                      hpull2lxg->GetPoint( ibin, x, y );
                      hpulllxg->Fill( y );
                    }
                  hpulllxg->Draw();
                  cpulllxg->Write();
                  gSystem->ProcessEvents();

                  delete cpulllxg;
                  delete hpulllxg;
                  delete framepulllxg;
                }

              //double correctionfactor,binxtemp;
              //datamc->GetPoint(i,binxtemp,correctionfactor);

              std::cout << "!!!!!====== slvar.getError() = " << slvar.getError() << std::endl;

              results[i].landmpv = ml.getVal()-0.22278298*sl.getVal(); ///correctionfactor; ///mcscale;
              results[i].landwidth = sl.getVal(); ///mcscale;
              results[i].gaussmean_l = mg.getVal(); ///mcscale;
              results[i].gausswidth_l = sg.getVal(); ///mcscale;
              //results[i].landmpverr = (mlvar.getError() /* / mcscale*/ )*mlscaleVal;
              results[i].landmpverr = sqrt(pow(mlvar.getError()*mlscaleVal,2)+pow(0.22278298,2)*pow(slvar.getError()*slscaleVal,2)+2*(-0.22278298)*fr->covarianceMatrix()[1][2]);
              results[i].landwidtherr = (slvar.getError() /* / mcscale */ )*slscaleVal;
              results[i].gaussmean_lerr = mg.getError(); ///mcscale;
              results[i].gausswidth_lerr = (sgvar.getError() /* / mcscale */ )*sgscaleVal;
              results[i].lxglow = lxglowval; ///mcscale;
              results[i].lxghigh = lxghighval; ///mcscale;
              //results[i].lxg_chi2ndf = lxgframe->chiSquare( "lxg", shortname );

              std::cout << "MPV=" << ml.getVal() << "  mlvar.getError()=" << mlvar.getError() << std::endl;

              if ( minMPVlxg > results[i].landmpv ) minMPVlxg = results[i].landmpv;
              if ( maxMPVlxg < results[i].landmpv ) maxMPVlxg = results[i].landmpv;

              results[i].fitsuccess = false;
              if ( fr->status() == 0 && fr->covQual() == 3 && fr->edm() < 1 && results[i].lxg_chi2ndf < 3)
                {
                  results[i].fitsuccess = true;
                  hbin_lxg->SetPoint( i, timeBinHist[i].bincenter, results[i].landmpv );
                  hbin_lxg->SetPointError( i, 0, results[i].landmpverr );
                }
              gw_lxg->SetPoint( i, timeBinHist[i].bincenter, results[i].gausswidth_l );
              gw_lxg->SetPointError( i, 0, results[i].gausswidth_lerr );
              lw->SetPoint( i, timeBinHist[i].bincenter, results[i].landwidth );
              lw->SetPointError( i, 0, results[i].landwidtherr );

              double ratiowidmpv = results[i].landwidth / results[i].landmpv;
              double ratiowiderr = results[i].landwidtherr / results[i].landwidth;
              double ratiompverr = results[i].landmpverr / results[i].landmpv;
              landwidmpv->SetPoint( i, timeBinHist[i].bincenter, ratiowidmpv);
              landwidmpv->SetPointError( i, 0, ratiowidmpv*sqrt(ratiowiderr*ratiowiderr+ratiompverr*ratiompverr));
              std::cout << "WIDTHDIVMPV = " << ratiowidmpv << "  WIDTHDIVMPVERR = " << ratiowidmpv*sqrt(ratiowiderr*ratiowiderr+ratiompverr*ratiompverr) << std::endl;

              lxgchi2->SetPoint( i, timeBinHist[i].bincenter, results[i].lxg_chi2ndf );

              delete canvlxg2;
              //delete imglxg;
              delete lxgframe;
            }

          if ( DoPG )
            {
              cpg.fitTo( *roodata, RooFit::Save(), RooFit::PrintLevel( -1 ), RooFit::Range( lowval, highval ), RooFit::SumW2Error( true ) );
              RooPlot * pgframe = charge.frame( RooFit::Title( TString::Format( "%s PeakGaussian", timeBinHist[i].hist->GetTitle() ) ) );
              roodata->plotOn( pgframe, RooFit::Binning( 200 ) );
              cpg.plotOn( pgframe, RooFit::Name( "pgb" ), RooFit::LineColor( kBlue ), RooFit::Range( lowval, highval ) );
              //results[i].pg_chi2ndf = pgframe->chiSquare();

              TPaveLabel * tpg = new TPaveLabel( 0.7, 0.83, 0.99, 0.9, Form( " #chi^{2}/NDF = %f", pgframe->chiSquare() ), "NDC" );
              tpg->SetFillColor( 0 );
              tpg->SetBorderSize( 1 );
              tpg->SetTextAlign( 12 );
              tpg->SetTextSize( 0.3 );
              pgframe->addObject( tpg );
              pgframe->getAttText()->SetTextSize( 0.3 );
              roodata->statOn( pgframe, RooFit::Layout( 0.7, 0.99, 0.8 ) );
              pgframe->getAttText()->SetTextSize( 0.019 );
              cpg.paramOn( pgframe, RooFit::Layout( 0.7, 0.99, 0.6 ) );
              pgframe->getAttText()->SetTextSize( 0.019 );

              char buff[100];
              sprintf( buff, "pg_bin%d", i );
              TCanvas* canvpg2 = new TCanvas( buff, buff, 1600, 1600 );
              gPad->SetLeftMargin( 0.15 );
              pgframe->GetYaxis()->SetTitleOffset( 1.4 );
              pgframe->Draw();
              canvpg2->Write();
              //canvpg2->SaveAs( TString::Format( "%s.png", buff ) );

              //gSystem->ProcessEvents();
              //TImage *imgpg = TImage::Create();
              //imgpg->FromPad( canvpg2 );
              //char newbuff[100];
              //sprintf( newbuff, "%s.png", buff );
              //imgpg->WriteImage( newbuff );


              if ( DoResidual )
                {
                  char buffresid[100];
                  sprintf( buffresid, "bin%d_resid_pg", i );
                  TCanvas* cresidpg = new TCanvas( buffresid, buffresid, 1600, 1600 );
                  gPad->SetLeftMargin( 0.15 );
                  RooPlot* frameresidpg = charge.frame( RooFit::Title( "pg_residual" ) );
                  RooHist* hresidpg = pgframe->residHist();
                  frameresidpg->addPlotable( hresidpg, "P" );
                  frameresidpg->GetYaxis()->SetTitleOffset( 1.4 );
                  frameresidpg->Draw();
                  cresidpg->Write();
                  gSystem->ProcessEvents();

                  delete cresidpg;
                  delete frameresidpg;

                }

              if ( DoPull )
                {
                  char buffpull[100];
                  sprintf( buffpull, "bin%d_pull_pg", i );
                  TCanvas* cpullpg = new TCanvas( buffpull, buffpull, 1600, 1600 );
                  gPad->SetLeftMargin( 0.15 );
                  RooPlot* framepullpg = charge.frame( RooFit::Title( "pg_pull" ) );
                  RooHist* hpull2pg = pgframe->pullHist();
                  TH1F* hpullpg = new TH1F( "hpull_pg", "pull_pg", 100, 0, 0 );
                  for ( Int_t ibin=0; ibin<hpull2pg->GetN(); ibin++ )
                    {
                      Double_t x, y;
                      hpull2pg->GetPoint( ibin, x, y );
                      hpullpg->Fill( y );
                    }
                  hpullpg->Draw();
                  cpullpg->Write();
                  gSystem->ProcessEvents();

                  delete cpullpg;
                  delete hpullpg;
                  delete framepullpg;
                }

              results[i].gaussmean = pgmean.getVal();
              results[i].gausswidth = pgwidth.getVal();
              results[i].gaussmeanerr = pgmean.getError();
              results[i].gausswidtherr = pgwidth.getError();
              results[i].pg_chi2ndf = pgframe->chiSquare( "cpg", shortname );

              if ( minMPVpg > results[i].gaussmean ) minMPVpg = results[i].gaussmean;
              if ( maxMPVpg < results[i].gaussmean ) maxMPVpg = results[i].gaussmean;

              //if ( results[i].pg_chi2ndf < 2 )
              {
                hbin_pg->SetPoint( i, timeBinHist[i].bincenter, results[i].gaussmean );
                hbin_pg->SetPointError( i, 0, results[i].gaussmeanerr );
              }
              gw_pg->SetPoint( i, timeBinHist[i].bincenter, results[i].gausswidth );
              gw_pg->SetPointError( i, 0, results[i].gausswidtherr );

              delete canvpg2;
              delete pgframe;
              //delete imgpg;
            }

          results[i].binnum = i;
          results[i].bincenter = timeBinHist[i].bincenter;
          results[i].binwidtherr = ( timeBinHist[i].binwidth )/sqrt( 12.0 );

        }
      if ( timeBinHist[i].histMCT->GetEntries() > 0 && analyse_MCTruth)
        {
          char shortname[100];
          sprintf( shortname, "dQdx_bin%dMCT", i );

          char cut[100];
          sprintf( cut, "tdriftMCT > %f && tdriftMCT < %f", timeBinHist[i].binlow, timeBinHist[i].binhigh );
          RooDataSet* roodataMCT = ( RooDataSet* ) chgdataMCT->reduce( RooFit::SelectVars( RooArgSet(chargeMCT) ), RooFit::Cut( cut ), RooFit::Name( shortname ), RooFit::Title( shortname ) );

          TH1D * clonehistMCT = (TH1D*)roodataMCT->createHistogram(TString::Format("roodata_clone_%sMCT",timeBinHist[i].histMCT->GetName()),chargeMCT,RooFit::Binning(200,ADCcutofflowMCT,ADCcutoffhighMCT));

          clonehistMCT->Write();

          Float_t maxbin = clonehistMCT->GetBinCenter( clonehistMCT->GetMaximumBin() );
          Float_t maxvalue = clonehistMCT->GetBinContent( clonehistMCT->GetMaximumBin() );
          Float_t lxgfitlow = maxvalue*0.2;
          Float_t lxgfithigh = maxvalue*0.3;
          Float_t lxglowval = clonehistMCT->GetBinCenter( clonehistMCT->FindFirstBinAbove( lxgfitlow ) );
          Float_t lxghighval = clonehistMCT->GetBinCenter( clonehistMCT->FindLastBinAbove( lxgfithigh ) );

          // Setup component pdfs
          // --------------------

          //setup observable
          //charge.setRange( shortname, maxbin-2000, maxbin+3000 );
          chargeMCT.setRange( shortname, ADCcutofflowMCT, ADCcutoffhighMCT );

          //setup landau( t, ml, sl )
          mlvarMCT.setVal( maxbin/mlscaleVal );
          //mlvar.setMin( (maxbin-400)/mlscaleVal ); //0 );
          //mlvar.setMax( (maxbin+400)/mlscaleVal ); //0 );
          //ml.setVal(maxbin);
          //ml.setMin(maxbin-2000);
          //ml.setMax(maxbin+2000);
          slvarMCT.setVal( ((i>0) ? resultsMCT[i-1].landwidth : 600)/slscaleVal );
          //sl.setMin( 10 );
          //sl.setMax( 2000 );

          //setup gauss( t, mg, sg )
          //mg.setVal( 1 );
          //mg.setMin( -100 );
          //mg.setMax( 100 );
          sgvarMCT.setVal( ((i>0) ? resultsMCT[i-1].gausswidth_l : 500)/slscaleVal );
          //sg.setMin( 10 );
          //sg.setMax( 2000 );

          // construct convolution pdf
          // -------------------------

          // set num bins to be used for FFT sampling
          chargeMCT.setBins( 20000, "cache" );

          // fit convoluted pdf to binHist
          // -----------------------------



          if ( DoLxG )
            {
              RooFitResult * frMCT = lxgMCT.fitTo( *roodataMCT,/* RooFit::ExternalConstraints(ratioconstraint),*/ RooFit::Save(), RooFit::PrintLevel( 1 ), RooFit::Range( lxglowval, lxghighval ), RooFit::Minimizer("Minuit2"), RooFit::Minos(kTRUE), RooFit::Offset(kTRUE));   //, RooFit::SumW2Error(kTRUE));
              RooPlot * lxgframeMCT = chargeMCT.frame( RooFit::Title("" /*TString::Format( "%s L( x )g", timeBinHist[i].hist->GetTitle()) */) );
              roodataMCT->plotOn( lxgframeMCT, RooFit::Binning( 200 ), RooFit::DataError(RooAbsData::SumW2), RooFit::LineColor(kBlack), RooFit::MarkerColor(kBlack) );
              lxgMCT.plotOn( lxgframeMCT, RooFit::Name( "lxgMCT" ), RooFit::LineColor( kRed ), RooFit::Range( lxglowval, lxghighval ) );
              resultsMCT[i].lxg_chi2ndf = lxgframeMCT->chiSquare();

              std::cout << "status = " << frMCT->status() << "  covQual = " << frMCT->covQual() << "  numInvalidNLL = " << frMCT->numInvalidNLL() << "  EDM = " << frMCT->edm() << "  minNll = " << frMCT->minNll() << std::endl;

              std::cout << "Correlation Matrix:" << std::endl;
              frMCT->correlationMatrix().Print();
              std::cout << "Covariance Matrix:" << std::endl;
              frMCT->covarianceMatrix().Print();


              TPaveLabel * tlxgMCT = new TPaveLabel( 0.7, 0.83, 0.99, 0.9, Form( " #chi^{2}/NDF = %f", lxgframeMCT->chiSquare() ), "NDC" );
              tlxgMCT->SetFillColor( 0 );
              tlxgMCT->SetBorderSize( 1 );
              tlxgMCT->SetTextAlign( 12 );
              lxgframeMCT->addObject( tlxgMCT );
              lxgframeMCT->getAttText()->SetTextSize( 0.30 );
              roodataMCT->statOn( lxgframeMCT, RooFit::Layout( 0.7, 0.99, 0.8 ) );
              lxgframeMCT->getAttText()->SetTextSize( 0.019 );
              lxgMCT.paramOn( lxgframeMCT, RooFit::Layout( 0.7, 0.99, 0.6 ) );
              lxgframeMCT->getAttText()->SetTextSize( 0.019 );
              TPaveLabel * tlxgmpvMCT = new TPaveLabel(0.7, 0.3,0.99,0.4,Form("LandMPV=%.1f +/- %.1f",mlMCT.getVal()-0.22278298*slMCT.getVal(),sqrt(pow(mlvarMCT.getError()*mlscaleVal,2)+pow(0.22278298,2)*pow(slvarMCT.getError()*slscaleVal,2)+2*(-0.22278298)*frMCT->covarianceMatrix()[1][2])),"NDC");
              tlxgmpvMCT->SetFillColor(0);
              tlxgmpvMCT->SetBorderSize(1);
              tlxgmpvMCT->SetTextAlign(12);
              lxgframeMCT->addObject(tlxgmpvMCT);
              //lxgframe->getAttText()->SetTextSize(0.2);

              char buff[100];
              sprintf( buff, "lxg_bin%dMCT", i );
              TCanvas* canvlxg2MCT = new TCanvas( buff, buff, 1600, 1600 );
              gPad->SetLeftMargin( 0.15 );
              lxgframeMCT->GetYaxis()->SetTitleOffset( 1.6 );
              lxgframeMCT->SetTitle("");
              lxgframeMCT->Draw();
              canvlxg2MCT->Write();
              canvlxg2MCT->SaveAs( TString::Format( "%s.png", buff ) );

              if ( DoResidual )
                {
                  char buffresid[100];
                  sprintf( buffresid, "bin%d_resid_lxg", i );
                  TCanvas* cresidlxg = new TCanvas( buffresid, buffresid, 1600, 1600 );
                  gPad->SetLeftMargin( 0.15 );
                  RooPlot* frameresidlxg = chargeMCT.frame( RooFit::Title( "lxg_residual" ) );
                  RooHist* hresidlxg = lxgframeMCT->residHist();
                  frameresidlxg->addPlotable( hresidlxg, "P" );
                  frameresidlxg->GetYaxis()->SetTitleOffset( 1.4 );
                  frameresidlxg->Draw();
                  cresidlxg->Write();
                  gSystem->ProcessEvents();

                  delete frameresidlxg;
                  delete cresidlxg;
                }

              if ( DoPull )
                {
                  char buffpull[100];
                  sprintf( buffpull, "bin%d_pull_lxg", i );
                  TCanvas* cpulllxg = new TCanvas( buffpull, buffpull, 1600, 1600 );
                  gPad->SetLeftMargin( 0.15 );
                  RooPlot* framepulllxg = chargeMCT.frame( RooFit::Title( "lxg_pull" ) );
                  RooHist* hpull2lxg = lxgframeMCT->pullHist();
                  TH1F* hpulllxg = new TH1F( "hpull_lxg", "pull_lxg", 100, 0, 0 );
                  for ( Int_t ibin=0; ibin<hpull2lxg->GetN(); ibin++ )
                    {
                      Double_t x, y;
                      hpull2lxg->GetPoint( ibin, x, y );
                      hpulllxg->Fill( y );
                    }
                  hpulllxg->Draw();
                  cpulllxg->Write();
                  gSystem->ProcessEvents();

                  delete cpulllxg;
                  delete hpulllxg;
                  delete framepulllxg;
                }

              std::cout << "!!!!!====== slvarMCT.getError() = " << slvarMCT.getError() << std::endl;

              resultsMCT[i].landmpv = mlMCT.getVal()-0.22278298*slMCT.getVal();   ///correctionfactor; ///mcscale;
              resultsMCT[i].landwidth = slMCT.getVal();   ///mcscale;
              resultsMCT[i].gaussmean_l = mg.getVal();   ///mcscale;
              resultsMCT[i].gausswidth_l = sgMCT.getVal();   ///mcscale;
              //resultsMCT[i].landmpverr = (mlvarMCT.getError() /* / mcscale*/ )*mlscaleVal;
              resultsMCT[i].landmpverr = sqrt(pow(mlvarMCT.getError()*mlscaleVal,2)+pow(0.22278298,2)*pow(slvarMCT.getError()*slscaleVal,2)+2*(-0.22278298)*frMCT->covarianceMatrix()[1][2]);
              resultsMCT[i].landwidtherr = (slvarMCT.getError() /* / mcscale */ )*slscaleVal;
              resultsMCT[i].gaussmean_lerr = mg.getError();   ///mcscale;
              resultsMCT[i].gausswidth_lerr = (sgvarMCT.getError() /* / mcscale */ )*sgscaleVal;
              resultsMCT[i].lxglow = lxglowval;   ///mcscale;
              resultsMCT[i].lxghigh = lxghighval;   ///mcscale;

              std::cout << "MCT: MPV=" << mlMCT.getVal() << "  mlvarMCT.getError()=" << mlvarMCT.getError() << std::endl;

              if ( minMPVlxgMCT > resultsMCT[i].landmpv ) minMPVlxgMCT = resultsMCT[i].landmpv;
              if ( maxMPVlxgMCT < resultsMCT[i].landmpv ) maxMPVlxgMCT = resultsMCT[i].landmpv;

              resultsMCT[i].fitsuccess = false;
              if ( frMCT->status() == 0 && frMCT->covQual() == 3 && frMCT->edm() < 1 && resultsMCT[i].lxg_chi2ndf < 3)
                {
                  resultsMCT[i].fitsuccess = true;
                  hbin_lxgMCT->SetPoint( i, timeBinHist[i].bincenter, resultsMCT[i].landmpv );
                  hbin_lxgMCT->SetPointError( i, 0, resultsMCT[i].landmpverr );
                }
              gw_lxgMCT->SetPoint( i, timeBinHist[i].bincenter, resultsMCT[i].gausswidth_l );
              gw_lxgMCT->SetPointError( i, 0, resultsMCT[i].gausswidth_lerr );
              lwMCT->SetPoint( i, timeBinHist[i].bincenter, resultsMCT[i].landwidth );
              lwMCT->SetPointError( i, 0, resultsMCT[i].landwidtherr );

              double ratiowidmpv = resultsMCT[i].landwidth / resultsMCT[i].landmpv;
              double ratiowiderr = resultsMCT[i].landwidtherr / resultsMCT[i].landwidth;
              double ratiompverr = resultsMCT[i].landmpverr / resultsMCT[i].landmpv;
              landwidmpvMCT->SetPoint( i, timeBinHist[i].bincenter, ratiowidmpv);
              landwidmpvMCT->SetPointError( i, 0, ratiowidmpv*sqrt(ratiowiderr*ratiowiderr+ratiompverr*ratiompverr));
              std::cout << "MCT: WIDTHDIVMPV = " << ratiowidmpv << "  WIDTHDIVMPVERR = " << ratiowidmpv*sqrt(ratiowiderr*ratiowiderr+ratiompverr*ratiompverr) << std::endl;

              lxgchi2MCT->SetPoint( i, timeBinHist[i].bincenter, resultsMCT[i].lxg_chi2ndf );

              delete canvlxg2MCT;
              //delete imglxg;
              delete lxgframeMCT;
            }

          resultsMCT[i].binnum = i;
          resultsMCT[i].bincenter = timeBinHist[i].bincenter;
          resultsMCT[i].binwidtherr = ( timeBinHist[i].binwidth )/sqrt( 12.0 );

        }
    }

  if ( DoLxG )
    {
      std::cout << "Landau( x )Gauss: " << std::endl;
      for ( UInt_t i = 0; i < Nbins; i++ )
        {
          std::cout << "     Bin " << i << ": ml=" << results[i].landmpv << " ( +/- " << results[i].landmpverr << " )  sl=" << results[i].landwidth << " ( +/- " << results[i].landwidtherr << " )  sg_l=" << results[i].gausswidth_l << " ( +/- " << results[i].gausswidth_lerr << " )  chi2ndf=" << results[i].lxg_chi2ndf << std::endl;
          std::cout << "MCT: Bin " << i << ": ml=" << resultsMCT[i].landmpv << " ( +/- " << resultsMCT[i].landmpverr << " )  sl=" << resultsMCT[i].landwidth << " ( +/- " << resultsMCT[i].landwidtherr << " )  sg_l=" << resultsMCT[i].gausswidth_l << " ( +/- " << resultsMCT[i].gausswidth_lerr << " )  chi2ndf=" << resultsMCT[i].lxg_chi2ndf << std::endl;
        }
    }
  if ( DoPG )
    {
      std::cout << "PeakFit Gauss: " << std::endl;
      for ( UInt_t i = 0; i < Nbins; i++ )
        {
          std::cout << "     Bin " << i << ": mean=" << results[i].gaussmean << " ( +/- " << results[i].gaussmeanerr << " )  width=" << results[i].gausswidth << " ( +/- " << results[i].gausswidtherr << " )  chi2ndf=" << results[i].pg_chi2ndf << std::endl;
        }
    }

  TF1 * expo = new TF1( "expo", "[0]*exp( -x/[1] )", 100, 1000 );

  expo->SetParNames( "dQdx0", "eLifetime" );
  expo->SetParameters( 3000, 3000 );
  expo->SetLineWidth(3);

  gStyle->SetOptStat( 0 );
  gStyle->SetOptFit( 0 );
  Int_t trans_red = TColor::GetColorTransparent( kRed, 0.35 );

  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");

  if ( DoLxG )
    {
      if (hbin_lxg->GetN() <2)
        {
          /*
             TCanvas * canvlxg = new TCanvas( "canvMPV", "LandauMPV", 1600, 900 );
             canvlxg->cd();
             hbin_lxg->Draw( "ap" );
             hbin_lxg->SetMarkerStyle( kFullDotLarge );
             hbin_lxg->SetMarkerSize( 2 );
             hbin_lxg->GetYaxis()->SetRangeUser( minMPVlxg-100, maxMPVlxg+100 );
             hbin_lxg->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
             hbin_lxg->GetYaxis()->SetTitle( "dQ/dx ( ADC/cm )" );
             canvlxg->Update();
             canvlxg->Write();
             canvlxg->SaveAs( "canvMPV.png" );
           */
        }
      else
        {
          TCanvas * canvlxg = new TCanvas( "canvMPV", "LandauMPV", 2400, 1600 );
          canvlxg->SetMargin(0.12,0.02,0.11,0.02); //l,r,b,t
          canvlxg->cd();
          canvlxg->SetLogy();
          hbin_lxg->Fit( "expo", "R" );
          hbin_lxg->Draw( "ap" );
          hbin_lxg->SetMarkerStyle( kFullTriangleUp );
          hbin_lxg->SetMarkerSize( 2 );
          hbin_lxg->SetLineWidth(3);


          TH1F * hint_lxg = new TH1F( "hint_lxg", "68#% confidence band", 100, 0, drifttimemax );
          ( TVirtualFitter::GetFitter() )->GetConfidenceIntervals( hint_lxg, 0.68 );
          hint_lxg->SetStats( kFALSE );
          hint_lxg->SetFillColor( trans_red );
          //hint_lxg->Draw( "e3 same" );

          hbin_lxg->SetTitle("");
          //hbin_lxg->GetYaxis()->SetMoreLogLabels(true);
          //hbin_lxg->GetYaxis()->SetNoExponent(true);
          hbin_lxg->GetYaxis()->SetNdivisions(0);
          hbin_lxg->GetXaxis()->SetTitle( "Drift Time [ #mus]" );
          hbin_lxg->GetYaxis()->SetTitle( "Most Probable dQ/dx  [ADC/cm]" );
          hbin_lxg->GetYaxis()->SetTitleOffset(1.2);
          hbin_lxg->GetYaxis()->SetRangeUser( minMPVlxg-100, maxMPVlxg+100 );


          TPaveText * pt_lxg = new TPaveText(0.5,0.7,0.97,0.9,"brNDC");
          pt_lxg->AddText(TString::Format("#tau_{lifetime} = %.0f +/- %.0f #mus",expo->GetParameter(1),expo->GetParError(1)));
          pt_lxg->AddText(TString::Format("dQ/dx0 = %.0f +/- %.0f ADC/cm",expo->GetParameter(0),expo->GetParError(0)));
          pt_lxg->SetTextSize(0.05);
          pt_lxg->SetTextFont(42);
          pt_lxg->SetTextColor(kRed);
          pt_lxg->Draw();

          canvlxg->Update();

          TLine tick;
          tick.SetLineWidth(1);
          tick.SetLineColor(1);
          TLatex latex;
          //latex.SetTextSize(0.035);
          latex.SetTextFont(42);
          latex.SetTextAlign(32);
          for (Int_t i=ceil((minMPVlxg-100)/100)*100; i<maxMPVlxg+100; i+=100)
            {
              tick.DrawLine(0,i,(i%100==0) ? 38 : 20,i);
              if (i%100==0 && i/100 % 2 == 0) latex.DrawLatex(-10,i,TString::Format("%i",i));
            }
          canvlxg->Update();


          canvlxg->Write();
          canvlxg->SaveAs( "canvMPV.png" );
          canvlxg->SaveAs("canvMPV.C");

          TGraphErrors * lifetimedrift = new TGraphErrors();
          lifetimedrift->SetTitle( "Fitted Electron Lifetime vs. Fiducial Drift Cut" );
          lifetimedrift->SetMarkerStyle( 20 );
          lifetimedrift->SetMarkerColor( kBlue );
          lifetimedrift->GetXaxis()->SetTitle( "Drift Time [ #mus]" );
          lifetimedrift->GetYaxis()->SetTitle( "eLifetime [ #mus]" );
          TF1 * expotrunc = new TF1( "expotrunc", "[0]*exp( -x/[1] )", 0, 1 );
          expotrunc->SetParNames( "dQdx0", "eLifetime" );
          expotrunc->SetParameters( 3000, 3000 );
          UInt_t nbinwidth = 5;
          int whichbin = 0;
          for ( UInt_t i_bin = 1; i_bin < Nbins-nbinwidth; i_bin += nbinwidth )
            {
              expotrunc->SetRange( timeBinHist[i_bin].binlow, timeBinHist[i_bin+nbinwidth].binhigh );
              if (hbin_lxg->GetN() != 0) hbin_lxg->Fit( "expotrunc", "RQ" );
              lifetimedrift->SetPoint( whichbin, timeBinHist[i_bin+nbinwidth].bincenter, expotrunc->GetParameter( 1 ) );
              lifetimedrift->SetPointError( whichbin, 0, expotrunc->GetParError( 1 ) );
              ++whichbin;
            }
          TCanvas * canvlifetimedrift_lxg = new TCanvas( "cltd_lxg", "cltd_lxg", 1600, 900 );
          canvlifetimedrift_lxg->cd();
          lifetimedrift->Draw( "alpe" );
          lifetimedrift->GetXaxis()->SetTitle( "Drift Distance Cut [cm]" );
          lifetimedrift->GetYaxis()->SetTitle( "Elifetime [ #mus]" );
          lifetimedrift->SetTitle("");
          canvlifetimedrift_lxg->Update();
          canvlifetimedrift_lxg->Write();
          canvlifetimedrift_lxg->SaveAs( "cltd_lxg.png" );

        }

      TCanvas * canvgw_lxg = new TCanvas( "canvgw_lxg", "GaussWidth of Lxg", 1000, 800 );
      canvgw_lxg->cd();
      gw_lxg->Draw( "ap" );
      gw_lxg->SetMarkerStyle( kFullDotLarge );
      gw_lxg->SetMarkerSize( 2 );
      gw_lxg->GetXaxis()->SetTitle( "Drift Time [ #mus]" );
      gw_lxg->GetYaxis()->SetTitle( "Gauss Width [ #mus]" );
      gw_lxg->SetTitle("");
      canvgw_lxg->Update();
      canvgw_lxg->Write();
      canvgw_lxg->SaveAs( "canvgw_lxg.png" );

      TCanvas * canvlw = new TCanvas( "canvlw", "LandauWidth", 1600, 900 );
      canvlw->cd();
      lw->Draw( "ap" );
      lw->SetMarkerStyle( kFullDotLarge );
      lw->SetMarkerSize( 2 );
      lw->GetXaxis()->SetTitle( "Drift Time [ #mus]" );
      lw->GetYaxis()->SetTitle( "Landau Width [ #mus]" );
      lw->SetTitle("");
      canvlw->Update();
      canvlw->Write();
      canvlw->SaveAs( "canvlw.png" );

      TCanvas * canvlxgchi2 = new TCanvas( "canvlxgchi2", "LxG Chi2/NDF", 2000, 1600 );
      canvlxgchi2->SetMargin(0.1,0.02,0.11,0.02); //l,r,b,t
      canvlxgchi2->cd();
      lxgchi2->Draw( "ap" );
      lxgchi2->SetMarkerStyle( kFullDotLarge );
      lxgchi2->SetMarkerSize( 3 );
      lxgchi2->GetXaxis()->SetTitle( "Bin Drift Time [ #mus]" );
      lxgchi2->GetYaxis()->SetTitle( "#chi^{2}/NDF of L#otimesG Fit" );
      lxgchi2->GetXaxis()->SetNdivisions(505,kTRUE);
      lxgchi2->SetTitle("");
      canvlxgchi2->Update();
      canvlxgchi2->Write();
      canvlxgchi2->SaveAs( "canvlxgchi2.png" );

      TCanvas * canvlandwidmpv = new TCanvas("canvlandwidmpv","LandWidth/LandMPV", 1600,900);
      canvlandwidmpv->cd();
      landwidmpv->Draw("ap");
      landwidmpv->SetMarkerStyle(kFullDotLarge);
      landwidmpv->SetMarkerSize(2);
      landwidmpv->GetXaxis()->SetTitle("Drift Time [ #mus]");
      landwidmpv->GetYaxis()->SetTitle("LandWidth / LandMPV");
      landwidmpv->SetTitle("");
      canvlandwidmpv->Update();
      canvlandwidmpv->Write();
      canvlandwidmpv->SaveAs("canvlandwidmpv.png");
    }

  if ( DoLxG && analyse_MCTruth)
    {
      if (hbin_lxgMCT->GetN() <2)
        {
          /*
             TCanvas * canvlxg = new TCanvas( "canvMPV", "LandauMPV", 1600, 900 );
             canvlxg->cd();
             hbin_lxg->Draw( "ap" );
             hbin_lxg->SetMarkerStyle( kFullDotLarge );
             hbin_lxg->SetMarkerSize( 2 );
             hbin_lxg->GetYaxis()->SetRangeUser( minMPVlxg-100, maxMPVlxg+100 );
             hbin_lxg->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
             hbin_lxg->GetYaxis()->SetTitle( "dQ/dx ( ADC/cm )" );
             canvlxg->Update();
             canvlxg->Write();
             canvlxg->SaveAs( "canvMPV.png" );
           */
        }
      else
        {
          TCanvas * canvlxgMCT = new TCanvas( "canvMPVMCT", "LandauMPVMCT", 1600, 900 );
          canvlxgMCT->cd();
          canvlxgMCT->SetLogy();
          hbin_lxgMCT->Fit( "expo", "R" );
          hbin_lxgMCT->Draw( "ap" );
          hbin_lxgMCT->SetMarkerStyle( kFullDotLarge );
          hbin_lxgMCT->SetMarkerSize( 1 );


          TH1F * hint_lxg = new TH1F( "hint_lxgMCT", "68#% confidence band", 100, 0, drifttimemax );
          ( TVirtualFitter::GetFitter() )->GetConfidenceIntervals( hint_lxg, 0.68 );
          hint_lxg->SetStats( kFALSE );
          hint_lxg->SetFillColor( trans_red );
          //hint_lxg->Draw( "e3 same" );

          hbin_lxgMCT->SetTitle("");
          //hbin_lxg->GetYaxis()->SetMoreLogLabels(true);
          //hbin_lxg->GetYaxis()->SetNoExponent(true);
          hbin_lxgMCT->GetYaxis()->SetNdivisions(0);
          hbin_lxgMCT->GetXaxis()->SetTitle( "Drift Time [ #mus]" );
          hbin_lxgMCT->GetYaxis()->SetTitle( "Most Probable dQ/dx  [ADC/cm]" );
          hbin_lxgMCT->GetYaxis()->SetTitleOffset(1.2);
          hbin_lxgMCT->GetYaxis()->SetRangeUser( minMPVlxgMCT-100, maxMPVlxgMCT+100 );


          TPaveText * pt_lxg = new TPaveText(0.6,0.74,0.9,0.9,"brNDC");
          pt_lxg->AddText(TString::Format("#tau_{lifetime} = %.0f +/- %.0f #mus",expo->GetParameter(1),expo->GetParError(1)));
          pt_lxg->AddText(TString::Format("dQdx0 = %.0f +/- %.0f ADC/cm",expo->GetParameter(0),expo->GetParError(0)));
          pt_lxg->SetTextSize(0.038);
          pt_lxg->SetTextFont(42);
          pt_lxg->SetTextColor(kRed);
          pt_lxg->Draw();

          canvlxgMCT->Update();

          TLine tick;
          tick.SetLineWidth(1);
          tick.SetLineColor(1);
          TLatex latex;
          latex.SetTextSize(0.035);
          latex.SetTextFont(42);
          latex.SetTextAlign(32);
          for (Int_t i=ceil((minMPVlxgMCT-100)/100)*100; i<maxMPVlxgMCT+100; i+=100)
            {
              tick.DrawLine(0,i,(i%100==0) ? 38 : 20,i);
              if (i%100==0) latex.DrawLatex(-10,i,TString::Format("%i",i));
            }
          canvlxgMCT->Update();


          canvlxgMCT->Write();
          canvlxgMCT->SaveAs( "canvMPVMCT.png" );
          canvlxgMCT->SaveAs("canvMPVMCT.C");

          TGraphErrors * lifetimedrift = new TGraphErrors();
          lifetimedrift->SetTitle( "Fitted Electron Lifetime vs. Fiducial Drift Cut" );
          lifetimedrift->SetMarkerStyle( 20 );
          lifetimedrift->SetMarkerColor( kBlue );
          lifetimedrift->GetXaxis()->SetTitle( "Drift Time [ #mus]" );
          lifetimedrift->GetYaxis()->SetTitle( "eLifetime [ #mus]" );
          TF1 * expotrunc = new TF1( "expotrunc", "[0]*exp( -x/[1] )", 0, 1 );
          expotrunc->SetParNames( "dQdx0", "eLifetime" );
          expotrunc->SetParameters( 3000, 3000 );
          UInt_t nbinwidth = 5;
          int whichbin = 0;
          for ( UInt_t i_bin = 1; i_bin < Nbins-nbinwidth; i_bin += nbinwidth )
            {
              expotrunc->SetRange( timeBinHist[i_bin].binlow, timeBinHist[i_bin+nbinwidth].binhigh );
              if (hbin_lxgMCT->GetN() != 0) hbin_lxgMCT->Fit( "expotrunc", "RQ" );
              lifetimedrift->SetPoint( whichbin, timeBinHist[i_bin+nbinwidth].bincenter, expotrunc->GetParameter( 1 ) );
              lifetimedrift->SetPointError( whichbin, 0, expotrunc->GetParError( 1 ) );
              ++whichbin;
            }
          TCanvas * canvlifetimedrift_lxgMCT = new TCanvas( "cltd_lxgMCT", "cltd_lxgMCT", 1600, 900 );
          canvlifetimedrift_lxgMCT->cd();
          lifetimedrift->Draw( "alpe" );
          lifetimedrift->GetXaxis()->SetTitle( "Drift Distance Cut [cm]" );
          lifetimedrift->GetYaxis()->SetTitle( "Elifetime [ #mus]" );
          lifetimedrift->SetTitle("");
          canvlifetimedrift_lxgMCT->Update();
          canvlifetimedrift_lxgMCT->Write();
          canvlifetimedrift_lxgMCT->SaveAs( "cltd_lxgMCT.png" );

        }

      TCanvas * canvgw_lxgMCT = new TCanvas( "canvgw_lxgMCT", "GaussWidth of Lxg MCT", 1000, 800 );
      canvgw_lxgMCT->cd();
      gw_lxgMCT->Draw( "ap" );
      gw_lxgMCT->SetMarkerStyle( kFullDotLarge );
      gw_lxgMCT->SetMarkerSize( 2 );
      gw_lxgMCT->GetXaxis()->SetTitle( "Drift Time [ #mus]" );
      gw_lxgMCT->GetYaxis()->SetTitle( "Gauss Width [ #mus]" );
      gw_lxgMCT->SetTitle("");
      canvgw_lxgMCT->Update();
      canvgw_lxgMCT->Write();
      canvgw_lxgMCT->SaveAs( "canvgw_lxgMCT.png" );

      TCanvas * canvlwMCT = new TCanvas( "canvlwMCT", "LandauWidthMCT", 1600, 900 );
      canvlwMCT->cd();
      lwMCT->Draw( "ap" );
      lwMCT->SetMarkerStyle( kFullDotLarge );
      lwMCT->SetMarkerSize( 2 );
      lwMCT->GetXaxis()->SetTitle( "Drift Time [ #mus]" );
      lwMCT->GetYaxis()->SetTitle( "Landau Width [ #mus]" );
      lwMCT->SetTitle("");
      canvlwMCT->Update();
      canvlwMCT->Write();
      canvlwMCT->SaveAs( "canvlwMCT.png" );

      TCanvas * canvlxgchi2MCT = new TCanvas( "canvlxgchi2MCT", "LxG Chi2/NDF MCT", 1000, 800 );
      canvlxgchi2MCT->cd();
      lxgchi2MCT->Draw( "ap" );
      lxgchi2MCT->SetMarkerStyle( kFullDotLarge );
      lxgchi2MCT->SetMarkerSize( 2 );
      lxgchi2MCT->GetXaxis()->SetTitle( "Drift Time [ #mus]" );
      lxgchi2MCT->GetYaxis()->SetTitle( "Chi^2/NDF of L(x)g Fit" );
      lxgchi2MCT->SetTitle("");
      canvlxgchi2MCT->Update();
      canvlxgchi2MCT->Write();
      canvlxgchi2MCT->SaveAs( "canvlxgchi2MCT.png" );

      TCanvas * canvlandwidmpvMCT = new TCanvas("canvlandwidmpvMCT","LandWidth/LandMPV MCT", 1600,900);
      canvlandwidmpvMCT->cd();
      landwidmpvMCT->Draw("ap");
      landwidmpvMCT->SetMarkerStyle(kFullDotLarge);
      landwidmpvMCT->SetMarkerSize(2);
      landwidmpvMCT->GetXaxis()->SetTitle("Drift Time [ #mus]");
      landwidmpvMCT->GetYaxis()->SetTitle("LandWidth / LandMPV");
      landwidmpvMCT->SetTitle("");
      canvlandwidmpvMCT->Update();
      canvlandwidmpvMCT->Write();
      canvlandwidmpvMCT->SaveAs("canvlandwidmpvMCT.png");
    }

  if ( DoPG )
    {
      TCanvas * canvpg = new TCanvas( "canvpg", "PeakGaussMean", 1600, 900 );
      canvpg->cd();
      hbin_pg->Fit( "expo", "R" );
      hbin_pg->Draw( "ap" );
      hbin_pg->SetMarkerStyle( kFullDotLarge );
      hbin_pg->SetMarkerSize( 2 );

      TH1F * hint_pg = new TH1F( "hint_pg", "95#% confidence band", 100, 0, drifttimemax );
      ( TVirtualFitter::GetFitter() )->GetConfidenceIntervals( hint_pg, 0.68 );
      hint_pg->SetStats( kFALSE );
      hint_pg->SetFillColor( trans_red );
      hint_pg->Draw( "e3 same" );

      hbin_pg->GetYaxis()->SetRangeUser( minMPVpg-100, maxMPVpg+100 );
      hbin_pg->GetXaxis()->SetTitle( "Drift Time [ #mus]" );
      hbin_pg->GetYaxis()->SetTitle( "dQ/dx [ADC/cm]" );
      canvpg->Update();
      canvpg->Write();
      //canvpg->SaveAs( "canvpg.png" );

      TCanvas * canvgw_pg = new TCanvas( "canvgw_pg", "PeakGaussWidth", 1600, 900 );
      canvgw_pg->cd();
      gw_pg->Draw( "ap" );
      gw_pg->SetMarkerStyle( kFullDotLarge );
      gw_pg->SetMarkerSize( 2 );
      gw_pg->GetXaxis()->SetTitle( "Drift Time ( #mu s )" );
      gw_pg->GetYaxis()->SetTitle( "Gauss Width ( ticks )" );
      canvgw_pg->Update();
      canvgw_pg->Write();
      //canvgw_pg->SaveAs( "canvgw_pg.png" );

      TGraphErrors * lifetimedrift = new TGraphErrors();
      lifetimedrift->SetTitle( "Fitted Electron Lifetime vs. Fiducial Drift Cut" );
      lifetimedrift->SetMarkerStyle( 20 );
      lifetimedrift->SetMarkerColor( kBlue );
      lifetimedrift->GetXaxis()->SetTitle( "Drift Time [ #mus]" );
      lifetimedrift->GetYaxis()->SetTitle( "eLifetime [ #mus]" );
      TF1 * expotrunc = new TF1( "expotrunc", "[0]*exp( -x/[1] )", 0, 1 );
      expotrunc->SetParNames( "dQdx0", "eLifetime" );
      expotrunc->SetParameters( 3000, 3000 );
      UInt_t nbinwidth=4;
      for ( UInt_t i_bin = 1; i_bin < Nbins-nbinwidth; ++i_bin )
        {
          expotrunc->SetRange( timeBinHist[i_bin].binlow, timeBinHist[i_bin+nbinwidth].binhigh );
          hbin_pg->Fit( "expotrunc", "RQ" );
          lifetimedrift->SetPoint( i_bin-1, timeBinHist[i_bin+nbinwidth].bincenter, expotrunc->GetParameter( 1 ) );
          lifetimedrift->SetPointError( i_bin-1, 0, expotrunc->GetParError( 1 ) );
        }
      TCanvas * canvlifetimedrift_pg = new TCanvas( "cltd_pg", "cltd_pg", 1600, 900 );
      canvlifetimedrift_pg->cd();
      lifetimedrift->Draw( "alpe" );
      lifetimedrift->GetXaxis()->SetTitle( "Drift Distance Cut [cm]" );
      lifetimedrift->GetYaxis()->SetTitle( "Elifetime [ #mus]" );
      canvlifetimedrift_pg->Update();
      canvlifetimedrift_pg->Write();
      //canvlifetimedrift_pg->SaveAs( "cltd_pg.png" );
    }

  TTree * resulttree = new TTree( "results", "Results of fits" );
  Float_t landmpv;
  Float_t landwidth;
  Float_t gaussmean_l;
  Float_t gausswidth_l;
  Float_t landmpverr;
  Float_t landwidtherr;
  Float_t gaussmean_lerr;
  Float_t gausswidth_lerr;
  Float_t lxg_chi2ndf;
  Float_t lxglow;
  Float_t lxghigh;
  Bool_t fitsuccess;
  resulttree->Branch( "landmpv", &landmpv, "landmpv/F" );
  resulttree->Branch( "landwidth", &landwidth, "landwidth/F" );
  resulttree->Branch( "gaussmean_l", &gaussmean_l, "gaussmean_l/F" );
  resulttree->Branch( "gausswidth_l", &gausswidth_l, "gausswidth_l/F" );
  resulttree->Branch( "landmpverr", &landmpverr, "landmpverr/F" );
  resulttree->Branch( "landwidtherr", &landwidtherr, "landwidtherr/F" );
  resulttree->Branch( "gaussmean_lerr", &gaussmean_lerr, "gaussmean_lerr/F" );
  resulttree->Branch( "gausswidth_lerr", &gausswidth_lerr, "gausswidth_lerr/F" );
  resulttree->Branch( "lxg_chi2ndf", &lxg_chi2ndf, "lxg_chi2ndf/F" );
  resulttree->Branch( "lxglow", &lxglow, "lxglow/F" );
  resulttree->Branch( "lxghigh", &lxghigh, "lxghigh/F" );
  resulttree->Branch( "fitsuccess", &fitsuccess, "fitsuccess/O");

  Float_t gaussmean;
  Float_t gausswidth;
  Float_t gaussmeanerr;
  Float_t gausswidtherr;
  Float_t bincenter;
  Float_t binwidtherr;
  Float_t pg_chi2ndf;
  resulttree->Branch( "binnum", &binnum, "binnum/I" );
  resulttree->Branch( "gaussmean", &gaussmean, "gaussmean/F" );
  resulttree->Branch( "gausswidth", &gausswidth, "gausswidth/F" );
  resulttree->Branch( "gaussmeanerr", &gaussmeanerr, "gaussmeanerr/F" );
  resulttree->Branch( "gausswidtherr", &gausswidtherr, "gausswidtherr/F" );
  resulttree->Branch( "bincenter", &bincenter, "bincenter/F" );
  resulttree->Branch( "binwidtherr", &binwidtherr, "binwidtherr/F" );
  resulttree->Branch( "pg_chi2ndf", &pg_chi2ndf, "pg_chi2ndf/F" );

  for ( auto const & r : results )
    {
      binnum = r.binnum;
      landmpv = r.landmpv;
      landwidth = r.landwidth;
      gaussmean = r.gaussmean;
      gausswidth = r.gausswidth;
      gaussmean_l = r.gaussmean_l;
      gausswidth_l = r.gausswidth_l;
      landmpverr = r.landmpverr;
      landwidtherr = r.landwidtherr;
      gaussmeanerr = r.gaussmeanerr;
      gausswidtherr = r.gausswidtherr;
      gaussmean_lerr = r.gaussmean_lerr;
      gausswidth_lerr = r.gausswidth_lerr;
      bincenter = r.bincenter;
      binwidtherr = r.binwidtherr;
      lxg_chi2ndf = r.lxg_chi2ndf;
      pg_chi2ndf = r.pg_chi2ndf;
      lxglow = r.lxglow;
      lxghigh = r.lxghigh;
      fitsuccess = r.fitsuccess;
      resulttree->Fill();
    }
  resulttree->Write();

  TTree * resulttreeMCT = new TTree( "resultsMCT", "Results of fits MCT" );
  resulttreeMCT->Branch( "landmpv", &landmpv, "landmpv/F" );
  resulttreeMCT->Branch( "landwidth", &landwidth, "landwidth/F" );
  resulttreeMCT->Branch( "gaussmean_l", &gaussmean_l, "gaussmean_l/F" );
  resulttreeMCT->Branch( "gausswidth_l", &gausswidth_l, "gausswidth_l/F" );
  resulttreeMCT->Branch( "landmpverr", &landmpverr, "landmpverr/F" );
  resulttreeMCT->Branch( "landwidtherr", &landwidtherr, "landwidtherr/F" );
  resulttreeMCT->Branch( "gaussmean_lerr", &gaussmean_lerr, "gaussmean_lerr/F" );
  resulttreeMCT->Branch( "gausswidth_lerr", &gausswidth_lerr, "gausswidth_lerr/F" );
  resulttreeMCT->Branch( "lxg_chi2ndf", &lxg_chi2ndf, "lxg_chi2ndf/F" );
  resulttreeMCT->Branch( "lxglow", &lxglow, "lxglow/F" );
  resulttreeMCT->Branch( "lxghigh", &lxghigh, "lxghigh/F" );
  resulttreeMCT->Branch( "fitsuccess", &fitsuccess, "fitsuccess/O");
  resulttreeMCT->Branch( "binnum", &binnum, "binnum/I" );
  resulttreeMCT->Branch( "bincenter", &bincenter, "bincenter/F" );
  resulttreeMCT->Branch( "binwidtherr", &binwidtherr, "binwidtherr/F" );

  for ( auto const & r : resultsMCT )
    {
      binnum = r.binnum;
      landmpv = r.landmpv;
      landwidth = r.landwidth;
      gaussmean = r.gaussmean;
      gausswidth = r.gausswidth;
      gaussmean_l = r.gaussmean_l;
      gausswidth_l = r.gausswidth_l;
      landmpverr = r.landmpverr;
      landwidtherr = r.landwidtherr;
      gaussmeanerr = r.gaussmeanerr;
      gausswidtherr = r.gausswidtherr;
      gaussmean_lerr = r.gaussmean_lerr;
      gausswidth_lerr = r.gausswidth_lerr;
      bincenter = r.bincenter;
      binwidtherr = r.binwidtherr;
      lxg_chi2ndf = r.lxg_chi2ndf;
      pg_chi2ndf = r.pg_chi2ndf;
      lxglow = r.lxglow;
      lxghigh = r.lxghigh;
      fitsuccess = r.fitsuccess;
      resulttreeMCT->Fill();
    }
  resulttreeMCT->Write();

  std::cout << "# Runs included = " << runseventshits.size() << std::endl;
  TH1F * hitsperevent = new TH1F("hpe","# Hits Per Event",100,0,340);
  int numberevents = 0;
  int numberhits = 0;
  for (auto run : runseventshits)
    {
      for (auto event : run.second)
        {
          numberevents++;
          numberhits += event.second;
          hitsperevent->Fill(event.second);
        }
    }
  std::cout << "# Events included = " << numberevents << std::endl;
  std::cout << "# Hits included = " << numberhits << std::endl;
  TCanvas * canvhitsperevent = new TCanvas("canvhitsperevent","",2000,1600);
  hitsperevent->Draw();
  canvhitsperevent->Update();
  canvhitsperevent->SaveAs("canvhitsperevent.png");

  std::cout << "Found Number = " << (int)(found->GetEntries()) << "  Assumed Number = " << (int)(assumed->GetEntries()) << std::endl;
  fileout->Close();
  return;
}
