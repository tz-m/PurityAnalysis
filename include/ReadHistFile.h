#ifndef READHISTFILE_H
#define READHISTFILE_H

#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <string>
#include <sstream>
#include <map>
#include <algorithm>

#include "UsefulTypes.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<double>+;
#pragma link C++ class vector<int>+;
#pragma link C++ class vector<float>+;
#endif

class ReadHistFile {

public:

  ReadHistFile();
  ReadHistFile(std::string filename);
  UInt_t ReadFile();

  types::HitMap * GetHitMap() {
    return &hits;
  }

  std::vector<Int_t> GetRunList() {
    return runs;
  }

private:
  types::HitMap hits;
  std::vector<Int_t> runs;
  UInt_t i;


};

ReadHistFile::ReadHistFile()
{
}

ReadHistFile::ReadHistFile(std::string filename)
{
  TFile * file = TFile::Open(filename.c_str(),"READ");
  if (!file || file->IsZombie())
    {
      std::stringstream ss;
      ss << "ReadHistFile::ReadFile() -- Input file is not read";
      throw std::runtime_error(ss.str());
    }

  TTree * tree = (TTree*)file->Get("robusthit/RobustHitFinder");
  //!!!!!!! For some reason, I used robusthit/RobustHitFinder for some, and robusthits/RobustHitFinder for others. Check both!

  TTreeReader reader;

  TTreeReaderValue<Int_t> *run;
  TTreeReaderValue<Int_t> *event;
  //TTreeReaderValue<Double_t> *t0;
  TTreeReaderValue<UInt_t> *c1;
  TTreeReaderValue<UInt_t> *c2;
  TTreeReaderValue<UInt_t> *trignum;
  TTreeReaderValue<Float_t> *c1x;
  TTreeReaderValue<Float_t> *c1y;
  TTreeReaderValue<Float_t> *c1z;
  TTreeReaderValue<Float_t> *c2x;
  TTreeReaderValue<Float_t> *c2y;
  TTreeReaderValue<Float_t> *c2z;
  //TTreeReaderValue<Float_t> *distancecut;
  TTreeReaderValue<Int_t> *channel;
  //TTreeReaderValue<Int_t> *wire;
  TTreeReaderValue<Int_t> *tpc;
  TTreeReaderValue<Int_t> *signalsize;
  //TTreeReaderArray<Float_t> * signal;
  //TTreeReaderArray<Float_t> * signalFilter;
  TTreeReaderValue<Float_t> *baseline;
  TTreeReaderValue<Float_t> *rms;
  TTreeReaderValue<Float_t> *baselineFilter;
  TTreeReaderValue<Float_t> *rmsFilter;
  //TTreeReaderValue<Float_t> *pedmean;
  //TTreeReaderValue<Float_t> *pedrms;
  TTreeReaderValue<Float_t> *integral;
  TTreeReaderValue<Float_t> *integralFilter;
  TTreeReaderValue<Float_t> *sumADC;
  TTreeReaderValue<Float_t> *sigmaintegral;
  TTreeReaderValue<Float_t> *sigmaintegralFilter;
  TTreeReaderValue<Float_t> *amplitude;
  TTreeReaderValue<Float_t> *amplitudeFilter;
  TTreeReaderValue<Float_t> *peaktick;
  TTreeReaderValue<Float_t> *peaktickFilter;
  TTreeReaderValue<Float_t> *peaktime;
  TTreeReaderValue<Float_t> *peaktimeFilter;
  TTreeReaderValue<Int_t> *begintick;
  TTreeReaderValue<Int_t> *endtick;
  TTreeReaderValue<Int_t> *width;
  TTreeReaderValue<Float_t> *hitx;
  TTreeReaderValue<Float_t> *hity;
  TTreeReaderValue<Float_t> *hitz;
  TTreeReaderValue<Float_t> *hiterrxlo;
  TTreeReaderValue<Float_t> *hiterrxhi;
  TTreeReaderValue<Float_t> *hiterrylo;
  TTreeReaderValue<Float_t> *hiterryhi;
  TTreeReaderValue<Float_t> *hiterrzlo;
  TTreeReaderValue<Float_t> *hiterrzhi;
  TTreeReaderValue<Float_t> *perpdist;
  TTreeReaderValue<Float_t> *hitt;
  TTreeReaderValue<Float_t> *driftdist;
  TTreeReaderValue<Bool_t> *countercut;
  TTreeReaderValue<Float_t> *fitconstant;
  TTreeReaderValue<Float_t> *fitconstanterr;
  TTreeReaderValue<Float_t> *fitlinear;
  TTreeReaderValue<Float_t> *fitlinearerr;
  TTreeReaderValue<Float_t> *fitquadratic;
  TTreeReaderValue<Float_t> *fitquadraticerr;
  TTreeReaderValue<Float_t> *fitchi2;
  TTreeReaderValue<Float_t> *fitsumsqrresidual;
  TTreeReaderValue<Float_t> *fitndf;
  TTreeReaderValue<Float_t> *fitmle;
  TTreeReaderValue<Bool_t> *fitsuccess;
  TTreeReaderValue<Bool_t> *fitrealhit;
  TTreeReaderValue<Float_t> *segmentlength;
  TTreeReaderValue<Bool_t> *assumedhit;
/*  TTreeReaderValue<Int_t> *numGoodHitsChan;
   TTreeReaderValue<Int_t> *nwiresTPC0;
   TTreeReaderValue<Int_t> *nwiresTPC1;
   TTreeReaderValue<Int_t> *nwiresTPC2;
   TTreeReaderValue<Int_t> *nwiresTPC3;
   TTreeReaderValue<Int_t> *nwiresTPC4;
   TTreeReaderValue<Int_t> *nwiresTPC5;
   TTreeReaderValue<Int_t> *nwiresTPC6;
   TTreeReaderValue<Int_t> *nwiresTPC7;
   TTreeReaderValue<Float_t> *prebaseline;
   TTreeReaderValue<Float_t> *postbaseline;
   TTreeReaderValue<Float_t> *prebaserms;
   TTreeReaderValue<Float_t> *postbaserms;
   TTreeReaderValue<Int_t> *trackid;
   TTreeReaderValue<Int_t> *numtrajpts;
   TTreeReaderValue<Double_t> *tracklength;
   TTreeReaderValue<Bool_t> *isOnTrack;
   TTreeReaderValue<Double_t> *dqdxatpt;
   TTreeReaderValue<Int_t> *prevStartTick;
   TTreeReaderValue<Int_t> *prevEndTick;
   TTreeReaderValue<Float_t> *prevPeakTime;
   TTreeReaderValue<Float_t> *prevSigmaPeakTime;
   TTreeReaderValue<Float_t> *prevRMS;
   TTreeReaderValue<Float_t> *prevPeakAmplitude;
   TTreeReaderValue<Float_t> *prevSigmaPeakAmplitude;
   TTreeReaderValue<Float_t> *prevSummedADC;
   TTreeReaderValue<Float_t> *prevIntegral;
   TTreeReaderValue<Float_t> *prevSigmaIntegral;
 */


  reader.SetTree(tree);

  if (tree->GetBranchStatus("run")) run = new TTreeReaderValue<Int_t>(reader,"run");
  if (tree->GetBranchStatus("event")) event = new TTreeReaderValue<Int_t>(reader,"event");
  //if (tree->GetBranchStatus("t0")) t0 = new TTreeReaderValue<Double_t>(reader,"t0");
  if (tree->GetBranchStatus("c1")) c1 = new TTreeReaderValue<UInt_t>(reader,"c1");
  if (tree->GetBranchStatus("c2")) c2 = new TTreeReaderValue<UInt_t>(reader,"c2");
  if (tree->GetBranchStatus("trignum")) trignum = new TTreeReaderValue<UInt_t>(reader,"trignum");
  if (tree->GetBranchStatus("c1x")) c1x = new TTreeReaderValue<Float_t>(reader,"c1x");
  if (tree->GetBranchStatus("c1y")) c1y = new TTreeReaderValue<Float_t>(reader,"c1y");
  if (tree->GetBranchStatus("c1z")) c1z = new TTreeReaderValue<Float_t>(reader,"c1z");
  if (tree->GetBranchStatus("c2x")) c2x = new TTreeReaderValue<Float_t>(reader,"c2x");
  if (tree->GetBranchStatus("c2y")) c2y = new TTreeReaderValue<Float_t>(reader,"c2y");
  if (tree->GetBranchStatus("c2z")) c2z = new TTreeReaderValue<Float_t>(reader,"c2z");
  //if (tree->GetBranchStatus("distancecut")) distancecut = new TTreeReaderValue<Float_t>(reader,"distancecut");
  if (tree->GetBranchStatus("channel")) channel = new TTreeReaderValue<Int_t>(reader,"channel");
  //if (tree->GetBranchStatus("wire")) wire = new TTreeReaderValue<Int_t>(reader,"wire");
  if (tree->GetBranchStatus("tpc")) tpc = new TTreeReaderValue<Int_t>(reader,"tpc");
  if (tree->GetBranchStatus("signalsize")) signalsize = new TTreeReaderValue<Int_t>(reader,"signalsize");
  //if (tree->GetBranchStatus("signal")) signal = new TTreeReaderArray<Float_t>(reader,"signal");
  //if (tree->GetBranchStatus("signalFilter")) signalFilter = new TTreeReaderArray<Float_t>(reader,"signalFilter");
  if (tree->GetBranchStatus("baseline")) baseline = new TTreeReaderValue<Float_t>(reader,"baseline");
  if (tree->GetBranchStatus("rms")) rms = new TTreeReaderValue<Float_t>(reader,"rms");
  if (tree->GetBranchStatus("baselineFilter")) baselineFilter = new TTreeReaderValue<Float_t>(reader,"baselineFilter");
  if (tree->GetBranchStatus("rmsFilter")) rmsFilter = new TTreeReaderValue<Float_t>(reader,"rmsFilter");
  //if (tree->GetBranchStatus("pedmean")) pedmean = new TTreeReaderValue<Float_t>(reader,"pedmean");
  //if (tree->GetBranchStatus("pedrms")) pedrms = new TTreeReaderValue<Float_t>(reader,"pedrms");
  if (tree->GetBranchStatus("integral")) integral = new TTreeReaderValue<Float_t>(reader,"integral");
  if (tree->GetBranchStatus("integralFilter")) integralFilter = new TTreeReaderValue<Float_t>(reader,"integralFilter");
  if (tree->GetBranchStatus("sumADC")) sumADC = new TTreeReaderValue<Float_t>(reader,"sumADC");
  if (tree->GetBranchStatus("sigmaintegral")) sigmaintegral = new TTreeReaderValue<Float_t>(reader,"sigmaintegral");
  if (tree->GetBranchStatus("sigmaintegralFilter")) sigmaintegralFilter = new TTreeReaderValue<Float_t>(reader,"sigmaintegralFilter");
  if (tree->GetBranchStatus("amplitude")) amplitude = new TTreeReaderValue<Float_t>(reader,"amplitude");
  if (tree->GetBranchStatus("amplitudeFilter")) amplitudeFilter = new TTreeReaderValue<Float_t>(reader,"amplitudeFilter");
  if (tree->GetBranchStatus("peaktick")) peaktick = new TTreeReaderValue<Float_t>(reader,"peaktick");
  if (tree->GetBranchStatus("peaktickFilter")) peaktickFilter = new TTreeReaderValue<Float_t>(reader,"peaktickFilter");
  if (tree->GetBranchStatus("peaktime")) peaktime = new TTreeReaderValue<Float_t>(reader,"peaktime");
  if (tree->GetBranchStatus("peaktimeFilter")) peaktimeFilter = new TTreeReaderValue<Float_t>(reader,"peaktimeFilter");
  if (tree->GetBranchStatus("begintick")) begintick = new TTreeReaderValue<Int_t>(reader,"begintick");
  if (tree->GetBranchStatus("endtick")) endtick = new TTreeReaderValue<Int_t>(reader,"endtick");
  if (tree->GetBranchStatus("width")) width = new TTreeReaderValue<Int_t>(reader,"width");
  if (tree->GetBranchStatus("hitx")) hitx = new TTreeReaderValue<Float_t>(reader,"hitx");
  if (tree->GetBranchStatus("hity")) hity = new TTreeReaderValue<Float_t>(reader,"hity");
  if (tree->GetBranchStatus("hitz")) hitz = new TTreeReaderValue<Float_t>(reader,"hitz");
  if (tree->GetBranchStatus("hiterrxlo")) hiterrxlo = new TTreeReaderValue<Float_t>(reader,"hiterrxlo");
  if (tree->GetBranchStatus("hiterrxhi")) hiterrxhi = new TTreeReaderValue<Float_t>(reader,"hiterrxhi");
  if (tree->GetBranchStatus("hiterrylo")) hiterrylo = new TTreeReaderValue<Float_t>(reader,"hiterrylo");
  if (tree->GetBranchStatus("hiterryhi")) hiterryhi = new TTreeReaderValue<Float_t>(reader,"hiterryhi");
  if (tree->GetBranchStatus("hiterrzlo")) hiterrzlo = new TTreeReaderValue<Float_t>(reader,"hiterrzlo");
  if (tree->GetBranchStatus("hiterrzhi")) hiterrzhi = new TTreeReaderValue<Float_t>(reader,"hiterrzhi");
  if (tree->GetBranchStatus("perpdist")) perpdist = new TTreeReaderValue<Float_t>(reader,"perpdist");
  if (tree->GetBranchStatus("hitt")) hitt = new TTreeReaderValue<Float_t>(reader,"hitt");
  if (tree->GetBranchStatus("driftdist")) driftdist = new TTreeReaderValue<Float_t>(reader,"driftdist");
  if (tree->GetBranchStatus("countercut")) countercut = new TTreeReaderValue<Bool_t>(reader,"countercut");
  if (tree->GetBranchStatus("fitconstant")) fitconstant = new TTreeReaderValue<Float_t>(reader,"fitconstant");
  if (tree->GetBranchStatus("fitconstanterr")) fitconstanterr = new TTreeReaderValue<Float_t>(reader,"fitconstanterr");
  if (tree->GetBranchStatus("fitlinear")) fitlinear = new TTreeReaderValue<Float_t>(reader,"fitlinear");
  if (tree->GetBranchStatus("fitlinearerr")) fitlinearerr = new TTreeReaderValue<Float_t>(reader,"fitlinearerr");
  if (tree->GetBranchStatus("fitquadratic")) fitquadratic = new TTreeReaderValue<Float_t>(reader,"fitquadratic");
  if (tree->GetBranchStatus("fitquadraticerr")) fitquadraticerr = new TTreeReaderValue<Float_t>(reader,"fitquadraticerr");
  if (tree->GetBranchStatus("fitchi2")) fitchi2 = new TTreeReaderValue<Float_t>(reader,"fitchi2");
  if (tree->GetBranchStatus("fitsumsqrresidual")) fitsumsqrresidual = new TTreeReaderValue<Float_t>(reader,"fitsumsqrresidual");
  if (tree->GetBranchStatus("fitndf")) fitndf = new TTreeReaderValue<Float_t>(reader,"fitndf");
  if (tree->GetBranchStatus("fitmle")) fitmle = new TTreeReaderValue<Float_t>(reader,"fitmle");
  if (tree->GetBranchStatus("fitsuccess")) fitsuccess = new TTreeReaderValue<Bool_t>(reader,"fitsuccess");
  if (tree->GetBranchStatus("fitrealhit")) fitrealhit = new TTreeReaderValue<Bool_t>(reader,"fitrealhit");
  if (tree->GetBranchStatus("segmentlength")) segmentlength = new TTreeReaderValue<Float_t>(reader,"segmentlength");
  if (tree->GetBranchStatus("assumedhit")) assumedhit = new TTreeReaderValue<Bool_t>(reader,"assumedhit");
  /* if (tree->GetBranchStatus("numGoodHitsChan")) numGoodHitsChan = new TTreeReaderValue<Int_t>(reader,"numGoodHitsChan");
     if (tree->GetBranchStatus("nwiresTPC0")) nwiresTPC0 = new TTreeReaderValue<Int_t>(reader,"nwiresTPC0");
     if (tree->GetBranchStatus("nwiresTPC1")) nwiresTPC1 = new TTreeReaderValue<Int_t>(reader,"nwiresTPC1");
     if (tree->GetBranchStatus("nwiresTPC2")) nwiresTPC2 = new TTreeReaderValue<Int_t>(reader,"nwiresTPC2");
     if (tree->GetBranchStatus("nwiresTPC3")) nwiresTPC3 = new TTreeReaderValue<Int_t>(reader,"nwiresTPC3");
     if (tree->GetBranchStatus("nwiresTPC4")) nwiresTPC4 = new TTreeReaderValue<Int_t>(reader,"nwiresTPC4");
     if (tree->GetBranchStatus("nwiresTPC5")) nwiresTPC5 = new TTreeReaderValue<Int_t>(reader,"nwiresTPC5");
     if (tree->GetBranchStatus("nwiresTPC6")) nwiresTPC6 = new TTreeReaderValue<Int_t>(reader,"nwiresTPC6");
     if (tree->GetBranchStatus("nwiresTPC7")) nwiresTPC7 = new TTreeReaderValue<Int_t>(reader,"nwiresTPC7");
     if (tree->GetBranchStatus("prebaseline")) prebaseline = new TTreeReaderValue<Float_t>(reader,"prebaseline");
     if (tree->GetBranchStatus("postbaseline")) postbaseline = new TTreeReaderValue<Float_t>(reader,"postbaseline");
     if (tree->GetBranchStatus("prebaserms")) prebaserms = new TTreeReaderValue<Float_t>(reader,"prebaserms");
     if (tree->GetBranchStatus("postbaserms")) postbaserms = new TTreeReaderValue<Float_t>(reader,"postbaserms");
     if (tree->GetBranchStatus("trackid")) trackid = new TTreeReaderValue<Int_t>(reader,"trackid");
     if (tree->GetBranchStatus("numtrajpts")) numtrajpts = new TTreeReaderValue<Int_t>(reader,"numtrajpts");
     if (tree->GetBranchStatus("tracklength")) tracklength = new TTreeReaderValue<Double_t>(reader,"tracklength");
     if (tree->GetBranchStatus("isOnTrack")) isOnTrack = new TTreeReaderValue<Bool_t>(reader,"isOnTrack");
     if (tree->GetBranchStatus("dqdxatpt")) dqdxatpt = new TTreeReaderValue<Double_t>(reader,"dqdxatpt");
     if (tree->GetBranchStatus("prevStartTick")) prevStartTick = new TTreeReaderValue<Int_t>(reader,"prevStartTick");
     if (tree->GetBranchStatus("prevEndTick")) prevEndTick = new TTreeReaderValue<Int_t>(reader,"prevEndTick");
     if (tree->GetBranchStatus("prevPeakTime")) prevPeakTime = new TTreeReaderValue<Float_t>(reader,"prevPeakTime");
     if (tree->GetBranchStatus("prevSigmaPeakTime")) prevSigmaPeakTime = new TTreeReaderValue<Float_t>(reader,"prevSigmaPeakTime");
     if (tree->GetBranchStatus("prevRMS")) prevRMS = new TTreeReaderValue<Float_t>(reader,"prevRMS");
     if (tree->GetBranchStatus("prevPeakAmplitude")) prevPeakAmplitude = new TTreeReaderValue<Float_t>(reader,"prevPeakAmplitude");
     if (tree->GetBranchStatus("prevSigmaPeakAmplitude")) prevSigmaPeakAmplitude = new TTreeReaderValue<Float_t>(reader,"prevSigmaPeakAmplitude");
     if (tree->GetBranchStatus("prevSummedADC")) prevSummedADC = new TTreeReaderValue<Float_t>(reader,"prevSummedADC");
     if (tree->GetBranchStatus("prevIntegral")) prevIntegral = new TTreeReaderValue<Float_t>(reader,"prevIntegral");
     if (tree->GetBranchStatus("prevSigmaIntegral")) prevSigmaIntegral = new TTreeReaderValue<Float_t>(reader,"prevSigmaIntegral");*/


  std::cout << "ReadHistFile::ReadFile() -- Reading file " << std::endl;

  i = 0;
  while (reader.Next())
    {
      types::HitInfo hi;
      if (run!=nullptr) hi.run = **run;
      if (event!=nullptr) hi.event = **event;
      //if (t0!=nullptr) hi.t0 = **t0;
      if (c1!=nullptr) hi.c1 = **c1;
      if (c2!=nullptr) hi.c2 = **c2;
      if (trignum!=nullptr) hi.trignum = **trignum;
      if (c1x!=nullptr) hi.c1x = **c1x;
      if (c1y!=nullptr) hi.c1y = **c1y;
      if (c1z!=nullptr) hi.c1z = **c1z;
      if (c2x!=nullptr) hi.c2x = **c2x;
      if (c2y!=nullptr) hi.c2y = **c2y;
      if (c2z!=nullptr) hi.c2z = **c2z;
      //if (distancecut!=nullptr) hi.distancecut = **distancecut;
      if (channel!=nullptr) hi.channel = **channel;
      //if (wire!=nullptr) hi.wire = **wire;
      if (tpc!=nullptr) hi.tpc = **tpc;
      if (signalsize!=nullptr) hi.signalsize = **signalsize;
      //if (signal!=nullptr) for (int j = 0; j < *signalsize; ++j) hi.signal.push_back(signal->At(j));
      //if (signalFilter!=nullptr) for (int j = 0; j < *signalsize; ++j) hi.signalFilter.push_back(signalFilter->At(j));
      if (baseline!=nullptr) hi.baseline = **baseline;
      if (rms!=nullptr) hi.rms = **rms;
      if (baselineFilter!=nullptr) hi.baselineFilter = **baselineFilter;
      if (rmsFilter!=nullptr) hi.rmsFilter = **rmsFilter;
      //if (pedmean!=nullptr) hi.pedmean = **pedmean;
      //if (pedrms!=nullptr) hi.pedrms = **pedrms;
      if (integral!=nullptr) hi.integral = **integral;
      if (integralFilter!=nullptr) hi.integralFilter = **integralFilter;
      if (sumADC!=nullptr) hi.sumADC = **sumADC;
      if (sigmaintegral!=nullptr) hi.sigmaintegral = **sigmaintegral;
      if (sigmaintegralFilter!=nullptr) hi.sigmaintegralFilter = **sigmaintegralFilter;
      if (amplitude!=nullptr) hi.amplitude = **amplitude;
      if (amplitudeFilter!=nullptr) hi.amplitudeFilter = **amplitudeFilter;
      if (peaktick!=nullptr) hi.peaktick = **peaktick;
      if (peaktickFilter!=nullptr) hi.peaktickFilter = **peaktickFilter;
      if (peaktime!=nullptr) hi.peaktime = **peaktime;
      if (peaktimeFilter!=nullptr) hi.peaktimeFilter = **peaktimeFilter;
      if (begintick!=nullptr) hi.begintick = **begintick;
      if (endtick!=nullptr) hi.endtick = **endtick;
      if (width!=nullptr) hi.width = **width;
      if (hitx!=nullptr) hi.hitx = **hitx;
      if (hity!=nullptr) hi.hity = **hity;
      if (hitz!=nullptr) hi.hitz = **hitz;
      if (hiterrxlo!=nullptr) hi.hiterrxlo = **hiterrxlo;
      if (hiterrxhi!=nullptr) hi.hiterrxhi = **hiterrxhi;
      if (hiterrylo!=nullptr) hi.hiterrylo = **hiterrylo;
      if (hiterryhi!=nullptr) hi.hiterryhi = **hiterryhi;
      if (hiterrzlo!=nullptr) hi.hiterrzlo = **hiterrzlo;
      if (hiterrzhi!=nullptr) hi.hiterrzhi = **hiterrzhi;
      if (perpdist!=nullptr) hi.perpdist = **perpdist;
      if (hitt!=nullptr) hi.hitt = **hitt;
      if (driftdist!=nullptr) hi.driftdist = **driftdist;
      if (countercut!=nullptr) hi.countercut = **countercut;
      if (fitconstant!=nullptr) hi.fitconstant = **fitconstant;
      if (fitconstanterr!=nullptr) hi.fitconstanterr = **fitconstanterr;
      if (fitlinear!=nullptr) hi.fitlinear = **fitlinear;
      if (fitlinearerr!=nullptr) hi.fitlinearerr = **fitlinearerr;
      if (fitquadratic!=nullptr) hi.fitquadratic = **fitquadratic;
      if (fitquadraticerr!=nullptr) hi.fitquadraticerr = **fitquadraticerr;
      if (fitchi2!=nullptr) hi.fitchi2 = **fitchi2;
      if (fitsumsqrresidual!=nullptr) hi.fitsumsqrresidual = **fitsumsqrresidual;
      if (fitndf!=nullptr) hi.fitndf = **fitndf;
      if (fitmle!=nullptr) hi.fitmle = **fitmle;
      if (fitsuccess!=nullptr) hi.fitsuccess = **fitsuccess;
      if (fitrealhit!=nullptr) hi.fitrealhit = **fitrealhit;
      if (segmentlength!=nullptr) hi.segmentlength = **segmentlength;
      if (assumedhit!=nullptr) hi.assumedhit = **assumedhit;
      /*   if (numGoodHitsChan!=nullptr) hi.numGoodHitsChan = **numGoodHitsChan;
         if (nwiresTPC0!=nullptr) hi.nwiresTPC0 = **nwiresTPC0;
         if (nwiresTPC1!=nullptr) hi.nwiresTPC1 = **nwiresTPC1;
         if (nwiresTPC2!=nullptr) hi.nwiresTPC2 = **nwiresTPC2;
         if (nwiresTPC3!=nullptr) hi.nwiresTPC3 = **nwiresTPC3;
         if (nwiresTPC4!=nullptr) hi.nwiresTPC4 = **nwiresTPC4;
         if (nwiresTPC5!=nullptr) hi.nwiresTPC5 = **nwiresTPC5;
         if (nwiresTPC6!=nullptr) hi.nwiresTPC6 = **nwiresTPC6;
         if (nwiresTPC7!=nullptr) hi.nwiresTPC7 = **nwiresTPC7;
         if (prebaseline!=nullptr) hi.prebaseline = **prebaseline;
         if (postbaseline!=nullptr) hi.postbaseline = **postbaseline;
         if (prebaserms!=nullptr) hi.prebaserms = **prebaserms;
         if (postbaserms!=nullptr) hi.postbaserms = **postbaserms;
         if (trackid!=nullptr) hi.trackid = **trackid;
         if (numtrajpts!=nullptr) hi.numtrajpts = **numtrajpts;
         if (tracklength!=nullptr) hi.tracklength = **tracklength;
         if (isOnTrack!=nullptr) hi.isOnTrack = **isOnTrack;
         if (dqdxatpt!=nullptr) hi.dqdxatpt = **dqdxatpt;
         if (prevStartTick!=nullptr) hi.prevStartTick = **prevStartTick;
         if (prevEndTick!=nullptr) hi.prevEndTick = **prevEndTick;
         if (prevPeakTime!=nullptr) hi.prevPeakTime = **prevPeakTime;
         if (prevSigmaPeakTime!=nullptr) hi.prevSigmaPeakTime = **prevSigmaPeakTime;
         if (prevRMS!=nullptr) hi.prevRMS = **prevRMS;
         if (prevPeakAmplitude!=nullptr) hi.prevPeakAmplitude = **prevPeakAmplitude;
         if (prevSigmaPeakAmplitude!=nullptr) hi.prevSigmaPeakAmplitude = **prevSigmaPeakAmplitude;
         if (prevSummedADC!=nullptr) hi.prevSummedADC = **prevSummedADC;
         if (prevIntegral!=nullptr) hi.prevIntegral = **prevIntegral;
         if (prevSigmaIntegral!=nullptr) hi.prevSigmaIntegral = **prevSigmaIntegral;*/
      hits.emplace(std::make_pair(i,hi));
      i++;
      //std::cout << "hi.run=" << hi.run << std::endl;
      if (std::find(runs.begin(),runs.end(),hi.run) == runs.end()) runs.push_back(hi.run);
      //if (i>1e5) break;
    }
  std::cout << "ReadHistFile::ReadFile() -- Finished reading file. " << i << " hits were found in " << runs.size() << " runs." << std::endl;

  file->Close();
}

UInt_t ReadHistFile::ReadFile()
{
  return i;
}

#endif
