#ifndef READHISTFILEALT_H
#define READHISTFILEALT_H

#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeReader.h"

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

class ReadHistFileAlt {

public:

  ReadHistFileAlt();
  UInt_t ReadFile(std::string filename);

  types::HitMap * GetHitMap() {
    return &hits;
  }

  std::vector<Int_t> GetRunList() {
    return runs;
  }

private:
  types::HitMap hits;
  std::vector<Int_t> runs;
};

ReadHistFileAlt::ReadHistFileAlt()
{
}

UInt_t ReadHistFileAlt::ReadFile(std::string filename)
{
  TFile * file = TFile::Open(filename.c_str(),"READ");
  if (!file || file->IsZombie())
    {
      std::stringstream ss;
      ss << "ReadHistFileAlt::ReadFile() -- Input file is not read";
      throw std::runtime_error(ss.str());
    }

  std::cout << "ReadHistFileAlt::ReadFile() -- Reading file " << filename << std::endl;

  TTreeReader reader("trackhit/trackhit",file);

  TTreeReaderValue<Int_t> run(reader,"run");
  TTreeReaderValue<Int_t> event(reader,"event");
  TTreeReaderValue<Int_t> channel(reader,"channel");
  TTreeReaderValue<Int_t> wire(reader,"wire");
  TTreeReaderValue<Int_t> tpc(reader,"tpc");
  TTreeReaderValue<Int_t> signalsize(reader,"signalsize");
  TTreeReaderValue<Float_t> prebaseline(reader,"prebaseline");
  TTreeReaderValue<Float_t> postbaseline(reader,"postbaseline");
  TTreeReaderValue<Float_t> prebaserms(reader,"prebaserms");
  TTreeReaderValue<Float_t> postbaserms(reader,"postbaserms");
  TTreeReaderValue<Float_t> integral(reader,"integral");
  TTreeReaderValue<Float_t> sigmaintegral(reader,"sigmaintegral");
  TTreeReaderValue<Float_t> amplitude(reader,"amplitude");
  TTreeReaderValue<Float_t> peaktick(reader,"peaktick");
  TTreeReaderValue<Int_t> begintick(reader,"begintick");
  TTreeReaderValue<Int_t> endtick(reader,"endtick");
  TTreeReaderValue<Float_t> segmentlength(reader,"segmentlength");
  TTreeReaderValue<Int_t> trackid(reader,"trackid");
  TTreeReaderValue<Int_t> numtrajpts(reader,"numtrajpts");
  TTreeReaderValue<Double_t> tracklength(reader,"tracklength");
  TTreeReaderValue<Bool_t> isOnTrack(reader,"isOnTrack");
  TTreeReaderValue<Double_t> dqdxatpt(reader,"dqdxatpt");
  TTreeReaderValue<Int_t> prevStartTick(reader,"prevStartTick");
  TTreeReaderValue<Int_t> prevEndTick(reader,"prevEndTick");
  TTreeReaderValue<Float_t> prevPeakTime(reader,"prevPeakTime");
  TTreeReaderValue<Float_t> prevSigmaPeakTime(reader,"prevSigmaPeakTime");
  TTreeReaderValue<Float_t> prevRMS(reader,"prevRMS");
  TTreeReaderValue<Float_t> prevPeakAmplitude(reader,"prevPeakAmplitude");
  TTreeReaderValue<Float_t> prevSigmaPeakAmplitude(reader,"prevSigmaPeakAmplitude");
  TTreeReaderValue<Float_t> prevSummedADC(reader,"prevSummedADC");
  TTreeReaderValue<Float_t> prevIntegral(reader,"prevIntegral");
  TTreeReaderValue<Float_t> prevSigmaIntegral(reader,"prevSigmaIntegral");

  UInt_t i = 0;
  while (reader.Next())
    {
      types::HitInfo hi;
      hi.run = *run;
      hi.event = *event;
      hi.channel = *channel;
      hi.wire = *wire;
      hi.tpc = *tpc;
      hi.signalsize = *signalsize;
      hi.prebaseline = *prebaseline;
      hi.postbaseline = *postbaseline;
      hi.prebaserms = *prebaserms;
      hi.postbaserms = *postbaserms;
      hi.integral = *integral;
      hi.sigmaintegral = *sigmaintegral;
      hi.amplitude = *amplitude;
      hi.peaktick = *peaktick;
      hi.begintick = *begintick;
      hi.endtick = *endtick;
      hi.segmentlength = *segmentlength;
      hi.trackid = *trackid;
      hi.numtrajpts = *numtrajpts;
      hi.tracklength = *tracklength;
      hi.isOnTrack = *isOnTrack;
      hi.dqdxatpt = *dqdxatpt;
      hi.prevStartTick = *prevStartTick;
      hi.prevEndTick = *prevEndTick;
      hi.prevPeakTime = *prevPeakTime;
      hi.prevSigmaPeakTime = *prevSigmaPeakTime;
      hi.prevRMS = *prevRMS;
      hi.prevPeakAmplitude = *prevPeakAmplitude;
      hi.prevSigmaPeakAmplitude = *prevSigmaPeakAmplitude;
      hi.prevSummedADC = *prevSummedADC;
      hi.prevIntegral = *prevIntegral;
      hi.prevSigmaIntegral = *prevSigmaIntegral;
      hits.emplace(std::make_pair(i,hi));
      i++;
      if (std::find(runs.begin(),runs.end(),hi.run) == runs.end()) runs.push_back(hi.run);
    }
  std::cout << "ReadHistFileAlt::ReadFile() -- Finished reading file. " << i << " hits were found in " << runs.size() << " runs." << std::endl;
  return i;
}

#endif
