#ifndef CUTS_H
#define CUTS_H

#include "PedCheck.h"

#include <fstream>
#include <string>
#include <algorithm>
#include <map>
#include <sstream>

class Cuts {
public:
  Cuts();

  Bool_t PrepareCuts(std::string configfile, std::vector<Int_t> runvec);

  Bool_t ChannelPass(const types::HitInfo & hit);
  Bool_t HitPass(const types::HitInfo & hit);
  Bool_t CounterPass(const types::HitInfo & hit);

  Bool_t GapWire(Int_t channel);
  Bool_t PedestalCut(Int_t run, Int_t channel);

  Float_t GetAnalysisCut(std::string ananame, std::string cutname);

private:
  Float_t minPedMean;
  Float_t maxPedMean;
  Float_t minPedRms;
  Float_t maxPedRms;

  Bool_t cutpedestals;
  Bool_t fitrealhit;
  Bool_t oppcounter;

  PedCheck pc;

  std::map<std::string,std::map<std::string,Float_t> > cutmapfloat;
};

Cuts::Cuts()
  : minPedMean(-99999),
    maxPedMean(-99999),
    minPedRms(-99999),
    maxPedRms(-99999),
    cutpedestals(false),
    fitrealhit(true),
    oppcounter(true)
{

}

Bool_t Cuts::PrepareCuts(std::string configfile, std::vector<Int_t> runvec)
{
  Bool_t valid = true;
  std::string TYPEOFPAR, NAMEOFANA, NAMEOFPAR, PAR;
  std::string line;
  std::ifstream file(configfile.c_str(),std::ifstream::in);
  if (file.is_open())
    {
      while (getline(file,line))
        {
          std::stringstream ss(line);
          ss >> TYPEOFPAR >> NAMEOFANA >> NAMEOFPAR >> PAR;

          if (TYPEOFPAR == "bool")
            {
              if (NAMEOFPAR == "CutPedestals")
                {
                  if (PAR == "true") cutpedestals = true;
                  else if (PAR == "false") cutpedestals = false;
                }
              if (NAMEOFPAR == "FitRealHit")
                {
                  if (PAR == "true") fitrealhit = true;
                  else if (PAR == "false") fitrealhit = false;
                }
              if (NAMEOFPAR == "EastWestOpposite")
                {
                  if (PAR == "true") oppcounter = true;
                  else if (PAR == "false") oppcounter = false;
                }
            }
          else if (TYPEOFPAR == "float")
            {
              if (NAMEOFPAR == "MinPedMean") minPedMean = stod(PAR);

              if (NAMEOFPAR == "MaxPedMean") maxPedMean = stod(PAR);

              if (NAMEOFPAR == "MinPedRms") minPedRms  = stod(PAR);

              if (NAMEOFPAR == "MaxPedRms") maxPedRms  = stod(PAR);

              if (NAMEOFPAR == "MakePlotFromTree")
                {

                }
            }
        }
    }
  else
    {
      std::cout << "Cuts::PrepareCuts() -- " << configfile << " not opened." << std::endl;
      valid=false;
    }

  if (cutpedestals)
    {
      valid = valid && pc.ReadPeds(runvec);
    }

  return valid;
}

Bool_t Cuts::ChannelPass(const types::HitInfo & hit)
{
  Bool_t pedestalpass = !PedestalCut(hit.run,hit.channel);
  Bool_t gapwirepass = !GapWire(hit.channel);
  return pedestalpass && gapwirepass;
}

Bool_t Cuts::HitPass(const types::HitInfo & hit)
{
  Bool_t realhitpass = (fitrealhit == hit.fitrealhit);
  return realhitpass;
}

Bool_t Cuts::CounterPass(const types::HitInfo & hit)
{
  Bool_t ewopposite = !(oppcounter && !((hit.c1>=6 && hit.c1<=15) || (hit.c1>=28 && hit.c1<=37)) && (hit.c1%22==hit.c2 || hit.c2%22==hit.c1));
  return ewopposite;
}

Bool_t Cuts::PedestalCut(Int_t run, Int_t channel)
{
  Double_t m = pc.PedMean(run,channel);
  Double_t r = pc.PedRms(run,channel);
  return cutpedestals && (r<minPedRms || r>maxPedRms || m<minPedMean || m>maxPedMean);
}

Bool_t Cuts::GapWire(Int_t channel)
{
  return (channel == 288 || channel == 399 ||
          channel == 400 || channel == 511 ||
          channel == 800 || channel == 911 ||
          channel == 912 || channel == 1023 ||
          channel == 1312 || channel == 1423 ||
          channel == 1424 || channel == 1535 ||
          channel == 1824 || channel == 1935 ||
          channel == 1936 || channel == 2047);
}

#endif
