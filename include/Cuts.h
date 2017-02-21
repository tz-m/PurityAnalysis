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

  Bool_t ChannelPass(const types::HitInfo * hit);
  Bool_t HitPass(const types::HitInfo * hit);
  Bool_t CounterPass(const types::HitInfo * hit);

  Bool_t GapWire(Int_t channel);
  Bool_t PedestalCut(Int_t run, Int_t channel);

  Float_t GetAnalysisCutFloat(std::string ananame, std::string parname);
  Bool_t GetAnalysisCutBool(std::string ananame, std::string parname);
  std::string GetAnalysisCutString(std::string ananame, std::string parname);

  PedCheck * GetPedCheckPtr() {
    return pc;
  }

private:
  PedCheck * pc;

  types::AnalysisCuts_F anacutsfloat;
  types::AnalysisCuts_B anacutsbool;
  types::AnalysisCuts_S anacutsstring;
};

Cuts::Cuts()
  : pc(new PedCheck())
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
              if (PAR == "true") anacutsbool[NAMEOFANA][NAMEOFPAR] = true;
              if (PAR == "false") anacutsbool[NAMEOFANA][NAMEOFPAR] = false;
            }
          if (TYPEOFPAR == "float")
            {
              anacutsfloat[NAMEOFANA][NAMEOFPAR] = stof(PAR);
            }
          if (TYPEOFPAR == "string")
            {
              anacutsstring[NAMEOFANA][NAMEOFPAR] = PAR;
            }
        }
    }
  else
    {
      std::cout << "Cuts::PrepareCuts() -- " << configfile << " not opened." << std::endl;
      valid=false;
    }

  //if (GetAnalysisCutBool("All","CutPedestals"))
  {
    valid = valid && pc->ReadPeds(runvec);
  }

  if (!valid) std::cout << "Problem with setting up cuts" << std::endl;
  return valid;
}

Float_t Cuts::GetAnalysisCutFloat(std::string ananame, std::string parname)
{
  if (anacutsfloat[ananame].find(parname) == anacutsfloat[ananame].end())
    {
      std::stringstream ss;
      ss << "Analysis \"" << ananame << "\" or parameter \"" << parname << "\" was not found";
      throw std::runtime_error(ss.str());
    }
  return anacutsfloat[ananame][parname];
}

Bool_t Cuts::GetAnalysisCutBool(std::string ananame, std::string parname)
{
  if (anacutsbool[ananame].find(parname) == anacutsbool[ananame].end())
    {
      std::stringstream ss;
      ss << "Analysis \"" << ananame << "\" or parameter \"" << parname << "\" was not found";
      throw std::runtime_error(ss.str());
    }
  return anacutsbool[ananame][parname];
}

std::string Cuts::GetAnalysisCutString(std::string ananame, std::string parname)
{
  if (anacutsstring[ananame].find(parname) == anacutsstring[ananame].end())
    {
      std::stringstream ss;
      ss << "Analysis \"" << ananame << "\" or parameter \"" << parname << "\" was not found";
      throw std::runtime_error(ss.str());
    }
  return anacutsstring[ananame][parname];
}

Bool_t Cuts::ChannelPass(const types::HitInfo * hit)
{
  Bool_t pedestalpass = !PedestalCut(hit->run,hit->channel);
  Bool_t gapwirepass = !GapWire(hit->channel);
  return pedestalpass && gapwirepass;
}

Bool_t Cuts::HitPass(const types::HitInfo * hit)
{
  Bool_t realhitpass = (GetAnalysisCutBool("All","FitRealHit") == hit->fitrealhit);
  Bool_t consistenttime = true; //fabs(hit->peaktime - hit->peaktimeFilter) < 5;
  return realhitpass && consistenttime;
}

Bool_t Cuts::CounterPass(const types::HitInfo * hit)
{
  Bool_t ewopposite = !(GetAnalysisCutBool("All","EastWestOpposite") && !(((hit->c1>=6 && hit->c1<=15) || (hit->c1>=28 && hit->c1<=37)) && (hit->c1%22==hit->c2 || hit->c2%22==hit->c1)));
  return ewopposite;
}

Bool_t Cuts::PedestalCut(Int_t run, Int_t channel)
{
  Float_t m = pc->PedMean(run,channel);
  Float_t r = pc->PedRms(run,channel);
  Float_t minRms = GetAnalysisCutFloat("All","MinPedRms");
  Float_t maxRms = GetAnalysisCutFloat("All","MaxPedRms");
  Float_t minMean = GetAnalysisCutFloat("All","MinPedMean");
  Float_t maxMean = GetAnalysisCutFloat("All","MaxPedMean");
  return GetAnalysisCutBool("All","CutPedestals") && (r<minRms || r>maxRms || m<minMean || m>maxMean);
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
