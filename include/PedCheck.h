#ifndef PEDCHECK_H
#define PEDCHECK_H

#include <fstream>
#include <glob.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream>

class PedCheck
{
public:
  PedCheck();

  Bool_t ReadPeds(std::vector<Int_t> runvec);

  Double_t PedRms(Int_t run, Int_t offlinechannel) {
    return ( (valid) ? runchanpedrms[run][offlinechannel].second : -99999 );
  }

  Double_t PedMean(Int_t run, Int_t offlinechannel) {
    return ( (valid) ? runchanpedrms[run][offlinechannel].first : -99999);
  }

private:
  Bool_t ReadPed(Int_t r);

  Bool_t valid;
  std::map<Int_t,Int_t> onlineoffline;
  std::map<Int_t,std::map<Int_t,std::pair<Double_t,Double_t> > > runchanpedrms;
};

PedCheck::PedCheck()
  : valid(true)
{
}

Bool_t PedCheck::ReadPeds(std::vector<Int_t> runvec)
{
  valid = true;
  for (auto const & run : runvec)
    {
      valid = valid && ReadPed(run);
    }
  return valid;
}

Bool_t PedCheck::ReadPed(Int_t run)
{
  Int_t pedRun = -1;
  std::string ONLINE,OFFLINE;
  std::string line;
  std::string oofilename = "/home/mthiesse/PurityAnalysis/BadChannels/online_offline.txt";
  std::ifstream oofile(oofilename.c_str(),std::ifstream::in);
  if (oofile.is_open())
    {
      while(getline(oofile,line))
        {
          std::stringstream ss(line);
          ss >> ONLINE >> OFFLINE;

          onlineoffline[stoi(ONLINE)] = stoi(OFFLINE);
        }
    }
  else
    {
      std::cout << "PedCheck::ReadPed() -- " << oofilename << " not found." << std::endl;
      return false;
    }

  glob_t globbuf;
  glob ("/home/mthiesse/PurityAnalysis/BadChannels/Pedestals/online_databaseRun_*.csv",GLOB_TILDE,NULL,&globbuf);

  std::vector<std::pair<Int_t,std::string> > pedestalRuns;
  for (UInt_t i = 0; i < globbuf.gl_pathc; ++i)
    {
      std::string file = globbuf.gl_pathv[i];
      pedestalRuns.push_back(std::pair<Int_t,std::string>(stoi(file.substr(71,5)),file));
    }
  pedestalRuns.push_back(std::pair<Int_t,std::string>(19999,"")); //not sure of the highest useful run number...
  std::sort(pedestalRuns.begin(),pedestalRuns.end());

  std::string pedFileName;
  for (UInt_t j = 0; j < pedestalRuns.size()-1; j++)
    {
      if (run >= pedestalRuns.at(j).first && run < pedestalRuns.at(j+1).first)
        {
          pedRun = pedestalRuns.at(j).first;
          pedFileName = pedestalRuns.at(j).second;
        }
    }
  if (pedRun != -1)
    {
      std::string ln;
      std::ifstream pedfile(pedFileName,std::ifstream::in);
      if (pedfile.is_open())
        {
          while(getline(pedfile,ln))
            {
              std::stringstream ss(ln);
              std::string lineagain;
              std::vector<std::string> linevec;
              while(getline(ss,lineagain,','))
                {
                  std::string VAL;
                  std::stringstream ssagain(lineagain);
                  ssagain >> VAL;
                  linevec.push_back(VAL);
                }
              runchanpedrms[run][onlineoffline[stoi(linevec[0])]] = std::pair<Double_t,Double_t>(stod(linevec[1]),stod(linevec[2]));
            }
        }
      else
        {
          std::cout << "PedCheck::ReadPed() -- " << pedFileName << " not found." << std::endl;
          return false;
        }
    }
  return true;
}

#endif
