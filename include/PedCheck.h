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
  PedCheck(int r);

  bool check(int ch);
  void setMinMaxRMS(double min, double max) {
    minRMS = min; maxRMS = max;
  }
  void setMinMaxPed(double min, double max) {
    minPed = min; maxPed = max;
  }

private:
  int run;
  int pedRun;
  std::string pedFileName;
  double minRMS;
  double maxRMS;
  double minPed;
  double maxPed;
  bool valid;
  std::map<int,int> onlineoffline;
  std::map<int,int> offlineonline;
  std::map<int,std::pair<double,double> > chanpedrms;
};

PedCheck::PedCheck()
  : run(0),
    pedRun(-1),
    valid(false)
{
}

PedCheck::PedCheck(int r)
  : run(r),
    pedRun(-1),
    valid(false)
{
  std::string ONLINE,OFFLINE;
  std::string line;
  std::ifstream oofile("/home/mthiesse/Documents/BadChannel/online_offline.txt",std::ifstream::in);
  if (oofile.is_open())
    {
      while(getline(oofile,line))
        {
          std::stringstream ss(line);
          ss >> ONLINE >> OFFLINE;

          onlineoffline[stoi(ONLINE)] = stoi(OFFLINE);
          offlineonline[stoi(OFFLINE)] = stoi(ONLINE);
        }
    }

  glob_t globbuf;
  glob ("/home/mthiesse/Documents/BadChannel/online_databaseRun_*.csv",GLOB_TILDE,NULL,&globbuf);

  std::vector<std::pair<int,std::string> > pedestalRuns;
  for (unsigned int i = 0; i < globbuf.gl_pathc; ++i)
    {
      std::string file = globbuf.gl_pathv[i];
      pedestalRuns.push_back(std::pair<int,std::string>(stoi(file.substr(55,5)),file));
    }
  pedestalRuns.push_back(std::pair<int,std::string>(19999,"")); //not sure of the highest useful run number...
  std::sort(pedestalRuns.begin(),pedestalRuns.end());

  for (unsigned int j = 0; j < pedestalRuns.size()-1; j++)
    {
      if (run >= pedestalRuns.at(j).first && run < pedestalRuns.at(j+1).first)
        {
          pedRun = pedestalRuns.at(j).first;
          pedFileName = pedestalRuns.at(j).second;
        }
    }
  if (pedRun != -1)
    {
      valid = true;
      //std::cout << "Run " << run << " uses pedestal file " << pedFileName << std::endl;

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
              chanpedrms[onlineoffline[stoi(linevec[0])]] = std::pair<double,double>(stod(linevec[1]),stod(linevec[2]));
            }
        }
    }
}

// OFFLINE channel number as argument
bool PedCheck::check(int ch)
{
  if (!valid) return false;

  //std::cout << "Online channel: " << offlineonline[ch] << "  Offline channel: " << ch << "  pedestal: " << chanpedrms[ch].first << "  rms: " << chanpedrms[ch].second << std::endl;

  if (chanpedrms[ch].first < minPed || chanpedrms[ch].first > maxPed || chanpedrms[ch].second < minRMS || chanpedrms[ch].second > maxRMS) return false;

  return true;
}


// requires 2 arguments:
//     1. current run number (program will find the correct pedestal run)
//     2. offline channel number
/*
   int main(int argc, char** argv)
   {
   PedCheck pc(stoi(argv[1]));
   pc.setMinMaxPed(200,1500);
   pc.setMinMaxRMS(6,40);

   int goodwires = 0;
   int totalwires = 0;
   std::string CHAN,TPC,PLANE,WIRE;
   std::string line;
   std::ifstream chanfile("/dune/app/users/mthiesse/PersistentFiles/CollectionWireList.txt",std::ifstream::in);
   if (chanfile.is_open())
    {
      while(getline(chanfile,line))
        {
          totalwires++;
          std::stringstream ss(line);
          ss >> CHAN >> TPC >> PLANE >> WIRE;
          if (pc.check(stoi(CHAN))) goodwires++;
        }
    }

   std::cout << "Good Wires: " << goodwires << "  Total Wires: " << totalwires << std::endl;

   return 0;
   }
 */
