#include "MakePlotFromTree.h"
#include "GraphChargeDistributions.h"
#include "GraphAvgHitsWires.h"

#include "InputParser.h"

#include "TApplication.h"
#include "TROOT.h"

int main(int argc, char** argv)
{
  TApplication app("hitvalidation",0,(char**)"");
  app.ExitOnException();

  InputParser input(argc,argv);
  if (input.cmdOptionExists("-h"))
    {
      std::cout << "Required options: " << std::endl;
      std::cout << "  -f [FILENAME] : Input data file" << std::endl;
      std::cout << "  -c [FILENAME] : Configuration file" << std::endl;
      std::cout << "  -a [ANALYSISNAME] : Name of analysis to run" << std::endl;
      std::cout << std::endl;
      std::cout << "Available analyses:" << std::endl;
      std::cout << std::endl;
      std::cout << "   MakePlotFromTree : Read the tree and" << std::endl;
      std::cout << "       plot a variable. default: hit amplitude" << std::endl;
      std::cout << "   GraphChargeDistributions : Make plots of" << std::endl;
      std::cout << "       charge distributions of \"real\" and" << std::endl;
      std::cout << "       \"fake\" hits" << std::endl;
      gApplication->Terminate(0);
      return 0;
    }

  const std::string & datafile = input.getCmdOption("-f");
  if (datafile.empty())
    {
      std::cout << "DoAnalysis -- Input file is required. Use \"-f [FILENAME]\"" << std::endl;
      gApplication->Terminate(1);
      return 1;
    }

  const std::string & configfile = input.getCmdOption("-c");
  if (configfile.empty())
    {
      std::cout << "DoAnalysis -- Configuration file is required. Use \"-c [FILENAME]\"" << std::endl;
      gApplication->Terminate(1);
      return 1;
    }

  const std::string & analysis = input.getCmdOption("-a");
  if (!analysis.empty())
    {
      if (analysis == "MakePlotFromTree")
        {
          std::cout << "DoAnalysis -- Making plot from TTree in file " << datafile << " using configuration in " << configfile << std::endl;
          MakePlotFromTree mpft;
          mpft.Setup(datafile,configfile);
          mpft.run();
        }
      else if (analysis == "GraphChargeDistributions")
        {
          std::cout << "DoAnalysis -- Graphing Charge Distributions using data file " << datafile << " and configuration in " << configfile << std::endl;
          GraphChargeDistributions gcds;
          gcds.Setup(datafile,configfile);
          gcds.run();
        }
      else if (analysis == "GraphAvgHitsWires")
        {
          std::cout << "DoAnalysis -- Graphing channel efficiency" << std::endl;
          GraphAvgHitsWires gahw;
          gahw.Setup(datafile,configfile);
          gahw.run();
        }
      else
        {
          std::cout << "DoAnalysis -- UNKNOWN ANALYSIS" << std::endl;
          gApplication->Terminate(1);
          return 1;
        }
    }
  else
    {
      std::cout << "DoAnalysis -- Nothing to do. Provide a directive with the \"-a [ANALYSISNAME]\" option" << std::endl;
      gApplication->Terminate(1);
      return 1;
    }

  app.Run();
  return 0;
}
