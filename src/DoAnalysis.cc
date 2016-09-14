#include "MakePlotFromTree.h"
#include "GraphChargeDistributions.h"

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
      std::cout << "  -f : Input data filename" << std::endl;
      std::cout << std::endl;
      std::cout << "To choose an analysis: -a [ANALYSIS NAME]" << std::endl;
      std::cout << std::endl;
      std::cout << "Possible analyses to run:" << std::endl;
      std::cout << "   MakePlotFromTree : Read the tree and" << std::endl;
      std::cout << "       plot a variable. default: hit amplitude" << std::endl;
      std::cout << "   GraphChargeDistributions : Make plots of" << std::endl;
      std::cout << "       charge distributions of \"real\" and" << std::endl;
      std::cout << "       \"fake\" hits" << std::endl;
      gApplication->Terminate(0);
      return 0;
    }

  const std::string & filename = input.getCmdOption("-f");
  if (filename.empty())
    {
      std::cout << "DoAnalysis -- Can't do anything without an input filename!" << std::endl;
      gApplication->Terminate(1);
      return 1;
    }

  const std::string & analysis = input.getCmdOption("-a");
  if (!analysis.empty())
    {
      if (analysis == "MakePlotFromTree")
        {
          std::cout << "DoAnalysis -- Making plot from TTree in file " << filename << std::endl;
          MakePlotFromTree mpft;
          mpft.Setup(filename);
          mpft.run();
        }
      if (analysis == "GraphChargeDistributions")
        {
          std::cout << "DoAnalysis -- Graphing Charge Distributions" << std::endl;
          GraphChargeDistributions gcds;
          gcds.Setup(filename);
          gcds.run();
        }
    }

  app.Run();
  return 0;
}
