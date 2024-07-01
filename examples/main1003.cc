// Authors: Marius Utheim <marius.m.utheim@jyu.fi>;
//          Chloé Gaudu <gaudu@uni-wuppertal.de>

// Keywords: analysis; rivet; hA collisions

// This example simulates events for specific particle beam, energy and frame settings.
// These events are analyzed using Rivet, set to ignore messages regarding beam 
// properties. The Rivet configuration includes adding a specific analysis and specifying
// the YODA output file. A log file is printed with Pythia statistics after the run.

#include "Pythia8/HeavyIons.h"
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/Pythia8Rivet.h"

using namespace Pythia8;

int main() {

  // Run settings.
  int nEvents = 1e6;
  constexpr int nCases = 1; //2; 
  vector<string> caseLabels = { "158" }; //, "350" };
  string out = "main1003_-211_1000060120_frameType2_158";
  
  for (int iCase = 0; iCase < nCases; ++iCase) {

    // Pythia configuration.
    Pythia pythia;
    pythia.readString("Beams:idA = -211");
    pythia.readString("Beams:idB = 1000060120");
    pythia.readString("Beams:frameType = 2");
    if (iCase == 0)
      pythia.readString("Beams:eA = 158.");
    else
    	pythia.readString("Beams:eA = 350.");
    pythia.readString("Beams:eB = "+std::to_string(pythia.particleData.m0(2212))); // mass_proton ~= mass_carbon/12

    // Parameters to reuse initialization.
    pythia.readString("MultipartonInteractions:reuseInit = 3");
    pythia.readString("HeavyIon:SasdMpiReuseInit = 3");
    pythia.readString("HeavyIon:SigFitReuseInit = 3");
    string initPrefix = "main1003-"+caseLabels[iCase];
    pythia.readString("MultipartonInteractions:initFile = "+initPrefix+".mpi");
    pythia.readString("HeavyIon:SasdMpiInitFile = "+initPrefix+".sasd.mpi");
    pythia.readString("HeavyIon:SigFitInitFile = "+initPrefix+".sigfit");

    if (!pythia.init()) {
      cout << "Pythia failed to initialize." << endl;
      return -1;
    }

    // Rivet initialization.
    Pythia8Rivet rivet(pythia, out + ".yoda");
    rivet.ignoreBeams(true);
    rivet.addAnalysis("NA61SHINE_2022_I2155140");

    // Loop over events.
    for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) { 
      if (!pythia.next()) continue;
      rivet();
    }
    rivet.done();
    pythia.stat();

  }

  return 0;
}
