// test496.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Analysis of target-remnant issues in PythiaCascade
// at variable beams and energies. (Requires main424 to be run first.)

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PythiaCascade.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events. Beams (p, pi+, pi-) + (14N, 16O, 40Ar). Energy range.
  int nEvent       = 900;
  int idBeams[3]   = {2212, 211, -211};
  int idTargets[3] = {1000070140, 1000080160, 1000180400};
  double eBeamMn = 1e6;
  double eBeamMx = 1e7;

  // Create PythiaCascade generator. Shorthand for event.
  PythiaCascade pythia;

  // Exit with error if Pythia(Cascade) fails to initialize.
  if (!pythia.init()) return 1;

  // Book histograms for momentum checks.
  Hist Arem( "remnant atomic number A", 50, -0.5, 49.5);
  Hist pTot( "momentum fraction full final state", 100, 0.8, 1.2);

  // Loop over events.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Variable incoming beam type and energy, and target type.
    int    idBeam = idBeams[(iEvent/3)%3]; 
    int    idTarg = idTargets[iEvent%3]; 
    double eBeam  = eBeamMn + pythia.rndm().flat() * (eBeamMx - eBeamMn);
    double mBeam  = pythia.particleData().m0(idBeam);
    Vec4  p4Beam( 0., 0., sqrtpos( pow2(eBeam) - pow2(mBeam)), eBeam); 

    // Simulate interaction.
    pythia.sigmaSetuphN( idBeam, p4Beam, mBeam);
    int Ztarg = (idTarg / 10000) % 1000;
    int Atarg = (idTarg / 10) % 1000;  
    Event& event = pythia.nextColl( Ztarg, Atarg);
    //if (iEvent < 9) event.list();

    // The event weight. Always unity in PythiaCascade.
    double weight = 1.;

    // Check final momentum. Count number of kicked-out nucleons.
    Vec4 p4Tot;
    int nTargKicked = 0;
    for (int i = 0; i < event.size(); ++i) {
      if (event[i].isFinal()) p4Tot += event[i].p();
      if (event[i].statusAbs() == 181 || event[i].statusAbs() == 182) 
        ++nTargKicked;
    }

    // A of remnant.
    int Aremn = Atarg - nTargKicked;
    Arem.fill( Aremn, weight);

    // Check how off the final-state momentum is.
    pTot.fill( p4Tot.pAbs() / eBeam, weight);

  // End of event loop.  
  }

  // Print statistics, mainly for errors. Also histograms.
  pythia.stat();
  cout << Arem << pTot;

}
