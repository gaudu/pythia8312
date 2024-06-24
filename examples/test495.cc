// test495.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Analysis of target-remnant issues in Angantyr,
// at variable beams and energies. (Requires main424 to be run first.)

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events. Beams (p, pi+, pi-) + (14N, 16O, 40Ar). Energy range.
  int nEvent       = 900;
  int idBeams[3]   = {2212, 211, -211};
  int idTargets[3] = {1000070140, 1000080160, 1000180400};
  double eBeamMn = 1e6;
  double eBeamMx = 1e7;

  // Create Angantyr generator. Shorthand for event.
  Pythia pythia;
  pythia.readString("HeavyIon:mode = 2"); 
  Event& event = pythia.event;

  // Beam type and energy values for initialization
  pythia.settings.mode("Beams:idA", idBeams[0]);
  pythia.settings.mode("Beams:idB", idTargets[2]);
  pythia.readString("Beams:frameType = 2");
  pythia.settings.parm("Beams:eA", eBeamMx);
  double mTarg = pythia.particleData.m0(2212);
  pythia.settings.parm("Beams:eB", mTarg);

  // Variable incoming beam type and energy.
  pythia.readString("Beams:allowVariableEnergy = on");
  pythia.readString("Beams:allowIDAswitch = on");

  // Reuse MPI and heavy-ion geometry initialization files. 
  // If they don't exist, you can generate them by running main424.
  pythia.readString("MultipartonInteractions:reuseInit = 2");
  pythia.readString("MultipartonInteractions:initFile = main424.mpi");
  pythia.readString("HeavyIon:SasdMpiReuseInit = 2");
  pythia.readString("HeavyIon:SasdMpiInitFile = main424.sasd.mpi");
  pythia.readString("HeavyIon:SigFitReuseInit = 2");
  pythia.readString("HeavyIon:SigFitInitFile = main424.sigfit");

  // Simplified generation for debug purposes.
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:all = off");

  // No check on energy-momentum conservation.
  pythia.readString("Check:event = off");

  // Exit with error if Pythia fails to initialize.
  if (!pythia.init()) return 1;

  // Book histograms for momentum checks.
  Hist Arem( "remnant atomic number A", 50, -0.5, 49.5);
  Hist pTot( "momentum fraction full final state", 100, 0.8, 1.2);
  Hist pRem( "momentum fraction target remnant",   100, 0.0, 1.0);
  Hist pHad( "momentum fraction hadronic system",  100, 0.8, 1.2);

  // Loop over events.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Variable incoming beam type and energy, and target type.
    double eBeam = eBeamMn + pythia.rndm.flat() * (eBeamMx - eBeamMn);
    int idBeam   = idBeams[(iEvent/3)%3]; 
    int idTarg   = idTargets[iEvent%3]; 
    pythia.setBeamIDs( idBeam, idTarg);
    pythia.setKinematics( eBeam, mTarg);

    // Simulate interaction.
    pythia.next();
    if (iEvent < 9) event.list();

    // The event weight. For test more useful to histogram by number.
    //double weight = pythia.info.weight();
    double weight = 1.;

    // Check final momentum and catch target remnant. 
    Vec4 p4Tot;
    int iTargRemn  = 0;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {
      p4Tot += event[i].p();
      if (event[i].id() > 1000000000) iTargRemn = i; 
    }

    // A of remnant. Skip if no nuclear remnant or if elastic scattering.
    int Atarg = (idTarg / 10) % 1000; 
    int Aremn = (iTargRemn > 0) ? (event[iTargRemn].id() / 10) % 1000 : 0;
    Arem.fill( Aremn, weight);
    if (iTargRemn == 0 || Aremn == Atarg) continue;

    // Check how off the final-state momentum is.
    Vec4 p4Rem = event[iTargRemn].p();
    Vec4 p4Had  = p4Tot - p4Rem;
    pTot.fill( p4Tot.pAbs() / eBeam, weight);
    pRem.fill( p4Rem.pAbs() / eBeam, weight);
    pHad.fill( p4Had.pAbs() / eBeam, weight);

  // End of event loop.  
  }

  // Print statistics, mainly for errors. Also histograms.
  pythia.stat();
  cout << Arem << pTot << pRem << pHad;

}
