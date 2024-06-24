// test494.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Analysis of target-remnant issues in Angantyr,
// at variable energies. (Requires main424 to be run first.)

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events. Beams p + 14N.
  int nEvent     = 1000;
  int idBeam     = 2212;
  int idTarget   = 1000070140;
  // Fixed-energy array and variable-energy range.
  double eBeamVal[5] = {1e6, 1e7, 1e8, 1e9, 1e10};
  double eBeamMn = 1e9;
  double eBeamMx = 1e10;

  // Begin loop over energy options.
  for (int iE = 0; iE < 6; ++iE) {
    bool varyE   = (iE == 5);
    double eBeam = (iE < 5) ? eBeamVal[iE] : eBeamMx; 

    // Create Angantyr generator. Shorthand for event.
    Pythia pythia;
    pythia.readString("HeavyIon:mode = 2"); 
    Event& event = pythia.event;

    // Setup the beams: p + 14N.
    pythia.settings.mode("Beams:idA", idBeam);
    pythia.settings.mode("Beams:idB", idTarget);
    pythia.readString("Beams:frameType = 2");
    pythia.settings.parm("Beams:eA", eBeam);
    double mTarg = pythia.particleData.m0(2212);
    pythia.settings.parm("Beams:eB", mTarg);

    // Variable incoming beam type and energy.
    if (varyE) {
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
    }

    // Simplified generation for debug purposes.
    pythia.readString("PartonLevel:ISR = off");
    pythia.readString("PartonLevel:FSR = off");
    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("HadronLevel:all = off");

    // No check on energy-momentum conservation.
    pythia.readString("Check:event = off");

    // Reduced printout.
    pythia.readString("Init:showProcesses  = off");
    pythia.readString("Init:showMultipartonInteractions = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showChangedParticleData = off");
    
    // Exit with error if Pythia fails to initialize.
    if (!pythia.init()) return 1;

    // Book histograms for momentum checks.
    double mHadMax = sqrt(2. * eBeam * 14. * mTarg);
    Hist Arem( "remnant atomic number A",            20, -0.5, 19.5);
    Hist pTot( "momentum fraction full final state", 100, 0.8, 1.2);
    Hist pRem( "momentum fraction target remnant",   100, 0.0, 1.0);
    Hist pHad( "momentum fraction hadronic system",  100, 0.8, 1.2);
    Hist mHad( "invariant mass hadronic system",     100, 0., mHadMax); 

    // Loop over events.
    int nSmallPT = 0;
    int nLargePT = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Variable incoming beam type and energy.
      if (varyE) {
        eBeam = eBeamMn + pythia.rndm.flat() * (eBeamMx - eBeamMn);

        // Set collision kinematics.
        pythia.setBeamIDs(idBeam, idTarget);
        pythia.setKinematics( eBeam, mTarg);
      }

      // Simulate interaction.
      pythia.next();

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
      int Aremn = (iTargRemn > 0) ? (event[iTargRemn].id() / 10) % 1000 : 0;
      Arem.fill( Aremn, weight);
      if (iTargRemn == 0 || Aremn == 14) continue;
      Vec4 p4Rem = event[iTargRemn].p();

      // Sum up hadronic four-momentum, evading target remnant.
      Vec4 p4Had;
      for (int i = 0; i < event.size(); ++i)
        if (event[i].isFinal() && i != iTargRemn) p4Had += event[i].p();

      // Check how off the final-state momentum is. Also mass.
      pTot.fill( p4Tot.pAbs() / eBeam, weight);
      pRem.fill( p4Rem.pAbs() / eBeam, weight);
      pHad.fill( p4Had.pAbs() / eBeam, weight);
      mHad.fill( p4Had.mCalc(), weight);

      // Identify events with small or large pT. Print first of each  kind.
      bool hasSmallPT = (event[iTargRemn].pT() < 1e-3);
      if (hasSmallPT) ++nSmallPT;
      bool hasLargePT = (event[iTargRemn].pT() > 1.);
      if (hasLargePT) ++nLargePT;
      if ((hasSmallPT && nSmallPT == 1) || (hasLargePT && nLargePT == 1)) 
        event.list();

      // Print key quantities for a few events, with rescaled pz and E.
      if ((hasSmallPT && nSmallPT < 20) || (hasLargePT && nLargePT < 20)) 
      cout << "\n Output (px, py, pz / eBeam, e / eBeam, m) for event no "
           << iEvent << " with eBeam = " << setprecision(5) << scientific
           << setw(14) << eBeam << " and Aremn = " << Aremn << endl  
           << " Total    final state " << fixed << setw(10) << p4Tot.px()
           << setw(10) << p4Tot.py() << setw(10) << p4Tot.pz() / eBeam
           << setw(10) << p4Tot.e() / eBeam << scientific << setw(14)
           << p4Tot.mCalc() << endl
           << " Remnant  final state " << fixed << setw(10) << p4Rem.px()
           << setw(10) << p4Rem.py() << setw(10) << p4Rem.pz() / eBeam
           << setw(10) << p4Rem.e() / eBeam << scientific << setw(14)
           << p4Rem.mCalc() << endl
           << " Hadronic final state " << fixed << setw(10) << p4Had.px()
           << setw(10) << p4Had.py() << setw(10) << p4Had.pz() / eBeam
           << setw(10) << p4Had.e() / eBeam << scientific << setw(14)
           << p4Had.mCalc() << endl;
    
    // End of event loop.  
    }

    // Print statistics, mainly for errors. Also histograms.
    pythia.stat();
    cout << Arem << pTot << pRem << pHad << mHad;

  // End of loop over energy options.
  }

}
