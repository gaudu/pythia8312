// test490.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim <marius.m.utheim@jyu.fi>
//          Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: cosmic ray cascade; switch beam; switch collision energy

// Test of afterburner to patch up Angantyr for unwanted beam remnant
// and poor momentum conservation.
// Note that this example requires main424 to be run first to produce 
// the necessary initialization files.
// Warning: this is currently a tryout code.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Event number and (variable) energy.
  int nEvent     = 100;
  bool varyE     = false;
  double pBeamMx = 1e10;
  double pBeamMn = 1e9;

  // Pythia object for performing individual collisions.
  Pythia pythiaColl;
  Event& eventColl = pythiaColl.event;

  // Beam and target properties.
  int idBeam   = 2212;
  int idTarg   = 1000070140;  
  int Atarg    = (idTarg / 10) % 1000;
  double mBeam = pythiaColl.particleData.m0(idBeam);
  double mTarg = pythiaColl.particleData.m0(2212);

  // Incoming kinematics.
  double pBeam = pBeamMx;
  double eBeam = sqrt(pow2(pBeam) + pow2(mBeam));
  Vec4 p4Beam( 0., 0., pBeam, eBeam);
  Vec4 p4Targ( 0., 0., 0., mTarg);
  double mCM   = ( p4Beam + p4Targ).mCalc();

  // Enable Angantyr.
  pythiaColl.readString("HeavyIon:mode = 2"); 
 
  // Variable incoming beam type and energy.
  pythiaColl.readString("Beams:allowVariableEnergy = on");
  pythiaColl.readString("Beams:allowIDAswitch = on");
  
  // Set up for collisions frame and energy.
  pythiaColl.readString("Beams:frameType = 2");
  pythiaColl.settings.parm("Beams:eA", eBeam);
  pythiaColl.settings.parm("Beams:eB", mTarg);
  
  // Simplified generation for debug purposes.
  pythiaColl.readString("PartonLevel:ISR = off");
  pythiaColl.readString("PartonLevel:FSR = off");
  pythiaColl.readString("PartonLevel:MPI = off");
  pythiaColl.readString("HadronLevel:Hadronize = off");

  // Decays should be handled separately.
  pythiaColl.readString("HadronLevel:Decay = off");
  
  // Reduced energy-momentum conservation check (hiding problems here!).
  pythiaColl.readString("Check:epTolErr = 0.1");
  
  // Reuse MPI and heavy-ion geometry initialization file. 
  // If they don't exist, you can generate them by running main424.
  pythiaColl.readString("MultipartonInteractions:reuseInit = 2");
  pythiaColl.readString("MultipartonInteractions:initFile = main424.mpi");
  pythiaColl.readString("HeavyIon:SasdMpiReuseInit = 2");
  pythiaColl.readString("HeavyIon:SasdMpiInitFile = main424.sasd.mpi");
  pythiaColl.readString("HeavyIon:SigFitReuseInit = 2");
  pythiaColl.readString("HeavyIon:SigFitInitFile = main424.sigfit");
  
  // Initialize.
  if (!pythiaColl.init()) return 1;

  // Book histograms for momentum check.
  Hist hAremn("remnant A number",       20, -0.5, 19.5);
  Hist hBef( "momentum fraction full final state",       100, 0.8, 1.2);
  Hist hRemn("momentum fraction target remnant",         100, 0.0, 1.0);
  Hist hAft( "momentum fraction without target remnant", 100, 0.8, 1.2);
  Hist hCorr("momentum fraction after correction",       100, 0.99, 1.01);

  // Loop over events.
  int nElastic = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Optionally vary incoming energy from event to event.
    if (varyE) {
      pBeam = pBeamMn + pythiaColl.rndm.flat() * (pBeamMx - pBeamMn);
      eBeam = sqrt(pow2(pBeam) + pow2(mBeam));
      p4Beam = Vec4( 0., 0., pBeam, eBeam);
      mCM   = ( p4Beam + p4Targ).mCalc();
    }

    // Insert primary particle in cleared main event record.
    //eventColl.clear();
    //eventColl.append(90,  -11, 0, 0, 1, 1, 0, 0, p4Beam + p4Targ, mCM);
    // Next line not needed (?).
    //eventColl.append(idBeam, 12, 0, 0, 0, 0, 0, 0, p4Beam, mBeam);

    // Set collision kinematics.
    pythiaColl.setBeamIDs(idBeam, idTarg);
    pythiaColl.setKinematics( eBeam, mTarg);

    // Simulate interaction.
    pythiaColl.next();
    if (iEvent == 0) eventColl.list();

    // Check momentum and catch target remnant.
    Vec4 p4Final;
    int iTargRemn  = 0;
    for (int i = 0; i < eventColl.size(); ++i) if (eventColl[i].isFinal()) {
      p4Final += eventColl[i].p();
      if (eventColl[i].id() > 1000000000) iTargRemn = i; 
    }

    // A of remnant. Catch elastic scatterings.
    int Aremn = (iTargRemn > 0) ? (eventColl[iTargRemn].id() / 10) % 1000 : 0;
    hAremn.fill( Aremn);
    if (Aremn == Atarg) {
      ++nElastic;

    // Check how off the (inelastic) final state is. Remove target remnant.
    } else {
      hBef.fill( p4Final.pAbs() / pBeam);
      if (iTargRemn > 0) {
        hRemn.fill( eventColl[iTargRemn].pAbs() / pBeam);
        eventColl[iTargRemn].statusNeg();
        p4Final -= eventColl[iTargRemn].p();
      }
      hAft.fill( p4Final.pAbs() / pBeam);

      // Do boost to bring final state to correct three-momentum.
      Vec4 p4Want( 0., 0., pBeam, sqrt( pow2(pBeam) + p4Final.m2Calc() ) );
      RotBstMatrix Mtransfer;
      Mtransfer.bst( p4Final, p4Want);
      eventColl.rotbst( Mtransfer);
    }
    
    // Check outcome.
    if (iEvent == 0) eventColl.list();
    Vec4 p4Corr;
    for (int i = 0; i < eventColl.size(); ++i) if (eventColl[i].isFinal()) 
      p4Corr += eventColl[i].p();
    hCorr.fill( p4Corr.pAbs() / pBeam);

  // End of event loop.  
  }

  // Print statistics, mainly for errors. Also histograms.
  pythiaColl.stat();
  cout << "\n number of elastic events = " << nElastic << endl;
  cout << hAremn << hBef << hRemn << hAft << hCorr;

}
