// Authors: Chloé Gaudu <gaudu@uni-wuppertal.de>
// Keywords: electron-positron annihilation; photon exchange; unpolarized tau decay channels; hepmc

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"

using namespace Pythia8;

int main() {

    Pythia8ToHepMC toHepMC("main1021.hepmc");
    Pythia pythia;

    // Electron-positron annihilation at 91.2 GeV
    pythia.readString("Beams:idA = 11"); 
    pythia.readString("Beams:idB = -11"); 
    pythia.readString("Beams:frameType = 1");
    pythia.readString("Beams:eCM = 91.2");

    /* Electroweak Processes

    WeakZ0:gmZmode
    ┌────────┬───────────────────────────────────────────────────────────────────────┐
    │ Option │ Description                                                           │
    ├────────┼───────────────────────────────────────────────────────────────────────┤
    │   0    │ Full γ^* or Z^0 structure, with interference included                 │
    │   1    │ Only pure γ^* contribution                                            │
    │   2    │ Only pure Z^0 contribution                                            │
    └────────┴───────────────────────────────────────────────────────────────────────┘
    Boson Exchange
    ┌───────────────────────────────────────────┬────────────────────────────────────┐
    │ Parameter                                 │ Description                        │
    ├───────────────────────────────────────────┼────────────────────────────────────┤
    │ WeakBosonExchange:ff2ff(t:gmZ)            │ f f' → f f' via γ^* or Z^0         │
    │ Code 211                                  │ t-channel exchange, with full      │
    │                                           │ interference between γ^* and Z^0   │
    └───────────────────────────────────────────┴────────────────────────────────────┘
    Single Boson
    ┌───────────────────────────────────────────┬────────────────────────────────────┐
    │ Parameter                                 │ Description                        │
    ├───────────────────────────────────────────┼────────────────────────────────────┤
    │ WeakSingleBoson:ffbar2gmZ                 │ f fbar → γ^* or Z^0, with full     │
    │ Code 221                                  │ interference between γ^* and Z^0   │
    ├───────────────────────────────────────────┼────────────────────────────────────┤
    │ WeakSingleBoson:ffbar2ffbar(s:gm)         │ f fbar → γ^* → f' fbar',           │
    │ Code 223                                  │ written as a 2 → 2 process,        │
    │                                           │ hardcoded for final states of      │
    │                                           │ 5 quark flavours or 3 lepton ones  │
    ├───────────────────────────────────────────┼────────────────────────────────────┤
    │ WeakSingleBoson:ffbar2ffbar(s:gmZ)        │ f fbar → γ^* or Z^0 → f'           │
    │ Code 224                                  │ fbar', written as a 2 → 2 process, │
    │                                           │ final-state flavour selection is   │
    │                                           │ based on Z^0 allowed decay modes   │
    └───────────────────────────────────────────┴────────────────────────────────────┘ */
    pythia.readString("WeakZ0:gmZmode = 0");
    //pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
    //pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("WeakSingleBoson:ffbar2ffbar(s:gm) = on");
    //pythia.readString("WeakSingleBoson:ffbar2ffbar(s:gmZ) = on");    

    // Decay channels for Z
    pythia.readString("23:onMode = off"); // Turn off Z decays
    pythia.readString("23:oneChannel = 15 -15"); // Z → τ+ τ−

    // Decay channels for τ
    //pythia.readString("15:onMode = on"); // Turn on τ decays
    pythia.readString("15:oneChannel = 1 0.17835 0 -16 11 -12");     // τ → ν_τ e- ν_e
    pythia.readString("15:addChannel = 1 0.17364 0 -16 13 -14");     // τ → ν_τ μ- ν_μ
    pythia.readString("15:addChannel = 1 0.11070 0 -16 211");        // τ → ν_τ π-
    pythia.readString("15:addChannel = 1 0.25490 0 -16 111 111");    // τ → ν_τ π0 π0 (CLEO model)
    pythia.readString("15:addChannel = 1 0.18600 0 -16 211 111 111");// τ → ν_τ π- π0 π0 (Kühn-Santamaria model)

    // Initialize
    pythia.init();

    // Event generation loop
    for (int iEvent = 0; iEvent < 1000000; ++iEvent) {
        if (!pythia.next()) continue;

        Pythia8::Event& event = pythia.event;

        // Write events with unpolarized taus to HepMC file
        bool keepEvent = false;
        for (int i = 0; i < event.size(); ++i) {
            if (abs(event[i].id()) == 15) { // Check for τ or τ-
                int motherId = event[event[i].mother1()].id();
                if (motherId == 22) { // Check for γ mother
                    keepEvent = true;
                    break;
                }
            }
        }

        if (keepEvent) {
            toHepMC.writeNextEvent(pythia);
        }

        //event.list();
    }

    // Print statistics to terminal
    pythia.stat();

    return 0;
}
