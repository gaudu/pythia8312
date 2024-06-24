#include <string_view>

#include <Pythia8/Pythia.h>
#include <corsika/corsika.hpp>

class Interaction {

public:
Interaction(std::string_view dataPath) {

    // is dataPath supposed to be equivalent to CORSIKA_DATA_DIR as defined in corsika.hpp?
    // cross-section tables are in CORSIKA_DATA_DIR/Pythia/pythia8311_xsec/.

    std::string mpi_path = CORSIKA_DATA_DIR+"/Pythia/pythia8312.mpi";
    std::string sasd_mpi_path = CORSIKA_DATA_DIR+"/Pythia/pythia8312.sasd.mpi";
    std::string sigfit_path = CORSIKA_DATA_DIR+"/Pythia/pythia8312.sigfit";

    int nPts = 11;

    vector<int> idA_list = {2212, -2212, 2112, -2112, 111, 211, -211, 311, 321, -321, 130, 310, 3122, 3212, 3222, 3112, 3322};
    vector<int> idB_list = {2212, 1000060120, 1000070140, 1000080160, 1000180400};
    double idA = idA_list[0]; // proton
    double idB = idB_list[3]; // nitrogen-14
    pythia.readString("Beams:allowIDASwitch = on");
    pythia.settings.mode("Beams:idA", idA);
    pythia.settings.mode("Beams:idB", idB);
    
    double pA_max = 1e12;
    double pA_min = 1e2; 
    double mA_max = pythia.particleData.m0(3322); // Xi0 heaviest
    double mA_min = pythia.particleData.m0(111); // pi0 lightest
    double mB = pythia.particleData.m0(2212); // Angantyr does a weird mass per nucleus operation, therefore mB is always the proton mass (even for eCM calculations I was told)
    double eA_max = sqrt(pow2(pA_max) + pow2(mA_max));  
    double eA_min = sqrt(pow2(pA_min) + pow2(mA_min)); 
    double eCM_max = sqrt(pow2(mA_max) + pow2(mB) + 2.* eA_max * mB);
    double eCM_min = sqrt(pow2(mA_min) + pow2(mB) + 2.* eA_min * mB);
    pythia8_.readString("Beams:allowVariableEnergy = on");
    pythia8_.settings.parm("HeavyIon:varECMMax", eCM_max);
    pythia8_.settings.parm("HeavyIon:varECMMin", eCM_min);
    pythia8_.settings.parm("HeavyIon:varECMSigFitNPts", nPts);  
    pythia8_.readString("Beams:frameType = 1");
    pythia8_.settings.parm("Beams:eCM", eCM_max);

    pythia8_.readString("SoftQCD:all = on"); // needed for regular Pythia pp run, but unused by Angantyr

    pythia8_.readString("MultipartonInteractions:reuseInit = 3");
    pythia8_.settings.parm("MultipartonInteractions:initFile", mpi_path);
    pythia8_.readString("HeavyIon:SasdMpiReuseInit = 3");
    pythia8_.settings.parm("HeavyIon:SasdMpiInitFile", sasd_mpi_path );
    pythia8_.readString("HeavyIon:SigFitReuseInit = 3");
    pythia8_.settings.parm("HeavyIon:SigFitInitFile", sigfit_path);
}


void generateEvent(double sqrtS, int projectilePDG, int targetPDG) {
    
    int nEvent = 1;

    double eCM_switch = sqrtS;

    vector<int> idA_list = {2212, -2212, 2112, -2112, 111, 211, -211, 311, 321, -321, 130, 310, 3122, 3212, 3222, 3112, 3322};
    auto itA = std::find(idA_list.begin(), idA_list.end(), projectilePDG);
    if (itA != idA_list.end()) {
        idA_switch = projectilePDG;
    } else {
        throw std::runtime_error("Error: This projectile is not (yet) in the pre-computed cross-section tables");
    }

    vector<int> idB_list = {2212, 1000060120, 1000070140, 1000080160, 1000180400};
    auto itB = std::find(idB_list.begin(), idB_list.end(), targetPDG);
    if (itB != idB_list.end()) {
        idB_switch = targetPDG;
    } else {
        throw std::runtime_error("This target is not (yet) in the pre-computed cross-section tables");
    }

    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    
        pythia8_.setBeamIDs(idA_switch, idB_switch);
        pythia8_.setKinematics(eCM_switch);

        pythia8_.next();
    }
    
    pythia8_.event.clear();
}

private:
  Pythia8::Pythia pythia8_;
};