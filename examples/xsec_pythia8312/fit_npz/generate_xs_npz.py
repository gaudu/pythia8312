import glob # search for files using patterns
import numpy as np
from particle import Particle # handling pdg id database

# define the structure of the data to be loaded with specific fields and types
dtype = np.dtype([
    ("plab", np.float32),
    ("sig_tot", np.float32),
    ("err_tot", np.float32),
    ("sig_el", np.float32),
    ("err_el", np.float32),
    ("sig_inel", np.float32),
    ("err_inel", np.float32)
])

# path for the input data files
path_xs = "/home/gaudu/PYTHIA/pythia8312/examples/xsec_pythia8312/fit_dat/main1011_{projectile}_{target}_pNow1e{logPlab}_*.dat"

# map of projectile particle names and pdg ids
projectiles = {
    "p": 2212,
    "ap": -2212,
    "n": 2112,
    "an": -2112,
    "pip": 211,
    "pim": -211,
    "pi0": 111,
    "kp": 321,
    "km": -321,
    "k0": 311,
    "k0s": 310,
    "k0l": 130,
    "lambda0": 3122,
    "sigmap": 3222,
    "sigmam": 3112,
    "xi0": 3322,
    "omegam": 3334,
    "dp": 411,
    "d0": 421,
    "b0": 511,
    "bp": 521,
    "4He": 1000020040,
    "14N": 1000070140,
    "56Fe": 1000260560
}
# map of target particle names and pdg ids
targets = {
    "1H": 2212,
    "12C": 1000060120,
    "14N": 1000070140,
    "16O": 1000080160,
    "40Ar": 1000180400
}

# reverse mapping from pdg id to particle names
proj_pdg_to_str = dict((v, k) for k, v in projectiles.items())
target_pdg_to_str = dict((v, k) for k, v in targets.items())

# generate a list of discrete logarithmic lab momentum values
momenta = ["{:02}".format(i) for i in range(2, 12)]

# format the file path based on projectile, target, logarithmic lab momentum
def format_path(proj, targ, logP):
    return path_xs.format(projectile=proj_pdg_to_str[proj], target=target_pdg_to_str[targ], logPlab=str(logP))

particle_data = {}
# loop over each projectile
for proj_pdg in projectiles.values():
    projectile = Particle.from_pdgid(proj_pdg)
    # loop over each target
    for target_pdg in targets.values():
        table = np.zeros(len(momenta), dtype=dtype)
        # loop over each momentum value
        for i, logPlab in enumerate(momenta):
            path = format_path(proj_pdg, target_pdg, logPlab) # format input file path
            gl = glob.glob(path) # search for files matching the path
            if not gl:
                print(f"No file found for {proj_pdg}, {target_pdg}, {logPlab}. Skipping.")
                continue
            # open and read the first matching file 
            with open(next(iter(gl)), "rt") as f:
                data = dict()
                 # loop over each line
                for l in f:
                    key, val = l.split()
                    try:
                        data[key] = int(val)
                    except ValueError:
                        data[key] = float(val)
            # validate the projectile and target ids from the file with filename
            assert(data["targ_id"] == target_pdg)
            assert(data["proj_id"] == proj_pdg)
            # fill the data table for valid keys
            for key, val in data.items():
                if key in dtype.fields:
                    table[i][key] = float(val)
        # store the table in particle_data 
        particle_data[(proj_pdg), (target_pdg)] = table
        # save to npz file
        d = {key: particle_data[(proj_pdg), (target_pdg)][key] for key in dtype.fields}
        np.savez(f"xs_{proj_pdg}_{target_pdg}.npz", **d)