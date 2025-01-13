import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_dat_file(filepath):
    data = {}
    with open(filepath, 'r') as file:
        for line in file:
            key, value = line.split()
            data[key] = float(value)
    return data

def extract_info_from_filename(filename):
    match = re.match(r"main1011_(\w+)_(\w+)_pNow([\d\.e]+)_nEvents", filename)
    if match:
        projectile, target, momentum = match.groups()
        return projectile, target, float(momentum)
    else:
        raise ValueError(f"Filename '{filename}' does not match the expected format.")

def gather_data(folder_path):
    data_dict = {}
    required_keys = ['sig_tot', 'err_tot', 'sig_el', 'err_el', 'sig_inel', 'err_inel']
    
    for filename in os.listdir(folder_path):
        if filename.endswith(".dat"):
            filepath = os.path.join(folder_path, filename)
            projectile, target, momentum = extract_info_from_filename(filename)
            parsed_data = parse_dat_file(filepath)

            if not all(key in parsed_data for key in required_keys):
                print(f"Skipping {filename}: Missing required keys.")
                continue

            key = (projectile, target)
            if key not in data_dict:
                data_dict[key] = []
            data_dict[key].append((momentum, parsed_data))
    
    for key in data_dict:
        data_dict[key].sort(key=lambda x: x[0])
    return data_dict

def organize_data(data_dict):
    target_order = ['1H', '12C', '14N', '16O', '40Ar']
    projectile_order = ['p', 'ap', 'n', 'an', 'pip', 'pim', 'pi0', 'kp', 'km', 'k0', 'k0s', 'k0l', 'lambda0', 'sigmap', 'sigmam', 'xi0', 'omegam', 'dp', 'd0', 'b0', 'bp', '4He', '14N', '56Fe']
    
    sorted_keys = sorted(
        data_dict.keys(),
        key=lambda x: (
            projectile_order.index(x[0]) if x[0] in projectile_order else float('inf'),
            target_order.index(x[1]) if x[1] in target_order else float('inf')
        )
    )

    organized_data = {key: data_dict[key] for key in sorted_keys}
    return organized_data

def plot_data_to_pdf(data_dict, output_pdf):
    """
    Plot data from data_dict dictionary into a single PDF.

    Parameters:
    - data_dict: dict
        Dictionary with keys (projectile, target) and values as lists of (momentum, data) tuples.
    - output_pdf: str
        Path to the output PDF file.
    """
    with PdfPages(output_pdf) as pdf:
        for (projectile, target), data_list in data_dict.items():
            momenta = []
            sig_inel = []
            err_inel = []

            for entry in data_list:
                momentum, data = entry
                momenta.append(momentum)
                sig_inel.append(data['sig_inel'])
                err_inel.append(data['err_inel'])
                
            plt.rcParams['text.usetex'] = True
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['font.serif'] = ['Computer Modern Roman']

            plt.figure(figsize=(5, 3))
            plt.errorbar(momenta, sig_inel, yerr=err_inel, fmt='-o', capsize=2)

            plt.xscale('log')
            plt.xlabel(r"P$_\mathrm{lab}$ / GeV", fontsize=13)
            plt.ylabel(r"$\sigma^\mathrm{inel}$ / mb", fontsize=13)
            plt.title(f"Cross-sections for {projectile} on {target}", fontsize=14)
            plt.grid(True, which="both", linestyle="--", linewidth=0.5)
            plt.tight_layout()

            pdf.savefig()
            plt.close()
            print(f"Added plot for {projectile} on {target} to PDF.")

def main():   
    folder_path = "xsec_pythia8312/fit_dat"
    output_pdf = "xsec_pythia8312.pdf"

    dict_dat = gather_data(folder_path)
    organized_dat = organize_data(dict_dat)
    plot_data_to_pdf(organized_dat, output_pdf)

    print(f"All plots have been saved to {output_pdf}.")

if __name__ == "__main__":
    main()
