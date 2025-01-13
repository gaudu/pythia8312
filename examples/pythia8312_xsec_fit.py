import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import linregress

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
            data_dict[key].append((momentum, parsed_data, filename))
    
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

def log_log_fit(data_dict, folder_path, projectiles_to_fit):
    fit_data = {}
    for (projectile, target), data_list in data_dict.items():
        if projectile not in projectiles_to_fit:
            continue

        # extract momenta and inelastic cross-sections
        momenta = np.array([entry[0] for entry in data_list])  
        sig_inel = np.array([entry[1]['sig_inel'] for entry in data_list])  

        # perform log-log linear regression
        log_momenta = np.log10(momenta) 
        log_sig_inel = np.log10(sig_inel)
        slope, intercept, _, _, _ = linregress(log_momenta, log_sig_inel)

        # generate fitted values
        sig_inel_fit = 10 ** (slope * log_momenta + intercept)
        
        output_path = "xsec_pythia8312/fit_data_2025"
        os.makedirs(output_path, exist_ok=True)

        # save fitted values in new files
        for i, entry in enumerate(data_list):
            momentum, data, filename = entry  # Unpack all three elements
            original_filepath = os.path.join(folder_path, filename)  # Use the correct filename
            fit_filename = filename.replace(".dat", "_fit.dat")  # Generate the fit filename
            fit_filepath = os.path.join(output_path, fit_filename)

            # write fitted file
            with open(original_filepath, 'r') as original_file, open(fit_filepath, 'w') as fit_file:
                for line in original_file:
                    if line.startswith("sig_inel"):
                        fit_file.write(f"sig_inel\t{sig_inel_fit[i]:.6f}\n")
                    else:
                        fit_file.write(line)

        # store fit data for plotting
        fit_data[(projectile, target)] = {
            'momenta': momenta,
            'sig_inel_fit': sig_inel_fit,
            'fit_values': sig_inel_fit,  
            'slope': slope,
            'intercept': intercept,
        }

    return fit_data

def plot_data_to_pdf(data_dict, fit_data, output_pdf):
    with PdfPages(output_pdf) as pdf:
        for (projectile, target), data_list in data_dict.items():
            momenta = np.array([entry[0] for entry in data_list])
            sig_inel = np.array([entry[1]['sig_inel'] for entry in data_list])
            err_inel = np.array([entry[1]['err_inel'] for entry in data_list])

            plt.figure(figsize=(10, 6))
            plt.errorbar(momenta, sig_inel, yerr=err_inel, label="σ_inel", fmt='^-', capsize=3)

            # add fit to the plot if available
            if (projectile, target) in fit_data:
                fit_momenta = fit_data[(projectile, target)]['momenta']
                fit_values = fit_data[(projectile, target)]['fit_values']
                plt.plot(fit_momenta, fit_values, label="σ_inel fit", linestyle='--', color='red')

            plt.xscale('log')
            plt.xlabel("Momentum (GeV/c)", fontsize=12)
            plt.ylabel("Cross-section (mb)", fontsize=12)
            plt.title(f"Cross-sections for {projectile} on {target}", fontsize=14)
            plt.legend()
            plt.grid(True, which="both", linestyle="--", linewidth=0.5)
            plt.tight_layout()

            pdf.savefig()
            plt.close()
            print(f"Added plot for {projectile} on {target} to PDF.")

def main():
    folder_path = "xsec_pythia8312"
    output_folder = "xsec_pythia8312/fit_dat"
    os.makedirs(output_folder, exist_ok=True)
    output_pdf = "xsec_pythia8312_fit.pdf"
    
    projectiles_to_fit = ['p', 'ap', 'n', 'an', 'pip', 'pim', 'pi0', 'kp', 'km', 'k0', 'k0s', 'k0l', 'lambda0', 'sigmap', 'sigmam', 'xi0', 'omegam', 'dp', 'd0', 'b0', 'bp', '4He', '14N', '56Fe']

    dict_dat = gather_data(folder_path)
    organized_dat = organize_data(dict_dat)
    fit = log_log_fit(organized_dat, folder_path, projectiles_to_fit)
    plot_data_to_pdf(organized_dat, fit, output_pdf)

    print(f"All plots have been saved to {output_pdf}.")

if __name__ == "__main__":
    main()
