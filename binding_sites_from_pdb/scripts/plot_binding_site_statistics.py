#!/usr/bin/env python3
'''Plot the ligand binding site statistics.
'''

import json

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt


def print_binding_sites_with_given_properties(bs_infos, min_lig_size, max_lig_size, min_site_size, max_site_size):
    '''Print the binding sites that have the given properties'''
   
    for bs in bs_infos:
        lig_size = bs['binding_site_info']['ligand_n_heavy_atoms']
        site_size = len(bs['binding_site_info']['binding_site_energies'])

        if lig_size > max_lig_size or lig_size < min_lig_size:
            continue

        if site_size > max_site_size or site_size < min_site_size:
            continue

        print(bs['pdb_file'], lig_size, site_size)

def print_most_common_ligand_types(bs_infos):
    '''Print most common ligand types.'''
    binding_sites_by_lig_names = {}

    for bs in bs_infos:
        ligand_name = (bs['binding_site_info']['ligand_name'].strip(), bs['binding_site_info']['ligand_n_heavy_atoms'])
        pdb_file = bs['pdb_file']

        if ligand_name in binding_sites_by_lig_names:
            binding_sites_by_lig_names[ligand_name].append(pdb_file)
        else:
            binding_sites_by_lig_names[ligand_name] = [pdb_file]
    
    unique_ligands = binding_sites_by_lig_names.keys()
    unique_ligands_sorted = sorted(unique_ligands, key=lambda x : len(binding_sites_by_lig_names[x]), reverse=True)

    print('The number of unique ligands is {0}'.format(len(unique_ligands)))
    print('The most common ligands are (lig_name, lig_size, N_sites, example_pdb):')
   
    for ligand_name in unique_ligands_sorted[:50]:
        print(ligand_name[0], ligand_name[1], len(binding_sites_by_lig_names[ligand_name]), binding_sites_by_lig_names[ligand_name][0] )

def get_n_ligand_heavy_atoms(bs_infos):
    '''Get a list of numbers of ligand heavy atoms.'''
    ligand_sizes = []
    
    for bs in bs_infos:
        ligand_sizes.append(bs['binding_site_info']['ligand_n_heavy_atoms'])

    return ligand_sizes

def get_binding_site_sizes(bs_infos):
    '''Get the number of residues in the ligand binding sizes.'''
    binding_site_sizes = []
    
    for bs in bs_infos:
        binding_site_sizes.append(len(bs['binding_site_info']['binding_site_energies']))

    return binding_site_sizes

def print_frequencies_of_n_lig_heavy_atoms_and_binding_site_sizes(bs_infos):
    ligand_sizes = get_n_ligand_heavy_atoms(bs_infos)
    binding_site_sizes = get_binding_site_sizes(bs_infos)

    N_total = len(ligand_sizes)

    count_lig_sizes = {}
    for n in ligand_sizes:
        if n in count_lig_sizes.keys():
            count_lig_sizes[n] += 1
        else:
            count_lig_sizes[n] = 1

    acculumated_count = 0
    for k in sorted(count_lig_sizes.keys()):
        acculumated_count += count_lig_sizes[k]
        print('N_heavy_atom = {0}, count = {1}, frequency = {2:.2f}%, acculumated_count = {3}, acculumated_freq = {4:.2f}%'.format(
            k, count_lig_sizes[k], count_lig_sizes[k] / N_total * 100, acculumated_count, acculumated_count / N_total * 100))

    count_site_sizes = {}
    for n in binding_site_sizes:
        if n in count_site_sizes.keys():
            count_site_sizes[n] += 1
        else:
            count_site_sizes[n] = 1

    acculumated_count = 0
    for k in sorted(count_site_sizes.keys()):
        acculumated_count += count_site_sizes[k]
        print('binding_site_size = {0}, count = {1}, frequency = {2:.2f}%, acculumated_count = {3}, acculumated_freq = {4:.2f}%'.format(
            k, count_site_sizes[k], count_site_sizes[k] / N_total * 100, acculumated_count, acculumated_count / N_total * 100))



def get_binding_site_energies(bs_infos, energy_type):
    '''Get a list of energies of a given type'''
    energies = []

    for bs in bs_infos:
        for r in bs['binding_site_info']['binding_site_energies']:
            e_list = r[3]
            energies.append(e_list[energy_type])

    print('{0} out of {1} structures have 0 {2} energy.'.format(len([e for e in energies if e==0]), len(energies), energy_type))
    return energies


def print_highest_and_lowest_energy_residues(bs_infos, energy_type):
    '''Print the residues that have highest and lowest energies
    for a given type.
    '''
    e_max = -float('inf')
    e_min = float('inf')

    high_pdb = ''
    high_ligand = ''
    high_res_id = ''
    low_pdb = ''
    low_ligand = ''
    low_res_id = ''
    
    for bs in bs_infos:
        for r in bs['binding_site_info']['binding_site_energies']:
            e_list = r[3]
            e = e_list[energy_type]

            if e > e_max:
                e_max = e
                high_pdb = bs['info_file']
                high_ligand = bs['binding_site_info']['ligand_number'] 
                high_res_id = r[1]

            if e < e_min:
                e_min = e
                low_pdb = bs['info_file']
                low_ligand = bs['binding_site_info']['ligand_number'] 
                low_res_id = r[1]

    print('The highest {0} energy is {1}. Structure: {2}. Ligand: {3}. Residue: {4}.'.format(energy_type, e_max, high_pdb, high_ligand, high_res_id))
    print('The lowest {0} energy is {1}. Structure: {2}. Ligand: {3}. Residue: {4}.'.format(energy_type, e_min, low_pdb, low_ligand, low_res_id))



def make_histogram(data, x_label, y_label='Count', upper_bound=None, lower_bound=None, num_bins=50, output_file=None, y_max=None):
    '''Make a histogram from a list of data.'''
    if upper_bound is None:
        upper_bound = max(data)
    if lower_bound is None:
        lower_bound = min(data)

    bin_width = (upper_bound - lower_bound) / (num_bins - 1)

    fig, ax1 = plt.subplots()

    bins = [lower_bound + (i - 0.5) * bin_width for i in range(num_bins + 1)]
    hist, bin_edges = np.histogram(data, bins=bins)

    ax1.bar(np.array(bin_edges[0:-1]) + 0.5 * bin_width, hist, width=bin_width, alpha=1)

    if not y_max is None:
        plt.ylim(top=y_max)

    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)

    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)

    plt.clf()
    
def plot_num_aa_types(bs_infos, output_file=None, ligand_sizes=None):
    # The Amino acid scale: Amino acid composition (%) in the UniProtKB/Swiss-Prot data bank.
    # Reference: Release notes for UniProtKB/Swiss-Prot release 2013_04 - April 2013.
    # https://web.expasy.org/protscale/pscale/A.A.Swiss-Prot.html

    backgroud_aa_freq = {
            'ALA':8.25, 'ARG':5.53, 'ASN':4.06, 'ASP':5.45, 'CYS':1.37, 'GLN':3.93, 'GLU':6.75, 'GLY':7.07, 'HIS':2.27, 'ILE':5.96, 'LEU':9.66,
            'LYS':5.84, 'MET':2.42, 'PHE':3.86, 'PRO':4.70, 'SER':6.56, 'THR':5.34, 'TRP':1.08, 'TYR':2.92, 'VAL':6.87} 

    # Get the AA composition in binding sites
    
    num_aas = {}
    
    for bs in bs_infos:
        lig_size = bs['binding_site_info']['ligand_n_heavy_atoms']
        
        if (ligand_sizes is None) or lig_size in ligand_sizes:
        
            for r in bs['binding_site_info']['binding_site_energies']:
                aa_type = r[2]


                if aa_type in num_aas:
                    num_aas[aa_type] += 1
                else:
                    num_aas[aa_type] = 0

    N_total = sum(num_aas[k] for k in num_aas.keys())
    
    enrichment_ratios = {}
    for k in num_aas.keys():
        enrichment_ratios[k] = 100 * num_aas[k] / N_total / backgroud_aa_freq[k] 

    enrichment_ratios_list = [(k, enrichment_ratios[k]) for k in enrichment_ratios.keys()]
    enrichment_ratios_list_sorted = sorted(enrichment_ratios_list, key=lambda x:x[1])
    print(enrichment_ratios_list_sorted)
    N_aa_types = len(enrichment_ratios_list_sorted)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx() 

    rec1 = ax1.bar([i - 0.15 for i in range(N_aa_types)], [x[1] for x in enrichment_ratios_list_sorted], 0.3, color='blue')
    rec2 = ax2.bar([i + 0.15 for i in range(N_aa_types)], [100 * num_aas[x[0]] / N_total for x in enrichment_ratios_list_sorted], 0.3, color='red')
   
    ax1.set_ylabel('Enrichment ratio', fontsize=15)
    ax2.set_ylabel('Frequency (%)', fontsize=15)
    ax1.set_xticks([i for i in range(N_aa_types)])
    ax1.set_xticklabels([x[0] for x in enrichment_ratios_list_sorted], rotation=-40)

    plt.legend([rec1, rec2], ['enrichment ratio', 'AA frequency'], prop={'size': 15})

    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)

    plt.clf()

def plot_lig_sizes_vs_site_sizes(ligand_sizes, binding_site_sizes, output_file=None):
   
    # Print the linear regression

    print('Linear regression for ligand_sizes and binding_site_sizes:')
    print(scipy.stats.linregress(ligand_sizes, binding_site_sizes))

    # Get the matrix for the heat map

    m = []
    max_lig_size = 70
    max_site_size = 35

    for i in range(max_lig_size):
        m.append([0] * (max_site_size))

    for i in range(len(ligand_sizes)):
        lig_size = ligand_sizes[i]
        site_size = binding_site_sizes[i]

        if lig_size <= max_lig_size and 0 < site_size <= max_site_size:
            m[lig_size - 1][site_size - 1] += 1

    fig, ax = plt.subplots()
    cax = ax.imshow(np.transpose(m), origin='lower', interpolation='nearest', cmap='cividis', vmin=0, vmax=30, 
            extent=(0.5, max_lig_size + 0.5, 0.5, max_site_size + 0.5), aspect='auto')

    #plt.scatter(ligand_sizes, binding_site_sizes)

    fig.colorbar(cax)
    
    ax.tick_params(axis='both', labelsize=15)
    plt.xlabel('N ligand heavy atoms', fontsize=15)
    plt.ylabel('N binding site residues', fontsize=15)
    plt.subplots_adjust(bottom=0.15)

    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)

    plt.clf()


if __name__ == '__main__':

    input_file = './binding_site_database_summary.json'
    #input_file = './obsolete/2020_08_07/binding_site_database_summary.json'
    output_path = './plots'

    with open(input_file, 'r') as f:
        bs_infos = json.load(f)

    print('The total number of binding sites is {0}.'.format(len(bs_infos)))

    print_binding_sites_with_given_properties(bs_infos, 1, 1, 40, 1000)

    # Print most common ligands

    print_most_common_ligand_types(bs_infos)

    # Print frequencies of ligand sizes and binding site sizes

    print_frequencies_of_n_lig_heavy_atoms_and_binding_site_sizes(bs_infos)

    # Plot the ligand sizes vs binding site sizes

    ligand_sizes = get_n_ligand_heavy_atoms(bs_infos)
    binding_site_sizes = get_binding_site_sizes(bs_infos)
    print(min(ligand_sizes), max(ligand_sizes))
    print(min(binding_site_sizes), max(binding_site_sizes))
    print('There are {0} binding sites with at least 3 protein residues'.format(len([x for x in binding_site_sizes if x > 2])))
   
    plot_lig_sizes_vs_site_sizes(ligand_sizes, binding_site_sizes, output_file='./plots/lig_sizes_vs_site_sizes.svg')
    
    make_histogram(ligand_sizes, 'number of ligand heavy atoms', upper_bound=100, lower_bound=1, num_bins=100, output_file='./plots/num_ligand_heavy_atoms.png')
    make_histogram(binding_site_sizes, 'Number of residues in ligand binding site', upper_bound=50, lower_bound=1, num_bins=50, output_file='./plots/binding_site_sizes.png')

    # Plot the histogram energies

    atr_energies = get_binding_site_energies(bs_infos, 'fa_atr')
    make_histogram(atr_energies, 'fa_atr (REU)', upper_bound=0, lower_bound=-15, output_file='./plots/hist_atr_energies.png')
    
    rep_energies = get_binding_site_energies(bs_infos, 'fa_rep')
    make_histogram(rep_energies, 'fa_rep (REU)', upper_bound=5, lower_bound=0, output_file='./plots/hist_rep_energies.png')
    
    sol_energies = get_binding_site_energies(bs_infos, 'fa_sol')
    make_histogram(sol_energies, 'fa_sol (REU)', upper_bound=20, lower_bound=-5, output_file='./plots/hist_sol_energies.png')
    
    elec_energies = get_binding_site_energies(bs_infos, 'fa_elec')
    make_histogram(elec_energies, 'fa_elec (REU)', upper_bound=5, lower_bound=-30, output_file='./plots/hist_elec_energies.png')
    
    lk_ball_energies = get_binding_site_energies(bs_infos, 'lk_ball_wtd')
    make_histogram(lk_ball_energies, 'lk_ball_wtd (REU)', upper_bound=2, lower_bound=-2, output_file='./plots/hist_lk_ball_energies.png')
    
    hbond_bb_sc_energies = get_binding_site_energies(bs_infos, 'hbond_bb_sc')
    make_histogram(hbond_bb_sc_energies, 'hbond_bb_sc (REU)', upper_bound=0, lower_bound=-4, output_file='./plots/hist_hbond_bb_sc_energies.png', y_max=2000)
    
    hbond_sc_energies = get_binding_site_energies(bs_infos, 'hbond_sc')
    make_histogram(hbond_sc_energies, 'hbond_sc (REU)', upper_bound=0, lower_bound=-4, output_file='./plots/hist_hbond_sc_energies.png', y_max=2000)
    
    total_energies = get_binding_site_energies(bs_infos, 'weighted_total')
    make_histogram(total_energies, 'weighted total 2 body ref2015 energy (REU)', upper_bound=20, lower_bound=-40, output_file='./plots/hist_total_energies.png')
    
    # Print the residues that give the highest and lowest energies
    
    print_highest_and_lowest_energy_residues(bs_infos, 'fa_atr')
    print_highest_and_lowest_energy_residues(bs_infos, 'fa_rep')
    print_highest_and_lowest_energy_residues(bs_infos, 'fa_sol')
    print_highest_and_lowest_energy_residues(bs_infos, 'fa_elec')
    print_highest_and_lowest_energy_residues(bs_infos, 'lk_ball_wtd')
    print_highest_and_lowest_energy_residues(bs_infos, 'hbond_bb_sc')
    print_highest_and_lowest_energy_residues(bs_infos, 'hbond_sc')
    print_highest_and_lowest_energy_residues(bs_infos, 'weighted_total')

    # Plot the AA composition of binding sites

    plot_num_aa_types(bs_infos, output_file='./plots/aa_composition.svg')
   
    plot_num_aa_types(bs_infos, ligand_sizes=[1], output_file='./plots/aa_composition_lig_size_1.svg')
    plot_num_aa_types(bs_infos, ligand_sizes=list(range(2, 101)), output_file='./plots/aa_composition_lig_size_2_to_100.svg')
