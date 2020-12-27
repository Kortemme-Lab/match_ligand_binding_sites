#!/usr/bin/env python3

import json

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def binding_site_summary_to_dict(bs_summary):
    '''Convert a binding site summary into a dictionary.'''
    binding_site_summary_dict = {}

    for bs_info in bs_summary:
        pdb_file = bs_info['pdb_file']
    
        d1 = pdb_file.split('/')[-3]  
        d2 = pdb_file.split('/')[-2]
        d3 = pdb_file.split('/')[-1].split('.')[0].split('_')[-1]

        binding_site_summary_dict['/'.join([d1, d2, d3])] = bs_info

    return binding_site_summary_dict

def aggregated_results_to_dict(aggregated_results):
    '''Convert a list of aggregated_results into a dictionary.'''
    aggregated_results_dict = {}

    for i in aggregated_results:
        if 'match_info_file' in i:
            key = '/'.join(i['match_info_file'].split('/')[-4:-1])
            aggregated_results_dict[key] = i

    return aggregated_results_dict

def calculate_mean_sequence_distance_for_one_binding_site(bs_info):
    binding_site_residues = [x[0] for x in bs_info['binding_site_info']['binding_site_energies']]

    sequence_distances = []

    for i in range(len(binding_site_residues)):
        for j in range(i + 1, len(binding_site_residues)):
            sequence_distances.append(np.absolute(binding_site_residues[i] - binding_site_residues[j]))
   
    return np.mean(sequence_distances)

def get_binding_ligand_size_and_match_results(list_aggregated_results_dict):
    '''Get the binding site ligand sizes and the match results into a list.
    If multiple aggregated_results_dict are provided, only matches that are successful
    in all results are counted as successful.
    '''
    ligand_size_and_match_results = []

    for k in list_aggregated_results_dict[0]:
        if not ('ligand_n_heavy_atoms' in list_aggregated_results_dict[0][k]):
            continue

        ligand_size = list_aggregated_results_dict[0][k]['ligand_n_heavy_atoms']

        fast_match_succeed = True
        rosetta_match_succeed = True
       
        for aggregated_results_dict in list_aggregated_results_dict:

            if not (k in aggregated_results_dict):
                fast_match_succeed = False
                rosetta_match_succeed = False
                break

            match_info = aggregated_results_dict[k]
            
            if not ('num_fast_match_succeed' in match_info and match_info['num_fast_match_succeed'] > 0):
                fast_match_succeed = False

            if not ('num_rosetta_match_succeed' in match_info and match_info['num_rosetta_match_succeed'] > 0):
                rosetta_match_succeed = False

        ligand_size_and_match_results.append((ligand_size, fast_match_succeed, rosetta_match_succeed))

    return ligand_size_and_match_results


def get_binding_site_mean_sequence_distances_and_match_results(binding_site_summary_dict, list_aggregated_results_dict):
    '''Get the binding site mean sequence distances and the match results into a list.
    If multiple aggregated_results_dict are provided, only matches that are successful
    in all results are counted as successful.
    '''
    seq_dist_and_match_results = []

    for k in binding_site_summary_dict:

        bs_info = binding_site_summary_dict[k]
        mean_sequence_distance = calculate_mean_sequence_distance_for_one_binding_site(bs_info)

        fast_match_succeed = True
        rosetta_match_succeed = True
       
        for aggregated_results_dict in list_aggregated_results_dict:

            if not (k in aggregated_results_dict):
                fast_match_succeed = False
                rosetta_match_succeed = False
                break

            match_info = aggregated_results_dict[k]
            
            if not ('num_fast_match_succeed' in match_info and match_info['num_fast_match_succeed'] > 0):
                fast_match_succeed = False

            if not ('num_rosetta_match_succeed' in match_info and match_info['num_rosetta_match_succeed'] > 0):
                rosetta_match_succeed = False

        seq_dist_and_match_results.append((mean_sequence_distance, fast_match_succeed, rosetta_match_succeed))

    return seq_dist_and_match_results

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

def plot_mean_sequence_distance_distribution(seq_dist_and_match_results):
    seq_dists = [np.log(x[0]) for x in seq_dist_and_match_results]

    make_histogram(seq_dists, 'Log (mean sequence distance)') 

def assort_success_matches_into_bins(metric_and_match_results, bin_centers, bin_width, use_log=True):
    N_bins = len(bin_centers)
    
    N_rosetta_succeeds = [0] * N_bins
    N_fast_succeeds = [0] * N_bins
    N_total_results = [0] * N_bins

    for x in metric_and_match_results:
        if use_log:
            metric = np.log(x[0])
        else:
            metric = x[0]
        
        bin_id = None

        for i in range(len(bin_centers)):
            if bin_centers[i] - 0.5 * bin_width <= metric < bin_centers[i] + 0.5 * bin_width:
                bin_id = i

        if bin_id is None:
            continue

        N_total_results[bin_id] += 1

        if x[1]:
            N_fast_succeeds[bin_id] += 1
        if x[2]:
            N_rosetta_succeeds[bin_id] += 1

    return N_rosetta_succeeds, N_fast_succeeds, N_total_results


def plot_num_matches_vs_mean_binding_site_dist(all_seq_dist_and_match_results, output_file=None):
    bin_width = 0.5
    N_bins = 14
    bin_centers = [bin_width * (i + 0.5) for i in range(N_bins)]
    linear_bin_centers = [np.exp(x) for x in bin_centers]

    fig, ax1 = plt.subplots()
    
    for k in all_seq_dist_and_match_results:
        N_rosetta_succeeds, N_fast_succeeds, N_total_results = assort_success_matches_into_bins(all_seq_dist_and_match_results[k], bin_centers, bin_width)

        if k == '3res_native_rossmann':
            color = 'yellow' 
            rec1 = ax1.plot(linear_bin_centers, N_rosetta_succeeds, color=color, label='native Rossmann')
        elif k == '3res_2lv8':
            color = 'blue'
            rec2 = ax1.plot(linear_bin_centers, N_rosetta_succeeds, color=color, label='de novo Rossmann')
        elif k == '3res_native_ntf2':
            color = 'cyan'
            rec3 = ax1.plot(linear_bin_centers, N_rosetta_succeeds, color=color, label='native NTF2')
        elif k == '3res_5tpj':
            color = 'red'
            rec4 = ax1.plot(linear_bin_centers, N_rosetta_succeeds, color=color, label='de novo NTF2')

    ax1.set_xlim(1, 1000)
    #ax1.set_ylim(0, 100)
    ax1.set_xscale('log')

    plt.legend(fontsize=15)
    ax1.tick_params(axis='both', labelsize=15)
    ax1.set_xlabel('Mean binding site primary sequence distance', fontsize=15)
    ax1.set_ylabel('N matches', fontsize=15)

    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)
    
    plt.clf()

def plot_success_rate_vs_ligand_sizes(all_ligand_size_and_match_results, output_file=None):
    bin_width = 2
    N_bins = 25
    bin_centers = [bin_width * (i + 0.5) for i in range(N_bins)]

    fig, ax1 = plt.subplots()
    
    for k in all_ligand_size_and_match_results:
        N_rosetta_succeeds, N_fast_succeeds, N_total_results = assort_success_matches_into_bins(
                all_ligand_size_and_match_results[k], bin_centers, bin_width, use_log=False)
        fast_success_rates = [100 * N_fast_succeeds[i] / (N_total_results[i] + 0.001) for i in range(N_bins)]
        rosetta_success_rates = [100 * N_rosetta_succeeds[i] / (N_total_results[i] + 0.001) for i in range(N_bins)]

        if k == '3res_native_rossmann':
            color = 'yellow' 
            rec1 = ax1.plot(bin_centers, rosetta_success_rates, color=color, label='native Rossmann')
        elif k == '3res_2lv8':
            color = 'blue'
            rec2 = ax1.plot(bin_centers, rosetta_success_rates, color=color, label='de novo Rossmann')
        elif k == '3res_native_ntf2':
            color = 'cyan'
            rec3 = ax1.plot(bin_centers, rosetta_success_rates, color=color, label='native NTF2')
        elif k == '3res_5tpj':
            color = 'red'
            rec4 = ax1.plot(bin_centers, rosetta_success_rates, color=color, label='de novo NTF2')

    ax1.set_xlim(0, 50)
    ax1.set_ylim(0, 33)

    plt.legend(fontsize=15)
    ax1.tick_params(axis='both', labelsize=15)
    ax1.set_xlabel('N ligand heavy atoms', fontsize=15)
    ax1.set_ylabel('Match success rate (%)', fontsize=15)

    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)
    
    plt.clf()

def plot_success_rate_vs_mean_binding_site_dist(all_seq_dist_and_match_results, output_file=None):
    bin_width = 0.5
    N_bins = 14
    bin_centers = [bin_width * (i + 0.5) for i in range(N_bins)]
    linear_bin_centers = [np.exp(x) for x in bin_centers]

    fig, ax1 = plt.subplots()
    
    for k in all_seq_dist_and_match_results:
        N_rosetta_succeeds, N_fast_succeeds, N_total_results = assort_success_matches_into_bins(all_seq_dist_and_match_results[k], bin_centers, bin_width)
        fast_success_rates = [100 * N_fast_succeeds[i] / N_total_results[i] for i in range(N_bins)]
        rosetta_success_rates = [100 * N_rosetta_succeeds[i] / N_total_results[i] for i in range(N_bins)]

        if k == '3res_native_rossmann':
            color = 'yellow' 
            rec1 = ax1.plot(linear_bin_centers, rosetta_success_rates, color=color, label='native Rossmann')
        elif k == '3res_2lv8':
            color = 'blue'
            rec2 = ax1.plot(linear_bin_centers, rosetta_success_rates, color=color, label='de novo Rossmann')
        elif k == '3res_native_ntf2':
            color = 'cyan'
            rec3 = ax1.plot(linear_bin_centers, rosetta_success_rates, color=color, label='native NTF2')
        elif k == '3res_5tpj':
            color = 'red'
            rec4 = ax1.plot(linear_bin_centers, rosetta_success_rates, color=color, label='de novo NTF2')

    ax1.set_xlim(1, 1000)
    ax1.set_ylim(0, 100)
    ax1.set_xscale('log')

    plt.legend(fontsize=15)
    ax1.tick_params(axis='both', labelsize=15)
    ax1.set_xlabel('Mean binding site primary sequence distance', fontsize=15)
    ax1.set_ylabel('Match success rate (%)', fontsize=15)

    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)
    
    plt.clf()

def print_overlapping_matches_for_two_datasets(binding_site_summary_dict, dict_aggregated_results, key1, key2):
    bin_width = 0.5
    N_bins = 14
    bin_centers = [bin_width * (i + 0.5) for i in range(N_bins)]
    linear_bin_centers = [np.exp(x) for x in bin_centers]

    seq_dist_and_match_results = get_binding_site_mean_sequence_distances_and_match_results(binding_site_summary_dict, 
            [aggregated_results_to_dict(dict_aggregated_results[key1])])
    N_rosetta_succeeds_1, N_fast_succeeds_1, N_total_results_1 = assort_success_matches_into_bins(seq_dist_and_match_results, bin_centers, bin_width)

    seq_dist_and_match_results = get_binding_site_mean_sequence_distances_and_match_results(binding_site_summary_dict, 
            [aggregated_results_to_dict(dict_aggregated_results[key2])])
    N_rosetta_succeeds_2, N_fast_succeeds_2, N_total_results_2 = assort_success_matches_into_bins(seq_dist_and_match_results, bin_centers, bin_width)

    seq_dist_and_match_results = get_binding_site_mean_sequence_distances_and_match_results(binding_site_summary_dict, 
            [aggregated_results_to_dict(dict_aggregated_results[key1]), 
             aggregated_results_to_dict(dict_aggregated_results[key2])])
    N_rosetta_succeeds_overlap, N_fast_succeeds_overlap, N_total_results_overlap = assort_success_matches_into_bins(seq_dist_and_match_results, bin_centers, bin_width)

    print('Bin centers', linear_bin_centers)
    print('Total sites', N_total_results_overlap)
    print(key1, N_rosetta_succeeds_1)
    print(key2, N_rosetta_succeeds_2)
    print('overlap', N_rosetta_succeeds_overlap)



if __name__ == '__main__':
    
    binding_site_summary_file = '../binding_sites_from_pdb/3residue_binding_site_database_summary.json'

    with open(binding_site_summary_file, 'r') as f:
        bs_summary = json.load(f)
    binding_site_summary_dict = binding_site_summary_to_dict(bs_summary)

    aggregated_results_files = {
            '3res_2lv8':'matching_data/aggregated_matching_3res_site_to_2lv8_two_LHL.json',
            '3res_5tpj':'matching_data/aggregated_matching_3res_site_to_5tpj.json',
            '3res_native_rossmann':'matching_data/aggregated_matching_3res_site_to_native_rossmann.json',
            '3res_native_ntf2':'matching_data/aggregated_matching_3res_site_to_native_ntf2.json',
            }
   
    dict_aggregated_results = {}

    for k in aggregated_results_files.keys():
        with open(aggregated_results_files[k], 'r') as f:
            dict_aggregated_results[k] = json.load(f)

    all_ligand_size_and_match_results = {}
    for k in dict_aggregated_results:
        aggregated_results_dict = aggregated_results_to_dict(dict_aggregated_results[k])
        ligand_size_and_match_results = get_binding_ligand_size_and_match_results([aggregated_results_dict])
        all_ligand_size_and_match_results[k] = ligand_size_and_match_results

    #plot_success_rate_vs_ligand_sizes(all_ligand_size_and_match_results)
    plot_success_rate_vs_ligand_sizes(all_ligand_size_and_match_results, output_file='./plots/success_rate_vs_ligand_sizes.svg')

    #all_seq_dist_and_match_results = {}
    #for k in dict_aggregated_results:
    #    aggregated_results_dict = aggregated_results_to_dict(dict_aggregated_results[k])
    #    seq_dist_and_match_results = get_binding_site_mean_sequence_distances_and_match_results(binding_site_summary_dict, [aggregated_results_dict])
    #    all_seq_dist_and_match_results[k] = seq_dist_and_match_results

    #plot_mean_sequence_distance_distribution(seq_dist_and_match_results)
    #plot_num_matches_vs_mean_binding_site_dist(all_seq_dist_and_match_results, output_file=None)
    
    #plot_success_rate_vs_mean_binding_site_dist(all_seq_dist_and_match_results, output_file=None)
    #plot_success_rate_vs_mean_binding_site_dist(all_seq_dist_and_match_results, output_file='./plots/success_rate_vs_mean_binding_site_dist.png')
    #plot_success_rate_vs_mean_binding_site_dist(all_seq_dist_and_match_results, output_file='./plots/success_rate_vs_mean_binding_site_dist.svg')

    #print_overlapping_matches_for_two_datasets(binding_site_summary_dict, dict_aggregated_results, '3res_2lv8', '3res_5tpj')

