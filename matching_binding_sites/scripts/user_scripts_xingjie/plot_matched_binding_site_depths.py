#!/usr/bin/env python3

import os
import json

import numpy as np
import matplotlib.pyplot as plt

def get_single_binding_site_layers(binding_site_residues, scaffold_layers):
    binding_site_layers = []

    for r in binding_site_residues:
        if r in scaffold_layers['core']:
            binding_site_layers.append('C')
        elif r in scaffold_layers['boundary']:
            binding_site_layers.append('B')
        elif r in scaffold_layers['surface']:
            binding_site_layers.append('S')

    return binding_site_layers


def get_matched_binding_site_layers(scaffold_layers_file, matching_result_dir):
    '''Load the binding site layers for all the matches found in the matching_result_dir.
    Return a list of fast matched binding site layers and a list of Rosetta matched binding
    site layers.
    Layers of each binding site is represented as a tuple of letters indicating layers for
    each binding site residue.
    '''
    with open(scaffold_layers_file, 'r') as f:
        all_scaffold_layers = json.load(f) 

    fast_matched_binding_site_layers = []
    rosetta_matched_binding_site_layers = []

    # Load the matching results

    for d1 in os.listdir(matching_result_dir):
        for d2 in os.listdir(os.path.join(matching_result_dir, d1)):
            for d3 in os.listdir(os.path.join(matching_result_dir, d1, d2)):
                match_info_file = os.path.join(matching_result_dir, d1, d2, d3, 'all_match_info.json')
   
                try:
                    with open(match_info_file, 'r') as f:
                        match_info = json.load(f)
                except:
                    continue

                for scaffold_info in match_info:
                    if scaffold_info['num_fast_matches'] > 0:
                        scaffold_layers = all_scaffold_layers[scaffold_info['scaffold']] 
                        
                        for bs in scaffold_info['fast_matched_positions']:
                            fast_matched_binding_site_layers.append(get_single_binding_site_layers(bs, scaffold_layers))

                        for bs in scaffold_info['rosetta_matched_positions']:
                            rosetta_matched_binding_site_layers.append(get_single_binding_site_layers(bs, scaffold_layers))

    return fast_matched_binding_site_layers, rosetta_matched_binding_site_layers

def dump_matched_binding_site_layers(output_file, fast_matched_binding_site_layers, rosetta_matched_binding_site_layers):
    with open(output_file, 'w') as f:
        json.dump({'fast_matched_binding_site_layers':fast_matched_binding_site_layers,
            'rosetta_matched_binding_site_layers':rosetta_matched_binding_site_layers}, f, indent=' ')

def load_matched_binding_site_layers(input_file):
    with open(input_file, 'r') as f:
        data = json.load(f)

    return data['fast_matched_binding_site_layers'], data['rosetta_matched_binding_site_layers']

def get_scaffold_layer_frequencies(scaffold_layers_file):
    with open(scaffold_layers_file, 'r') as f:
        all_scaffold_layers = json.load(f) 

    layer_count = {'core': 0, 'boundary': 0, 'surface': 0}

    for k in all_scaffold_layers:
        layer_count['core'] += len(all_scaffold_layers[k]['core'])
        layer_count['boundary'] += len(all_scaffold_layers[k]['boundary'])
        layer_count['surface'] += len(all_scaffold_layers[k]['surface'])

    N_total = layer_count['core'] + layer_count['boundary'] + layer_count['surface']

    layer_frequencies = {'core': layer_count['core'] / N_total,
            'boundary': layer_count['boundary'] / N_total,
            'surface': layer_count['surface'] / N_total}

    return layer_frequencies

def get_binding_site_layer_frequencies(layers):
    layer_count = {'core': 0, 'boundary': 0, 'surface': 0}

    for l in layers:
        for ll in l:
            if ll == 'C':
                layer_count['core'] += 1
            elif ll == 'B':
                layer_count['boundary'] += 1
            elif ll == 'S':
                layer_count['surface'] += 1

    N_total = layer_count['core'] + layer_count['boundary'] + layer_count['surface']

    layer_frequencies = {'core': layer_count['core'] / N_total,
            'boundary': layer_count['boundary'] / N_total,
            'surface': layer_count['surface'] / N_total}
    
    return layer_frequencies

def get_binding_site_depth_scores(layers):
    depth_scores = []

    for l in layers:
        depth_score = 0
        for ll in l:
            if ll == 'C':
                depth_score += 2
            elif ll == 'B':
                depth_score += 1

        depth_scores.append(depth_score)

    return depth_scores

def plot_scaffold_layer_distributions(scaffold_layers_files, output_file=None):

    bins = np.array([1, 2, 3])
    bin_width = 0.15

    fig, ax1 = plt.subplots()
    
    for k in scaffold_layers_files:
        layer_frequencies = get_scaffold_layer_frequencies(scaffold_layers_files[k])
        list_layer_frequencies = [100 * layer_frequencies['surface'], 100 * layer_frequencies['boundary'], 100 * layer_frequencies['core']]
       
        print(k, layer_frequencies)

        if k == 'native_rossmann':
            shift = -1.5
            color = 'yellow' 
            rec1 = ax1.bar(np.array(bins) + 0.5 + shift * bin_width, list_layer_frequencies, width=bin_width, color=color, alpha=1)
        elif k == '2lv8':
            shift = -0.5
            color = 'blue'
            rec2 = ax1.bar(np.array(bins) + 0.5 + shift * bin_width, list_layer_frequencies, width=bin_width, color=color, alpha=1)
        elif k == 'native_ntf2':
            shift = 0.5
            color = 'cyan'
            rec3 = ax1.bar(np.array(bins) + 0.5 + shift * bin_width, list_layer_frequencies, width=bin_width, color=color, alpha=1)
        elif k == '5tpj':
            shift = 1.5
            color = 'red'
            rec4 = ax1.bar(np.array(bins) + 0.5 + shift * bin_width, list_layer_frequencies, width=bin_width, color=color, alpha=1)

    plt.legend([rec1, rec2, rec3, rec4], ['native Rossmann', 'de novo Rossmann', 'native NTF2', 'de novo NTF2'], prop={'size': 15})
    plt.xticks([1.5, 2.5, 3.5], ['surface', 'boundary', 'core'])
   
    ax1.set_ylim(0, 70)
    ax1.tick_params(axis='both', labelsize=15)
    ax1.set_xlabel('Scaffold residue layer', fontsize=15)
    ax1.set_ylabel('Frequency (%)', fontsize=15)

    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)
    
    plt.clf()

def plot_binding_site_layer_distribution(scaffold_layers_files, binding_site_layer_files, output_file=None):

    bins = np.array([1, 2, 3])
    bin_width = 0.15

    fig, ax1 = plt.subplots()
    
    for k in scaffold_layers_files:
        fast_matched_binding_site_layers, rosetta_matched_binding_site_layers = load_matched_binding_site_layers(binding_site_layer_files[k])
        
        layer_frequencies = get_binding_site_layer_frequencies(rosetta_matched_binding_site_layers)
        list_layer_frequencies = [100 * layer_frequencies['surface'], 100 * layer_frequencies['boundary'], 100 * layer_frequencies['core']]
        
        print(k, layer_frequencies)
        
        if k == 'native_rossmann':
            shift = -1.5
            color = 'yellow' 
            rec1 = ax1.bar(np.array(bins) + 0.5 + shift * bin_width, list_layer_frequencies, width=bin_width, color=color, alpha=1)
        elif k == '2lv8':
            shift = -0.5
            color = 'blue'
            rec2 = ax1.bar(np.array(bins) + 0.5 + shift * bin_width, list_layer_frequencies, width=bin_width, color=color, alpha=1)
        elif k == 'native_ntf2':
            shift = 0.5
            color = 'cyan'
            rec3 = ax1.bar(np.array(bins) + 0.5 + shift * bin_width, list_layer_frequencies, width=bin_width, color=color, alpha=1)
        elif k == '5tpj':
            shift = 1.5
            color = 'red'
            rec4 = ax1.bar(np.array(bins) + 0.5 + shift * bin_width, list_layer_frequencies, width=bin_width, color=color, alpha=1)

    plt.legend([rec1, rec2, rec3, rec4], ['native Rossmann', 'de novo Rossmann', 'native NTF2', 'de novo NTF2'], prop={'size': 15})
    plt.xticks([1.5, 2.5, 3.5], ['surface', 'boundary', 'core'])
    
    ax1.set_ylim(0, 70)
    ax1.tick_params(axis='both', labelsize=15)
    ax1.set_xlabel('Matched residue layer', fontsize=15)
    ax1.set_ylabel('Frequency (%)', fontsize=15)

    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)
    
    plt.clf()


def plot_binding_site_depth_distribution(all_depth_scores, output_file=None):
    all_freq_hists = {}

    bins = np.array(list(range(8))) - 0.5
    bin_width = 0.1
    
    for k in all_depth_scores:
        hist, bin_edges = np.histogram(all_depth_scores[k], bins=bins)
        N = np.sum(hist)
        freqs = 100 * np.array(hist) / N
        all_freq_hists[k] = freqs

    fig, ax1 = plt.subplots()
    

    for k in all_depth_scores:
        print(k, all_freq_hists[k])
        
        if k == 'native_rossmann':
            shift = -1.5
            color = 'yellow' 
            rec1 = ax1.bar(np.array(bins[:-1]) + 0.5 + shift * bin_width, all_freq_hists[k], width=bin_width, color=color, alpha=1)
        elif k == '2lv8':
            shift = -0.5
            color = 'blue'
            rec2 = ax1.bar(np.array(bins[:-1]) + 0.5 + shift * bin_width, all_freq_hists[k], width=bin_width, color=color, alpha=1)
        elif k == 'native_ntf2':
            shift = 0.5
            color = 'cyan'
            rec3 = ax1.bar(np.array(bins[:-1]) + 0.5 + shift * bin_width, all_freq_hists[k], width=bin_width, color=color, alpha=1)
        elif k == '5tpj':
            shift = 1.5
            color = 'red'
            rec4 = ax1.bar(np.array(bins[:-1]) + 0.5 + shift * bin_width, all_freq_hists[k], width=bin_width, color=color, alpha=1)


    plt.legend([rec1, rec2, rec3, rec4], ['native Rossmann', 'de novo Rossmann', 'native NTF2', 'de novo NTF2'], prop={'size': 15})

    ax1.tick_params(axis='both', labelsize=15)
    ax1.set_xlabel('Match depth score', fontsize=15)
    ax1.set_ylabel('Frequency (%)', fontsize=15)

    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)
    
    plt.clf()



if __name__ == '__main__':

    scaffold_layers_files = {
            '2lv8':'./scaffolds/designable_2lv8_two_LHL_1000_random_residue_layers.json',
            '5tpj':'./scaffolds/designable_5tpj_no_c_term_1000_random_residue_layers.json',
            'native_rossmann':'./scaffolds/native_rossmann_scaffolds_residue_layers.json',
            'native_ntf2':'./scaffolds/native_ntf2_scaffolds_residue_layers.json'
        }
    matching_result_dirs = {
            '2lv8':'./matching_data/matching_to_2lv8_two_LHL_no_dump',
            '5tpj':'./matching_data/matching_site_to_5tpj_no_dump',
            'native_rossmann':'./matching_data/matching_site_to_native_rossmann_no_dump',
            'native_ntf2':'./matching_data/matching_site_to_native_ntf2_no_dump'
        }
    matching_3res_result_dirs = {
            '2lv8':'./matching_data/matching_3res_site_to_2lv8_two_LHL_no_dump',
            '5tpj':'./matching_data/matching_3res_site_to_5tpj_no_dump',
            'native_rossmann':'./matching_data/matching_3res_site_to_native_rossmann_no_dump',
            'native_ntf2':'./matching_data/matching_3res_site_to_native_ntf2_no_dump'
        }
    binding_site_layer_files = {
            '2lv8':'./matching_data/matching_to_2lv8_two_LHL_no_dump_layers.json',
            '5tpj':'./matching_data/matching_site_to_5tpj_no_dump_layers.json',
            'native_rossmann':'./matching_data/matching_site_to_native_rossmann_no_dump_layers.json',
            'native_ntf2':'./matching_data/matching_site_to_native_ntf2_no_dump_layers.json'
        }
    binding_3res_site_layer_files = {
            '2lv8':'./matching_data/matching_3res_site_to_2lv8_two_LHL_no_dump_layers.json',
            '5tpj':'./matching_data/matching_3res_site_to_5tpj_no_dump_layers.json',
            'native_rossmann':'./matching_data/matching_3res_site_to_native_rossmann_no_dump_layers.json',
            'native_ntf2':'./matching_data/matching_3res_site_to_native_ntf2_no_dump_layers.json'
        }

    # Dump all the layer files

#    for k in scaffold_layers_files:
#        fast_matched_binding_site_layers, rosetta_matched_binding_site_layers = get_matched_binding_site_layers(
#                scaffold_layers_files[k], matching_result_dirs[k])
#        dump_matched_binding_site_layers(binding_site_layer_files[k], fast_matched_binding_site_layers, rosetta_matched_binding_site_layers)
     
#    for k in scaffold_layers_files:
#        fast_matched_binding_site_layers, rosetta_matched_binding_site_layers = get_matched_binding_site_layers(
#                scaffold_layers_files[k], matching_3res_result_dirs[k])
#        dump_matched_binding_site_layers(binding_3res_site_layer_files[k], fast_matched_binding_site_layers, rosetta_matched_binding_site_layers)
     
    # Plot scaffold layer distributions
   
#    plot_scaffold_layer_distributions(scaffold_layers_files, output_file='./plots/scaffold_layer_distributions.svg')

    # Plot matching layer distributions

#    plot_binding_site_layer_distribution(scaffold_layers_files, binding_3res_site_layer_files, output_file='./plots/match_layer_distributions.svg')

    # Plot matching depth distributions

#    all_depth_scores = {}
#
#    for k in scaffold_layers_files:
#        fast_matched_binding_site_layers, rosetta_matched_binding_site_layers = load_matched_binding_site_layers(binding_site_layer_files[k])
#        
#        #get_scaffold_layer_frequencies(scaffold_layers_files[k])
#        #get_binding_site_layer_frequencies(fast_matched_binding_site_layers)
#        #get_binding_site_layer_frequencies(rosetta_matched_binding_site_layers)
#
#        depth_scores = get_binding_site_depth_scores(rosetta_matched_binding_site_layers)
#   
#        all_depth_scores[k] = depth_scores
#
#    plot_binding_site_depth_distribution(all_depth_scores, output_file='./plots/match_depth_distribution.svg')

    all_depth_scores = {}

    for k in scaffold_layers_files:
        fast_matched_binding_site_layers, rosetta_matched_binding_site_layers = load_matched_binding_site_layers(binding_3res_site_layer_files[k])
        depth_scores = get_binding_site_depth_scores(rosetta_matched_binding_site_layers)
        all_depth_scores[k] = depth_scores

    plot_binding_site_depth_distribution(all_depth_scores, output_file='./plots/match_depth_3res_distribution.svg')

