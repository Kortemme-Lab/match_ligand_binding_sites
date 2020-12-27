#!/usr/bin/env python3
'''Plot the statistics of a matching result.
Usage:
    ./plot_matching_statistics.py aggregated_results_file
'''

import json

import numpy as np


def print_binding_sites(aggregated_results, fast_success, rosetta_success, binding_site_range=(1, 100), 
        ligand_name=None, lig_heavy_atm_range=(1,float('inf')), max_num_output=10):
    '''Get binding sites that satisfies given
    properties.
    '''
    selected_matches = []

    for match_info in aggregated_results:
        if 'num_rosetta_match_succeed' in match_info:
            binding_site_size = len(match_info['binding_site_residues'])

            if match_info['num_fast_match_succeed'] < 1 and fast_success:
                continue

            if match_info['num_rosetta_match_succeed'] < 1 and rosetta_success:
                continue

            if not(binding_site_range[0] <= binding_site_size <= binding_site_range[1]):
                continue

            if not(ligand_name is None) and match_info['ligand_name'] != ligand_name:
                continue

            if not(lig_heavy_atm_range[0] <= match_info['ligand_n_heavy_atoms'] <= lig_heavy_atm_range[1]):
                continue

            selected_matches.append(match_info)

    
    np.random.shuffle(selected_matches)

    if max_num_output is None:
        max_num_output = len(selected_matches)

    for x in selected_matches[:max_num_output]:
        print(x)


if __name__ == '__main__':
    
    aggregated_results_file = 'matching_data/aggregated_matching_3res_site_to_5tpj.json'
    #aggregated_results_file = 'matching_data/aggregated_matching_to_2lv8_two_LHL.json'
    #aggregated_results_file = 'matching_data/aggregated_matching_3res_site_to_native_ntf2.json'

    with open(aggregated_results_file, 'r') as f:
        aggregated_results = json.load(f)

    #print_binding_sites(aggregated_results, True, True, binding_site_range=(3,3), ligand_name='MES', lig_heavy_atm_range=(10,20), max_num_output=100)
    #print_binding_sites(aggregated_results, True, True, binding_site_range=(3,10), lig_heavy_atm_range=(1,1), ligand_name=' ZN', max_num_output=10)
    #print_binding_sites(aggregated_results, True, True, binding_site_range=(4,10), lig_heavy_atm_range=(10,100), max_num_output=100)
    print_binding_sites(aggregated_results, True, True, binding_site_range=(3,3), lig_heavy_atm_range=(30,100), max_num_output=100)
