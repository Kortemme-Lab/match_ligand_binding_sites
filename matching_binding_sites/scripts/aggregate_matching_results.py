#!/usr/bin/env python3
'''Aggregate the matching results into a json file and a csv file.
Usage:
    ./aggregate_matching_results.py matching_result_dir output_file_prefix [binding_site_size] 
'''

import sys
import os

import json

import numpy as np
import pandas as pd


def aggregate_results_for_one_binding_site(match_info_file, binding_site_db, binding_site_size=None):
    '''Aggregate matching results for one binding site.'''
    if not os.path.exists(match_info_file):
        print('Cannot find file {0}'.format(match_info_file))
        return {}
  
    try:
        with open(match_info_file, 'r') as f:
            match_info = json.load(f)
    except:
        return {}

    aggregated_match_info = {
            'match_info_file' : match_info_file,
            'num_scaffold_tested' : 0,
            'num_fast_match_succeed' : 0,
            'num_rosetta_match_succeed' : 0,
            'num_rosetta_match_failed' : 0,
            'job_id' : 0,
            'success_rosetta_match_scaffolds':[],
            'n_scaffolds_to_first_fast_match':-1,
            'n_scaffolds_to_first_rosetta_match':-1,
            'average_success_fast_match_time':-1,
            'average_success_rosetta_match_time':-1,
            }

    success_fast_match_times = []
    success_rosetta_match_times = []

    for scaffold_info in match_info:
        aggregated_match_info['job_id'] = scaffold_info['job_id']
        aggregated_match_info['num_scaffold_tested'] += 1
       
        if scaffold_info['num_fast_matches'] > 0:
            aggregated_match_info['num_fast_match_succeed'] += 1
            success_fast_match_times.append(scaffold_info['fast_match_time'])

            if aggregated_match_info['n_scaffolds_to_first_fast_match'] == -1:
                aggregated_match_info['n_scaffolds_to_first_fast_match'] = aggregated_match_info['num_scaffold_tested']
        else:
            continue

        if scaffold_info['num_rosetta_matches'] > 0:
            aggregated_match_info['num_rosetta_match_succeed'] += 1
            success_rosetta_match_times.append(scaffold_info['rosetta_match_time'])
            aggregated_match_info['success_rosetta_match_scaffolds'].append(scaffold_info['scaffold'])
            
            if aggregated_match_info['n_scaffolds_to_first_rosetta_match'] == -1:
                aggregated_match_info['n_scaffolds_to_first_rosetta_match'] = aggregated_match_info['num_scaffold_tested']

        if scaffold_info['num_rosetta_matches'] == -1:
            aggregated_match_info['num_rosetta_match_failed'] += 1

    # Get the average times

    if len(success_fast_match_times) > 0:
        aggregated_match_info['average_success_fast_match_time'] = np.mean(success_fast_match_times)
    if len(success_rosetta_match_times) > 0:
        aggregated_match_info['average_success_rosetta_match_time'] = np.mean(success_rosetta_match_times)

    # Get the information from the binding site database
   
    s_match_info_file = match_info_file.split('/')
    ligand_id = s_match_info_file[-2]
    binding_site_info_file = os.path.join(binding_site_db, '/'.join(s_match_info_file[-4:-2]), 'binding_site_info_{0}.json'.format(ligand_id))

    with open(binding_site_info_file, 'r') as f:
        binding_site_info = json.load(f)

    aggregated_match_info['ligand_n_heavy_atoms'] = binding_site_info['ligand_n_heavy_atoms']
    aggregated_match_info['ligand_name'] = binding_site_info['ligand_name']

    # Get the binding site residues

    binding_site_residues = []

    if binding_site_size is None:
        for b_e in binding_site_info['binding_site_energies']:
            binding_site_residues.append(b_e[2])

    else:
        r_n_e = []
        for r in binding_site_info['binding_site_energies']:
            r_n_e.append((r[2], r[3]['weighted_total']))
    
        r_n_e_sorted = sorted(r_n_e, key=lambda x : x[1])

        binding_site_residues = [r_n_e_sorted[i][0] for i in range(binding_site_size)]

    aggregated_match_info['binding_site_residues'] = binding_site_residues

    return aggregated_match_info
    

if __name__ == '__main__':

    matching_result_dir = sys.argv[1] 
    output_file_prefix = sys.argv[2]
    
    if len(sys.argv) > 3:
        binding_site_size = int(sys.argv[3])
    else:
        binding_site_size = None

    binding_site_db = '../binding_sites_from_pdb/binding_site_database' 


    all_aggregated_match_info = []

    for d1 in os.listdir(matching_result_dir):
        for d2 in os.listdir(os.path.join(matching_result_dir, d1)):
            for d3 in os.listdir(os.path.join(matching_result_dir, d1, d2)):
                match_info_file = os.path.join(matching_result_dir, d1, d2, d3, 'all_match_info.json')

                aggregated_match_info = aggregate_results_for_one_binding_site(match_info_file, binding_site_db, binding_site_size) 
               
                all_aggregated_match_info.append(aggregated_match_info)

                if 'num_rosetta_match_failed' in aggregated_match_info and aggregated_match_info['num_rosetta_match_failed']  > 0:
                    print(aggregated_match_info)

    # Write the json file

    with open(output_file_prefix + '.json', 'w') as f:
        json.dump(all_aggregated_match_info, f, indent=2)

    # Write the table

    df = pd.DataFrame(all_aggregated_match_info)
    with open(output_file_prefix + '.tsv', 'w') as f:
        df.to_csv(f, sep='\t', index=False, header=True)

