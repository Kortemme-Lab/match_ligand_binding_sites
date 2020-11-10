#!/usr/bin/env python3
'''Aggregate the binding site informations into a single json file.
'''

import os
import json


if __name__ == '__main__':

    input_path = './binding_site_database'
    output_file = 'binding_site_database_summary.json'

    aggregated_binding_sites = []

    for d1 in os.listdir(input_path):
        for d2 in os.listdir(os.path.join(input_path, d1)):
        
            for bf in os.listdir(os.path.join(input_path, d1, d2)):
                
                if bf.endswith('.json'):
                    info_file = os.path.join(input_path, d1, d2, bf)
                    bs_id = bf.split('.')[0].split('_')[-1]
                    
                    pdb_file = os.path.join(input_path, d1, d2, 'binding_site_{0}.pdb.gz'.format(bs_id))
                    ligand_params_file = os.path.join(input_path, d1, d2, 'ligand_{0}.params'.format(bs_id)) 

                    if not os.path.exists(ligand_params_file):
                        ligand_params_file = 'NA'

                    with open(info_file, 'r') as f:
                        b_info = json.load(f)

                    aggregated_binding_sites.append({'info_file': info_file,
                        'pdb_file' : pdb_file,
                        'ligand_params_file' : ligand_params_file,
                        'binding_site_info':b_info})


    with open(output_file, 'w') as f:
        json.dump(aggregated_binding_sites, f, indent='  ')
