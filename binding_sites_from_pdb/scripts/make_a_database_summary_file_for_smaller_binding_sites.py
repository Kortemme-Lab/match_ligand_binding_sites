#!/usr/bin/env python3
'''
Make a new database summary table such
that each binding site has a user specified
number of protein residues.
'''

import copy
import json

def get_sub_binding_site_info(bs_info, n_residue_to_keep):
    '''Return a binding site info dictionary
    of n_residue_to_keep protein residues with lowest energies.
    '''

    # Return None if there are not enough protein residues

    if len(bs_info['binding_site_info']['binding_site_energies']) < n_residue_to_keep:
        return None

    # Find the lowest energy residues to keep

    energies_sorted = sorted(bs_info['binding_site_info']['binding_site_energies'], key=lambda x : x[3]['weighted_total'])
    
    sub_bs_info = copy.deepcopy(bs_info)
    sub_bs_info['binding_site_info']['binding_site_energies'] = copy.deepcopy(energies_sorted[:n_residue_to_keep])

    return sub_bs_info


if __name__ == '__main__':

    binding_site_summary_file = 'binding_site_database_summary.json'
    new_binding_site_summary_file = '3residue_binding_site_database_summary.json'

    with open(binding_site_summary_file, 'r') as f:
       bs_summary = json.load(f)

    print('Found {0} binding sites in the file {1}'.format(len(bs_summary), binding_site_summary_file))
    
    new_bs_summary = []

    for bs_info in bs_summary:
        sub_bs_info = get_sub_binding_site_info(bs_info, 3)
        if not (sub_bs_info is None):
            new_bs_summary.append(sub_bs_info)

    with open(new_binding_site_summary_file, 'w') as f:
        json.dump(new_bs_summary, f, indent='  ')

    print('Wrote {0} binding sites to {1}'.format(len(new_bs_summary), new_binding_site_summary_file))



