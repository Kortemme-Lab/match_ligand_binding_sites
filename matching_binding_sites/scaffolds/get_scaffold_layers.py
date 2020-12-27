#!/usr/bin/env python3
'''Get the layers of scaffold residues. 
'''

import os
import json

import pyrosetta
from pyrosetta import rosetta



def print_pymol_selection_string(residue_list):
    print('sele resi ' + '+'.join([str(x) for x in residue_list]))


def residue_subset_to_list(residue_subset):
    '''Convert a residue subset to a list'''
    residue_list = []

    for i in range(1, len(residue_subset) + 1):
        
        if residue_subset[i]:
            residue_list.append(i)

    return residue_list


def init_residue_selectors():
    '''Initialize the residue selectors for core, boundary and surface layers.'''
    selector_core = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_residue_selector(
            '<Layer name="layer" select_core="true"/>')

    selector_boundary = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_residue_selector(
            '<Layer name="layer" select_boundary="true"/>')
    
    selector_surface = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_residue_selector(
            '<Layer name="surface" select_surface="true"/>')
   
    return selector_core, selector_boundary, selector_surface


def get_residue_layers_of_pose(pose, selector_core, selector_boundary, selector_surface):
    '''Get the sequence independent residue layers of a pose.
    Return a dictionary of different layers.
    '''
    residue_layers = {'core':[], 'boundary':[], 'surface':[]}
  
    # Select the core layer

    selected_residues = selector_core.apply(pose)
    residue_layers['core'] = residue_subset_to_list(selected_residues)

    # Select the boundary layer

    selected_residues = selector_boundary.apply(pose)
    residue_layers['boundary'] = residue_subset_to_list(selected_residues)

    # Select the core layer

    selected_residues = selector_surface.apply(pose)
    residue_layers['surface'] = residue_subset_to_list(selected_residues)

    return residue_layers

def dump_residue_layers_of_a_library(pdb_path, output_file):
   
    all_residue_layers = {}

    selector_core, selector_boundary, selector_surface = init_residue_selectors()

    for f in os.listdir(pdb_path):
        if f.endswith('.pdb') or f.endswith('.pdb.gz'):
            
            pose = pyrosetta.pose_from_file(os.path.join(pdb_path,f))
            residue_layers = get_residue_layers_of_pose(pose, selector_core, selector_boundary, selector_surface)

            all_residue_layers[f] = residue_layers

    with open(output_file, 'w') as f:
        json.dump(all_residue_layers, f, indent=' ')


if __name__ == '__main__':
    pyrosetta.init()

    libraries_and_output_files = [
            ('./native_rossmann_scaffolds', 'native_rossmann_scaffolds_residue_layers.json'),
            ('./native_ntf2_scaffolds', 'native_ntf2_scaffolds_residue_layers.json'),
            ('./designable_2lv8_two_LHL_1000_random', 'designable_2lv8_two_LHL_1000_random_residue_layers.json'),
            ('./designable_5tpj_no_c_term_1000_random', 'designable_5tpj_no_c_term_1000_random_residue_layers.json')
            ]

    for lao in libraries_and_output_files:
        dump_residue_layers_of_a_library(lao[0], lao[1])

