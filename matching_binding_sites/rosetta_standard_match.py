#!/usr/bin/env python3
'''Match binding sites to scaffolds
using standard Rosetta matcher.
'''
import os
import subprocess
import shutil

import numpy as np

import pyrosetta
from pyrosetta import rosetta

def angle(x1, x2, x3):
    '''Return the angle x1-x2-x3 in degrees.
    x1, x2 and x3 are xyzVectors.
    '''
    v1 = (x1 - x2).normalized()
    v2 = (x3 - x2).normalized()

    return np.arccos(v1.dot(v2)) * 180 / np.pi

def dihedral(x1, x2, x3, x4):
    '''Return the dihedral x1-x2-x3-x4 in degrees.
    x1, x2, x3 and x4 are xyzVectors.
    '''
    v1 = (x1 - x2).normalized()
    v2 = (x3 - x2).normalized()
    v3 = (x4 - x3).normalized()

    v1_parallel = v2.normalized()
    v1_parallel *= v1.dot(v2)
    v3_parallel = v2.normalized()
    v3_parallel *= v3.dot(v2)
    
    v1_n = (v1 - v1_parallel).normalized()
    v3_n = (v3 - v3_parallel).normalized()

    cos = v1_n.dot(v3_n)
    sin = v1_n.cross(v3_n).dot(v2)

    return np.arctan2(sin, cos) * 180 / np.pi


def atom_residue_distance(pose, res, atom_id, res2):
    '''Return the minimum heavy atom distance
    between an atom and a residue.
    '''
    min_distance = float('inf')

    for a2 in range(1, pose.residue(res2).nheavyatoms() + 1):
        dist = pose.residue(res).xyz(atom_id).distance(pose.residue(res2).xyz(a2))

        if dist < min_distance:
            min_distance = dist

    return min_distance

def find_contact_heavy_atoms(pose, res1, res2):
    '''Find 3 contact heavy atoms for a residue 
     i.e. the heavy atoms that is closest to the opposing residue
     and the 2 heavy atoms that are closest to 
     this residue
    '''
    # Get the atoms that are can be considered as contact atoms for residue 1

    potential_contact_atom_ids = []

    for a in range(1, pose.residue(res1).nheavyatoms() + 1):
        # Exclude protein backbone atoms except for N, CA and C
        # If backbon variants are included, for example OXT, the
        # matcher application would crash
        
        if pose.residue(res1).is_protein():
            if pose.residue(res1).atom_is_backbone(a):
                if pose.residue(res1).atom_name(a).strip() in ['N', 'CA', 'C']:
                    potential_contact_atom_ids.append(a)
            elif not pose.residue(res1).is_virtual(a): # Ignore virtual atoms for protein residues
                potential_contact_atom_ids.append(a)

        else:
            potential_contact_atom_ids.append(a)

    # Calculate distances for all heavy atoms
    
    distances = []

    for a in potential_contact_atom_ids:
        distances.append((a, atom_residue_distance(pose, res1, a, res2)))

    # Get the names of the closest

    distances_sorted = sorted(distances, key=lambda x : x[1])

    contact_atoms = []
    contact_atoms.append(pose.residue(res1).atom_name(distances_sorted[0][0]))

    # Get the two closest intra-residue atoms

    intra_distances = []

    for a in potential_contact_atom_ids:
        if a == distances_sorted[0][0]:
            continue

        intra_distances.append((a, pose.residue(res1).xyz(distances_sorted[0][0]).distance(
            pose.residue(res1).xyz(a))))

    intra_distances_sorted = sorted(intra_distances, key=lambda x:x[1])

    for i in range(2):
        contact_atoms.append(pose.residue(res1).atom_name(intra_distances_sorted[i][0]))

    return contact_atoms

def generate_cst_file(site_pose, output_file):
    '''Generate the constraint file for Rosetta matching.'''
    # Find the ligand and protein residues

    ligand_residue = 0
    protein_residues = []
    for i in range(1, site_pose.size() + 1):
        if site_pose.residue(i).is_ligand():
            ligand_residue = i
        elif site_pose.residue(i).is_protein():
            protein_residues.append(i)

    assert(ligand_residue != 0)

    # Generate constraints between the ligand and binding site protein residues 

    template_sections = \
'''CST::BEGIN
  TEMPLATE::   ATOM_MAP: 1 atom_name: {a1_1} {a1_2} {a1_3}
  TEMPLATE::   ATOM_MAP: 1 residue3: {r1}

  TEMPLATE::   ATOM_MAP: 2 atom_name: {a2_1} {a2_2} {a2_3}
  TEMPLATE::   ATOM_MAP: 2 residue3: {r2}

  CONSTRAINT:: distanceAB: {dist:7.2f}   5.00 100.00  0        0
  CONSTRAINT::    angle_A: {angA:7.2f}  10.00 100.00  360.00   1
  CONSTRAINT::    angle_B: {angB:7.2f}  10.00 100.00  360.00   1
  CONSTRAINT::  torsion_A: {torA:7.2f}  10.00 100.00  360.00   1
  CONSTRAINT::  torsion_B: {torB:7.2f}  10.00 100.00  360.00   1
  CONSTRAINT:: torsion_AB: {toAB:7.2f}  10.00 100.00  360.00   1
CST::END'''

    constraint_sections = []

    for pro_res in protein_residues:
        ligand_contact_atoms = find_contact_heavy_atoms(site_pose, ligand_residue, pro_res)
        protein_contact_atoms = find_contact_heavy_atoms(site_pose, pro_res, ligand_residue)

        # Get the positions of contact atoms

        v_r1a1 = site_pose.residue(ligand_residue).xyz(ligand_contact_atoms[0])
        v_r1a2 = site_pose.residue(ligand_residue).xyz(ligand_contact_atoms[1])
        v_r1a3 = site_pose.residue(ligand_residue).xyz(ligand_contact_atoms[2])

        v_r2a1 = site_pose.residue(pro_res).xyz(protein_contact_atoms[0])
        v_r2a2 = site_pose.residue(pro_res).xyz(protein_contact_atoms[1])
        v_r2a3 = site_pose.residue(pro_res).xyz(protein_contact_atoms[2])

        # Calculate the constraint distance, angles and torsions

        dAB = v_r1a1.distance(v_r2a1)
        aA = angle(v_r1a2, v_r1a1, v_r2a1)
        aB = angle(v_r1a1, v_r2a1, v_r2a2)
        tA = dihedral(v_r1a3, v_r1a2, v_r1a1, v_r2a1)
        tB = dihedral(v_r1a1, v_r2a1, v_r2a2, v_r2a3)
        tAB = dihedral(v_r1a2, v_r1a1, v_r2a1, v_r2a2)

        # Generate the constraint file section

        constraint_section = template_sections.format(
                a1_1=ligand_contact_atoms[0], a1_2=ligand_contact_atoms[1], a1_3=ligand_contact_atoms[2], r1=site_pose.residue(ligand_residue).name3(),
                a2_1=protein_contact_atoms[0], a2_2=protein_contact_atoms[1], a2_3=protein_contact_atoms[2], r2=site_pose.residue(pro_res).name3(),
                dist=dAB, angA=aA, angB=aB, torA=tA, torB=tB, toAB=tAB)

        constraint_sections.append(constraint_section)

    # Write the constraint file

    with open(output_file, 'w') as f:
        f.write('\n\n'.join(constraint_sections))

def write_params_file_for_ligand(binding_site_pdb, output_file):
    '''Write a params file for the ligand residue in the binding site.
    Carefull!!! The default output is not aligned and does not work for matcher.
    '''
    site_pose = pyrosetta.pose_from_file(binding_site_pdb) 
    
    # Find the ligand and protein residues

    ligand_residue = 0
    for i in range(1, site_pose.size() + 1):
        if site_pose.residue(i).is_ligand():
            ligand_residue = i

    assert(ligand_residue != 0)

    # Write the .params file

    rosetta.core.chemical.write_topology_file(site_pose.residue(ligand_residue).type(), output_file)

def standard_rosetta_match(scratching_path, scaffold_pdb, constraint_file, ligand_name, params_file=None, keep_outputs=False):
    '''Run standard Rosetta matching.
    Return the number of matches.
    Return -1 if matching failed.
    '''
    abs_scaffold_pdb = os.path.abspath(scaffold_pdb)
    abs_constraint_file = os.path.abspath(constraint_file)
    
    if not (params_file is None):
        abs_params_file = os.path.abspath(params_file)
    
    matcher_app = 'match.linuxgccrelease'

    # Go to the output path

    cwd = os.getcwd()

    os.makedirs(scratching_path, exist_ok=True)
    os.chdir(scratching_path)

    # Generate pos file

    scaffold_pose = pyrosetta.pose_from_file(abs_scaffold_pdb)
    scaffold_protein_residues = []

    for i in range(1, scaffold_pose.size() + 1):
        if scaffold_pose.residue(i).is_protein():
            scaffold_protein_residues.append(str(i))

    pos_file = 'scaffold.pos'
    with open(pos_file, 'w') as f:
        f.write(' '.join(scaffold_protein_residues))

    # Run matcher
        
    matcher_cmd = [matcher_app,
        '-match:output_format', 'PDB',
        '-match:match_grouper', 'SameSequenceGrouper',
        '-match:consolidate_matches',
        '-match:output_matches_per_group', '1',
        '-use_input_sc',
        '-in:ignore_unrecognized_res',
        '-ex1', '-ex2',
        '-enumerate_ligand_rotamers', 'false',
        '-match::lig_name', ligand_name,
        '-match:geometric_constraint_file', abs_constraint_file,
        '-s', abs_scaffold_pdb,
        '-match::scaffold_active_site_residues', pos_file
        ]

    if not (params_file is None):
        matcher_cmd += ['-extra_res_fa', abs_params_file]

    return_code = subprocess.call(matcher_cmd)

    # Count the number of matches

    num_matches = 0

    for f in os.listdir('.'):
        if f.startswith('UM') and f.endswith('.pdb'):
            num_matches += 1

    # Clear the scratch path

    os.chdir(cwd)

    if not keep_outputs:
        shutil.rmtree(scratching_path)

    if 0 == return_code:
        return num_matches

    else:
        return -1

if __name__ == '__main__':
    pyrosetta.init()

    binding_site_pdb = './test_inputs/test_site.pdb'
    site_pose = pyrosetta.pose_from_file(binding_site_pdb) 

    generate_cst_file(site_pose, 'test.cst')
    ##write_params_file_for_ligand(binding_site_pdb, 'test.params')

    
    num_matches = standard_rosetta_match('test_rosetta_match', './test_inputs/2src_scaffold.pdb', './test_inputs/test_short.cst', 'ANP', 
            './test_inputs/debug_params/params_file_from_mol2/ANP.params',
            #keep_outputs=False)
            keep_outputs=True)

    print('Found {0} matches.'.format(num_matches))
