#!/usr/bin/env python3
'''Match binding sites from pdb to de novo
scaffolds.
Usage:
    ./match_binding_sites_to_scaffolds.py scaffold_path output_path job_id
'''

import sys
import os
import time
import json
from optparse import OptionParser

import pyrosetta
from pyrosetta import rosetta

import fast_match
import rosetta_standard_match


def extract_subpose(pose, residues):
    '''Extract a sub-pose from a pose.
    Return the new pose.
    '''
    seqposes = rosetta.utility.vector1_unsigned_long()
    for seqpos in residues:
        seqposes.append(seqpos)

    new_pose = rosetta.core.pose.Pose()
    rosetta.core.pose.pdbslice(new_pose, pose, seqposes)

    return new_pose

def get_sub_binding_site_pose(pdb_file, bs_info, n_residue_to_keep):
    '''Get a pose of a binding site that has
    only n_residue_to_keep low energy protein residues.

    Return: (pose, protein_residues)

    If there are less than n_residue_to_keep protein
    residues in the pdb file, return None.
    '''

    r_n_e = []
    for r in bs_info['binding_site_info']['binding_site_energies']:
        r_n_e.append((r[1], r[3]['weighted_total']))

    # Return None if there are not enough protein residues

    if len(r_n_e) < n_residue_to_keep:
        return (None, None)

    # Extract the sub-pose

    r_n_e_sorted = sorted(r_n_e, key=lambda x : x[1])
    ligand = bs_info['binding_site_info']['ligand_number']
    
    pose = pyrosetta.pose_from_file(pdb_file)
  
    residues_to_keep = [ligand] + [r_n_e_sorted[i][0] for i in range(n_residue_to_keep)]
    residues_to_keep_pose = []

    for r in residues_to_keep:
        chain = r.split()[1]
        res = int(r.split()[0])

        residues_to_keep_pose.append(pose.pdb_info().pdb2pose(chain, res))
  
    new_pose = extract_subpose(pose, residues_to_keep_pose)
    
    # Get the protein residues of the new pose

    protein_residues = [r_n_e_sorted[i][0] for i in range(n_residue_to_keep)]
    new_protein_residues_pose = []

    for r in protein_residues:
        chain = r.split()[1]
        res = int(r.split()[0])

        new_protein_residues_pose.append(new_pose.pdb_info().pdb2pose(chain, res))

    return new_pose, new_protein_residues_pose

if __name__ == '__main__':

    # Parse the options

    parser = OptionParser()
    parser.add_option("-n", "--n_residue_to_keep", dest="n_residue_to_keep", action="store", default=None,
            help="The number of binding sites to keep")
    parser.add_option('-d', '--dump_matches', dest='dump_matches', action='store_true', default=False,
            help='Dump the matches')

    (options, args) = parser.parse_args()

    scaffold_path = args[0]
    output_path = args[1]
    job_id = int(args[2])

    binding_site_summary_file = '../binding_sites_from_pdb/binding_site_database_summary.json'
    binding_site_base = '../binding_sites_from_pdb'

    try:
        main_scratch_path = os.environ['TMPDIR']
    except:
        main_scratch_path = None

    pyrosetta.init()

    with open(binding_site_summary_file, 'r') as f:
        bs_summary = json.load(f)

    # Run matching for the binding site specified by the job_id

    bs_info = bs_summary[job_id - 1]

    # Get the input files for matching 

    pdb_file = os.path.join(binding_site_base, bs_info['pdb_file'])
    ligand_name = bs_info['binding_site_info']['ligand_name'].strip()

    if bs_info['ligand_params_file'] == 'NA':
       ligand_params_file = None
    else:
        ligand_params_file = os.path.join(binding_site_base, bs_info['ligand_params_file'])

    # Load the binding site pose

    if options.n_residue_to_keep is None:
        n_residue_to_keep = max(2, len(bs_info['binding_site_info']['binding_site_energies']))
    else:
        n_residue_to_keep = int(options.n_residue_to_keep)

    site_pose, site_protein_residues = get_sub_binding_site_pose(pdb_file, bs_info, n_residue_to_keep)
    
    if site_pose is None: exit()
    
    d1 = pdb_file.split('/')[-3]  
    d2 = pdb_file.split('/')[-2]
    d3 = pdb_file.split('/')[-1].split('.')[0].split('_')[-1]

    # Generate the constraint file

    binding_site_constraint_file = os.path.join(output_path, d1, d2, d3, 'binding_site.cst')
    os.makedirs(os.path.join(output_path, d1, d2, d3), exist_ok=True)
    rosetta_standard_match.generate_cst_file(site_pose, binding_site_constraint_file)  
    site_pose.dump_pdb(os.path.join(output_path, d1, d2, d3, 'binding_site.pdb.gz'))

    # Iterate through all scaffold_path

    all_match_info = []

    for f in os.listdir(scaffold_path):
        if f.endswith('.pdb') or f.endswith('.pdb.gz'):
            match_start_time = time.time()
            
            match_info = {}
            match_info['job_id'] = job_id
            match_info['scaffold'] = f
            match_info['scaffold'] = f

            d4 = 'match_scaffold_{0}'.format(f.split('.')[0])
            
            match_output_path = os.path.join(output_path, d1, d2, d3, d4)
            scaffold_pose = pyrosetta.pose_from_file(os.path.join(scaffold_path, f))
       
            # Run fast match

            scaffold_protein_residues = [i for i in range(1, scaffold_pose.size() + 1) if scaffold_pose.residue(i).is_protein()]

            num_fast_matches, fast_matched_positions = fast_match.fast_match(site_pose, site_protein_residues, scaffold_pose, scaffold_protein_residues,
                    match_output_path, dump_matches=options.dump_matches)
            match_info['num_fast_matches'] = num_fast_matches
            match_info['fast_matched_positions'] = fast_matched_positions

            fast_match_finish_time = time.time()

            # Run standard Rosetta match if fast match found hits 

            if num_fast_matches > 0:
               
                if main_scratch_path is None:
                    rosetta_match_scratch_path = os.path.join(output_path, d1, d2, d3, 'rosetta_match_scaffold_{0}'.format(f.split('.')[0])) 
                else:
                    rosetta_match_scratch_path = os.path.join(main_scratch_path, 'rosetta_match_scaffold_{0}'.format(f.split('.')[0]))

                num_rosetta_matches, rosetta_matched_positions = rosetta_standard_match.standard_rosetta_match(rosetta_match_scratch_path, 
                        os.path.join(scaffold_path, f), binding_site_constraint_file, ligand_name, ligand_params_file, keep_outputs=options.dump_matches)
                
                # Try a different ligand name if the match failed

                if ligand_params_file is None and num_rosetta_matches < 0:

                    num_rosetta_matches, rosetta_matched_positions = rosetta_standard_match.standard_rosetta_match(rosetta_match_scratch_path, 
                            os.path.join(scaffold_path, f), binding_site_constraint_file, 'pdb_' + ligand_name, ligand_params_file, keep_outputs=options.dump_matches)

                match_info['num_rosetta_matches'] = num_rosetta_matches
                match_info['rosetta_matched_positions'] = rosetta_matched_positions

            rosetta_standard_match_finish_time = time.time()

            match_info['fast_match_time'] = int(fast_match_finish_time - match_start_time)
            match_info['rosetta_match_time'] = int(rosetta_standard_match_finish_time - fast_match_finish_time)

            all_match_info.append(match_info)

            # Break if found rosetta matches

            if num_fast_matches > 0 and num_rosetta_matches > 0:
                break


    os.makedirs(os.path.join(output_path, d1, d2, d3), exist_ok=True)
    with open(os.path.join(output_path, d1, d2, d3, 'all_match_info.json'), 'w') as f:
        json.dump(all_match_info, f, indent=' '*4)
