#!/usr/bin/env python3

import sys
import os
import json
import subprocess

import pyrosetta
from pyrosetta import rosetta

def real_heavy_atoms(pose, residue):
    '''Return the indices of heavy atoms that are not virtual.'''
    rh_atoms = []

    for i in range(1, pose.residue(residue).nheavyatoms() + 1):
        if not pose.residue(residue).is_virtual(i):
            rh_atoms.append(i)

    return rh_atoms

def min_heavy_atom_b_factor(pose, residue):
    '''Return the minimal heavy atom b factor of a residue.'''
    rh_atoms = real_heavy_atoms(pose, residue)
    b_factors = [pose.pdb_info().bfactor(residue, i) for i in rh_atoms]

    return min(b_factors)

def average_heavy_atom_b_factor(pose, residue):
    '''Return the heavy atom b factor of a residue.'''
    rh_atoms = real_heavy_atoms(pose, residue)

    b_factors = [pose.pdb_info().bfactor(residue, i) for i in rh_atoms]

    return sum(b_factors) / len(b_factors)

def minimize_pose(pose):
    '''Minimize a pose.'''
    mm = rosetta.core.kinematics.MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)
    mm.set_jump(True)
    min_opts = rosetta.core.optimization.MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, True )
 
    min_mover = rosetta.protocols.minimization_packing.MinMover()
    min_mover = rosetta.protocols.minimization_packing.MinMover()
    min_mover.movemap(mm)

    min_mover.apply(pose)

def fast_relax(pose, fast_relax_rounds=1, use_ex_rotamers=False):
    '''Relax a pose.'''
    rosetta.basic.options.set_boolean_option('relax:constrain_relax_to_start_coords', True)
    rosetta.basic.options.set_boolean_option('relax:ramp_constraints', True)

    task_factory = rosetta.core.pack.task.TaskFactory()
    
    if use_ex_rotamers:
        ers = rosetta.core.pack.task.operation.ExtraRotamersGeneric()
        ers.ex1(True)
        ers.ex2(True)
        ers.extrachi_cutoff(18)
        task_factory.push_back(ers)
    
    task_factory.push_back(rosetta.core.pack.task.operation.RestrictToRepacking())
    lac = rosetta.protocols.task_operations.LimitAromaChi2Operation()
    task_factory.push_back(lac)

    sfxn = rosetta.core.scoring.get_score_function()
    fast_relax_mover = rosetta.protocols.relax.FastRelax(sfxn, fast_relax_rounds)

    fast_relax_mover.apply(pose)

def find_ligand_residues(pose, size_min=0, size_max=float('inf')):
    '''Find the ligand residues in a pose. The number of heavy atoms of the
    ligand should be within the range [size_min, size_max]. 
    Return a list of residue numbers.
    '''

    ligand_residues = []

    for i in range(1, pose.size() + 1):
        if pose.residue(i).is_ligand():
            if size_min <= len(real_heavy_atoms(pose, i)) <= size_max:
                ligand_residues.append(i)

    return ligand_residues

def min_heavy_atom_distance(pose, residue1, residue2):
    '''Return the minimum heavy atom distance
    between two residues.
    '''
    min_distance = float('inf')

    for a1 in real_heavy_atoms(pose, residue1):
        for a2 in real_heavy_atoms(pose, residue2):
            dist = pose.residue(residue1).xyz(a1).distance(pose.residue(residue2).xyz(a2))

            if dist < min_distance:
                min_distance = dist

    return min_distance

def residue_residue_energies(pose, res1, res2):
    '''Get the interaction energies between two residues.'''
    e_edge = pose.energies().energy_graph().find_energy_edge(res1, res2)

    energies = {}

    if e_edge:
        energies['fa_atr'] = e_edge[rosetta.core.scoring.fa_atr]
        energies['fa_rep'] = e_edge[rosetta.core.scoring.fa_rep]
        energies['fa_sol'] = e_edge[rosetta.core.scoring.fa_sol]
        energies['fa_elec'] = e_edge[rosetta.core.scoring.fa_elec]
        energies['lk_ball_wtd'] = e_edge[rosetta.core.scoring.lk_ball_wtd]
        energies['hbond_bb_sc'] = e_edge[rosetta.core.scoring.hbond_bb_sc]
        energies['hbond_sc'] = e_edge[rosetta.core.scoring.hbond_sc]
        energies['weighted_total'] = e_edge.dot(pose.energies().weights())
       
        return energies
    else:
        return None

def get_interaction_residues(pose, target_residue, distance_cutoff, a_e_h_energy_cutoff, rep_cutoff=5, residue_total_energy_cutoff=50):
    '''Get the residues interacting with a target residue.
    The pose should be scored.
    Only return the residues that satisfy the distance and energy cutoffs.
    The only attraction, electrostatic and hydrogne bond energies are considered.
    Return a list of [residue_rosetta_number, residue_pdb_number, residue_type, interaction_energies].
    '''
    interaction_residues = []

    canonical_aa = ['ALA', 'PRO', 'VAL', 'LEU', 'ILE', 'MET',
                    'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS',
                    'LYS', 'ARG', 'HIS', 'ASP', 'GLU', 'ASN',
                    'GLN', 'GLY']


    for i in range(1, pose.size() + 1):
        if i == target_residue: continue

        # Only allow canonical protein residues
        if not pose.residue(i).is_protein(): continue
        if not pose.residue(i).name3() in canonical_aa: continue

        d = min_heavy_atom_distance(pose, i, target_residue)
        
        if d > distance_cutoff: continue
        
        energies = residue_residue_energies(pose, i, target_residue)
        energies['residue_total_energy'] = pose.energies().residue_total_energy(i)

        if energies is None:
            continue

        a_e_h_energies = sum(energies[k] for k in ['fa_atr', 'fa_elec', 'hbond_bb_sc', 'hbond_sc'])

        if a_e_h_energies < a_e_h_energy_cutoff and energies['fa_rep'] < rep_cutoff \
                and energies['residue_total_energy'] < residue_total_energy_cutoff and min_heavy_atom_b_factor(pose, i) > 0:

            ir = [i, pose.pdb_info().pose2pdb(i), pose.residue(i).name3(), energies]
            interaction_residues.append(ir)
            
    return interaction_residues

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

def generate_params_file_for_ligand(ligand_pdb, ligand_name, ligand_id, output_path):
    '''Generate params file for a ligand.'''
    openbabel = '/home/xingjie/Softwares/OpenBabel/openbabel-3.0.0/build/bin/obabel'
    molfile_to_params_script = '/home/xingjie/Softwares/Rosetta/githubRepo/main/source/scripts/python/public/molfile_to_params.py'
    abs_ligand_pdb = os.path.abspath(ligand_pdb)

    # Go to the output path

    cwd = os.getcwd()
    os.chdir(output_path)

    # Run OpenBabel to convert the ligand PDB file to mol2 file

    subprocess.call([openbabel, '-ipdb', abs_ligand_pdb, '-omol2', '-O', 'ligand_{0}.mol2'.format(ligand_id)])

    # Generate params file using the molfile_to_params.py script 

    generate_params_cmd =[
        'python2', molfile_to_params_script,
        '--clobber', # Overwrite existing files
        '--no-pdb',
        '-n', ligand_name,
        '-p', 'ligand_{0}'.format(ligand_id),
        '--conformers-in-one-file',
        '--keep-names',
        'ligand_{0}.mol2'.format(ligand_id)]
   
    subprocess.call(generate_params_cmd)

    # Go back to the current working directory
    
    os.chdir(cwd)
    

def get_binding_sites_for_a_structure(pdb_file, output_path, b_factor_cutoff=60, min_binding_site_size=2, minimize=False, relax=False):
    ''''''
    
    # Load and score the structure
    
    pose = pyrosetta.pose_from_file(pdb_file)
    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose)

    # Relax the pose if required
    
    if minimize:
        minimize_pose(pose)
    if relax:
        fast_relax(pose)

    # Find the ligand residues that have no more than 100 heavy atoms
    
    ligand_residues = find_ligand_residues(pose, size_min=1, size_max=100)

    ligand_id = 0
    ligand_names = set()

    for ligand in ligand_residues:

        # Skip ligands that have high b factors

        if average_heavy_atom_b_factor(pose, ligand) > b_factor_cutoff: continue
      
        # Skip the ligand if the same type of ligand has already been extracted
        # This is to avoid over counting for ligands in structures that have multiple chains

        ligand_name = pose.residue(ligand).name3()
        print(ligand_name)
        if ligand_name in ligand_names:
            continue

        # Skip the ligand if the ligand has bad energy (probably due to clashing)
       
        ligand_total_energy = pose.energies().residue_total_energy(ligand)
        if ligand_total_energy > 10:
            continue

        interaction_residues = get_interaction_residues(pose, ligand, distance_cutoff=5, a_e_h_energy_cutoff=-1)
        binding_pose = extract_subpose(pose, [ligand] + list(i[0] for i in interaction_residues))
        ligand_pose = extract_subpose(pose, [ligand])

        # Get the binding site information
        
        binding_site_info = {}
        binding_site_info['pdb_file'] = pdb_file
        binding_site_info['ligand_name'] = pose.residue(ligand).name3()
        binding_site_info['ligand_n_heavy_atoms'] = len(real_heavy_atoms(pose, ligand))                                                                           
        binding_site_info['ligand_number'] = pose.pdb_info().pose2pdb(ligand)
        binding_site_info['ligand_total_energy'] = ligand_total_energy
        binding_site_info['binding_site_energies'] = interaction_residues

        # Skip the ligand if there are too few binding site residues

        if len(interaction_residues) < min_binding_site_size:
            continue

        # Make the output dir
        
        os.makedirs(output_path, exist_ok=True)
        
        # Write the outputs
        
        ligand_names.add(ligand_name) # record the ligand name
        
        with open(os.path.join(output_path, 'binding_site_info_{0}.json'.format(ligand_id)), 'w') as f:
            json.dump(binding_site_info, f, indent=' ')

        binding_pose.dump_pdb(os.path.join(output_path, 'binding_site_{0}.pdb.gz'.format(ligand_id)))
        ligand_pdb = os.path.join(output_path, 'ligand_{0}.pdb.gz'.format(ligand_id))
        ligand_pose.dump_pdb(ligand_pdb)
       
        # Generate params files for ligands

        generate_params_file_for_ligand(ligand_pdb, ligand_pose.residue(1).name3(), ligand_id, output_path)
        
        ligand_id += 1



if __name__ == '__main__':

    job_id = int(sys.argv[1])

    input_path = 'high_resolution_structures'
    output_path = 'binding_site_database'
    #output_path = 'test_output'

    pyrosetta.init()

    pdb_paths = sorted(os.listdir(input_path))

    pdb_path = pdb_paths[job_id - 1] 
    
    for f in os.listdir(os.path.join(input_path, pdb_path)):
        pdb_code = f[3:7]

        pdb_file = os.path.join(input_path, pdb_path, f)

        output_path_single = os.path.join(output_path, pdb_path, pdb_code)
        
        get_binding_sites_for_a_structure(pdb_file, output_path_single)
