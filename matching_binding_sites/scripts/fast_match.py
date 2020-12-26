#!/usr/bin/env python3
'''A method to quickly match a active site to a 
protein scaffold.
'''

import os
import json

import numpy as np

import pyrosetta
from pyrosetta import rosetta


def xyzV_to_np_array(xyz):
    return np.array([xyz.x, xyz.y, xyz.z])

def np_array_to_xyzV(a):
    return rosetta.numeric.xyzVector_double_t(a[0], a[1], a[2])

def np_array_to_xyzM(a):
    return rosetta.numeric.xyzMatrix_double_t.rows(
            a[0][0], a[0][1], a[0][2],
            a[1][0], a[1][1], a[1][2],
            a[2][0], a[2][1], a[2][2])

def get_backbone_points(pose, residues):
    '''Get backbone points for residues in a pose.'''
    points = []

    for res in residues:
        for atom in ['N', 'CA', 'C']:
            points.append(xyzV_to_np_array(pose.residue(res).xyz(atom)))

    return points

def get_superimpose_transformation(P1, P2):
    '''Get the superimpose transformation that transfoms a list of
    points P1 to another list of points P2.'''
    if len(P1) != len(P2):
        raise Exception("Sets to be superimposed must have same number of points.")

    com1 = np.mean(P1, axis=0)
    com2 = np.mean(P2, axis=0)

    R = np.dot(np.transpose(np.array(P1) - com1), np.array(P2) - com2)
    V, S, W = np.linalg.svd(R)

    if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
        V[:, -1] = -V[:, -1]

    M = np.transpose(np.array(np.dot(V, W)))

    return M, com2 - np.dot(M, com1)

def superimpose_poses_by_residues(pose_source, residues_source, pose_target, residues_target):
    '''Superimpose residues in a source pose into residues in a target pose.
    Only backbone atoms are used for the superimposition.
    '''
    assert(len(residues_source) == len(residues_target))

    # Get the points to be superimposed

    points_source = get_backbone_points(pose_source, residues_source)
    points_target = get_backbone_points(pose_target, residues_target)

    # Get the rigid body transformation

    M, t = get_superimpose_transformation(points_source, points_target)

    # Transform the source pose

    pose_source.apply_transform_Rx_plus_v(np_array_to_xyzM(M), 
            np_array_to_xyzV(t))


def calc_backbone_RMSD(pose1, residues1, pose2, residues2):
    '''Calculate backbone RMSD between two poses for specific positions.'''
    assert(len(residues1) == len(residues2))

    def RMSD(points1, poinsts2):
        '''Calcualte RMSD between two lists of numpy points.'''
        diff = [points1[i] - poinsts2[i] for i in range(len(points1))]
        return np.sqrt(sum(np.dot(d, d) for d in diff) / len(diff))

    points1 = get_backbone_points(pose1, residues1)
    points2 = get_backbone_points(pose2, residues2)

    return RMSD(points1, points2)

def calc_res_direction_cos(pose1, residues1, pose2, residues2):
    '''Calculate the cos between CA-CB directions.
    And the cos between N-C directions.
    '''
    coses = []

    for i in range(len(residues1)):
        res1 = pose1.residue(residues1[i])
        res2 = pose2.residue(residues2[i])
        if 'GLY' == res1.name3() or 'GLY' == res2.name3():
            coses.append(1)

        else:
            ca1 = res1.xyz('CA')
            cb1 = res1.xyz('CB')
            ca2 = res2.xyz('CA')
            cb2 = res2.xyz('CB')

            coses.append((cb1 - ca1).normalize().dot((cb2 - ca1).normalize()))

        n1 = res1.xyz('N')
        c1 = res1.xyz('C')
        n2 = res2.xyz('N')
        c2 = res2.xyz('C')
        
        coses.append((c1 - n1).normalize().dot((c2 - n2).normalize()))

    return coses


def match_anchor_position(target_pose, target_anchor_seqpos, movable_pose, movable_anchor_seqpos):
    '''Match the movable pose to the anchor
    position of the target pose.
    '''
    res_t = target_pose.residue(target_anchor_seqpos)
    res_m = movable_pose.residue(movable_anchor_seqpos)

    current_frame = rosetta.numeric.xyzTransform_double_t(res_m.xyz('N'), res_m.xyz('CA'), res_m.xyz('C'))
    inv_cf = current_frame.inverse()
    target_frame = rosetta.numeric.xyzTransform_double_t(res_t.xyz('N'), res_t.xyz('CA'), res_t.xyz('C'))

    R = target_frame.R * inv_cf.R
    v = target_frame.R * inv_cf.t + target_frame.t

    movable_pose.apply_transform_Rx_plus_v(R, v) 

def num_clash_between_ligand_pose_and_scaffold_bb(ligand_pose, scaffold_pose, scaffold_pose_matched_residues, scale_factor=0.36):
    '''Return the number of clashes between the 
    ligand pose and the scaffold backbone.
    The matched residues in the scaffold pose are excluded
    for calculation.
    Only heavy atoms are considered.
    '''
    num_clashes = 0

    # Find the residues in the scaffold to calculate clashes

    scaffold_residues_to_calc = []
    for seqpos in range(1, scaffold_pose.size() + 1):
        if not scaffold_pose.residue(seqpos).is_protein(): continue
        if seqpos in scaffold_pose_matched_residues: continue

        scaffold_residues_to_calc.append(seqpos)

    # Calculate clashes

    # Iterate through scaffold backbone heavy atoms

    for seqpos1 in scaffold_residues_to_calc:
        res1 = scaffold_pose.residue(seqpos1)
        
        for i in range(1, res1.last_backbone_atom() + 1):
            vi = res1.xyz(i)
            ri = res1.atom_type(i).lj_radius()

            # Iterate through ligand pose heavy atoms.
            # The protein backbone atoms are ignored

            for seqpos2 in range(1, ligand_pose.size() + 1):
                res2 = ligand_pose.residue(seqpos2)

                if res2.is_ligand():
                    a_start = 1
                else:
                    a_start = res2.first_sidechain_atom()

                a_stop = res2.nheavyatoms()
            
                for j in range(a_start, a_stop + 1):
                    vj = res2.xyz(j)
                    rj = res2.atom_type(j).lj_radius()

                    # Check clash

                    if (vi - vj).length_squared() < scale_factor * ((ri + rj) ** 2):
                        num_clashes += 1

    return num_clashes
    



def fast_match(site_pose, site_protein_residues, scaffold_pose, scaffold_protein_residues, output_path, 
        cut_off_max_dist=2, cut_off_max_rmsd=1, cut_off_min_direction_cos=0.7, dump_matches=False, dump_trajectory=False):
    '''Match a binding site pose to scaffold_pose.
    Dump the matched results to the output_path.
    Return the number of matches.
    '''
    # Use the first protein residue as the anchoring residue

    anchor_residue = site_protein_residues[0]

    num_matches = 0
    matched_positions = []

    # Iterate through scaffold stubs

    for j in scaffold_protein_residues:

        # Anchor the site pose to the scaffold

        match_anchor_position(scaffold_pose, j, site_pose, anchor_residue)
    
        # Find the best matched residues for each site residue

        best_matched_residues = []
        best_matched_max_dists = []

        for ii in site_protein_residues:
            best_max_dist = float('inf')

            for jj in scaffold_protein_residues:

                dist_N = site_pose.residue(ii).xyz('N').distance(scaffold_pose.residue(jj).xyz('N'))
                dist_CA = site_pose.residue(ii).xyz('CA').distance(scaffold_pose.residue(jj).xyz('CA'))
                dist_C = site_pose.residue(ii).xyz('C').distance(scaffold_pose.residue(jj).xyz('C'))

                max_dist = max([dist_N, dist_CA, dist_C])
       
                if max_dist < best_max_dist:
                    best_max_dist = max_dist
                    best_residue = jj

            best_matched_residues.append(best_residue)
            best_matched_max_dists.append(best_max_dist)

        # Discard the match if multiple site residues are matched to one scaffold residue

        if len(set(best_matched_residues)) != len(best_matched_residues):
            continue

        # Discard the match if the worst distance of best matches are greater than the cut_off_max_dist

        if max(best_matched_max_dists) > cut_off_max_dist:
            continue

        if dump_trajectory:
            os.makedirs(output_path, exist_ok=True)
            site_pose.dump_pdb(os.path.join(output_path, 'match_{0}_step1.pdb.gz'.format(num_matches)))

        # Discard the match if the RMSD of the aligned site is greater than the cut_off_max_rmsd

        superimpose_poses_by_residues(site_pose, site_protein_residues, scaffold_pose, best_matched_residues)
        bb_rmsd = calc_backbone_RMSD(site_pose, site_protein_residues, scaffold_pose, best_matched_residues)

        if bb_rmsd > cut_off_max_rmsd:
            continue
        
        # Discard the match if the side chain directions don't match

        res_dir_coses = calc_res_direction_cos(site_pose, site_protein_residues, scaffold_pose, best_matched_residues)

        if min(res_dir_coses) < cut_off_min_direction_cos:
            continue

        # Check clashing between the ligand and the backbone atoms of the scaffold pose

        num_clashes = num_clash_between_ligand_pose_and_scaffold_bb(site_pose, scaffold_pose, best_matched_residues)
        if num_clashes > 0:
            continue

        if dump_matches:

            # Make the folder for output

            os.makedirs(output_path, exist_ok=True)

            # Dump the match

            scaffold_pose.dump_pdb(os.path.join(output_path, 'scaffold.pdb.gz'))
            site_pose.dump_pdb(os.path.join(output_path, 'match_{0}.pdb.gz'.format(num_matches)))
            
            match_info = {'binding_site_residues': site_protein_residues,
                          'matched_residues': best_matched_residues}

            with open(os.path.join(output_path, 'match_info_{0}.json'.format(num_matches)), 'w') as f:
                json.dump(match_info, f)

        matched_positions.append(best_matched_residues)
        num_matches += 1

    return num_matches, matched_positions


if __name__ == '__main__':

    pyrosetta.init()

    site_pose = pyrosetta.pose_from_file('./test_inputs/test_site.pdb')
    scaffold_pose = pyrosetta.pose_from_file('./test_inputs/2src_scaffold.pdb')

    site_protein_residues = [i for i in range(1, site_pose.size() + 1) if site_pose.residue(i).is_protein()] 
    scaffold_protein_residues = [i for i in range(1, scaffold_pose.size() + 1) if scaffold_pose.residue(i).is_protein()]

    print(fast_match(site_pose, site_protein_residues, scaffold_pose, scaffold_protein_residues, 'test_output'))
