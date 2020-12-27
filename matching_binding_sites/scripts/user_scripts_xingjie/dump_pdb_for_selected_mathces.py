#!/usr/bin/env python3

import os
import shutil
import subprocess


def dump_pdb_for_selected_match(scaffold_pdb, output_path, job_id, binding_site_size=None):
    os.makedirs(output_path, exist_ok=True)
    new_scaffold_path = os.path.join(output_path, 'scaffold')
    
    os.makedirs(new_scaffold_path, exist_ok=True) 

    shutil.copy(scaffold_pdb, new_scaffold_path)
   
    if None is binding_site_size:
        subprocess.call(['./scripts/match_binding_sites_to_scaffolds.py',
                new_scaffold_path,
                output_path,
                str(job_id),
                '-d'
                ])
    else:
        subprocess.call(['./scripts/match_binding_sites_to_scaffolds.py',
                new_scaffold_path,
                output_path,
                str(job_id),
                '-d',
                '-n', str(binding_site_size)
                ])



if __name__ == '__main__':

    
    #dump_pdb_for_selected_match('scaffolds/designable_5tpj_no_c_term_1000_random/185196.pdb.gz', 'test_output', 1644, binding_site_size=3) # 3 res match
    #dump_pdb_for_selected_match('scaffolds/designable_5tpj_no_c_term_1000_random/238211.pdb.gz', 'test_output', 6891, binding_site_size=3) # 3 res match
    #dump_pdb_for_selected_match('scaffolds/designable_5tpj_no_c_term_1000_random/168841.pdb.gz', 'test_output', 17649, binding_site_size=3) # 3 res match
    #dump_pdb_for_selected_match('scaffolds/designable_5tpj_no_c_term_1000_random/224740.pdb.gz', 'test_output', 6447, binding_site_size=3) # 3 res match
    #dump_pdb_for_selected_match('scaffolds/designable_5tpj_no_c_term_1000_random/304367.pdb.gz', 'test_output', 6797, binding_site_size=3) # 3 res match
    
    #dump_pdb_for_selected_match('scaffolds/designable_2lv8_two_LHL_1000_random/model_301404.pdb.gz', 'test_output', 16307) # 5 res match
    
    #dump_pdb_for_selected_match('scaffolds/native_ntf2_scaffolds/3f9sB00.pdb.gz', 'test_output', 608, binding_site_size=3) # 5 res match
    
    #dump_pdb_for_selected_match('scaffolds/designable_5tpj_no_c_term_1000_random/388152.pdb.gz', 'test_output', 20747, binding_site_size=3) # 3 res large ligand match
    dump_pdb_for_selected_match('scaffolds/designable_5tpj_no_c_term_1000_random/366802.pdb.gz', 'test_output', 2990, binding_site_size=3) # 3 res match
