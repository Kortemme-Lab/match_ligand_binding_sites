#!/usr/bin/env python3
'''Get the high quality crystal structures
from PDB. 
'''

import os
import shutil

from prody import *

if __name__ == '__main__':
  pdb_path = '/wynton/home/database/pdb/remediated/pdb'
  pdb95_list_file = '/wynton/home/database/pdb/remediated/pdb_95.cod'
  output_path = 'high_resolution_structures'

  with open(pdb95_list_file, 'r') as f:
    pdb95_list = f.readlines()

  os.makedirs(output_path, exist_ok=True)

  n_high_res_structures = set()

  for c in pdb95_list:
    d = c[1:3]
    f = 'pdb{0}.ent.gz'.format(c[0:4])
    pdb_file = os.path.join(pdb_path, d, f)
   
    if os.path.exists(pdb_file): 
      print(f)
       
      resolution = parsePDBHeader(pdb_file, 'resolution')
      
      if resolution is None: continue
      
      if resolution < 2:
        os.makedirs(os.path.join(output_path, d), exist_ok=True)
        #shutil.copy(pdb_file, os.path.join(output_path, d, f))
        
        n_high_res_structures.add(f)

  print('Found {0} high resolution structures.'.format(len(n_high_res_structures)))
