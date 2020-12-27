This directory contains scripts to extract binding sites from the PDB95 database.

If you have the PDB95 database, you can extract the binding sites using the `./scripts/get_high_quality_structures.py` and the `./scripts/extract_ligand_binding_sites.py` scripts.

Use can also use a precompiled binding site database for matching. Concatenate the two parts of the compressed binding site database
```
cat binding_site_database.tar.gz.part* > binding_site_database.tar.gz
```
Decompress the binding site database
```
tar -xzf binding_site_database.tar.gz
```

Finally, you need to generate a binding site summary file for matching jobs by running
```
./scripts/aggregate_binding_site_informations.py
```
