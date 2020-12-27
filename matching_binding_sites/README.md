This directory contains scripts to match ligand binding site libraries to scaffold libraries.

In order to run the matching scripts, the `../binding_sites_from_pdb/binding_site_database` folder and the `../binding_sites_from_pdb/binding_site_database_summary.json` file must be generated in advance.

The scaffold libraries are stored in the `scaffolds` directory.

The fast matching and standard rosetta matching protocols are in the `./scripts/fast_match.py` and the `./scripts/rosetta_standard_match.py` files respectively. The `./scripts/match_binding_sites_to_scaffolds.py` script matches the binding site database to the user defined scaffold libraries.

Here is the command to run an example matching job
```
./scripts/run_test_match.sh
```
You can open the file with a text editor to see what this script does.

Matching all the binding sites to a library of scaffolds is computationally intense and requires a cluster. Refer to the `./scripts/user_scripts_xingjie/submit_matching_5tpj_no_c_term.sh` for submitting matching jobs to the UCSF Wynton cluster.
