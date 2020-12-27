#!/bin/bash

# Match the 28th ligand binding site in the binding site library to the scaffolds in the scaffolds/test_scaffolds directory. Write the results to the ./test_output directory
# The -n 3 flag tells the script to only keep top 3 protein residues in the binding site that have lowest two-body energy with the ligand
./scripts/match_binding_sites_to_scaffolds.py scaffolds/test_scaffolds ./test_output 9444 -n 3

# The -d flag asks the scripts to dump the matched structures. Using this flag can take a lot of disk space when running libray-to-libray matching.
./scripts/match_binding_sites_to_scaffolds.py scaffolds/test_scaffolds ./test_output 2990 -d -n 3

# Aggregate the matching results to a json file
./scripts/aggregate_matching_results.py test_output aggregated_test_output 
