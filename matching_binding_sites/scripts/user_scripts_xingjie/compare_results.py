#!/usr/bin/env python3
'''Compare two aggregated results.
'''

import json

import numpy as np
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt


def aggregated_results_to_dict(aggregated_results):
    '''Convert a list of aggregated_results into a dictionary.'''
    aggregated_results_dict = {}

    for i in aggregated_results:
        if 'match_info_file' in i:
            key = '/'.join(i['match_info_file'].split('/')[-4:-1])
            aggregated_results_dict[key] = i

    return aggregated_results_dict

def get_over_lapping_sites(aggregated_results_file1, aggregated_results_file2, ligand_size=None, N_total_sites=None):
    with open(aggregated_results_file1, 'r') as f:
        aggregated_results1 = json.load(f)
   
    with open(aggregated_results_file2, 'r') as f:
        aggregated_results2 = json.load(f)

    aggregated_results_dict1 = aggregated_results_to_dict(aggregated_results1)
    aggregated_results_dict2 = aggregated_results_to_dict(aggregated_results2)

    keys1 = set(aggregated_results_dict1)   
    keys2 = set(aggregated_results_dict2)

    #keys = keys1.intersection(keys2)
    keys = keys1.union(keys2)

    overlap_fast_match = {'one_only':0, 'overlap':0, 'two_only':0, 'none':0}     
    overlap_rosetta_match = {'one_only':0, 'overlap':0, 'two_only':0, 'none':0}     

    for k in keys:
        if not k in keys1:
            aggregated_results_dict1[k] = {'num_fast_match_succeed':0, 'num_rosetta_match_succeed':0, 'ligand_n_heavy_atoms':0}
        if not k in keys2:
            aggregated_results_dict2[k] = {'num_fast_match_succeed':0, 'num_rosetta_match_succeed':0, 'ligand_n_heavy_atoms':0}

        if not (ligand_size is None) and aggregated_results_dict1[k]['ligand_n_heavy_atoms'] != ligand_size:
            continue

        if aggregated_results_dict1[k]['num_fast_match_succeed'] > 0 and aggregated_results_dict2[k]['num_fast_match_succeed'] < 1:
            overlap_fast_match['one_only'] += 1
        elif aggregated_results_dict1[k]['num_fast_match_succeed'] > 0 and aggregated_results_dict2[k]['num_fast_match_succeed'] > 0:
            overlap_fast_match['overlap'] += 1
        elif aggregated_results_dict1[k]['num_fast_match_succeed'] < 1 and aggregated_results_dict2[k]['num_fast_match_succeed'] > 0:
            overlap_fast_match['two_only'] += 1
        else:
            overlap_fast_match['none'] += 1

        if aggregated_results_dict1[k]['num_rosetta_match_succeed'] > 0 and aggregated_results_dict2[k]['num_rosetta_match_succeed'] < 1:
            overlap_rosetta_match['one_only'] += 1
        elif aggregated_results_dict1[k]['num_rosetta_match_succeed'] > 0 and aggregated_results_dict2[k]['num_rosetta_match_succeed'] > 0:
            overlap_rosetta_match['overlap'] += 1
        elif aggregated_results_dict1[k]['num_rosetta_match_succeed'] < 1 and aggregated_results_dict2[k]['num_rosetta_match_succeed'] > 0:
            overlap_rosetta_match['two_only'] += 1
        else:
            overlap_rosetta_match['none'] += 1

    if not (N_total_sites is None):
        overlap_fast_match['none'] = N_total_sites - overlap_fast_match['overlap'] - overlap_fast_match['one_only'] - overlap_fast_match['two_only']
        overlap_rosetta_match['none'] = N_total_sites - overlap_rosetta_match['overlap'] - overlap_rosetta_match['one_only'] - overlap_rosetta_match['two_only']

    # Run chi-square test 

    fast_match_table = np.array([[overlap_fast_match['overlap'], overlap_fast_match['one_only']], 
        [overlap_fast_match['two_only'], overlap_fast_match['none']]])

    rosetta_match_table = np.array([[overlap_rosetta_match['overlap'], overlap_rosetta_match['one_only']], 
        [overlap_rosetta_match['two_only'], overlap_rosetta_match['none']]])

    chi_square_test_fast = chi2_contingency(fast_match_table)
    chi_square_test_rosetta = chi2_contingency(rosetta_match_table)

    print('{0} VS {1}'.format(aggregated_results_file1, aggregated_results_file2))
    print(overlap_fast_match)
    print(chi_square_test_fast)
    print(overlap_rosetta_match, sum(sum(rosetta_match_table)))
    print(chi_square_test_rosetta)

    one = rosetta_match_table[0][0] + rosetta_match_table[0][1]
    two = rosetta_match_table[0][0] + rosetta_match_table[1][0]
    total = sum(sum(rosetta_match_table))
    print('The expected number of overlaps is {0}'.format(one / total * two))
    print('')

def plot_over_lapping_vs_lig_size(aggregated_results_file1, aggregated_results_file2, output_file=None, fast_match=False):
    # Load the results
    
    with open(aggregated_results_file1, 'r') as f:
        aggregated_results1 = json.load(f)
   
    with open(aggregated_results_file2, 'r') as f:
        aggregated_results2 = json.load(f)

    aggregated_results_dict1 = aggregated_results_to_dict(aggregated_results1)
    aggregated_results_dict2 = aggregated_results_to_dict(aggregated_results2)

    keys1 = set(aggregated_results_dict1)   
    keys2 = set(aggregated_results_dict2)

    keys = keys1.intersection(keys2)

    dict_overlapping_results = {}

    for k in keys:
        lig_size = aggregated_results_dict1[k]['ligand_n_heavy_atoms']

        if not (lig_size in dict_overlapping_results):
            dict_overlapping_results[lig_size] = {'one_only':0, 'overlap':0, 'two_only':0, 'none':0}     

        if fast_match:
            if aggregated_results_dict1[k]['num_fast_match_succeed'] > 0 and aggregated_results_dict2[k]['num_fast_match_succeed'] < 1:
                dict_overlapping_results[lig_size]['one_only'] += 1
            elif aggregated_results_dict1[k]['num_fast_match_succeed'] > 0 and aggregated_results_dict2[k]['num_fast_match_succeed'] > 0:
                dict_overlapping_results[lig_size]['overlap'] += 1
            elif aggregated_results_dict1[k]['num_fast_match_succeed'] < 1 and aggregated_results_dict2[k]['num_fast_match_succeed'] > 0:
                dict_overlapping_results[lig_size]['two_only'] += 1
            else:
                dict_overlapping_results[lig_size]['none'] += 1
        else:
            if aggregated_results_dict1[k]['num_rosetta_match_succeed'] > 0 and aggregated_results_dict2[k]['num_rosetta_match_succeed'] < 1:
                dict_overlapping_results[lig_size]['one_only'] += 1
            elif aggregated_results_dict1[k]['num_rosetta_match_succeed'] > 0 and aggregated_results_dict2[k]['num_rosetta_match_succeed'] > 0:
                dict_overlapping_results[lig_size]['overlap'] += 1
            elif aggregated_results_dict1[k]['num_rosetta_match_succeed'] < 1 and aggregated_results_dict2[k]['num_rosetta_match_succeed'] > 0:
                dict_overlapping_results[lig_size]['two_only'] += 1
            else:
                dict_overlapping_results[lig_size]['none'] += 1

    # Calculate the difference between observed and expected overlaps  

    dict_overlapping_freqs = {}

    for lig_size in dict_overlapping_results.keys():
        overlap = dict_overlapping_results[lig_size]['overlap']
        one_all = dict_overlapping_results[lig_size]['one_only'] + overlap
        two_all = dict_overlapping_results[lig_size]['two_only'] + overlap
        none = dict_overlapping_results[lig_size]['none']
        count_all = overlap + one_all + two_all + none

        try:
            chi_square_test = chi2_contingency(np.array([[overlap, one_all - overlap], [two_all - overlap, none]]))
            freq_one = one_all / count_all
            freq_two = two_all / count_all
            freq_overlap = overlap / count_all
            expected_freq = freq_one * freq_two
       
            dict_overlapping_freqs[lig_size] = {'freq_overlap': freq_overlap, 'expected_freq': expected_freq, 
                    'count_overlap': overlap, 'p': chi_square_test[1]}

        except:
            continue

    # Plot

    X = []
    Y = []
    counts = []
    S = []
    logP = []

    for x in sorted(dict_overlapping_freqs.keys()):
        print(x, dict_overlapping_freqs[x])

    for x in dict_overlapping_freqs.keys():
        X.append(x)
        Y.append(dict_overlapping_freqs[x]['freq_overlap'] - dict_overlapping_freqs[x]['expected_freq'])
        counts.append(dict_overlapping_freqs[x]['count_overlap'])
        S.append(10 * dict_overlapping_freqs[x]['count_overlap'] ** 0.5)
        if dict_overlapping_freqs[x]['p'] == 0:
            logP.append(-9999)
        else:
            logP.append(np.log10(dict_overlapping_freqs[x]['p']))

    fig, ax = plt.subplots()
    scatter = ax.scatter(X, Y, s=S, c=logP, cmap='coolwarm', vmin=-30, vmax=0)
    #scatter = ax.scatter(X, Y, s=S, c=logP, cmap='coolwarm')
    cb = fig.colorbar(scatter)
    cb.set_label('log(P-value)', fontsize=20)

    handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6, num=[1, 10, 100, 500, 1000], func=lambda x:(x / 10)**2)
    legend = ax.legend(handles, labels, title="N overlaps", fontsize=12)
    
    plt.xlabel('N ligand heavy atoms', fontsize=20)
    plt.ylabel('$F_{observed} - F_{expected}$', fontsize=20)
    
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.subplots_adjust(left=0.15, bottom=0.15)
    
    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)
    
    plt.clf()




if __name__ == '__main__':

    aggregated_results_files = {
            '3res_2lv8':'matching_data/aggregated_matching_3res_site_to_2lv8_two_LHL.json', 
            '3res_5tpj':'matching_data/aggregated_matching_3res_site_to_5tpj.json', 
            '3res_native_rossmann':'matching_data/aggregated_matching_3res_site_to_native_rossmann.json', 
            '3res_native_ntf2':'matching_data/aggregated_matching_3res_site_to_native_ntf2.json', 
            
            '2lv8':'matching_data/aggregated_matching_to_2lv8_two_LHL.json', 
            '5tpj':'matching_data/aggregated_matching_to_5tpj.json', 
            'native_rossmann':'matching_data/aggregated_matching_to_native_rossmann.json', 
            'native_ntf2':'matching_data/aggregated_matching_to_native_ntf2.json', 
            }

    get_over_lapping_sites(aggregated_results_files['3res_native_ntf2'], aggregated_results_files['3res_5tpj'], N_total_sites=20102)
    get_over_lapping_sites(aggregated_results_files['3res_2lv8'], aggregated_results_files['3res_5tpj'], N_total_sites=20102)
    get_over_lapping_sites(aggregated_results_files['3res_2lv8'], aggregated_results_files['3res_native_rossmann'], N_total_sites=20102)
    get_over_lapping_sites(aggregated_results_files['3res_native_rossmann'], aggregated_results_files['3res_native_ntf2'], N_total_sites=20102)

    #plot_over_lapping_vs_lig_size(aggregated_results_files['3res_native_ntf2'], aggregated_results_files['3res_5tpj'], fast_match=False, 
    #        output_file='./plots/overlap_vs_lig_size_3res_native_ntf2_and_3res_5tpj.svg')
    #plot_over_lapping_vs_lig_size(aggregated_results_files['3res_2lv8'], aggregated_results_files['3res_5tpj'], fast_match=False, 
    #        output_file='./plots/overlap_vs_lig_size_3res_2lv8_and_3res_5tpj.svg')
