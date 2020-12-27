#!/usr/bin/env python3
'''Plot the statistics of a matching result.
Usage:
    ./plot_matching_statistics.py
'''

import sys
import json

import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt


def print_average_successful_running_times(aggregated_results):
    fast_match_times = []
    rosetta_match_times = []
    
    for match_info in aggregated_results:
        if 'num_fast_match_succeed' in match_info:
            if match_info['num_fast_match_succeed'] > 0:
                fast_match_times += [match_info['average_success_fast_match_time']] * match_info['num_fast_match_succeed']

                if match_info['num_rosetta_match_succeed'] > 0:
                    rosetta_match_times += [match_info['average_success_rosetta_match_time']] * match_info['num_rosetta_match_succeed']
        
    print('average_fast_match_time = {0}, average_rosetta_match_time = {1}'.format(np.mean(fast_match_times), np.mean(rosetta_match_times)))


def parse_fast_match_results(aggregated_results):

    num_success = 0
    num_failure = 0
    
    for match_info in aggregated_results:
        if 'num_fast_match_succeed' in match_info:
            if match_info['num_fast_match_succeed'] > 0:
                num_success += 1
            else:
                num_failure += 1

    return num_success, num_failure

def parse_rosetta_match_results(aggregated_results):

    num_success = 0
    num_failure = 0
    
    for match_info in aggregated_results:
        if 'num_rosetta_match_succeed' in match_info:
            if match_info['num_fast_match_succeed'] > 0:
                if match_info['num_rosetta_match_succeed'] > 0:
                    num_success += 1
                elif match_info['num_rosetta_match_succeed'] == 0:
                    num_failure += 1

    return num_success, num_failure

def print_number_of_matches(dict_aggregated_results):
    
    for k in dict_aggregated_results.keys():
        n_fast_match_success = parse_fast_match_results(dict_aggregated_results[k])[0]
        n_rosetta_match_success = parse_rosetta_match_results(dict_aggregated_results[k])[0]
        
        print('{0}, n_fast_match_success={1}, n_rosetta_match_success={2}'.format(k, n_fast_match_success, n_rosetta_match_success))



def plot_match_results_vs_lig_and_site_size(aggregated_results):
    # Initialize the matrices
    
    m_rosetta_success = []
    m_fast_match_success = []
    m_total = []
    m_fast_success_rate =[]
    m_rosetta_success_rate = []

    max_lig_size = 70
    max_site_size = 35

    for i in range(max_lig_size):
        m_rosetta_success.append([0] * (max_site_size))
        m_fast_match_success.append([0] * (max_site_size))
        m_total.append([0] * (max_site_size))
        m_fast_success_rate.append([0] * (max_site_size))
        m_rosetta_success_rate.append([0] * (max_site_size))

    # Parse the results

    for match_info in aggregated_results:
        
        # Skip the bad entries

        if 'ligand_n_heavy_atoms' in match_info and 'binding_site_residues' in match_info:

            lig_size = match_info['ligand_n_heavy_atoms']
            site_size = len(match_info['binding_site_residues'])
            
            if lig_size <= max_lig_size and 0 < site_size <= max_site_size:

                if 'num_fast_match_succeed' in match_info and match_info['num_fast_match_succeed'] > -1:
                    m_total[lig_size - 1][site_size - 1] += 1
 
                    if match_info['num_fast_match_succeed'] > 0:
                        m_fast_match_success[lig_size - 1][site_size - 1] += 1
                        
                        if match_info['num_rosetta_match_succeed'] > 0:
                            m_rosetta_success[lig_size - 1][site_size - 1] += 1

    # Calculate the success rates

    for i in range(max_lig_size):
        for j in range(max_site_size):
            if m_total[i][j] > 0:
                m_fast_success_rate[i][j] = m_fast_match_success[i][j] / m_total[i][j]
            else:
                m_fast_success_rate[i][j] = np.nan

            if m_fast_match_success[i][j] > 0:
                m_rosetta_success_rate[i][j] = m_rosetta_success[i][j] / m_fast_match_success[i][j]
            else:
                m_rosetta_success_rate[i][j] = np.nan

    # Make the plots
    
    fig, ax = plt.subplots()

    #cax = ax.imshow(np.transpose(m_rosetta_success), origin='lower', interpolation='nearest', cmap='cividis', vmin=0, vmax=30, 
    #        extent=(0.5, max_lig_size + 0.5, 0.5, max_site_size + 0.5), aspect='auto')
    cax = ax.imshow(np.transpose(m_fast_success_rate), origin='lower', interpolation='nearest', cmap='cividis', vmin=0, vmax=1, 
            extent=(0.5, max_lig_size + 0.5, 0.5, max_site_size + 0.5), aspect='auto')
    #cax = ax.imshow(np.transpose(m_rosetta_success_rate), origin='lower', interpolation='nearest', cmap='cividis', vmin=0, vmax=1, 
    #        extent=(0.5, max_lig_size + 0.5, 0.5, max_site_size + 0.5), aspect='auto')

    plt.show()



def plot_match_results_vs_ligand_size(aggregated_results, x_max=None):
    results_by_ligand_size = {}
    
    for match_info in aggregated_results:
        if 'num_rosetta_match_succeed' in match_info:
            ligand_size = match_info['ligand_n_heavy_atoms']
            
            if not ligand_size in results_by_ligand_size:
                results_by_ligand_size[ligand_size] = {
                        'fast_match_failed': 0,
                        'fast_match_succeed': 0,
                        'rosetta_match_failed': 0,
                        'rosetta_match_succeed': 0
                        }

            if match_info['num_fast_match_succeed'] > 0:
                results_by_ligand_size[ligand_size]['fast_match_succeed'] += 1
                
                if match_info['num_rosetta_match_succeed'] > 0:
                    results_by_ligand_size[ligand_size]['rosetta_match_succeed'] += 1
                elif match_info['num_rosetta_match_succeed'] == 0:
                    results_by_ligand_size[ligand_size]['rosetta_match_failed'] += 1

            else:
                results_by_ligand_size[ligand_size]['fast_match_failed'] += 1

    keys = sorted(results_by_ligand_size.keys())

    # Append the missing keys

    max_key = max(results_by_ligand_size.keys())

    for i in range(1, max_key):
        if not (i in results_by_ligand_size):
           results_by_ligand_size[i] = {
                   'fast_match_failed': 0,
                   'fast_match_succeed': 0,
                   'rosetta_match_failed': 0,
                   'rosetta_match_succeed': 0
                   }

    keys = sorted(results_by_ligand_size.keys())

    if not (x_max is None):
        keys = [k for k in keys if k <= x_max]


    no_match_values = []
    only_fast_values = []
    rosetta_match_values = []
    fast_success_rate = []
    rosetta_success_rate = []
    x_ticks = []

    for k in keys:
        print(k, results_by_ligand_size[k])
        no_match_values.append(results_by_ligand_size[k]['fast_match_failed']) 
        only_fast_values.append(results_by_ligand_size[k]['fast_match_succeed'] - results_by_ligand_size[k]['rosetta_match_succeed']) 
        rosetta_match_values.append(results_by_ligand_size[k]['rosetta_match_succeed']) 
        x_ticks.append(k)
        
        n_total = no_match_values[-1] + only_fast_values[-1] + rosetta_match_values[-1]
        if n_total > 0:
            fast_success_rate.append((only_fast_values[-1] + rosetta_match_values[-1]) / n_total)
            rosetta_success_rate.append((rosetta_match_values[-1]) / n_total)
        else:
            fast_success_rate.append(0)
            rosetta_success_rate.append(0)


    #make_stacked_bar_plots(no_match_values, only_fast_values, rosetta_match_values, 'Binding site size')

    positions = list(range(len(no_match_values)))
    plt.bar(positions, fast_success_rate)
    plt.bar(positions, rosetta_success_rate)
    
    plt.xticks(positions, x_ticks)
    plt.show()


def plot_match_results_vs_site_size(aggregated_results, x_max=None, output_file=None):
    results_by_site_size = {}
    
    for match_info in aggregated_results:
        if 'num_rosetta_match_succeed' in match_info:
            binding_site_size = len(match_info['binding_site_residues'])
            
            if not binding_site_size in results_by_site_size:
                results_by_site_size[binding_site_size] = {
                        'fast_match_failed': 0,
                        'fast_match_succeed': 0,
                        'rosetta_match_failed': 0,
                        'rosetta_match_succeed': 0
                        }

            if match_info['num_fast_match_succeed'] > 0:
                results_by_site_size[binding_site_size]['fast_match_succeed'] += 1
                
                if match_info['num_rosetta_match_succeed'] > 0:
                    results_by_site_size[binding_site_size]['rosetta_match_succeed'] += 1
                elif match_info['num_rosetta_match_succeed'] == 0:
                    results_by_site_size[binding_site_size]['rosetta_match_failed'] += 1

            else:
                results_by_site_size[binding_site_size]['fast_match_failed'] += 1

    # Append the missing keys

    max_key = max(results_by_site_size.keys())

    for i in range(2, max_key):
        if not (i in results_by_site_size):
           results_by_site_size[i] = {
                   'fast_match_failed': 0,
                   'fast_match_succeed': 0,
                   'rosetta_match_failed': 0,
                   'rosetta_match_succeed': 0
                   }

    keys = sorted(results_by_site_size.keys())

    if not (x_max is None):
        keys = [k for k in keys if k <= x_max]

    no_match_values = []
    only_fast_values = []
    rosetta_match_values = []
    x_ticks = []

    N_fast_total_matches = 0
    N_rosetta_total_matches = 0
    for k in keys:
        
        print(k, results_by_site_size[k], 'success_rate = {0}%'.format( 
                results_by_site_size[k]['rosetta_match_succeed'] / (results_by_site_size[k]['fast_match_succeed'] + results_by_site_size[k]['fast_match_failed']) * 100))
        no_match_values.append(results_by_site_size[k]['fast_match_failed']) 
        only_fast_values.append(results_by_site_size[k]['fast_match_succeed'] - results_by_site_size[k]['rosetta_match_succeed']) 
        rosetta_match_values.append(results_by_site_size[k]['rosetta_match_succeed']) 
       
        N_fast_total_matches += results_by_site_size[k]['fast_match_succeed']
        N_rosetta_total_matches += results_by_site_size[k]['rosetta_match_succeed']

        x_ticks.append(k)

    print('N_fast_total_matches = {0}, N_rosetta_total_matches = {1}'.format(N_fast_total_matches, N_rosetta_total_matches))

    make_stacked_bar_plots(no_match_values, only_fast_values, rosetta_match_values, 'Binding site size', 
            positions=keys, output_file=output_file)



def make_stacked_bar_plots(no_match_values, only_fast_values, rosetta_match_values, xlabel, positions=None, x_ticks=None, output_file=None):
    no_match_values = np.array(no_match_values)
    only_fast_values = np.array(only_fast_values)
    rosetta_match_values = np.array(rosetta_match_values)
  
    if None is positions:
        positions = list(range(len(no_match_values)))

    rects1 = plt.bar(positions, no_match_values, bottom=rosetta_match_values+only_fast_values)
    rects2 = plt.bar(positions, only_fast_values, bottom=rosetta_match_values)
    rects3 = plt.bar(positions, rosetta_match_values)
  
    if not (x_ticks is None):
        plt.xticks(positions, x_ticks)

    plt.legend([rects1, rects2, rects3], ['No match', 'Only fast match', 'Success Rosetta match'], fontsize=15)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel('Number of binding sites', fontsize=20)

    plt.subplots_adjust(left=0.15, bottom=0.15)

    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)
    
    plt.clf()

def plot_success_for_different_data_sets(list_aggregated_results, x_ticks):
   
    no_match_values = []
    only_fast_values = []
    rosetta_match_values = []

    for aggregated_results in list_aggregated_results:
        fast_match_succeed, fast_match_failed = parse_fast_match_results(aggregated_results)
        rosetta_match_succeed, rosetta_match_failed = parse_rosetta_match_results(aggregated_results)

        no_match_values.append(fast_match_failed)
        only_fast_values.append(fast_match_succeed - rosetta_match_succeed)
        rosetta_match_values.append(rosetta_match_succeed)

    make_stacked_bar_plots(no_match_values, only_fast_values, rosetta_match_values, '', x_ticks=x_ticks)


def plot_match_vs_n_scaffold(aggregated_results, total_n_scaffold=1000, output_file=None):
    
    def make_plot(list_n_first_match, ax, label, color, text_shift=0):
        list_n_first_match_sorted = sorted(list_n_first_match)

        X = list(range(1, total_n_scaffold + 1))
        Y = [0] * total_n_scaffold

        current_position = 0

        for x in X:
            while current_position < len(list_n_first_match_sorted) and list_n_first_match_sorted[current_position] <= x:
                current_position += 1

            Y[x - 1] = current_position

        ll_plot = ax.loglog(X, Y, label=label, color=color)
       
        linreg = stats.linregress(np.log(np.array(X)), np.log(np.array(Y)))

        x_mid = np.exp(np.log(total_n_scaffold) / 2)

        ax.text(x_mid, np.min(Y) + text_shift, 'r={0:.3f} slope={1:.3f}'.format(linreg.rvalue, linreg.slope), fontsize=15, color=color)
       
        print(linreg)
        n_scaffold = 23238
        print('The extraploated value when n_scaffold = {0} is {1}'.format(n_scaffold, np.exp(linreg.slope * np.log(n_scaffold) + linreg.intercept )))


    list_n_scaffolds_to_first_fast_match = []
    list_n_scaffolds_to_first_rosetta_match = []

    for match_info in aggregated_results:
        if 'n_scaffolds_to_first_fast_match' in match_info and match_info['n_scaffolds_to_first_fast_match'] > 0:
            list_n_scaffolds_to_first_fast_match.append(match_info['n_scaffolds_to_first_fast_match'])
        if 'n_scaffolds_to_first_rosetta_match' in match_info and match_info['n_scaffolds_to_first_rosetta_match'] > 0:
            list_n_scaffolds_to_first_rosetta_match.append(match_info['n_scaffolds_to_first_rosetta_match'])
   
    fig, ax = plt.subplots()

    make_plot(list_n_scaffolds_to_first_fast_match, ax, 'fast match', '#ff7f0e', text_shift=-300)
    make_plot(list_n_scaffolds_to_first_rosetta_match, ax, 'rosetta match', '#2ca02c', text_shift=0)

    plt.legend(fontsize=15, loc='upper left')

    ax.tick_params(axis='both', which='both', labelsize=12)
    plt.xlabel('N scaffolds', fontsize=20)
    plt.ylabel('N matched binding sites', fontsize=20)
    
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())

    plt.subplots_adjust(left=0.15, bottom=0.15)
    
    if output_file is None:
        plt.show()
    else:
        plt.savefig(output_file)
    
    plt.clf()

def aggregated_results_to_dict(aggregated_results):
    '''Convert a list of aggregated_results into a dictionary.'''
    aggregated_results_dict = {}

    for i in aggregated_results:
        if 'match_info_file' in i:
            key = '/'.join(i['match_info_file'].split('/')[-4:-1])
            aggregated_results_dict[key] = i

    return aggregated_results_dict

def plot_new_match_vs_n_scaffold(aggregated_results, ref_aggregated_results, total_n_scaffold=1000, output_file=None):
    
    def make_plot(list_n_first_match, ax, label, color, text_shift=0):
        list_n_first_match_sorted = sorted(list_n_first_match)

        X = list(range(1, total_n_scaffold + 1))
        Y = [0] * total_n_scaffold

        current_position = 0

        for x in X:
            while current_position < len(list_n_first_match_sorted) and list_n_first_match_sorted[current_position] <= x:
                current_position += 1

            Y[x - 1] = current_position

        ll_plot = ax.loglog(X, Y, label=label, color=color)
        #ll_plot = ax.plot(X, Y, label=label, color=color)
       
        linreg = stats.linregress(np.log(np.array(X)), np.log(np.array(Y)))

        x_mid = np.exp(np.log(total_n_scaffold) / 2)

        ax.text(x_mid, np.min(Y) + text_shift, 'r={0:.3f} slope={1:.3f}'.format(linreg.rvalue, linreg.slope), fontsize=15, color=color)
        
        print(linreg)

    aggregated_results_dict = aggregated_results_to_dict(aggregated_results)
    ref_aggregated_results_dict = aggregated_results_to_dict(ref_aggregated_results)
   
    list_n_scaffolds_to_first_fast_match = []
    list_n_scaffolds_to_first_rosetta_match = []
    
    for match_info in aggregated_results:
        if 'match_info_file' in match_info:
            key = '/'.join(match_info['match_info_file'].split('/')[-4:-1])
        else:
            key = None

        if key in ref_aggregated_results_dict and ref_aggregated_results_dict[key]['num_rosetta_match_succeed'] > 0:
            continue
        if 'n_scaffolds_to_first_rosetta_match' in match_info and match_info['n_scaffolds_to_first_rosetta_match'] > 0:
            list_n_scaffolds_to_first_rosetta_match.append(match_info['n_scaffolds_to_first_rosetta_match'])
        
        if key in ref_aggregated_results_dict and ref_aggregated_results_dict[key]['num_fast_match_succeed'] > 0:
            continue
        if 'n_scaffolds_to_first_fast_match' in match_info and match_info['n_scaffolds_to_first_fast_match'] > 0:
            list_n_scaffolds_to_first_fast_match.append(match_info['n_scaffolds_to_first_fast_match'])

    fig, ax = plt.subplots()

    make_plot(list_n_scaffolds_to_first_fast_match, ax, 'novel fast match', 'magenta', text_shift=0)
    make_plot(list_n_scaffolds_to_first_rosetta_match, ax, 'novel rosetta match', 'cyan', text_shift=0)

    plt.legend(fontsize=15, loc='upper left')

    ax.tick_params(axis='both', which='both', labelsize=12)
    plt.xlabel('N scaffolds', fontsize=20)
    plt.ylabel('N matched binding sites', fontsize=20)
    
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax.get_yaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())

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


    dict_aggregated_results = {}

    for k in aggregated_results_files.keys():
        with open(aggregated_results_files[k], 'r') as f:
            dict_aggregated_results[k] = json.load(f)

    


#    print_average_successful_running_times(dict_aggregated_results['native_rossmann'])
#    print_average_successful_running_times(dict_aggregated_results['native_ntf2'])
#    print_average_successful_running_times(dict_aggregated_results['2lv8'])
#    print_average_successful_running_times(dict_aggregated_results['5tpj'])

    print_number_of_matches(dict_aggregated_results)

#    plot_match_results_vs_lig_and_site_size(dict_aggregated_results['native_rossmann'])

#    plot_match_results_vs_site_size(dict_aggregated_results['2lv8'], x_max=30, output_file='./plots/match_results_vs_site_size_2lv8.svg')
#    plot_match_results_vs_site_size(dict_aggregated_results['5tpj'], x_max=30, output_file='./plots/match_results_vs_site_size_5tpj.svg')
#    plot_match_results_vs_site_size(dict_aggregated_results['native_rossmann'], x_max=30, output_file='./plots/match_results_vs_site_size_native_rossmann.svg')
#    plot_match_results_vs_site_size(dict_aggregated_results['native_ntf2'], x_max=30, output_file='./plots/match_results_vs_site_size_native_ntf2.svg')

#    plot_match_results_vs_ligand_size(dict_aggregated_results['3res_5tpj'])
    
#    plot_success_for_different_data_sets([dict_aggregated_results[k] for k in ['2lv8', '5tpj', 'native_rossmann', 'native_ntf2']], ['a', 'b', 'c', 'd'])
#    plot_success_for_different_data_sets([aggregated_results[k] for k in ['3res_2lv8', '3res_5tpj', '3res_native_rossmann', '3res_native_ntf2']], ['a', 'b', 'c', 'd'])

#    plot_match_vs_n_scaffold(dict_aggregated_results['3res_2lv8'], output_file='./plots/match_vs_n_scaffolds_3res_2lv8.svg')
#    plot_match_vs_n_scaffold(dict_aggregated_results['3res_5tpj'], output_file='./plots/match_vs_n_scaffolds_3res_5tpj.svg')
#    plot_match_vs_n_scaffold(dict_aggregated_results['3res_native_rossmann'], total_n_scaffold=20, output_file='./plots/match_vs_n_scaffolds_3res_native_rossmann.svg')
#    plot_match_vs_n_scaffold(dict_aggregated_results['3res_native_ntf2'], total_n_scaffold=103, output_file='./plots/match_vs_n_scaffolds_3res_native_ntf2.svg')

#    plot_match_vs_n_scaffold(dict_aggregated_results['2lv8'])
#    plot_match_vs_n_scaffold(dict_aggregated_results['5tpj'])
#    plot_match_vs_n_scaffold(dict_aggregated_results['native_rossmann'], total_n_scaffold=20)
#    plot_match_vs_n_scaffold(dict_aggregated_results['native_ntf2'], total_n_scaffold=103)

#    plot_new_match_vs_n_scaffold(dict_aggregated_results['3res_2lv8'], dict_aggregated_results['3res_native_rossmann'], 
#            output_file='./plots/new_match_vs_n_scaffolds_3res_2lv8.svg')
#    plot_new_match_vs_n_scaffold(dict_aggregated_results['3res_5tpj'], dict_aggregated_results['3res_native_ntf2'], 
#            output_file='./plots/new_match_vs_n_scaffolds_3res_5tpj.svg')
    
