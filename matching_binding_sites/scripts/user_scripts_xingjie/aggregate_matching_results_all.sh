#!/bin/bash

./scripts/aggregate_matching_results.py matching_data/matching_site_to_native_rossmann_no_dump matching_data/aggregated_matching_to_native_rossmann
./scripts/aggregate_matching_results.py matching_data/matching_site_to_native_ntf2_no_dump matching_data/aggregated_matching_to_native_ntf2
./scripts/aggregate_matching_results.py matching_data/matching_to_2lv8_two_LHL_no_dump matching_data/aggregated_matching_to_2lv8_two_LHL
./scripts/aggregate_matching_results.py matching_data/matching_site_to_5tpj_no_dump matching_data/aggregated_matching_to_5tpj

./scripts/aggregate_matching_results.py matching_data/matching_3res_site_to_native_rossmann_no_dump matching_data/aggregated_matching_3res_site_to_native_rossmann 3
./scripts/aggregate_matching_results.py matching_data/matching_3res_site_to_native_ntf2_no_dump matching_data/aggregated_matching_3res_site_to_native_ntf2 3
./scripts/aggregate_matching_results.py matching_data/matching_3res_site_to_2lv8_two_LHL_no_dump matching_data/aggregated_matching_3res_site_to_2lv8_two_LHL 3
./scripts/aggregate_matching_results.py matching_data/matching_3res_site_to_5tpj_no_dump matching_data/aggregated_matching_3res_site_to_5tpj 3
