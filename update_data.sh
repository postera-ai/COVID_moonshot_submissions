#!/bin/bash

python ./orders_data/get_all_ordered_mols.py
python ./shipments_data/get_all_received_mols.py
python ./shipments_data/create_diamond_files.py
python ./data_for_CDD/get_data_for_CDD.py
python ./compound_tracking/get_status_info.py
python ./compound_tracking/create_tracking_plot.py

