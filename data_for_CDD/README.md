This folder contains the data needed for upload to the CDD Vault. Data in this folder is automatically generated using the script `get_data_for_CDD.py` using the data in other folders.

The `compounds/` folder contains three separate files
- `Compounds_Virtual.csv`: The superset of compounds. This file contains all designs submitted to covid.postera.ai/covid
- `Compounds_for_Synthesis.csv`: The compounds that haven been sent for synthesis. This file is a subset of `Compounds_Virtual.csv` and superset of `Compounds_Made.csv`
- `Compounds_Made.csv`: The compounds that have been succesfully made. This file is a subset of `Compounds_for_Sythesis.csv`

The `experimental_results/` folder cotains the experimental data for upload to CDD
- `protease_assay.csv` contains the data from all compounds that have undergone protease assays.

