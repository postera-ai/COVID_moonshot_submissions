This folder contains all of the compound shipments made to COVID Moonshot consortium by others who are synthesizing the compounds.

The `all_received_mols.csv` file contains all of the molecules received. The file is created by the `get_all_received_mols.py` script.

The format for folders containing order csvs is as follows:
`{YYYYMMDD}_{vendor}_shipment/`. 

Each folder contains files with similar naming conventions:
`{YYYYMMDD}_{vendor}_to_{receiver}.csv`

In this field, the following are required:
- **SMILES**: the standardized and canonicalized smiles: This is done using RDKit for Canonicalization following standardization by ChEMBL structure pipeline https://github.com/chembl/ChEMBL_Structure_Pipeline
- **CID**: Compound ID according to covid.postera.ai/covid
