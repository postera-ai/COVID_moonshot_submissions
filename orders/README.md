This folder contains all of the compound orders made by the COVID Moonshot consortium and others.

The `all_ordered_mols.csv` file contains all of the molecules ordered across all orders. The file is created by the `get_all_ordered_mols.py` script.

The format for folders containing order csvs is as follows:
`{YYYYMMDD}_{vendor}_order/`. 

Each folder contains files with similar naming conventions:
`{YYYYMMDD}_{orderer}_to_{vendor}.csv`

In this field, the following are required:
- **SMILES**: the standardized and canonicalized smiles: This is done using RDKit for Canonicalization following standardization by ChEMBL structure pipeline https://github.com/chembl/ChEMBL_Structure_Pipeline
- **CID**: Compound ID according to covid.postera.ai/covid
