# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

from chembl_structure_pipeline import standardizer

import requests
import sys
import os

# get parent path of file
from pathlib import Path

dir_path = Path(__file__).parent.absolute()

### GET CURRENT CDD COMPOUNDS
current_cdd_df = pd.read_csv(
    dir_path / "current_vault_data" / "current_vault_data.csv"
)
cdd_virtual_external_id_list = list(
    current_cdd_df.loc[current_cdd_df["virtual_project"] == True][
        "external_ID"
    ]
)
cdd_synthesis_external_id_list = list(
    current_cdd_df.loc[current_cdd_df["for_synthesis_project"] == True][
        "external_ID"
    ]
)
cdd_made_external_id_list = list(
    current_cdd_df.loc[current_cdd_df["made_project"] == True]["external_ID"]
)

### GET CURRENT MOONSHOT SUBMISSIONS
virtual_df = pd.read_csv(dir_path / "compounds" / "Compounds_Virtual.csv")
synthesis_df = pd.read_csv(
    dir_path / "compounds" / "Compounds_for_Synthesis.csv"
)
made_df = pd.read_csv(dir_path / "compounds" / "Compounds_Made.csv")

### FILTER
add_to_virtual_df = virtual_df.loc[
    ~virtual_df["old_CID"].isin(cdd_virtual_external_id_list)
]
add_to_virtual_df = add_to_virtual_df.reset_index(drop=True)
add_to_virtual_df = add_to_virtual_df[
    ["SMILES", "old_CID", "creator", "fragments", "covalent_warhead"]
]
add_to_virtual_df.to_csv(
    dir_path / "vault_updates" / "add_to_virtual_df.csv", index=False
)

add_to_synthesis_df = synthesis_df.loc[
    ~synthesis_df["old_CID"].isin(cdd_synthesis_external_id_list)
]
add_to_synthesis_df = add_to_synthesis_df.reset_index(drop=True)
add_to_synthesis_df = add_to_synthesis_df[
    ["SMILES", "old_CID", "creator", "fragments", "covalent_warhead"]
]
add_to_synthesis_df.to_csv(
    dir_path / "vault_updates" / "add_to_synthesis_df.csv", index=False
)

add_to_made_df = made_df.loc[
    ~made_df["old_CID"].isin(cdd_made_external_id_list)
]
add_to_made_df = add_to_made_df.reset_index(drop=True)
add_to_made_df = add_to_made_df[
    ["SMILES", "old_CID", "creator", "fragments", "covalent_warhead"]
]
add_to_made_df.to_csv(
    dir_path / "vault_updates" / "add_to_made_df.csv", index=False
)
