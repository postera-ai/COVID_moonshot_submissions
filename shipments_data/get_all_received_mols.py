# Matt Robinson, matthew.robinson@postera.ai
# Compile all of the shipped mols from `shipments_data/` into `all_received_mols.csv`

# general imports
import numpy as np
import pandas as pd

from rdkit import Chem

# get parent path of file
from pathlib import Path

dir_path = Path(__file__).parent.absolute()

import sys

sys.path.append(dir_path.parent.absolute())

from utils import strip_and_standardize_smi, get_CID, get_comments

# get all csvs from folders
received_csv_files = [
    f
    for f in dir_path.glob("**/*.csv")
    if (
        ("all_received_mols.csv" not in str(f)) and ("annotated" not in str(f))
    )
]

smiles_dict = {}
for csv_file in received_csv_files:
    try:
        received_df = pd.read_csv(csv_file)
        received_df["SMILES"] = received_df["SMILES"].apply(
            lambda x: strip_and_standardize_smi(x)
        )

        received_smi = list(received_df["SMILES"])
        for smi in received_smi:
            if smi not in smiles_dict:
                smiles_dict[smi] = str(csv_file).split("/")[-1]
            else:
                smiles_dict[smi] = (
                    smiles_dict[smi] + ", " + str(csv_file).split("/")[-1]
                )

    except Exception as e:
        print(f"FAILED ON {csv_file}")
        print(e)
        pass

# write out final csv
all_smiles = list(smiles_dict.keys())
all_iks = [Chem.MolToInchiKey(Chem.MolFromSmiles(x)) for x in all_smiles]
all_shipments = [smiles_dict[x] for x in all_smiles]

all_received_df = pd.DataFrame(
    {"SMILES": all_smiles, "inchikey": all_iks, "shipments": all_shipments}
)


all_received_df["CID"] = all_received_df["inchikey"].apply(
    lambda x: get_CID(x)
)
all_received_df["comments"] = all_received_df["inchikey"].apply(
    lambda x: get_comments(x)
)

all_received_df = all_received_df[["SMILES", "CID", "shipments", "comments"]]
all_received_df.to_csv(dir_path / "all_received_mols.csv", index=False)
