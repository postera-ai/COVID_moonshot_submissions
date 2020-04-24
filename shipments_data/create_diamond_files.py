# Matt Robinson, matthew.robinson@postera.ai
# Compile all of the shipped mols from `shipments_data/` into `all_received_mols.csv`

# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

from chembl_structure_pipeline import standardizer

# get parent path of file
from pathlib import Path

dir_path = Path(__file__).parent.absolute()
all_df = pd.read_csv(dir_path / "../covid_submissions_all_info.csv")

# get all csvs from folders
received_csv_files = [
    f
    for f in dir_path.glob("**/*.csv")
    if "all_received_mols.csv" not in str(f)
]

smiles_dict = {}
for csv_file in received_csv_files:
    if "xchem" in str(csv_file):
        try:
            received_df = pd.read_csv(csv_file)
            # received_df["SMILES"] = received_df["SMILES"].apply(
            #     lambda x: Chem.MolToSmiles(
            #         Chem.MolFromSmiles(
            #             Chem.MolToSmiles(
            #                 standardizer.standardize_mol(
            #                     standardizer.get_parent_mol(Chem.MolFromSmiles(x))[
            #                         0
            #                     ]
            #                 )
            #             )
            #         )
            #     )
            # )
            new_filename = str(csv_file).split("/")[-1].replace("xchem", "diamond")

            received_df["Library Name"] = str(csv_file).split("/")[-1]
            received_df = received_df[
                [
                    "plate_ID",
                    "well",
                    "Library Name",
                    "SMILES",
                    "shipment_ID",
                    "volume(uL)",
                    "concentration(mM)",
                ]
            ]
            received_df = received_df.rename(
                columns={
                    "plate_ID": "Library Plate",
                    "well": "Source well",
                    "shipment_ID": "Code",
                    "volume(uL)": "Volume (ul)",
                    "concentration(mM)": "Concentration (mM)",
                }
            )

            received_df.to_csv(
                dir_path / "diamond_files" / new_filename, index=False
            )

        except Exception as e:
            print(f"FAILED ON {csv_file}")
            print(e)
            pass
