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
id_df = pd.read_csv(dir_path / "../covid_moonshot_ids.csv")
cdd_df = pd.read_csv(dir_path / "../data_for_CDD/current_vault_data/current_vault_data.csv")

def get_CID(ik):
    short_ik = ik.split("-")[0]
    if ik in list(id_df["inchikey"]):
        return list(id_df.loc[id_df["inchikey"] == ik]["canonical_CID"])[0]
    elif short_ik in list(id_df["short_inchikey"]):
        return list(
            id_df.loc[id_df["short_inchikey"] == short_ik]["canonical_CID"]
        )[0]
    else:
        print("NOT FOUND")
        return np.nan


def get_CDD_ID(external_id):
    if external_id in list(cdd_df['external_ID']):
        return list(cdd_df.loc[cdd_df['external_ID']==external_id]['CDD_name'])[0]
    else:
        print("NOT FOUND")
        return np.nan


def get_comments(ik):
    short_ik = ik.split("-")[0]
    if ik in list(id_df["inchikey"]):
        return ""
    elif short_ik in list(id_df["short_inchikey"]):
        return "imperfect match"
    else:
        return "not found"

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
            received_df["standard_SMILES"] = received_df["SMILES"].apply(
                lambda x: Chem.MolToSmiles(
                    Chem.MolFromSmiles(
                        Chem.MolToSmiles(
                            standardizer.standardize_mol(
                                standardizer.get_parent_mol(
                                    Chem.MolFromSmiles(x)
                                )[0]
                            )
                        )
                    )
                )
            )
            received_df["inchikey"] = received_df["standard_SMILES"].apply(
                lambda x: Chem.MolToInchiKey(Chem.MolFromSmiles(x))
            )
            received_df["external_ID"] = received_df["inchikey"].apply(
                lambda x: get_CID(x)
            )
            received_df["CDD_ID"] = received_df["external_ID"].apply(
                lambda x: get_CDD_ID(x)
            )
            received_df["PostEra_comments"] = received_df["inchikey"].apply(
                lambda x: get_comments(x)
            )
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
                    "external_ID",
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
