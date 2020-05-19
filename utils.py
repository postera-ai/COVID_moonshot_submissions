### util functions for parsing all the moonshot data
### matthew.robinson@postera.ai

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
# all_df = pd.read_csv(dir_path / "../covid_submissions_all_info.csv")
id_df = pd.read_csv(dir_path / "covid_moonshot_ids.csv")
cdd_df = pd.read_csv(
    dir_path / "data_for_CDD/current_vault_data/current_vault_data.csv"
)
CID_df = pd.read_csv("https://covid.postera.ai/covid/submissions.csv")


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
    if external_id in list(cdd_df["external_ID"]):
        return list(
            cdd_df.loc[cdd_df["external_ID"] == external_id]["CDD_name"]
        )[0]
    else:
        print("NOT FOUND")
        return np.nan


def get_comments(ik):
    short_ik = ik.split("-")[0]
    if ik in list(id_df["inchikey"]):
        return ""
    elif short_ik in list(id_df["short_inchikey"]):
        return f"imperfect match for {list(id_df.loc[id_df['short_inchikey']==short_ik]['canonical_CID'])[0]}"
    else:
        return "not found"


def strip_and_standardize_smi(smi):
    return Chem.MolToSmiles(
        Chem.MolFromSmiles(
            Chem.MolToSmiles(
                standardizer.standardize_mol(
                    standardizer.get_parent_mol(Chem.MolFromSmiles(smi))[0]
                )
            )
        )
    )


# code to retrieve new and old CIDS
new_CID_list = list(CID_df.CID)
old_CID_list = list(CID_df.old_CID)
old_to_new_CID_dict = {}
for old_CID, new_CID in zip(old_CID_list, new_CID_list):
    if "None" in old_CID:
        old_to_new_CID_dict[new_CID] = new_CID
    else:
        old_to_new_CID_dict[old_CID] = new_CID
new_to_old_CID_dict = {v: k for k, v in old_to_new_CID_dict.items()}


def get_new_CID_from_old(old_CID):
    return old_to_new_CID_dict[old_CID]


def get_old_CID_from_new(new_CID):
    return new_to_old_CID_dict[new_CID]
