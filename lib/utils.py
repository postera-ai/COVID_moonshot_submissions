### util functions for parsing all the moonshot data
### matthew.robinson@postera.ai

# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

from chembl_structure_pipeline import standardizer

# get parent path of file
from pathlib import Path

dir_path = Path(__file__).parent.absolute()
# all_df = pd.read_csv(dir_path / "../covid_submissions_all_info.csv")
id_df = pd.read_csv(dir_path / "../covid_moonshot_ids.csv")
cdd_df = pd.read_csv(
    dir_path / "../data_for_CDD/current_vault_data/current_vault_data.csv"
)
CID_df = pd.read_csv("https://covid.postera.ai/covid/submissions.csv")


def get_CID(ik):
    short_ik = ik.split("-")[0]
    if ik in list(id_df["inchikey"]):
        return list(id_df.loc[id_df["inchikey"] == ik]["canonical_CID"])[0]
    elif short_ik in list(id_df["short_inchikey"]):
        return list(
            id_df.loc[id_df["short_inchikey"] == short_ik]["canonical_CID"]
        )[0] # this will pick up the first one, which is what we want when enantiopures are separated
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
        return f"imperfect stereochemical match for {list(id_df.loc[id_df['short_inchikey']==short_ik]['canonical_CID'])[0]}"
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
new_CID_list = list(CID_df["CID"])
old_CID_list = [str(x) for x in list(CID_df["CID (old format)"])]
old_to_new_CID_dict = {}
for old_CID, new_CID in zip(old_CID_list, new_CID_list):
    if "nan" in old_CID:
        old_to_new_CID_dict[new_CID] = new_CID
    else:
        old_to_new_CID_dict[old_CID] = new_CID
new_to_old_CID_dict = {v: k for k, v in old_to_new_CID_dict.items()}


def get_new_CID_from_old(old_CID):
    return old_to_new_CID_dict[old_CID]


def get_old_CID_from_new(new_CID):
    return new_to_old_CID_dict[new_CID]


def get_series(smi):
    series_SMARTS_dict = {
        # "3-aminopyridine": "[R1][C,N;R0;!$(NC(=O)CN)]C(=O)[C,N;R0;!$(NC(=O)CN)][c]1cnccc1",
        "Ugi": "[c,C:1][C](=[O])[N]([c,C,#1:2])[C]([c,C,#1:3])([c,C,#1:4])[C](=[O])[NH1][c,C:5]",
        "Isatins": "O=C1Nc2ccccc2C1=O",
        "3-aminopyridine-like": "[cR1,cR2]-[C,N]C(=O)[C,N]!@[R1]",
        "quinolones": "NC(=O)c1cc(=O)[nH]c2ccccc12",
        "piperazine-chloroacetamide": "O=C(CCl)N1CCNCC1",
        "activated-ester": "O=C(Oc1cncc(Cl)c1)c1cccc2[nH]ccc12"
    }

    def check_if_smi_in_series(
        smi, SMARTS, MW_cutoff=550, num_atoms_cutoff=70, num_rings_cutoff=10
    ):
        mol = Chem.MolFromSmiles(smi)
        MW = Chem.Descriptors.MolWt(mol)
        num_heavy_atoms = mol.GetNumHeavyAtoms()
        num_rings = Chem.rdMolDescriptors.CalcNumRings(mol)
        patt = Chem.MolFromSmarts(SMARTS)
        if (
            (
                len(
                    Chem.AddHs(Chem.MolFromSmiles(smi)).GetSubstructMatches(
                        patt
                    )
                )
                > 0
            )
            and (MW <= MW_cutoff)
            and (num_heavy_atoms <= num_atoms_cutoff)
            and (num_rings <= num_rings_cutoff)
        ):
            return True
        else:
            return False

    for series in series_SMARTS_dict:
        series_SMARTS = series_SMARTS_dict[series]
        if series == "3-amonipyridine-like":
            if check_if_smi_in_series(
                smi,
                series_SMARTS,
                MW_cutoff=450,
                num_rings_cutoff=4,
                num_atoms_cutoff=35,
            ):
                return series
        else:
            if check_if_smi_in_series(smi, series_SMARTS):
                return series
    return None
