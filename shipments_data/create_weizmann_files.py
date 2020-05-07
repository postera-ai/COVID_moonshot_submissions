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
all_df["SMILES"] = all_df["SMILES"].apply(
    lambda x: Chem.MolToSmiles(
        Chem.MolFromSmiles(
            Chem.MolToSmiles(
                standardizer.standardize_mol(
                    standardizer.get_parent_mol(Chem.MolFromSmiles(x))[0]
                )
            )
        )
    )
)

achiral_all_df = all_df.copy()
achiral_all_df["SMILES"] = all_df["SMILES"].apply(
    lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x), isomericSmiles=False)
)

# code
CID_df = pd.read_csv("https://covid.postera.ai/covid/submissions.csv")
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

# get all csvs from folders
received_csv_files = [
    f
    for f in dir_path.glob("**/*.csv")
    if (("all_received_mols.csv" not in str(f)) and ("weizmann_annotated" not in str(f)))
]

smiles_dict = {}
for csv_file in received_csv_files:
    if "weizmann" in str(csv_file):
        try:
            received_df = pd.read_csv(csv_file)
            received_df["standard_SMILES"] = received_df["SMILES"].apply(
                lambda x: Chem.MolToSmiles(
                    Chem.MolFromSmiles(
                        Chem.MolToSmiles(
                            standardizer.standardize_mol(
                                standardizer.get_parent_mol(Chem.MolFromSmiles(x))[
                                    0
                                ]
                            )
                        )
                    )
                )
            )
            new_filename = str(csv_file).split("/")[-1].replace("weizmann", "weizmann_annotated")

            CID_dict = {}
            for smi in list(all_df.SMILES):
                inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
                CID_dict[inchikey] = list(all_df.loc[all_df["SMILES"] == smi]["CID"])[
                    0
                ]
            for smi in list(achiral_all_df.SMILES):
                inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
                CID_dict[inchikey] = list(
                    achiral_all_df.loc[achiral_all_df["SMILES"] == smi]["CID"]
                )[0]

            def get_CID(smi):
                inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
                if inchikey in CID_dict:
                    return CID_dict[inchikey]
                no_stereo_inchikey = Chem.MolToInchiKey(
                    Chem.MolFromSmiles(
                        Chem.MolToSmiles(Chem.MolFromSmiles(smi), isomericSmiles=False)
                    )
                )
                if no_stereo_inchikey in CID_dict:
                    return CID_dict[no_stereo_inchikey]
                else:
                    print(smi)
                    return None

            received_df["CID"] = received_df["standard_SMILES"].apply(lambda x: get_CID(x))
            received_df["old_CID"] = received_df.loc[:, "CID"].apply(
                lambda x: get_old_CID_from_new(x)
            )
            received_df = received_df.rename(columns={"old_CID": "external_ID"})

            received_df = received_df[[
                "SMILES",
                "external_ID",
                "shipment_ID",
                "catalog_ID",
                "sample_MW",
                "purity",
                "volume(uL)",
                "concentration(mM)",
                "plate_ID",
                "well",
                "stereochemistry",
                "CID"
            ]]

            received_df.to_csv(
                dir_path / "weizmann_files" / new_filename, index=False
            )

        except Exception as e:
            print(f"FAILED ON {csv_file}")
            print(e)
            pass