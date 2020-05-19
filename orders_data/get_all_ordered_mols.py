# Matt Robinson, matthew.robinson@postera.ai
# Compile all of the ordered mols from `orders_data/` into `all_ordered_mols.csv`

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
id_df = pd.read_csv(dir_path / "../covid_moonshot_ids.csv")

# get all csvs from folders
order_csv_files = [
    f
    for f in dir_path.glob("**/*.csv")
    if "all_ordered_mols.csv" not in str(f)
]

smiles_dict = {}
for csv_file in order_csv_files:
    try:
        order_df = pd.read_csv(csv_file)
        order_df["SMILES"] = order_df["SMILES"].apply(
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

        order_smi = list(order_df["SMILES"])
        for smi in order_smi:
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
all_orders = [smiles_dict[x] for x in all_smiles]

all_ordered_df = pd.DataFrame(
    {"SMILES": all_smiles, "inchikey": all_iks, "orders": all_orders}
)


def get_CID(ik):
    short_ik = ik.split("-")[0]
    if ik in list(id_df["inchikey"]):
        return list(id_df.loc[id_df["inchikey"] == ik]["canonical_CID"])[0]
    elif short_ik in list(id_df["short_inchikey"]):
        return list(id_df.loc[id_df["short_inchikey"] == short_ik]["canonical_CID"])[
            0
        ]
    else:
        print("NOT FOUND")
        return np.nan


def get_comments(ik):
    short_ik = ik.split("-")[0]
    if ik in list(id_df["inchikey"]):
        return ''
    elif short_ik in list(id_df["short_inchikey"]):
        return 'imperfect match'
    else:
        return 'not found'


all_ordered_df["CID"] = all_ordered_df["inchikey"].apply(lambda x: get_CID(x))
all_ordered_df['comments'] = all_ordered_df["inchikey"].apply(lambda x: get_comments(x))

all_ordered_df = all_ordered_df[["SMILES", "CID", "orders", "comments"]]
all_ordered_df.to_csv(dir_path / "all_ordered_mols.csv", index=False)
