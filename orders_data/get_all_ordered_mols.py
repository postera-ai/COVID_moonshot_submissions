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
all_df = pd.read_csv(dir_path / "../covid_submissions_all_info.csv")

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
            smiles_dict[smi] = str(csv_file).split("/")[-1]

    except Exception as e:
        print(f"FAILED ON {csv_file}")
        print(e)
        pass


# write out final csv
all_smiles = list(smiles_dict.keys())
all_orders = [smiles_dict[x] for x in all_smiles]

all_ordered_df = pd.DataFrame({"SMILES": all_smiles, "order": all_orders})
all_ordered_df["CID"] = all_ordered_df["SMILES"].apply(
    lambda x: list(all_df.loc[all_df["SMILES"] == x]["CID"])[0]
    if x in list(all_df.SMILES)
    else np.nan
)

achiral_all_df = all_df
achiral_all_df["SMILES"] = all_df["SMILES"].apply(
    lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x), isomericSmiles=False)
)

# this is inefficient, should not check back for list element 0
all_ordered_df["CID"] = all_ordered_df["SMILES"].apply(
    lambda x: list(achiral_all_df.loc[achiral_all_df["SMILES"] == x]["CID"])[0]
    if list(all_ordered_df.loc[all_ordered_df["SMILES"] == x]["CID"])[0] is np.nan
    else list(all_ordered_df.loc[all_ordered_df["SMILES"] == x]["CID"])[0]
)

all_ordered_df = all_ordered_df[["SMILES", "CID", "order"]]
all_ordered_df.to_csv(dir_path / "all_ordered_mols.csv", index=False)
