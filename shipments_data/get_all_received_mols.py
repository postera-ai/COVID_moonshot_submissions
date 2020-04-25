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
    try:
        received_df = pd.read_csv(csv_file)
        received_df["SMILES"] = received_df["SMILES"].apply(
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

        received_smi = list(received_df["SMILES"])
        for smi in received_smi:
            smiles_dict[smi] = str(csv_file).split("/")[-1]

    except Exception as e:
        print(f"FAILED ON {csv_file}")
        print(e)
        pass


# write out final csv
all_smiles = list(smiles_dict.keys())
all_shipments = [smiles_dict[x] for x in all_smiles]

all_received_df = pd.DataFrame(
    {"SMILES": all_smiles, "shipment": all_shipments}
)
all_received_df["CID"] = all_received_df["SMILES"].apply(
    lambda x: list(all_df.loc[all_df["SMILES"] == x]["CID"])[0]
    if x in list(all_df.SMILES)
    else np.nan
)

achiral_all_df = all_df
achiral_all_df["SMILES"] = all_df["SMILES"].apply(
    lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x), isomericSmiles=False)
)

all_received_df["CID"] = all_received_df["SMILES"].apply(
    lambda x: list(achiral_all_df.loc[achiral_all_df["SMILES"] == x]["CID"])[0]
    if x == np.nan
    else list(all_received_df.loc[all_received_df["SMILES"] == x]["CID"])[0]
)

all_received_df = all_received_df[["SMILES", "CID", "shipment"]]
all_received_df.to_csv(dir_path / "all_received_mols.csv", index=False)
