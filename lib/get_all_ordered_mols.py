# Matt Robinson, matthew.robinson@postera.ai
# Compile all of the ordered mols from `orders_data/` into `all_ordered_mols.csv`

# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

from .utils import strip_and_standardize_smi, get_CID, get_comments


def update_orders_data(order_csv_files):

    smiles_dict = {}
    for csv_file in order_csv_files:
        try:
            order_df = pd.read_csv(csv_file)
            order_df["SMILES"] = order_df["SMILES"].apply(
                lambda x: strip_and_standardize_smi(x)
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
            raise e
            pass

    # write out final csv
    all_smiles = list(smiles_dict.keys())
    all_iks = [Chem.MolToInchiKey(Chem.MolFromSmiles(x)) for x in all_smiles]
    all_orders = [smiles_dict[x] for x in all_smiles]

    all_ordered_df = pd.DataFrame(
        {"SMILES": all_smiles, "inchikey": all_iks, "orders": all_orders}
    )

    all_ordered_df["CID"] = all_ordered_df["inchikey"].apply(
        lambda x: get_CID(x)
    )
    all_ordered_df["comments"] = all_ordered_df["inchikey"].apply(
        lambda x: get_comments(x)
    )

    all_ordered_df = all_ordered_df[["SMILES", "CID", "orders", "comments"]]
    return all_ordered_df


if __name__ == "__main__":

    # get parent path of file
    from pathlib import Path
    dir_path = Path(__file__).parent.parent.absolute()
    # get all csvs from folders
    order_csv_files = [
        f
        for f in dir_path.glob("orders/**/*.csv")
        if "all_ordered_mols.csv" not in str(f)
    ]
    orders_df = update_orders_data(order_csv_files)
    orders_df.to_csv(dir_path / "orders" / "all_ordered_mols.csv", index=False)
