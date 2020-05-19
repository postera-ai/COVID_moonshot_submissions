# Matt Robinson, matthew.robinson@postera.ai
# master script to update data for COVID Moonshot project

# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

from utils import strip_and_standardize_smi, get_CID, get_CDD_ID, get_comments

# get parent path of file
from pathlib import Path

dir_path = Path(__file__).parent.absolute()

### First Get the data ###
all_df = pd.read_csv("https://covid.postera.ai/covid/submissions.csv")
all_df["SMILES"] = all_df["SMILES"].apply(
    lambda x: strip_and_standardize_smi(x)
)
all_df.to_csv(dir_path / "covid_submissions_all_info.csv", index=False)
id_df = pd.read_csv(dir_path / "covid_moonshot_ids.csv")

### Update the Moonshot IDs
from update_moonshot_ids import update_CIDs

id_df = update_CIDs(all_df, id_df)
id_df.to_csv(dir_path / "covid_moonshot_ids.csv", index=False)

### Update the orders data
from get_all_ordered_mols import update_orders_data

order_csv_files = [
    f
    for f in dir_path.glob("orders/**/*.csv")
    if "all_ordered_mols.csv" not in str(f)
]
orders_df = update_orders_data(order_csv_files)
orders_df.to_csv(dir_path / "orders" / "all_ordered_mols.csv", index=False)

### Update shipments data
from get_all_received_mols import (
    update_shipments_data,
    create_diamond_files,
    create_weizmann_files,
)

received_csv_files = [
    f
    for f in dir_path.glob("shipments/**/*.csv")
    if (
        ("all_received_mols.csv" not in str(f))
        and ("annotated" not in str(f))
        and ("diamond" not in str(f))
    )
]
shipments_df = update_shipments_data(received_csv_files)
shipments_df.to_csv(
    dir_path / "shipments" / "all_received_mols.csv", index=False
)

diamond_dfs = create_diamond_files(received_csv_files)
for diamond_df, diamond_fn in diamond_dfs:
    diamond_df.to_csv(
        dir_path / "shipments" / "diamond_files" / diamond_fn, index=False
    )

weizmann_dfs = create_weizmann_files(received_csv_files)
for weizmann_df, weizmann_fn in weizmann_dfs:
    weizmann_df.to_csv(
        dir_path / "shipments" / "weizmann_files" / weizmann_fn, index=False
    )
