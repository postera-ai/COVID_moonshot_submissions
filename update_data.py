# Matt Robinson, matthew.robinson@postera.ai
# master script to update data for COVID Moonshot project

# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

from lib.utils import (
    strip_and_standardize_smi,
    get_CID,
    get_CDD_ID,
    get_comments,
)

# get parent path of file
from pathlib import Path

dir_path = Path(__file__).parent.absolute()

### First Get the data ###
all_df = pd.read_csv("https://covid.postera.ai/covid/submissions.csv")
all_df["SMILES"] = all_df["SMILES"].apply(
    lambda x: strip_and_standardize_smi(x)
)
all_df = all_df.rename(
    columns={
        "Submission Creator": "creator",
        "Submission Rationale": "rationale",
        "Submission Fragments": "fragments",
        "Covalent Warhead": "covalent_warhead",
    }
)
all_df.to_csv(dir_path / "covid_submissions_all_info.csv", index=False)


def create_old_cid(x):
    if x["old_CID"] is np.nan:
        return x["CID"]
    else:
        return x["old_CID"]


id_df = all_df.copy()[["SMILES", "CID", "CID (canonical)", "CID (old format)"]]
id_df = id_df.rename(
    columns={"CID (old format)": "old_CID", "CID (canonical)": "canonical_CID"}
)
id_df["inchikey"] = id_df["SMILES"].apply(
    lambda x: Chem.MolToInchiKey(Chem.MolFromSmiles(x))
)
id_df["short_inchikey"] = id_df["inchikey"].apply(lambda x: x.split("-")[0])
id_df["old_CID"] = id_df.apply(create_old_cid, axis=1)
id_df.to_csv(dir_path / "covid_moonshot_ids.csv", index=False)

# ### Update the Moonshot IDs (replaced by above)
# from lib.update_moonshot_ids import update_CIDs

# id_df = update_CIDs(all_df, id_df)
# id_df.to_csv(dir_path / "covid_moonshot_ids.csv", index=False)

### Update the orders data
from lib.get_all_ordered_mols import update_orders_data

order_csv_files = [
    f
    for f in dir_path.glob("orders/**/*.csv")
    if "all_ordered_mols.csv" not in str(f)
]
orders_df = update_orders_data(order_csv_files)
orders_df.to_csv(dir_path / "orders" / "all_ordered_mols.csv", index=False)

### Update shipments data
from lib.get_all_received_mols import (
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

### Update CDD info

# first get the necessary data
received_df = pd.read_csv(dir_path / "shipments" / "all_received_mols.csv")
made_df = received_df.copy()
made_df.to_csv(
    dir_path / "data_for_CDD" / "compounds" / "Compounds_Made.csv", index=False
)

ordered_df = pd.read_csv(dir_path / "orders" / "all_ordered_mols.csv")
synthesis_df = ordered_df.copy()
synthesis_df.to_csv(
    dir_path / "data_for_CDD" / "compounds" / "Compounds_for_Synthesis.csv",
    index=False,
)

virtual_df = all_df.copy()
virtual_df.to_csv(
    dir_path / "data_for_CDD" / "compounds" / "Compounds_Virtual.csv",
    index=False,
)

# get current data in the vault
from lib.get_current_vault_data import get_current_vault_data

current_cdd_df = get_current_vault_data()
current_cdd_df.to_csv(
    dir_path
    / "data_for_CDD"
    / "current_vault_data"
    / "current_vault_data.csv",
    index=False,
)

# update master file
all_df["CDD_name"] = all_df["CID"].apply(
    lambda x: current_cdd_df.loc[current_cdd_df["external_ID"] == x][
        "CDD_name"
    ].item()
    if (x in list(current_cdd_df["external_ID"]))
    else np.nan
)
all_df["CDD_mol_ID"] = all_df["CID"].apply(
    lambda x: current_cdd_df.loc[current_cdd_df["external_ID"] == x][
        "molecule_ID"
    ].item()
    if (x in list(current_cdd_df["external_ID"]))
    else np.nan
)

# get the necessary updates to CDD
from lib.get_CDD_updates import get_CDD_updates

add_to_virtual_df, add_to_synthesis_df, add_to_made_df = get_CDD_updates(
    all_df, current_cdd_df, virtual_df, synthesis_df, made_df
)

add_to_virtual_df.to_csv(
    dir_path / "data_for_CDD" / "vault_updates" / "add_to_virtual_df.csv",
    index=False,
)
add_to_synthesis_df.to_csv(
    dir_path / "data_for_CDD" / "vault_updates" / "add_to_synthesis_df.csv",
    index=False,
)
add_to_made_df.to_csv(
    dir_path / "data_for_CDD" / "vault_updates" / "add_to_made_df.csv",
    index=False,
)

# get assay data
from lib.get_experimental_data import (
    get_rapidfire_inhibition_data,
    get_rapidfire_IC50_data,
    get_fluorescense_inhibition_data,
    get_fluorescense_IC50_data,
    get_solubility_data,
    get_trypsin_data,
)

rapidfire_inhibition_df = get_rapidfire_inhibition_data()
rapidfire_IC50_df = get_rapidfire_IC50_data()

fluorescence_inhibition_df = get_fluorescense_inhibition_data()
fluorescence_IC50_df = get_fluorescense_IC50_data()

solubility_df = get_solubility_data()
trypsin_df = get_trypsin_data()

all_df = pd.merge(
    all_df, rapidfire_inhibition_df, how="left", on=["CDD_mol_ID"]
)
all_df = pd.merge(all_df, rapidfire_IC50_df, how="left", on=["CDD_mol_ID"])
all_df = pd.merge(
    all_df, fluorescence_inhibition_df, how="left", on=["CDD_mol_ID"]
)
all_df = pd.merge(all_df, fluorescence_IC50_df, how="left", on=["CDD_mol_ID"])
all_df = pd.merge(all_df, solubility_df, how="left", on=["CDD_mol_ID"])
all_df = pd.merge(all_df, trypsin_df, how="left", on=["CDD_mol_ID"])

### only list things with at least one inhibition value as assayed
ordered_iks = [
    Chem.MolToInchiKey(Chem.MolFromSmiles(x))
    for x in list(synthesis_df["SMILES"])
]
made_iks = [
    Chem.MolToInchiKey(Chem.MolFromSmiles(x)) for x in list(made_df["SMILES"])
]

assayed_df = all_df.loc[
    (
        (all_df["r_inhibition_at_20_uM"].notnull())
        | (all_df["r_inhibition_at_50_uM"].notnull())
        | (all_df["f_inhibition_at_20_uM"].notnull())
        | (all_df["f_inhibition_at_50_uM"].notnull())
    )
]
assayed_iks = [
    Chem.MolToInchiKey(Chem.MolFromSmiles(x))
    for x in list(assayed_df["SMILES"])
]

all_df["ORDERED"] = all_df["InChIKey"].apply(
    lambda x: "TRUE"
    if ((x in ordered_iks) or (x in made_iks) or (x in assayed_iks))
    else "FALSE"
)
all_df["MADE"] = all_df["InChIKey"].apply(
    lambda x: "TRUE" if ((x in made_iks) or (x in assayed_iks)) else "FALSE"
)
all_df["ASSAYED"] = all_df["InChIKey"].apply(
    lambda x: "TRUE" if x in assayed_iks else "FALSE"
)
all_df.to_csv(dir_path / "covid_submissions_all_info.csv", index=False)
