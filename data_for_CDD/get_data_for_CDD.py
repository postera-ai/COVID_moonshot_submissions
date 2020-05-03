# Matt Robinson, matthew.robinson@postera.ai
# Compile all of data for csv

# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

from chembl_structure_pipeline import standardizer

# get parent path of file
from pathlib import Path


# code
CID_df = pd.read_csv("https://covid.postera.ai/covid/submissions.csv")
new_CID_list = list(CID_df.CID)
old_CID_list = list(CID_df.old_CID)
old_to_new_CID_dict = {}
for old_CID, new_CID in zip(old_CID_list, new_CID_list):
    if 'None' in old_CID:
        old_to_new_CID_dict[new_CID] = new_CID
    else:
        old_to_new_CID_dict[old_CID] = new_CID
new_to_old_CID_dict = {v: k for k, v in old_to_new_CID_dict.items()}


def get_new_CID_from_old(old_CID):
    return old_to_new_CID_dict[old_CID]


def get_old_CID_from_new(new_CID):
    return new_to_old_CID_dict[new_CID]


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

### GET MADE MOLS ###
received_df = pd.read_csv(dir_path / "../shipments_data/all_received_mols.csv")
made_df = all_df.loc[all_df["SMILES"].isin(list(received_df.SMILES))]
made_df["old_CID"] = made_df.loc[:, "CID"].apply(
    lambda x: get_old_CID_from_new(x)
)
made_df.to_csv(dir_path / "compounds" / "Compounds_Made.csv", index=False)

### GET ORDERED MOLS ###
ordered_df = pd.read_csv(dir_path / "../orders_data/all_ordered_mols.csv")
# synthesis_df = all_df.loc[
#     (all_df["SMILES"].isin(list(ordered_df.SMILES)))
#     & (~all_df["SMILES"].isin(list(received_df.SMILES)))
# ]
synthesis_df = all_df.loc[all_df["SMILES"].isin(list(ordered_df.SMILES))]
synthesis_df["old_CID"] = synthesis_df.loc[:, "CID"].apply(
    lambda x: get_old_CID_from_new(x)
)
synthesis_df.to_csv(
    dir_path / "compounds" / "Compounds_for_Synthesis.csv", index=False
)

### GET VIRTUAL MOLS ###
virtual_df = all_df.copy(deep=True)
# virtual_df = all_df.loc[
#     (~all_df["SMILES"].isin(list(ordered_df.SMILES)))
#     & (~all_df["SMILES"].isin(list(received_df.SMILES)))
# ]
virtual_df["old_CID"] = virtual_df.loc[:, "CID"].apply(
    lambda x: get_old_CID_from_new(x)
)
virtual_df.to_csv(
    dir_path / "compounds" / "Compounds_Virtual.csv", index=False
)

### GET EXPERIMENTAL RESULTS ###
all_assay_df = pd.DataFrame()
all_assay_csvs = (dir_path / "../experimental_data/protease_assay").glob(
    "**/*.csv"
)
for assay_csv in all_assay_csvs:
    assay_df = pd.read_csv(assay_csv)
    assay_df["SMILES"] = assay_df["SMILES"].apply(
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
    assay_df = assay_df[
        [
            "SMILES",
            "purity",
            "volume(uL)",
            "concentration(mM)",
            "% Inhibition at 20 mM (N=1)",
            "% Inhibition at 20 mM (N=2)",
            "% Inhibition at 100 mM (N=1)",
            "% Inhibition at 100 mM (N=2)",
        ]
    ]
    assay_df["CID"] = assay_df["SMILES"].apply(
        lambda x: list(all_df.loc[all_df["SMILES"] == x]["CID"])[0]
        if x in list(all_df.SMILES)
        else np.nan
    )

    achiral_all_df = all_df
    achiral_all_df["SMILES"] = all_df["SMILES"].apply(
        lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x), isomericSmiles=False)
    )

    assay_df["CID"] = assay_df["SMILES"].apply(
        lambda x: list(
            achiral_all_df.loc[achiral_all_df["SMILES"] == x]["CID"]
        )[0]
        if list(assay_df.loc[assay_df["SMILES"] == x]["CID"])[0] is np.nan
        else list(assay_df.loc[assay_df["SMILES"] == x]["CID"])[0]
    )

    assay_df = assay_df[
        [
            "SMILES",
            "CID",
            "purity",
            "volume(uL)",
            "concentration(mM)",
            "% Inhibition at 20 mM (N=1)",
            "% Inhibition at 20 mM (N=2)",
            "% Inhibition at 100 mM (N=1)",
            "% Inhibition at 100 mM (N=2)",
        ]
    ]
    all_assay_df = pd.concat([all_assay_df, assay_df], axis=0)

all_assay_df["old_CID"] = all_assay_df.loc[:, "CID"].apply(
    lambda x: get_old_CID_from_new(x)
)
all_assay_df.to_csv(
    dir_path / "experimental_results" / "protease_assay.csv", index=False
)

### GET SOLUBILITY DATA ###

# I have no idea why I have to run this again, but it makes it work
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

all_sol_df = pd.DataFrame()
all_sol_csvs = (dir_path / "../experimental_data/solubility").glob("**/*.csv")
for sol_csv in all_sol_csvs:
    sol_df = pd.read_csv(sol_csv)
    sol_df["SMILES"] = sol_df["SMILES"].apply(
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
    sol_df = sol_df[
        [
            "SMILES",
            "Raw signal @20 µM",
            "Relative solubility @20 µM",
            "Raw signal @100 µM",
            "Relative solubility @100 µM",
            "100 µM / 20 µM",
        ]
    ]
    sol_df["CID"] = sol_df["SMILES"].apply(
        lambda x: list(all_df.loc[all_df["SMILES"] == x]["CID"])[0]
        if x in list(all_df.SMILES)
        else np.nan
    )

    achiral_all_df = all_df
    achiral_all_df["SMILES"] = all_df["SMILES"].apply(
        lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x), isomericSmiles=False)
    )

    sol_df["CID"] = sol_df["SMILES"].apply(
        lambda x: list(
            achiral_all_df.loc[achiral_all_df["SMILES"] == x]["CID"]
        )[0]
        if list(sol_df.loc[sol_df["SMILES"] == x]["CID"])[0] is np.nan
        else list(sol_df.loc[sol_df["SMILES"] == x]["CID"])[0]
    )

    sol_df = sol_df[
        [
            "SMILES",
            "CID",
            "Raw signal @20 µM",
            "Relative solubility @20 µM",
            "Raw signal @100 µM",
            "Relative solubility @100 µM",
            "100 µM / 20 µM",
        ]
    ]
    all_sol_df = pd.concat([all_sol_df, sol_df], axis=0)


all_sol_df["old_CID"] = all_sol_df.loc[:, "CID"].apply(
    lambda x: get_old_CID_from_new(x)
)
all_sol_df.to_csv(
    dir_path / "experimental_results" / "solubility.csv", index=False
)
