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
made_df.to_csv(dir_path / "compounds" / "made.csv", index=False)

### GET ORDERED MOLS ###
ordered_df = pd.read_csv(dir_path / "../orders_data/all_ordered_mols.csv")
synthesis_df = all_df.loc[
    (all_df["SMILES"].isin(list(ordered_df.SMILES)))
    & (~all_df["SMILES"].isin(list(received_df.SMILES)))
]
synthesis_df.to_csv(dir_path / "compounds" / "for_synthesis.csv", index=False)

### GET VIRTUAL MOLS ###
virtual_df = all_df
virtual_df = all_df.loc[
    (~all_df["SMILES"].isin(list(ordered_df.SMILES)))
    & (~all_df["SMILES"].isin(list(received_df.SMILES)))
]
virtual_df.to_csv(dir_path / "compounds" / "virtual.csv", index=False)

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

all_assay_df.to_csv(
    dir_path / "experimental_results" / "protease_assay.csv", index=False
)
