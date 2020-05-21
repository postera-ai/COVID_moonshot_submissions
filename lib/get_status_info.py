# Matt Robinson, matthew.robinson@postera.ai
# Compile all of data for csv

# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

from .utils import strip_and_standardize_smi, get_CID, get_CDD_ID, get_comments

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
made_df.to_csv(dir_path / "made.csv", index=False)

### GET ORDERED MOLS ###
ordered_df = pd.read_csv(dir_path / "../orders_data/all_ordered_mols.csv")
synthesis_df = all_df.loc[
    (all_df["SMILES"].isin(list(ordered_df.SMILES)))
    & (~all_df["SMILES"].isin(list(received_df.SMILES)))
]
synthesis_df.to_csv(dir_path / "ordered_not_made.csv", index=False)

### GET VIRTUAL MOLS ###
virtual_df = all_df
virtual_df = all_df.loc[
    (~all_df["SMILES"].isin(list(ordered_df.SMILES)))
    & (~all_df["SMILES"].isin(list(received_df.SMILES)))
]
virtual_df.to_csv(dir_path / "designed_not_ordered_nor_made.csv", index=False)

# ### GET EXPERIMENTAL RESULTS ###
# all_assay_df = pd.DataFrame()
# all_assay_csvs = (dir_path / "../experimental_data/protease_assay").glob(
#     "**/*.csv"
# )
# for assay_csv in all_assay_csvs:
#     assay_df = pd.read_csv(assay_csv)
#     assay_df["SMILES"] = assay_df["SMILES"].apply(
#         lambda x: Chem.MolToSmiles(
#             Chem.MolFromSmiles(
#                 Chem.MolToSmiles(
#                     standardizer.standardize_mol(
#                         standardizer.get_parent_mol(Chem.MolFromSmiles(x))[0]
#                     )
#                 )
#             )
#         )
#     )
#     assay_df = assay_df[
#         [
#             "SMILES",
#             "purity",
#             "volume(uL)",
#             "concentration(mM)",
#             "% Inhibition at 20 mM (N=1)",
#             "% Inhibition at 20 mM (N=2)",
#             "% Inhibition at 100 mM (N=1)",
#             "% Inhibition at 100 mM (N=2)",
#         ]
#     ]
#     assay_df["CID"] = assay_df["SMILES"].apply(
#         lambda x: list(all_df.loc[all_df["SMILES"] == x]["CID"])[0]
#         if x in list(all_df.SMILES)
#         else np.nan
#     )
#     assay_df = assay_df[
#         [
#             "SMILES",
#             "CID",
#             "purity",
#             "volume(uL)",
#             "concentration(mM)",
#             "% Inhibition at 20 mM (N=1)",
#             "% Inhibition at 20 mM (N=2)",
#             "% Inhibition at 100 mM (N=1)",
#             "% Inhibition at 100 mM (N=2)",
#         ]
#     ]
#     all_assay_df = pd.concat([all_assay_df, assay_df], axis=0)

all_assay_df = pd.read_csv(
    dir_path / "../experimental_data/protease_assay/all_inhibition_protease_assay_weizmann.csv"
)
all_assay_df["SMILES"] = all_assay_df["SMILES"].apply(
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


all_assay_smi = list(all_assay_df.SMILES)
made_not_assayed_df = made_df.loc[~made_df["SMILES"].isin(all_assay_smi)]
made_not_assayed_df.to_csv(dir_path / "made_not_assayed.csv", index=False)

made_smi = list(made_df.SMILES)
ordered_smi = list(ordered_df.SMILES)

made_cid = list(made_df.CID)
ordered_cid = list(ordered_df.CID)
assayed_cid = list(all_assay_df.CID)

made_short_ik = [Chem.MolToInchiKey(Chem.MolFromSmiles(x)).split('-')[0] for x in made_smi]
ordered_short_ik = [Chem.MolToInchiKey(Chem.MolFromSmiles(x)).split('-')[0] for x in ordered_smi]
assayed_short_ik = [Chem.MolToInchiKey(Chem.MolFromSmiles(x)).split('-')[0] for x in all_assay_smi]

all_df['short_ik'] = all_df['SMILES'].apply(lambda x: Chem.MolToInchiKey(Chem.MolFromSmiles(x)).split('-')[0])



### NOTE THAT switching to inchikey from CID means we will overcount the duplicate designs when we count orders and such

all_df["ORDERED"] = all_df["short_ik"].apply(
    lambda x: "TRUE"
    if ((x in ordered_short_ik) or (x in made_short_ik) or (x in assayed_short_ik))
    else "FALSE"
)
all_df["MADE"] = all_df["short_ik"].apply(
    lambda x: "TRUE" if ((x in made_short_ik) or (x in assayed_short_ik)) else "FALSE"
)
all_df["ASSAYED"] = all_df["short_ik"].apply(
    lambda x: "TRUE" if x in assayed_short_ik else "FALSE"
)

# all_df["ORDERED"] = all_df["CID"].apply(
#     lambda x: "TRUE"
#     if ((x in ordered_cid) or (x in made_cid) or (x in assayed_cid))
#     else "FALSE"
# )
# all_df["MADE"] = all_df["CID"].apply(
#     lambda x: "TRUE" if ((x in made_cid) or (x in assayed_cid)) else "FALSE"
# )
# all_df["ASSAYED"] = all_df["CID"].apply(
#     lambda x: "TRUE" if x in assayed_cid else "FALSE"
# )

all_df.to_csv(dir_path / "../covid_submissions_all_info.csv", index=False)
