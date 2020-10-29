# Matt Robinson, matthew.robinson@postera.ai
# Compile all of the shipped mols from `shipments_data/` into `all_received_mols.csv`

# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

from .utils import strip_and_standardize_smi, get_CID, get_CDD_ID, get_comments


def update_shipments_data(received_csv_files):

    smiles_dict = {}
    shipment_info_dict = {}
    for csv_file in received_csv_files:
        try:
            received_df = pd.read_csv(csv_file)
            received_df["SMILES"] = received_df["SMILES"].apply(
                lambda x: strip_and_standardize_smi(x)
            )

            received_smi = list(received_df["SMILES"])
            for smi in received_smi:
                catalog_ID = list(
                    received_df.loc[received_df["SMILES"] == smi]["catalog_ID"]
                )[0]
                receiver = str(csv_file).split("/")[-1].split("to_")[-1][:-4]
                shipment_ID = list(
                    received_df.loc[received_df["SMILES"] == smi]["shipment_ID"]
                )[0]
                stereochemistry = list(
                    received_df.loc[received_df["SMILES"] == smi]["stereochemistry"]
                )[0]
                plate = list(received_df.loc[received_df["SMILES"] == smi]["plate_ID"])[
                    0
                ]
                well = list(received_df.loc[received_df["SMILES"] == smi]["well"])[0]
                if "Batch ID" in list(received_df.columns):
                    batch = list(received_df.loc[received_df["SMILES"] == smi]["Batch ID"])[0]
                else:
                    batch = "batch-unknown"
                shipment_str = f"{str(receiver)}_{str(catalog_ID)}_{str(batch)}_{str(shipment_ID)}_{str(stereochemistry)}_{str(plate)}_{str(well)}"
                if smi not in smiles_dict:
                    smiles_dict[smi] = str(csv_file).split("/")[-1]
                else:
                    smiles_dict[smi] = (
                        smiles_dict[smi] + ", " + str(csv_file).split("/")[-1]
                    )
                if smi not in shipment_info_dict:
                    shipment_info_dict[smi] = shipment_str
                else:
                    shipment_info_dict[smi] = (
                        shipment_info_dict[smi] + ", " + shipment_str
                    )

        except Exception as e:
            print(f"FAILED ON {csv_file}")
            raise e
            pass

    # write out final csv
    all_smiles = list(smiles_dict.keys())
    all_iks = [Chem.MolToInchiKey(Chem.MolFromSmiles(x)) for x in all_smiles]
    all_shipments = [smiles_dict[x] for x in all_smiles]
    all_shipment_strs = [shipment_info_dict[x] for x in all_smiles]

    all_received_df = pd.DataFrame(
        {"SMILES": all_smiles, "inchikey": all_iks,
         "shipments": all_shipments, "shipment_info": all_shipment_strs}
    )
    all_received_df["CID"] = all_received_df["inchikey"].apply(lambda x: get_CID(x))
    all_received_df["comments"] = all_received_df["inchikey"].apply(
        lambda x: get_comments(x)
    )

    all_received_df = all_received_df[["SMILES", "CID", "shipments","shipment_info", "comments"]]
    return all_received_df


def create_diamond_files(received_csv_files):

    diamond_dfs = []

    for csv_file in received_csv_files:
        if "xchem" in str(csv_file):
            try:
                received_df = pd.read_csv(csv_file)
                received_df["standard_SMILES"] = received_df["SMILES"].apply(
                    strip_and_standardize_smi
                )
                received_df["inchikey"] = received_df["standard_SMILES"].apply(
                    lambda x: Chem.MolToInchiKey(Chem.MolFromSmiles(x))
                )
                received_df["external_ID"] = received_df["inchikey"].apply(
                    lambda x: get_CID(x)
                )
                received_df["CDD_ID"] = received_df["external_ID"].apply(
                    lambda x: get_CDD_ID(x)
                )
                received_df["PostEra_comments"] = received_df["inchikey"].apply(
                    lambda x: get_comments(x)
                )
                new_filename = str(csv_file).split("/")[-1].replace("xchem", "diamond")

                if "Batch ID" not in list(received_df.columns):
                    received_df["Batch ID"] = 'UNK'

                received_df["Library Name"] = str(csv_file).split("/")[-1]
                received_df = received_df[
                    [
                        "plate_ID",
                        "well",
                        "Library Name",
                        "SMILES",
                        "shipment_ID",
                        "volume(uL)",
                        "concentration(mM)",
                        "external_ID",
                        "stereochemistry",
                        "catalog_ID",
                        "Batch ID"
                    ]
                ]
                received_df = received_df.rename(
                    columns={
                        "plate_ID": "Library Plate",
                        "well": "Source well",
                        "shipment_ID": "Code",
                        "volume(uL)": "Volume (ul)",
                        "concentration(mM)": "Concentration (mM)",
                    }
                )
                # get rid of chemaxon style smiles
                received_df["SMILES"] = received_df["SMILES"].apply(
                    lambda x: x.split(" ")[0]
                )

                diamond_dfs.append([received_df, new_filename])

            except Exception as e:
                print(f"FAILED ON {csv_file}")
                raise e
                pass

    return diamond_dfs


def create_weizmann_files(received_csv_files):

    weizmann_dfs = []

    for csv_file in received_csv_files:
        if "weizmann" in str(csv_file):
            try:
                received_df = pd.read_csv(csv_file)
                received_df["standard_SMILES"] = received_df["SMILES"].apply(
                    lambda x: strip_and_standardize_smi(x)
                )
                received_df["inchikey"] = received_df["standard_SMILES"].apply(
                    lambda x: Chem.MolToInchiKey(Chem.MolFromSmiles(x))
                )
                received_df["external_ID"] = received_df["inchikey"].apply(
                    lambda x: get_CID(x)
                )
                received_df["CDD_ID"] = received_df["external_ID"].apply(
                    lambda x: get_CDD_ID(x)
                )
                received_df["PostEra_comments"] = received_df["inchikey"].apply(
                    lambda x: get_comments(x)
                )
                if "Batch ID" not in list(received_df.columns):
                    received_df["Batch ID"] = 'UNK'
                new_filename = (
                    str(csv_file)
                    .split("/")[-1]
                    .replace("weizmann", "weizmann_annotated")
                )

                received_df = received_df[
                    [
                        "SMILES",
                        "external_ID",
                        "CDD_ID",
                        "shipment_ID",
                        "catalog_ID",
                        "purity",
                        "volume(uL)",
                        "concentration(mM)",
                        "plate_ID",
                        "well",
                        "stereochemistry",
                        "PostEra_comments",
                        "Batch ID"
                    ]
                ]

                weizmann_dfs.append([received_df, new_filename])

            except Exception as e:
                print(f"FAILED ON {csv_file}")
                raise e
                pass

    return weizmann_dfs


def create_oxford_files(received_csv_files):

    oxford_dfs = []

    for csv_file in received_csv_files:
        if "oxford" in str(csv_file):
            try:
                received_df = pd.read_csv(csv_file)
                received_df["standard_SMILES"] = received_df["SMILES"].apply(
                    lambda x: strip_and_standardize_smi(x)
                )
                received_df["inchikey"] = received_df["standard_SMILES"].apply(
                    lambda x: Chem.MolToInchiKey(Chem.MolFromSmiles(x))
                )
                received_df["external_ID"] = received_df["inchikey"].apply(
                    lambda x: get_CID(x)
                )
                received_df["CDD_ID"] = received_df["external_ID"].apply(
                    lambda x: get_CDD_ID(x)
                )
                received_df["PostEra_comments"] = received_df["inchikey"].apply(
                    lambda x: get_comments(x)
                )
                new_filename = (
                    str(csv_file).split("/")[-1].replace("oxford", "oxford_annotated")
                )
                if "Batch ID" not in list(received_df.columns):
                    received_df["Batch ID"] = 'UNK'

                received_df = received_df[
                    [
                        "SMILES",
                        "external_ID",
                        "CDD_ID",
                        "shipment_ID",
                        "catalog_ID",
                        "purity",
                        "volume(uL)",
                        "concentration(mM)",
                        "plate_ID",
                        "well",
                        "stereochemistry",
                        "PostEra_comments",
                        "Batch ID"
                    ]
                ]

                oxford_dfs.append([received_df, new_filename])

            except Exception as e:
                print(f"FAILED ON {csv_file}")
                raise e
                pass

    return oxford_dfs


if __name__ == "__main__":

    # get parent path of file
    from pathlib import Path

    dir_path = Path(__file__).parent.parent.absolute()

    # get all csvs from folders
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
    shipments_df.to_csv(dir_path / "shipments" / "all_received_mols.csv", index=False)

    diamond_dfs = create_diamond_files(received_csv_files)
    for diamond_df, diamond_fn in diamond_dfs:
        diamond_df.to_csv(
            dir_path / "shipments" / "diamond_files" / diamond_fn, index=False
        )

    weizmann_dfs = create_weizmann_files(received_csv_files)
    for weizmann_df, weizmann_fn in weizmann_dfs:
        weizmann_df.to_csv(
            dir_path / "shipments" / "weizmann_files" / weizmann_fn, index=False,
        )
