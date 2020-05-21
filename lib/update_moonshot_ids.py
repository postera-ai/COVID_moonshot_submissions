### update the ids so that they have a canonical CID

# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

from .utils import (
    strip_and_standardize_smi,
    get_new_CID_from_old,
    get_old_CID_from_new,
)


def update_CIDs(all_df, id_df):

    existing_inchikeys_list = list(id_df["inchikey"])
    existing_cids_list = list(id_df["CID"])

    new_smi_list = []
    new_cids_list = []
    old_cids_list = []
    canonical_cids_list = []
    inchikeys_list = []
    short_inchikeys_list = []

    for smi, cid in zip(list(all_df["SMILES"]), list(all_df["CID"])):
        if cid in existing_cids_list:
            continue

        inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
        if inchikey is None:
            print(f"Could not produce inchikey for {cid} with SMILES {smi}")
            continue
        if inchikey in existing_inchikeys_list:
            new_cids_list.append(cid)
            inchikeys_list.append(inchikey)
            short_inchikeys_list.append(inchikey.split("-")[0])
            new_smi_list.append(smi)
            old_cids_list.append(get_old_CID_from_new(cid))
            canonical_cids_list.append(
                list(
                    id_df.loc[id_df["inchikey"] == inchikey]["canonical_CID"]
                )[0]
            )

        else:
            new_cids_list.append(cid)
            inchikeys_list.append(inchikey)
            short_inchikeys_list.append(inchikey.split("-")[0])
            new_smi_list.append(smi)
            old_cids_list.append(get_old_CID_from_new(cid))
            canonical_cids_list.append(
                list(all_df.loc[all_df["CID"] == cid]["CID (canonical)"])[0]
            )

    new_id_df = pd.DataFrame(
        {
            "SMILES": new_smi_list,
            "CID": new_cids_list,
            "old_CID": old_cids_list,
            "canonical_CID": canonical_cids_list,
            "inchikey": inchikeys_list,
            "short_inchikey": short_inchikeys_list,
        }
    )
    new_id_df = new_id_df.loc[~(new_id_df["inchikey"] == "")]

    id_df = pd.concat([id_df.copy(), new_id_df], axis=0)
    id_df = id_df.sort_values(by=["CID"]).reset_index(drop=True)
    return id_df


if __name__ == "__main__":

    # need to update all_df
    # get parent path of file
    from pathlib import Path

    dir_path = Path(__file__).parent.parent.absolute()

    all_df = pd.read_csv(dir_path / "covid_submissions_all_info.csv")
    all_df["SMILES"] = all_df["SMILES"].apply(
        lambda x: strip_and_standardize_smi(x)
    )
    all_df.to_csv(dir_path / "covid_submissions_all_info.csv", index=False)
    id_df = pd.read_csv(dir_path / "covid_moonshot_ids.csv")
    update_CIDs(all_df, id_df)
