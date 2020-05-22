# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem


def get_CDD_updates(all_df, current_cdd_df, virtual_df, synthesis_df, made_df):

    cdd_virtual_external_id_list = list(
        current_cdd_df.loc[current_cdd_df["virtual_project"] == True][
            "external_ID"
        ]
    )
    cdd_synthesis_external_id_list = list(
        current_cdd_df.loc[current_cdd_df["for_synthesis_project"] == True][
            "external_ID"
        ]
    )
    cdd_made_external_id_list = list(
        current_cdd_df.loc[current_cdd_df["made_project"] == True][
            "external_ID"
        ]
    )

    ### FILTER
    add_to_virtual_df = virtual_df.loc[
        ~virtual_df["CID"].isin(cdd_virtual_external_id_list)
    ]
    add_to_virtual_df = add_to_virtual_df.reset_index(drop=True)
    add_to_virtual_df = add_to_virtual_df[
        [
            "SMILES",
            "CID",
            "creator",
            "fragments",
            "covalent_warhead",
            "CID (canonical)",
        ]
    ]

    add_to_synthesis_df = synthesis_df.loc[
        ~synthesis_df["CID"].isin(cdd_synthesis_external_id_list)
    ]
    add_to_synthesis_df = add_to_synthesis_df.reset_index(drop=True)
    add_to_synthesis_df = add_to_synthesis_df[["SMILES", "CID"]]

    add_to_made_df = made_df.loc[
        ~made_df["CID"].isin(cdd_made_external_id_list)
    ]
    add_to_made_df = add_to_made_df.reset_index(drop=True)
    add_to_made_df = add_to_made_df[["SMILES", "CID", "comments"]]

    extra_cids_to_add_to_synthesis = []
    extra_cids_to_add_to_made = []

    for cid in list(add_to_synthesis_df["CID"]):
        CDD_name = all_df.loc[all_df["CID"] == cid]["CDD_name"].item()
        if len(list(all_df.loc[all_df["CDD_name"] == CDD_name])) > 1:
            cids_to_add = list(
                all_df.loc[all_df["CDD_name"] == CDD_name]["CID"]
            )
            cids_to_add.remove(cid)
            extra_cids_to_add_to_synthesis.extend(cids_to_add)

    for cid in list(add_to_made_df["CID"]):
        CDD_name = all_df.loc[all_df["CID"] == cid]["CDD_name"].item()
        if len(list(all_df.loc[all_df["CDD_name"] == CDD_name])) > 1:
            cids_to_add = list(
                all_df.loc[all_df["CDD_name"] == CDD_name]["CID"]
            )
            cids_to_add.remove(cid)
            extra_cids_to_add_to_made.extend(cids_to_add)

    # # mapping all the extra ids, which are batches of those mols,
    # # slight delay of one round for these
    # add_to_synthesis = []
    # add_to_made = []
    # for cdd_id in list(current_cdd_df["CDD_name"]):
    #     for_syn_list = list(
    #         current_cdd_df[current_cdd_df["CDD_name"] == cdd_id][
    #             "for_synthesis_project"
    #         ]
    #     )
    #     made_list = list(
    #         current_cdd_df[current_cdd_df["CDD_name"] == cdd_id][
    #             "made_project"
    #         ]
    #     )
    #     if (
    #         (len(for_syn_list) > 1)
    #         and (True in for_syn_list)
    #         and (False in for_syn_list)
    #     ):
    #         add_to_synthesis.extend(
    #             list(
    #                 current_cdd_df[
    #                     (current_cdd_df["CDD_name"] == cdd_id)
    #                     & (current_cdd_df["for_synthesis_project"] == False)
    #                 ]["external_ID"]
    #             )
    #         )
    #     if (
    #         (len(made_list) > 1)
    #         and (True in made_list)
    #         and (False in made_list)
    #     ):
    #         add_to_made.extend(
    #             list(
    #                 current_cdd_df[
    #                     (current_cdd_df["CDD_name"] == cdd_id)
    #                     & (current_cdd_df["made_project"] == False)
    #                 ]["external_ID"]
    #             )
    #         )
    add_to_made = list(set(extra_cids_to_add_to_made))
    add_to_synthesis = list(set(extra_cids_to_add_to_synthesis))
    add_to_synthesis_df = pd.concat(
        [
            add_to_synthesis_df,
            pd.DataFrame(
                {
                    "SMILES": [""] * len(add_to_synthesis),
                    "CID": add_to_synthesis,
                }
            ),
        ],
        axis=0,
    )
    add_to_made_df = pd.concat(
        [
            add_to_made_df,
            pd.DataFrame(
                {
                    "SMILES": [""] * len(add_to_made),
                    "CID": add_to_made,
                    "comments": [""] * len(add_to_made),
                }
            ),
        ],
        axis=0,
    )

    return add_to_virtual_df, add_to_synthesis_df, add_to_made_df
