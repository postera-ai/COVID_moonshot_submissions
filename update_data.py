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


def update_data(
    fetch_submissions=True,
    fetch_orders=True,
    fetch_shipments=True,
    fetch_CDD=True,
    get_CDD_updates=True,
    fetch_assays=True,
    fetch_structures=True,
    update_tracking_status=True,
    update_plots=True,
):

    ### First Get the data ###
    if fetch_submissions:
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
        all_df = all_df.drop(columns=[
            "r_inhibition_at_20_uM",
            "r_inhibition_at_50_uM",
            "r_IC50",
            "f_inhibition_at_20_uM",
            "f_inhibition_at_50_uM",
            "f_avg_IC50",
            "f_avg_pIC50",
            "f_max_inhibition_reading",
            "f_min_inhibition_reading",
            "f_hill_slope",
            "f_R2",
            "relative_solubility_at_20_uM",
            "relative_solubility_at_100_uM",
            "trypsin_IC50"
        ])

        def create_old_cid(x):
            if x["old_CID"] is np.nan:
                return x["CID"]
            else:
                return x["old_CID"]

        id_df = all_df.copy()[
            ["SMILES", "CID", "CID (canonical)", "CID (old format)"]
        ]
        id_df = id_df.rename(
            columns={
                "CID (old format)": "old_CID",
                "CID (canonical)": "canonical_CID",
            }
        )
        id_df["inchikey"] = id_df["SMILES"].apply(
            lambda x: Chem.MolToInchiKey(Chem.MolFromSmiles(x))
        )
        id_df["short_inchikey"] = id_df["inchikey"].apply(
            lambda x: x.split("-")[0]
        )
        id_df["old_CID"] = id_df.apply(create_old_cid, axis=1)
        id_df.to_csv(dir_path / "covid_moonshot_ids.csv", index=False)
    else:
        all_df = pd.read_csv(dir_path / "covid_submissions_all_info.csv")

    ### Update the orders data
    from lib.get_all_ordered_mols import update_orders_data

    if fetch_orders:
        order_csv_files = [
            f
            for f in dir_path.glob("orders/**/*.csv")
            if "all_ordered_mols.csv" not in str(f)
        ]
        orders_df = update_orders_data(order_csv_files)
        orders_df.to_csv(
            dir_path / "orders" / "all_ordered_mols.csv", index=False
        )
    else:
        orders_df = pd.read_csv(dir_path / "orders" / "all_ordered_mols.csv")

    ### Update shipments data
    from lib.get_all_received_mols import (
        update_shipments_data,
        create_diamond_files,
        create_weizmann_files,
        create_oxford_files,
    )

    if fetch_shipments:
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
                dir_path / "shipments" / "diamond_files" / diamond_fn,
                index=False,
            )

        weizmann_dfs = create_weizmann_files(received_csv_files)
        for weizmann_df, weizmann_fn in weizmann_dfs:
            weizmann_df.to_csv(
                dir_path / "shipments" / "weizmann_files" / weizmann_fn,
                index=False,
            )

        oxford_dfs = create_oxford_files(received_csv_files)
        for oxford_df, oxford_fn in oxford_dfs:
            oxford_df.to_csv(
                dir_path / "shipments" / "oxford_files" / oxford_fn,
                index=False,
            )
    else:
        shipments_df = pd.read_csv(
            dir_path / "shipments" / "all_received_mols.csv"
        )

    ### Update CDD info

    # first get the necessary data
    received_df = pd.read_csv(dir_path / "shipments" / "all_received_mols.csv")
    made_df = received_df.copy()
    made_df.to_csv(
        dir_path / "data_for_CDD" / "compounds" / "Compounds_Made.csv",
        index=False,
    )

    ordered_df = pd.read_csv(dir_path / "orders" / "all_ordered_mols.csv")
    synthesis_df = ordered_df.copy()
    synthesis_df.to_csv(
        dir_path
        / "data_for_CDD"
        / "compounds"
        / "Compounds_for_Synthesis.csv",
        index=False,
    )

    virtual_df = all_df.copy()
    virtual_df.to_csv(
        dir_path / "data_for_CDD" / "compounds" / "Compounds_Virtual.csv",
        index=False,
    )

    # get current data in the vault
    from lib.get_current_vault_data import get_current_vault_data

    if fetch_CDD:
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
    else:
        current_cdd_df = pd.read_csv(
            dir_path
            / "data_for_CDD"
            / "current_vault_data"
            / "current_vault_data.csv"
        )

    # get the necessary updates to CDD
    from lib.get_CDD_updates import get_CDD_updates

    if get_CDD_updates:
        add_to_virtual_df, add_to_synthesis_df, add_to_made_df = get_CDD_updates(
            all_df, current_cdd_df, virtual_df, synthesis_df, made_df
        )

        add_to_virtual_df.to_csv(
            dir_path
            / "data_for_CDD"
            / "vault_updates"
            / "add_to_virtual_df.csv",
            index=False,
        )
        add_to_synthesis_df.to_csv(
            dir_path
            / "data_for_CDD"
            / "vault_updates"
            / "add_to_synthesis_df.csv",
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

    if fetch_assays:
        rapidfire_inhibition_df = get_rapidfire_inhibition_data()
        rapidfire_IC50_df = get_rapidfire_IC50_data()

        fluorescence_inhibition_df = get_fluorescense_inhibition_data()
        fluorescence_IC50_df = get_fluorescense_IC50_data()

        solubility_df = get_solubility_data()
        trypsin_df = get_trypsin_data()

        all_df = pd.merge(
            all_df, rapidfire_inhibition_df, how="left", on=["CDD_mol_ID"]
        )
        all_df = pd.merge(
            all_df, rapidfire_IC50_df, how="left", on=["CDD_mol_ID"]
        )
        all_df = pd.merge(
            all_df, fluorescence_inhibition_df, how="left", on=["CDD_mol_ID"]
        )
        all_df = pd.merge(
            all_df, fluorescence_IC50_df, how="left", on=["CDD_mol_ID"]
        )
        all_df = pd.merge(all_df, solubility_df, how="left", on=["CDD_mol_ID"])
        all_df = pd.merge(all_df, trypsin_df, how="left", on=["CDD_mol_ID"])

    if fetch_structures:
        # update structural info
        structures_df = pd.read_csv(
            dir_path / "structures" / "fragalysis_structures.csv"
        )
        structures_df["InChIKey"] = structures_df["SMILES"].apply(
            lambda x: Chem.MolToInchiKey(Chem.MolFromSmiles(x))
        )
        structures_df = structures_df[
            ["InChIKey", "structure_ID", "structure_LINK"]
        ]
        all_df = pd.merge(all_df, structures_df, how="left", on=["InChIKey"])

    if update_tracking_status:
        ### only list things with at least one inhibition value as assayed
        ordered_iks = [
            Chem.MolToInchiKey(Chem.MolFromSmiles(x))
            for x in list(synthesis_df["SMILES"])
        ]
        made_iks = [
            Chem.MolToInchiKey(Chem.MolFromSmiles(x))
            for x in list(made_df["SMILES"])
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
            lambda x: "TRUE"
            if ((x in made_iks) or (x in assayed_iks))
            else "FALSE"
        )
        all_df["ASSAYED"] = all_df["InChIKey"].apply(
            lambda x: "TRUE" if x in assayed_iks else "FALSE"
        )

    all_df.to_csv(dir_path / "covid_submissions_all_info.csv", index=False)

    if update_plots:
        from lib.create_tracking_plot import create_tracking_plot_spec

        tracking_plot_spec_data = create_tracking_plot_spec(all_df)
        with open(
            dir_path / "plots" / "tracking_plot_vega_spec.json", "w"
        ) as f:
            f.writelines(tracking_plot_spec_data)

        from lib.create_pIC50_plot import create_pIC50_html

        pIC50_html_data = create_pIC50_html(all_df)
        with open(dir_path / "plots" / "pIC50_plot.html", "w") as f:
            f.writelines(pIC50_html_data)

        from lib.create_dose_response_plot import (
            create_dose_response_spec,
            create_dose_response_spec_chloroacetamides,
        )

        dose_response_spec_data = create_dose_response_spec(all_df)
        with open(
            dir_path / "plots" / "dose_response_vega_spec.json", "w"
        ) as f:
            f.writelines(dose_response_spec_data)

        chloroacetamides_dose_response_spec_data = create_dose_response_spec_chloroacetamides(
            all_df
        )
        with open(
            dir_path
            / "plots"
            / "chloroacetamides_dose_response_vega_spec.json",
            "w",
        ) as f:
            f.writelines(chloroacetamides_dose_response_spec_data)


if __name__ == "__main__":
    update_data(
        fetch_submissions=True,
        fetch_orders=True,
        fetch_shipments=True,
        fetch_CDD=True,
        get_CDD_updates=True,
        fetch_assays=True,
        fetch_structures=True,
        update_tracking_status=True,
        update_plots=True,
    )