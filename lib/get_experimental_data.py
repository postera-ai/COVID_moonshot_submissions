# Matt Robinson, matthew.robinson@postera.ai
# Compile all experimental info from CDD

# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

import requests
import sys
import time
import json

# get parent path of file
from pathlib import Path

lib_path = Path(__file__).parent.absolute()
with open(lib_path / "CDD.env", "r") as env_file:
    env_vars_dict = dict(
        tuple(line.split("="))
        for line in [x.strip("\n") for x in env_file.readlines()]
        if not line.startswith("#")
    )

vault_num = env_vars_dict["VAULT_num"]
vault_token = env_vars_dict["VAULT_token"]

virtual_project_id = "12336"
synthesis_project_id = "12334"
made_project_id = "12335"

rapidfire_IC50_protocol_id = "49700"
rapidfire_inhibition_protocol_id = "49192"
fluorescence_IC50_protocol_id = "49439"
fluorescence_inhibition_protocol_id = "49412"
solubility_protocol_id = "49275"
trypsin_protocol_id = "49443"


def get_async_export(async_url):
    headers = {"X-CDD-token": vault_token}
    response = requests.get(async_url, headers=headers)
    print("BEGINNING EXPORT")
    print(response)
    export_info = response.json()

    export_id = export_info["id"]

    # CHECK STATUS of Export
    status = None
    seconds_waiting = 0

    while status != "finished":
        headers = {"X-CDD-token": vault_token}
        url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_num}/export_progress/{export_id}"

        response = requests.get(url, headers=headers)

        print("CHECKING STATUS of EXPORT:")
        # to view the status, use:
        print(response)
        status = response.json()["status"]
        print(status)

        time.sleep(5)
        seconds_waiting += 5
        if seconds_waiting > 200:
            print("Export Never Finished")
            break

    if status != "finished":
        sys.exit("EXPORT IS BROKEN")

    headers = {"X-CDD-token": vault_token}
    url = (
        f"https://app.collaborativedrug.com/api/v1/vaults/{vault_num}/exports/{export_id}"
    )

    print("RETRIEVING FINISHED EXPORT")
    response = requests.get(url, headers=headers)
    print(response)
    return response


def get_rapidfire_inhibition_data():

    url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_num}/protocols/{rapidfire_inhibition_protocol_id}/data?async=True"
    response = get_async_export(url)
    inhibition_response_dict = response.json()["objects"]

    with open(lib_path / "scr" / "rapidfire_inhibition_data.json", "w") as f:
        json.dump(inhibition_response_dict, f)

    inhibition_data_dict = {}
    for mol_dict in inhibition_response_dict:
        if "molecule" not in mol_dict:
            continue
        mol_id = mol_dict["molecule"]
        if mol_dict["readouts"]["553839"] == 20.0:
            if mol_id not in inhibition_data_dict:
                inhibition_data_dict[mol_id] = {
                    "20_uM": mol_dict["readouts"]["553894"]["value"]
                }
            else:
                inhibition_data_dict[mol_id]["20_uM"] = mol_dict["readouts"][
                    "553894"
                ]["value"]

        elif mol_dict["readouts"]["553839"] == 50.0:
            if mol_id not in inhibition_data_dict:
                inhibition_data_dict[mol_id] = {
                    "50_uM": mol_dict["readouts"]["553894"]["value"]
                }
            else:
                inhibition_data_dict[mol_id]["50_uM"] = mol_dict["readouts"][
                    "553894"
                ]["value"]

    mol_id_list = [float(x) for x in inhibition_data_dict.keys()]
    for mol_id in mol_id_list:
        if "20_uM" not in inhibition_data_dict[mol_id]:
            inhibition_data_dict[mol_id]["20_uM"] = np.nan
        if "50_uM" not in inhibition_data_dict[mol_id]:
            inhibition_data_dict[mol_id]["50_uM"] = np.nan
    inhibition_at_20_uM_list = [
        inhibition_data_dict[mol_id]["20_uM"] for mol_id in mol_id_list
    ]
    inhibition_at_50_uM_list = [
        inhibition_data_dict[mol_id]["50_uM"] for mol_id in mol_id_list
    ]

    inhibition_df = pd.DataFrame(
        {
            "CDD_mol_ID": mol_id_list,
            "r_inhibition_at_20_uM": inhibition_at_20_uM_list,
            "r_inhibition_at_50_uM": inhibition_at_50_uM_list,
        }
    )
    inhibition_df = inhibition_df.drop_duplicates(subset="CDD_mol_ID")
    return inhibition_df


def get_rapidfire_IC50_data():

    url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_num}/protocols/{rapidfire_IC50_protocol_id}/data?async=True"
    response = get_async_export(url)
    rapid_fire_dose_response_dict = response.json()["objects"]

    with open(lib_path / "scr" / "rapidfire_IC50_data.json", "w") as f:
        json.dump(rapid_fire_dose_response_dict, f)

    mol_id_list = []
    ic50_list = []

    for mol_dict in rapid_fire_dose_response_dict:
        if "molecule" not in mol_dict:
            continue
        mol_id = mol_dict["molecule"]
        if "564286" in mol_dict["readouts"]:
            if type(mol_dict["readouts"]["564286"]) == dict:
                ic50 = 99
            else:
                ic50 = mol_dict["readouts"]["564286"]
        else:
            ic50 = np.nan

        mol_id_list.append(float(mol_id))
        ic50_list.append(ic50)

    rapidfire_df = pd.DataFrame(
        {"CDD_mol_ID": mol_id_list, "r_IC50": ic50_list}
    )
    rapidfire_df = rapidfire_df.drop_duplicates(subset="CDD_mol_ID")
    return rapidfire_df


def get_fluorescense_inhibition_data():

    url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_num}/protocols/{fluorescence_inhibition_protocol_id}/data?async=True"
    response = get_async_export(url)
    inhibition_response_dict = response.json()["objects"]

    with open(lib_path / "scr" / "fluorescense_inhibition_data.json", "w") as f:
        json.dump(inhibition_response_dict, f)

    inhibition_data_dict = {}
    for mol_dict in inhibition_response_dict:
        if "molecule" not in mol_dict:
            continue
        mol_id = mol_dict["molecule"]
        if mol_dict["readouts"]["556717"] == 20.0:
            if mol_id not in inhibition_data_dict:
                inhibition_data_dict[mol_id] = {
                    "20_uM": mol_dict["readouts"]["556718"]["value"]
                }
            else:
                inhibition_data_dict[mol_id]["20_uM"] = mol_dict["readouts"][
                    "556718"
                ]["value"]

        elif mol_dict["readouts"]["556717"] == 50.0:
            if mol_id not in inhibition_data_dict:
                inhibition_data_dict[mol_id] = {
                    "50_uM": mol_dict["readouts"]["556718"]["value"]
                }
            else:
                inhibition_data_dict[mol_id]["50_uM"] = mol_dict["readouts"][
                    "556718"
                ]["value"]

    mol_id_list = [float(x) for x in inhibition_data_dict.keys()]
    for mol_id in mol_id_list:
        if "20_uM" not in inhibition_data_dict[mol_id]:
            inhibition_data_dict[mol_id]["20_uM"] = np.nan
        if "50_uM" not in inhibition_data_dict[mol_id]:
            inhibition_data_dict[mol_id]["50_uM"] = np.nan

    inhibition_at_20_uM_list = [
        inhibition_data_dict[mol_id]["20_uM"] for mol_id in mol_id_list
    ]
    inhibition_at_50_uM_list = [
        inhibition_data_dict[mol_id]["50_uM"] for mol_id in mol_id_list
    ]

    inhibition_df = pd.DataFrame(
        {
            "CDD_mol_ID": mol_id_list,
            "f_inhibition_at_20_uM": inhibition_at_20_uM_list,
            "f_inhibition_at_50_uM": inhibition_at_50_uM_list,
        }
    )
    inhibition_df = inhibition_df.drop_duplicates(subset="CDD_mol_ID")
    return inhibition_df


def get_fluorescense_IC50_data():

    url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_num}/protocols/{fluorescence_IC50_protocol_id}/data?async=True"
    response = get_async_export(url)
    fluorescence_response_dict = response.json()["objects"]

    with open(lib_path / "scr" / "fluorescense_IC50_data.json", "w") as f:
        json.dump(fluorescence_response_dict, f)

    mol_id_list = []
    avg_ic50_list = []
    avg_pic50_list = []
    max_reading_list = []
    min_reading_list = []
    hill_slope_list = []
    r2_list = []

    for mol_dict in fluorescence_response_dict:
        if "molecule" not in mol_dict:
            continue

        elif mol_dict["molecule"] in mol_id_list:
            continue

        mol_id = mol_dict["molecule"]
        if "557736" in mol_dict["readouts"]:
            avg_ic50 = mol_dict["readouts"]["557736"]["value"]
        else:
            avg_ic50 = np.nan

        if "557738" in mol_dict["readouts"]:
            if type(mol_dict["readouts"]["557738"]) == dict:
                avg_pic50 = np.nan
            else:
                avg_pic50 = mol_dict["readouts"]["557738"]
        else:
            avg_pic50 = np.nan

        if "557085" in mol_dict["readouts"]:
            min_reading = mol_dict["readouts"]["557085"]["value"]
        else:
            min_reading = np.nan

        if "557086" in mol_dict["readouts"]:
            max_reading = mol_dict["readouts"]["557086"]["value"]
        else:
            max_reading = np.nan

        if "557078" in mol_dict["readouts"]:
            hill_slope = mol_dict["readouts"]["557078"]
        else:
            hill_slope = np.nan

        if "557082" in mol_dict["readouts"]:
            r2 = mol_dict["readouts"]["557082"]
        else:
            r2 = np.nan

        mol_id_list.append(float(mol_id))
        avg_ic50_list.append(avg_ic50)
        avg_pic50_list.append(avg_pic50)
        max_reading_list.append(max_reading)
        min_reading_list.append(min_reading)
        hill_slope_list.append(hill_slope)
        r2_list.append(r2)

    fluorescence_df = pd.DataFrame(
        {
            "CDD_mol_ID": mol_id_list,
            "f_avg_IC50": avg_ic50_list,
            "f_avg_pIC50": avg_pic50_list,
            "f_max_inhibition_reading": max_reading_list,
            "f_min_inhibition_reading": min_reading_list,
            "f_hill_slope": hill_slope_list,
            "f_R2": r2_list,
        }
    )
    fluorescence_df = fluorescence_df.drop_duplicates(subset="CDD_mol_ID")
    return fluorescence_df


def get_solubility_data():

    url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_num}/protocols/{solubility_protocol_id}/data?async=True"
    response = get_async_export(url)
    solubility_response_dict = response.json()["objects"]

    with open(lib_path / "scr" / "solubility_data.json", "w") as f:
        json.dump(solubility_response_dict, f)

    solubility_data_dict = {}
    for mol_dict in solubility_response_dict:
        mol_id = mol_dict["molecule"]
        if mol_dict["readouts"]["554984"] == 20.0:
            if mol_id not in solubility_data_dict:
                solubility_data_dict[mol_id] = {
                    "20_uM": mol_dict["readouts"]["555388"]["value"]
                }
            else:
                solubility_data_dict[mol_id]["20_uM"] = mol_dict["readouts"][
                    "555388"
                ]["value"]

        elif mol_dict["readouts"]["554984"] == 100.0:
            if mol_id not in solubility_data_dict:
                solubility_data_dict[mol_id] = {
                    "100_uM": mol_dict["readouts"]["555388"]["value"]
                }
            else:
                solubility_data_dict[mol_id]["100_uM"] = mol_dict["readouts"][
                    "555388"
                ]["value"]

    mol_id_list = [float(x) for x in solubility_data_dict.keys()]
    relative_solubility_at_20_uM_list = [
        solubility_data_dict[mol_id]["20_uM"] for mol_id in mol_id_list
    ]
    relative_solubility_at_100_uM_list = [
        solubility_data_dict[mol_id]["100_uM"] for mol_id in mol_id_list
    ]

    solubility_df = pd.DataFrame(
        {
            "CDD_mol_ID": mol_id_list,
            "relative_solubility_at_20_uM": relative_solubility_at_20_uM_list,
            "relative_solubility_at_100_uM": relative_solubility_at_100_uM_list,
        }
    )
    solubility_df = solubility_df.drop_duplicates(subset="CDD_mol_ID")
    return solubility_df


def get_trypsin_data():

    url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_num}/protocols/{trypsin_protocol_id}/data?async=True"
    response = get_async_export(url)
    trypsin_response_dict = response.json()["objects"]

    with open(lib_path / "scr" / "trypsin_data.json", "w") as f:
        json.dump(trypsin_response_dict, f)

    mol_id_list = []
    ic50_list = []

    for mol_dict in trypsin_response_dict:
        if "molecule" not in mol_dict:
            continue
        mol_id = mol_dict["molecule"]
        if mol_id in mol_id_list:
            continue

        if "557122" in mol_dict["readouts"]:
            if type(mol_dict["readouts"]["557122"]) == float:
                ic50 = mol_dict["readouts"]["557122"]
            else:
                if "value" not in mol_dict["readouts"]["557122"]:
                    ic50 = np.nan
                else:
                    ic50 = 99  # mol_dict["readouts"]["557122"]["value"]
        else:
            ic50 = np.nan

        mol_id_list.append(float(mol_id))
        ic50_list.append(ic50)

    trypsin_df = pd.DataFrame(
        {"CDD_mol_ID": mol_id_list, "trypsin_IC50": ic50_list}
    )
    trypsin_df = trypsin_df.drop_duplicates(subset="CDD_mol_ID")
    return trypsin_df
