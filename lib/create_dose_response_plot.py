# Create plot showing dose-response curves of Weizmann data

# general imports
import numpy as np
import pandas as pd
import altair as alt

from rdkit import Chem
from rdkit.Chem import AllChem

from pathlib import Path

dir_path = Path(__file__).parent.absolute()


def create_dose_response_spec(all_df):
    def get_sigmoid(x):
        return float(x["min"]) + (
            (float(x["max"]) - float(x["min"]))
            / (
                1
                + (float(x["IC50 (µM)"]) / x["Concentration (µM)"])
                ** x["Hill slope"]
            )
        )

    w_ic50_df = all_df.loc[(all_df["f_avg_IC50"].notnull())].drop_duplicates(
        subset="InChIKey"
    )
    w_ic50_df = w_ic50_df.loc[w_ic50_df["f_avg_IC50"] != 99]

    smi_list = []
    ic_50_list = []
    hill_slope_list = []
    conc_list = []
    min_list = []
    max_list = []

    for smi in np.unique(list(w_ic50_df.SMILES)):
        for x in np.logspace(-1, 3, 100):
            smi_list.append(smi)
            ic_50_list.append(
                round(
                    list(w_ic50_df.loc[w_ic50_df.SMILES == smi]["f_avg_IC50"])[
                        0
                    ],
                    3,
                )
            )
            hill_slope_list.append(
                round(
                    list(
                        w_ic50_df.loc[w_ic50_df.SMILES == smi]["f_hill_slope"]
                    )[0],
                    3,
                )
            )
            conc_list.append(x)
            min_list.append(
                min(
                    0,
                    list(
                        w_ic50_df.loc[w_ic50_df.SMILES == smi][
                            "f_min_inhibition_reading"
                        ]
                    )[0],
                )
            )
            max_list.append(
                max(
                    100,
                    list(
                        w_ic50_df.loc[w_ic50_df.SMILES == smi][
                            "f_max_inhibition_reading"
                        ]
                    )[0],
                )
            )

    sigmoid_df = pd.DataFrame(
        {
            "SMILES": smi_list,
            "IC50 (µM)": ic_50_list,
            "Hill slope": hill_slope_list,
            "Concentration (µM)": conc_list,
            "min": min_list,
            "max": max_list,
        }
    )

    sigmoid_df["chloroacetamide"] = sigmoid_df["SMILES"].apply(
        lambda x: True
        if len(
            Chem.MolFromSmiles(x).GetSubstructMatch(
                Chem.MolFromSmarts("Cl[C;H2:1]C(N)=O")
            )
        )
        > 0
        else False
    )
    sigmoid_df = sigmoid_df.loc[sigmoid_df["chloroacetamide"]==False]

    sigmoid_df["R squared"] = sigmoid_df["SMILES"].apply(
        lambda x: list(w_ic50_df.loc[w_ic50_df["SMILES"] == x]["f_R2"])[0]
    )
    sigmoid_df["CID"] = sigmoid_df["SMILES"].apply(
        lambda x: list(
            w_ic50_df.loc[w_ic50_df["SMILES"] == x]["CID (canonical)"]
        )[0]
    )
    sigmoid_df["calc_inh"] = sigmoid_df.apply(get_sigmoid, axis=1)

    selection = alt.selection_multi(fields=["CID"])
    color = alt.condition(
        selection, alt.Color("CID:N", legend=None), alt.value("lightgray")
    )
    line_size = alt.condition(selection, alt.value(5), alt.value(1))

    line_source = sigmoid_df
    lines = (
        alt.Chart(line_source)
        .mark_line(opacity=0.5)
        .encode(
            alt.X(
                "Concentration (µM)",
                scale=alt.Scale(type="log", domain=[0.1, 125]),
            ),
            alt.Y(
                "calc_inh",
                scale=alt.Scale(domain=[-25, 125]),
                axis=alt.Axis(title="Inhibition (%)"),
            ),
            color=color,
            tooltip=["Hill slope", "IC50 (µM)", "R squared", "CID"],
            size=line_size,
        )
    )

    legend = (
        alt.Chart(line_source)
        .mark_point()
        .encode(
            y=alt.Y("CID:N", axis=alt.Axis(orient="right", title="")),
            color=color,
        )
        .add_selection(selection)
    )

    chart = (lines).properties(width=800, height=600).interactive() | legend

    chart = chart.configure_axis(labelFontSize=10, titleFontSize=20)

    chart.save(str(dir_path / "scr" / "weizmann_dose_response_curves.html"))
    chart.save(str(dir_path / "scr" / "weizmann_dose_response_curves.json"))
    with open(
        str(dir_path / "scr" / "weizmann_dose_response_curves.json"), "r"
    ) as f:
        json_data = f.readlines()

    return json_data


def create_dose_response_spec_chloroacetamides(all_df):
    def get_sigmoid(x):
        return float(x["min"]) + (
            (float(x["max"]) - float(x["min"]))
            / (
                1
                + (float(x["IC50 (µM)"]) / x["Concentration (µM)"])
                ** x["Hill slope"]
            )
        )

    w_ic50_df = all_df.loc[(all_df["f_avg_IC50"].notnull())].drop_duplicates(
        subset="InChIKey"
    )
    w_ic50_df = w_ic50_df.loc[w_ic50_df["f_avg_IC50"] != 99]

    smi_list = []
    ic_50_list = []
    hill_slope_list = []
    conc_list = []
    min_list = []
    max_list = []

    for smi in np.unique(list(w_ic50_df.SMILES)):
        for x in np.logspace(-1, 3, 100):
            smi_list.append(smi)
            ic_50_list.append(
                round(
                    list(w_ic50_df.loc[w_ic50_df.SMILES == smi]["f_avg_IC50"])[
                        0
                    ],
                    3,
                )
            )
            hill_slope_list.append(
                round(
                    list(
                        w_ic50_df.loc[w_ic50_df.SMILES == smi]["f_hill_slope"]
                    )[0],
                    3,
                )
            )
            conc_list.append(x)
            min_list.append(
                min(
                    0,
                    list(
                        w_ic50_df.loc[w_ic50_df.SMILES == smi][
                            "f_min_inhibition_reading"
                        ]
                    )[0],
                )
            )
            max_list.append(
                max(
                    100,
                    list(
                        w_ic50_df.loc[w_ic50_df.SMILES == smi][
                            "f_max_inhibition_reading"
                        ]
                    )[0],
                )
            )

    sigmoid_df = pd.DataFrame(
        {
            "SMILES": smi_list,
            "IC50 (µM)": ic_50_list,
            "Hill slope": hill_slope_list,
            "Concentration (µM)": conc_list,
            "min": min_list,
            "max": max_list,
        }
    )
    sigmoid_df["chloroacetamide"] = sigmoid_df["SMILES"].apply(
        lambda x: True
        if len(
            Chem.MolFromSmiles(x).GetSubstructMatch(
                Chem.MolFromSmarts("Cl[C;H2:1]C(N)=O")
            )
        )
        > 0
        else False
    )
    sigmoid_df = sigmoid_df.loc[sigmoid_df["chloroacetamide"]==True]

    sigmoid_df["R squared"] = sigmoid_df["SMILES"].apply(
        lambda x: list(w_ic50_df.loc[w_ic50_df["SMILES"] == x]["f_R2"])[0]
    )
    sigmoid_df["CID"] = sigmoid_df["SMILES"].apply(
        lambda x: list(
            w_ic50_df.loc[w_ic50_df["SMILES"] == x]["CID (canonical)"]
        )[0]
    )
    sigmoid_df["calc_inh"] = sigmoid_df.apply(get_sigmoid, axis=1)

    selection = alt.selection_multi(fields=["CID"])
    color = alt.condition(
        selection, alt.Color("CID:N", legend=None), alt.value("lightgray")
    )
    line_size = alt.condition(selection, alt.value(5), alt.value(1))

    line_source = sigmoid_df
    lines = (
        alt.Chart(line_source)
        .mark_line(opacity=0.5)
        .encode(
            alt.X(
                "Concentration (µM)",
                scale=alt.Scale(type="log", domain=[0.1, 125]),
            ),
            alt.Y(
                "calc_inh",
                scale=alt.Scale(domain=[-25, 125]),
                axis=alt.Axis(title="Inhibition (%)"),
            ),
            color=color,
            tooltip=["Hill slope", "IC50 (µM)", "R squared", "CID"],
            size=line_size,
        )
    )

    legend = (
        alt.Chart(line_source)
        .mark_point()
        .encode(
            y=alt.Y("CID:N", axis=alt.Axis(orient="right", title="")),
            color=color,
        )
        .add_selection(selection)
    )

    chart = (lines).properties(width=800, height=600).interactive() | legend

    chart = chart.configure_axis(labelFontSize=10, titleFontSize=20)

    chart.save(str(dir_path / "scr" / "weizmann_chloroacetamide_dose_response_curves.html"))
    chart.save(str(dir_path / "scr" / "weizmann_chloroacetamide_dose_response_curves.json"))
    with open(
        str(dir_path / "scr" / "weizmann_chloroacetamide_dose_response_curves.json"), "r"
    ) as f:
        json_data = f.readlines()

    return json_data
