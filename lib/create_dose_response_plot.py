# Create plot showing dose-response curves of Weizmann data

# general imports
import numpy as np
import pandas as pd
import altair as alt

from rdkit import Chem
from rdkit.Chem import AllChem

from pathlib import Path

dir_path = Path(__file__).parent.absolute()


def create_fluorescence_dose_response_specs(fluorescense_df):
    def get_sigmoid(x):
        return float(x["min"]) + (
            (float(x["max"]) - float(x["min"]))
            / (1 + (float(x["IC50 (µM)"]) / x["Concentration (µM)"]) ** x["Hill slope"])
        )

    for j in range(fluorescense_df.shape[0]):
        f_IC50_dict = fluorescense_df.iloc[j].to_dict()
        CID = f_IC50_dict["CID (canonical)"]
        num_curves = len(f_IC50_dict["f_R2"])

        sigmoids_df = pd.DataFrame({})
        points_df = pd.DataFrame({})
        for i in range(num_curves):
            ic_50_list = []
            hill_slope_list = []
            conc_list = []
            min_list = []
            max_list = []
            r_squared_list = []

            for x in np.logspace(-3, 3, 100):
                ic_50_list.append(round(f_IC50_dict["f_curve_IC50"][i], 3))
                hill_slope_list.append(round(f_IC50_dict["f_hill_slope"][i], 3))
                conc_list.append(x)
                min_list.append(min(0, f_IC50_dict["f_min_inhibition_reading"][i]))
                max_list.append(max(100, f_IC50_dict["f_max_inhibition_reading"][i]))
                r_squared_list.append(f_IC50_dict["f_R2"][i])

            sigmoid_df = pd.DataFrame(
                {
                    "IC50 (µM)": ic_50_list,
                    "Hill slope": hill_slope_list,
                    "Concentration (µM)": conc_list,
                    "min": min_list,
                    "max": max_list,
                    "R squared": r_squared_list,
                    "batch": [(i + 1)] * len(ic_50_list),
                }
            )
            sigmoids_df = pd.concat([sigmoids_df, sigmoid_df], axis=0)

            conc_points = f_IC50_dict["f_concentration_uM"][i]
            inh_points = f_IC50_dict["f_inhibition_list"][i]
            inh_df = pd.DataFrame(
                {
                    "concentration_uM": conc_points,
                    "percent_inhibition": inh_points,
                    "batch": [(i + 1)] * len(conc_points),
                }
            )
            points_df = pd.concat([points_df, inh_df], axis=0)

        sigmoids_df["calc_inh"] = sigmoids_df.apply(get_sigmoid, axis=1)

        selection = alt.selection_multi(fields=["batch"])
        color = alt.condition(
            selection, alt.Color("batch:N", legend=None), alt.value("lightgray")
        )
        line_size = alt.condition(selection, alt.value(5), alt.value(1))

        point_source = points_df
        points = (
            alt.Chart(point_source)
            .mark_circle(opacity=0.5)
            .encode(
                alt.X(
                    "concentration_uM",
                    scale=alt.Scale(type="log", domain=[0.001, 125]),
                    axis=alt.Axis(title="Concentration (µM)"),
                ),
                alt.Y(
                    "percent_inhibition",
                    scale=alt.Scale(domain=[-25, 125]),
                    axis=alt.Axis(title="Inhibition (%)"),
                ),
                color=color,
                size=alt.condition(selection, alt.value(50), alt.value(50)),
            )
        )

        line_source = sigmoids_df
        lines = (
            alt.Chart(line_source)
            .mark_line(opacity=0.5)
            .encode(
                alt.X(
                    "Concentration (µM)",
                    scale=alt.Scale(type="log", domain=[0.001, 125]),
                ),
                alt.Y(
                    "calc_inh",
                    scale=alt.Scale(domain=[-25, 125]),
                    axis=alt.Axis(title="Inhibition (%)"),
                ),
                color=color,
                tooltip=["Hill slope", "IC50 (µM)", "R squared"],
                size=line_size,
            )
        )

        legend = (
            alt.Chart(line_source)
            .mark_point()
            .encode(
                y=alt.Y("batch:N", axis=alt.Axis(orient="right", title="")),
                color=color,
            )
            .add_selection(selection)
        )

        chart = (lines + points).properties(
            width=800, height=600
        ).interactive() | legend

        chart = chart.configure_axis(labelFontSize=10, titleFontSize=20).properties(
            title=f"Dose-Response Curve for {CID} with AVG IC50={round(f_IC50_dict['f_avg_IC50'],3)} µM"
        )
        chart = chart.configure_title(
            fontSize=20, font="Courier", anchor="start", color="gray"
        )

        chart.save(
            str(
                dir_path
                / "scr"
                / "dose_response_curves"
                / f"{CID}_f_dose_response_curve.html"
            )
        )
        chart.save(
            str(
                dir_path
                / "scr"
                / "dose_response_curves"
                / f"{CID}_f_dose_response_curve.json"
            )
        )


def create_rapidfire_dose_response_specs(rapidfire_df):
    def get_sigmoid(x):
        return float(x["min"]) + (
            (float(x["max"]) - float(x["min"]))
            / (1 + (float(x["IC50 (µM)"]) / x["Concentration (µM)"]) ** x["Hill slope"])
        )

    for j in range(rapidfire_df.shape[0]):
        r_IC50_dict = rapidfire_df.iloc[j].to_dict()
        CID = r_IC50_dict["CID (canonical)"]
        num_curves = len(r_IC50_dict["r_R2"])

        sigmoids_df = pd.DataFrame({})
        points_df = pd.DataFrame({})
        for i in range(num_curves):
            ic_50_list = []
            hill_slope_list = []
            conc_list = []
            min_list = []
            max_list = []
            r_squared_list = []

            for x in np.logspace(-3, 3, 100):
                ic_50_list.append(round(r_IC50_dict["r_curve_IC50"][i], 3))
                hill_slope_list.append(round(r_IC50_dict["r_hill_slope"][i], 3))
                conc_list.append(x)
                min_list.append(min(0, r_IC50_dict["r_min_inhibition_reading"][i]))
                max_list.append(max(100, r_IC50_dict["r_max_inhibition_reading"][i]))
                r_squared_list.append(r_IC50_dict["r_R2"][i])

            sigmoid_df = pd.DataFrame(
                {
                    "IC50 (µM)": ic_50_list,
                    "Hill slope": hill_slope_list,
                    "Concentration (µM)": conc_list,
                    "min": min_list,
                    "max": max_list,
                    "R squared": r_squared_list,
                    "batch": [(i + 1)] * len(ic_50_list),
                }
            )
            sigmoids_df = pd.concat([sigmoids_df, sigmoid_df], axis=0)

            conc_points = r_IC50_dict["r_concentration_uM"][i]
            inh_points = r_IC50_dict["r_inhibition_list"][i]
            inh_df = pd.DataFrame(
                {
                    "concentration_uM": conc_points,
                    "percent_inhibition": inh_points,
                    "batch": [(i + 1)] * len(conc_points),
                }
            )
            points_df = pd.concat([points_df, inh_df], axis=0)

        sigmoids_df["calc_inh"] = sigmoids_df.apply(get_sigmoid, axis=1)

        selection = alt.selection_multi(fields=["batch"])
        color = alt.condition(
            selection, alt.Color("batch:N", legend=None), alt.value("lightgray")
        )
        line_size = alt.condition(selection, alt.value(5), alt.value(1))

        point_source = points_df
        points = (
            alt.Chart(point_source)
            .mark_circle(opacity=0.5)
            .encode(
                alt.X(
                    "concentration_uM",
                    scale=alt.Scale(type="log", domain=[0.001, 125]),
                    axis=alt.Axis(title="Concentration (µM)"),
                ),
                alt.Y(
                    "percent_inhibition",
                    scale=alt.Scale(domain=[-25, 125]),
                    axis=alt.Axis(title="Inhibition (%)"),
                ),
                color=color,
                size=alt.condition(selection, alt.value(50), alt.value(50)),
            )
        )

        line_source = sigmoids_df
        lines = (
            alt.Chart(line_source)
            .mark_line(opacity=0.5)
            .encode(
                alt.X(
                    "Concentration (µM)",
                    scale=alt.Scale(type="log", domain=[0.001, 125]),
                ),
                alt.Y(
                    "calc_inh",
                    scale=alt.Scale(domain=[-25, 125]),
                    axis=alt.Axis(title="Inhibition (%)"),
                ),
                color=color,
                tooltip=["Hill slope", "IC50 (µM)", "R squared"],
                size=line_size,
            )
        )

        legend = (
            alt.Chart(line_source)
            .mark_point()
            .encode(
                y=alt.Y("batch:N", axis=alt.Axis(orient="right", title="")),
                color=color,
            )
            .add_selection(selection)
        )

        chart = (lines + points).properties(
            width=800, height=600
        ).interactive() | legend

        chart = chart.configure_axis(labelFontSize=10, titleFontSize=20).properties(
            title=f"Dose-Response Curve for {CID} with AVG IC50={round(r_IC50_dict['r_avg_IC50'],3)} µM"
        )
        chart = chart.configure_title(
            fontSize=20, font="Courier", anchor="start", color="gray"
        )

        chart.save(
            str(
                dir_path
                / "scr"
                / "dose_response_curves"
                / f"{CID}_r_dose_response_curve.html"
            )
        )
        chart.save(
            str(
                dir_path
                / "scr"
                / "dose_response_curves"
                / f"{CID}_r_dose_response_curve.json"
            )
        )

