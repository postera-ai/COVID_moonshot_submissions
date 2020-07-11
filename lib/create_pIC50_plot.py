# Create plot showing comparison of IC50 values
# general imports
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

import bokeh
from bokeh.models import HoverTool
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.transform import factor_cmap
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.embed import json_item

import json
from pathlib import Path

dir_path = Path(__file__).parent.absolute()


def create_pIC50_html_and_json(all_df):

    ic50_df = all_df.loc[
        (all_df["f_avg_IC50"].notnull()) & (all_df["r_avg_IC50"].notnull())
    ].drop_duplicates(subset="InChIKey")

    mol_imgs = []
    for idx, smi in enumerate(list(ic50_df["SMILES"])):
        png_file_name = str(
            dir_path
            / "scr"
            / "mol_images"
            / (Chem.MolToInchiKey(Chem.MolFromSmiles(smi)) + ".png")
        )
        Draw.MolToFile(Chem.MolFromSmiles(smi), png_file_name)
        mol_imgs.append(png_file_name)
    ic50_df["mol_imgs"] = mol_imgs

    ic50_df["r_pIC50"] = ic50_df["r_avg_IC50"].apply(
        lambda x: -1.0 * np.log10(x * 1.0e-6)
    )
    ic50_df["f_pIC50"] = ic50_df["f_avg_IC50"].apply(
        lambda x: -1.0 * np.log10(x * 1.0e-6)
    )

    source = ColumnDataSource(
        data=dict(
            x=ic50_df["f_pIC50"],
            y=ic50_df["r_pIC50"],
            desc1=ic50_df["CID (canonical)"],
            imgs=ic50_df["mol_imgs"],
        )
    )

    TOOLTIPS = """
        <div>
            <div>
                <img
                    src="@imgs" height="150" alt="@imgs" width="150"
                    style="float: left; margin: 0px 15px 15px 0px;"
                    border="1"
                ></img>
                </br>
            </div>
            <div>
                <span style="font-size: 12px;">@desc1</span>
            </div>
            <div>
                <span style="font-size: 12px;">Values</span></br>
                <span style="font-size: 12px; color: #696;">(Fluorescense pIC50: $x)</span></br>
                <span style="font-size: 12px; color: #696;">(RapidFire pIC50: $y)
            </div>
        </div>
    """

    p = figure(
        plot_width=800,
        plot_height=620,
        tooltips=TOOLTIPS,
        x_axis_label="Fluorescense pIC50",
        y_axis_label="RapidFire pIC50",
        title="Comparing pIC50 values between assays",
    )

    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"
    p.circle("x", "y", size=10, source=source, color="red", alpha=0.5)

    html = file_html(p, CDN, "Comparing pIC50 values between assays")
    # return html
    return html, json.dumps(json_item(p, "Comparing pIC50 values between assays"))

