# Create plot tracking status of compounds using Altair.

# general imports
import numpy as np
import pandas as pd
import altair as alt

# get parent path of file
from pathlib import Path

dir_path = Path(__file__).parent.absolute()
all_df = pd.read_csv(dir_path / "../covid_submissions_all_info.csv")

num_designed = all_df.shape[0]
num_ordered = all_df[all_df["ORDERED"] == True].shape[0]
num_made = all_df[all_df["MADE"] == True].shape[0]
num_assayed = all_df[all_df["ASSAYED"] == True].shape[0]

tracking_df = pd.DataFrame(
    {
        "stage": ["Designed", "Ordered", "Made", "Assayed"],
        "color": ["Designed", "Ordered", "Made", "Assayed"],
        "num_mols": [num_designed, num_ordered, num_made, num_assayed],
    }
)

# create Altair plot
source = tracking_df
tracking_plot = (
    alt.Chart(source)
    .mark_bar(opacity=0.7)
    .encode(
        y=alt.Y(
            "stage:O",
            sort=["Designed", "Ordered", "Made", "Assayed"],
            title="Stage",
        ),
        x=alt.X("num_mols:Q", stack=None, title="Number of Molecules"),
        color=alt.Color("color", legend=None),
    )
)

tracking_plot.save(
    str(dir_path / "tracking_plot.html"),
    format="html",
    embed_options={"renderer": "svg"},
)
