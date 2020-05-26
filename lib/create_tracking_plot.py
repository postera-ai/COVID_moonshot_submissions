# Create plot tracking status of compounds using Altair.

# general imports
import numpy as np
import pandas as pd
import altair as alt

from pathlib import Path
dir_path = Path(__file__).parent.absolute()


def create_tracking_plot_spec(all_df):

    num_designed = all_df['InChIKey'].unique().shape[0]
    num_ordered = all_df[all_df["ORDERED"] == "TRUE"]['InChIKey'].unique().shape[0]
    num_made = all_df[all_df["MADE"] == "TRUE"]['InChIKey'].unique().shape[0]
    num_assayed = all_df[all_df["ASSAYED"] == "TRUE"]['InChIKey'].unique().shape[0]

    # tracking_df = pd.DataFrame(
    #     {
    #         "stage": ["Designed", "Ordered", "Made", "Assayed"],
    #         "color": ["Designed", "Ordered", "Made", "Assayed"],
    #         "num_mols": [num_designed, num_ordered, num_made, num_assayed],
    #     }
    # )

    tracking_df = pd.DataFrame(
        {
            "stage": ["Ordered", "Made", "Assayed"],
            "color": ["Ordered", "Made", "Assayed"],
            "num_mols": [num_ordered, num_made, num_assayed],
        }
    )

    print(tracking_df)

    # create Altair plot
    source = tracking_df
    tracking_plot = (
        alt.Chart(source)
        .mark_bar(opacity=0.7)
        .encode(
            y=alt.Y("stage:O", sort=["Ordered", "Made", "Assayed"], title="Stage"),
            x=alt.X(
                "num_mols:Q",
                stack=None,
                title=f"Number of Molecules (from {num_designed} unique designs)",
            ),
            color=alt.Color("color", legend=None),
        )
    )

    tracking_plot.save(
        str(dir_path / "scr" / "tracking_plot.json")
    )
    with open(str(dir_path / "scr" / "tracking_plot.json"), "r") as f:
        json_data = f.readlines()

    return json_data

