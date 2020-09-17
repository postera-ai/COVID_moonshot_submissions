# Create plot showing comparison of map of world and location
# general imports
import numpy as np
import pandas as pd

import altair as alt
from vega_datasets import data

from rdkit import Chem
from rdkit.Chem import AllChem

from pathlib import Path

dir_path = Path(__file__).parent.absolute()


def create_map_plot_spec(all_df):

    num_designed = all_df["InChIKey"].unique().shape[0]
    num_ordered = (
        all_df[all_df["ORDERED"] == "TRUE"]["InChIKey"].unique().shape[0]
    )
    num_made = all_df[all_df["MADE"] == "TRUE"]["InChIKey"].unique().shape[0]
    num_assayed = (
        all_df[all_df["ASSAYED"] == "TRUE"]["InChIKey"].unique().shape[0]
    )

    tracking_df = pd.DataFrame(
        {
            "stage": ["Ordered", "Made", "Assayed"],
            "color": ["Ordered", "Made", "Assayed"],
            "num_mols": [num_ordered, num_made, num_assayed],
        }
    )

    making_df = all_df.loc[
        (all_df["ORDERED"] == "TRUE") & (all_df["MADE"] == "FALSE")
    ]
    making_df = making_df.drop_duplicates(subset="InChIKey")
    num_compounds_dict = making_df["MAKER"].value_counts().to_dict()
    print(num_compounds_dict)

    world = data.world_110m.url
    world_topo = data.world_110m()

    loc_data = pd.DataFrame(
        {
            "location": [
                "Enamine",
                "Sai",
                "Mcule",
                "WuXi",
                "Molport",
                "Weizmann",
                "Diamond",
            ],
            "city": [
                "Kiev",
                "Hyderabad",
                "Budapest",
                "Shanghai",
                "Riga",
                "Tel Aviv",
                "Oxford",
            ],
            "country": [
                "Ukraine",
                "India",
                "Hungary",
                "China",
                "Latvia",
                "Israel",
                "UK",
            ],
            "latitude": [
                50.4501,
                17.3850,
                47.4979,
                31.2304,
                56.9496,
                32.0853,
                51.7520,
            ],
            "longitude": [
                30.5234,
                78.4867,
                19.0402,
                121.4737,
                24.1052,
                34.7818,
                1.2577,
            ],
            "num_compounds": [
                num_compounds_dict["enamine"],
                num_compounds_dict["sai"],
                num_compounds_dict["mcule"],
                num_compounds_dict["wuxi"],
                0,  # num_compounds_dict["molport"],
                num_made,
                num_made,
            ],
        }
    )

    routes = pd.DataFrame(
        {
            "origin": [
                "Enamine",
                "Enamine",
                "Sai",
                "Mcule",
                "WuXi",
                "Molport",
                "Weizmann",
                "Diamond",
                "Enamine",
                "Enamine",
                "Enamine",
                "Enamine",
            ],
            "destination": [
                "Weizmann",
                "Diamond",
                "Enamine",
                "Enamine",
                "Enamine",
                "Enamine",
                "Enamine",
                "Enamine",
                "Sai",
                "Mcule",
                "WuXi",
                "Molport",
            ],
        }
    )

    # interactive selection for origin
    # select nearest origin to mouse cursor
    origin = alt.selection_single(
        on="mouseover", nearest=True, fields=["origin"], empty="none"
    )

    # shared data reference for lookup transforms
    foreign = alt.LookupData(
        data=loc_data, key="location", fields=["latitude", "longitude"]
    )

    m = (
        alt.layer(
            # use the sphere of the Earth as the base layer
            alt.Chart({"sphere": True}).mark_geoshape(fill="#e6f3ff"),
            # add a graticule for geographic reference lines
            alt.Chart({"graticule": True}).mark_geoshape(
                stroke="#ffffff", strokeWidth=1
            ),
            # and then the countries of the world
            alt.Chart(alt.topo_feature(world, "countries")).mark_geoshape(
                fill="#ddd", stroke="#fff", strokeWidth=0.5
            ),
            alt.Chart(routes)
            .mark_rule(color="#000", opacity=0.35)
            .transform_filter(origin)  # filter to selected origin only
            .transform_lookup(lookup="origin", from_=foreign)  # origin lat/lon
            .transform_lookup(
                lookup="destination",
                from_=foreign,
                as_=["lat2", "lon2"],  # dest lat/lon
            )
            .encode(
                latitude="latitude:Q",
                longitude="longitude:Q",
                latitude2="lat2",
                longitude2="lon2",
            ),
            alt.Chart(routes)
            .mark_circle()
            .transform_aggregate(groupby=["origin"])
            .transform_lookup(
                lookup="origin",
                from_=alt.LookupData(
                    data=loc_data,
                    key="location",
                    fields=[
                        "location",
                        "city",
                        "country",
                        "latitude",
                        "longitude",
                        "num_compounds",
                    ],
                ),
            )
            .add_selection(origin)
            .encode(
                latitude="latitude:Q",
                longitude="longitude:Q",
                tooltip=["num_compounds:Q", "location:N"],
                size=alt.Size(
                    "num_compounds:Q",
                    scale=alt.Scale(range=[100, 500]),
                    legend=None,
                ),
            ),
        )
        .project(type="naturalEarth1", scale=500, translate=[45, 800])
        .configure_view(stroke=None)
        .properties(width=1200, height=800)
    )

    m.save(str(dir_path / "scr" / "world_map.html"))
    m.save(str(dir_path / "scr" / "world_map.json"))
    with open(str(dir_path / "scr" / "world_map.json"), "r") as f:
        json_data = f.readlines()

    return json_data
