import argparse
import json
import re
from collections import defaultdict

import altair as alt
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--input-file",
    default="benchmark.json",
    help="Path to the collected JSON results",
)
parser.add_argument(
    "-o",
    "--output-file",
    default="benchmark.html",
    help="Path to save the output HTML file",
)
args = parser.parse_args()

with open(args.input_file, "r") as f:
    data = json.load(f)

df = defaultdict(list)
for entry in data:
    identifier = entry["identifier"]
    for group in entry["results"]:
        group_name = group["name"].split(" - ")[0]
        for row in group["results"]:
            df["url"].append(f"https://github.com/dfm/commit/{identifier}")
            df["identifier"].append(identifier)
            df["hash"].append(identifier[:7])
            df["group"].append(group_name)
            e, n = re.match(r"e=(.+); n=(.+)", row["name"].strip(": ")).groups()
            df["eccentricity"].append(float(e))
            for k in ["mean", "standardDeviation"]:
                df[k].append(row[k]["value"] / int(n))

df = pd.DataFrame(df)
print(df.head())

plot1 = (
    alt.Chart(
        df[(df["identifier"] == identifier) & (df["group"] == "baselined")],
    )
    .mark_rule(color="black", size=2)
    .encode(y="mean:Q")
)

m = df["identifier"] == identifier
m &= ~df["group"].str.startswith("baseline")
for mp, shape in [
    (df["group"].str.endswith("d"), "circle"),
    (df["group"].str.endswith("dv"), "square"),
    (df["group"].str.startswith("ref:"), "diamond"),
]:
    plot1 += (
        alt.Chart(df[m & mp], width=400)
        .mark_line(point=alt.OverlayMarkDef(shape=shape, size=50))
        .encode(
            x=alt.X("eccentricity:Q"),
            y=alt.Y(
                "mean:Q",
                scale=alt.Scale(type="log", nice=False, padding=1),
                title="Time (ns)",
            ),
            color=alt.Color("group:N", scale=alt.Scale(scheme="tableau20")),
        )
    )

plot2 = (
    alt.Chart(
        df[(df["identifier"] == identifier) & (df["group"] == "baselinef")],
    )
    .mark_rule(color="black", size=2)
    .encode(y="mean:Q")
)

m = df["identifier"] == identifier
m &= ~df["group"].str.startswith("baseline")
for mp, shape in [
    (df["group"].str.endswith("f"), "circle"),
    (df["group"].str.endswith("fv"), "square"),
]:
    plot2 += (
        alt.Chart(df[m & mp], width=400)
        .mark_line(point=alt.OverlayMarkDef(shape=shape, size=50))
        .encode(
            x=alt.X("eccentricity:Q"),
            y=alt.Y(
                "mean:Q",
                scale=alt.Scale(type="log", nice=False, padding=1),
                title="Time (ns)",
            ),
            color=alt.Color("group:N", scale=alt.Scale(scheme="tableau20")),
        )
    )

plot = plot1
# plot = plot1 | plot2
plot.save(args.output_file)


# print(df[df["group"] == "baselined"])

# m = df["eccentricity"] == 0.2
# plot1 = (
#     alt.Chart(df[m & df["group"].str.endswith("d")], width=40)
#     .mark_line(point=alt.OverlayMarkDef(shape="circle"))
#     .encode(
#         x=alt.X("hash:N"),
#         y=alt.Y(
#             "mean:Q",
#             scale=alt.Scale(type="log", nice=False, padding=1),
#             title="Time (ns)",
#         ),
#         color=alt.Color("group:N", scale=alt.Scale(scheme="tableau20")),
#         href="url",
#         tooltip=["identifier", "url"],
#     )
# )
# plot2 = (
#     alt.Chart(df[m & df["group"].str.endswith("dv")], width=40)
#     .mark_line(point=alt.OverlayMarkDef(shape="square"))
#     .encode(
#         x=alt.X("hash:N"),
#         y=alt.Y(
#             "mean:Q",
#             scale=alt.Scale(type="log", nice=False, padding=1),
#             title="Time (ns)",
#         ),
#         color=alt.Color("group:N", scale=alt.Scale(scheme="tableau20")),
#         href="url",
#         tooltip=["identifier", "url"],
#     )
# )
# plot = plot1 + plot2

# plot.save("plot.html")
