import json
import argparse
import subprocess
from pathlib import Path
from xml.etree import ElementTree

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--input-file",
    default="benchmark.xml",
    help="Path to the XML benchmark results",
)
parser.add_argument(
    "-o",
    "--output-file",
    default="benchmark.json",
    help="Path to collect the results into",
)
parser.add_argument("--identifier", help="A unique identifier for this set of results")
parser.add_argument(
    "--overwrite",
    action="store_true",
    help="Overwrite the previous results for this identifier if they exist",
)
args = parser.parse_args()

if args.identifier:
    identifier = args.identifier
else:
    identifier = (
        subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()
    )
print(f"Collecting results for: {identifier}")

out_path = Path(args.output_file)
if out_path.exists():
    with open(out_path, "r") as f:
        previous = json.load(f)
else:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    previous = []
for result in previous:
    if result["identifier"] == identifier:
        if args.overwrite:
            previous.remove(result)
        else:
            raise RuntimeError(f"Results already collected for {identifier}")

print(f"Collecting from: {args.input_file}")
results = []
tree = ElementTree.parse(args.input_file)
root = tree.getroot()
for test_case in root:
    if test_case.tag != "TestCase":
        continue
    info = dict(test_case.attrib)
    info["results"] = []
    for test_result in test_case:
        if test_result.tag != "BenchmarkResults":
            continue
        data = dict(test_result.attrib)
        for test in test_result:
            data[test.tag] = {k: float(v) for k, v in test.attrib.items()}
        info["results"].append(data)
    results.append(info)

print(f"Writing to: {out_path}")
previous.append({"identifier": identifier, "results": results})
with open(out_path, "w") as f:
    json.dump(previous, f, indent=4)
