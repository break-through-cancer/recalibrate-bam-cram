"""
preprocess.py -- Cirro dataset preprocessing hook for the BQSR pipeline.
Injects params.bqsr_runs directly, same pattern as the MuTect1 and
CRAM-to-BAM pipelines' preprocess.py hooks.

Only CRAMs from the MarkDuplicates stage are considered -- this is the
duplicate-marked, NOT-yet-recalibrated stage BQSR expects as input.
Filters primarily on the dataset's own `bamType` column (seen as
"markduplicates" / "mapped" / "recalibrated" in this project's real
ds.files output) since that's more robust than pattern-matching the path;
falls back to a path-based regex if that column isn't present.
"""
import json
import re
from cirro.helpers.preprocess_dataset import PreprocessDataset

MD_FOLDER_PATTERN = r"/markduplicates/"


def extract_crams(ds):
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    if "bamType" in df.columns:
        df = df[df["bamType"] == "markduplicates"]
    else:
        df = df[df["file"].str.contains(MD_FOLDER_PATTERN, regex=True)]

    df = df[df["file"].str.endswith(".cram") | df["file"].str.endswith(".cram.crai")]

    cram_map = {}
    for sample, group in df.groupby("sample"):
        cram, crai = "", ""
        for f in group["file"]:
            if f.endswith(".cram") and not f.endswith(".cram.crai"):
                if cram:
                    raise ValueError(
                        f"Multiple markduplicates CRAMs found for sample {sample!r}: "
                        f"{cram!r} and {f!r} -- expected exactly one."
                    )
                cram = f
            elif f.endswith(".cram.crai"):
                crai = f
        if cram:
            cram_map[str(sample)] = {"cram": cram, "crai": crai}

    if not cram_map:
        raise ValueError(
            "No markduplicates-stage CRAMs found in this dataset "
            "(checked bamType == 'markduplicates', or path matching "
            f"{MD_FOLDER_PATTERN!r} as a fallback)"
        )
    return cram_map


def main():
    ds = PreprocessDataset.from_running()

    print("=== ds.files preview ===")
    print(ds.files.head(20).to_string(index=False))

    cram_map = extract_crams(ds)

    bqsr_runs = []
    for sample, files in cram_map.items():
        if not files["crai"]:
            raise ValueError(f"Sample {sample!r} has a CRAM but no matching .crai index")
        bqsr_runs.append({
            "sample_id": sample,
            "cram": files["cram"],
            "crai": files["crai"],
        })

    ds.add_param("bqsr_runs", bqsr_runs)

    print("\nFinal parameters:")
    print(json.dumps(ds.params, indent=2, default=str))


if __name__ == "__main__":
    main()