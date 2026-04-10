#!/usr/bin/env python3
"""
Download bacterial genomes and Bakta annotations from AllTheBacteria (ATB).

Usage:
  python download_atb_bacteria.py \
      --file-list file_list.all.202505.tsv \
      --bakta-status atb.bakta.r0.2.status.tsv \
      --batch-urls batch_urls.tsv \
      --outdir atb_bacteria \
      --max-per-species 10

Input files:
  file_list.all.202505.tsv   — ATB assembly manifest (one row per genome,
                               includes sample ID, sylph_species, tar path)
  atb.bakta.r0.2.status.tsv — Bakta annotation status per sample; only
                               samples with status=PASS are kept
  batch_urls.tsv             — Two-column TSV: tar_xz filename, download URL
                               for each Bakta batch archive

Output layout (inside --outdir):
  manifest.tsv               — Final selected samples with all relevant paths
  asm_batch_index.tsv        — Index for SLURM assembly array: tar_xz, url
  bkt_batch_index.tsv        — Index for SLURM bakta array: tar_xz, url
  logs/                      — Created for SLURM log output
"""

import os
import re
import argparse
import pandas as pd

# Excluded because codon usage analyses assume the standard code.
NON_STANDARD_CODE_GENERA = {
    "Mycoplasma", "Spiroplasma", "Ureaplasma", "Phytoplasma",
    "Acholeplasma", "Mesoplasma", "Entomoplasma",
}


def safe_dirname(name):
    """Convert species name to a safe directory name."""
    return name.strip().replace(" ", "_").replace("/", "_").replace("'", "")


def strip_gtdb_suffix(genus):
    """Remove GTDB-style suffixes (e.g. _D, _A) appended to genus names."""
    return re.sub(r"_[A-Z]$", "", genus)


def load_file_list(path):
    """Load the ATB assembly manifest TSV."""
    print(f"Loading assembly file list: {path}")
    df = pd.read_csv(path, sep="\t", low_memory=False)
    print(f"  Total entries: {len(df)}")
    return df


def load_bakta_status(path):
    """Load the Bakta annotation status TSV."""
    print(f"Loading Bakta status: {path}")
    df = pd.read_csv(path, sep="\t", low_memory=False)
    print(f"  Total entries: {len(df)}")
    return df


def load_batch_urls(path):
    """Load the Bakta batch download URL table."""
    print(f"Loading batch URLs: {path}")
    df = pd.read_csv(path, sep="\t", header=None, names=["tar_xz", "url"])
    print(f"  Total batches: {len(df)}")
    return df


def filter_and_select(file_list, bakta_status, batch_urls, max_per_species):
    """Apply quality filters and select up to max_per_species genomes per species.

      1. Drop unknown species.
      2. Drop non-standard genetic code genera.
      3. Inner-join with Bakta status.
      4. Attach Bakta batch download URLs.
      5. Take the first max_per_species rows per species.
    """

    df = file_list.copy()
    n_start = len(df)

    # 1. Remove samples where sylph could not assign a species
    df = df[df["sylph_species"] != "unknown"].copy()
    print(f"  After removing unknown species: {len(df)} / {n_start}")

    # 2. Remove non-standard genetic code genera.
    def is_nonstandard(species_name):
        genus = str(species_name).split()[0]
        return strip_gtdb_suffix(genus) in NON_STANDARD_CODE_GENERA

    mask = df["sylph_species"].apply(is_nonstandard)
    df = df[~mask].copy()
    print(f"After removing non-std code genera: {len(df)} (removed {mask.sum()})")

    # 3. Keep only samples with a successful Bakta annotation (status == PASS).
    # The join also attaches the bakta tar_xz filename for each sample.
    bakta_pass = bakta_status[bakta_status["status"] == "PASS"][["sample", "tar_xz"]].copy()
    bakta_pass = bakta_pass.rename(columns={"tar_xz": "bakta_tar_xz"})
    df = df.merge(bakta_pass, on="sample", how="inner")
    print(f"After joining with PASS bakta: {len(df)}")

    # 4. Map each sample's bakta tar_xz to its download URL.
    # Samples whose batch tar has no URL in batch_urls.tsv are dropped.
    url_map = dict(zip(batch_urls["tar_xz"], batch_urls["url"]))
    df["bakta_url"] = df["bakta_tar_xz"].map(url_map)
    n_no_url = df["bakta_url"].isna().sum()
    if n_no_url > 0:
        print(f"WARNING: {n_no_url} samples have no bakta batch URL — dropping them")
        df = df.dropna(subset=["bakta_url"])

    # 5. Take N rows per species.
    df = (
        df.groupby("sylph_species", sort=False)
          .head(max_per_species)
          .reset_index(drop=True)
    )
    n_species = df["sylph_species"].nunique()
    print(f"  After {max_per_species} per species: {len(df)} samples across {n_species} species")

    return df


def write_manifest(df, outdir):
    """Write the full selected-sample manifest to outdir/manifest.tsv.

    Columns: sample, sylph_species, filename_in_tar_xz, tar_xz (assembly), tar_xz_url (assembly URL), bakta_tar_xz (bakta), bakta_url.

    """

    manifest_path = os.path.join(outdir, "manifest.tsv")
    cols = ["sample", "sylph_species", "filename_in_tar_xz", "tar_xz", "tar_xz_url",
            "bakta_tar_xz", "bakta_url"]
    df[cols].to_csv(manifest_path, sep="\t", index=False)
    print(f"\nManifest written: {manifest_path} ({len(df)} samples)")
    return manifest_path


def write_batch_index(df, outdir):
    """Write SLURM array index files: one line per tar.xz, with format: tar_xz   url

    SLURM worker scripts filter manifest.tsv by tar_xz to find their samples.
    """
    asm_index_path = os.path.join(outdir, "asm_batch_index.tsv")
    with open(asm_index_path, "w") as f:
        for row in df[["tar_xz", "tar_xz_url"]].drop_duplicates().itertuples(index=False):
            f.write(f"{row.tar_xz}\t{row.tar_xz_url}\n")

    bkt_index_path = os.path.join(outdir, "bkt_batch_index.tsv")
    with open(bkt_index_path, "w") as f:
        for row in df[["bakta_tar_xz", "bakta_url"]].drop_duplicates().itertuples(index=False):
            f.write(f"{row.bakta_tar_xz}\t{row.bakta_url}\n")

    n_asm = df["tar_xz"].nunique()
    n_bkt = df["bakta_tar_xz"].nunique()
    print(f"Assembly batches: {n_asm}\nBakta batches: {n_bkt}")

    if n_asm == n_bkt:
        print("Assembly and Bakta batch counts match.")


def main():
    parser = argparse.ArgumentParser(
        description="Select and prepare ATB genome + annotation downloads."
    )
    parser.add_argument(
        "--outdir", default="atb_bacteria",
        help="Root output directory (default: atb_bacteria)"
    )
    parser.add_argument(
        "--max-per-species", type=int, default=10,
        help="Maximum genomes per species (default: 10)"
    )
    parser.add_argument(
        "--file-list", required=True,
        help="Path to file_list.all.202505.tsv (assembly list with species)"
    )
    parser.add_argument(
        "--bakta-status", required=True,
        help="Path to atb.bakta.r0.2.status.tsv (bakta annotation status)"
    )
    parser.add_argument(
        "--batch-urls", required=True,
        help="Path to batch_urls.tsv (bakta batch download URLs)"
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "logs"), exist_ok=True)

    file_list = load_file_list(args.file_list)
    bakta_status = load_bakta_status(args.bakta_status)
    batch_urls = load_batch_urls(args.batch_urls)

    print("\nFiltering and selecting samples:")
    selected = filter_and_select(file_list, bakta_status, batch_urls, args.max_per_species)

    write_manifest(selected, args.outdir)
    write_batch_index(selected, args.outdir)

if __name__ == "__main__":
    main()
