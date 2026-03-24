#!/usr/bin/env python3
"""
Download ~N genomic FASTA files per bacterial species from the NCBI RefSeq FTP.

Usage:
    python download_ncbi_bacteria.py [--outdir ncbi_bacteria] [--max-per-species 10]

Output layout:
    ncbi_bacteria/genomes/{species_name}/{accession}.fna

Filters applied before downloading:
  - Drop assemblies with no FTP path
  - Drop assemblies marked as excluded from RefSeq
  - Drop lineages known to use non-standard genetic codes (table != 11)
  - Within each species (species_taxid), keep at most --max-per-species genomes,
    preferring higher assembly levels (Complete Genome > Chromosome > Scaffold > Contig)
"""

import os
import gzip
import shutil
import argparse
import subprocess
import pandas as pd


ASSEMBLY_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
GENBANK_SUMMARY_URL  = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"

ASSEMBLY_LEVEL_PRIORITY = {
    "Complete Genome": 0,
    "Chromosome": 1,
    "Scaffold": 2,
    "Contig": 3,
}

# Genera known to use non-standard translation tables.
# Annotating with table 11 corrupts codon-bias features, so we exclude them.
# Table 4: Mycoplasma, Spiroplasma, Ureaplasma, Phytoplasma (UGA = Trp)
# Table 25: some CPR / Gracilibacteria (UGA = Gly)
NON_STANDARD_CODE_GENERA = {
    "Mycoplasma", "Spiroplasma", "Ureaplasma", "Phytoplasma",
    "Acholeplasma", "Mesoplasma", "Entomoplasma",  # other Mollicutes with table 4
}


def download_file(url, dest_path):
    """Download a URL to dest_path via wget; skips if already present."""
    if os.path.exists(dest_path):
        return
    result = subprocess.run(['wget', '-q', '-O', dest_path, url])
    if result.returncode != 0:
        raise RuntimeError(f"wget failed for {url} (exit code {result.returncode})")


def load_assembly_summary(local_path):
    """Download (if needed) and parse assembly_summary.txt into a DataFrame."""
    if not os.path.exists(local_path):
        print("Downloading assembly summary from NCBI FTP …")
        download_file(ASSEMBLY_SUMMARY_URL, local_path)
        print(f"  Saved to {local_path}")
    else:
        print(f"Using cached assembly summary: {local_path}")

    df = pd.read_csv(local_path, sep="\t", skiprows=1, low_memory=False)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    return df


def filter_assemblies(df, max_per_species):
    """
    Apply all pre-download filters and subsample to max_per_species per taxid.

    Steps
    -----
    1. Drop rows with no FTP path.
    2. Drop assemblies the NCBI has flagged to exclude from RefSeq.
    3. Drop genera that use non-standard genetic codes.
    4. Sort by assembly level quality, then subsample to max_per_species per species_taxid.
    """
    n_start = len(df)

    # 1. Must have a valid FTP path
    df = df[df["ftp_path"].notna() & (df["ftp_path"] != "na")].copy()
    print(f"  After removing missing FTP paths:  {len(df):>7,} / {n_start:,}")

    # 2. Drop assemblies excluded from RefSeq (contamination, low quality, etc.)
    if "excluded_from_refseq" in df.columns:
        df = df[df["excluded_from_refseq"].isna() | (df["excluded_from_refseq"].str.strip() == "")].copy()
        print(f"  After removing excluded assemblies:{len(df):>7,}")

    # 3. Drop non-standard-code lineages (first word of organism_name = genus)
    def is_nonstandard(name):
        genus = str(name).split()[0]
        return genus in NON_STANDARD_CODE_GENERA

    mask = df["organism_name"].apply(is_nonstandard)
    df = df[~mask].copy()
    print(f"  After removing non-std code genera: {len(df):>7,} (removed {mask.sum():,})")

    # 4. Sort by assembly quality and subsample per species
    df["_level_rank"] = df["assembly_level"].map(ASSEMBLY_LEVEL_PRIORITY).fillna(99)
    df = df.sort_values("_level_rank")

    df = (
        df.groupby("species_taxid", sort=False)
          .head(max_per_species)
          .drop(columns=["_level_rank"])
          .reset_index(drop=True)
    )
    print(f"  After ≤{max_per_species} per species:         {len(df):>7,}")
    return df


def safe_dirname(name):
    """Convert an organism name to a filesystem-safe directory name."""
    return name.strip().replace(" ", "_").replace("/", "_").replace("'", "")


def download_genomes(df, outdir):
    """
    Download genomic FASTA files from NCBI FTP.
    Output: outdir/genomes/{species_name}/{accession}.fna
    """
    genomes_dir = os.path.join(outdir, "genomes")
    os.makedirs(genomes_dir, exist_ok=True)

    n_total  = len(df)
    n_done   = 0
    n_skip   = 0
    n_failed = 0

    for _, row in df.iterrows():
        accession = str(row["assembly_accession"])
        organism  = str(row["organism_name"])
        ftp_base  = str(row["ftp_path"]).rstrip("/")

        fname_base  = ftp_base.split("/")[-1]
        fna_gz_url  = f"{ftp_base}/{fname_base}_genomic.fna.gz"

        species_dir = os.path.join(genomes_dir, safe_dirname(organism))
        os.makedirs(species_dir, exist_ok=True)

        dest_fna = os.path.join(species_dir, f"{accession}.fna")

        if os.path.exists(dest_fna):
            n_skip += 1
            continue

        gz_tmp = dest_fna + ".gz"
        try:
            print(f"[{n_done + n_skip + 1}/{n_total}] {accession}  {organism}")
            download_file(fna_gz_url, gz_tmp)

            with gzip.open(gz_tmp, "rb") as f_in, open(dest_fna, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gz_tmp)
            n_done += 1

        except Exception as e:
            print(f"  FAILED: {accession} — {e}")
            if os.path.exists(gz_tmp):
                os.remove(gz_tmp)
            n_failed += 1

    print(f"\nDone. {n_done} downloaded, {n_skip} already present, {n_failed} failed.")


def load_genbank_mags(local_path):
    """Download (if needed) and parse GenBank assembly_summary.txt, filtering to MAGs only."""
    if not os.path.exists(local_path):
        print("Downloading GenBank assembly summary …")
        download_file(GENBANK_SUMMARY_URL, local_path)
        print(f"  Saved to {local_path}")
    else:
        print(f"Using cached GenBank assembly summary: {local_path}")

    df = pd.read_csv(local_path, sep="\t", skiprows=1, low_memory=False)
    df.columns = [c.lstrip("#").strip() for c in df.columns]

    # Keep only metagenome-assembled genomes
    if "assembly_type" in df.columns:
        df = df[df["assembly_type"].str.contains("metagenome", case=False, na=False)].copy()
    else:
        print("  Warning: 'assembly_type' column not found; cannot filter for MAGs.")
        return pd.DataFrame()

    print(f"  MAGs found in GenBank: {len(df):,}")
    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir", default="ncbi_bacteria",
        help="Root output directory (default: ncbi_bacteria)"
    )
    parser.add_argument(
        "--max-per-species", type=int, default=10,
        help="Maximum genomes per species taxid (default: 10)"
    )
    parser.add_argument(
        "--assembly-summary", default="assembly_summary.txt",
        help="Local cache path for RefSeq assembly_summary.txt (downloaded if absent)"
    )
    parser.add_argument(
        "--include-mags", action="store_true",
        help="Also download MAGs from GenBank (assembly_type = metagenome-assembled)"
    )
    args = parser.parse_args()

    # RefSeq isolates
    df = load_assembly_summary(args.assembly_summary)
    print(f"Total RefSeq assemblies in summary: {len(df):,}")

    print("\nApplying filters …")
    df_filtered = filter_assemblies(df, args.max_per_species)

    # GenBank MAGs (optional)
    if args.include_mags:
        genbank_cache = args.assembly_summary.replace("assembly_summary", "genbank_assembly_summary")
        if genbank_cache == args.assembly_summary:
            genbank_cache = "genbank_assembly_summary.txt"
        df_mags = load_genbank_mags(genbank_cache)
        if not df_mags.empty:
            print("\nApplying filters to MAGs …")
            df_mags_filtered = filter_assemblies(df_mags, args.max_per_species)
            df_filtered = pd.concat([df_filtered, df_mags_filtered], ignore_index=True)
            print(f"\nCombined total (RefSeq + MAGs): {len(df_filtered):,}")

    print(f"\nDownloading {len(df_filtered):,} genomes to {args.outdir}/genomes/ …\n")
    download_genomes(df_filtered, outdir=args.outdir)


if __name__ == "__main__":
    main()
