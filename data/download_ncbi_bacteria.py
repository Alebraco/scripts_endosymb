#!/usr/bin/env python3

import os
import argparse
import subprocess
import pandas as pd


ASSEMBLY_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"   # RefSeq bacteria assembly summary URL

GENBANK_SUMMARY_URL  = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"  # GenBank bacteria assembly summary URL (used for MAGs)

NON_STANDARD_CODE_GENERA = {
    "Mycoplasma", "Spiroplasma", "Ureaplasma", "Phytoplasma",
    "Acholeplasma", "Mesoplasma", "Entomoplasma",  
}


def download_file(url, dest_path):
    if os.path.exists(dest_path):
        return
    try:
        subprocess.run(['curl', '-o', dest_path, url])
    except Exception as e:
        raise RuntimeError(f"Failed to download {url}: {e}")
    

def load_assembly_summary(local_path):
    if not os.path.exists(local_path):
        print("Downloading assembly summary from NCBI FTP.")
        download_file(ASSEMBLY_SUMMARY_URL, local_path)
        print(f"  Saved to {local_path}")
    else:
        print(f"Using cached assembly summary: {local_path}")

    df = pd.read_csv(local_path, sep="\t", skiprows=1, low_memory=False)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    return df


def filter_assemblies(df, max_per_species):
    '''Apply filters to the assembly summary DataFrame:
       1. Keep only rows with valid FTP paths
       2. Drop assemblies marked as excluded from RefSeq
       3. Drop lineages known to use non-standard genetic codes
       4. Within each species_taxid, keep at most max_per_species assemblies
    '''

    n_start = len(df)

    # 1. Valid FTP path
    df = df[df["ftp_path"].notna()].copy()
    print(f"  After removing missing FTP paths:  {len(df)} / {n_start}")

    # 2. Drop assemblies excluded from RefSeq
    if "excluded_from_refseq" in df.columns:
        df = df[df["excluded_from_refseq"].isna() | (df["excluded_from_refseq"].str.strip().isin(["", "na"]))].copy()
        print(f"  After removing excluded assemblies:{len(df)}")

    # 3. Drop non-standard-code lineages
    def is_nonstandard(name):
        genus = str(name).split()[0]
        return genus in NON_STANDARD_CODE_GENERA  

    mask = df["organism_name"].apply(is_nonstandard)
    df = df[~mask].copy()
    print(f"  After removing non-std code genera: {len(df)} (removed {mask.sum()})")

    df = (
        df.groupby("species_taxid", sort=False)
          .head(max_per_species)
          .reset_index(drop=True)
    )
    print(f"  After {max_per_species} per species: {len(df)}")
    return df


def safe_dirname(name):
    return name.strip().replace(" ", "_").replace("/", "_").replace("'", "")

def download_genomes(df, outdir):
    """
    Download genomic FASTA files from NCBI FTP.
    Output: outdir/genomes/{species_name}/{accession}.fna
    """

    genomes_dir = os.path.join(outdir, "genomes")
    os.makedirs(genomes_dir, exist_ok=True)

    n_total  = len(df)   # total genomes to download
    n_done   = 0         # count of newly downloaded files
    n_skip   = 0         # count of files already present
    n_failed = 0         # count of download failures

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

        try:
            print(f"[{n_done + n_skip + 1}/{n_total}] {accession}  {organism}")
            download_file(fna_gz_url, dest_fna + ".gz")
            subprocess.run(['gzip', '-d', dest_fna])
            n_done += 1

        except Exception as e:
            print(f"  FAILED: {accession} — {e}")
            n_failed += 1

    print(f"\nDone. {n_done} downloaded, {n_skip} already present, {n_failed} failed.")


def load_genbank_mags(local_path):
    """Download and parse GenBank assembly_summary.txt, filtering to MAGs only."""
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
        print("'assembly_type' column not found, ignoring MAG filtering.")
        return pd.DataFrame()

    print(f"  MAGs found in GenBank: {len(df):,}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Download bacterial genomes from NCBI RefSeq."
    )
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
        help="Local path for RefSeq assembly_summary.txt (downloaded if absent)"
    )
    parser.add_argument(
        "--include-mags", action="store_true",
        help="Also download MAGs from GenBank (assembly_type = metagenome-assembled)"
    )
    args = parser.parse_args()

    # RefSeq genomes
    df = load_assembly_summary(args.assembly_summary)
    print(f"Total RefSeq assemblies in summary: {len(df)}")

    df_filtered = filter_assemblies(df, args.max_per_species)

    # GenBank MAGs
    if args.include_mags:
        if "assembly_summary" in args.assembly_summary:
            genbank_file = args.assembly_summary.replace("assembly_summary", "genbank_assembly_summary")
        else:
            genbank_file = "genbank_assembly_summary.txt"
        df_mags = load_genbank_mags(genbank_file)

        if not df_mags.empty:
            df_mags_filtered = filter_assemblies(df_mags, args.max_per_species)
            df_filtered = pd.concat([df_filtered, df_mags_filtered], ignore_index=True)
            print(f"\nCombined total (RefSeq + MAGs): {len(df_filtered)}")

    print(f"\nDownloading {len(df_filtered):,} genomes to {args.outdir}/genomes/ …\n")
    download_genomes(df_filtered, outdir=args.outdir)   # download all selected assemblies

if __name__ == "__main__":
    main()
