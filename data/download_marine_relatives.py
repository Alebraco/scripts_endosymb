#!/usr/bin/env python3
"""Download marine free-living relative genomes by TAXON NAME (NCBI datasets CLI).

These genomes become Group=relatives_only rows for the Model B training set.

Output structure:
    <outdir>/genomes/<Species_dir>/<ACCESSION>.fna

Workflow:
    - Verification mode: Print the assemblies that would be fetched for verification.
    - Download mode: fetch each selected genome, unzip, and place the .fna.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

# (taxon, max genomes), group=relatives_only.
TARGETS = [
    ("Pelagibacter ubique HTCC1062", 1),
    ("Pelagibacter sp. HTCC7211", 1),
    ("Pelagibacter sp. IMCC9063", 1),
    ("Prochlorococcus marinus", 1),
    ("Polaribacter", 3),
    ("Dokdonia", 3),
    ("Croceibacter", 2),
    ("Formosa", 3),
]

# Prefer the most finished assemblies first.
LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
FASTA_EXTS = (".fna", ".fa", ".fasta")


def species_dir(org, strain):
    """Sanitize organism (+strain if not already present) into a directory name."""
    label = org or "Unknown"
    if strain and strain.lower() not in label.lower():
        label = f"{label} {strain}"
    label = re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_")
    return label or "Unknown"


def query_and_select(taxon, max_n):
    """Query RefSeq for taxon and return (all_reports, top_n_chosen)."""
    try:
        proc = subprocess.run(
            ["datasets", "summary", "genome", "taxon", taxon,
             "--assembly-source", "RefSeq", "--as-json-lines"],
            capture_output=True, text=True, check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"  Summary failed for '{taxon}': {exc}", file=sys.stderr)
        return [], []
    reports = []
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            reports.append(json.loads(line))
        except json.JSONDecodeError:
            continue

    def rank(r):
        lvl = (r.get("assembly_info") or {}).get("assembly_level", "")
        size = (r.get("assembly_stats") or {}).get("total_sequence_length")
        return (LEVEL_RANK.get(lvl, 9), -(size if isinstance(size, int) else 0))

    seen, chosen = set(), []
    for r in sorted(reports, key=rank):
        acc = r.get("accession")
        if not acc or acc in seen:
            continue
        seen.add(acc)
        chosen.append(r)
        if len(chosen) >= max_n:
            break
    return reports, chosen


def download_one(accession, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        zip_path = os.path.join(tmp, f"{accession}.zip")
        subprocess.run(
            ["datasets", "download", "genome", "accession", accession,
             "--include", "genome", "--filename", zip_path],
            check=True,
        )
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmp)
        fna_files = []
        for root, _, files in os.walk(os.path.join(tmp, "ncbi_dataset", "data")):
            for f in files:
                if f.endswith(FASTA_EXTS):
                    fna_files.append(os.path.join(root, f))
        if not fna_files:
            print(f"  No FASTA found in download for {accession}", file=sys.stderr)
            return False
        dest = os.path.join(dest_dir, f"{accession}.fna")
        with open(dest, "wb") as out:
            for fp in sorted(fna_files):
                with open(fp, "rb") as src:
                    shutil.copyfileobj(src, out)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default="marine_free_livers",
                        help="Base dir, genomes go to <outdir>/genomes/<Species>/.")
    parser.add_argument("--download", action="store_true",
                        help="Download mode. Without this flag, only resolve and "
                             "print selections for verification.")
    args = parser.parse_args()

    selections = []
    for taxon, max_n in TARGETS:
        reports, chosen = query_and_select(taxon, max_n)
        print(f"\n# {taxon}  (found {len(reports)} RefSeq assemblies, keeping {len(chosen)})")
        if not chosen:
            print("  Nothing selected, verify taxon name")
            continue
        for r in chosen:
            org = (r.get("organism") or {}).get("organism_name", "?")
            strain = ((r.get("organism") or {}).get("infraspecific_names") or {}).get("strain")
            stats = r.get("assembly_stats") or {}
            size = stats.get("total_sequence_length", "?")
            gc = stats.get("gc_percent", "?")
            level = (r.get("assembly_info") or {}).get("assembly_level", "?")
            sdir = species_dir(org, strain)
            acc = r["accession"]
            size_str = f"{size:,}" if isinstance(size, int) else str(size)
            print(f"  {acc} {size_str} bp  GC={gc}  {level} "
                  f"{org}{(' / ' + strain) if strain else ''}  -> {sdir}")
            selections.append((sdir, acc))

    print(f"\nTotal selected: {len(selections)} genomes")

    if not args.download:
        print("\nVerify the organisms, then re-run with --download.")
        return

    if not selections:
        sys.exit("Nothing to download.")

    print("\nDownloading...")
    genomes_root = os.path.join(args.outdir, "genomes")
    n = 0
    for sdir, acc in selections:
        dest_dir = os.path.join(genomes_root, sdir)
        print(f"  {acc} -> {dest_dir}")
        if download_one(acc, dest_dir):
            n += 1
    print(f"\nDownloaded {n}/{len(selections)} genomes into {genomes_root}/")


if __name__ == "__main__":
    main()
