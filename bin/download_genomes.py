#!/usr/bin/env python3
"""
Created: Connor S. Murray
Date: 2.14.2025

Download NCBI assembly files (proteome FASTA, genomic/cds FASTA, and GTF) 
given a list of assembly entries, and unzip the downloaded files.

The script expects each line of the input file to contain three columns:
  1. The full assembly identifier (e.g., GCF_000001405.39_GRCh38.p13)
  2. The renamed basename to use for the downloaded files
  3. A group status value (ignored in this script)

It then builds the FTP path using the standard NCBI directory structure:
  https://ftp.ncbi.nih.gov/genomes/all/GC[A/F]/###/###/###/AssemblyFolder

And downloads the following files:
  - Proteome FASTA: <basename>.protein.faa.gz
  - Genomic FASTA:  <basename>.genomic.fna.gz
  - CDS FASTA:      <basename>.cds_from_genomic.fna.gz
  - GTF:            <basename>.genomic.gtf.gz

After downloading, each file is automatically unzipped.

Usage:
    python download_ncbi_assemblies.py --assembly_list assemblies.txt --out_dir ncbi_downloads
"""

import argparse
import os
import sys
import urllib.request
import gzip
import shutil

def construct_ftp_path(assembly):
    # The assembly is assumed to start with GCF_ or GCA_
    # and be in the form: GCF_XXXXXXXXX.<version>_<assembly_name>
    # Extract the prefix and the nine-digit number.
    try:
        # Split at underscores.
        parts = assembly.split("_")
        if len(parts) < 2:
            raise ValueError("Assembly string is not in expected format.")
        prefix = parts[0]  # e.g. GCF or GCA
        num_and_version = parts[1]  # e.g. 000001405.39
        # Get the nine-digit number (ignore version part)
        nine_digit = num_and_version.split('.')[0]  # e.g. 000001405
        if len(nine_digit) != 9:
            raise ValueError("Expected nine digits in assembly number.")
        # Break into 3-digit folders.
        folder1 = nine_digit[0:3]
        folder2 = nine_digit[3:6]
        folder3 = nine_digit[6:9]
        ftp_path = f"https://ftp.ncbi.nih.gov/genomes/all/{prefix}/{folder1}/{folder2}/{folder3}/{assembly}"
        return ftp_path
    except Exception as e:
        print(f"[ERROR] Failed to construct FTP path for {assembly}: {e}", file=sys.stderr)
        return None

def download_file(url, out_path):
    if os.path.exists(out_path):
        print(f"File already exists: {out_path} (skipping)")
        return
    print(f"Downloading:\n  {url}\n  -> {out_path}")
    try:
        urllib.request.urlretrieve(url, out_path)
    except Exception as e:
        print(f"[ERROR] Could not download {url}: {e}", file=sys.stderr)

def unzip_file(gz_path):
    if not os.path.exists(gz_path):
        print(f"File does not exist: {gz_path} (skipping unzip)")
        return
    if not gz_path.endswith('.gz'):
        print(f"File {gz_path} is not a gzip file (skipping unzip)")
        return
    out_path = gz_path[:-3]  # Remove the '.gz'
    print(f"Unzipping {gz_path} to {out_path}")
    try:
        with gzip.open(gz_path, 'rb') as f_in, open(out_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        # Optionally, remove the original .gz file by uncommenting the next line:
        # os.remove(gz_path)
    except Exception as e:
        print(f"[ERROR] Failed to unzip {gz_path}: {e}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description="Download NCBI assembly proteome, genomic FASTA, CDS FASTA, and GTF files given a list of assembly entries, and unzip them."
    )
    parser.add_argument("--assembly_list", required=True,
                        help="File with three columns per line: [AssemblyID] [RenamedBasename] [GroupStatus]")
    parser.add_argument("--out_dir", required=True,
                        help="Output directory for downloads")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.assembly_list, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    for line in lines:
        parts = line.split()
        if len(parts) < 2:
            print(f"[WARN] Skipping line with insufficient columns: {line}", file=sys.stderr)
            continue

        assembly_id = parts[0]
        renamed_basename = parts[1]
        # The third column (group status) is ignored.
        
        print(f"\nProcessing assembly: {assembly_id} (renamed as {renamed_basename})")
        ftp_path = construct_ftp_path(assembly_id)
        if ftp_path is None:
            print(f"[WARN] Skipping assembly {assembly_id} due to FTP path construction error.", file=sys.stderr)
            continue

        # Construct file names using the renamed basename, using periods instead of underscores before the file type suffix.
        proteome_file = f"{renamed_basename}.protein.faa.gz"
        genomic_file  = f"{renamed_basename}.genomic.fna.gz"
        cds_file      = f"{renamed_basename}.cds_from_genomic.fna.gz"
        gtf_file      = f"{renamed_basename}.genomic.gtf.gz"

        # Note: The URLs remain based on the original assembly naming convention.
        proteome_url = f"{ftp_path}/{assembly_id}_protein.faa.gz"
        genomic_url  = f"{ftp_path}/{assembly_id}_genomic.fna.gz"
        cds_url      = f"{ftp_path}/{assembly_id}_cds_from_genomic.fna.gz"
        gtf_url      = f"{ftp_path}/{assembly_id}_genomic.gtf.gz"

        # Create a subdirectory for the renamed assembly.
        assembly_dir = os.path.join(args.out_dir, renamed_basename)
        os.makedirs(assembly_dir, exist_ok=True)

        # Download files
        proteome_path = os.path.join(assembly_dir, proteome_file)
        genomic_path = os.path.join(assembly_dir, genomic_file)
        cds_path = os.path.join(assembly_dir, cds_file)
        gtf_path = os.path.join(assembly_dir, gtf_file)

        download_file(proteome_url, proteome_path)
        download_file(genomic_url, genomic_path)
        download_file(cds_url, cds_path)
        download_file(gtf_url, gtf_path)

        # Unzip the downloaded files
        unzip_file(proteome_path)
        unzip_file(genomic_path)
        unzip_file(cds_path)
        unzip_file(gtf_path)

if __name__ == "__main__":
    main()
