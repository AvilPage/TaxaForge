#!/usr/bin/env python3

import hashlib
import json
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
from pathlib import Path

import click
import ncbi_genome_download

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


NCBI_SERVER = "https://ftp.ncbi.nlm.nih.gov"


DB_TYPE_CONFIG = {
    'standard': ("archaea", "bacteria", "viral", "plasmid", "human", "UniVec_Core")
}


def hash_file(filename, buf_size=8192):
    md5 = hashlib.md5()
    with open(filename, "rb") as in_file:
        while True:
            data = in_file.read(buf_size)
            if not data:
                break
            md5.update(data)
    digest = md5.hexdigest()
    return digest


def run_basic_checks():
    if not shutil.which("kraken2-build"):
        logger.error("kraken2-build not found in PATH. Exiting.")
        sys.exit(1)

    if not shutil.which("ncbi-genome-download"):
        logger.error("ncbi-genome-download not found in PATH. Exiting.")
        sys.exit(1)


def create_cache_dir():
    # Unix ~/.cache/kdb
    # macOS ~/Library/Caches/kdb
    if sys.platform == "darwin":
        cache_dir = Path.home() / "Library" / "Caches" / "kdb"
    if sys.platform == "linux":
        cache_dir = Path.home() / ".cache" / "kdb"

    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def download_taxanomy(cache_dir, skip_maps=None, protein=None):
    taxonomy_path = os.path.join(cache_dir, "taxonomy")
    os.makedirs(taxonomy_path, exist_ok=True)
    os.chdir(taxonomy_path)

    if not skip_maps:
        if not protein:
            # Define URLs for nucleotide accession to taxon map
            urls = [
                f"{NCBI_SERVER}/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
                f"{NCBI_SERVER}/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz"
            ]
        else:
            # Define URL for protein accession to taxon map
            urls = ["ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"]
    else:
        logger.info("Skipping maps download")

    # Download taxonomy tree data
    urls.append(f"{NCBI_SERVER}/pub/taxonomy/taxdump.tar.gz")

    cmd = f"echo {' '.join(urls)} | xargs -n 1 -P 4 wget -c"
    subprocess.run(cmd, shell=True, check=True)

    logger.info("Untarring taxonomy tree data")
    cmd = f"tar -k -xvf taxdump.tar.gz"
    run_cmd(cmd)

    logger.info("Decompressing taxonomy data")
    cmd = f"find {cache_dir}/taxonomy -name '*.gz' | xargs -n 1 gunzip -k"
    run_cmd(cmd)
    logger.info("Finished downloading taxonomy data")


def run_cmd(cmd):
    logger.info(f"Running command: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        pass


def download_genomes(cache_dir, cwd, db_type, db_name, threads, force=False):
    organisms = DB_TYPE_CONFIG.get(db_type, [db_type])
    if force:
        shutil.rmtree(cwd / db_name, ignore_errors=True)

    os.makedirs(cwd / db_name, exist_ok=True)
    os.chdir(cwd)

    for organism in organisms:
        logger.info(f"Downloading genomes for {organism}")
        os.chdir(cache_dir)
        ncbi_genome_download.download(
            section='refseq', groups=organism, file_formats='fasta',
            progress_bar=True, parallel=threads,
            assembly_levels=['complete'],
            # retries=3
        )

        cmd = f"find {cache_dir}/refseq/{organism} -name '*.gz' | xargs -n 1 -P {threads} gunzip -k"
        run_cmd(cmd)
        logger.info(f"Finished downloading {organism} genomes")

        # os.chdir(cwd)
        # cmd = f"find {cache_dir}/refseq/{organism} -name '*.fna' | xargs -n 1 -P {threads} kraken2-build --db {db_name} --add-to-library"
        # cmd = f"find {cache_dir}/refseq/{organism} -name '*.fna' | xargs -n 1 -P {threads} k2 add-to-library --db {db_name} --files"
        # run_cmd(cmd)

    logger.info("Finished downloading all genomes")


def build_db(
        cache_dir, cwd, db_type, db_name, threads, kmer_len,
        fast_build, rebuild, load_factor
):
    os.chdir(cwd)

    if not os.path.exists(f"{db_name}/taxonomy"):
        cmd = f"ln -s {cache_dir}/taxonomy {db_name}/taxonomy"
        run_cmd(cmd)

    if rebuild:
        cmd = f"rm -rf {db_name}/*.k2d"
        run_cmd(cmd)

    # TODO: Fix issue with macos threads
    if sys.platform == "darwin":
        threads = 1

    cmd = f"kraken2-build --build --db {db_name} --threads {threads} --load-factor {load_factor} --kmer-len {kmer_len}"

    if fast_build:
        cmd += " --fast-build"

    run_cmd(cmd)

    cmd = f"du -sh {db_name}/*.k2d"
    run_cmd(cmd)


def add_to_library(cache_dir, cwd, genomes_dir, db_type, db_name, threads):
    files_to_add = []

    hashes = {}
    md5_file = cwd / db_name / "library" / "added.md5"
    if os.path.exists(md5_file):
        with open(md5_file, "r") as in_file:
            hashes = json.load(in_file)

    if genomes_dir:
        files = subprocess.check_output(f"find {genomes_dir} -name '*.fna'", shell=True).decode("utf-8").split("\n")
    else:
        organisms = DB_TYPE_CONFIG.get(db_type, [db_type])
        files = []
        for organism in organisms:
            org_files = subprocess.check_output(f"find {cache_dir}/refseq/{organism} -name '*.fna'", shell=True).decode("utf-8").split("\n")
            files.extend(org_files)

    for file in files:
        if not file:
            continue
        md5sum = hash_file(file)
        if md5sum in hashes:
            continue
        files_to_add.append(file)
        hashes[md5sum] = file

    if not files_to_add:
        logger.info(f"No new genomes to add for")
        return

    os.chdir(cwd)
    batch_size = threads * 10
    for i in range(0, len(files_to_add), batch_size):
        batch = files_to_add[i:i + batch_size]
        cmd = f"echo {' '.join(batch)} | xargs -n 1 -P {threads} kraken2-build --db {db_name} --add-to-library"
        run_cmd(cmd)

        with open(md5_file, "w") as out_file:
            json.dump(hashes, out_file)

    logger.info(f"Added downloaded genomes to library")


@click.command()
@click.option('--db-type', default=None, help='database type to build', required=True)
@click.option('--db-name', default=None, help='database name to build')
@click.option('--genomes-dir', default=None, help='Directory containing genomes')
@click.option('--cache-dir', default=create_cache_dir(), help='Cache directory')
@click.option('--threads', default=multiprocessing.cpu_count(), help='Number of threads to use')
@click.option('--load-factor', default=0.7, help='Proportion of the hash table to be populated')
@click.option('--kmer-len', default=31, help='Kmer length in bp/aa. Used only in build task')
@click.option('--force', is_flag=True, help='Force download and build')
@click.option('--rebuild', is_flag=True, help='Clean existing build files and re-build')
@click.option('--fast-build', is_flag=True, help='Non deterministic but faster build')
@click.pass_context
def main(
        context,
        db_type: str, db_name, cache_dir, genomes_dir,
        threads, load_factor, kmer_len: int,
        force: bool, rebuild, fast_build: bool
):
    logger.info(f"Building Kraken2 database of type {db_type}")
    run_basic_checks()
    cwd = Path(os.getcwd())

    if cache_dir == '.':
        cache_dir = cwd

    if not db_name:
        db_name = f"k2_{context.params['db_type']}"

    logger.info(f"Using cache directory {cache_dir}")

    download_taxanomy(cache_dir)
    if not genomes_dir:
        download_genomes(cache_dir, cwd, db_type, db_name, threads, force)
    add_to_library(cache_dir, cwd, genomes_dir, db_type, db_name, threads)
    build_db(
        cache_dir, cwd, db_type, db_name, threads, kmer_len,
        fast_build, rebuild, load_factor
    )


if __name__ == '__main__':
    main()
