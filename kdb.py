#!/usr/bin/env python3

import glob
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
import tarfile
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from pathlib import Path

import click

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


NCBI_SERVER = "https://ftp.ncbi.nlm.nih.gov"


DB_TYPE_CONFIG = {
    'standard': ("archaea", "bacteria", "viral", "plasmid", "human", "UniVec_Core")
}


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


def run_cmd(cmd):
    logger.info(f"Running command: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        pass


def download_genomes(cache_dir, db_type, db_name, threads, force=False):
    organisms = DB_TYPE_CONFIG.get(db_type, [db_type])
    if force:
        shutil.rmtree(db_name, ignore_errors=True)

    os.makedirs(db_name, exist_ok=True)
    os.chdir(cache_dir)

    for organism in organisms:
        logger.info(f"Downloading genomes for {organism}")

        cmd = f"""ncbi-genome-download --section refseq --format fasta --assembly-level complete --retries 3 --parallel {threads} --progress-bar {organism}"""
        logger.info(f"Running command: {cmd}")
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        cmd = f"find {cache_dir}/refseq/{organism} -name '*.gz' | xargs -n 1 -P {threads} gunzip -k"
        run_cmd(cmd)

        cmd = f"find {cache_dir}/refseq/{organism} -name '*.fna' | xargs -n 1 -P {threads} kraken2-build --db {k2_db_dir} --add-to-library"
        run_cmd(cmd)
        logger.info(f"Finished downloading {organism} genomes")

    logger.info("Finished downloading all genomes")


def build_db(cache_dir, db_name, threads):
    if not os.path.exists(f"{db_name}/taxonomy"):
        cmd = f"ln -s {cache_dir}/taxonomy {db_name}/taxonomy"
        run_cmd(cmd)
    
    os.chdir(cache_dir)
    cmd = f"kraken2-build --db {db_name} --build --threads {threads}"
    logger.info(f"Running command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


@click.command()
@click.option('--db-type', default='standard', help='database type to build')
@click.option('--threads', default=multiprocessing.cpu_count(), help='Number of threads to use')
@click.option('--cache-dir', default=create_cache_dir(), help='Cache directory')
@click.option('--force', is_flag=True, help='Force download and build')
@click.pass_context
def main(context, db_type: str, threads: int, cache_dir, force: bool):
    logger.info(f"Building Kraken2 database of type {db_type}")
    run_basic_checks()

    if cache_dir == '.':
        cache_dir = os.getcwd()
        print(cache_dir)

    db_name = f"k2_{context.params['db_type']}"

    logger.info(f"Using cache directory {cache_dir}")

    with ProcessPoolExecutor(max_workers=2) as executor:
        future1 = executor.submit(download_taxanomy, cache_dir)
        future2 = executor.submit(download_genomes, cache_dir, db_type, db_name, threads, force)

        future1.result()
        future2.result()

    build_db(cache_dir, db_name, threads)


if __name__ == '__main__':
    main()
