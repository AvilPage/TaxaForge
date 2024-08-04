#!/usr/bin/env python3
import hashlib
import json
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import click

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
    os.chdir(cache_dir)

    hashes = {}
    md5_file = cwd / db_name / "library" / "added.md5"
    if os.path.exists(md5_file):
        with open(md5_file, "r") as in_file:
            hashes = json.load(in_file)

    for organism in organisms:
        logger.info(f"Downloading genomes for {organism}")

        cmd = f"""ncbi-genome-download --section refseq --format fasta --assembly-level complete --retries 3 --parallel {threads} --progress-bar {organism}"""
        logger.info(f"Running command: {cmd}")
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        cmd = f"find {cache_dir}/refseq/{organism} -name '*.gz' | xargs -n 1 -P {threads} gunzip -k"
        run_cmd(cmd)
        logger.info(f"Finished downloading {organism} genomes")

        cmd = f"find {cache_dir}/refseq/{organism} -name '*.fna' | xargs -n 1 -P {threads} kraken2-build --db {db_name} --add-to-library"

        files_to_add = []

        files = subprocess.check_output(f"find {cache_dir}/refseq/{organism} -name '*.fna'", shell=True).decode("utf-8").split("\n")

        for file in files:
            if not file:
                continue
            md5sum = hash_file(file)
            if md5sum in hashes:
                continue
            files_to_add.append(file)
            hashes[md5sum] = file

        os.chdir(cwd)
        cmd = f"echo {' '.join(files_to_add)} | xargs -n 1 -P {threads} kraken2-build --db {db_name} --add-to-library"
        run_cmd(cmd)

        with open(md5_file, "w") as out_file:
            json.dump(hashes, out_file)

        logger.info(f"Added downloaded {organism} genomes to library")

    logger.info("Finished downloading all genomes")


def build_db(cache_dir, cwd, db_name, threads):
    if not os.path.exists(f"{db_name}/taxonomy"):
        cmd = f"ln -s {cache_dir}/taxonomy {db_name}/taxonomy"
        run_cmd(cmd)
    
    cmd = f"export OMP_NUM_THREADS={threads}; export KRAKEN2_NUM_THREADS={threads}; kraken2-build --build --db {db_name} --threads {threads}"
    run_cmd(cmd)


@click.command()
@click.option('--db-type', default=None, help='database type to build', required=True)
@click.option('--threads', default=multiprocessing.cpu_count(), help='Number of threads to use')
@click.option('--cache-dir', default=create_cache_dir(), help='Cache directory')
@click.option('--force', is_flag=True, help='Force download and build')
@click.pass_context
def main(context, db_type: str, threads: int, cache_dir, force: bool):
    logger.info(f"Building Kraken2 database of type {db_type}")
    run_basic_checks()
    cwd = Path(os.getcwd())

    if cache_dir == '.':
        cache_dir = cwd
        print(cache_dir)

    db_name = f"k2_{context.params['db_type']}"

    logger.info(f"Using cache directory {cache_dir}")

    # with ProcessPoolExecutor(max_workers=2) as executor:
    #     future1 = executor.submit(download_taxanomy, cache_dir)
    #     future2 = executor.submit(download_genomes, cache_dir, db_type, db_name, threads, force)

        # future1.result()
        # future2.result()

    download_taxanomy(cache_dir)
    download_genomes(cache_dir, cwd, db_type, db_name, threads, force)
    build_db(cache_dir, cwd, db_name, threads)


if __name__ == '__main__':
    main()
