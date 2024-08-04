#!/usr/bin/env python3

import glob
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
import tarfile
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


def download_taxanomy(context):
    cache_dir = context.params.get('cache_dir')
    taxonomy_path = os.path.join(cache_dir, "taxonomy")
    os.makedirs(taxonomy_path, exist_ok=True)
    os.chdir(taxonomy_path)

    if not context.params.get('skip_maps'):
        if not context.params.get('protein'):
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
    cmd = f"find . -name '*.gz' | xargs -n 1 -P {context.params['threads']} gunzip -k"
    run_cmd(cmd)


def run_cmd(cmd):
    logger.info(f"Running command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def download_genomes(context, db_type):
    organisms = DB_TYPE_CONFIG.get(db_type, [db_type])
    k2_db_dir = f"k2_{context.params['db_type']}"
    if context.params.get('force'):
        shutil.rmtree(k2_db_dir, ignore_errors=True)

    os.makedirs(k2_db_dir, exist_ok=True)

    cache_dir = context.cache_dir

    for organism in organisms:
        logger.info(f"Downloading genomes for {organism}")
        os.chdir(cache_dir)

        cmd = f"""ncbi-genome-download --section refseq --format fasta --assembly-level complete --retries 3 --parallel {context.params['threads']} --progress-bar {organism}"""
        logger.info(f"Running command: {cmd}")
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # import ncbi_genome_download as ngd
        # ngd.group = 'fungi'
        # ngd.progess_bar = True
        # ngd.format = 'fasta'
        # ngd.download()

        # check if gunziped exits, if not unzip all files parallely but keep the original
        cmd = f"find refseq/{organism} -name '*.gz' | xargs -n 1 -P {context.params['threads']} gunzip -k"
        # cmd = f"find refseq/{organism} -name '*.gz' | xargs -n 1 -P {context.params['threads']} gunzip -k"
        logger.info(f"Running command: {cmd}")
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Files already present, skipping gunzip")

        #     kraken2-build --db $DB --add-to-library "$file"
        # os.chdir(k2_db_dir)
        cmd = f"find {cache_dir}/refseq/{organism} -name '*.fna' | xargs -n 1 -P {context.params['threads']} kraken2-build --db {k2_db_dir} --add-to-library"
        logger.info(f"Running command: {cmd}")
        subprocess.run(cmd, shell=True, check=True)


def build_db(context):
    k2_db_dir = f"k2_{context.params['db_type']}"
    cache_dir = context.cache_dir

    if not os.path.exists(f"{k2_db_dir}/taxonomy"):
        cmd = f"ln -s {cache_dir}/taxonomy {k2_db_dir}/taxonomy"
        run_cmd(cmd)
    
    os.chdir(cache_dir)
    cmd = f"kraken2-build --db {k2_db_dir} --build --threads {context.params['threads']}"
    logger.info(f"Running command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def create_assembly(context):
    logger.info("Creating assembly")
    k2_db_dir = f"k2_{context.params['db_type']}"
    os.makedirs(k2_db_dir, exist_ok=True)
    # add all fasta files
    cmd = f"find "
    pass


@click.command()
@click.option('--db-type', default='standard', help='database type to build')
@click.option('--threads', default=multiprocessing.cpu_count(), help='Number of threads to use')
@click.option('--cache-dir', default=create_cache_dir(), help='Cache directory')
@click.option('--force', is_flag=True, help='Force download and build')
@click.pass_context
def main(context, db_type: str, threads: int, cache_dir, force: bool):
    logger.info(f"Building Kraken2 database of type {db_type}")
    run_basic_checks()

    logger.info(f"Using cache directory {context.params.get('cache_dir')}")

    download_taxanomy(context)
    download_genomes(context, db_type)
    build_db(context)

    # db_types = [c.strip() for c in db_type.split(',')]



if __name__ == '__main__':
    main()
