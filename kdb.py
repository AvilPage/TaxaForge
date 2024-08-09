#!/usr/bin/env python3
import atexit
import hashlib
import logging
import multiprocessing
import os
import shutil
import signal
import subprocess
import sys
from pathlib import Path

import click
import ncbi_genome_download
from tqdm import tqdm

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


NCBI_SERVER = "https://ftp.ncbi.nlm.nih.gov"


DB_TYPE_CONFIG = {
    'standard': ("archaea", "bacteria", "viral", "plasmid", "human", "UniVec_Core")
}
hashes = set()
md5_file = None


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


def run_cmd(cmd, return_output=False, no_output=False):
    if not no_output:
        logger.info(f"Running command: {cmd}")

    if return_output:
        return subprocess.check_output(cmd, shell=True).decode("utf-8").strip().split("\n")

    try:
        if no_output:
            subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
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
        )

        cmd = f"find {cache_dir}/refseq/{organism} -name '*.gz' | xargs -n 1 -P {threads} gunzip -k"
        run_cmd(cmd)
        logger.info(f"Finished downloading {organism} genomes")

    logger.info("Finished downloading all genomes")


def build_db(
        cache_dir, cwd, db_type, db_name, threads, kmer_len,
        fast_build, rebuild, load_factor
):
    run_cmd(f"cd {cwd}")

    if not os.path.exists(f"{db_name}/taxonomy"):
        cmd = f"ln -s {cache_dir}/taxonomy {db_name}/"
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


def get_files(genomes_dir, cache_dir, db_type):
    if genomes_dir:
        logger.info(f"Adding {genomes_dir} genomes to library")
        cmd = f"find {genomes_dir} -type f -name '*.fna'"
        files = run_cmd(cmd, return_output=True)
        logger.info(f"Found {len(files)} genomes to add")
    else:
        organisms = DB_TYPE_CONFIG.get(db_type, [db_type])
        files = []
        for organism in organisms:
            cmd = f"find {cache_dir}/refseq/{organism} -name '*.fna'"
            org_files = run_cmd(cmd, return_output=True)
            logger.info(f"Found {len(org_files)} genomes for {organism}")
            files.extend(org_files)

    return files


def save_md5_file(*args, **kwargs):
    global md5_file
    # save set of md5 hashes to file
    with open(md5_file, "w") as out_file:
        for line in hashes:
            out_file.write(line + "\n")
    logger.info(f"Saved {len(hashes)} md5 hashes")

def add_to_library(cache_dir, cwd, genomes_dir, db_type, db_name, threads):
    # atexit.register(save_md5_file)
    # signal.signal(signal.SIGTERM, save_md5_file)
    # signal.signal(signal.SIGINT, save_md5_file)

    run_cmd(f"cd {cwd}")

    global hashes
    global md5_file
    md5_file = cwd / db_name / "library" / "added.md5"
    os.makedirs(cwd / db_name / "library", exist_ok=True)

    if os.path.exists(md5_file):
        with open(md5_file, "r") as in_file:
            hashes = {line.strip() for line in in_file}

        logger.info(f"Found {len(hashes)} md5 hashes in {md5_file}")

    files = get_files(genomes_dir, cache_dir, db_type)

    for file in tqdm(files):
        if not os.path.exists(f"{file}.md5"):
            md5sum = hash_file(file)
            with open(f"{file}.md5", "w") as fh:
                fh.write(md5sum)
        else:
            with open(f"{file}.md5", "r") as in_file:
                md5sum = in_file.read()

        if md5sum in hashes:
            continue

        cmd = f"kraken2-build --db {db_name} --add-to-library {file}"
        run_cmd(cmd, no_output=True)

        # append md5sum to file
        with open(md5_file, "a") as out_file:
            out_file.write(md5sum + "\n")

        hashes.add(md5sum)

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

    if force:
        run_cmd(f"rm -rf {db_name}")
        run_cmd(f"mkdir -p {db_name}")

    logger.info(f"Using cache directory {cache_dir}")

    download_taxanomy(cache_dir)

    if not genomes_dir:
        download_genomes(cache_dir, cwd, db_type, db_name, threads, force)
    try:
        add_to_library(cache_dir, cwd, genomes_dir, db_type, db_name, threads)
    except KeyboardInterrupt:
        print("Exiting")
        # save_md5_file()
        sys.exit(1)
    build_db(
        cache_dir, cwd, db_type, db_name, threads, kmer_len,
        fast_build, rebuild, load_factor
    )


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Exiting")
        save_md5_file()
        sys.exit(1)
