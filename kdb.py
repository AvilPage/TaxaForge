#!/usr/bin/env python3
import datetime
import hashlib
import logging
import multiprocessing
import os
import shutil
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

    cmd = f"echo {' '.join(urls)} | xargs -n 1 -P 4 wget -q -c"
    run_cmd(cmd, no_output=True)

    cmd = f"tar -k -xvf taxdump.tar.gz"
    run_cmd(cmd, no_output=True)

    logger.info("Decompressing taxonomy data")
    cmd = f"find {cache_dir}/taxonomy -name '*.gz' | xargs -n 1 gunzip -k"
    run_cmd(cmd, no_output=True)

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

    for organism in organisms:
        logger.info(f"Downloading genomes for {organism}")
        os.chdir(cache_dir)
        ncbi_genome_download.download(
            section='refseq', groups=organism, file_formats='fasta',
            progress_bar=True, parallel=threads,
            assembly_levels=['complete'],
            output=cache_dir
        )

        cmd = f"find {cache_dir}/refseq/{organism} -name '*.gz' | xargs -n 1 -P {threads} gunzip -k"
        run_cmd(cmd)
        logger.info(f"Finished downloading {organism} genomes")

    os.chdir(cwd)
    logger.info("Finished downloading all genomes")


def build_db(
        cache_dir, cwd, db_type, db_name, threads, kmer_len,
        fast_build, rebuild, load_factor, use_k2
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

    if use_k2:
        cmd = f"k2 build"
    else:
        cmd = f"kraken2-build --build"

    cmd += f" --db {db_name} --threads {threads} --kmer-len {kmer_len} --load-factor {load_factor}"
    if fast_build:
        cmd += " --fast-build"

    run_cmd(cmd)

    cmd = f"du -sh {db_name}/*.k2d"
    run_cmd(cmd)


def get_files(genomes_dir, cache_dir, db_type, db_name, threads):
    if genomes_dir:
        logger.info(f"Adding {genomes_dir} genomes to library")

        # cmd = f"find {genomes_dir} -name '*.gz' | xargs -n 1 -P {threads} gunzip -k"
        # run_cmd(cmd)

        cmd = f"find {genomes_dir} -type f -name '*.fna'"
        files = run_cmd(cmd, return_output=True)
        logger.info(f"Found {len(files)} genomes to add to {db_name} library")
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
    with open(md5_file, "w") as out_file:
        for line in hashes:
            out_file.write(line + "\n")
    logger.info(f"Saved {len(hashes)} md5 hashes")


def add_to_library(
        cache_dir, cwd, genomes_dir, db_type, db_name,
        limit, batch_size, threads, use_k2
):
    os.chdir(cwd)
    os.makedirs(cwd / db_name / "library", exist_ok=True)

    files = get_files(genomes_dir, cache_dir, db_type, db_name, threads)
    if limit:
        logger.info(f"Limiting number of genomes to {limit}")
        files = files[:limit]

    step = batch_size
    dynamic_step = len(files) // 10
    step = min(step, dynamic_step)
    if step == 0:
        step = 1

    logger.info(f"Using step size of {step}")

    file_count = len(files)
    start = datetime.datetime.now()

    if use_k2:
        for index, file in enumerate(files, start=1):
            if index % step == 0:
                duration = datetime.datetime.now() - start
                average_speed = duration / step
                eta = (file_count - index) * average_speed
                logger.info(f"{datetime.datetime.now()}: Added {index} genomes in {duration}. ETA: {eta}")
                start = datetime.datetime.now()

            cmd = f"k2 add-to-library --db {db_name} --files {file}"
            run_cmd(cmd, no_output=True)

        logger.info(f"Added downloaded genomes to library")
        end = datetime.datetime.now()
        print(f"Time taken: {end - start}")
        return

    global hashes
    global md5_file
    md5_file = cwd / db_name / "library" / "added.md5"

    if os.path.exists(md5_file):
        with open(md5_file, "r") as in_file:
            hashes = {line.strip() for line in in_file}

        logger.info(f"Found {len(hashes)} md5 hashes in {md5_file}")

    for index, file in enumerate(files, start=1):
        if index % step == 0:
            duration = datetime.datetime.now() - start
            average_speed = duration / step
            eta = (file_count - index) * average_speed
            logger.info(f"{datetime.datetime.now()}: Added {index} genomes in {duration}. ETA: {eta}")
            start = datetime.datetime.now()

        if not os.path.exists(f"{file}.md5"):
            md5sum = hash_file(file)
            with open(f"{file}.md5", "w") as fh:
                fh.write(md5sum)
        else:
            with open(f"{file}.md5", "r") as in_file:
                md5sum = in_file.read()

        if md5sum in hashes:
            continue

        cmd = f"kraken2-build --db {db_name} --add-to-library {file} --threads {threads}"
        run_cmd(cmd, no_output=True)

        with open(md5_file, "a") as out_file:
            out_file.write(md5sum + "\n")

        hashes.add(md5sum)

    end = datetime.datetime.now()
    print(f"Time taken: {end - start}")

    logger.info(f"Added downloaded genomes to library")


@click.command()
@click.option('--db-type', default=None, help='database type to build')
@click.option('--db-name', default=None, help='database name to build')
@click.option('--genomes-dir', default=None, help='Directory containing genomes')
@click.option('--cache-dir', default=create_cache_dir(), help='Cache directory')
@click.option('--threads', default=multiprocessing.cpu_count(), help='Number of threads to use', type=int)
@click.option('--load-factor', default=0.7, help='Proportion of the hash table to be populated')
@click.option('--kmer-len', default=35, help='Kmer length in bp/aa. Used only in build task', type=int)
@click.option('--limit', default=None, help='Limit number of genomes to use', type=int)
@click.option('--batch-size', default=1000, help='Number of genomes to add to library at a time', type=int)
@click.option('--force', is_flag=True, help='Force download and build')
@click.option('--rebuild', is_flag=True, help='Clean existing build files and re-build')
@click.option('--fast-build', is_flag=True, help='Non deterministic but faster build')
@click.option('--use-k2', is_flag=True, help='Non deterministic but faster build')
@click.pass_context
def main(
        context,
        db_type: str, db_name, cache_dir, genomes_dir,
        threads, load_factor, kmer_len: int, limit: int, batch_size: int,
        force: bool, rebuild, fast_build: bool, use_k2: bool
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

    add_to_library(
        cache_dir, cwd, genomes_dir, db_type, db_name,
        limit, batch_size, threads, use_k2
    )
    build_db(
        cache_dir, cwd, db_type, db_name, threads, kmer_len,
        fast_build, rebuild, load_factor, use_k2
    )


if __name__ == '__main__':
    main()
