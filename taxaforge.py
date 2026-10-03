#!/usr/bin/env python3
import concurrent.futures
import configparser
import csv
import datetime
import gzip
import hashlib
import importlib.metadata
import itertools
import logging
import multiprocessing
import os
import re
import shlex
import shutil
import socket
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.error
import urllib.request
import zlib
from collections import deque
from pathlib import Path

import click
import ncbi_genome_download
from tqdm import tqdm

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


NCBI_SERVER = "https://ftp.ncbi.nlm.nih.gov"
DOWNLOAD_ATTEMPTS = 3
MAX_DOWNLOAD_THREADS = 4
KRAKEN2_DEFAULT_MINIMIZER_SPACES = 7
DEFAULT_EXTRACTION_WORKERS = max(1, min(8, multiprocessing.cpu_count()))
DEFAULT_LIBRARY_REPORT_WORKERS = max(
    1, min(8, multiprocessing.cpu_count() - 4)
)
EXTRACTION_CHUNK_SIZE = 1000
LIBRARY_REPORT_BATCH_SIZE = 500
ACCESSION_MAP_PARALLEL_MIN_BYTES = 16 * 1024 * 1024


DB_TYPE_CONFIG = {
    'standard': ("archaea", "bacteria", "viral", "plasmid", "human", "UniVec_Core")
}
REQUIRED_BINS = {
    'kraken2': "kraken2-build",
    'ganon2': "ganon",
    'ganon': "ganon",
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


def parse_max_db_size(context, parameter, value):
    if value is None:
        return None

    value = value.upper()
    if not re.fullmatch(r"[1-9]\d*(?:[KMGTP]?)", value):
        raise click.BadParameter(
            "must be a positive byte count optionally suffixed with K, M, G, T, or P"
        )
    return value


def run_basic_checks(tool, use_k2=False):
    if not shutil.which("ncbi-genome-download"):
        logger.error("ncbi-genome-download not found in PATH. Exiting.")
        sys.exit(1)

    if tool not in REQUIRED_BINS:
        logger.error(f"Unknown tool: {tool}. Supported tools: {', '.join(REQUIRED_BINS)}")
        sys.exit(1)

    binary = "k2" if (tool == 'kraken2' and use_k2) else REQUIRED_BINS[tool]
    if not shutil.which(binary):
        logger.error(f"{binary} not found in PATH. Exiting.")
        sys.exit(1)


def create_cache_dir():
    # Unix ~/.cache/taxaforge
    # macOS ~/Library/Caches/taxaforge
    if sys.platform == "darwin":
        cache_dir = Path.home() / "Library" / "Caches" / "taxaforge"
    if sys.platform == "linux":
        cache_dir = Path.home() / ".cache" / "taxaforge"

    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def create_config_dir():
    # Unix ~/.config/taxaforge
    # macOS ~/Library/Application Support/taxaforge
    if sys.platform == "darwin":
        config_dir = Path.home() / "Library" / "Application Support" / "taxaforge"
    if sys.platform == "linux":
        config_dir = Path.home() / ".config" / "taxaforge"

    config_dir.mkdir(parents=True, exist_ok=True)
    return config_dir


CONFIG_SECTION = "taxaforge"


def get_config_path():
    return create_config_dir() / "config.ini"


def load_config():
    parser = configparser.ConfigParser()
    parser.read(get_config_path())
    if not parser.has_section(CONFIG_SECTION):
        parser.add_section(CONFIG_SECTION)
    return parser


def save_config(parser):
    with open(get_config_path(), "w") as out_file:
        parser.write(out_file)


def download_file(url, position):
    filename = url.rsplit("/", 1)[-1]
    existing = os.path.getsize(filename) if os.path.exists(filename) else 0

    request = urllib.request.Request(url)
    if existing:
        request.add_header("Range", f"bytes={existing}-")

    for attempt in range(DOWNLOAD_ATTEMPTS):
        try:
            response = urllib.request.urlopen(request, timeout=30)
            break
        except urllib.error.HTTPError as error:
            if existing and error.code == 416:
                return
            raise
        except urllib.error.URLError as error:
            if not isinstance(
                error.reason, (socket.gaierror, TimeoutError, ConnectionError)
            ):
                raise
            if attempt + 1 == DOWNLOAD_ATTEMPTS:
                raise click.ClickException(
                    f"Unable to download {url} after {DOWNLOAD_ATTEMPTS} attempts: "
                    f"{error.reason}"
                ) from error
            time.sleep(2**attempt)

    with response:
        resumed = response.status == 206
        total = int(response.headers.get("Content-Length", 0)) + (existing if resumed else 0)

        with open(filename, "ab" if resumed else "wb") as out_file, tqdm(
            total=total, initial=existing if resumed else 0, unit="B", unit_scale=True,
            desc=filename, position=position, leave=True
        ) as bar:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                out_file.write(chunk)
                bar.update(len(chunk))


def download_files(urls, max_workers=4):
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(download_file, url, index % max_workers)
            for index, url in enumerate(urls)
        ]
        for future in concurrent.futures.as_completed(futures):
            future.result()


def validate_gzip_file(compressed_path):
    with gzip.open(compressed_path, "rb") as compressed_file:
        while compressed_file.read(1024 * 1024):
            pass


def decompress_gzip_file_atomically(compressed_path):
    decompressed_path = compressed_path.with_suffix("")
    temporary_path = decompressed_path.with_name(
        f"{decompressed_path.name}.partial"
    )
    try:
        with gzip.open(compressed_path, "rb") as compressed_file, temporary_path.open(
            "wb"
        ) as decompressed_file:
            shutil.copyfileobj(compressed_file, decompressed_file)
        temporary_path.replace(decompressed_path)
    except (EOFError, gzip.BadGzipFile, zlib.error):
        temporary_path.unlink(missing_ok=True)
        raise


def download_taxanomy(cache_dir, skip_maps=None, protein=None):
    taxonomy_path = os.path.join(cache_dir, "taxonomy")
    os.makedirs(taxonomy_path, exist_ok=True)
    original_cwd = os.getcwd()
    os.chdir(taxonomy_path)

    try:
        urls = []
        map_urls = []
        if not skip_maps:
            if not protein:
                # Define URLs for nucleotide accession to taxon map
                map_urls = [
                    (
                        f"{NCBI_SERVER}/pub/taxonomy/accession2taxid/"
                        "nucl_gb.accession2taxid.gz"
                    ),
                    (
                        f"{NCBI_SERVER}/pub/taxonomy/accession2taxid/"
                        "nucl_wgs.accession2taxid.gz"
                    ),
                ]
            else:
                # Define URL for protein accession to taxon map
                map_urls = [
                    "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"
                    "prot.accession2taxid.gz"
                ]
            urls.extend(
                url for url in map_urls
                if not os.path.exists(url.rsplit("/", 1)[-1])
            )
        else:
            logger.info("Skipping maps download")

        taxonomy_files = ("nodes.dmp", "names.dmp")
        taxonomy_ready = all(
            os.path.isfile(os.path.join(taxonomy_path, filename))
            for filename in taxonomy_files
        )
        taxdump_path = os.path.join(taxonomy_path, "taxdump.tar.gz")
        if not taxonomy_ready and not os.path.isfile(taxdump_path):
            urls.append(f"{NCBI_SERVER}/pub/taxonomy/taxdump.tar.gz")

        if urls:
            logger.info(f"Downloading {len(urls)} taxonomy files")
            download_files(urls)

        if not taxonomy_ready:
            logger.info("Extracting taxdump.tar.gz")
            with tarfile.open(taxdump_path, "r:gz") as archive:
                archive.extractall(path=taxonomy_path)

        logger.info("Decompressing taxonomy data")
        for compressed_path in Path(taxonomy_path).glob("*.gz"):
            if compressed_path.name == "taxdump.tar.gz":
                continue
            decompressed_path = compressed_path.with_suffix("")
            completion_marker = decompressed_path.with_name(
                f"{decompressed_path.name}.complete"
            )
            if decompressed_path.exists() and completion_marker.exists():
                continue

            try:
                if decompressed_path.exists():
                    validate_gzip_file(compressed_path)
                else:
                    decompress_gzip_file_atomically(compressed_path)
            except (EOFError, gzip.BadGzipFile, zlib.error):
                decompressed_path.unlink(missing_ok=True)
                completion_marker.unlink(missing_ok=True)
                source_url = next(
                    (url for url in map_urls if url.endswith(compressed_path.name)),
                    None,
                )
                if source_url is None:
                    raise click.ClickException(
                        f"Compressed taxonomy file is corrupt: {compressed_path}"
                    )
                logger.warning(
                    f"Incomplete taxonomy map {compressed_path.name}; "
                    "resuming its download"
                )
                download_file(source_url, 0)
                try:
                    decompress_gzip_file_atomically(compressed_path)
                except (EOFError, gzip.BadGzipFile, zlib.error) as error:
                    compressed_path.unlink(missing_ok=True)
                    logger.warning(
                        f"Resumed taxonomy map {compressed_path.name} is still "
                        "corrupt; downloading a fresh copy"
                    )
                    try:
                        download_file(source_url, 0)
                        decompress_gzip_file_atomically(compressed_path)
                    except (
                        EOFError, gzip.BadGzipFile, zlib.error
                    ) as retry_error:
                        compressed_path.unlink(missing_ok=True)
                        raise click.ClickException(
                            f"Downloaded taxonomy map is corrupt: {compressed_path}"
                        ) from retry_error
            completion_marker.touch()

        logger.info("Finished downloading taxonomy data")
    finally:
        os.chdir(original_cwd)


def run_cmd(cmd, return_output=False, no_output=False):
    if not no_output:
        logger.info(f"Running command: {cmd}")

    if return_output:
        output = subprocess.check_output(cmd, shell=True).decode("utf-8").splitlines()
        return [line for line in output if line]

    try:
        if no_output:
            subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as error:
        raise click.ClickException(
            f"Command failed with exit code {error.returncode}: {cmd}"
        ) from error


def import_seed_genomes(dest_dir, seed_dirs):
    """Symlink existing *_genomic.fna.gz files from seed_dirs into dest_dir/{accession}/,
    matching ncbi_genome_download's own layout so it skips accessions already present
    (it checks filename + md5 against dest_dir/{accession}/, see has_file_changed)."""
    if not seed_dirs:
        return
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    accession_re = re.compile(r'(GC[AF]_\d+\.\d+)')
    linked = 0
    for seed_dir in seed_dirs:
        seed_path = Path(seed_dir).expanduser()
        if not seed_path.is_dir():
            logger.warning(f"Seed dir {seed_path} not found, skipping")
            continue
        for genome_file in seed_path.rglob('*_genomic.fna.gz'):
            match = accession_re.search(genome_file.name)
            if not match:
                continue
            accession_dir = dest_dir / match.group(1)
            dest_file = accession_dir / genome_file.name
            if dest_file.exists():
                continue
            accession_dir.mkdir(parents=True, exist_ok=True)
            dest_file.symlink_to(genome_file.resolve())
            linked += 1
    if linked:
        logger.info(f"Linked {linked} existing genome files from seed dirs into {dest_dir}")


def _decompress_gzip_file(compressed_path):
    decompressed_path = compressed_path.with_suffix("")
    if decompressed_path.exists():
        return None

    temporary_path = decompressed_path.with_name(
        f".{decompressed_path.name}.{os.getpid()}.tmp"
    )
    try:
        with gzip.open(compressed_path, "rb") as compressed_file:
            with temporary_path.open("wb") as decompressed_file:
                shutil.copyfileobj(compressed_file, decompressed_file)
        temporary_path.replace(decompressed_path)
    except (gzip.BadGzipFile, EOFError, zlib.error) as error:
        temporary_path.unlink(missing_ok=True)
        return compressed_path, error
    return None


def decompress_gzip_files(directory, threads, retry_download=None):
    """Decompress cached files, retrying once if a compressed file is corrupt."""
    directory = Path(directory)
    compressed_paths = list(directory.rglob("*.gz"))
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, threads)) as executor:
        corrupt_files = [
            result
            for result in executor.map(_decompress_gzip_file, compressed_paths)
            if result is not None
        ]

    non_gzip_responses = []
    for compressed_path, _ in corrupt_files:
        with compressed_path.open("rb") as compressed_file:
            prefix = compressed_file.read(256).lstrip(b"\xef\xbb\xbf \t\r\n")
        if prefix.startswith(b"<"):
            non_gzip_responses.append(compressed_path)

    if corrupt_files and retry_download is not None and not non_gzip_responses:
        for compressed_path, _ in corrupt_files:
            compressed_path.unlink(missing_ok=True)
            compressed_path.with_suffix("").unlink(missing_ok=True)
        logger.warning(
            f"Found {len(corrupt_files)} corrupt compressed genome file(s); "
            "removing them and retrying the download once"
        )
        retry_download()
        return decompress_gzip_files(directory, threads)

    if corrupt_files:
        if non_gzip_responses:
            if retry_download is not None:
                for compressed_path in non_gzip_responses:
                    compressed_path.unlink(missing_ok=True)
            examples = ", ".join(str(path) for path in non_gzip_responses[:3])
            remaining = len(non_gzip_responses) - min(len(non_gzip_responses), 3)
            suffix = f" (and {remaining} more)" if remaining else ""
            raise click.ClickException(
                f"Received non-gzip XML/HTML content for "
                f"{len(non_gzip_responses)} genome file(s), for example: "
                f"{examples}{suffix}. NCBI returned an error page instead of "
                "genome data; check NCBI access, proxy, or network connectivity "
                "before retrying. Invalid cached response files were removed."
            )
        details = "; ".join(
            f"{path}: {error}" for path, error in corrupt_files[:5]
        )
        remaining = len(corrupt_files) - min(len(corrupt_files), 5)
        suffix = f"; and {remaining} more" if remaining else ""
        raise click.ClickException(
            f"Unable to decompress {len(corrupt_files)} corrupt gzip file(s): "
            f"{details}{suffix}"
        )


def refresh_incomplete_assembly_checksums(directory):
    """Remove checksum manifests for assemblies that have no downloaded FASTA."""
    refreshed = 0
    for checksum_path in Path(directory).rglob("MD5SUMS"):
        assembly_dir = checksum_path.parent
        if not any(assembly_dir.glob("*.fna")) and not any(
            assembly_dir.glob("*.fna.gz")
        ):
            checksum_path.unlink()
            refreshed += 1
    return refreshed


def ensure_organism_genomes(cache_dir, organism, threads, seed_dirs=None):
    # ncbi_genome_download skips assemblies already present under cache_dir/refseq/{organism},
    # so re-running for the same organism only fetches what's missing
    import_seed_genomes(Path(cache_dir) / "refseq" / organism, seed_dirs)
    logger.info(f"Downloading genomes for {organism}")
    original_cwd = os.getcwd()

    def download():
        return ncbi_genome_download.download(
            section='refseq', groups=organism, file_formats='fasta',
            progress_bar=True, parallel=threads,
            assembly_levels=['complete'],
            output=cache_dir, uri=f"{NCBI_SERVER}/genomes"
        )

    try:
        os.chdir(cache_dir)
        download_result = download()
        organism_dir = Path(cache_dir) / "refseq" / organism
        decompress_gzip_files(organism_dir, threads, retry_download=download)
        genome_files = [
            path for path in organism_dir.rglob("*.fna") if path.is_file()
        ]
        if not genome_files and download_result in (None, 0):
            refreshed = refresh_incomplete_assembly_checksums(organism_dir)
            if refreshed:
                logger.warning(
                    f"Refreshing {refreshed} incomplete NCBI checksum manifest(s) "
                    "and retrying the genome download once"
                )
                download_result = download()
                decompress_gzip_files(organism_dir, threads)
                genome_files = [
                    path for path in organism_dir.rglob("*.fna") if path.is_file()
                ]
        if download_result == 75 and not genome_files:
            raise click.ClickException(
                f"Unable to download {organism} genomes from NCBI, and no cached "
                "FASTA files are available. Check the network/DNS connection and "
                "retry."
            )
        if download_result not in (None, 0, 1, 75):
            raise click.ClickException(
                f"NCBI genome download for {organism} failed with exit code "
                f"{download_result}"
            )
        if not genome_files:
            raise click.ClickException(
                f"No FASTA files were found for {organism} after the download. "
                "Check the NCBI download output and retry."
            )
        if download_result == 75:
            logger.warning(
                f"NCBI is unreachable; continuing with {len(genome_files)} cached "
                f"{organism} FASTA file(s)"
            )
    finally:
        os.chdir(original_cwd)
    logger.info(f"Finished downloading {organism} genomes")


def read_assembly_accessions(accession_file, limit=None):
    """Read exact NCBI assembly accessions, also accepting FASTA filenames."""
    if limit is not None and limit < 1:
        raise ValueError("limit must be a positive integer")

    accession_pattern = re.compile(r'GC[AF]_\d+\.\d+')
    accessions = []
    seen = set()
    with open(accession_file, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            matches = accession_pattern.findall(line)
            filename = Path(line).name
            accession = next(
                (match for match in matches if filename.startswith(match)), None
            )
            if accession is None and len(matches) == 1:
                accession = matches[0]
            if accession is None:
                line_preview = line[:120]
                if len(line) > len(line_preview):
                    line_preview += "..."
                raise ValueError(
                    f"Expected one NCBI assembly accession on line {line_number} "
                    f"of {accession_file}; got {line_preview!r}"
                )
            if accession not in seen:
                accessions.append(accession)
                seen.add(accession)
                if limit is not None and len(accessions) >= limit:
                    break

    if not accessions:
        raise ValueError(f"No NCBI assembly accessions found in {accession_file}")
    return accessions


def download_assembly_accessions(
    cache_dir, accession_file, threads, limit=None, genomes_cache_dir=None
):
    """Download only the assemblies listed in an accession/FASTA list."""
    accessions = read_assembly_accessions(accession_file, limit)
    digest = hashlib.sha256("\n".join(accessions).encode("utf-8")).hexdigest()[:16]
    genome_dir = Path(cache_dir) / "genomes" / f"accessions-{digest}"
    genome_dir.mkdir(parents=True, exist_ok=True)

    cached_accessions = set()
    if genomes_cache_dir:
        cache_path = Path(genomes_cache_dir).expanduser()
        if not cache_path.is_dir():
            raise ValueError(f"Genome cache directory does not exist: {cache_path}")
        accession_set = set(accessions)
        fasta_suffixes = (".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")
        resolved_genome_dir = genome_dir.resolve()
        for genome_path in cache_path.rglob("*"):
            if not genome_path.is_file() or not genome_path.name.endswith(fasta_suffixes):
                continue
            if resolved_genome_dir in genome_path.resolve().parents:
                continue
            match = re.search(r"(GC[AF]_\d+\.\d+)", genome_path.name)
            if not match or match.group(1) not in accession_set:
                continue
            accession = match.group(1)
            target = genome_dir / genome_path.name
            if target.is_symlink() and not target.exists():
                target.unlink()
            if not target.exists():
                target.symlink_to(genome_path.resolve())
            cached_accessions.add(accession)

        if cached_accessions:
            logger.info(
                f"Reusing {len(cached_accessions)} requested assemblies from "
                f"genome cache {cache_path}"
            )

    missing_accessions = [
        accession for accession in accessions if accession not in cached_accessions
    ]
    for section, prefix in (("refseq", "GCF_"), ("genbank", "GCA_")):
        section_accessions = [
            accession for accession in missing_accessions
            if accession.startswith(prefix)
        ]
        if not section_accessions:
            continue

        logger.info(
            f"Downloading {len(section_accessions)} selected {section} assemblies"
        )
        result = ncbi_genome_download.download(
            section=section,
            groups="all",
            assembly_accessions=section_accessions,
            file_formats="fasta",
            assembly_levels="all",
            flat_output=True,
            progress_bar=True,
            parallel=threads,
            output=str(genome_dir),
            uri=f"{NCBI_SERVER}/genomes",
        )
        if result not in (0, 1):
            raise RuntimeError(
                f"NCBI genome download failed for {section} accessions "
                f"(exit code {result})"
            )
        if result == 1:
            logger.info(
                f"No current {section} assemblies matched; checking NCBI summaries"
            )

    downloaded = {
        match.group(1)
        for genome_path in genome_dir.iterdir()
        if (match := re.search(r"(GC[AF]_\d+\.\d+)", genome_path.name))
    }
    missing = [
        accession for accession in missing_accessions if accession not in downloaded
    ]
    if missing:
        summary_sections = {
            "refseq" if accession.startswith("GCF_") else "genbank"
            for accession in missing
        }
        summaries = ensure_assembly_summaries(cache_dir, summary_sections)
        summary_columns = None
        ftp_paths = {}
        missing_set = set(missing)
        for summary_path in summaries:
            with open(summary_path, encoding="utf-8") as summary:
                for line in summary:
                    if line.startswith("#assembly_accession") or line.startswith(
                        "# assembly_accession"
                    ):
                        summary_columns = line.lstrip("#").strip().split("\t")
                        accession_column = summary_columns.index("assembly_accession")
                        ftp_column = summary_columns.index("ftp_path")
                        continue
                    if line.startswith("#") or not line.strip() or summary_columns is None:
                        continue
                    columns = line.rstrip("\n").split("\t")
                    accession = columns[accession_column]
                    if accession in missing_set:
                        ftp_path = columns[ftp_column]
                        if ftp_path != "na":
                            ftp_paths[accession] = ftp_path

        fallback_urls = []
        for accession in missing:
            ftp_path = ftp_paths.get(accession)
            if not ftp_path:
                continue
            ftp_path = ftp_path.replace(
                "ftp://ftp.ncbi.nlm.nih.gov",
                "https://ftp.ncbi.nlm.nih.gov",
            )
            assembly_dir = ftp_path.rstrip("/").rsplit("/", 1)[-1]
            fallback_urls.append(
                f"{ftp_path}/{assembly_dir}_genomic.fna.gz"
            )

        if fallback_urls:
            logger.info(
                f"Downloading {len(fallback_urls)} selected historical assemblies "
                "from NCBI assembly summaries"
            )
            original_cwd = os.getcwd()
            os.chdir(genome_dir)
            try:
                download_files(fallback_urls, max_workers=threads)
            finally:
                os.chdir(original_cwd)

        downloaded = {
            match.group(1)
            for genome_path in genome_dir.iterdir()
            if (match := re.search(r"(GC[AF]_\d+\.\d+)", genome_path.name))
        }
        missing = [
            accession for accession in missing_accessions if accession not in downloaded
        ]

    if missing:
        preview = ", ".join(missing[:10])
        remaining = len(missing) - min(len(missing), 10)
        suffix = f" (and {remaining} more)" if remaining else ""
        raise RuntimeError(
            f"NCBI did not provide {len(missing)} requested assemblies: "
            f"{preview}{suffix}"
        )

    logger.info(
        f"Prepared all {len(accessions)} requested assemblies "
        f"({len(cached_accessions)} reused from cache) in {genome_dir}"
    )
    return str(genome_dir)


def build_accession_taxid_map(summary_paths, accessions=None):
    """Read accession-to-taxid mappings from NCBI assembly summary files."""
    wanted = None if accessions is None else set(accessions)
    accession_to_taxid = {}
    for summary_path in summary_paths:
        columns = None
        with open(summary_path, encoding="utf-8") as summary:
            for line in summary:
                if line.startswith("#assembly_accession") or line.startswith(
                    "# assembly_accession"
                ):
                    columns = [column.strip() for column in line.lstrip("#").strip().split("\t")]
                    required_columns = {"assembly_accession", "taxid"}
                    if not required_columns.issubset(columns):
                        raise ValueError(
                            f"Assembly summary {summary_path} is missing required "
                            "assembly_accession or taxid columns"
                        )
                    accession_column = columns.index("assembly_accession")
                    taxid_column = columns.index("taxid")
                    continue
                if line.startswith("#") or not line.strip() or columns is None:
                    continue
                values = line.rstrip("\n").split("\t")
                if len(values) <= max(accession_column, taxid_column):
                    continue
                accession = values[accession_column]
                taxid = values[taxid_column]
                if (wanted is None or accession in wanted) and taxid.isdigit():
                    accession_to_taxid[accession] = taxid
    return accession_to_taxid


def create_createtaxdb_samplesheet(
    accession_file, genomes_dir, cache_dir, output_file, limit=None
):
    """Create a complete nf-core/createtaxdb samplesheet or fail on missing inputs."""
    accessions = read_assembly_accessions(accession_file, limit)
    genome_dir = Path(genomes_dir).expanduser().resolve()
    if not genome_dir.is_dir():
        raise ValueError(f"Genome directory does not exist: {genome_dir}")

    fasta_suffixes = (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta")
    fasta_by_accession = {}
    for genome_path in genome_dir.rglob("*"):
        if not genome_path.is_file() or not genome_path.name.endswith(fasta_suffixes):
            continue
        match = re.match(r"(GC[AF]_\d+\.\d+)", genome_path.name)
        if match and match.group(1) in accessions:
            fasta_by_accession.setdefault(match.group(1), genome_path.resolve())

    missing_fastas = [accession for accession in accessions if accession not in fasta_by_accession]
    if missing_fastas:
        raise ValueError(
            f"Missing FASTA files for {len(missing_fastas)} requested assemblies: "
            f"{', '.join(missing_fastas[:10])}"
        )

    sections = {
        "refseq" if accession.startswith("GCF_") else "genbank"
        for accession in accessions
    }
    summary_paths = ensure_assembly_summaries(cache_dir, sections)
    accession_to_taxid = build_accession_taxid_map(summary_paths, accessions)
    missing_taxids = [accession for accession in accessions if accession not in accession_to_taxid]
    if missing_taxids:
        raise ValueError(
            f"Missing NCBI taxids for {len(missing_taxids)} requested assemblies: "
            f"{', '.join(missing_taxids[:10])}"
        )

    output_path = Path(output_file).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as samplesheet:
        writer = csv.DictWriter(
            samplesheet,
            fieldnames=["id", "taxid", "fasta_dna", "fasta_aa"],
        )
        writer.writeheader()
        for accession in accessions:
            writer.writerow(
                {
                    "id": accession,
                    "taxid": accession_to_taxid[accession],
                    "fasta_dna": str(fasta_by_accession[accession]),
                    "fasta_aa": "",
                }
            )
    logger.info(f"Wrote {len(accessions)} samples to {output_path}")
    return str(output_path)


def run_createtaxdb(
    nextflow_bin, pipeline, pipeline_revision, output_dir, db_name,
    cache_dir, accession_file, genomes_dir, limit, threads, kmer_len, min_len,
    level, resume,
):
    """Build a ganon index with nf-core/createtaxdb."""
    nextflow_bin = os.path.expanduser(nextflow_bin)
    executable = (
        os.path.isfile(nextflow_bin) and os.access(nextflow_bin, os.X_OK)
        if os.path.sep in nextflow_bin
        else shutil.which(nextflow_bin) is not None
    )
    if not executable:
        raise click.ClickException(f"Nextflow executable not found: {nextflow_bin}")

    taxonomy_dir = Path(cache_dir).expanduser().resolve() / "taxonomy"
    for filename in ("nodes.dmp", "names.dmp"):
        if not (taxonomy_dir / filename).is_file():
            raise click.ClickException(f"Required taxonomy file is missing: {taxonomy_dir / filename}")

    samplesheet_dir = Path(cache_dir).expanduser().resolve() / "samplesheets"
    samplesheet = samplesheet_dir / f"{db_name}.csv"
    try:
        create_createtaxdb_samplesheet(
            accession_file, genomes_dir, cache_dir, samplesheet, limit
        )
    except ValueError as error:
        raise click.ClickException(str(error)) from error

    command = [
        nextflow_bin,
        "run",
        pipeline,
    ]
    if pipeline_revision:
        command.extend(["-r", pipeline_revision])
    command.extend(
        [
            "-profile",
            "docker",
            "--input",
            str(samplesheet),
            "--nodesdmp",
            str(taxonomy_dir / "nodes.dmp"),
            "--namesdmp",
            str(taxonomy_dir / "names.dmp"),
            "--dbname",
            db_name,
            "--build_ganon",
            "--outdir",
            str(Path(output_dir).resolve()),
        ]
    )

    ganon_options = []
    if threads is not None:
        ganon_options.extend(["--threads", str(threads)])
    if kmer_len is not None:
        ganon_options.extend(["--kmer-size", str(kmer_len)])
    if min_len is not None:
        ganon_options.extend(["--window-size", str(min_len)])
    if level and level != "file":
        ganon_options.extend(["--level", level])
    if ganon_options:
        command.append("--ganon_build_options=" + " ".join(ganon_options))
    if resume:
        command.append("-resume")

    logger.info(f"Running nf-core/createtaxdb: {' '.join(shlex.quote(arg) for arg in command)}")
    try:
        subprocess.run(command, check=True, cwd=output_dir)
    except OSError as error:
        raise click.ClickException(f"Nextflow executable not found: {nextflow_bin}") from error
    except subprocess.CalledProcessError as error:
        raise click.ClickException(
            f"nf-core/createtaxdb failed with exit code {error.returncode}"
        ) from error


def download_genomes(cache_dir, cwd, db_type, db_name, threads, force=False, seed_dirs=None):
    organisms = DB_TYPE_CONFIG.get(db_type, [db_type])
    if force:
        shutil.rmtree(cwd / db_name, ignore_errors=True)

    os.makedirs(cwd / db_name, exist_ok=True)

    for organism in organisms:
        ensure_organism_genomes(cache_dir, organism, threads, seed_dirs)

    os.chdir(cwd)
    logger.info("Finished downloading all genomes")


def build_db(
        cache_dir, cwd, db_type, db_name, threads, kmer_len, min_len,
        fast_build, rebuild, load_factor, use_k2,
        minimizer_spaces=KRAKEN2_DEFAULT_MINIMIZER_SPACES,
        max_db_size=None,
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

    cmd += (
        f" --db {db_name} --threads {threads} --kmer-len {kmer_len}"
        f" --minimizer-len {min_len} --minimizer-spaces {minimizer_spaces}"
        f" --load-factor {load_factor}"
    )
    if max_db_size:
        cmd += f" --max-db-size {max_db_size}"
    if fast_build:
        cmd += " --fast-build"

    run_cmd(cmd)

    cmd = f"du -sh {db_name}/*.k2d"
    run_cmd(cmd)


def detect_taxonomy_flag():
    # ganon CLI flag for local taxonomy files has changed across versions; ask the
    # installed binary rather than hardcoding one, so build-custom doesn't fail on --help mismatch
    result = subprocess.run("ganon build-custom --help", shell=True, capture_output=True, text=True)
    help_text = result.stdout + result.stderr
    if "--taxdump-files" in help_text:
        return "--taxdump-files"
    if "--taxonomy-files" in help_text:
        return "--taxonomy-files"
    if "--taxdump-file" in help_text:
        return "--taxdump-file"
    return "--taxonomy-files"


def ensure_assembly_summaries(cache_dir, sections=None):
    # Cached so `ganon build-custom --ncbi-file-info` doesn't re-download these on every build
    taxonomy_path = os.path.join(cache_dir, "taxonomy")
    os.makedirs(taxonomy_path, exist_ok=True)

    summary_urls = {
        "refseq": [
            f"{NCBI_SERVER}/genomes/refseq/assembly_summary_refseq.txt",
            f"{NCBI_SERVER}/genomes/refseq/assembly_summary_refseq_historical.txt",
        ],
        "genbank": [
            f"{NCBI_SERVER}/genomes/genbank/assembly_summary_genbank.txt",
            f"{NCBI_SERVER}/genomes/genbank/assembly_summary_genbank_historical.txt",
        ],
    }
    sections = set(summary_urls) if sections is None else set(sections)
    unknown_sections = sections - summary_urls.keys()
    if unknown_sections:
        raise ValueError(f"Unsupported NCBI assembly summary section(s): {unknown_sections}")

    urls = [url for section in ("refseq", "genbank") if section in sections for url in summary_urls[section]]
    paths = [os.path.join(taxonomy_path, url.rsplit("/", 1)[-1]) for url in urls]
    missing_urls = [url for url, path in zip(urls, paths) if not os.path.exists(path)]

    if missing_urls:
        logger.info(f"Downloading {len(missing_urls)} assembly summary files")
        original_cwd = os.getcwd()
        os.chdir(taxonomy_path)
        try:
            download_files(missing_urls)
        finally:
            os.chdir(original_cwd)

    return paths


def ensure_genome_size_file(cache_dir):
    # Cached so `ganon build-custom --genome-size-files` doesn't try to (re-)download it, which
    # fails outright on hosts without internet access
    taxonomy_path = os.path.join(cache_dir, "taxonomy")
    os.makedirs(taxonomy_path, exist_ok=True)

    path = os.path.join(taxonomy_path, "species_genome_size.txt.gz")
    if not os.path.exists(path):
        logger.info("Downloading species genome size file")
        original_cwd = os.getcwd()
        os.chdir(taxonomy_path)
        try:
            download_files([f"{NCBI_SERVER}/genomes/ASSEMBLY_REPORTS/species_genome_size.txt.gz"])
        finally:
            os.chdir(original_cwd)

    return path


def setup_raptor_shim(shim_dir):
    # Older raptor builds expect --input, but ganon build-custom calls it with --input-file;
    # this shim rewrites the flag so builds don't fail on that mismatch.
    # Only rewrite if the installed raptor actually needs the old flag - newer
    # raptor builds want --input-file, and rewriting there would break them.
    real_raptor = shutil.which("raptor") or "/usr/local/bin/raptor"
    help_text = subprocess.run(
        [real_raptor, "layout", "--help"], capture_output=True, text=True
    ).stdout
    needs_old_flag = "--input-file" not in help_text
    shim_script = shim_dir / "raptor"
    if needs_old_flag:
        content = f"""#!/usr/bin/env bash
REAL_BIN="{real_raptor}"
args=()
for arg in "$@"; do
    if [[ "$arg" == "--input-file" ]]; then
        args+=("--input")
    else
        args+=("$arg")
    fi
done
exec "$REAL_BIN" "${{args[@]}}"
"""
    else:
        content = f"""#!/usr/bin/env bash
exec "{real_raptor}" "$@"
"""
    shim_script.write_text(content)
    shim_script.chmod(0o755)
    return shim_dir


INPUT_EXTENSION_CANDIDATES = ["fna.gz", "fna", "fa.gz", "fa", "fasta.gz", "fasta"]


def detect_input_extension(input_dirs):
    # ganon defaults to fna.gz; genomes-dir may hold uncompressed/differently-named files.
    # Mixed dirs (a few stray .gz among mostly .fna) pick by count, not first match, or the
    # majority format gets silently dropped.
    counts = {ext: 0 for ext in INPUT_EXTENSION_CANDIDATES}
    for ext in INPUT_EXTENSION_CANDIDATES:
        for input_dir in input_dirs:
            path = Path(input_dir)
            if path.is_dir():
                counts[ext] += sum(1 for _ in path.rglob(f"*.{ext}"))
    best_ext = max(counts, key=lambda ext: counts[ext])
    return best_ext if counts[best_ext] else INPUT_EXTENSION_CANDIDATES[0]


def detect_assembly_summary_sections(input_dirs, input_extension):
    sections = set()
    file_count = 0
    for input_dir in input_dirs:
        for path in Path(input_dir).rglob(f"*.{input_extension}"):
            file_count += 1
            match = re.match(r"(GCF|GCA)_\d+\.\d+", path.name)
            if not match:
                return None
            sections.add("refseq" if match.group(1) == "GCF" else "genbank")
    return sections if file_count else None

ACCESSION_RE = re.compile(r"^(GC[AF]_\d+\.\d+)")
RAPTOR_MAX_KMER = 32


def build_accession_metadata_map(assembly_summary_paths):
    # Same as build_accession_taxid_map but also keeps ftp_path (col 20), which points at
    # the NCBI HTTPS/FTP directory for that assembly. nf-core/createtaxdb's fasta_dna column
    # accepts remote URLs directly (Nextflow downloads them), so this lets us build a
    # samplesheet in --remote mode without ever needing the genome files on local disk.
    accession_to_info = {}
    for path in assembly_summary_paths:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) > 19 and cols[0] and cols[5] and cols[19] and cols[19] != "na":
                    accession_to_info[cols[0]] = (cols[5], cols[19])
    return accession_to_info


def ftp_path_to_fasta_url(ftp_path):
    basename = ftp_path.rstrip("/").rsplit("/", 1)[-1]
    return f"{ftp_path.rstrip('/')}/{basename}_genomic.fna.gz"


def read_genome_name_list(input_file):
    names = [
        line.strip()
        for line in Path(input_file).read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    logger.info(f"Read {len(names)} entries from --input-file list {input_file}")
    return names


def write_ganon_input_tsv(files, assembly_summary_paths):
    # NOTE: --input-file does NOT auto-run --ncbi-file-info matching the way --input
    # does - it requires the target/node (taxid) columns to be supplied manually, and
    # silently invalidates every entry otherwise (ganon reports "Unable to match
    # taxonomy to targets"). Passing thousands of files as bare --input args instead
    # would work for small slices but blow past ARG_MAX for large ones (e.g. 600k
    # genomes), so resolve accession -> taxid ourselves and write target/node into the
    # tsv, same info --ncbi-file-info would have extracted.
    accession_to_taxid = build_accession_taxid_map(assembly_summary_paths)
    input_file_list = tempfile.NamedTemporaryFile(mode="w", prefix="ganon_input_", suffix=".tsv", delete=False)
    matched, unmatched = 0, 0
    for f in files:
        match = ACCESSION_RE.match(Path(f).name)
        accession = match.group(1) if match else None
        taxid = accession_to_taxid.get(accession) if accession else None
        if taxid:
            matched += 1
            input_file_list.write(f"{f}\t{accession}\t{taxid}\n")
        else:
            unmatched += 1
            input_file_list.write(f"{f}\n")
    input_file_list.close()
    if unmatched:
        logger.warning(f"{unmatched} of {len(files)} files had no accession/taxid match; they will be skipped by ganon")
    else:
        logger.info(f"Resolved taxid for all {matched} files")
    return input_file_list.name


def resolve_input_file_list(input_file, input_dirs, input_extension):
    # input_file may list bare filenames (e.g. GCF_..._genomic.fna) with no directory
    # component - resolve each against input_dirs so ganon gets real paths. Entries that
    # are already absolute/relative paths pointing at existing files are used as-is.
    names = read_genome_name_list(input_file)

    index = {}
    for input_dir in input_dirs:
        for path in Path(input_dir).rglob(f"*.{input_extension}"):
            index.setdefault(path.name, path)

    resolved, missing = [], []
    for name in names:
        candidate = Path(name)
        if candidate.is_file():
            resolved.append(str(candidate))
            continue
        found = index.get(Path(name).name)
        if found:
            resolved.append(str(found))
        else:
            missing.append(name)
    if missing:
        logger.warning(f"{len(missing)} of {len(names)} entries in {input_file} not found under {input_dirs}")
    logger.info(f"Resolved {len(resolved)} genome files from {input_file}")
    return resolved, missing


def download_missing_genomes(missing_names, genomes_dir, threads):
    # Missing entries are grouped by accession prefix - GCF_ lives under NCBI's refseq
    # section, GCA_ only under genbank - and downloaded straight into genomes_dir so the
    # caller's directory rglob picks them up on re-resolve without any extra copy/move step.
    accessions_by_section = {"refseq": [], "genbank": []}
    unmatched = 0
    for name in missing_names:
        match = ACCESSION_RE.match(Path(name).name)
        if not match:
            unmatched += 1
            continue
        accession = match.group(1)
        section = "refseq" if accession.startswith("GCF_") else "genbank"
        accessions_by_section[section].append(accession)

    if unmatched:
        logger.warning(f"{unmatched} missing entries have no recognizable GCA_/GCF_ accession; cannot auto-download them")

    genomes_dir = Path(genomes_dir)
    genomes_dir.mkdir(parents=True, exist_ok=True)
    for section, accessions in accessions_by_section.items():
        if not accessions:
            continue
        logger.info(f"Downloading {len(accessions)} missing genomes from NCBI {section}")
        result = ncbi_genome_download.download(
            section=section, groups="all", file_formats="fasta",
            assembly_accessions=accessions,
            progress_bar=True, parallel=threads,
            output=str(genomes_dir), uri=f"{NCBI_SERVER}/genomes"
        )
        if result not in (0, 1):
            raise click.ClickException(
                f"NCBI genome download failed for {section} (exit code {result})"
            )
    decompress_gzip_files(genomes_dir, threads)

def run_ganon_build_custom(
    cache_dir, input_dirs, db_name, threads, kmer_len, min_len, level,
    from_idx=None, to_idx=None, restart=False, custom_args=None,
    input_file=None, download_missing=False, download_threads=None,
):
    if kmer_len is not None and int(kmer_len) > RAPTOR_MAX_KMER:
        logger.warning(
            f"--kmer-len {kmer_len} exceeds raptor's max of {RAPTOR_MAX_KMER} "
            f"for ganon2 builds; clamping to {RAPTOR_MAX_KMER}"
        )
        kmer_len = RAPTOR_MAX_KMER
    if (
        kmer_len is not None
        and min_len is not None
        and int(min_len) < int(kmer_len)
    ):
        logger.warning(
            f"--min-len {min_len} is smaller than --kmer-len {kmer_len}; "
            f"raising --min-len to {kmer_len}"
        )
        min_len = kmer_len

    taxa_flag = detect_taxonomy_flag()
    if taxa_flag == "--taxonomy-files":
        tax_args = f"{taxa_flag} {cache_dir}/taxonomy/nodes.dmp {cache_dir}/taxonomy/names.dmp"
    else:
        tax_args = f"{taxa_flag} {cache_dir}/taxonomy/taxdump.tar.gz"

    input_extension = detect_input_extension(input_dirs)
    summary_sections = detect_assembly_summary_sections(input_dirs, input_extension)
    assembly_summary_paths = ensure_assembly_summaries(cache_dir, summary_sections)
    genome_size_file = ensure_genome_size_file(cache_dir)

    # ganon defaults to file-level targets; omit --level for that default.
    level_args = f" --level {level}" if level and level != "file" else ""
    kmer_args = f" --kmer-size {kmer_len}" if kmer_len is not None else ""
    window_args = f" --window-size {min_len}" if min_len is not None else ""
    custom_args_str = " ".join(shlex.quote(arg) for arg in (custom_args or []))

    input_file_list_path = None
    ncbi_file_info_args = f"--ncbi-file-info {' '.join(assembly_summary_paths)} "
    if input_file is not None:
        files, missing = resolve_input_file_list(input_file, input_dirs, input_extension)
        if missing and download_missing:
            download_missing_genomes(
                missing, input_dirs[0], download_threads or threads
            )
            files, missing = resolve_input_file_list(input_file, input_dirs, input_extension)
            if missing:
                logger.warning(f"{len(missing)} entries still missing after download attempt; they will be skipped")
        input_file_list_path = write_ganon_input_tsv(files, assembly_summary_paths)
        input_args = f"--input-file {shlex.quote(input_file_list_path)}"
        # target/node already supplied per-row above; --ncbi-file-info would be ignored
        # (and only slows things down re-parsing assembly summaries) for matched rows.
        ncbi_file_info_args = ""
    elif from_idx is not None or to_idx is not None:
        files = sorted(
            str(path) for input_dir in input_dirs
            for path in Path(input_dir).rglob(f"*.{input_extension}")
        )
        logger.info(f"Using genome files range [{from_idx}:{to_idx}] out of {len(files)}")
        files = files[from_idx:to_idx]

        input_file_list_path = write_ganon_input_tsv(files, assembly_summary_paths)
        input_args = f"--input-file {shlex.quote(input_file_list_path)}"
        # target/node already supplied per-row above; --ncbi-file-info would be ignored
        # (and only slows things down re-parsing assembly summaries) for matched rows.
        ncbi_file_info_args = ""
    else:
        input_args = f"--input {' '.join(input_dirs)} --input-recursive --input-extension {input_extension}"

    cmd = (
        f"ganon build-custom {input_args} "
        f"{tax_args} --taxonomy ncbi "
        f"{ncbi_file_info_args}"
        f"--genome-size-files {genome_size_file} "
        f"--db-prefix {db_name} --threads {threads}"
        f"{kmer_args}{window_args}{level_args}"
        f"{' ' + custom_args_str if custom_args_str else ''}"
        f"{' --restart' if restart else ''}"
    )

    with tempfile.TemporaryDirectory(prefix="raptor_shim_") as shim_tmp:
        setup_raptor_shim(Path(shim_tmp))
        old_path = os.environ.get("PATH", "")
        os.environ["PATH"] = f"{shim_tmp}:{old_path}"
        try:
            run_cmd(cmd)
        finally:
            os.environ["PATH"] = old_path
            if input_file_list_path:
                os.unlink(input_file_list_path)

    cmd = f"du -sh {db_name}.*"
    run_cmd(cmd)


def build_ganon2(
    cache_dir, cwd, genomes_dir, db_type, db_name, threads, kmer_len, min_len,
    level, rebuild, from_idx=None, to_idx=None, custom_args=None,
    input_file=None, download_missing=False, download_threads=None,
):
    os.chdir(cwd)

    if genomes_dir:
        input_dirs = [str(genomes_dir)]
    else:
        organisms = DB_TYPE_CONFIG.get(db_type, [db_type])
        input_dirs = [f"{cache_dir}/refseq/{organism}" for organism in organisms]

    if rebuild:
        cmd = f"rm -f {db_name}.*"
        run_cmd(cmd)

    os.chdir(cwd)
    run_ganon_build_custom(
        cache_dir, input_dirs, db_name, threads, kmer_len, min_len, level,
        from_idx, to_idx, restart=rebuild, custom_args=custom_args,
        input_file=input_file, download_missing=download_missing,
        download_threads=download_threads,
    )


FILTER_DOWNLOAD_PASSES = {
    'cgrg': [{'assembly_levels': ['complete']}, {'refseq_categories': ['reference']}],
}


def ensure_filtered_genomes(cache_dir, organism, threads, filter_code, seed_dirs=None):
    # Union of assembly-level/refseq-category filters via two ncbi-genome-download passes into
    # the same dir, instead of ganon's own --genome-updater -F filter (which returns 0 files for
    # some groups because ganon ANDs --complete-genomes/--reference-genomes together)
    filtered_root = Path(cache_dir) / "refseq_filtered" / filter_code
    filtered_root.mkdir(parents=True, exist_ok=True)
    import_seed_genomes(filtered_root / "refseq" / organism, seed_dirs)
    def download():
        for kwargs in FILTER_DOWNLOAD_PASSES[filter_code]:
            logger.info(f"Downloading {organism} genomes ({filter_code}) with {kwargs}")
            ncbi_genome_download.download(
                section='refseq', groups=organism, file_formats='fasta',
                progress_bar=True, parallel=threads,
                output=str(filtered_root), uri=f"{NCBI_SERVER}/genomes", **kwargs
            )

    download()
    organism_dir = filtered_root / "refseq" / organism
    decompress_gzip_files(organism_dir, threads, retry_download=download)
    return organism_dir


def build_filtered_ganon(
    cache_dir, db_prefix, organism_groups, filter_code, ganon_args, seed_dirs=None,
    download_only=False, from_idx=None, to_idx=None, download_threads=None,
    build_threads=None,
):
    # cgrg-style combined filter: download via ncbi-genome-download instead of ganon's own
    # --genome-updater -F filter, then hand the files to build-custom
    logger.info(f"Filtered ({filter_code}) build for organism groups: {', '.join(organism_groups)} via ncbi-genome-download + ganon build-custom")
    ganon_threads = int((get_arg_values(ganon_args, '-t', '--threads') or ['4'])[0])
    build_threads = build_threads or ganon_threads
    download_threads = download_threads or build_threads
    kmer_len = (get_arg_values(ganon_args, '-k', '--kmer-size') or ['19'])[0]
    min_len = (get_arg_values(ganon_args, '-w', '--window-size') or ['31'])[0]
    level = (get_arg_values(ganon_args, '-l', '--level') or ['species'])[0]

    input_dirs = [
        str(
            ensure_filtered_genomes(
                cache_dir, organism, download_threads, filter_code, seed_dirs
            )
        )
        for organism in organism_groups
    ]
    if download_only:
        logger.info("Downloaded filtered genomes, skipping build (--download-only)")
        return
    run_ganon_build_custom(
        cache_dir, input_dirs, db_prefix, build_threads, kmer_len, min_len,
        level, from_idx, to_idx,
    )


def detect_filter_code(ganon_args):
    """cgrg (genome-updater) builds route through build_filtered_ganon; plain rg/cg builds
    keep using native `ganon build`, which already works fine for those on their own."""
    has_gu = get_option_value(ganon_args, '--genome-updater') is not None
    has_rg = any(arg in ('-r', '--reference-genomes') for arg in ganon_args)
    has_cg = any(arg in ('-c', '--complete-genomes') for arg in ganon_args)
    if has_gu or (has_rg and has_cg):
        return 'cgrg'
    return None


def build_combo_ganon(
    cache_dir, db_prefix, organism_groups, ganon_args, seed_dirs=None,
    download_only=False, from_idx=None, to_idx=None, download_threads=None,
    build_threads=None,
):
    # Combo db-prefix (e.g. archaea+viral): download each organism group into its own
    # persistent cache dir and reuse ganon build-custom, instead of native `ganon build`
    # which would re-download the whole combo as a single genome_updater job every time
    logger.info(f"Combo build for organism groups: {', '.join(organism_groups)} (downloaded/reused separately)")
    ganon_threads = int((get_arg_values(ganon_args, '-t', '--threads') or ['4'])[0])
    build_threads = build_threads or ganon_threads
    download_threads = download_threads or build_threads
    kmer_len = (get_arg_values(ganon_args, '-k', '--kmer-size') or ['19'])[0]
    min_len = (get_arg_values(ganon_args, '-w', '--window-size') or ['31'])[0]
    level = (get_arg_values(ganon_args, '-l', '--level') or ['species'])[0]

    for organism in organism_groups:
        ensure_organism_genomes(cache_dir, organism, download_threads, seed_dirs)

    if download_only:
        logger.info("Downloaded genomes, skipping build (--download-only)")
        return

    input_dirs = [f"{cache_dir}/refseq/{organism}" for organism in organism_groups]
    run_ganon_build_custom(
        cache_dir, input_dirs, db_prefix, build_threads, kmer_len, min_len,
        level, from_idx, to_idx,
    )


def get_files(genomes_dir, cache_dir, db_type, db_name, threads):
    if genomes_dir:
        logger.info(f"Adding {genomes_dir} genomes to library")

        decompress_gzip_files(genomes_dir, threads)

        cmd = f"find {genomes_dir} -name '*.gbff'"
        files = run_cmd(cmd, return_output=True)
        for file in files:
            if os.path.exists(f"{file}.fna"):
                continue
            cmd = f"any2fasta -u {file} > {file}.fna"
            run_cmd(cmd)

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
        limit, from_idx, to_idx, batch_size, threads, use_k2
):
    os.chdir(cwd)
    os.makedirs(cwd / db_name / "library", exist_ok=True)

    files = get_files(genomes_dir, cache_dir, db_type, db_name, threads)
    if from_idx is not None or to_idx is not None:
        logger.info(f"Using genome files range [{from_idx}:{to_idx}] out of {len(files)}")
        files = files[from_idx:to_idx]
    elif limit:
        logger.info(f"Limiting number of genomes to {limit}")
        files = files[:limit]
    if not files:
        raise click.ClickException(
            "No genome FASTA files are available to add to the database. "
            "Check the genome download and selected file range."
        )

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


ORGANISM_GROUPS = [
    'archaea', 'bacteria', 'fungi', 'human', 'invertebrate', 'metagenomes',
    'other', 'plant', 'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral',
]
LIBRARY_GROUP_TAXIDS = (
    ("human", "9606"),
    ("vertebrate_mammalian", "40674"),
    ("vertebrate_other", "7742"),
    ("invertebrate", "33208"),
    ("fungi", "4751"),
    ("plant", "33090"),
    ("protozoa", "5794"),
    ("viral", "10239"),
    ("archaea", "2157"),
    ("bacteria", "2"),
)
SOURCE_ALIASES = {'rs': 'refseq', 'gb': 'genbank'}
BOOL_FLAG_ALIASES = {'rg': 'reference-genomes', 'cg': 'complete-genomes'}
FILTER_CODE_EXPR = {'rg': '$5 == "reference genome"', 'cg': '$12 == "Complete Genome"'}
CATEGORY_ARG_ALIASES = {
    'source': ('-b', '--source'),
    'organism-group': ('-g', '--organism-group'),
    'reference-genomes': ('-r', '--reference-genomes'),
    'complete-genomes': ('-c', '--complete-genomes'),
    'genome-updater': ('--genome-updater',),
}


def get_arg_values(ganon_args, short, long):
    """Collect values passed to a repeatable ganon CLI option (e.g. -g archaea bacteria)."""
    values = []
    for index, arg in enumerate(ganon_args):
        if arg in (short, long):
            for v in ganon_args[index + 1:]:
                if v.startswith('-'):
                    break
                values.append(v)
        elif arg.startswith(f"{long}="):
            values.append(arg.split('=', 1)[1])
    return values


def get_option_value(ganon_args, *names):
    """Value of a single-value option, e.g. get_option_value(args, '--genome-updater')."""
    for index, arg in enumerate(ganon_args):
        if arg in names and index + 1 < len(ganon_args):
            return ganon_args[index + 1]
        for name in names:
            if arg.startswith(f"{name}="):
                return arg.split('=', 1)[1]
    return None


def genome_cache_key(ganon_args):
    """Stable cache key for the genome set a ganon build downloads, independent of --db-prefix."""
    source = get_arg_values(ganon_args, '-b', '--source') or ['refseq']
    organism_group = get_arg_values(ganon_args, '-g', '--organism-group')
    taxid = get_arg_values(ganon_args, '-a', '--taxid')
    flags = [f for flag_short, flag_long in (('-r', '--reference-genomes'), ('-c', '--complete-genomes'))
             for f in (flag_long.lstrip('-'),) if flag_short in ganon_args or flag_long in ganon_args]
    genome_updater_value = get_option_value(ganon_args, '--genome-updater')
    if genome_updater_value:
        flags.append('gu-' + hashlib.md5(genome_updater_value.encode()).hexdigest()[:8])
    parts = sorted(source) + sorted(organism_group or taxid) + sorted(flags) or ['default']
    return "_".join(parts)


def link_genome_cache(cache_dir, db_prefix, ganon_args):
    """Point ganon's {db_prefix}_files/ download dir at a shared, reusable genome cache."""
    genomes_cache_dir = Path(cache_dir) / "genomes" / genome_cache_key(ganon_args)
    genomes_cache_dir.mkdir(parents=True, exist_ok=True)
    files_link = Path.cwd() / f"{db_prefix}_files"
    if not files_link.exists():
        files_link.symlink_to(genomes_cache_dir)
        logger.info(f"Linked {files_link} -> shared genome cache {genomes_cache_dir}")
    elif not (files_link.is_symlink() and files_link.resolve() == genomes_cache_dir.resolve()):
        logger.info(f"{files_link} already exists, leaving as-is (not using shared genome cache)")


ORGANISM_GROUP_CODES = {}
for _group in ORGANISM_GROUPS:
    ORGANISM_GROUP_CODES.setdefault(''.join(w[0] for w in _group.split('_')), _group)
del _group


def word_break(segment, code_map):
    """Word-break a segment into codes from code_map, e.g. 'av' -> ['a', 'v'] for a 1-char code map."""
    max_code_len = max(len(code) for code in code_map)
    decoded = [[]] + [None] * len(segment)
    for end in range(1, len(segment) + 1):
        for length in range(min(max_code_len, end), 0, -1):
            start = end - length
            if decoded[start] is None:
                continue
            if segment[start:end] in code_map:
                decoded[end] = decoded[start] + [segment[start:end]]
                break
    return decoded[len(segment)]


def decode_organism_codes(segment):
    """Word-break a segment into organism-group codes, e.g. 'av' -> ['archaea', 'viral']."""
    codes = word_break(segment, ORGANISM_GROUP_CODES)
    return [ORGANISM_GROUP_CODES[code] for code in codes] if codes else None


def infer_ganon_options(db_prefix):
    """Infer --source/--organism-group/--reference-genomes/--complete-genomes from db-prefix segments (e.g. 'rs_a_rg')."""
    inferred = {}
    for segment in db_prefix.split('_'):
        if segment in SOURCE_ALIASES:
            inferred.setdefault('source', ['--source', SOURCE_ALIASES[segment]])
            continue
        if segment in BOOL_FLAG_ALIASES:
            category = BOOL_FLAG_ALIASES[segment]
            inferred.setdefault(category, [f'--{category}'])
            continue
        filter_codes = word_break(segment, FILTER_CODE_EXPR)
        if filter_codes and len(filter_codes) > 1:
            expr = ' || '.join(FILTER_CODE_EXPR[code] for code in filter_codes)
            inferred.setdefault('genome-updater', ['--genome-updater', f'-F {expr}'])
            continue
        groups = decode_organism_codes(segment)
        if groups:
            inferred.setdefault('organism-group', ['--organism-group'] + groups)
    return inferred


@click.group(no_args_is_help=True, epilog=f"Config file: {get_config_path()}")
@click.version_option(importlib.metadata.version('taxaforge'), '--version', '-v')
def cli():
    pass


@cli.command(no_args_is_help=True, context_settings={"ignore_unknown_options": True})
@click.option('--tool', default=lambda: load_config().get(CONFIG_SECTION, 'tool', fallback='kraken2'), type=click.Choice(list(REQUIRED_BINS)), help='Classifier to build the database for')
@click.option('--backend', default='direct', type=click.Choice(['direct', 'createtaxdb']), help='Build ganon directly or through nf-core/createtaxdb.')
@click.option('--db-type', default=None, help='database type to build')
@click.option('--db-name', default=None, help='database name to build')
@click.option('--genomes-dir', default=None, help='Directory containing genomes')
@click.option('--genomes-cache-dir', type=click.Path(file_okay=False, path_type=Path), default=None, help='Reuse requested FASTA files from this directory before downloading missing assemblies')
@click.option('--assembly-accessions-file', type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path), default=None, help='Download only NCBI assemblies listed in this file (accessions or FASTA filenames, one per line). Used only for ganon2')
@click.option('--nextflow-bin', default=lambda: os.environ.get('NEXTFLOW_BIN', '~/tools/nextflow/nextflow'), help='Nextflow executable for the createtaxdb backend.')
@click.option('--createtaxdb-pipeline', default=lambda: os.environ.get('CREATETAXDB_PIPELINE', 'nf-core/createtaxdb'), help='Nextflow pipeline name or path.')
@click.option('--pipeline-revision', default=lambda: os.environ.get('CREATETAXDB_REVISION'), help='Optional nf-core/createtaxdb pipeline revision.')
@click.option('--resume', is_flag=True, help='Resume a previous Nextflow run (createtaxdb backend only).')
@click.option('--input-file', default=None, help='Path to a text file listing genome filenames (one per line, matched against --genomes-dir) to build a ganon/ganon2 index from a specific subset, bypassing full directory scan')
@click.option('--download-missing', is_flag=True, help='With --input-file: auto-download any listed genomes missing from --genomes-dir via ncbi-genome-download (matched by GCA_/GCF_ accession in filename), then retry the build')
@click.option('--seed-dir', 'seed_dirs', multiple=True, default=lambda: tuple(d for d in load_config().get(CONFIG_SECTION, 'seed-dirs', fallback='').split(',') if d), help='Existing folder(s) with genomes to reuse before downloading; missing ones still get downloaded (repeatable, config key: seed-dirs)')
@click.option('--cache-dir', default=lambda: load_config().get(CONFIG_SECTION, 'cache-dir', fallback=str(create_cache_dir())), help='Cache directory for downloaded genomes/taxonomy (config key: cache-dir)')
@click.option('--output-dir', default=lambda: load_config().get(CONFIG_SECTION, 'output-dir', fallback='.'), help='Directory to build the database in, instead of cwd (config key: output-dir)')
@click.option('--threads', default=None, type=click.IntRange(min=1), help='Set both download and build threads (legacy alias).')
@click.option('--download-threads', default=None, type=click.IntRange(min=1, max=MAX_DOWNLOAD_THREADS), help=f'Number of parallel workers for downloading genomes (maximum {MAX_DOWNLOAD_THREADS}).')
@click.option('--build-threads', default=None, type=click.IntRange(min=1), help='Number of threads for adding genomes and building the index.')
@click.option('--load-factor', default=0.7, help='Proportion of the hash table to be populated. Used only for kraken2')
@click.option('--max-db-size', default=None, callback=parse_max_db_size, help='Maximum Kraken2 hash table size (for example, 50G).')
@click.option('--kmer-len', default=None, help='Override ganon k-mer size; otherwise use its default. Kraken2 defaults to 35.', type=int)
@click.option('--min-len', '--minimizer-len', 'min_len', default=None, help='Override ganon window size or Kraken2 minimizer length. Kraken2 defaults to 31.', type=int)
@click.option('--minimizer-spaces', default=None, type=click.IntRange(min=0), help='Kraken2 minimizer spaces (default 7, automatically reduced if the minimizer is too short).')
@click.option('--level', default=None, type=click.Choice(['leaves', 'species', 'genus', 'assembly', 'file']), help='Override ganon taxonomic level; otherwise use its default.')
@click.option('--limit', default=None, help='Limit number of genomes to use', type=int)
@click.option('--from', 'from_idx', default=None, help='Start index of genome files range, for testing a small slice first. Used only for kraken2', type=int)
@click.option('--to', 'to_idx', default=None, help='End index of genome files range, for testing a small slice first. Used only for kraken2', type=int)
@click.option('--batch-size', default=1000, help='Number of genomes to add to library at a time. Used only for kraken2', type=int)
@click.option('--force', is_flag=True, help='Force download and build')
@click.option('--rebuild', is_flag=True, help='Clean existing build files and re-build')
@click.option('--fast-build', is_flag=True, help='Non deterministic but faster build. Used only for kraken2')
@click.option('--use-k2', is_flag=True, help='Use k2 CLI instead of kraken2-build. Used only for kraken2')
@click.option('--download-only', is_flag=True, help='Download genomes/taxonomy only, skip the db build step')
@click.argument('ganon_args', nargs=-1, type=click.UNPROCESSED)
@click.pass_context
def build(
        context,
        tool: str, backend, db_type: str, db_name, cache_dir, output_dir, genomes_dir, genomes_cache_dir, assembly_accessions_file,
        nextflow_bin, createtaxdb_pipeline, pipeline_revision, resume, seed_dirs,
        threads, download_threads, build_threads, load_factor, max_db_size, kmer_len: int, min_len, minimizer_spaces, level: str, limit: int, from_idx: int, to_idx: int, batch_size: int,
        input_file, download_missing,
        force: bool, rebuild, fast_build: bool, use_k2: bool, download_only: bool, ganon_args
):
    if backend == "createtaxdb" and tool != "ganon2":
        raise click.UsageError("--backend createtaxdb requires --tool ganon2")
    if backend == "createtaxdb" and not assembly_accessions_file:
        raise click.UsageError(
            "--backend createtaxdb requires --assembly-accessions-file"
        )
    if backend == "createtaxdb" and genomes_dir:
        raise click.UsageError(
            "--backend createtaxdb uses --assembly-accessions-file, not --genomes-dir"
        )
    if resume and backend != "createtaxdb":
        raise click.UsageError("--resume is available only with --backend createtaxdb")
    if resume and (force or rebuild):
        raise click.UsageError("--resume cannot be combined with --force or --rebuild")
    if backend == "createtaxdb" and (force or rebuild):
        raise click.UsageError(
            "Use a fresh run without --resume to rebuild with the createtaxdb backend"
        )
    if backend == "createtaxdb" and ganon_args:
        raise click.UsageError(
            "Pass ganon options with TaxaForge options such as --kmer-len and --min-len"
        )
    if assembly_accessions_file and genomes_dir:
        raise click.UsageError("--assembly-accessions-file cannot be combined with --genomes-dir")
    if assembly_accessions_file and tool != "ganon2":
        raise click.UsageError("--assembly-accessions-file is currently supported only with --tool ganon2")
    if input_file and not genomes_dir:
        raise click.UsageError("--input-file requires --genomes-dir")
    if download_missing and not input_file:
        raise click.UsageError("--download-missing requires --input-file")
    if input_file and tool not in ("ganon", "ganon2"):
        raise click.UsageError("--input-file is supported only with --tool ganon or --tool ganon2")
    if input_file and backend != "direct":
        raise click.UsageError("--input-file cannot be used with --backend createtaxdb")
    if max_db_size and tool != "kraken2":
        raise click.UsageError("--max-db-size is only available with --tool kraken2")
    if max_db_size and use_k2:
        raise click.UsageError("--max-db-size cannot be combined with --use-k2")
    if assembly_accessions_file:
        try:
            read_assembly_accessions(assembly_accessions_file, limit)
        except ValueError as error:
            raise click.BadParameter(str(error), param_hint="--assembly-accessions-file") from error

    default_threads = max(1, int(multiprocessing.cpu_count() * 0.8))
    download_threads = min(
        download_threads or threads or default_threads, MAX_DOWNLOAD_THREADS
    )
    build_threads = build_threads or threads or default_threads
    threads = build_threads

    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(output_dir)
    logger.info(f"Building in {output_dir}")

    if kmer_len is None:
        # ganon's own default k-mer size is 19; kraken2-build's is 35. Use whichever
        # matches the selected tool instead of forcing one default on both.
        kmer_len = 19 if tool in ('ganon', 'ganon2') else 35

    if tool == 'kraken2':
        min_len = 31 if min_len is None else min_len
        if min_len > kmer_len:
            logger.warning(
                f"--minimizer-len {min_len} exceeds --kmer-len {kmer_len}; "
                f"reducing --minimizer-len to {kmer_len} for Kraken2"
            )
            min_len = kmer_len
        max_minimizer_spaces = min_len // 4
        if minimizer_spaces is None:
            minimizer_spaces = min(
                KRAKEN2_DEFAULT_MINIMIZER_SPACES, max_minimizer_spaces
            )
            if minimizer_spaces < KRAKEN2_DEFAULT_MINIMIZER_SPACES:
                logger.warning(
                    f"Kraken2 allows at most {max_minimizer_spaces} minimizer "
                    f"spaces for --minimizer-len {min_len}; reducing the default "
                    f"from {KRAKEN2_DEFAULT_MINIMIZER_SPACES}"
                )
        elif minimizer_spaces > max_minimizer_spaces:
            raise click.UsageError(
                f"--minimizer-spaces {minimizer_spaces} exceeds the maximum "
                f"{max_minimizer_spaces} for --minimizer-len {min_len}"
            )

    if tool in ('ganon', 'ganon2') and ganon_args and not assembly_accessions_file:
        logger.info(f"Passing through unrecognized options to native ganon build: {' '.join(ganon_args)}")
        run_basic_checks(tool, use_k2)
        ganon_args = list(ganon_args)

        db_prefix_value = None
        for index, arg in enumerate(ganon_args):
            if arg in ('-d', '--db-prefix') and index + 1 < len(ganon_args):
                db_prefix_value = ganon_args[index + 1]
            elif arg.startswith('--db-prefix='):
                db_prefix_value = arg.split('=', 1)[1]

        if not db_prefix_value and db_name:
            db_prefix_value = db_name
            ganon_args += ['-d', db_name]

        if db_prefix_value and genomes_dir:
            if download_only:
                logger.info(f"Genomes dir {genomes_dir} provided, nothing to download (--download-only)")
                return
            run_ganon_build_custom(
                cache_dir, [str(genomes_dir)], db_prefix_value, threads,
                kmer_len, min_len, level, from_idx, to_idx,
                custom_args=list(ganon_args),
                input_file=input_file, download_missing=download_missing,
                download_threads=download_threads,
            )
            return

        if db_prefix_value:
            for category, tokens in infer_ganon_options(db_prefix_value).items():
                aliases = CATEGORY_ARG_ALIASES[category]
                already_set = any(arg in aliases or arg.startswith(f"{aliases[-1]}=") for arg in ganon_args)
                if not already_set:
                    logger.info(f"Inferred {' '.join(tokens)} from db-prefix '{db_prefix_value}'")
                    ganon_args += tokens

        organism_groups = get_arg_values(ganon_args, '-g', '--organism-group')
        filter_code = detect_filter_code(ganon_args)
        if db_prefix_value and organism_groups and filter_code:
            build_filtered_ganon(
                cache_dir, db_prefix_value, organism_groups, filter_code,
                ganon_args, seed_dirs, download_only, from_idx, to_idx,
                download_threads, build_threads,
            )
            return
        if db_prefix_value and len(organism_groups) > 1:
            build_combo_ganon(
                cache_dir, db_prefix_value, organism_groups, ganon_args,
                seed_dirs, download_only, from_idx, to_idx, download_threads,
                build_threads,
            )
            return

        if download_only:
            logger.warning("Native `ganon build` downloads and builds in one step, can't split; skipping (--download-only)")
            return

        taxdump = Path(cache_dir) / "taxonomy" / "taxdump.tar.gz"
        if taxdump.exists() and not any(arg in ('-m', '--taxonomy-files') for arg in ganon_args):
            logger.info(f"Reusing cached taxonomy files at {taxdump}")
            ganon_args += ["-m", str(taxdump)]
        if db_prefix_value and '--restart' not in ganon_args:
            link_genome_cache(cache_dir, db_prefix_value, ganon_args)
        cmd = "ganon build " + " ".join(shlex.quote(arg) for arg in ganon_args)
        run_cmd(cmd)
        return

    logger.info(f"Building {tool} database of type {db_type}")
    if backend == "direct":
        run_basic_checks(tool, use_k2)
    cwd = output_dir

    if cache_dir == '.':
        cache_dir = cwd

    if not db_name:
        db_name = f"{tool}_{context.params['db_type']}"

    if force:
        if tool == 'kraken2':
            run_cmd(f"rm -rf {db_name}")
            run_cmd(f"mkdir -p {db_name}")
        else:
            run_cmd(f"rm -f {db_name}.*")

    logger.info(f"Using cache directory {cache_dir}")

    download_taxanomy(cache_dir, skip_maps=(tool == 'ganon2'))

    if assembly_accessions_file:
        try:
            genomes_dir = download_assembly_accessions(
                cache_dir, assembly_accessions_file, download_threads, limit, genomes_cache_dir
            )
        except ValueError as error:
            raise click.BadParameter(str(error), param_hint="--assembly-accessions-file") from error
        except RuntimeError as error:
            raise click.ClickException(str(error)) from error
    elif not genomes_dir:
        download_genomes(
            cache_dir, cwd, db_type, db_name, download_threads, force, seed_dirs
        )

    if download_only:
        logger.info("Downloaded genomes/taxonomy, skipping build (--download-only)")
        return

    if tool == 'kraken2':
        kmer_len = 35 if kmer_len is None else kmer_len
        add_to_library(
            cache_dir, cwd, genomes_dir, db_type, db_name,
            limit, from_idx, to_idx, batch_size, threads, use_k2
        )
        build_db(
            cache_dir, cwd, db_type, db_name, threads, kmer_len, min_len,
            fast_build, rebuild, load_factor, use_k2, minimizer_spaces,
            max_db_size=max_db_size,
        )
    elif tool == 'ganon2':
        if backend == "createtaxdb":
            try:
                run_createtaxdb(
                    nextflow_bin, createtaxdb_pipeline, pipeline_revision,
                    cwd, db_name, cache_dir, assembly_accessions_file,
                    genomes_dir, limit, threads, kmer_len, min_len, level, resume,
                )
            except ValueError as error:
                raise click.ClickException(str(error)) from error
        else:
            build_ganon2(
                cache_dir, cwd, genomes_dir, db_type, db_name, threads,
                kmer_len, min_len, level, rebuild, from_idx, to_idx,
                custom_args=list(ganon_args) or None,
                input_file=input_file, download_missing=download_missing,
                download_threads=download_threads,
            )


@cli.group(name='config', invoke_without_command=True)
@click.pass_context
def config(context):
    if context.invoked_subcommand is not None:
        return

    parser = load_config()
    items = parser.items(CONFIG_SECTION)
    if not items:
        logger.info(f"No config set. Config file: {get_config_path()}")
        return

    for key, value in items:
        print(f"{key} = {value}")


@config.command(name='get')
@click.argument('key')
def config_get(key):
    parser = load_config()
    if not parser.has_option(CONFIG_SECTION, key):
        logger.error(f"{key} not set")
        sys.exit(1)
    print(parser.get(CONFIG_SECTION, key))


@config.command(name='set')
@click.argument('key')
@click.argument('value')
def config_set(key, value):
    parser = load_config()
    parser.set(CONFIG_SECTION, key, value)
    save_config(parser)
    logger.info(f"Set {key} = {value}")


@cli.command(name="download-taxonomy")
@click.option(
    "--cache-dir",
    default=lambda: load_config().get(
        CONFIG_SECTION, "cache-dir", fallback=str(create_cache_dir())
    ),
    show_default="configured cache directory",
    type=click.Path(file_okay=False, path_type=Path),
    help="TaxaForge cache directory containing taxonomy data.",
)
def download_taxonomy(cache_dir):
    """Ensure Kraken taxonomy and nucleotide accession maps are downloaded."""
    download_taxanomy(cache_dir)
    logger.info(f"Taxonomy is ready in {cache_dir / 'taxonomy'}")


@cli.command(name='createtaxdb-input')
@click.option('--input-file', required=True, help='Path to a text file listing genome filenames/accessions (one per line) to include in the samplesheet')
@click.option('--genomes-dir', default=None, help='Directory containing the genome fasta files to resolve --input-file entries against. Not required with --remote')
@click.option('--output', default='samplesheet.csv', help='Path to write the nf-core/createtaxdb samplesheet CSV to (default: samplesheet.csv)')
@click.option('--remote', is_flag=True, help="Don't require genomes on local disk: use each accession's NCBI HTTPS ftp_path (from assembly_summary) as fasta_dna instead, since nf-core/createtaxdb/Nextflow can download those URLs itself. Skips --genomes-dir and --download-missing entirely")
@click.option('--download-missing', is_flag=True, help='With local (non --remote) mode: auto-download any listed genomes missing from --genomes-dir via ncbi-genome-download (matched by GCA_/GCF_ accession in filename), then retry resolution')
@click.option('--taxonomy-source', type=click.Choice(['auto', 'refseq', 'genbank', 'both']), default='auto', help='Which NCBI assembly_summary file(s) to download for accession->taxid/ftp_path lookup. "auto" (default) inspects --input-file and only fetches genbank\'s (~1.75GB) if GCA_ accessions are present, skipping it for refseq-only (GCF_) lists')
@click.option('--cache-dir', default=lambda: load_config().get(CONFIG_SECTION, 'cache-dir', fallback=str(create_cache_dir())), help='Cache directory for downloaded genomes/taxonomy (config key: cache-dir)')
@click.option('--threads', default=max(1, min(MAX_DOWNLOAD_THREADS, int(multiprocessing.cpu_count() * 0.8))), help=f'Number of download workers for --download-missing (maximum {MAX_DOWNLOAD_THREADS})', type=click.IntRange(min=1, max=MAX_DOWNLOAD_THREADS))
def createtaxdb_input(input_file, genomes_dir, output, remote, download_missing, taxonomy_source, cache_dir, threads):
    """Turn a plain genome filename list into an nf-core/createtaxdb samplesheet (id,taxid,fasta_dna)."""
    if not remote and not genomes_dir:
        logger.error("--genomes-dir is required unless --remote is set")
        sys.exit(1)

    if remote:
        names = read_genome_name_list(input_file)
        missing = []
    else:
        input_dirs = [genomes_dir]
        input_extension = detect_input_extension(input_dirs)
        files, missing = resolve_input_file_list(input_file, input_dirs, input_extension)
        if missing and download_missing:
            download_missing_genomes(missing, genomes_dir, threads)
            files, missing = resolve_input_file_list(input_file, input_dirs, input_extension)

    scan_names = names if remote else [Path(f).name for f in files] + missing
    if taxonomy_source == 'auto':
        sections = tuple(sorted({
            "refseq" if name.startswith("GCF_") else "genbank"
            for name in scan_names if ACCESSION_RE.match(name)
        })) or ("refseq", "genbank")
    elif taxonomy_source == 'both':
        sections = ("refseq", "genbank")
    else:
        sections = (taxonomy_source,)
    logger.info(f"Fetching assembly_summary for section(s): {', '.join(sections)}")
    assembly_summary_paths = ensure_assembly_summaries(cache_dir, sections=sections)

    rows, unmatched = [], 0
    if remote:
        accession_to_info = build_accession_metadata_map(assembly_summary_paths)
        for name in names:
            match = ACCESSION_RE.match(Path(name).name)
            accession = match.group(1) if match else None
            info = accession_to_info.get(accession) if accession else None
            if info:
                taxid, ftp_path = info
                rows.append((accession, taxid, ftp_path_to_fasta_url(ftp_path)))
            else:
                unmatched += 1
    else:
        accession_to_taxid = build_accession_taxid_map(assembly_summary_paths)
        for f in files:
            match = ACCESSION_RE.match(Path(f).name)
            accession = match.group(1) if match else None
            taxid = accession_to_taxid.get(accession) if accession else None
            if taxid:
                rows.append((accession, taxid, str(Path(f).resolve())))
            else:
                unmatched += 1

    with open(output, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["id", "taxid", "fasta_dna"])
        writer.writerows(rows)

    total = len(names) if remote else len(files)
    if unmatched:
        logger.warning(f"{unmatched} of {total} entries had no accession/taxid match; they were skipped")
    if missing and not remote and not download_missing:
        logger.warning(f"{len(missing)} entries from {input_file} could not be found under {genomes_dir}; re-run with --download-missing to fetch them locally, or use --remote to reference NCBI URLs instead")
    logger.info(f"Wrote {len(rows)} entries to {output}")


DOCTOR_BINS = ["ncbi-genome-download", "kraken2-build", "k2", "ganon", "any2fasta", "wget", "tar", "gunzip"]


def format_size(num_bytes):
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if num_bytes < 1024:
            return f"{num_bytes:.1f}{unit}"
        num_bytes /= 1024
    return f"{num_bytes:.1f}PB"


def genome_cache_status(cache_dir):
    """Count downloaded genome files and their size in the latest version of each cached organism group."""
    status = {}
    genomes_dir = Path(cache_dir) / "genomes"
    if not genomes_dir.is_dir():
        return status
    for group_dir in sorted(genomes_dir.iterdir()):
        if not group_dir.is_dir():
            continue
        versions = sorted(d for d in group_dir.iterdir() if d.is_dir())
        latest = versions[-1] if versions else None
        files_dir = latest / "files" if latest else None
        sizes = [f.stat().st_size for f in files_dir.rglob('*') if f.is_file()] if files_dir and files_dir.is_dir() else []
        total_size = sum(sizes)
        avg_size = total_size / len(sizes) if sizes else 0
        status[group_dir.name] = (len(sizes), total_size, avg_size)
    return status


FASTA_FILE_SUFFIXES = (
    ".fa", ".fasta", ".fna", ".fas", ".ffn", ".fsa",
    ".fa.gz", ".fasta.gz", ".fna.gz", ".fas.gz", ".ffn.gz", ".fsa.gz",
    ".zip",
)
FASTA_REPORT_SUFFIXES = FASTA_FILE_SUFFIXES[:-1]
GENOME_GROUP_ALIASES = {"viruses": "viral"}


def _genome_group_from_path(path_parts):
    for part in reversed(path_parts):
        group = part.casefold().replace("-", "_")
        group = GENOME_GROUP_ALIASES.get(group, group)
        if group in ORGANISM_GROUPS:
            return group
    return "unclassified"


def analyze_genome_directory(
    input_dir, progress_callback=None, progress_interval=1000
):
    """Count FASTA files by group inferred from group-named directories."""
    if progress_interval < 1:
        raise ValueError("progress_interval must be at least 1")

    root = Path(input_dir).expanduser().resolve()
    counts = {}
    scanned = 0
    for current_dir, _, filenames in os.walk(root):
        relative_parts = Path(current_dir).relative_to(root).parts
        group = _genome_group_from_path(root.parts + relative_parts)
        for filename in filenames:
            if filename.casefold().endswith(FASTA_FILE_SUFFIXES):
                counts[group] = counts.get(group, 0) + 1
            scanned += 1
            if progress_callback and scanned % progress_interval == 0:
                progress_callback(scanned, counts.copy())
    if progress_callback and scanned % progress_interval:
        progress_callback(scanned, counts.copy())
    return counts


def fasta_report_files(input_dir):
    root = Path(input_dir).expanduser().resolve()
    for current_dir, directory_names, filenames in os.walk(root):
        directory_names.sort()
        for filename in sorted(filenames):
            if filename.casefold().endswith(FASTA_REPORT_SUFFIXES):
                yield Path(current_dir) / filename


def fasta_headers(fasta_path, limit=None):
    opener = gzip.open if fasta_path.name.casefold().endswith(".gz") else open
    headers = []
    saw_content = False
    with opener(fasta_path, "rt", encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            content = line.strip()
            if not content:
                continue
            saw_content = True
            if content.startswith(">"):
                header = content[1:].strip()
                if not header:
                    raise ValueError(
                        f"Empty FASTA header on line {line_number} of {fasta_path}"
                    )
                headers.append(header)
                if limit is not None and len(headers) >= limit:
                    return headers
            elif not headers:
                raise ValueError(
                    f"Expected a FASTA header on line {line_number} of {fasta_path}"
                )
    if saw_content and not headers:
        raise ValueError(f"No FASTA headers found in {fasta_path}")
    return headers


def library_report_rows_for_batch(fasta_paths, root, one_per_file=False):
    rows = []
    for fasta_path in fasta_paths:
        group = _genome_group_from_path(
            fasta_path.relative_to(root).parts[:-1]
        )
        headers = fasta_headers(
            fasta_path, limit=1 if one_per_file else None
        )
        rows.extend(
            (group, header, fasta_path.as_uri())
            for header in headers
        )
    return rows


def library_group_from_taxid(taxid, nodes, cache):
    if taxid in cache:
        return cache[taxid]

    lineage = []
    seen = set()
    current = taxid
    while current not in seen:
        seen.add(current)
        lineage.append(current)
        node = nodes.get(current)
        if node is None or node[0] == current:
            break
        current = node[0]

    lineage_taxids = set(lineage)
    group = next(
        (
            group
            for group, ancestor_taxid in LIBRARY_GROUP_TAXIDS
            if ancestor_taxid in lineage_taxids
        ),
        "unclassified",
    )
    for descendant in lineage:
        cache[descendant] = group
    return group


def library_assembly_accession(url):
    match = re.search(r"(GC[AF]_\d+\.\d+)", url)
    return match.group(1) if match else None


def resolve_library_assembly_accessions(rows, taxonomy_dir):
    accessions = {
        accession
        for row in rows
        if (accession := library_assembly_accession(row[2]))
    }
    summary_paths = sorted(taxonomy_dir.glob("assembly_summary_*.txt"))
    if not accessions or not summary_paths:
        return {}
    return build_accession_taxid_map(summary_paths, accessions)


def classify_library_report_rows(rows, taxonomy_dir, workers=1):
    """Classify unclassified report rows from their sequence accessions."""
    accession_to_taxid = resolve_library_accessions(
        rows, taxonomy_dir, workers
    )
    unresolved_rows = [
        row
        for row in rows
        if row[0] == "unclassified"
        and library_record_accession(row[1]) not in accession_to_taxid
    ]
    assembly_accession_to_taxid = resolve_library_assembly_accessions(
        unresolved_rows, taxonomy_dir
    )
    nodes = parse_taxonomy_nodes(taxonomy_dir)
    cache = {}
    classified = []
    unresolved = 0

    for library, sequence_name, url in rows:
        if library != "unclassified":
            classified.append((library, sequence_name, url))
            continue

        accession = library_record_accession(sequence_name)
        taxid = accession_to_taxid.get(accession) or (
            assembly_accession_to_taxid.get(
                library_assembly_accession(url)
            )
        )
        group = (
            library_group_from_taxid(taxid, nodes, cache)
            if taxid is not None
            else "unclassified"
        )
        if group == "unclassified":
            unresolved += 1
        classified.append((group, sequence_name, url))

    return classified, unresolved


def fasta_report_file_batches(input_dir, limit=None):
    files = fasta_report_files(input_dir)
    if limit is not None:
        files = itertools.islice(files, limit)
    files = iter(files)
    while batch := list(itertools.islice(files, LIBRARY_REPORT_BATCH_SIZE)):
        yield batch


def iter_library_report_rows(
    input_dir, workers=1, limit=None, one_per_file=False
):
    """Yield library-report rows without retaining the entire report in memory."""
    root = Path(input_dir).expanduser().resolve()
    batches = iter(fasta_report_file_batches(root, limit))
    if workers == 1:
        for batch in batches:
            yield from library_report_rows_for_batch(
                batch, root, one_per_file
            )
        return

    max_pending_batches = workers * 4
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        pending = deque()

        def submit_batches():
            for batch in itertools.islice(
                batches, max_pending_batches - len(pending)
            ):
                pending.append(
                    executor.submit(
                        library_report_rows_for_batch,
                        batch,
                        root,
                        one_per_file,
                    )
                )

        submit_batches()
        while pending:
            yield from pending.popleft().result()
            submit_batches()


def create_library_report(
    input_dir, workers=1, limit=None, one_per_file=False
):
    """Return library-report rows for all FASTA records under input_dir."""
    return list(
        iter_library_report_rows(input_dir, workers, limit, one_per_file)
    )


def extract_kraken_species(inspect_file, include_s1=False):
    """Return distinct species taxids and names from Kraken inspect output."""
    accepted_ranks = {"S", "S1"} if include_s1 else {"S"}
    taxa = {}

    with open(inspect_file, encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t", 5)
            if len(fields) != 6:
                raise ValueError(
                    f"Expected 6 tab-separated fields on line {line_number}"
                )

            rank = fields[3]
            taxid = fields[4].strip()
            if not taxid.isdigit():
                raise ValueError(
                    f"Invalid taxid {taxid!r} on line {line_number}"
                )
            if rank in accepted_ranks:
                taxa[taxid] = (rank, fields[5].strip())

    return [
        (taxid, rank, name)
        for taxid, (rank, name) in sorted(taxa.items(), key=lambda item: int(item[0]))
    ]


def normalize_taxon_name(name):
    return " ".join(name.casefold().replace("[", "").replace("]", "").split())


def library_species_prefixes(sequence_name, max_words=12):
    """Yield normalized leading organism-name candidates from a library label."""
    parts = sequence_name.strip().split(maxsplit=1)
    if len(parts) < 2:
        return

    description = parts[1]
    if description.startswith("MAG: "):
        description = description[5:]
    words = description.split()
    for word_count in range(1, min(max_words, len(words)) + 1):
        yield normalize_taxon_name(
            " ".join(words[:word_count]).rstrip(",.;")
        )


def read_library_report(library_report):
    with open(library_report, newline="", encoding="utf-8") as source:
        sample = source.read(4096)
        source.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        except csv.Error:
            dialect = csv.excel_tab
        reader = csv.reader(source, dialect)
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("Library report is empty")
        if header != ["#Library", "Sequence Name", "URL"]:
            raise ValueError(
                "Expected header '#Library\\tSequence Name\\tURL'"
            )

        rows = []
        for line_number, row in enumerate(reader, start=2):
            if len(row) != 3:
                raise ValueError(
                    f"Expected 3 tab-separated fields on line {line_number}"
                )
            rows.append(row)
    return rows


def count_library_species_rows(rows, taxa):
    counts = {}
    unmatched = 0

    for row in rows:
        match = None
        for candidate in library_species_prefixes(row[1]):
            if candidate in taxa:
                match = taxa[candidate]
        if match is None:
            unmatched += 1
            continue

        counts[match] = counts.get(match, 0) + 1

    return counts, unmatched


library_taxa = None


def initialize_library_taxa(taxa):
    global library_taxa
    library_taxa = taxa


def count_library_species_chunk(rows):
    return count_library_species_rows(rows, library_taxa)


def combine_library_species_counts(results):
    counts = {}
    unmatched = 0
    for chunk_counts, chunk_unmatched in results:
        unmatched += chunk_unmatched
        for taxon, count in chunk_counts.items():
            counts[taxon] = counts.get(taxon, 0) + count
    return counts, unmatched


def parallel_library_species_counts(rows, taxa, workers):
    if workers == 1 or len(rows) < EXTRACTION_CHUNK_SIZE:
        return count_library_species_rows(rows, taxa)

    chunks = [
        rows[index : index + EXTRACTION_CHUNK_SIZE]
        for index in range(0, len(rows), EXTRACTION_CHUNK_SIZE)
    ]
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=min(workers, len(chunks)),
        initializer=initialize_library_taxa,
        initargs=(taxa,),
    ) as executor:
        return combine_library_species_counts(
            executor.map(count_library_species_chunk, chunks)
        )


def extract_library_species(
    library_report, inspect_file, include_s1=False, workers=1
):
    """Resolve library records to distinct Kraken species taxids by name."""
    taxa = {
        normalize_taxon_name(name): (taxid, rank, name)
        for taxid, rank, name in extract_kraken_species(inspect_file, include_s1)
    }
    counts, unmatched = parallel_library_species_counts(
        read_library_report(library_report), taxa, workers
    )

    rows = [
        (taxid, rank, name, count)
        for (taxid, rank, name), count in counts.items()
    ]
    rows.sort(key=lambda row: int(row[0]))
    return rows, unmatched


def library_record_accession(sequence_name):
    parts = sequence_name.strip().split(maxsplit=1)
    return parts[0] if len(parts) == 2 else None


def parse_taxonomy_nodes(taxonomy_dir):
    nodes_path = taxonomy_dir / "nodes.dmp"
    if not nodes_path.is_file():
        raise ValueError(f"Missing taxonomy nodes file: {nodes_path}")

    nodes = {}
    with open(nodes_path, encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            fields = line.rstrip("\t|\n").split("\t|\t")
            if len(fields) < 3 or not fields[0].isdigit() or not fields[1].isdigit():
                raise ValueError(
                    f"Invalid taxonomy node on line {line_number} of {nodes_path}"
                )
            nodes[fields[0]] = (fields[1], fields[2])
    return nodes


def parse_taxonomy_names(taxonomy_dir, taxids):
    names_path = taxonomy_dir / "names.dmp"
    if not names_path.is_file():
        raise ValueError(f"Missing taxonomy names file: {names_path}")

    names = {}
    with open(names_path, encoding="utf-8") as source:
        for line in source:
            fields = line.rstrip("\t|\n").split("\t|\t")
            if (
                len(fields) >= 4
                and fields[0] in taxids
                and fields[3] == "scientific name"
            ):
                names[fields[0]] = fields[1]
    return names


def resolve_accession_taxids_range(
    map_path, start, end, wanted, base_accessions
):
    accession_to_taxid = {}
    with open(map_path, "rb") as source:
        source.seek(start)
        if start:
            source.seek(start - 1)
            if source.read(1) != b"\n":
                source.readline()
            else:
                source.seek(start)
        while end is None or source.tell() < end:
            line = source.readline()
            if not line:
                break
            fields = line.rstrip(b"\n").split(b"\t")
            if len(fields) < 3 or fields[0] == b"accession":
                continue
            accession, accession_version, taxid = fields[:3]
            if not taxid.isdigit():
                continue
            if accession_version in wanted:
                accession_to_taxid[wanted[accession_version]] = taxid.decode(
                    "ascii"
                )
            elif accession in base_accessions:
                for library_accession in base_accessions[accession]:
                    accession_to_taxid[library_accession] = taxid.decode("ascii")
    return accession_to_taxid


accession_wanted = None
accession_base_accessions = None


def initialize_accession_maps(wanted, base_accessions):
    global accession_wanted
    global accession_base_accessions
    accession_wanted = wanted
    accession_base_accessions = base_accessions


def resolve_accession_taxids_chunk(task):
    map_path, start, end = task
    return resolve_accession_taxids_range(
        map_path, start, end, accession_wanted, accession_base_accessions
    )


def accession_map_ranges(map_path, parts):
    file_size = map_path.stat().st_size
    if parts == 1 or file_size == 0:
        return [(map_path, 0, None)]

    chunk_size = (file_size + parts - 1) // parts
    return [
        (map_path, start, min(start + chunk_size, file_size))
        for start in range(0, file_size, chunk_size)
    ]


def accession_map_workers(map_paths, workers):
    sizes = {map_path: map_path.stat().st_size for map_path in map_paths}
    total_size = sum(sizes.values())
    if workers == 1 or total_size < ACCESSION_MAP_PARALLEL_MIN_BYTES:
        return {map_path: 1 for map_path in map_paths}

    worker_counts = {map_path: 1 for map_path in map_paths}
    for _ in range(workers - len(map_paths)):
        map_path = max(
            map_paths,
            key=lambda path: sizes[path] / worker_counts[path],
        )
        worker_counts[map_path] += 1
    return worker_counts


def resolve_library_accessions(rows, taxonomy_dir, workers=1):
    wanted = {
        accession
        for row in rows
        if (accession := library_record_accession(row[1]))
    }
    base_accessions = {}
    for accession in wanted:
        base_accessions.setdefault(accession.rsplit(".", 1)[0], []).append(
            accession
        )
    encoded_wanted = {
        accession.encode("utf-8"): accession for accession in wanted
    }
    encoded_base_accessions = {
        accession.encode("utf-8"): library_accessions
        for accession, library_accessions in base_accessions.items()
    }

    map_paths = [
        taxonomy_dir / "nucl_gb.accession2taxid",
        taxonomy_dir / "nucl_wgs.accession2taxid",
    ]
    missing_paths = [path for path in map_paths if not path.is_file()]
    if missing_paths:
        paths = ", ".join(str(path) for path in missing_paths)
        raise ValueError(
            "Missing nucleotide accession-to-taxid map(s): "
            f"{paths}. Run a TaxaForge build without --skip-maps first."
        )

    worker_counts = accession_map_workers(map_paths, workers)
    tasks = [
        task
        for map_path in map_paths
        for task in accession_map_ranges(
            map_path, worker_counts[map_path]
        )
    ]
    if len(tasks) == 1:
        map_results = [
            resolve_accession_taxids_range(
                *tasks[0], encoded_wanted, encoded_base_accessions
            )
        ]
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=min(workers, len(tasks)),
            initializer=initialize_accession_maps,
            initargs=(encoded_wanted, encoded_base_accessions),
        ) as executor:
            map_results = list(executor.map(resolve_accession_taxids_chunk, tasks))

    accession_to_taxid = {}
    for map_path in map_paths:
        for task, map_result in zip(tasks, map_results):
            if task[0] == map_path:
                accession_to_taxid.update(map_result)
    return accession_to_taxid


def species_ancestor(taxid, nodes, cache):
    if taxid in cache:
        return cache[taxid]

    lineage = []
    seen = set()
    current = taxid
    while current not in seen:
        seen.add(current)
        lineage.append(current)
        node = nodes.get(current)
        if node is None:
            break
        parent, rank = node
        if rank == "species":
            for descendant in lineage:
                cache[descendant] = current
            return current
        if parent == current:
            break
        current = parent

    for descendant in lineage:
        cache[descendant] = None
    return None


taxonomy_nodes = None


def initialize_taxonomy_nodes(nodes):
    global taxonomy_nodes
    taxonomy_nodes = nodes


def resolve_species_taxids(taxids):
    cache = {}
    return {
        taxid: species_ancestor(taxid, taxonomy_nodes, cache)
        for taxid in taxids
    }


def parallel_species_ancestors(taxids, nodes, workers):
    if workers == 1 or len(taxids) < EXTRACTION_CHUNK_SIZE:
        cache = {}
        return {
            taxid: species_ancestor(taxid, nodes, cache) for taxid in taxids
        }

    chunks = [
        taxids[index : index + EXTRACTION_CHUNK_SIZE]
        for index in range(0, len(taxids), EXTRACTION_CHUNK_SIZE)
    ]
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=min(workers, len(chunks)),
        initializer=initialize_taxonomy_nodes,
        initargs=(nodes,),
    ) as executor:
        resolved_chunks = executor.map(resolve_species_taxids, chunks)
        return {
            taxid: species_taxid
            for resolved_chunk in resolved_chunks
            for taxid, species_taxid in resolved_chunk.items()
        }


def extract_library_species_from_taxonomy(
    library_report, taxonomy_dir, workers=1
):
    """Resolve library accessions to species with TaxaForge taxonomy maps."""
    rows = read_library_report(library_report)
    accession_to_taxid = resolve_library_accessions(
        rows, taxonomy_dir, workers
    )
    nodes = parse_taxonomy_nodes(taxonomy_dir)
    taxids = sorted(
        {
            taxid
            for row in rows
            if (accession := library_record_accession(row[1]))
            if (taxid := accession_to_taxid.get(accession))
        },
        key=int,
    )
    species_taxids = parallel_species_ancestors(taxids, nodes, workers)
    resolved = []
    unmatched = 0

    for row in rows:
        accession = library_record_accession(row[1])
        taxid = accession_to_taxid.get(accession)
        species_taxid = species_taxids.get(taxid)
        if species_taxid is None:
            unmatched += 1
            continue
        resolved.append(species_taxid)

    names = parse_taxonomy_names(taxonomy_dir, set(resolved))
    counts = {}
    for taxid in resolved:
        if taxid in names:
            counts[taxid] = counts.get(taxid, 0) + 1
        else:
            unmatched += 1

    return [
        (taxid, "species", names[taxid], count)
        for taxid, count in sorted(counts.items(), key=lambda item: int(item[0]))
    ], unmatched


@cli.command(name="extract-kraken-species")
@click.argument(
    "inspect_file",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
)
@click.option(
    "--output",
    default="unique-species.tsv",
    show_default=True,
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    help="TSV file to write (taxid, rank, name).",
)
@click.option(
    "--include-s1",
    is_flag=True,
    help="Include Kraken S1 taxa, which contain many viral species and strains.",
)
def extract_kraken_species_command(inspect_file, output, include_s1):
    """Extract unique species taxids from a Kraken inspect report."""
    try:
        taxa = extract_kraken_species(inspect_file, include_s1)
    except ValueError as error:
        raise click.ClickException(
            f"Unable to parse Kraken inspect report: {error}"
        ) from error

    with open(output, "w", newline="", encoding="utf-8") as destination:
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow(("taxid", "rank", "name"))
        writer.writerows(taxa)

    logger.info(f"Wrote {len(taxa)} unique taxa to {output}")


@cli.command(name="extract-library-species")
@click.argument(
    "library_report",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
)
@click.option(
    "--inspect-file",
    default=None,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    help="Kraken inspect report that supplies canonical taxids and names.",
)
@click.option(
    "--taxonomy-dir",
    default=None,
    type=click.Path(exists=True, file_okay=False, readable=True, path_type=Path),
    help="TaxaForge taxonomy directory with nodes, names, and nucleotide accession maps.",
)
@click.option(
    "--output",
    default="unique-library-species.tsv",
    show_default=True,
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    help="TSV file to write (taxid, rank, name, library_entries).",
)
@click.option(
    "--include-s1",
    is_flag=True,
    help="Include Kraken S1 taxa, which contain many viral species and strains.",
)
@click.option(
    "--workers",
    default=DEFAULT_EXTRACTION_WORKERS,
    show_default=True,
    type=click.IntRange(min=1),
    help="Number of CPU workers for library-record resolution.",
)
def extract_library_species_command(
    library_report, inspect_file, taxonomy_dir, output, include_s1, workers
):
    """Extract unique library species using Kraken or TaxaForge taxonomy."""
    if inspect_file and taxonomy_dir:
        raise click.UsageError(
            "Use either --inspect-file or --taxonomy-dir, not both"
        )
    if not inspect_file and taxonomy_dir is None:
        cache_dir = load_config().get(
            CONFIG_SECTION, "cache-dir", fallback=str(create_cache_dir())
        )
        taxonomy_dir = Path(cache_dir) / "taxonomy"
    if taxonomy_dir and include_s1:
        raise click.UsageError(
            "--include-s1 is only available with --inspect-file"
        )
    try:
        if inspect_file:
            taxa, unmatched = extract_library_species(
                library_report, inspect_file, include_s1, workers
            )
        else:
            taxa, unmatched = extract_library_species_from_taxonomy(
                library_report, taxonomy_dir, workers
            )
    except ValueError as error:
        raise click.ClickException(
            f"Unable to parse library report: {error}"
        ) from error

    with open(output, "w", newline="", encoding="utf-8") as destination:
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow(("taxid", "rank", "name", "library_entries"))
        writer.writerows(taxa)

    logger.info(f"Wrote {len(taxa)} unique taxa to {output}")
    if unmatched:
        if inspect_file:
            logger.warning(
                f"{unmatched} library records did not match an inspect taxon "
                "name; they may use a synonym, a renamed taxon, or a "
                "non-species label."
            )
        else:
            logger.warning(
                f"{unmatched} library records could not be resolved through "
                "the accession maps to a species-level taxon."
            )


@cli.command(name="analyze-genomes")
@click.argument(
    "input_dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
def analyze_genomes(input_dir):
    """Count FASTA and ZIP genome files by organism-group directory."""
    progress_shown = False

    def show_progress(scanned, counts):
        nonlocal progress_shown
        progress_shown = True
        groups = ", ".join(
            f"{group}: {count}"
            for group, count in sorted(counts.items())
        ) or "none found yet"
        sys.stdout.write(
            f"\rScanned {scanned:,} files; genome files found: "
            f"{sum(counts.values()):,} ({groups})"
        )
        sys.stdout.flush()

    counts = analyze_genome_directory(input_dir, progress_callback=show_progress)
    if progress_shown:
        print()
    if not counts:
        logger.info(f"No FASTA or ZIP genome files found under {input_dir}")
        return

    rows = [
        (group, str(count))
        for group, count in sorted(counts.items(), key=lambda item: (item[0] == "unclassified", item[0]))
    ]
    headers = ("GROUP", "FILES")
    widths = [
        max(len(row[index]) for row in (headers, *rows))
        for index in range(len(headers))
    ]
    print("  " + "  ".join(value.ljust(width) for value, width in zip(headers, widths)))
    for row in rows:
        print("  " + "  ".join(value.ljust(width) for value, width in zip(row, widths)))
    print(f"\nTotal genome files: {sum(counts.values())}")
    if counts.get("unclassified"):
        logger.warning(
            "Some files were not under a recognized organism-group directory; "
            "they are counted as 'unclassified'."
        )


@cli.command(name="create-library-report")
@click.argument(
    "input_dir",
    type=click.Path(exists=True, file_okay=False, readable=True, path_type=Path),
)
@click.option(
    "--output",
    default="library_report.csv",
    show_default=True,
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    help="CSV file to write (#Library, Sequence Name, URL).",
)
@click.option(
    "--workers",
    default=DEFAULT_LIBRARY_REPORT_WORKERS,
    show_default=True,
    type=click.IntRange(min=1),
    help="Number of CPU workers used to read FASTA files.",
)
@click.option(
    "--limit",
    default=None,
    type=click.IntRange(min=1),
    help="Maximum number of FASTA files to include.",
)
@click.option(
    "--one-per-file",
    is_flag=True,
    help="Write only the first FASTA record from each file.",
)
@click.option(
    "--taxonomy-dir",
    default=None,
    type=click.Path(exists=True, file_okay=False, readable=True, path_type=Path),
    help=(
        "TaxaForge taxonomy directory used to classify files outside "
        "organism-group directories from their sequence accessions."
    ),
)
def create_library_report_command(
    input_dir, output, workers, limit, one_per_file, taxonomy_dir
):
    """Create a library report from FASTA records in a directory."""
    try:
        rows = iter_library_report_rows(
            input_dir, workers, limit, one_per_file
        )
        unresolved = 0
        if taxonomy_dir:
            rows, unresolved = classify_library_report_rows(
                list(rows), taxonomy_dir, workers
            )
        with open(output, "w", newline="", encoding="utf-8") as destination:
            writer = csv.writer(destination, lineterminator="\n")
            writer.writerow(("#Library", "Sequence Name", "URL"))
            records_written = 0
            for row in rows:
                writer.writerow(row)
                records_written += 1
    except (OSError, UnicodeDecodeError, ValueError) as error:
        raise click.ClickException(f"Unable to read FASTA files: {error}") from error

    logger.info(f"Wrote {records_written} FASTA records to {output}")
    if taxonomy_dir and unresolved:
        logger.warning(
            f"{unresolved} FASTA records could not be classified from "
            "their sequence accessions."
        )


@cli.command(name='doctor')
def doctor():
    for binary in DOCTOR_BINS:
        path = shutil.which(binary)
        status = path if path else "MISSING"
        print(f"{binary:<20} {status}")

    print()
    cache_dir = load_config().get(CONFIG_SECTION, 'cache-dir', fallback=str(create_cache_dir()))
    print(f"Cache dir:   {cache_dir}")
    print(f"Config file: {get_config_path()}")

    cache_status = genome_cache_status(cache_dir)
    if cache_status:
        print()
        print("Genome cache:")
        rows = [
            (group, str(count), format_size(total_size), format_size(avg_size))
            for group, (count, total_size, avg_size) in cache_status.items()
        ]
        headers = ("GROUP", "FILES", "TOTAL", "AVG")
        widths = [max(len(row[i]) for row in (headers, *rows)) for i in range(4)]
        for row in (headers, *rows):
            print("  " + "  ".join(val.ljust(w) for val, w in zip(row, widths)))


if __name__ == '__main__':
    cli()
