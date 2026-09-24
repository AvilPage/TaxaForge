#!/usr/bin/env python3
import concurrent.futures
import configparser
import datetime
import hashlib
import importlib.metadata
import logging
import multiprocessing
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import urllib.request
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

    try:
        response = urllib.request.urlopen(request)
    except urllib.error.HTTPError as error:
        if existing and error.code == 416:
            return
        raise

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


def download_taxanomy(cache_dir, skip_maps=None, protein=None):
    taxonomy_path = os.path.join(cache_dir, "taxonomy")
    os.makedirs(taxonomy_path, exist_ok=True)
    os.chdir(taxonomy_path)

    urls = []
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

    logger.info(f"Downloading {len(urls)} taxonomy files")
    download_files(urls)

    logger.info("Extracting taxdump.tar.gz")
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


def ensure_organism_genomes(cache_dir, organism, threads, seed_dirs=None):
    # ncbi_genome_download skips assemblies already present under cache_dir/refseq/{organism},
    # so re-running for the same organism only fetches what's missing
    import_seed_genomes(Path(cache_dir) / "refseq" / organism, seed_dirs)
    logger.info(f"Downloading genomes for {organism}")
    cwd = os.getcwd()
    os.chdir(cache_dir)
    ncbi_genome_download.download(
        section='refseq', groups=organism, file_formats='fasta',
        progress_bar=True, parallel=threads,
        assembly_levels=['complete'],
        output=cache_dir, uri=f"{NCBI_SERVER}/genomes"
    )
    cmd = f"find {cache_dir}/refseq/{organism} -name '*.gz' | xargs -n 1 -P {threads} gunzip -k"
    run_cmd(cmd)
    os.chdir(cwd)
    logger.info(f"Finished downloading {organism} genomes")


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

    cmd += f" --db {db_name} --threads {threads} --kmer-len {kmer_len} --minimizer-len {min_len} --load-factor {load_factor}"
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


def ensure_assembly_summaries(cache_dir):
    # Cached so `ganon build-custom --ncbi-file-info` doesn't re-download these on every build
    taxonomy_path = os.path.join(cache_dir, "taxonomy")
    os.makedirs(taxonomy_path, exist_ok=True)

    urls = [
        f"{NCBI_SERVER}/genomes/refseq/assembly_summary_refseq.txt",
        f"{NCBI_SERVER}/genomes/genbank/assembly_summary_genbank.txt",
        f"{NCBI_SERVER}/genomes/refseq/assembly_summary_refseq_historical.txt",
        f"{NCBI_SERVER}/genomes/genbank/assembly_summary_genbank_historical.txt",
    ]
    paths = [os.path.join(taxonomy_path, url.rsplit("/", 1)[-1]) for url in urls]
    missing_urls = [url for url, path in zip(urls, paths) if not os.path.exists(path)]

    if missing_urls:
        logger.info(f"Downloading {len(missing_urls)} assembly summary files")
        os.chdir(taxonomy_path)
        download_files(missing_urls)

    return paths


def ensure_genome_size_file(cache_dir):
    # Cached so `ganon build-custom --genome-size-files` doesn't try to (re-)download it, which
    # fails outright on hosts without internet access
    taxonomy_path = os.path.join(cache_dir, "taxonomy")
    os.makedirs(taxonomy_path, exist_ok=True)

    path = os.path.join(taxonomy_path, "species_genome_size.txt.gz")
    if not os.path.exists(path):
        logger.info("Downloading species genome size file")
        os.chdir(taxonomy_path)
        download_files([f"{NCBI_SERVER}/genomes/ASSEMBLY_REPORTS/species_genome_size.txt.gz"])

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


ACCESSION_RE = re.compile(r"^(GC[AF]_\d+\.\d+)")


def build_accession_taxid_map(assembly_summary_paths):
    # Maps assembly_accession (col 1) -> taxid (col 6), same info ganon's --ncbi-file-info
    # extracts internally. Needed so we can pre-resolve target/node ourselves when writing
    # an --input-file tsv for a from/to slice (see run_ganon_build_custom for why).
    accession_to_taxid = {}
    for path in assembly_summary_paths:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) > 5 and cols[0] and cols[5]:
                    accession_to_taxid[cols[0]] = cols[5]
    return accession_to_taxid


RAPTOR_MAX_KMER = 32  # raptor's (h)ibf backend packs k-mers into a 64-bit word (2 bits/base for DNA)


def run_ganon_build_custom(cache_dir, input_dirs, db_name, threads, kmer_len, min_len, level, from_idx=None, to_idx=None):
    if kmer_len > RAPTOR_MAX_KMER:
        logger.warning(
            f"--kmer-len {kmer_len} exceeds raptor's max of {RAPTOR_MAX_KMER} for ganon2 builds; "
            f"clamping to {RAPTOR_MAX_KMER} (raptor prepare would otherwise fail with "
            f"'Value {kmer_len} is not in range [1,{RAPTOR_MAX_KMER}]')"
        )
        kmer_len = RAPTOR_MAX_KMER

    # raptor also requires window-size (minimizer length) >= kmer-size
    if min_len < kmer_len:
        logger.warning(
            f"--min-len {min_len} (window-size) is smaller than --kmer-len {kmer_len}; "
            f"raptor requires window-size >= kmer-size, raising --min-len to {kmer_len}"
        )
        min_len = kmer_len

    taxa_flag = detect_taxonomy_flag()
    if taxa_flag == "--taxonomy-files":
        tax_args = f"{taxa_flag} {cache_dir}/taxonomy/nodes.dmp {cache_dir}/taxonomy/names.dmp"
    else:
        tax_args = f"{taxa_flag} {cache_dir}/taxonomy/taxdump.tar.gz"

    assembly_summary_paths = ensure_assembly_summaries(cache_dir)
    genome_size_file = ensure_genome_size_file(cache_dir)
    input_extension = detect_input_extension(input_dirs)

    # ganon has no 'file' taxonomic rank - 'file' means per-file bins, achieved by
    # omitting --level so it defaults to --input-target (itself defaulting to 'file').
    # Passing --level file gets silently rejected by ganon and falls back to 'leaves'.
    level_args = "" if level == "file" else f" --level {level}"

    input_file_list = None
    ncbi_file_info_args = f"--ncbi-file-info {' '.join(assembly_summary_paths)} "
    if from_idx is not None or to_idx is not None:
        files = sorted(
            str(path) for input_dir in input_dirs
            for path in Path(input_dir).rglob(f"*.{input_extension}")
        )
        logger.info(f"Using genome files range [{from_idx}:{to_idx}] out of {len(files)}")
        files = files[from_idx:to_idx]

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
        input_args = f"--input-file {shlex.quote(input_file_list.name)}"
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
        f"--db-prefix {db_name} --threads {threads} --kmer-size {kmer_len} "
        f"--window-size {min_len}{level_args}"
    )

    with tempfile.TemporaryDirectory(prefix="raptor_shim_") as shim_tmp:
        setup_raptor_shim(Path(shim_tmp))
        old_path = os.environ.get("PATH", "")
        os.environ["PATH"] = f"{shim_tmp}:{old_path}"
        try:
            run_cmd(cmd)
        finally:
            os.environ["PATH"] = old_path
            if input_file_list:
                os.unlink(input_file_list.name)

    cmd = f"du -sh {db_name}.*"
    run_cmd(cmd)


def build_ganon2(cache_dir, cwd, genomes_dir, db_type, db_name, threads, kmer_len, min_len, level, rebuild, from_idx=None, to_idx=None):
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
    run_ganon_build_custom(cache_dir, input_dirs, db_name, threads, kmer_len, min_len, level, from_idx, to_idx)


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
    for kwargs in FILTER_DOWNLOAD_PASSES[filter_code]:
        logger.info(f"Downloading {organism} genomes ({filter_code}) with {kwargs}")
        ncbi_genome_download.download(
            section='refseq', groups=organism, file_formats='fasta',
            progress_bar=True, parallel=threads,
            output=str(filtered_root), uri=f"{NCBI_SERVER}/genomes", **kwargs
        )
    organism_dir = filtered_root / "refseq" / organism
    cmd = f"find {organism_dir} -name '*.gz' | xargs -n 1 -P {threads} gunzip -k"
    run_cmd(cmd)
    return organism_dir


def build_filtered_ganon(cache_dir, db_prefix, organism_groups, filter_code, ganon_args, seed_dirs=None, download_only=False, from_idx=None, to_idx=None):
    # cgrg-style combined filter: download via ncbi-genome-download instead of ganon's own
    # --genome-updater -F filter, then hand the files to build-custom
    logger.info(f"Filtered ({filter_code}) build for organism groups: {', '.join(organism_groups)} via ncbi-genome-download + ganon build-custom")
    threads = int((get_arg_values(ganon_args, '-t', '--threads') or ['4'])[0])
    kmer_len = (get_arg_values(ganon_args, '-k', '--kmer-size') or ['19'])[0]
    min_len = (get_arg_values(ganon_args, '-w', '--window-size') or ['31'])[0]
    level = (get_arg_values(ganon_args, '-l', '--level') or ['species'])[0]

    input_dirs = [str(ensure_filtered_genomes(cache_dir, organism, threads, filter_code, seed_dirs)) for organism in organism_groups]
    if download_only:
        logger.info("Downloaded filtered genomes, skipping build (--download-only)")
        return
    run_ganon_build_custom(cache_dir, input_dirs, db_prefix, threads, kmer_len, min_len, level, from_idx, to_idx)


def detect_filter_code(ganon_args):
    """cgrg (genome-updater) builds route through build_filtered_ganon; plain rg/cg builds
    keep using native `ganon build`, which already works fine for those on their own."""
    has_gu = get_option_value(ganon_args, '--genome-updater') is not None
    has_rg = any(arg in ('-r', '--reference-genomes') for arg in ganon_args)
    has_cg = any(arg in ('-c', '--complete-genomes') for arg in ganon_args)
    if has_gu or (has_rg and has_cg):
        return 'cgrg'
    return None


def build_combo_ganon(cache_dir, db_prefix, organism_groups, ganon_args, seed_dirs=None, download_only=False, from_idx=None, to_idx=None):
    # Combo db-prefix (e.g. archaea+viral): download each organism group into its own
    # persistent cache dir and reuse ganon build-custom, instead of native `ganon build`
    # which would re-download the whole combo as a single genome_updater job every time
    logger.info(f"Combo build for organism groups: {', '.join(organism_groups)} (downloaded/reused separately)")
    threads = int((get_arg_values(ganon_args, '-t', '--threads') or ['4'])[0])
    kmer_len = (get_arg_values(ganon_args, '-k', '--kmer-size') or ['19'])[0]
    min_len = (get_arg_values(ganon_args, '-w', '--window-size') or ['31'])[0]
    level = (get_arg_values(ganon_args, '-l', '--level') or ['species'])[0]

    for organism in organism_groups:
        ensure_organism_genomes(cache_dir, organism, threads, seed_dirs)

    if download_only:
        logger.info("Downloaded genomes, skipping build (--download-only)")
        return

    input_dirs = [f"{cache_dir}/refseq/{organism}" for organism in organism_groups]
    run_ganon_build_custom(cache_dir, input_dirs, db_prefix, threads, kmer_len, min_len, level, from_idx, to_idx)


def get_files(genomes_dir, cache_dir, db_type, db_name, threads):
    if genomes_dir:
        logger.info(f"Adding {genomes_dir} genomes to library")

        cmd = f"find {genomes_dir} -name '*.gz' | xargs -n 1 -P {threads} gunzip -k"
        run_cmd(cmd)

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
@click.option('--db-type', default=None, help='database type to build')
@click.option('--db-name', default=None, help='database name to build')
@click.option('--genomes-dir', default=None, help='Directory containing genomes')
@click.option('--seed-dir', 'seed_dirs', multiple=True, default=lambda: tuple(d for d in load_config().get(CONFIG_SECTION, 'seed-dirs', fallback='').split(',') if d), help='Existing folder(s) with genomes to reuse before downloading; missing ones still get downloaded (repeatable, config key: seed-dirs)')
@click.option('--cache-dir', default=lambda: load_config().get(CONFIG_SECTION, 'cache-dir', fallback=str(create_cache_dir())), help='Cache directory for downloaded genomes/taxonomy (config key: cache-dir)')
@click.option('--output-dir', default=lambda: load_config().get(CONFIG_SECTION, 'output-dir', fallback='.'), help='Directory to build the database in, instead of cwd (config key: output-dir)')
@click.option('--threads', default=max(1, int(multiprocessing.cpu_count() * 0.8)), help='Number of threads to use (default: 80% of CPU threads)', type=int)
@click.option('--load-factor', default=0.7, help='Proportion of the hash table to be populated. Used only for kraken2')
@click.option('--kmer-len', default=None, help='Kmer length in bp/aa. Used only in build task (default: 35 for kraken2, 19 for ganon/ganon2, matching each tool\'s own default)', type=int)
@click.option('--min-len', default=31, help='Minimizer/window length in bp/aa. Used only in build task', type=int)
@click.option('--level', default='leaves', type=click.Choice(['leaves', 'species', 'genus', 'assembly', 'file']), help='Taxonomic level to group sequences by. Used only for ganon2')
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
        tool: str, db_type: str, db_name, cache_dir, output_dir, genomes_dir, seed_dirs,
        threads, load_factor, kmer_len: int, min_len, level: str, limit: int, from_idx: int, to_idx: int, batch_size: int,
        force: bool, rebuild, fast_build: bool, use_k2: bool, download_only: bool, ganon_args
):
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(output_dir)
    logger.info(f"Building in {output_dir}")

    if kmer_len is None:
        # ganon's own default k-mer size is 19; kraken2-build's is 35. Use whichever
        # matches the selected tool instead of forcing one default on both.
        kmer_len = 19 if tool in ('ganon', 'ganon2') else 35

    if tool in ('ganon', 'ganon2') and ganon_args:
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
            run_ganon_build_custom(cache_dir, [str(genomes_dir)], db_prefix_value, threads, kmer_len, min_len, level, from_idx, to_idx)
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
            build_filtered_ganon(cache_dir, db_prefix_value, organism_groups, filter_code, ganon_args, seed_dirs, download_only, from_idx, to_idx)
            return
        if db_prefix_value and len(organism_groups) > 1:
            build_combo_ganon(cache_dir, db_prefix_value, organism_groups, ganon_args, seed_dirs, download_only, from_idx, to_idx)
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

    if not genomes_dir:
        download_genomes(cache_dir, cwd, db_type, db_name, threads, force, seed_dirs)

    if download_only:
        logger.info("Downloaded genomes/taxonomy, skipping build (--download-only)")
        return

    if tool == 'kraken2':
        add_to_library(
            cache_dir, cwd, genomes_dir, db_type, db_name,
            limit, from_idx, to_idx, batch_size, threads, use_k2
        )
        build_db(
            cache_dir, cwd, db_type, db_name, threads, kmer_len, min_len,
            fast_build, rebuild, load_factor, use_k2
        )
    elif tool == 'ganon2':
        build_ganon2(
            cache_dir, cwd, genomes_dir, db_type, db_name, threads,
            kmer_len, min_len, level, rebuild, from_idx, to_idx
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
