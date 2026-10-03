Installation
============

```bash
pip install taxaforge
```

Usage
=====

```bash
taxaforge --help
```

To create standard Kraken2 database

```bash
taxaforge build --db-type standard
```

Before creating a standard database, you can try a smaller database like fungi.

```bash
taxaforge build --db-type fungi
```

To build a ganon2 database instead, pass `--tool ganon2`

```bash
taxaforge build --tool ganon2 --db-type fungi
```

To use locally downloaded files, run the following command

```bash
taxaforge build --db-name k2_test --genomes-dir /path/to/genomes --taxonomy-dir /path/to/taxonomy
```

To limit the number of genomes in the database, use the `--limit` option

```bash
taxaforge build --db-name k2_test_100 --genomes-dir /path/to/genomes --limit 1000
```

To download and build only the assemblies in an accession list with ganon, use
`--assembly-accessions-file`. The file may contain assembly accessions or FASTA
filenames containing one accession per line; `--limit` selects the first
distinct accessions. Use `--genomes-cache-dir` to reuse matching FASTA files
already present in a directory tree; only requested assemblies missing from
that cache are downloaded. For ganon2 builds, omitted `--kmer-len`, `--min-len`,
and `--level` options are left unset so ganon applies its own defaults.
Additional ganon `build-custom` options can be passed through, for example
`--max-fp 0.01`.

Use `--download-threads` and `--build-threads` to tune download parallelism
separately from index building. Downloads default to at most 4 workers, and
`--download-threads` accepts values from 1 to 4. For example, reduce NCBI
request concurrency while keeping the index build parallel:

```bash
taxaforge build --db-type archaea --download-threads 2 --build-threads 12
```

The legacy `--threads` option remains available and sets both values unless
overridden by a specific option; its download value is capped at 4. Build
threads default to 80% of available CPU threads.

For Kraken2, `--minimizer-len` is an alias for `--min-len` and is passed as
`kraken2-build --minimizer-len`. The minimizer length cannot exceed the k-mer
length; if omitted and a shorter `--kmer-len` is selected, TaxaForge lowers
the default minimizer length to match. Kraken2's default 7 minimizer spaces
also cannot exceed one quarter of the minimizer length. TaxaForge lowers that
default when needed (for example, a minimizer length of 26 uses 6 spaces).
Override it with `--minimizer-spaces`; values above Kraken2's limit are rejected
before genome downloads begin.

Use `--max-db-size` to impose a Kraken2 hash-table budget. Kraken2 downsamples
minimizers to fit the budget, trading sensitivity for a bounded index size:

```bash
taxaforge build --tool kraken2 --db-type bacteria --db-name bacteria-50g \
  --max-db-size 50G
```

The value must be a positive byte count, optionally suffixed with `K`, `M`,
`G`, `T`, or `P`. This option is unavailable for ganon and `--use-k2`.

```bash
taxaforge build --tool ganon2 \
  --assembly-accessions-file /path/to/250k_genomes.txt \
  --genomes-cache-dir /path/to/existing/genomes \
  --limit 1000 --db-name g1000 --threads 4 --max-fp 0.01
```

For a workflow-managed build with nf-core/createtaxdb, install Nextflow, Java
17+, and Docker, then select the `createtaxdb` backend:

```bash
taxaforge build --tool ganon2 --backend createtaxdb \
  --assembly-accessions-file /path/to/250k_genomes.txt \
  --genomes-cache-dir /path/to/existing/genomes \
  --limit 1000 --db-name g1000 --threads 4 --resume
```

This creates a complete samplesheet from the requested FASTA files and NCBI
assembly summaries, then invokes nf-core/createtaxdb with Docker. Missing FASTA
files or taxids are errors rather than silently omitted rows. `--resume` uses
Nextflow's task cache; omit it for a fresh run. Set `NEXTFLOW_BIN`,
`CREATETAXDB_PIPELINE`, and `CREATETAXDB_REVISION` to override the executable,
pipeline, or revision. `--kmer-len`, `--min-len`, and `--level` are passed
through only when specified.

To turn a plain genome filename list into an [nf-core/createtaxdb](https://github.com/nf-core/createtaxdb) samplesheet (`id,taxid,fasta_dna`), use `createtaxdb-input`

```bash
taxaforge createtaxdb-input --input-file 250k_genomes.txt --genomes-dir /path/to/genomes --output samplesheet.csv
```

Add `--download-missing` to auto-fetch any listed genomes not found under `--genomes-dir`.

By default `createtaxdb-input` auto-detects from the accessions in `--input-file` whether it needs NCBI's refseq and/or genbank `assembly_summary` file (genbank's is ~1.75GB, so it's skipped for refseq-only lists). Override with `--taxonomy-source refseq|genbank|both`.

To build the samplesheet without downloading or even having any genomes locally, use `--remote` — it uses each accession's NCBI HTTPS URL (from `assembly_summary`'s `ftp_path`) as `fasta_dna`, since Nextflow/nf-core/createtaxdb can fetch those URLs itself

```bash
taxaforge createtaxdb-input --input-file 250k_genomes.txt --remote --output samplesheet.csv
```

Config
======

Read/write config, stored in the OS default config location (`~/.config/taxaforge/config.ini` on Linux, `~/Library/Application Support/taxaforge/config.ini` on macOS).

```bash
taxaforge config set threads 8
taxaforge config get threads
taxaforge config
```

Download taxonomy
=================

Download or reuse TaxaForge's taxonomy dump and nucleotide accession-to-taxid
maps. These files are required when `extract-library-species` uses its default
taxonomy mode:

```bash
taxaforge download-taxonomy
```

Use `--cache-dir` to store or reuse taxonomy data outside the configured cache:

```bash
taxaforge download-taxonomy --cache-dir /path/to/cache
```

Analyze a genome directory
==========================

Count supported FASTA files by organism-group folder without reading their contents. The scan supports `.fa`, `.fasta`, `.fna` (including gzip-compressed variants), and `.zip` files, and is suitable for large collections:

```bash
taxaforge analyze-genomes /path/to/genomes
```

The command refreshes a live scan count and group totals every 1,000 filesystem files, so large scans show ongoing progress. It looks for recognized group names such as `viral`, `bacteria`, or `archaea` in the scanned directory path and its descendants, so an assembly folder such as `/refseq/archaea/GCF_...` is classified as `archaea`. Files with no recognized group in their path are counted as `unclassified`; the command does not infer taxonomy from FASTA contents or filenames.

Create a library report
=======================

Create a CSV library report from every record in plain or gzip-compressed
FASTA files beneath a directory:

```bash
taxaforge create-library-report /path/to/genomes \
  --output library_report.csv --workers 8
```

The report contains `#Library`, `Sequence Name`, and `URL` columns. TaxaForge
uses the FASTA header as `Sequence Name`, writes the local file URI as `URL`,
and infers `#Library` from a recognized organism-group directory in the file's
path; files outside one are marked `unclassified`. The generated CSV can be
used directly with `extract-library-species`. Without taxonomy resolution it
streams rows to disk and uses the smaller of eight workers or the detected CPU
count minus four (at least one), keeping memory bounded and avoiding excess
disk contention on large machines. Set `--workers` to match the CPU and
storage bandwidth available to the job. Use `--limit N` to process only the
first `N` FASTA files (each selected file contributes all of its sequence
records). For genome-level comparisons, use `--one-per-file` to write only
the first FASTA record from each genome file. Pair that report with
`extract-library-species --taxonomy-dir ...` to produce a compact,
canonical-species summary for comparison with
`data/standard/unique-library-species.tsv`.

For a flat collection such as `250K_Genomes`, pass the TaxaForge taxonomy
directory to classify otherwise-unclassified records from their sequence
accessions:

```bash
taxaforge create-library-report /media/HD/anand/250K_Genomes \
  --taxonomy-dir /path/to/cache/taxonomy \
  --output data/250k/library_report.csv --one-per-file --workers 8
```

The taxonomy directory must contain `nodes.dmp`, `nucl_gb.accession2taxid`,
and `nucl_wgs.accession2taxid`; run `taxaforge download-taxonomy` first if
the maps are absent. If a sequence accession is absent from those maps,
TaxaForge also resolves a `GCF_` or `GCA_` accession in the genome filename
through any cached `assembly_summary_*.txt` files in that directory.

Extract unique Kraken species
==============================

Extract the distinct taxids at Kraken's `S` species rank from an inspect report:

```bash
taxaforge extract-kraken-species data/plusPF/inspect-ppf.txt \
  --output data/plusPF/unique-species.tsv
```

The output is a TSV with `taxid`, `rank`, and `name` columns. Kraken places
many viral species and strain-level taxa at `S1`; include those explicitly when
needed:

```bash
taxaforge extract-kraken-species data/plusPF/inspect-ppf.txt \
  --output data/plusPF/unique-species-and-s1.tsv --include-s1
```

Extract species represented by a Kraken library
===============================================

`library_report.tsv` contains sequence records and does not have taxids, so
TaxaForge resolves its leading organism names against an inspect report and
counts the matched records for each species:

```bash
taxaforge extract-library-species data/plusPF/library_report.tsv \
  --inspect-file data/plusPF/inspect-ppf.txt \
  --output data/plusPF/unique-library-species.tsv --include-s1 \
  --workers 8
```

The output contains `taxid`, `rank`, `name`, and `library_entries`. Records
whose labels use taxonomy synonyms, renamed taxa, or non-species descriptions
are reported as unmatched rather than silently assigned to an incorrect taxon.
TaxaForge uses up to eight CPU workers by default. Use `--workers` to reduce
CPU use or to match the CPUs allocated to the job. For large taxonomy
accession maps, workers are distributed across byte ranges in both map files,
so a higher value can use more than two CPU cores.

Alternatively, use TaxaForge's downloaded taxonomy files. This uses the
nucleotide accession-to-taxid maps and taxonomy lineage, rather than display
names, and does not require an inspect report:

```bash
taxaforge extract-library-species data/plusPF/library_report.tsv \
  --taxonomy-dir /path/to/cache/taxonomy \
  --output data/plusPF/unique-library-species.tsv --workers 8
```

The taxonomy directory must contain `nodes.dmp`, `names.dmp`,
`nucl_gb.accession2taxid`, and `nucl_wgs.accession2taxid`. The accession maps
are unavailable when the original build used `--skip-maps`. When neither
`--taxonomy-dir` nor `--inspect-file` is supplied, TaxaForge defaults to the
`taxonomy` directory in its configured cache.


Why TaxaForge?
==============

TaxaForge aims to provide a simple and easy to use tool to build wide variety of taxonomic classifier databases with a single command.

Why not kraken2-build/ganon directly?

kraken2-build and ganon each build databases for their own tool only. TaxaForge wraps the shared download/taxonomy pipeline once and dispatches to the right tool via `--tool`, so adding support for more classifiers is a matter of plugging in a new build step.


Documentation
=============

- [Kraken2 Database Builder](https://avilpage.com/kdb.html)
- [Mastering Kraken2](https://avilpage.com/tags/kraken2.html)
