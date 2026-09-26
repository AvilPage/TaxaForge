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

Config
======

Read/write config, stored in the OS default config location (`~/.config/taxaforge/config.ini` on Linux, `~/Library/Application Support/taxaforge/config.ini` on macOS).

```bash
taxaforge config set threads 8
taxaforge config get threads
taxaforge config
```


Why TaxaForge?
==============

TaxaForge aims to provide a simple and easy to use tool to build wide variety of taxonomic classifier databases with a single command.

Why not kraken2-build/ganon directly?

kraken2-build and ganon each build databases for their own tool only. TaxaForge wraps the shared download/taxonomy pipeline once and dispatches to the right tool via `--tool`, so adding support for more classifiers is a matter of plugging in a new build step.


Documentation
=============

- [Kraken2 Database Builder](https://avilpage.com/kdb.html)
- [Mastering Kraken2](https://avilpage.com/tags/kraken2.html)
