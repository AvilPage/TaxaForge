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
