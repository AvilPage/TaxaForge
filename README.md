Installation
============

```bash
pip install kraken-db-builder
```

Usage
=====

```bash
kdb --help

kraken-db-builder --help
```

To create standard Kraken2 database

```bash
kdb --db-type standard
```

Before creating a standard database, you can try a smaller database like fungi.

```bash
kdb --db-type fungi
```

To use locally downloaded files, run the following command

```bash
kdb --db-name k2_test --genomes-dir /path/to/genomes --taxonomy-dir /path/to/taxonomy
```

To limit the number of genomes in the database, use the `--limit` option

```bash
kdb --db-name k2_test_100 --genomes-dir /path/to/genomes --limit 1000
```


Why kdb(kraken-db-builder)?
============================

kdb was aimed to created to provide a simple and easy to use tool to build wide variety of databases with a single command.

Why not kraken2-build?

kraken2-build is a tool provided by Kraken2 to build Kraken2 database. It is a great tool but it is limited to build only Kraken2 database. kdb is a more generic tool which can be used to build databases for Kraken2, Centrifuge, Bracken, etc.


Documentation
=============

- [Kraken2 Database Builder](https://avilpage.com/kdb.html)
- [Mastering Kraken2](https://avilpage.com/tags/kraken2.html)
