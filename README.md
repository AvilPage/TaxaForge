Installation
============

```bash
pip install kraken-db-builder
```

Usage
=====

```bash
kdb --help
```

To create standard Kraken2 database

```bash
kdb --db-type standard
```

Before creating standard database, you can try a smaller database like fungi.

```bash
kdb --db-type fungi
```


Why kdb(kraken-db-builder)?
============================

kdb was aimed to created to provide a simple and easy to use tool to build wide variety of databases with a single command.

Why not kraken2-build?

kraken2-build is a tool provided by Kraken2 to build Kraken2 database. It is a great tool but it is limited to build only Kraken2 database. kdb is a more generic tool which can be used to build databases for Kraken2, Centrifuge, Bracken, etc.
