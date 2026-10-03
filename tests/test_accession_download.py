import csv
import gzip
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import taxaforge
from click.testing import CliRunner


class AssemblyAccessionDownloadTests(unittest.TestCase):
    def test_download_taxonomy_uses_configured_cache_directory(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "cache"
            config = taxaforge.configparser.ConfigParser()
            config.add_section(taxaforge.CONFIG_SECTION)
            config.set(taxaforge.CONFIG_SECTION, "cache-dir", str(cache_dir))

            with patch.object(
                taxaforge, "load_config", return_value=config
            ), patch.object(taxaforge, "download_taxanomy") as download:
                result = CliRunner().invoke(
                    taxaforge.cli, ["download-taxonomy"]
                )

            self.assertEqual(result.exit_code, 0, result.output)
            download.assert_called_once_with(cache_dir)

    def test_extract_kraken_species_writes_unique_s_rank_taxa(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            inspect_file = root / "inspect.txt"
            output_file = root / "species.tsv"
            inspect_file.write_text(
                "# Database options: nucleotide db\n"
                "100.00\t100\t10\tR\t1\troot\n"
                " 10.00\t10\t10\tS\t20\t  Species beta\n"
                "  5.00\t5\t5\tS1\t21\t    Species beta strain A\n"
                "  4.00\t4\t4\tS\t10\t  Species alpha\n"
                "  3.00\t3\t3\tS\t20\t  Species beta\n",
                encoding="utf-8",
            )

            result = CliRunner().invoke(
                taxaforge.cli,
                [
                    "extract-kraken-species",
                    str(inspect_file),
                    "--output",
                    str(output_file),
                ],
            )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(
                output_file.read_text(encoding="utf-8"),
                "taxid\trank\tname\n"
                "10\tS\tSpecies alpha\n"
                "20\tS\tSpecies beta\n",
            )

    def test_extract_kraken_species_can_include_s1_taxa(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            inspect_file = root / "inspect.txt"
            inspect_file.write_text(
                " 10.00\t10\t10\tS\t20\t  Species beta\n"
                "  5.00\t5\t5\tS1\t21\t    Species beta strain A\n",
                encoding="utf-8",
            )

            result = CliRunner().invoke(
                taxaforge.cli,
                [
                    "extract-kraken-species",
                    str(inspect_file),
                    "--output",
                    str(root / "species.tsv"),
                    "--include-s1",
                ],
            )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(
                (root / "species.tsv").read_text(encoding="utf-8"),
                "taxid\trank\tname\n"
                "20\tS\tSpecies beta\n"
                "21\tS1\tSpecies beta strain A\n",
            )

    def test_extract_kraken_species_rejects_malformed_records(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            inspect_file = Path(temp_dir) / "inspect.txt"
            inspect_file.write_text("invalid record\n", encoding="utf-8")

            result = CliRunner().invoke(
                taxaforge.cli,
                ["extract-kraken-species", str(inspect_file)],
            )

            self.assertNotEqual(result.exit_code, 0)
            self.assertIn("Expected 6 tab-separated fields on line 1", result.output)

    def test_extract_library_species_resolves_and_counts_taxa(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            inspect_file = root / "inspect.txt"
            library_report = root / "library.tsv"
            output_file = root / "species.tsv"
            inspect_file.write_text(
                " 10.00\t10\t10\tS\t20\t  Species beta\n"
                "  5.00\t5\t5\tS1\t21\t    Species beta strain A\n"
                "  4.00\t4\t4\tS\t10\t  Species alpha\n",
                encoding="utf-8",
            )
            library_report.write_text(
                "#Library\tSequence Name\tURL\n"
                "bacteria\tNC_000001.1 Species beta strain X, complete genome\thttps://example.test/1\n"
                "plasmid\tNC_000002.1 Species beta plasmid p1\thttps://example.test/2\n"
                "bacteria\tNC_000003.1 MAG: Species alpha isolate 1\thttps://example.test/3\n"
                "UniVec_Core\tgnl|uv|123 Vector pTest\thttps://example.test/4\n",
                encoding="utf-8",
            )

            result = CliRunner().invoke(
                taxaforge.cli,
                [
                    "extract-library-species",
                    str(library_report),
                    "--inspect-file",
                    str(inspect_file),
                    "--output",
                    str(output_file),
                ],
            )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(
                output_file.read_text(encoding="utf-8"),
                "taxid\trank\tname\tlibrary_entries\n"
                "10\tS\tSpecies alpha\t1\n"
                "20\tS\tSpecies beta\t2\n",
            )

    def test_extract_library_species_parallelizes_large_reports(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            inspect_file = root / "inspect.txt"
            library_report = root / "library.tsv"
            inspect_file.write_text(
                " 10.00\t10\t10\tS\t20\t  Species beta\n",
                encoding="utf-8",
            )
            library_report.write_text(
                "#Library\tSequence Name\tURL\n"
                + "".join(
                    "bacteria\tNC_000001.1 Species beta strain "
                    f"{index}\thttps://example.test/{index}\n"
                    for index in range(2_000)
                ),
                encoding="utf-8",
            )

            taxa, unmatched = taxaforge.extract_library_species(
                library_report, inspect_file, workers=2
            )

            self.assertEqual(taxa, [("20", "S", "Species beta", 2_000)])
            self.assertEqual(unmatched, 0)

    def test_extract_library_species_rejects_invalid_header(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            inspect_file = root / "inspect.txt"
            library_report = root / "library.tsv"
            inspect_file.write_text(
                " 10.00\t10\t10\tS\t20\t  Species beta\n",
                encoding="utf-8",
            )
            library_report.write_text(
                "Library\tSequence Name\tURL\n",
                encoding="utf-8",
            )

            result = CliRunner().invoke(
                taxaforge.cli,
                [
                    "extract-library-species",
                    str(library_report),
                    "--inspect-file",
                    str(inspect_file),
                ],
            )

            self.assertNotEqual(result.exit_code, 0)
            self.assertIn("Expected header", result.output)

    def test_extract_library_species_resolves_taxonomy_accession_maps(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            taxonomy_dir = root / "taxonomy"
            taxonomy_dir.mkdir()
            library_report = root / "library.tsv"
            output_file = root / "species.tsv"
            (taxonomy_dir / "nodes.dmp").write_text(
                "1\t|\t1\t|\tno rank\t|\n"
                "10\t|\t1\t|\tspecies\t|\n"
                "11\t|\t10\t|\tsubspecies\t|\n",
                encoding="utf-8",
            )
            (taxonomy_dir / "names.dmp").write_text(
                "10\t|\tSpecies alpha\t|\t\t|\tscientific name\t|\n",
                encoding="utf-8",
            )
            (taxonomy_dir / "nucl_gb.accession2taxid").write_text(
                "accession\taccession.version\ttaxid\tgi\n"
                "NC_000001\tNC_000001.1\t11\t0\n",
                encoding="utf-8",
            )
            (taxonomy_dir / "nucl_wgs.accession2taxid").write_text(
                "accession\taccession.version\ttaxid\tgi\n",
                encoding="utf-8",
            )
            library_report.write_text(
                "#Library\tSequence Name\tURL\n"
                "bacteria\tNC_000001.1 Species alpha strain A\thttps://example.test/1\n"
                "plasmid\tNC_000001.1 Species alpha plasmid p1\thttps://example.test/2\n",
                encoding="utf-8",
            )

            config = taxaforge.configparser.ConfigParser()
            config.add_section(taxaforge.CONFIG_SECTION)
            config.set(taxaforge.CONFIG_SECTION, "cache-dir", str(root))
            with patch.object(taxaforge, "load_config", return_value=config):
                result = CliRunner().invoke(
                    taxaforge.cli,
                    [
                        "extract-library-species",
                        str(library_report),
                        "--output",
                        str(output_file),
                    ],
                )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(
                output_file.read_text(encoding="utf-8"),
                "taxid\trank\tname\tlibrary_entries\n"
                "10\tspecies\tSpecies alpha\t2\n",
            )

    def test_accession_map_ranges_use_requested_workers_without_losing_rows(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            map_path = root / "nucl_gb.accession2taxid"
            wgs_map_path = root / "nucl_wgs.accession2taxid"
            map_path.write_text(
                "accession\taccession.version\ttaxid\tgi\n"
                "NC_000001\tNC_000001.1\t10\t0\n"
                "NC_000002\tNC_000002.1\t20\t0\n"
                "NC_000003\tNC_000003.1\t30\t0\n",
                encoding="utf-8",
            )
            wgs_map_path.write_text(
                "accession\taccession.version\ttaxid\tgi\n",
                encoding="utf-8",
            )

            wanted = {
                b"NC_000001.1": "NC_000001.1",
                b"NC_000002.1": "NC_000002.1",
                b"NC_000003.1": "NC_000003.1",
            }
            resolved = {}
            for task in taxaforge.accession_map_ranges(map_path, 4):
                resolved.update(
                    taxaforge.resolve_accession_taxids_range(
                        *task, wanted, {}
                    )
                )

            with patch.object(
                taxaforge, "ACCESSION_MAP_PARALLEL_MIN_BYTES", 1
            ):
                worker_counts = taxaforge.accession_map_workers(
                    [map_path, wgs_map_path], 70
                )

            self.assertEqual(
                resolved,
                {
                    "NC_000001.1": "10",
                    "NC_000002.1": "20",
                    "NC_000003.1": "30",
                },
            )
            self.assertEqual(sum(worker_counts.values()), 70)

    def test_extract_library_species_rejects_both_taxonomy_sources(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            library_report = root / "library.tsv"
            inspect_file = root / "inspect.txt"
            taxonomy_dir = root / "taxonomy"
            taxonomy_dir.mkdir()
            library_report.write_text(
                "#Library\tSequence Name\tURL\n",
                encoding="utf-8",
            )
            inspect_file.write_text("", encoding="utf-8")

            result = CliRunner().invoke(
                taxaforge.cli,
                [
                    "extract-library-species",
                    str(library_report),
                    "--inspect-file",
                    str(inspect_file),
                    "--taxonomy-dir",
                    str(taxonomy_dir),
                ],
            )

            self.assertNotEqual(result.exit_code, 0)
            self.assertIn(
                "Use either --inspect-file or --taxonomy-dir, not both",
                result.output,
            )

    def test_analyze_genome_directory_counts_supported_files_by_group(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "refseq" / "viral").mkdir(parents=True)
            (root / "genbank" / "bacteria").mkdir(parents=True)
            (root / "GCF_000000001.1_genomic.fa").touch()
            (root / "refseq" / "viral" / "GCF_000000002.1_genomic.fna.gz").touch()
            (root / "refseq" / "viral" / "GCF_000000003.1.zip").touch()
            (root / "genbank" / "bacteria" / "GCA_000000004.1.fa").touch()
            (root / "ignored.txt").touch()

            self.assertEqual(
                taxaforge.analyze_genome_directory(root),
                {"unclassified": 1, "viral": 2, "bacteria": 1},
            )

    def test_analyze_genome_directory_uses_group_root_and_normalizes_aliases(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            viral_dir = Path(temp_dir) / "Viruses"
            viral_dir.mkdir()
            (viral_dir / "genome.FA").touch()

            self.assertEqual(
                taxaforge.analyze_genome_directory(viral_dir), {"viral": 1}
            )

    def test_analyze_genome_directory_uses_group_ancestor_above_scan_root(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            assembly_dir = Path(temp_dir) / "refseq" / "archaea" / "GCF_000980155.1"
            assembly_dir.mkdir(parents=True)
            (assembly_dir / "genomic.fna").touch()

            self.assertEqual(
                taxaforge.analyze_genome_directory(assembly_dir), {"archaea": 1}
            )

    def test_analyze_genome_directory_reports_progress_every_thousand_files(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "bacteria").mkdir()
            for index in range(2001):
                (root / "bacteria" / f"genome-{index}.fa").touch()
            updates = []

            result = taxaforge.analyze_genome_directory(
                root,
                progress_callback=lambda scanned, counts: updates.append(
                    (scanned, counts["bacteria"])
                ),
            )

            self.assertEqual(updates, [(1000, 1000), (2000, 2000), (2001, 2001)])
            self.assertEqual(result, {"bacteria": 2001})

    def test_analyze_genome_directory_rejects_invalid_progress_interval(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            with self.assertRaisesRegex(ValueError, "at least 1"):
                taxaforge.analyze_genome_directory(
                    temp_dir, progress_interval=0
                )

    def test_create_library_report_reads_plain_and_gzip_fasta_files(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            bacteria = root / "bacteria"
            viral = root / "viruses"
            bacteria.mkdir()
            viral.mkdir()
            (bacteria / "genome.fna").write_text(
                ">NC_000001.1 Species alpha chromosome\nACGT\n"
                ">NC_000002.1 Species alpha plasmid\nTGCA\n",
                encoding="utf-8",
            )
            with gzip.open(viral / "genome.fa.gz", "wt", encoding="utf-8") as source:
                source.write(">NC_000003.1 Virus beta complete genome\nACGT\n")
            (root / "ignored.txt").write_text(">ignored\n", encoding="utf-8")
            output = root / "library_report.csv"

            result = CliRunner().invoke(
                taxaforge.cli,
                [
                    "create-library-report",
                    str(root),
                    "--output",
                    str(output),
                    "--workers",
                    "2",
                ],
            )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(
                output.read_text(encoding="utf-8"),
                "#Library,Sequence Name,URL\n"
                "bacteria,NC_000001.1 Species alpha chromosome,"
                f"{(bacteria / 'genome.fna').resolve().as_uri()}\n"
                "bacteria,NC_000002.1 Species alpha plasmid,"
                f"{(bacteria / 'genome.fna').resolve().as_uri()}\n"
                "viral,NC_000003.1 Virus beta complete genome,"
                f"{(viral / 'genome.fa.gz').resolve().as_uri()}\n",
            )
            self.assertEqual(
                taxaforge.read_library_report(output),
                [
                    [
                        "bacteria",
                        "NC_000001.1 Species alpha chromosome",
                        (bacteria / "genome.fna").resolve().as_uri(),
                    ],
                    [
                        "bacteria",
                        "NC_000002.1 Species alpha plasmid",
                        (bacteria / "genome.fna").resolve().as_uri(),
                    ],
                    [
                        "viral",
                        "NC_000003.1 Virus beta complete genome",
                        (viral / "genome.fa.gz").resolve().as_uri(),
                    ],
                ],
            )

    def test_create_library_report_classifies_flat_files_from_taxonomy(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            taxonomy_dir = root / "taxonomy"
            taxonomy_dir.mkdir()
            (root / "bacterial.fna").write_text(
                ">NC_000001.1 Bacterium alpha chromosome\nACGT\n",
                encoding="utf-8",
            )
            (root / "human.fna").write_text(
                ">NC_000002.1 Homo sapiens chromosome\nACGT\n",
                encoding="utf-8",
            )
            (root / "GCF_000003.1_fungal.fna").write_text(
                ">NW_000003.1 Fungus alpha scaffold\nACGT\n",
                encoding="utf-8",
            )
            (taxonomy_dir / "nodes.dmp").write_text(
                "1\t|\t1\t|\tno rank\t|\n"
                "2\t|\t1\t|\tsuperkingdom\t|\n"
                "9606\t|\t1\t|\tspecies\t|\n"
                "4751\t|\t1\t|\tkingdom\t|\n"
                "12345\t|\t4751\t|\tspecies\t|\n",
                encoding="utf-8",
            )
            (taxonomy_dir / "assembly_summary_refseq.txt").write_text(
                "#assembly_accession\ttaxid\n"
                "GCF_000003.1\t12345\n",
                encoding="utf-8",
            )
            for filename in (
                "nucl_gb.accession2taxid",
                "nucl_wgs.accession2taxid",
            ):
                (taxonomy_dir / filename).write_text(
                    "accession\taccession.version\ttaxid\tgi\n"
                    + (
                        "NC_000001\tNC_000001.1\t2\t0\n"
                        "NC_000002\tNC_000002.1\t9606\t0\n"
                        if filename == "nucl_gb.accession2taxid"
                        else ""
                    ),
                    encoding="utf-8",
                )
            output = root / "library_report.csv"

            result = CliRunner().invoke(
                taxaforge.cli,
                [
                    "create-library-report",
                    str(root),
                    "--output",
                    str(output),
                    "--taxonomy-dir",
                    str(taxonomy_dir),
                    "--workers",
                    "1",
                ],
            )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(
                taxaforge.read_library_report(output),
                [
                    [
                        "fungi",
                        "NW_000003.1 Fungus alpha scaffold",
                        (root / "GCF_000003.1_fungal.fna").resolve().as_uri(),
                    ],
                    [
                        "bacteria",
                        "NC_000001.1 Bacterium alpha chromosome",
                        (root / "bacterial.fna").resolve().as_uri(),
                    ],
                    [
                        "human",
                        "NC_000002.1 Homo sapiens chromosome",
                        (root / "human.fna").resolve().as_uri(),
                    ],
                ],
            )

    def test_create_library_report_rejects_invalid_fasta(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "invalid.fasta").write_text("ACGT\n", encoding="utf-8")

            result = CliRunner().invoke(
                taxaforge.cli, ["create-library-report", str(root)]
            )

            self.assertNotEqual(result.exit_code, 0)
            self.assertIn("Expected a FASTA header", result.output)

    def test_create_library_report_limits_fasta_files(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "first.fasta").write_text(
                ">first-a\nACGT\n>first-b\nTGCA\n", encoding="utf-8"
            )
            (root / "second.fasta").write_text(
                ">second\nACGT\n", encoding="utf-8"
            )
            output = root / "library_report.csv"

            result = CliRunner().invoke(
                taxaforge.cli,
                [
                    "create-library-report",
                    str(root),
                    "--output",
                    str(output),
                    "--workers",
                    "1",
                    "--limit",
                    "1",
                ],
            )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(
                taxaforge.read_library_report(output),
                [
                    [
                        "unclassified",
                        "first-a",
                        (root / "first.fasta").resolve().as_uri(),
                    ],
                    [
                        "unclassified",
                        "first-b",
                        (root / "first.fasta").resolve().as_uri(),
                    ],
                ],
            )

    def test_create_library_report_can_write_one_record_per_file(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "first.fasta").write_text(
                ">first-a\nACGT\n>first-b\nTGCA\n", encoding="utf-8"
            )
            (root / "second.fasta").write_text(
                ">second\nACGT\n", encoding="utf-8"
            )
            output = root / "library_report.csv"

            result = CliRunner().invoke(
                taxaforge.cli,
                [
                    "create-library-report",
                    str(root),
                    "--output",
                    str(output),
                    "--workers",
                    "1",
                    "--one-per-file",
                ],
            )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(
                taxaforge.read_library_report(output),
                [
                    [
                        "unclassified",
                        "first-a",
                        (root / "first.fasta").resolve().as_uri(),
                    ],
                    [
                        "unclassified",
                        "second",
                        (root / "second.fasta").resolve().as_uri(),
                    ],
                ],
            )

    def test_decompression_redownloads_corrupt_cached_gzip_once(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            genome_dir = Path(temp_dir)
            compressed_path = genome_dir / "GCF_000001985.1_genomic.fna.gz"
            compressed_path.write_bytes(b"truncated gzip data")
            downloads = []

            def retry_download():
                downloads.append(True)
                with gzip.open(compressed_path, "wb") as compressed_file:
                    compressed_file.write(b">assembly\nACGT\n")

            taxaforge.decompress_gzip_files(
                genome_dir, threads=2, retry_download=retry_download
            )

            self.assertEqual(downloads, [True])
            self.assertEqual(
                compressed_path.with_suffix("").read_bytes(), b">assembly\nACGT\n"
            )

    def test_decompression_reports_corrupt_gzip_when_redownload_is_unavailable(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            compressed_path = Path(temp_dir) / "GCF_000001985.1_genomic.fna.gz"
            compressed_path.write_bytes(b"truncated gzip data")

            with self.assertRaisesRegex(
                taxaforge.click.ClickException,
                r"corrupt gzip file.*GCF_000001985\.1_genomic\.fna\.gz",
            ):
                taxaforge.decompress_gzip_files(Path(temp_dir), threads=1)

    def test_decompression_identifies_ncbi_error_pages_and_skips_retries(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            compressed_path = Path(temp_dir) / "GCF_000001985.1_genomic.fna.gz"
            compressed_path.write_bytes(b"<?xml version='1.0'?><Error>blocked</Error>")
            retry_download = MagicMock()

            with self.assertRaisesRegex(
                taxaforge.click.ClickException,
                r"non-gzip XML/HTML content.*NCBI returned an error page",
            ):
                taxaforge.decompress_gzip_files(
                    Path(temp_dir), threads=1, retry_download=retry_download
                )

            retry_download.assert_not_called()
            self.assertFalse(compressed_path.exists())

    def test_organism_download_retries_corrupt_cache_and_restores_working_directory(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "cache"
            genome_dir = cache_dir / "refseq" / "archaea"
            genome_dir.mkdir(parents=True)
            compressed_path = genome_dir / "GCF_000001985.1_genomic.fna.gz"
            download_count = []

            def fake_download(**kwargs):
                download_count.append(True)
                if len(download_count) == 1:
                    compressed_path.write_bytes(b"truncated gzip data")
                else:
                    with gzip.open(compressed_path, "wb") as compressed_file:
                        compressed_file.write(b">assembly\nACGT\n")

            original_cwd = os.getcwd()
            with patch.object(
                taxaforge.ncbi_genome_download, "download", side_effect=fake_download
            ):
                taxaforge.ensure_organism_genomes(cache_dir, "archaea", threads=1)

            self.assertEqual(os.getcwd(), original_cwd)
            self.assertEqual(len(download_count), 2)
            self.assertEqual(
                compressed_path.with_suffix("").read_bytes(), b">assembly\nACGT\n"
            )

    def test_download_file_retries_transient_dns_failures(self):
        response = MagicMock()
        response.status = 200
        response.headers.get.return_value = "0"
        response.read.return_value = b""
        response.__enter__.return_value = response
        dns_error = taxaforge.urllib.error.URLError(
            taxaforge.socket.gaierror(-3, "Temporary failure in name resolution")
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            original_cwd = os.getcwd()
            os.chdir(temp_dir)
            try:
                with patch.object(
                    taxaforge.urllib.request,
                    "urlopen",
                    side_effect=[dns_error, response],
                ) as urlopen, patch.object(taxaforge.time, "sleep") as sleep:
                    taxaforge.download_file("https://example.test/taxdump.tar.gz", 0)
            finally:
                os.chdir(original_cwd)

        self.assertEqual(urlopen.call_count, 2)
        sleep.assert_called_once_with(1)

    def test_download_file_resumes_partial_file_with_http_range(self):
        response = MagicMock()
        response.status = 206
        response.headers.get.return_value = "3"
        response.read.side_effect = [b"def", b""]
        response.__enter__.return_value = response

        with tempfile.TemporaryDirectory() as temp_dir:
            original_cwd = os.getcwd()
            os.chdir(temp_dir)
            try:
                Path("taxonomy.gz").write_bytes(b"abc")
                with patch.object(
                    taxaforge.urllib.request, "urlopen", return_value=response
                ) as urlopen:
                    taxaforge.download_file(
                        "https://example.test/taxonomy.gz", 0
                    )
            finally:
                os.chdir(original_cwd)

            self.assertEqual(Path(temp_dir, "taxonomy.gz").read_bytes(), b"abcdef")
        request = urlopen.call_args.args[0]
        self.assertEqual(request.get_header("Range"), "bytes=3-")

    def test_download_file_reports_dns_failure_after_retry_limit(self):
        dns_error = taxaforge.urllib.error.URLError(
            taxaforge.socket.gaierror(-3, "Temporary failure in name resolution")
        )

        with patch.object(
            taxaforge.urllib.request, "urlopen", side_effect=dns_error
        ) as urlopen, patch.object(taxaforge.time, "sleep"):
            with self.assertRaisesRegex(
                taxaforge.click.ClickException,
                r"Unable to download https://example\.test/taxdump\.tar\.gz "
                r"after 3 attempts",
            ):
                taxaforge.download_file(
                    "https://example.test/taxdump.tar.gz", 0
                )

        self.assertEqual(urlopen.call_count, 3)

    def test_cached_taxonomy_files_are_not_downloaded_or_extracted_again(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir)
            taxonomy_dir = cache_dir / "taxonomy"
            taxonomy_dir.mkdir()
            (taxonomy_dir / "nodes.dmp").write_bytes(b"nodes")
            (taxonomy_dir / "names.dmp").write_bytes(b"names")

            accession_map = taxonomy_dir / "nucl_gb.accession2taxid.gz"
            with gzip.open(accession_map, "wb") as compressed_file:
                compressed_file.write(b"accession map")

            with patch.object(taxaforge, "download_files") as download_files:
                taxaforge.download_taxanomy(cache_dir, skip_maps=True)
                taxaforge.download_taxanomy(cache_dir, skip_maps=True)

            download_files.assert_not_called()
            self.assertEqual((taxonomy_dir / "nodes.dmp").read_bytes(), b"nodes")
            self.assertEqual((taxonomy_dir / "names.dmp").read_bytes(), b"names")
            self.assertEqual(
                (taxonomy_dir / "nucl_gb.accession2taxid").read_bytes(),
                b"accession map",
            )
            self.assertFalse((taxonomy_dir / "taxdump.tar.gz").exists())

    def test_taxonomy_download_replaces_corrupt_accession_map_and_partial_output(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir)
            taxonomy_dir = cache_dir / "taxonomy"
            taxonomy_dir.mkdir()
            (taxonomy_dir / "nodes.dmp").write_text("nodes", encoding="utf-8")
            (taxonomy_dir / "names.dmp").write_text("names", encoding="utf-8")
            corrupt_map = taxonomy_dir / "nucl_gb.accession2taxid.gz"
            corrupt_map.write_bytes(b"truncated gzip data")
            corrupt_map.with_suffix("").write_bytes(b"partial output")
            wgs_map = taxonomy_dir / "nucl_wgs.accession2taxid.gz"
            with gzip.open(wgs_map, "wb") as compressed_file:
                compressed_file.write(b"wgs map")

            def redownload(url, _position):
                self.assertTrue(url.endswith("nucl_gb.accession2taxid.gz"))
                with gzip.open(corrupt_map, "wb") as compressed_file:
                    compressed_file.write(b"replacement map")

            with patch.object(taxaforge, "download_file", side_effect=redownload):
                taxaforge.download_taxanomy(cache_dir)

            self.assertEqual(
                corrupt_map.with_suffix("").read_bytes(), b"replacement map"
            )
            self.assertTrue(
                corrupt_map.with_suffix("").with_name(
                    "nucl_gb.accession2taxid.complete"
                ).is_file()
            )

    def test_run_cmd_propagates_command_failures(self):
        error = subprocess.CalledProcessError(2, "ganon build-custom")
        with patch.object(taxaforge.subprocess, "run", side_effect=error):
            with self.assertRaisesRegex(taxaforge.click.ClickException, "exit code 2"):
                taxaforge.run_cmd("ganon build-custom --help")

    def test_run_cmd_returns_no_entries_for_empty_output(self):
        with patch.object(taxaforge.subprocess, "check_output", return_value=b""):
            self.assertEqual(taxaforge.run_cmd("find missing-dir", return_output=True), [])

    def test_organism_download_fails_clearly_when_ncbi_is_unreachable_and_cache_is_empty(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "cache"
            (cache_dir / "refseq" / "archaea").mkdir(parents=True)

            with patch.object(taxaforge.ncbi_genome_download, "download", return_value=75):
                with self.assertRaisesRegex(
                    taxaforge.click.ClickException,
                    "Unable to download archaea genomes from NCBI.*no cached FASTA",
                ):
                    taxaforge.ensure_organism_genomes(cache_dir, "archaea", threads=1)

    def test_organism_download_uses_existing_genomes_when_ncbi_is_unreachable(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "cache"
            genome_dir = cache_dir / "refseq" / "archaea"
            genome_dir.mkdir(parents=True)
            genome_path = genome_dir / "cached_genome.fna"
            genome_path.write_text(">assembly\nACGT\n", encoding="utf-8")

            with patch.object(taxaforge.ncbi_genome_download, "download", return_value=75):
                taxaforge.ensure_organism_genomes(cache_dir, "archaea", threads=1)

            self.assertTrue(genome_path.is_file())

    def test_organism_download_refreshes_stale_checksum_manifests_once(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "cache"
            genome_dir = cache_dir / "refseq" / "archaea"
            assembly_dir = genome_dir / "GCF_000001985.1"
            assembly_dir.mkdir(parents=True)
            checksum_path = assembly_dir / "MD5SUMS"
            checksum_path.write_text("stale checksums\n", encoding="utf-8")
            genome_path = assembly_dir / "GCF_000001985.1_genomic.fna.gz"
            downloads = []

            def fake_download(**_kwargs):
                downloads.append(True)
                if len(downloads) == 2:
                    with gzip.open(genome_path, "wb") as compressed_file:
                        compressed_file.write(b">assembly\nACGT\n")

            with patch.object(
                taxaforge.ncbi_genome_download, "download", side_effect=fake_download
            ):
                taxaforge.ensure_organism_genomes(cache_dir, "archaea", threads=1)

            self.assertEqual(len(downloads), 2)
            self.assertFalse(checksum_path.exists())
            self.assertEqual(
                genome_path.with_suffix("").read_bytes(), b">assembly\nACGT\n"
            )

    def test_add_to_library_fails_clearly_when_no_genomes_are_found(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            original_cwd = os.getcwd()
            try:
                with patch.object(taxaforge, "get_files", return_value=[]):
                    with self.assertRaisesRegex(
                        taxaforge.click.ClickException, "No genome FASTA files"
                    ):
                        taxaforge.add_to_library(
                            root, root, None, "archaea", "db", None, None, None,
                            100, 1, False,
                        )
            finally:
                os.chdir(original_cwd)

    def test_cli_separates_download_threads_from_index_build_threads(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            original_cwd = os.getcwd()
            try:
                with patch.object(taxaforge, "run_basic_checks"), patch.object(
                    taxaforge, "download_taxanomy"
                ), patch.object(taxaforge, "download_genomes") as download_genomes, patch.object(
                    taxaforge, "add_to_library"
                ) as add_to_library, patch.object(taxaforge, "build_db") as build_db:
                    result = CliRunner().invoke(
                        taxaforge.cli,
                        [
                            "build",
                            "--db-type", "archaea",
                            "--db-name", "test",
                            "--cache-dir", str(root / "cache"),
                            "--output-dir", str(root / "output"),
                            "--download-threads", "2",
                            "--build-threads", "12",
                        ],
                    )
            finally:
                os.chdir(original_cwd)

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(download_genomes.call_args.args[4], 2)
            self.assertEqual(add_to_library.call_args.args[9], 12)
            self.assertEqual(build_db.call_args.args[4], 12)

    def test_cli_minimizer_len_alias_is_passed_to_kraken2_build(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            original_cwd = os.getcwd()
            try:
                with patch.object(taxaforge, "run_basic_checks"), patch.object(
                    taxaforge, "download_taxanomy"
                ), patch.object(taxaforge, "download_genomes"), patch.object(
                    taxaforge, "add_to_library"
                ), patch.object(taxaforge, "build_db") as build_db:
                    result = CliRunner().invoke(
                        taxaforge.cli,
                        [
                            "build",
                            "--db-type", "archaea",
                            "--db-name", "test",
                            "--cache-dir", str(root / "cache"),
                            "--output-dir", str(root / "output"),
                            "--kmer-len", "30",
                            "--minimizer-len", "29",
                        ],
                    )
            finally:
                os.chdir(original_cwd)

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(build_db.call_args.args[5:7], (30, 29))

    def test_cli_caps_default_kraken2_minimizer_to_kmer_length(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            original_cwd = os.getcwd()
            try:
                with patch.object(taxaforge, "run_basic_checks"), patch.object(
                    taxaforge, "download_taxanomy"
                ), patch.object(taxaforge, "download_genomes"), patch.object(
                    taxaforge, "add_to_library"
                ), patch.object(taxaforge, "build_db") as build_db, self.assertLogs(
                    taxaforge.logger, level="WARNING"
                ) as logs:
                    result = CliRunner().invoke(
                        taxaforge.cli,
                        [
                            "build",
                            "--db-type", "archaea",
                            "--db-name", "test",
                            "--cache-dir", str(root / "cache"),
                            "--output-dir", str(root / "output"),
                            "--kmer-len", "30",
                        ],
                    )
            finally:
                os.chdir(original_cwd)

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(build_db.call_args.args[5:7], (30, 30))
            self.assertIn("reducing --minimizer-len to 30", logs.output[0])

    def test_cli_reduces_default_minimizer_spaces_for_short_minimizer(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            original_cwd = os.getcwd()
            try:
                with patch.object(taxaforge, "run_basic_checks"), patch.object(
                    taxaforge, "download_taxanomy"
                ), patch.object(taxaforge, "download_genomes"), patch.object(
                    taxaforge, "add_to_library"
                ), patch.object(taxaforge, "build_db") as build_db, self.assertLogs(
                    taxaforge.logger, level="WARNING"
                ) as logs:
                    result = CliRunner().invoke(
                        taxaforge.cli,
                        [
                            "build",
                            "--db-type", "archaea",
                            "--db-name", "test",
                            "--cache-dir", str(root / "cache"),
                            "--output-dir", str(root / "output"),
                            "--kmer-len", "30",
                            "--minimizer-len", "26",
                        ],
                    )
            finally:
                os.chdir(original_cwd)

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(build_db.call_args.args[5:8], (30, 26, False))
            self.assertEqual(build_db.call_args.args[-1], 6)
            self.assertIn("at most 6 minimizer spaces", logs.output[0])

    def test_build_db_passes_minimizer_spaces_to_kraken2(self):
        commands = []
        with tempfile.TemporaryDirectory() as temp_dir, patch.object(
            taxaforge, "run_cmd", side_effect=commands.append
        ):
            taxaforge.build_db(
                "/cache",
                Path(temp_dir),
                "archaea",
                "test",
                4,
                30,
                26,
                False,
                False,
                0.7,
                False,
                minimizer_spaces=6,
                max_db_size="50G",
            )

        self.assertTrue(
            any(
                "kraken2-build --build" in command
                and "--minimizer-len 26 --minimizer-spaces 6" in command
                and "--max-db-size 50G" in command
                for command in commands
            ),
            commands,
        )

    def test_cli_passes_max_db_size_to_kraken2_build(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            original_cwd = os.getcwd()
            try:
                with patch.object(taxaforge, "run_basic_checks"), patch.object(
                    taxaforge, "download_taxanomy"
                ), patch.object(taxaforge, "download_genomes"), patch.object(
                    taxaforge, "add_to_library"
                ), patch.object(taxaforge, "build_db") as build_db:
                    result = CliRunner().invoke(
                        taxaforge.cli,
                        [
                            "build",
                            "--db-type", "archaea",
                            "--db-name", "test",
                            "--cache-dir", str(root / "cache"),
                            "--output-dir", str(root / "output"),
                            "--max-db-size", "50g",
                        ],
                    )
            finally:
                os.chdir(original_cwd)

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(build_db.call_args.kwargs["max_db_size"], "50G")

    def test_cli_rejects_invalid_max_db_size(self):
        result = CliRunner().invoke(
            taxaforge.cli,
            ["build", "--max-db-size", "50 GB"],
        )

        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("positive byte count", result.output)

    def test_cli_rejects_max_db_size_for_ganon2(self):
        result = CliRunner().invoke(
            taxaforge.cli,
            ["build", "--tool", "ganon2", "--max-db-size", "50G"],
        )

        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("only available with --tool kraken2", result.output)

    def test_cli_rejects_excessive_explicit_minimizer_spaces(self):
        result = CliRunner().invoke(
            taxaforge.cli,
            [
                "build",
                "--db-type", "archaea",
                "--kmer-len", "30",
                "--minimizer-len", "26",
                "--minimizer-spaces", "7",
            ],
        )

        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("exceeds the maximum 6", result.output)

    def test_cli_legacy_threads_option_sets_both_thread_counts(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            original_cwd = os.getcwd()
            try:
                with patch.object(taxaforge, "run_basic_checks"), patch.object(
                    taxaforge, "download_taxanomy"
                ), patch.object(taxaforge, "download_genomes") as download_genomes, patch.object(
                    taxaforge, "add_to_library"
                ) as add_to_library, patch.object(taxaforge, "build_db") as build_db:
                    result = CliRunner().invoke(
                        taxaforge.cli,
                        [
                            "build",
                            "--db-type", "archaea",
                            "--db-name", "test",
                            "--cache-dir", str(root / "cache"),
                            "--output-dir", str(root / "output"),
                            "--threads", "5",
                        ],
                    )
            finally:
                os.chdir(original_cwd)

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(download_genomes.call_args.args[4], 4)
            self.assertEqual(add_to_library.call_args.args[9], 5)
            self.assertEqual(build_db.call_args.args[4], 5)

    def test_cli_rejects_download_threads_above_safe_maximum(self):
        result = CliRunner().invoke(
            taxaforge.cli,
            ["build", "--download-threads", str(taxaforge.MAX_DOWNLOAD_THREADS + 1)],
        )

        self.assertNotEqual(result.exit_code, 0)
        self.assertIn(f"x<={taxaforge.MAX_DOWNLOAD_THREADS}", result.output)

    def test_createtaxdb_samplesheet_contains_every_accession_and_taxid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            accession_file = root / "assemblies.txt"
            accession_file.write_text(
                "GCF_000001985.1_genomic.fna\nGCA_000002435.2_genomic.fna\n",
                encoding="utf-8",
            )
            genome_dir = root / "genomes"
            genome_dir.mkdir()
            for name in (
                "GCF_000001985.1_genomic.fna.gz",
                "GCA_000002435.2_genomic.fna",
            ):
                (genome_dir / name).touch()

            refseq_summary = root / "assembly_summary_refseq.txt"
            refseq_summary.write_text(
                "# assembly_accession\ttaxid\nGCF_000001985.1\t562\n",
                encoding="utf-8",
            )
            genbank_summary = root / "assembly_summary_genbank.txt"
            genbank_summary.write_text(
                "# assembly_accession\ttaxid\nGCA_000002435.2\t1280\n",
                encoding="utf-8",
            )
            output = root / "samplesheet.csv"
            with patch.object(
                taxaforge,
                "ensure_assembly_summaries",
                return_value=[refseq_summary, genbank_summary],
            ) as ensure_summaries:
                taxaforge.create_createtaxdb_samplesheet(
                    accession_file, genome_dir, root / "cache", output
                )

            self.assertEqual(
                ensure_summaries.call_args.args[1],
                {"refseq", "genbank"},
            )
            with output.open(newline="", encoding="utf-8") as samplesheet:
                rows = list(csv.DictReader(samplesheet))
            self.assertEqual([row["id"] for row in rows], [
                "GCF_000001985.1",
                "GCA_000002435.2",
            ])
            self.assertEqual([row["taxid"] for row in rows], ["562", "1280"])
            self.assertTrue(Path(rows[0]["fasta_dna"]).is_file())
            self.assertTrue(Path(rows[1]["fasta_dna"]).is_file())
            self.assertEqual([row["fasta_aa"] for row in rows], ["", ""])

    def test_createtaxdb_samplesheet_fails_instead_of_skipping_missing_genomes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            accession_file = root / "assemblies.txt"
            accession_file.write_text("GCF_000001985.1_genomic.fna\n", encoding="utf-8")
            genome_dir = root / "genomes"
            genome_dir.mkdir()

            with self.assertRaisesRegex(ValueError, "Missing FASTA files"):
                taxaforge.create_createtaxdb_samplesheet(
                    accession_file, genome_dir, root / "cache", root / "samplesheet.csv"
                )

    def test_createtaxdb_samplesheet_fails_instead_of_skipping_missing_taxids(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            accession_file = root / "assemblies.txt"
            accession_file.write_text("GCF_000001985.1_genomic.fna\n", encoding="utf-8")
            genome_dir = root / "genomes"
            genome_dir.mkdir()
            (genome_dir / "GCF_000001985.1_genomic.fna.gz").touch()
            summary = root / "assembly_summary_refseq.txt"
            summary.write_text(
                "# assembly_accession\ttaxid\nGCF_000001985.1\tna\n",
                encoding="utf-8",
            )

            with patch.object(taxaforge, "ensure_assembly_summaries", return_value=[summary]):
                with self.assertRaisesRegex(ValueError, "Missing NCBI taxids"):
                    taxaforge.create_createtaxdb_samplesheet(
                        accession_file, genome_dir, root / "cache", root / "samplesheet.csv"
                    )

    def test_createtaxdb_runner_uses_default_options_unless_overridden(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            cache = root / "cache"
            taxonomy_dir = cache / "taxonomy"
            taxonomy_dir.mkdir(parents=True)
            (taxonomy_dir / "nodes.dmp").touch()
            (taxonomy_dir / "names.dmp").touch()
            accession_file = root / "assemblies.txt"
            accession_file.write_text("GCF_000001985.1_genomic.fna\n", encoding="utf-8")
            genome_dir = root / "genomes"
            genome_dir.mkdir()
            (genome_dir / "GCF_000001985.1_genomic.fna.gz").touch()
            summary = root / "assembly_summary_refseq.txt"
            summary.write_text(
                "# assembly_accession\ttaxid\nGCF_000001985.1\t562\n",
                encoding="utf-8",
            )
            observed = []

            def fake_run(command, **kwargs):
                observed.append((command, kwargs))

            with patch.object(
                taxaforge, "ensure_assembly_summaries", return_value=[summary]
            ), patch.object(taxaforge.subprocess, "run", side_effect=fake_run), patch.object(
                taxaforge.shutil, "which", return_value="/usr/bin/nextflow"
            ):
                taxaforge.run_createtaxdb(
                    "nextflow", "nf-core/createtaxdb", "1.0.0", root, "gtest",
                    cache, accession_file, genome_dir, None, None, None, None, None, True,
                )

            command, run_kwargs = observed[0]
            self.assertEqual(command[:5], [
                "nextflow", "run", "nf-core/createtaxdb", "-r", "1.0.0",
            ])
            self.assertIn("-profile", command)
            self.assertIn("docker", command)
            self.assertIn("-resume", command)
            self.assertFalse(any(arg.startswith("--ganon_build_options=") for arg in command))
            self.assertEqual(run_kwargs["cwd"], root)

            observed.clear()
            with patch.object(
                taxaforge, "ensure_assembly_summaries", return_value=[summary]
            ), patch.object(taxaforge.subprocess, "run", side_effect=fake_run), patch.object(
                taxaforge.shutil, "which", return_value="/usr/bin/nextflow"
            ):
                taxaforge.run_createtaxdb(
                    "nextflow", "nf-core/createtaxdb", None, root, "gtest",
                    cache, accession_file, genome_dir, None, 4, 31, 50,
                    "species", False,
                )
            command, _ = observed[0]
            self.assertIn(
                "--ganon_build_options=--threads 4 --kmer-size 31 "
                "--window-size 50 --level species",
                command,
            )

    def test_reads_accessions_from_fasta_filenames_and_applies_limit(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            accession_file = Path(temp_dir) / "assemblies.txt"
            accession_file.write_text(
                "\n".join(
                    [
                        "# selected genomes",
                        "GCF_000001985.1_JCVI-PMFA1-2.0_genomic.fna",
                        "GCF_000002415.2_ASM241v2_genomic.fna",
                        "GCA_000002435.2_ASM243v2_genomic.fna",
                    ]
                ),
                encoding="utf-8",
            )

            accessions = taxaforge.read_assembly_accessions(accession_file, limit=2)

        self.assertEqual(
            accessions,
            ["GCF_000001985.1", "GCF_000002415.2"],
        )

    def test_prefers_filename_primary_when_refseq_and_genbank_aliases_are_present(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            accession_file = Path(temp_dir) / "assemblies.txt"
            accession_file.write_text(
                "./GCF_000006355.2_GCA_000006355.2_genomic.fna\n",
                encoding="utf-8",
            )

            accessions = taxaforge.read_assembly_accessions(accession_file)

        self.assertEqual(accessions, ["GCF_000006355.2"])

    def test_downloads_only_listed_refseq_and_genbank_accessions(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            accession_file = Path(temp_dir) / "assemblies.txt"
            accession_file.write_text(
                "GCF_000001985.1_genomic.fna\nGCA_000002435.2_genomic.fna\n",
                encoding="utf-8",
            )

            def fake_download(**kwargs):
                for accession in kwargs["assembly_accessions"]:
                    genome = Path(kwargs["output"]) / f"{accession}_ASM_genomic.fna.gz"
                    genome.write_bytes(b"compressed fasta placeholder")
                return 0

            with patch.object(
                taxaforge.ncbi_genome_download, "download", side_effect=fake_download
            ) as download:
                genome_dir = taxaforge.download_assembly_accessions(
                    Path(temp_dir) / "cache", accession_file, threads=3
                )
                calls = [call.kwargs for call in download.call_args_list]
                downloaded_files = sorted(
                    path.name for path in Path(genome_dir).glob("*_genomic.fna.gz")
                )

            self.assertEqual([call["section"] for call in calls], ["refseq", "genbank"])
            self.assertEqual(
                [call["assembly_accessions"] for call in calls],
                [["GCF_000001985.1"], ["GCA_000002435.2"]],
            )
            self.assertTrue(all(call["flat_output"] for call in calls))
            self.assertEqual(
                downloaded_files,
                [
                    "GCA_000002435.2_ASM_genomic.fna.gz",
                    "GCF_000001985.1_ASM_genomic.fna.gz",
                ],
            )

    def test_reuses_genome_cache_without_downloading_cached_accessions(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            accession_file = root / "assemblies.txt"
            accession_file.write_text("GCF_000001985.1_genomic.fna\n", encoding="utf-8")
            cached_genome = root / "genomes-cache" / "refseq" / "GCF_000001985.1"
            cached_genome.mkdir(parents=True)
            cached_file = cached_genome / "GCF_000001985.1_ASM_genomic.fna.gz"
            cached_file.write_bytes(b"cached genome")

            with patch.object(
                taxaforge.ncbi_genome_download, "download"
            ) as download:
                genome_dir = taxaforge.download_assembly_accessions(
                    root / "build-cache", accession_file, threads=1,
                    genomes_cache_dir=root / "genomes-cache",
                )

            download.assert_not_called()
            staged_file = Path(genome_dir) / cached_file.name
            self.assertTrue(staged_file.is_symlink())
            self.assertEqual(staged_file.resolve(), cached_file.resolve())

    def test_downloads_only_accessions_missing_from_genome_cache(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            accession_file = root / "assemblies.txt"
            accession_file.write_text(
                "GCF_000001985.1_genomic.fna\nGCA_000002435.2_genomic.fna\n",
                encoding="utf-8",
            )
            cached_dir = root / "genomes-cache"
            cached_dir.mkdir()
            cached_file = cached_dir / "GCF_000001985.1_ASM_genomic.fna"
            cached_file.write_bytes(b"cached genome")

            def fake_download(**kwargs):
                self.assertEqual(kwargs["assembly_accessions"], ["GCA_000002435.2"])
                (Path(kwargs["output"]) / "GCA_000002435.2_ASM_genomic.fna.gz").write_bytes(
                    b"downloaded genome"
                )
                return 0

            with patch.object(
                taxaforge.ncbi_genome_download, "download", side_effect=fake_download
            ) as download:
                genome_dir = taxaforge.download_assembly_accessions(
                    root / "build-cache", accession_file, threads=1,
                    genomes_cache_dir=cached_dir,
                )

            self.assertEqual(download.call_count, 1)
            self.assertTrue(
                (Path(genome_dir) / cached_file.name).is_symlink()
            )
            self.assertTrue(
                (Path(genome_dir) / "GCA_000002435.2_ASM_genomic.fna.gz").is_file()
            )

    def test_fails_if_a_requested_accession_was_not_downloaded(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            accession_file = Path(temp_dir) / "assemblies.txt"
            accession_file.write_text("GCF_000001985.1_genomic.fna\n", encoding="utf-8")

            with patch.object(taxaforge.ncbi_genome_download, "download", return_value=1), \
                    patch.object(taxaforge, "ensure_assembly_summaries", return_value=[]):
                with self.assertRaisesRegex(RuntimeError, "GCF_000001985.1"):
                    taxaforge.download_assembly_accessions(
                        Path(temp_dir) / "cache", accession_file, threads=1
                    )

    def test_downloads_historical_assemblies_from_summary_paths(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            accession_file = temp_path / "assemblies.txt"
            accession_file.write_text("GCF_000001985.1_genomic.fna\n", encoding="utf-8")
            summary_file = temp_path / "assembly_summary.txt"
            summary_file.write_text(
                "#assembly_accession\tftp_path\n"
                "GCF_000001985.1\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/985/GCF_000001985.1_JCVI-PMFA1-2.0/\n",
                encoding="utf-8",
            )

            def fake_download_files(urls, max_workers):
                for url in urls:
                    (Path.cwd() / url.rsplit("/", 1)[-1]).write_bytes(b"genome")

            with patch.object(taxaforge.ncbi_genome_download, "download", return_value=1), \
                    patch.object(taxaforge, "ensure_assembly_summaries", return_value=[summary_file]), \
                    patch.object(taxaforge, "download_files", side_effect=fake_download_files):
                genome_dir = taxaforge.download_assembly_accessions(
                    temp_path / "cache", accession_file, threads=1
                )

            self.assertTrue(
                (Path(genome_dir) / "GCF_000001985.1_JCVI-PMFA1-2.0_genomic.fna.gz").is_file()
            )

    def test_rejects_lines_without_exactly_one_accession(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            accession_file = Path(temp_dir) / "assemblies.txt"
            accession_file.write_text("not-an-accession\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "line 1.*not-an-accession"):
                taxaforge.read_assembly_accessions(accession_file)

    def test_detects_only_needed_assembly_summary_sections(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            refseq_dir = root / "refseq"
            genbank_dir = root / "genbank"
            refseq_dir.mkdir()
            genbank_dir.mkdir()
            (refseq_dir / "GCF_000001985.1_genomic.fna.gz").touch()
            (genbank_dir / "GCA_000002435.2_genomic.fna.gz").touch()

            self.assertEqual(
                taxaforge.detect_assembly_summary_sections(
                    [refseq_dir, genbank_dir], "fna.gz"
                ),
                {"refseq", "genbank"},
            )
            self.assertEqual(
                taxaforge.detect_assembly_summary_sections([refseq_dir], "fna.gz"),
                {"refseq"},
            )

    def test_assembly_summary_download_preserves_cwd_and_selects_section(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            original_cwd = Path.cwd()
            downloaded_urls = []

            def fake_download_files(urls):
                downloaded_urls.extend(urls)
                for url in urls:
                    (Path.cwd() / url.rsplit("/", 1)[-1]).touch()

            with patch.object(taxaforge, "download_files", side_effect=fake_download_files):
                paths = taxaforge.ensure_assembly_summaries(
                    temp_dir, sections={"refseq"}
                )

            self.assertEqual(Path.cwd(), original_cwd)
            self.assertEqual(len(paths), 2)
            self.assertTrue(all("/refseq/" in url for url in downloaded_urls))

    def test_ganon_build_uses_output_directory_and_only_needed_summaries(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            output_dir = root / "output"
            genomes_dir = root / "genomes"
            output_dir.mkdir()
            genomes_dir.mkdir()
            (genomes_dir / "GCF_000001985.1_genomic.fna.gz").touch()
            observed = []

            def fake_run_cmd(command):
                observed.append((command, Path.cwd()))

            def fake_summaries(_cache, sections):
                self.assertEqual(sections, {"refseq"})
                return ["summary.txt"]

            with patch.object(taxaforge, "detect_taxonomy_flag", return_value="--taxonomy-files"), \
                    patch.object(
                        taxaforge, "ensure_assembly_summaries", side_effect=fake_summaries
                    ), \
                    patch.object(taxaforge, "ensure_genome_size_file", return_value="genome-size.txt"), \
                    patch.object(taxaforge, "setup_raptor_shim", side_effect=lambda path: path), \
                    patch.object(taxaforge, "run_cmd", side_effect=fake_run_cmd):
                original_cwd = Path.cwd()
                try:
                    os.chdir(output_dir)
                    taxaforge.run_ganon_build_custom(
                        str(root / "cache"), [str(genomes_dir)], "g1000", 4, 31, 50,
                        "file", restart=True,
                    )
                    first_command = observed[0][0]
                    first_cwd = observed[0][1]
                    observed.clear()
                    taxaforge.run_ganon_build_custom(
                        str(root / "cache"), [str(genomes_dir)], "g1000", 4,
                        None, None, None, restart=True,
                    )
                    default_command = observed[0][0]
                finally:
                    os.chdir(original_cwd)

            self.assertTrue(first_command.startswith("ganon build-custom"))
            self.assertTrue(first_command.endswith("--restart"))
            self.assertEqual(first_cwd, output_dir.resolve())
            self.assertNotIn("--kmer-size", default_command)
            self.assertNotIn("--window-size", default_command)
            self.assertNotIn("--level", default_command)
            self.assertTrue(default_command.endswith("--restart"))

    def test_accession_list_forwards_ganon_options_to_custom_build(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            accession_file = root / "assemblies.txt"
            accession_file.write_text("GCF_000001985.1\n", encoding="utf-8")
            genomes_dir = root / "genomes"
            genomes_dir.mkdir()

            with patch.object(taxaforge, "run_basic_checks"), \
                    patch.object(taxaforge, "download_taxanomy"), \
                    patch.object(
                        taxaforge, "download_assembly_accessions",
                        return_value=str(genomes_dir),
                    ), \
                    patch.object(taxaforge, "build_ganon2") as build_ganon2:
                result = CliRunner().invoke(
                    taxaforge.cli,
                    [
                        "build", "--tool", "ganon2",
                        "--assembly-accessions-file", str(accession_file),
                        "--db-name", "gtest", "--max-fp", "0.01",
                    ],
                )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(
                build_ganon2.call_args.kwargs["custom_args"],
                ["--max-fp", "0.01"],
            )

    def test_custom_build_appends_ganon_options(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            genomes_dir = root / "genomes"
            genomes_dir.mkdir()
            (genomes_dir / "GCF_000001985.1_genomic.fna.gz").touch()
            observed = []

            with patch.object(
                taxaforge, "detect_taxonomy_flag", return_value="--taxonomy-files"
            ), patch.object(
                taxaforge, "ensure_assembly_summaries", return_value=["summary.txt"]
            ), patch.object(
                taxaforge, "ensure_genome_size_file", return_value="genome-size.txt"
            ), patch.object(
                taxaforge, "setup_raptor_shim", side_effect=lambda path: path
            ), patch.object(
                taxaforge, "run_cmd", side_effect=observed.append
            ):
                taxaforge.run_ganon_build_custom(
                    str(root / "cache"), [str(genomes_dir)], "gtest", 1,
                    None, None, None, custom_args=["--max-fp", "0.01"],
                )

            build_command = observed[0]
            self.assertTrue(build_command.startswith("ganon build-custom"))
            self.assertIn("--max-fp 0.01", build_command)
            self.assertNotIn("ganon build --max-fp", build_command)


if __name__ == "__main__":
    unittest.main()
