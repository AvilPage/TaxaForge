import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import taxaforge


class AssemblyAccessionDownloadTests(unittest.TestCase):
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
                "GCF_000006355.2_GCA_000006355.2_genomic.fna\n",
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

            with self.assertRaisesRegex(ValueError, "line 1"):
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
                        str(root / "cache"), [str(genomes_dir)], "g1000", 4, 31, 50, "file"
                    )
                finally:
                    os.chdir(original_cwd)

            self.assertTrue(observed[0][0].startswith("ganon build-custom"))
            self.assertEqual(observed[0][1], output_dir.resolve())


if __name__ == "__main__":
    unittest.main()
