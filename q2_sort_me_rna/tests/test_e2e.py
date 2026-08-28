# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase


@unittest.skipUnless(
    shutil.which("sortmerna"),
    "SortMeRNA is not installed; skipping end-to-end tests.",
)
class SortMeRNAEndToEndTests(TestPluginBase):
    package = "q2_sort_me_rna.tests"

    def setUp(self):
        super().setUp()
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.workdir = Path(self.temporary_directory.name)
        self.output_dir = self.workdir / "qiime-output"

    def tearDown(self):
        self.temporary_directory.cleanup()
        super().tearDown()

    def test_align_sequences(self):
        command = [
            "qiime",
            "sort-me-rna",
            "align-sequences",
            "--i-references",
            self.get_data_path("rrna_references.qza"),
            "--i-reads",
            self.get_data_path("paired_raw_sequence.qza"),
            "--p-num-alignments",
            "1",
            "--p-no-best",
            "--p-score-split",
            "--p-max-read-len",
            "30000",
            "--output-dir",
            str(self.output_dir),
        ]
        subprocess.run(command, check=True)

        self._assert_artifact("blast_aligned_seq.qza", "FeatureData[BLAST6]")
        self._assert_artifact(
            "fastx_aligned_seq.qza", "SampleData[SequencesWithQuality]"
        )
        self._assert_artifact("alignment_map.qza", "SampleData[AlignmentMap]")

    def test_sortmerna_major_version(self):
        completed = subprocess.run(
            ["sortmerna", "--version"],
            check=True,
            capture_output=True,
            text=True,
        )
        version_output = completed.stdout + completed.stderr
        self.assertRegex(version_output, r"SortMeRNA version 7\.")

    def test_otu_mapping(self):
        command = [
            "qiime",
            "sort-me-rna",
            "otu-mapping",
            "--i-references",
            self.get_data_path("rrna_references.qza"),
            "--i-reads",
            self.get_data_path("raw_sequence.qza"),
            "--p-id",
            "0.12",
            "--p-coverage",
            "0.12",
            "--output-dir",
            str(self.output_dir),
        ]
        subprocess.run(command, check=True)

        self._assert_artifact("otu_mapping.qza", "FeatureTable[Frequency]")

    def test_denovo_otu_mapping(self):
        command = [
            "qiime",
            "sort-me-rna",
            "denovo-otu-mapping",
            "--i-references",
            self.get_data_path("rrna_references.qza"),
            "--i-reads",
            self.get_data_path("raw_sequence.qza"),
            "--p-id",
            "0.7",
            "--p-coverage",
            "0.7",
            "--output-dir",
            str(self.output_dir),
        ]
        subprocess.run(command, check=True)

        self._assert_artifact(
            "denovo_aligned_seq.qza",
            "SampleData[SequencesWithQuality]",
        )

    def _assert_artifact(self, filename, expected_type):
        artifact_path = self.output_dir / filename
        self.assertTrue(artifact_path.exists())
        artifact = Artifact.load(artifact_path)
        self.assertEqual(str(artifact.type), expected_type)
