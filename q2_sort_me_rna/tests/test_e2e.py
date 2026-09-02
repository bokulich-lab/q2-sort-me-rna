# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import shutil
import subprocess
import unittest

from qiime2 import Artifact
from qiime2.plugins.sort_me_rna.methods import (
    align_sequences,
    denovo_otu_mapping,
    otu_mapping,
)
from qiime2.plugin.testing import TestPluginBase


@unittest.skipUnless(
    shutil.which("sortmerna"),
    "SortMeRNA is not installed; skipping end-to-end tests.",
)
class SortMeRNAEndToEndTests(TestPluginBase):
    package = "q2_sort_me_rna.tests"

    def setUp(self):
        super().setUp()
        self.references = Artifact.load(
            self.get_data_path("rrna_references.qza")
        )
        self.single_reads = Artifact.load(
            self.get_data_path("raw_sequence.qza")
        )
        self.paired_reads = Artifact.load(
            self.get_data_path("paired_raw_sequence.qza")
        )

    def test_align_sequences(self):
        results = align_sequences(
            references=self.references,
            reads=self.paired_reads,
            num_alignments=1,
            no_best=True,
            score_split=True,
            max_read_len=30000,
        )

        self._assert_artifact(
            results.blast_aligned_seq, "FeatureData[BLAST6]"
        )
        self._assert_artifact(
            results.fastx_aligned_seq,
            "SampleData[SequencesWithQuality]",
        )
        self._assert_artifact(
            results.alignment_map, "SampleData[AlignmentMap]"
        )

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
        results = otu_mapping(
            references=self.references,
            reads=self.single_reads,
            id=0.12,
            coverage=0.12,
        )

        self._assert_artifact(
            results.otu_mapping, "FeatureTable[Frequency]"
        )

    def test_denovo_otu_mapping(self):
        results = denovo_otu_mapping(
            references=self.references,
            reads=self.single_reads,
            id=0.7,
            coverage=0.7,
        )

        self._assert_artifact(
            results.denovo_aligned_seq,
            "SampleData[SequencesWithQuality]",
        )

    def _assert_artifact(self, artifact, expected_type):
        self.assertIsInstance(artifact, Artifact)
        self.assertEqual(str(artifact.type), expected_type)
