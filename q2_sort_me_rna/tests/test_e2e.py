# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from collections import Counter
import gzip
import importlib.util
import re
import shutil
import subprocess
import unittest
from pathlib import Path

import pandas as pd
import pandas.testing as pdt

from qiime2 import Artifact
from qiime2.plugins.sort_me_rna.pipelines import (
    align_sequences,
    denovo_otu_mapping,
    otu_mapping,
)
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import BLAST6Format
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    CasavaOneEightSingleLanePerSampleDirFmt,
)


@unittest.skipUnless(
    shutil.which("sortmerna")
    and importlib.util.find_spec("q2_fondue")
    and importlib.util.find_spec("q2_feature_table"),
    "SortMeRNA, q2-fondue, and q2-feature-table are required.",
)
class SortMeRNAEndToEndTests(TestPluginBase):
    package = "q2_sort_me_rna.tests"

    def setUp(self):
        super().setUp()
        self.references = Artifact.load(self.get_data_path("rrna_references.qza"))
        self.single_reads = Artifact.load(self.get_data_path("raw_sequence.qza"))
        self.paired_reads = Artifact.load(self.get_data_path("paired_raw_sequence.qza"))
        self.expected_aligned_sequences = {
            "Seq1": "AGCTGATCGATCGATCGATCGATCGATCGATC",
            "Seq2": "TACGTACGTACGTACGTACGTACGTACGTACG",
            "Seq3Perturbed2": "TTTTTACGTACGTACGTACGTACGTACGTACG",
            "Seq5CopyOfTwo": "TACGTACGTACGTACGTACGTACGTACGTACG",
        }
        self.expected_reference_assignments = {
            "Seq1": "16S_rRNA_1",
            "Seq2": "16S_rRNA_2",
            "Seq3Perturbed2": "16S_rRNA_2",
            "Seq5CopyOfTwo": "16S_rRNA_2",
        }

    def test_align_sequences(self):
        results = align_sequences(
            references=self.references,
            reads=self.paired_reads,
            num_alignments=1,
            no_best=True,
            score_split=True,
            max_read_len=30000,
        )

        self._assert_artifact(results.blast_aligned_seq, "FeatureData[BLAST6]")
        self._assert_artifact(
            results.fastx_aligned_seq,
            "SampleData[SequencesWithQuality]",
        )
        self._assert_artifact(results.alignment_map, "SampleData[AlignmentMap]")

        expected_sequences = {
            f"{read_id}\\{mate}": sequence
            for read_id, sequence in self.expected_aligned_sequences.items()
            for mate in (1, 2)
        }
        expected_assignments = Counter(
            (f"{read_id}\\{mate}", reference)
            for read_id, reference in self.expected_reference_assignments.items()
            for mate in (1, 2)
        )
        self.assertEqual(
            self._read_fastq_records(results.fastx_aligned_seq),
            {"raw_0_L001_R1_001.fastq.gz": expected_sequences},
        )
        self.assertEqual(
            self._read_blast_assignments(results.blast_aligned_seq),
            expected_assignments,
        )
        self.assertEqual(
            self._read_bam_assignments(results.alignment_map),
            {"raw.bam": expected_assignments},
        )

    def test_sortmerna_supported_version(self):
        completed = subprocess.run(
            ["sortmerna", "--version"],
            check=True,
            capture_output=True,
            text=True,
        )
        version_output = completed.stdout + completed.stderr
        match = re.search(r"SortMeRNA version (\d+)\.(\d+)\.(\d+)", version_output)
        self.assertIsNotNone(match, version_output)
        version = tuple(map(int, match.groups()))
        self.assertGreaterEqual(version, (8, 0, 1))
        self.assertLess(version, (9, 0, 0))

    def test_otu_mapping(self):
        results = otu_mapping(
            references=self.references,
            reads=self.single_reads,
            id=0.12,
            coverage=0.12,
        )

        self._assert_otu_table(results.otu_mapping, ["raw"])

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
        self._assert_otu_table(results.otu_mapping, ["raw"])
        self.assertEqual(
            self._read_fastq_records(results.denovo_aligned_seq),
            {
                "raw_0_L001_R1_001.fastq.gz": {
                    "Seq3Perturbed2": self.expected_aligned_sequences["Seq3Perturbed2"]
                }
            },
        )

    def test_parallel_align_sequences_combines_samples(self):
        source_dir = self.single_reads.view(CasavaOneEightSingleLanePerSampleDirFmt)
        source = next(Path(str(source_dir)).glob("*.fastq.gz"))
        multi = CasavaOneEightSingleLanePerSampleDirFmt()
        for sample_id in ("raw", "second"):
            destination = Path(str(multi)) / f"{sample_id}_S1_L001_R1_001.fastq.gz"
            shutil.copyfile(source, destination)
        reads = Artifact.import_data("SampleData[SequencesWithQuality]", multi)

        with self.test_config:
            results = align_sequences.parallel(self.references, reads)._result()

        self._assert_artifact(results.blast_aligned_seq, "FeatureData[BLAST6]")
        self._assert_artifact(
            results.fastx_aligned_seq,
            "SampleData[SequencesWithQuality]",
        )
        self._assert_artifact(results.alignment_map, "SampleData[AlignmentMap]")

        expected_assignments = Counter(self.expected_reference_assignments.items())
        observed_bams = self._read_bam_assignments(results.alignment_map)
        self.assertEqual(
            set(observed_bams),
            {"raw.bam", "second.bam"},
        )
        for assignments in observed_bams.values():
            self.assertEqual(assignments, expected_assignments)

        observed_fastqs = self._read_fastq_records(results.fastx_aligned_seq)
        self.assertEqual(
            set(observed_fastqs),
            {
                "raw_0_L001_R1_001.fastq.gz",
                "second_0_L001_R1_001.fastq.gz",
            },
        )
        for records in observed_fastqs.values():
            self.assertEqual(records, self.expected_aligned_sequences)

        self.assertEqual(
            self._read_blast_assignments(results.blast_aligned_seq),
            expected_assignments + expected_assignments,
        )

    def test_parallel_otu_mapping(self):
        with self.test_config:
            results = otu_mapping.parallel(
                self.references, self.single_reads, id=0.12, coverage=0.12
            )._result()
        self._assert_otu_table(results.otu_mapping, ["raw"])

    def test_parallel_denovo_otu_mapping(self):
        with self.test_config:
            results = denovo_otu_mapping.parallel(
                self.references, self.single_reads, id=0.7, coverage=0.7
            )._result()
        self._assert_artifact(
            results.denovo_aligned_seq, "SampleData[SequencesWithQuality]"
        )
        self._assert_otu_table(results.otu_mapping, ["raw"])
        self.assertEqual(
            self._read_fastq_records(results.denovo_aligned_seq),
            {
                "raw_0_L001_R1_001.fastq.gz": {
                    "Seq3Perturbed2": self.expected_aligned_sequences["Seq3Perturbed2"]
                }
            },
        )

    def _assert_artifact(self, artifact, expected_type):
        self.assertIsInstance(artifact, Artifact)
        self.assertEqual(str(artifact.type), expected_type)
        artifact.validate(level="max")

    def _assert_otu_table(self, artifact, sample_ids):
        self._assert_artifact(artifact, "FeatureTable[Frequency]")
        observed = artifact.view(pd.DataFrame).sort_index().sort_index(axis=1)
        expected = pd.DataFrame(
            {
                "16S_rRNA_1": [1.0] * len(sample_ids),
                "16S_rRNA_2": [2.0] * len(sample_ids),
            },
            index=sample_ids,
        ).sort_index()
        pdt.assert_frame_equal(observed, expected, check_names=False)

    def _read_blast_assignments(self, artifact):
        blast = artifact.view(BLAST6Format)
        assignments = Counter()
        for line in Path(str(blast)).read_text().splitlines():
            fields = line.split("\t")
            self.assertEqual(len(fields), 12)
            assignments[(fields[0], fields[1])] += 1
        return assignments

    def _read_fastq_records(self, artifact):
        fastqs = artifact.view(CasavaOneEightSingleLanePerSampleDirFmt)
        observed = {}
        for path in sorted(Path(str(fastqs)).glob("*.fastq.gz")):
            with gzip.open(path, "rt") as handle:
                lines = handle.read().splitlines()
            self.assertEqual(len(lines) % 4, 0)
            records = {}
            for index in range(0, len(lines), 4):
                name, sequence, separator, quality = lines[index : index + 4]
                self.assertTrue(name.startswith("@"))
                self.assertEqual(separator, "+")
                self.assertEqual(len(sequence), len(quality))
                self.assertNotIn(name[1:], records)
                records[name[1:]] = sequence
            observed[path.name] = records
        return observed

    def _read_bam_assignments(self, artifact):
        bams = artifact.view(BAMDirFmt)
        observed = {}
        for path in sorted(Path(str(bams)).glob("*.bam")):
            completed = subprocess.run(
                ["samtools", "view", str(path)],
                check=True,
                capture_output=True,
                text=True,
            )
            assignments = Counter()
            for line in completed.stdout.splitlines():
                fields = line.split("\t")
                self.assertGreaterEqual(len(fields), 11)
                assignments[(fields[0], fields[2])] += 1
            observed[path.name] = assignments
        return observed
