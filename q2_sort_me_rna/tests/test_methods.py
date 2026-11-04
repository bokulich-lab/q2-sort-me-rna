# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import base64
import inspect
from pathlib import Path
import shutil
import subprocess
from unittest.mock import patch

import pandas as pd
import pandas.testing as pdt

from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import BLAST6Format, DNAFASTAFormat
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    CasavaOneEightSingleLanePerSampleDirFmt,
)

from q2_sort_me_rna._methods import (
    align_sequences,
    denovo_otu_mapping,
    otu_mapping,
)
from q2_sort_me_rna._rna_sorter import _get_read_sets, _parse_parameters


class SortMeRNAMethodTests(TestPluginBase):
    package = "q2_sort_me_rna.tests"

    def setUp(self):
        super().setUp()
        self.reference = Artifact.load(self.get_data_path("rrna_references.qza")).view(
            DNAFASTAFormat
        )
        self.reference_artifact = Artifact.load(
            self.get_data_path("rrna_references.qza")
        )
        self.single_reads = Artifact.load(self.get_data_path("raw_sequence.qza")).view(
            CasavaOneEightSingleLanePerSampleDirFmt
        )
        self.single_reads_artifact = Artifact.load(
            self.get_data_path("raw_sequence.qza")
        )
        self.multi_reads = CasavaOneEightSingleLanePerSampleDirFmt()
        source_read = next(Path(str(self.single_reads)).glob("*.fastq.gz"))
        for sample_id in ("raw", "second"):
            destination = (
                Path(str(self.multi_reads)) / f"{sample_id}_S1_L001_R1_001.fastq.gz"
            )
            shutil.copyfile(source_read, destination)
        self.paired_reads = Artifact.load(
            self.get_data_path("paired_raw_sequence.qza")
        ).view(CasavaOneEightSingleLanePerSampleDirFmt)
        self.multi_paired_reads = CasavaOneEightSingleLanePerSampleDirFmt()
        for source_read in Path(str(self.paired_reads)).glob("*.fastq.gz"):
            direction = "R1" if "_R1_" in source_read.name else "R2"
            for sample_id in ("raw", "second"):
                destination = (
                    Path(str(self.multi_paired_reads))
                    / f"{sample_id}_S1_L001_{direction}_001.fastq.gz"
                )
                shutil.copyfile(source_read, destination)
        self.mock_output_dir = Path(self.get_data_path("sortmerna-output"))

    def _mock_sortmerna(self, command, **kwargs):
        if command[0] == "samtools":
            destination = Path(command[command.index("-o") + 1])
            encoded_bam = (self.mock_output_dir / "aligned.bam.b64").read_bytes()
            destination.write_bytes(base64.b64decode(encoded_bam))
            return subprocess.CompletedProcess(command, returncode=0)

        workdir = Path(command[command.index("--workdir") + 1])
        output_dir = workdir / "out"
        output_dir.mkdir(parents=True)

        for filename in ("aligned.blast", "aligned.fastq", "aligned.sam"):
            shutil.copyfile(self.mock_output_dir / filename, output_dir / filename)

        if "--otu_map" in command:
            shutil.copyfile(
                self.mock_output_dir / "otu_map.txt",
                output_dir / "otu_map.txt",
            )
        if "--de_novo_otu" in command:
            shutil.copyfile(
                self.mock_output_dir / "aligned_denovo.fastq",
                output_dir / "aligned_denovo.fastq",
            )

        return subprocess.CompletedProcess(command, returncode=0)

    def test_align_sequences_success(self):
        with patch(
            "q2_sort_me_rna._rna_sorter.run_command",
            side_effect=self._mock_sortmerna,
        ) as mock_run:
            observed = align_sequences(
                self.reference,
                self.single_reads,
                num_alignments=2,
                no_best=True,
                score_split=True,
                max_read_len=40000,
                threads=2,
            )

        self.assertIsInstance(observed[0], BLAST6Format)
        self.assertIsInstance(observed[1], CasavaOneEightSingleLanePerSampleDirFmt)
        self.assertIsInstance(observed[2], BAMDirFmt)
        observed[2].validate(level="max")

        commands = [call.args[0] for call in mock_run.call_args_list]
        command = next(command for command in commands if command[0] == "sortmerna")
        self.assertEqual(command[0], "sortmerna")
        self.assertIn("--ref", command)
        self.assertIn("--reads", command)
        self.assertIn("--fastx", command)
        self.assertIn("--sam", command)
        self.assertIn("--SQ", command)
        self.assertIn("--no-best", command)
        self.assertIn("--score_split", command)
        self.assertNotIn("True", command)
        self.assertNotIn("--best", command)
        self.assertEqual(command[command.index("--num_alignments") + 1], "2")
        self.assertEqual(command[command.index("--max_read_len") + 1], "40000")
        self.assertEqual(command[command.index("--threads") + 1], "2")
        temporary_workdir = Path(command[command.index("--workdir") + 1])
        self.assertFalse(temporary_workdir.exists())

    def test_otu_mapping_success(self):
        with patch(
            "q2_sort_me_rna._rna_sorter.run_command",
            side_effect=self._mock_sortmerna,
        ):
            observed = otu_mapping(self.reference, self.single_reads)

        expected = pd.DataFrame(
            {"ref1": [2]}, index=pd.Index(["raw"], name="sample-id")
        )
        pdt.assert_frame_equal(observed[3], expected)

    def test_multi_sample_outputs_are_merged(self):
        with patch(
            "q2_sort_me_rna._rna_sorter.run_command",
            side_effect=self._mock_sortmerna,
        ) as mock_run:
            observed = otu_mapping(self.reference, self.multi_reads)

        sortmerna_commands = [
            call.args[0]
            for call in mock_run.call_args_list
            if call.args[0][0] == "sortmerna"
        ]
        self.assertEqual(len(sortmerna_commands), 2)

        blast_lines = Path(str(observed[0])).read_text().splitlines()
        self.assertEqual(len(blast_lines), 2)
        self.assertEqual(
            {path.name for path in Path(str(observed[1])).iterdir()},
            {
                "raw_0_L001_R1_001.fastq.gz",
                "second_0_L001_R1_001.fastq.gz",
            },
        )
        self.assertEqual(
            {path.name for path in Path(str(observed[2])).iterdir()},
            {"raw.bam", "second.bam"},
        )
        expected = pd.DataFrame(
            {"ref1": [2, 2]},
            index=pd.Index(["raw", "second"], name="sample-id"),
        ).sort_index()
        pdt.assert_frame_equal(observed[3].sort_index(), expected)

        observed[0].validate("max")
        observed[1].validate(level="max")
        observed[2].validate(level="max")

    def test_multi_sample_denovo_outputs_are_merged(self):
        with patch(
            "q2_sort_me_rna._rna_sorter.run_command",
            side_effect=self._mock_sortmerna,
        ):
            observed = denovo_otu_mapping(self.reference, self.multi_reads)

        self.assertEqual(
            {path.name for path in Path(str(observed[4])).iterdir()},
            {
                "raw_0_L001_R1_001.fastq.gz",
                "second_0_L001_R1_001.fastq.gz",
            },
        )
        observed[4].validate(level="max")

    def test_sample_without_otu_assignments_is_retained(self):
        sortmerna_calls = 0

        def mock_one_empty_otu_map(command, **kwargs):
            nonlocal sortmerna_calls
            result = self._mock_sortmerna(command, **kwargs)
            if command[0] == "sortmerna":
                sortmerna_calls += 1
                if sortmerna_calls == 2:
                    workdir = Path(command[command.index("--workdir") + 1])
                    (workdir / "out" / "otu_map.txt").write_text("")
            return result

        with patch(
            "q2_sort_me_rna._rna_sorter.run_command",
            side_effect=mock_one_empty_otu_map,
        ):
            observed = otu_mapping(self.reference, self.multi_reads)

        expected = pd.DataFrame(
            {"ref1": [2, 0]},
            index=pd.Index(["raw", "second"], name="sample-id"),
        ).sort_index()
        pdt.assert_frame_equal(observed[3].sort_index(), expected)

    def test_denovo_otu_mapping_success(self):
        with patch(
            "q2_sort_me_rna._rna_sorter.run_command",
            side_effect=self._mock_sortmerna,
        ):
            observed = denovo_otu_mapping(self.reference, self.single_reads)

        self.assertEqual(len(observed), 5)
        self.assertIsInstance(observed[4], CasavaOneEightSingleLanePerSampleDirFmt)

    def test_sortmerna_failure_is_actionable(self):
        error = subprocess.CalledProcessError(
            returncode=2,
            cmd=["sortmerna"],
            stderr="invalid option",
        )
        with patch(
            "q2_sort_me_rna._rna_sorter.run_command",
            side_effect=error,
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "SortMeRNA failed with exit code 2 for sample " "'raw': invalid option",
            ):
                align_sequences(self.reference, self.single_reads)

    def test_samtools_failure_is_actionable(self):
        def fail_samtools(command, **kwargs):
            if command[0] == "samtools":
                raise subprocess.CalledProcessError(
                    returncode=2,
                    cmd=command,
                    stderr="invalid SAM",
                )
            return self._mock_sortmerna(command, **kwargs)

        with patch(
            "q2_sort_me_rna._rna_sorter.run_command",
            side_effect=fail_samtools,
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "samtools failed with exit code 2: invalid SAM",
            ):
                align_sequences(self.reference, self.single_reads)

    def test_parameter_parsing_returns_safe_argument_list(self):
        observed = _parse_parameters(
            {
                "ref": "/tmp/reference with spaces.fasta",
                "fastx": True,
                "threads": 0,
                "other": False,
            }
        )
        self.assertEqual(
            observed,
            [
                "--ref",
                "/tmp/reference with spaces.fasta",
                "--fastx",
                "--threads",
                "0",
            ],
        )

    def test_parameter_parsing_uses_sortmerna_v7_option_names(self):
        observed = _parse_parameters(
            {
                "no_best": True,
                "score_split": True,
                "max_read_len": 1000,
                "f": True,
                "n": -4,
                "r": True,
            }
        )
        self.assertEqual(
            observed,
            [
                "--no-best",
                "--score_split",
                "--max_read_len",
                "1000",
                "-F",
                "-N",
                "-4",
                "-R",
            ],
        )

    def test_parameter_parsing_resolves_default_alignment_limit(self):
        observed = _parse_parameters(
            {
                "num_alignments": 1,
                "min_lis": 2,
                "match": 2,
                "threads": 1,
            }
        )
        self.assertEqual(observed, ["--match", "2", "--threads", "1"])

    def test_parameter_parsing_rejects_two_custom_alignment_limits(self):
        with self.assertRaisesRegex(
            ValueError,
            "SortMeRNA accepts either num_alignments or min_lis, not both",
        ):
            _parse_parameters({"num_alignments": 2, "min_lis": 3})

    def test_paired_reads_are_ordered(self):
        read_sets = _get_read_sets(self.paired_reads)
        self.assertEqual(len(read_sets), 1)
        forward, reverse, sample_id = read_sets[0]
        self.assertIn("_R1_", forward.name)
        self.assertIn("_R2_", reverse.name)
        self.assertEqual(sample_id, "raw")

    def test_multiple_samples_are_grouped(self):
        read_sets = _get_read_sets(self.multi_reads)
        self.assertEqual(
            [sample_id for _, _, sample_id in read_sets],
            ["raw", "second"],
        )
        self.assertTrue(all(reverse is None for _, reverse, _ in read_sets))

    def test_multiple_paired_samples_are_grouped(self):
        read_sets = _get_read_sets(self.multi_paired_reads)
        self.assertEqual(
            [sample_id for _, _, sample_id in read_sets],
            ["raw", "second"],
        )
        for forward, reverse, _ in read_sets:
            self.assertIn("_R1_", forward.name)
            self.assertIn("_R2_", reverse.name)

    def test_actions_are_registered(self):
        self.assertIn("fetch_db", self.plugin.methods)
        self.assertIn("align_sequences", self.plugin.methods)
        self.assertIn("otu_mapping", self.plugin.methods)
        self.assertIn("denovo_otu_mapping", self.plugin.methods)

        signature = self.plugin.methods["align_sequences"].signature
        self.assertIn("references", signature.inputs)
        self.assertNotIn("ref", signature.inputs)
        parameters = signature.parameters
        self.assertNotIn("best", parameters)
        for name in (
            "no_best",
            "score_split",
            "max_read_len",
        ):
            self.assertIn(name, parameters)

        for action_name in (
            "align_sequences",
            "otu_mapping",
            "denovo_otu_mapping",
        ):
            action_signature = self.plugin.methods[action_name].signature
            self.assertIn("references", action_signature.inputs)
            self.assertNotIn("ref", action_signature.inputs)
            action_parameters = action_signature.parameters
            for name in (
                "workdir",
                "kvdb",
                "idx_dir",
                "readb",
                "readfeed",
                "sq",
                "fastx",
                "sam",
                "blast",
                "aligned",
                "other",
                "zip_out",
                "out2",
                "sout",
                "paired_out",
                "pid",
                "index",
                "v",
                "otu_map",
                "de_novo_otu",
                "h",
                "version",
                "cmd",
                "task",
                "dbg_level",
                "flush_delay",
            ):
                self.assertNotIn(name, action_parameters)

        otu_parameters = self.plugin.methods["otu_mapping"].signature.parameters
        self.assertNotIn("no_best", otu_parameters)

    def test_method_defaults_match_sortmerna(self):
        expected_defaults = {
            "num_alignments": 1,
            "no_best": False,
            "min_lis": 2,
            "print_all_reads": False,
            "paired_in": False,
            "match": 2,
            "mismatch": -3,
            "gap_open": 5,
            "gap_ext": 2,
            "e": 1e-5,
            "f": False,
            "n": -3,
            "r": False,
            "score_split": False,
            "max_read_len": 30000,
            "id": 0.97,
            "coverage": 0.97,
            "passes": "18,9,3",
            "edges": 4,
            "num_seeds": 2,
            "full_search": False,
            "threads": 1,
            "l": 18,
            "m": 3072.0,
            "interval": 1,
            "max_pos": 10000,
        }
        methods = {
            "align_sequences": align_sequences,
            "otu_mapping": otu_mapping,
            "denovo_otu_mapping": denovo_otu_mapping,
        }

        for action_name, method in methods.items():
            expected = expected_defaults.copy()
            if action_name == "align_sequences":
                expected.pop("id")
                expected.pop("coverage")
            else:
                expected.pop("no_best")

            method_parameters = inspect.signature(method).parameters
            action_parameters = self.plugin.methods[action_name].signature.parameters
            for name, value in expected.items():
                self.assertEqual(method_parameters[name].default, value)
                self.assertEqual(action_parameters[name].default, value)

            self.assertEqual(action_parameters["threads"].default, 1)

    def test_registered_action_materializes_artifacts(self):
        action = self.plugin.methods["align_sequences"]
        with patch(
            "q2_sort_me_rna._rna_sorter.run_command",
            side_effect=self._mock_sortmerna,
        ) as mock_run:
            observed = action(
                references=self.reference_artifact,
                reads=self.single_reads_artifact,
            )

        commands = [call.args[0] for call in mock_run.call_args_list]
        command = next(command for command in commands if command[0] == "sortmerna")
        self.assertEqual(command[command.index("--threads") + 1], "1")
        for option in ("--num_alignments", "--min_lis"):
            self.assertNotIn(option, command)
        for option, value in (
            ("--match", "2"),
            ("--mismatch", "-3"),
            ("--passes", "18,9,3"),
            ("--max_pos", "10000"),
        ):
            self.assertEqual(command[command.index(option) + 1], value)

        self.assertEqual(str(observed.blast_aligned_seq.type), "FeatureData[BLAST6]")
        self.assertEqual(
            str(observed.fastx_aligned_seq.type), "SampleData[SequencesWithQuality]"
        )
        self.assertEqual(str(observed.alignment_map.type), "SampleData[AlignmentMap]")
