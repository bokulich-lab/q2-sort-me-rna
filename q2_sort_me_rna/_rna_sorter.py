# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
from pathlib import Path
import shutil
import subprocess
import tempfile

import pandas as pd

from q2_types.feature_data import BLAST6Format
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    CasavaOneEightSingleLanePerSampleDirFmt,
)

from q2_sort_me_rna._utils import run_command


def sort_rna(
    references,
    reads,
    fastx: bool = None,
    sam: bool = None,
    sq: bool = None,
    blast: str = None,
    num_alignments: int = None,
    no_best: bool = None,
    min_lis: int = None,
    print_all_reads: bool = None,
    paired_in: bool = None,
    match: int = None,
    mismatch: int = None,
    gap_open: int = None,
    gap_ext: int = None,
    e: float = None,
    f: bool = None,
    n: int = None,
    r: bool = None,
    score_split: bool = None,
    max_read_len: int = None,
    id: float = None,
    coverage: float = None,
    de_novo_otu: bool = None,
    otu_map: bool = None,
    passes: str = None,
    edges: int = None,
    num_seeds: int = None,
    full_search: bool = None,
    threads: int = None,
    l: int = None,
    m: float = None,
    interval: int = None,
    max_pos: int = None,
):
    arguments = locals().copy()
    with tempfile.TemporaryDirectory(
        prefix="q2-sort-me-rna-"
    ) as temporary_workdir:
        arguments["workdir"] = temporary_workdir
        return _run_sort_rna(arguments)


def _run_sort_rna(arguments):
    arguments = arguments.copy()
    workdir = Path(arguments["workdir"])
    expected = ["blast", "fastx", "alignment_map"]
    if arguments.get("otu_map"):
        expected.append("otu_map")
    if arguments.get("de_novo_otu"):
        expected.append("denovo")

    sample_outputs = []
    for sample_number, (forward_read, reverse_read, sample_id) in enumerate(
        _get_read_sets(arguments["reads"]), start=1
    ):
        sample_workdir = workdir / f"sample-{sample_number:06d}"
        sample_workdir.mkdir()
        sample_arguments = arguments.copy()
        sample_arguments["workdir"] = str(sample_workdir)
        sample_arguments["ref"] = str(arguments["references"])
        del sample_arguments["references"]
        sample_arguments["reads"] = _un_gzip_file(
            forward_read, sample_workdir / "reads.fastq"
        )
        if reverse_read is not None:
            sample_arguments["reads_reverse"] = _un_gzip_file(
                reverse_read, sample_workdir / "reads-reverse.fastq"
            )

        command = ["sortmerna", *_parse_parameters(sample_arguments)]
        try:
            run_command(command, capture_output=True, text=True)
        except subprocess.CalledProcessError as error:
            detail = error.stderr or error.stdout or ""
            if isinstance(detail, bytes):
                detail = detail.decode(errors="replace")
            detail = detail.strip()
            message = (
                f"SortMeRNA failed with exit code {error.returncode} for "
                f"sample {sample_id!r}"
            )
            if detail:
                message = f"{message}: {detail}"
            raise RuntimeError(message) from error
        except OSError as error:
            raise RuntimeError(
                f"Unable to execute SortMeRNA for sample {sample_id!r}: "
                f"{error}"
            ) from error

        output_dir = sample_workdir / "out"
        if not output_dir.is_dir():
            raise RuntimeError(
                f"SortMeRNA did not create its expected output directory for "
                f"sample {sample_id!r}: {output_dir}"
            )

        outputs = _collect_outputs(output_dir, sample_id)
        missing = [name for name in expected if name not in outputs]
        if missing:
            raise RuntimeError(
                f"SortMeRNA completed for sample {sample_id!r} without "
                "producing expected output(s): " + ", ".join(missing)
            )
        sample_outputs.append(outputs)

    outputs = _merge_sample_outputs(sample_outputs, expected)
    return tuple(outputs[name] for name in expected)


def _merge_sample_outputs(sample_outputs, expected):
    outputs = {}
    for name in expected:
        values = [sample[name] for sample in sample_outputs]
        if name == "blast":
            outputs[name] = _merge_blast_formats(values)
        elif name in {"fastx", "denovo"}:
            outputs[name] = _merge_directory_formats(
                values, CasavaOneEightSingleLanePerSampleDirFmt
            )
        elif name == "alignment_map":
            outputs[name] = _merge_directory_formats(values, BAMDirFmt)
        elif name == "otu_map":
            mapping = pd.concat(values, axis="index").fillna(0).astype(int)
            if mapping.shape[1] == 0:
                raise ValueError(
                    "SortMeRNA produced no OTU assignments for any sample."
                )
            outputs[name] = mapping

    return outputs


def _merge_blast_formats(formats):
    merged = BLAST6Format()
    with Path(str(merged)).open("wb") as output_file:
        for blast_format in formats:
            with Path(str(blast_format)).open("rb") as input_file:
                shutil.copyfileobj(input_file, output_file)
    return merged


def _merge_directory_formats(formats, format_type):
    merged = format_type()
    destination_dir = Path(str(merged))
    for directory_format in formats:
        for input_path in Path(str(directory_format)).iterdir():
            destination = destination_dir / input_path.name
            if destination.exists():
                raise ValueError(
                    f"Cannot merge duplicate per-sample output: "
                    f"{input_path.name}"
                )
            shutil.copyfile(input_path, destination)
    return merged


def _parse_parameters(arguments):
    arguments = arguments.copy()
    _resolve_alignment_limit(arguments)
    uppercase_args = {"sq", "f", "n", "r", "l"}
    hyphenated_args = {"no_best"}
    duplicate_args = {"reads_reverse": "reads"}
    parameters = []

    for name, value in arguments.items():
        if value is None or value is False:
            continue
        name = duplicate_args.get(name, name)
        if name in hyphenated_args:
            name = name.replace("_", "-")
        if name in uppercase_args:
            name = name.upper()

        prefix = "-" if len(name) == 1 else "--"
        parameters.append(f"{prefix}{name}")
        if value is not True:
            parameters.append(str(value))

    return parameters


def _resolve_alignment_limit(arguments):
    num_alignments = arguments.get("num_alignments")
    min_lis = arguments.get("min_lis")
    if num_alignments is None or min_lis is None:
        return
    if num_alignments != 1 and min_lis != 2:
        raise ValueError(
            "SortMeRNA accepts either num_alignments or min_lis, not both."
        )
    if num_alignments == 1:
        del arguments["num_alignments"]
    if min_lis == 2:
        del arguments["min_lis"]


def _collect_outputs(output_dir, sample_id):
    outputs = {}
    for output_file in sorted(output_dir.iterdir()):
        name = output_file.name
        extension = _effective_extension(name)

        if extension == ".blast":
            outputs["blast"] = _construct_blast_fmt(output_file)
        elif _is_fastx(extension):
            key = "denovo" if "_denovo." in name else "fastx"
            outputs[key] = _construct_fastx_fmt(output_file, sample_id)
        elif extension == ".sam":
            outputs["alignment_map"] = _construct_bam_fmt(
                output_file, sample_id
            )
        elif name == "otu_map.txt":
            outputs["otu_map"] = _construct_otu_mapping(
                output_file, sample_id
            )

    return outputs


def _effective_extension(filename):
    path = Path(filename)
    if path.suffix == ".gz":
        return Path(path.stem).suffix
    return path.suffix


def _construct_blast_fmt(input_path):
    blast_fmt = BLAST6Format()
    shutil.copyfile(input_path, str(blast_fmt))
    return blast_fmt


def _construct_fastx_fmt(input_path, sample_id):
    input_path = Path(input_path)
    extension = _effective_extension(input_path.name)
    if not _is_fastq(extension):
        raise ValueError(
            f"Unsupported aligned sequence file type: {extension}"
        )

    if input_path.suffix != ".gz":
        compressed_path = input_path.with_suffix(f"{input_path.suffix}.gz")
        _gzip_file(input_path, compressed_path)
        input_path = compressed_path

    fastx_fmt = CasavaOneEightSingleLanePerSampleDirFmt()
    destination = (
        Path(str(fastx_fmt)) / f"{sample_id}_0_L001_R1_001.fastq.gz"
    )
    shutil.copyfile(input_path, destination)
    return fastx_fmt


def _construct_bam_fmt(input_path, sample_id):
    input_path = Path(input_path)
    if input_path.suffix == ".gz":
        uncompressed_path = input_path.with_suffix("")
        _un_gzip_file(input_path, uncompressed_path)
        input_path = uncompressed_path

    bam_fmt = BAMDirFmt()
    destination = Path(str(bam_fmt)) / f"{sample_id}.bam"
    command = [
        "samtools",
        "view",
        "-b",
        "-o",
        str(destination),
        str(input_path),
    ]
    try:
        run_command(command, capture_output=True, text=True)
    except subprocess.CalledProcessError as error:
        detail = error.stderr or error.stdout or ""
        if isinstance(detail, bytes):
            detail = detail.decode(errors="replace")
        detail = detail.strip()
        message = f"samtools failed with exit code {error.returncode}"
        if detail:
            message = f"{message}: {detail}"
        raise RuntimeError(message) from error
    except OSError as error:
        raise RuntimeError(f"Unable to execute samtools: {error}") from error

    return bam_fmt


def _construct_otu_mapping(input_path, sample_id):
    counts = {}
    with Path(input_path).open() as fh:
        for line_number, line in enumerate(fh, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if not fields or not fields[0]:
                raise ValueError(
                    f"Invalid OTU map record at line {line_number}."
                )
            counts[fields[0]] = len(fields) - 1

    mapping = pd.DataFrame([counts], index=[sample_id], dtype=int)
    mapping.index.name = "sample-id"
    return mapping


def _gzip_file(input_path, output_path):
    with Path(input_path).open("rb") as input_file:
        with gzip.open(output_path, "wb") as output_file:
            shutil.copyfileobj(input_file, output_file)
    return str(output_path)


def _un_gzip_file(input_path, output_path):
    with gzip.open(input_path, "rb") as input_file:
        with Path(output_path).open("wb") as output_file:
            shutil.copyfileobj(input_file, output_file)
    return str(output_path)


def _is_fastx(extension):
    return _is_fastq(extension) or _is_fasta(extension)


def _is_fastq(extension):
    return extension in {".fq", ".fastq"}


def _is_fasta(extension):
    return extension in {
        ".fasta",
        ".fas",
        ".fa",
        ".fna",
        ".ffn",
        ".faa",
        ".mpfa",
        ".frn",
    }


def _get_read_sets(reads_artifact):
    manifest = reads_artifact.manifest
    if manifest.empty:
        raise ValueError(
            "The reads artifact does not contain a gzipped FASTQ file."
        )

    read_sets = []
    layouts = set()
    for sample_id, row in manifest.sort_index().iterrows():
        if pd.isna(row["forward"]):
            raise ValueError(
                f"Sample {sample_id!r} does not contain a forward read."
            )

        forward_read = Path(row["forward"])
        if pd.isna(row["reverse"]):
            reverse_read = None
            layouts.add("single")
        else:
            reverse_read = Path(row["reverse"])
            layouts.add("paired")
        read_sets.append((forward_read, reverse_read, str(sample_id)))

    if len(layouts) > 1:
        raise ValueError(
            "A reads artifact cannot mix single-end and paired-end samples."
        )

    return read_sets
