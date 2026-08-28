# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from q2_types.feature_data import BLAST6Format, DNAFASTAFormat
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    CasavaOneEightSingleLanePerSampleDirFmt,
)

from q2_sort_me_rna._rna_sorter import sort_rna


def align_sequences(
    references: DNAFASTAFormat,
    reads: CasavaOneEightSingleLanePerSampleDirFmt,
    num_alignments: int = 1,
    no_best: bool = False,
    min_lis: int = 2,
    print_all_reads: bool = False,
    paired_in: bool = False,
    match: int = 2,
    mismatch: int = -3,
    gap_open: int = 5,
    gap_ext: int = 2,
    e: float = 1e-5,
    f: bool = False,
    n: int = -3,
    r: bool = False,
    score_split: bool = False,
    max_read_len: int = 30000,
    passes: str = "18,9,3",
    edges: int = 4,
    num_seeds: int = 2,
    full_search: bool = False,
    threads: int = 1,
    l: int = 18,
    m: float = 3072.0,
    interval: int = 1,
    max_pos: int = 10000,
) -> (
    BLAST6Format,
    CasavaOneEightSingleLanePerSampleDirFmt,
    BAMDirFmt,
):
    arguments = locals()
    return sort_rna(blast="1", fastx=True, sam=True, sq=True, **arguments)


def otu_mapping(
    references: DNAFASTAFormat,
    reads: CasavaOneEightSingleLanePerSampleDirFmt,
    num_alignments: int = 1,
    min_lis: int = 2,
    print_all_reads: bool = False,
    paired_in: bool = False,
    match: int = 2,
    mismatch: int = -3,
    gap_open: int = 5,
    gap_ext: int = 2,
    e: float = 1e-5,
    f: bool = False,
    n: int = -3,
    r: bool = False,
    score_split: bool = False,
    max_read_len: int = 30000,
    id: float = 0.97,
    coverage: float = 0.97,
    passes: str = "18,9,3",
    edges: int = 4,
    num_seeds: int = 2,
    full_search: bool = False,
    threads: int = 1,
    l: int = 18,
    m: float = 3072.0,
    interval: int = 1,
    max_pos: int = 10000,
) -> (
    BLAST6Format,
    CasavaOneEightSingleLanePerSampleDirFmt,
    BAMDirFmt,
    pd.DataFrame,
):
    arguments = locals()
    return sort_rna(blast="1", fastx=True, sam=True, sq=True, otu_map=True, **arguments)


def denovo_otu_mapping(
    references: DNAFASTAFormat,
    reads: CasavaOneEightSingleLanePerSampleDirFmt,
    num_alignments: int = 1,
    min_lis: int = 2,
    print_all_reads: bool = False,
    paired_in: bool = False,
    match: int = 2,
    mismatch: int = -3,
    gap_open: int = 5,
    gap_ext: int = 2,
    e: float = 1e-5,
    f: bool = False,
    n: int = -3,
    r: bool = False,
    score_split: bool = False,
    max_read_len: int = 30000,
    id: float = 0.97,
    coverage: float = 0.97,
    passes: str = "18,9,3",
    edges: int = 4,
    num_seeds: int = 2,
    full_search: bool = False,
    threads: int = 1,
    l: int = 18,
    m: float = 3072.0,
    interval: int = 1,
    max_pos: int = 10000,
) -> (
    BLAST6Format,
    CasavaOneEightSingleLanePerSampleDirFmt,
    BAMDirFmt,
    pd.DataFrame,
    CasavaOneEightSingleLanePerSampleDirFmt,
):
    arguments = locals()
    return sort_rna(
        blast="1",
        fastx=True,
        sam=True,
        sq=True,
        otu_map=True,
        de_novo_otu=True,
        **arguments,
    )
