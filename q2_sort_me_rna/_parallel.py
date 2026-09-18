# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData


def _run_parallel(ctx, action_name, references, reads, num_partitions, arguments):
    if reads.type <= SampleData[SequencesWithQuality]:
        partition = ctx.get_action("types", "partition_samples_single")
    elif reads.type <= SampleData[PairedEndSequencesWithQuality]:
        partition = ctx.get_action("types", "partition_samples_paired")
    else:
        raise ValueError(f"Unsupported read type: {reads.type}")

    worker = ctx.get_action("sort_me_rna", f"_{action_name}")
    collate_blasts = ctx.get_action("sort_me_rna", "_collate_blasts")
    collate_alignment_maps = ctx.get_action("sort_me_rna", "_collate_alignment_maps")
    combine_seqs = ctx.get_action("fondue", "combine_seqs")
    (partitions,) = partition(reads, num_partitions)

    output_count = {
        "align_sequences": 3,
        "otu_mapping": 4,
        "denovo_otu_mapping": 5,
    }[action_name]
    outputs = [[] for _ in range(output_count)]
    for sample_reads in partitions.values():
        for collection, result in zip(
            outputs, worker(references=references, reads=sample_reads, **arguments)
        ):
            collection.append(result)

    (blasts,) = collate_blasts(outputs[0])
    (fastxs,) = combine_seqs(outputs[1])
    (bams,) = collate_alignment_maps(outputs[2])
    merged = [blasts, fastxs, bams]
    if output_count >= 4:
        merge_tables = ctx.get_action("feature_table", "merge")
        (otu_map,) = merge_tables(outputs[3])
        merged.append(otu_map)
    if output_count == 5:
        (denovo_reads,) = combine_seqs(outputs[4])
        merged.append(denovo_reads)
    return tuple(merged)


def align_sequences(
    ctx,
    references,
    reads,
    num_alignments=1,
    no_best=False,
    min_lis=2,
    print_all_reads=False,
    paired_in=False,
    match=2,
    mismatch=-3,
    gap_open=5,
    gap_ext=2,
    e=1e-5,
    f=False,
    n=-3,
    r=False,
    score_split=False,
    max_read_len=30000,
    passes="18,9,3",
    edges=4,
    num_seeds=2,
    full_search=False,
    threads=1,
    l=18,
    m=3072.0,
    interval=1,
    max_pos=10000,
    num_partitions=None,
):
    arguments = locals().copy()
    for name in ("ctx", "references", "reads", "num_partitions"):
        del arguments[name]
    return _run_parallel(
        ctx, "align_sequences", references, reads, num_partitions, arguments
    )


def otu_mapping(
    ctx,
    references,
    reads,
    num_alignments=1,
    min_lis=2,
    print_all_reads=False,
    paired_in=False,
    match=2,
    mismatch=-3,
    gap_open=5,
    gap_ext=2,
    e=1e-5,
    f=False,
    n=-3,
    r=False,
    score_split=False,
    max_read_len=30000,
    id=0.97,
    coverage=0.97,
    passes="18,9,3",
    edges=4,
    num_seeds=2,
    full_search=False,
    threads=1,
    l=18,
    m=3072.0,
    interval=1,
    max_pos=10000,
    num_partitions=None,
):
    arguments = locals().copy()
    for name in ("ctx", "references", "reads", "num_partitions"):
        del arguments[name]
    return _run_parallel(
        ctx, "otu_mapping", references, reads, num_partitions, arguments
    )


def denovo_otu_mapping(
    ctx,
    references,
    reads,
    num_alignments=1,
    min_lis=2,
    print_all_reads=False,
    paired_in=False,
    match=2,
    mismatch=-3,
    gap_open=5,
    gap_ext=2,
    e=1e-5,
    f=False,
    n=-3,
    r=False,
    score_split=False,
    max_read_len=30000,
    id=0.97,
    coverage=0.97,
    passes="18,9,3",
    edges=4,
    num_seeds=2,
    full_search=False,
    threads=1,
    l=18,
    m=3072.0,
    interval=1,
    max_pos=10000,
    num_partitions=None,
):
    arguments = locals().copy()
    for name in ("ctx", "references", "reads", "num_partitions"):
        del arguments[name]
    return _run_parallel(
        ctx, "denovo_otu_mapping", references, reads, num_partitions, arguments
    )
