# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (
    Bool,
    Choices,
    Citations,
    Float,
    Int,
    Plugin,
    Range,
    Str,
    TypeMatch,
)
from q2_types.feature_data import BLAST6, FeatureData, Sequence
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.per_sample_sequences import (
    AlignmentMap,
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData

from q2_sort_me_rna import __version__
from q2_sort_me_rna._methods import (
    align_sequences,
    denovo_otu_mapping,
    otu_mapping,
)
from q2_sort_me_rna._database import DATABASE_CHOICES, fetch_db


citations = Citations.load("citations.bib", package="q2_sort_me_rna")
sortmerna_citation = citations["sortmerna-2012"]

plugin = Plugin(
    name="sort-me-rna",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-sort-me-rna",
    package="q2_sort_me_rna",
    description=(
        "A QIIME 2 wrapper for SortMeRNA, a sequence-alignment tool for "
        "filtering ribosomal RNA reads."
    ),
    short_description="Filter and align ribosomal RNA reads with SortMeRNA.",
    citations=[sortmerna_citation],
)

nonnegative_int = Int % Range(0, None)
positive_int = Int % Range(1, None)
nonnegative_float = Float % Range(0, None)
proportion = Float % Range(0, 1, inclusive_end=True)
edge_count = Int % Range(1, 11)
seed_length = Int % Range(8, 27)

# These v7 options change data represented by the actions' declared outputs.
# Filesystem, serialization, indexing-task, and developer controls are owned
# by the wrapper and are intentionally absent from the QIIME 2 interface.
sort_me_rna_parameters = {
    "num_alignments": nonnegative_int,
    "no_best": Bool,
    "min_lis": nonnegative_int,
    "print_all_reads": Bool,
    "paired_in": Bool,
    "match": Int,
    "mismatch": Int,
    "gap_open": Int,
    "gap_ext": Int,
    "e": nonnegative_float,
    "f": Bool,
    "n": Int,
    "r": Bool,
    "score_split": Bool,
    "max_read_len": positive_int,
    "id": proportion,
    "coverage": proportion,
    "passes": Str,
    "edges": edge_count,
    "num_seeds": positive_int,
    "full_search": Bool,
    "threads": positive_int,
    "l": seed_length,
    "m": Float,
    "interval": nonnegative_int,
    "max_pos": nonnegative_int,
}

sort_me_rna_parameter_descriptions = {
    "num_alignments": (
        "Number of alignments to report per read; zero reports all "
        "alignments that meet the E-value threshold."
    ),
    "no_best": (
        "Disable the exhaustive best-alignment search and report the first "
        "num-alignments hits that meet the E-value threshold."
    ),
    "min_lis": "Minimum longest-increasing-subsequence score.",
    "print_all_reads": "Include every read in alignment reports.",
    "paired_in": "Place both mates in aligned output when either aligns.",
    "match": "Smith-Waterman match score.",
    "mismatch": "Smith-Waterman mismatch penalty.",
    "gap_open": "Smith-Waterman gap-opening penalty.",
    "gap_ext": "Smith-Waterman gap-extension penalty.",
    "e": "Maximum E-value for reported alignments.",
    "f": "Search the forward strand only.",
    "n": (
        "Smith-Waterman score for ambiguous N bases. By default this follows "
        "the mismatch penalty."
    ),
    "r": "Search the reverse-complement strand only.",
    "score_split": (
        "Calculate the minimum Smith-Waterman score per read split instead "
        "of across all reads."
    ),
    "max_read_len": "Maximum accepted read length in nucleotides.",
    "id": "Minimum alignment identity as a fraction from zero to one.",
    "coverage": "Minimum query coverage as a fraction from zero to one.",
    "passes": (
        "Comma-delimited seed-search intervals for three passes. The default "
        "is derived from the seed length as L,L/2,3."
    ),
    "edges": "Bases to extend on each edge of an alignment region.",
    "num_seeds": "Minimum number of seeds required by the seed-search filter.",
    "full_search": "Search one-error seed matches after exact matches.",
    "threads": "Number of worker threads; the QIIME default is one.",
    "l": "Reference-index seed length; must be an even integer.",
    "m": "Reference-index memory limit.",
    "interval": "Reference-index seed interval.",
    "max_pos": (
        "Maximum reference positions stored per seed. SortMeRNA v7.0.0 uses "
        "10000 at runtime despite displaying 1000 in its help text."
    ),
}

alignment_parameters = {
    name: parameter
    for name, parameter in sort_me_rna_parameters.items()
    if name not in {"id", "coverage"}
}
otu_parameters = {
    name: parameter
    for name, parameter in sort_me_rna_parameters.items()
    if name != "no_best"
}
alignment_parameter_descriptions = {
    name: sort_me_rna_parameter_descriptions[name]
    for name in alignment_parameters
}
otu_parameter_descriptions = {
    name: sort_me_rna_parameter_descriptions[name]
    for name in otu_parameters
}

read_type = TypeMatch(
    [SequencesWithQuality, PairedEndSequencesWithQuality]
)
common_inputs = {
    "references": FeatureData[Sequence],
    "reads": SampleData[read_type],
}
common_input_descriptions = {
    "references": "Reference ribosomal RNA sequences in FASTA format.",
    "reads": (
        "Single-end or paired-end demultiplexed reads. Each sample is "
        "processed independently and represented in the combined outputs."
    ),
}
common_output_descriptions = {
    "blast_aligned_seq": "Aligned reads in BLAST tabular format.",
    "fastx_aligned_seq": "Aligned reads in FASTQ format.",
    "alignment_map": "Per-sample aligned reads in BAM format.",
}

plugin.methods.register_function(
    function=fetch_db,
    inputs={},
    parameters={"database": Str % Choices(DATABASE_CHOICES)},
    outputs=[("reference_db", FeatureData[Sequence])],
    input_descriptions={},
    parameter_descriptions={
        "database": (
            "Prepared upstream database variant. 'default' is recommended; "
            "'fast' is slightly less sensitive, while 'sensitive' and "
            "'sensitive-rfam-seeds' are more comprehensive but slower."
        )
    },
    output_descriptions={
        "reference_db": (
            "Prepared rRNA reference sequences from the SortMeRNA v7.0.0 "
            "release, represented with thymine for compatibility with the "
            "QIIME 2 DNA sequence type."
        )
    },
    name="Fetch a prepared SortMeRNA database",
    description=(
        "Download, verify, and decompress one of the prepared rRNA reference "
        "databases distributed with SortMeRNA v7.0.0. Upstream retains the "
        "v4.3 label in these database filenames. Sequences are normalized "
        "to uppercase DNA (uracil becomes thymine) so the database can be "
        "used directly by the alignment actions as FeatureData[Sequence]."
    ),
    citations=[sortmerna_citation],
)

plugin.methods.register_function(
    function=align_sequences,
    inputs=common_inputs,
    parameters=alignment_parameters,
    outputs=[
        ("blast_aligned_seq", FeatureData[BLAST6]),
        ("fastx_aligned_seq", SampleData[SequencesWithQuality]),
        ("alignment_map", SampleData[AlignmentMap]),
    ],
    input_descriptions=common_input_descriptions,
    parameter_descriptions=alignment_parameter_descriptions,
    output_descriptions=common_output_descriptions,
    name="Align sequences",
    description=(
        "Align single-end or paired-end reads against an rRNA reference "
        "database and return BLAST, FASTQ, and BAM alignment representations."
    ),
    citations=[sortmerna_citation],
)

plugin.methods.register_function(
    function=otu_mapping,
    inputs=common_inputs,
    parameters=otu_parameters,
    outputs=[
        ("blast_aligned_seq", FeatureData[BLAST6]),
        ("fastx_aligned_seq", SampleData[SequencesWithQuality]),
        ("alignment_map", SampleData[AlignmentMap]),
        ("otu_mapping", FeatureTable[Frequency]),
    ],
    input_descriptions=common_input_descriptions,
    parameter_descriptions=otu_parameter_descriptions,
    output_descriptions={
        **common_output_descriptions,
        "otu_mapping": "Counts of reads assigned to each reference sequence.",
    },
    name="Align sequences and create an OTU map",
    description=(
        "Align reads against an rRNA reference database and also return a "
        "feature table summarizing reference assignments."
    ),
    citations=[sortmerna_citation],
)

plugin.methods.register_function(
    function=denovo_otu_mapping,
    inputs=common_inputs,
    parameters=otu_parameters,
    outputs=[
        ("blast_aligned_seq", FeatureData[BLAST6]),
        ("fastx_aligned_seq", SampleData[SequencesWithQuality]),
        ("alignment_map", SampleData[AlignmentMap]),
        ("otu_mapping", FeatureTable[Frequency]),
        ("denovo_aligned_seq", SampleData[SequencesWithQuality]),
    ],
    input_descriptions=common_input_descriptions,
    parameter_descriptions=otu_parameter_descriptions,
    output_descriptions={
        **common_output_descriptions,
        "otu_mapping": "Counts of reads assigned to each reference sequence.",
        "denovo_aligned_seq": "Reads assigned during de novo OTU mapping.",
    },
    name="Align sequences with de novo OTU mapping",
    description=(
        "Align reads against an rRNA reference database, build an OTU map, "
        "and return reads assigned during de novo OTU mapping."
    ),
    citations=[sortmerna_citation],
)
