# q2-sort-me-rna

`q2-sort-me-rna` is a [QIIME 2](https://qiime2.org) community plugin that
wraps [SortMeRNA](https://github.com/sortmerna/sortmerna) for filtering and
aligning ribosomal RNA reads.

The plugin provides four methods:

- `fetch-db` downloads a prepared upstream rRNA reference database and
  returns it as `FeatureData[Sequence]`.
- `align-sequences` returns BLAST, FASTQ, and BAM alignment-map
  representations of aligned reads.
- `otu-mapping` additionally returns a feature table of reference assignments.
- `denovo-otu-mapping` additionally returns reads assigned during de novo OTU
  mapping.

## Installation

Create a development environment from the concise environment specification:

```shell
conda env create \
  -n q2-sort-me-rna-qiime2-tiny-dev \
  --file environment-files/q2-sort-me-rna-qiime2-tiny-dev.yml
conda activate q2-sort-me-rna-qiime2-tiny-dev
make install
qiime dev refresh-cache
```

Confirm that QIIME 2 discovered the plugin:

```shell
qiime sort-me-rna --help
```

## Prepared reference databases

Fetch the upstream recommended database as a provenance-tracked QIIME 2
artifact:

```shell
qiime sort-me-rna fetch-db \
  --p-database default \
  --o-reference-db smr-default-reference.qza
```

The available variants are `fast`, `default`, `sensitive`, and
`sensitive-rfam-seeds`. The action downloads the corresponding asset from the
SortMeRNA 7.0.0 GitHub release and verifies its upstream SHA-256 digest before
decompression. Although these assets are distributed with SortMeRNA 7.0.0,
their filenames retain the upstream `smr_v4.3` database label. SortMeRNA
recommends `default`; the two sensitive variants offer a small accuracy gain
but run at least twice as slowly. The upstream references use RNA bases and
mixed case; the action uppercases sequence records and normalizes `U` to `T`
so the result is a standard QIIME 2 `FeatureData[Sequence]` artifact accepted
directly by this plugin's alignment actions. FASTA identifiers and descriptions
are left unchanged.

The environment requires SortMeRNA 7.x. The plugin exposes alignment and
scoring controls that change data represented by the declared QIIME 2 outputs.
Their defaults match SortMeRNA 7.0.0, except that QIIME actions intentionally
default to one processing thread.
Filesystem paths, caches, compression, execution phases, and output switches
are managed by the wrapper so an action cannot omit or redirect a required
result. SortMeRNA's `no-best` mode is excluded from the two OTU actions because
v7 declares it incompatible with OTU-map generation.

## Example

```shell
qiime sort-me-rna align-sequences \
  --i-reads q2_sort_me_rna/tests/data/raw_sequence.qza \
  --i-references q2_sort_me_rna/tests/data/rrna_references.qza \
  --output-dir output \
  --verbose
```

Each action uses an isolated temporary SortMeRNA workspace and removes it after
collecting the declared QIIME 2 results.

## Input processing

SortMeRNA accepts one single-end file or one paired-end file pair per run. For
multi-sample `SampleData` artifacts, the plugin runs SortMeRNA independently for
each sample and combines the results into the declared QIIME 2 outputs.

Reads must be FASTQ; reference sequences must be FASTA. An empty OTU map cannot
be represented as a QIIME 2 `FeatureTable[Frequency]`, so OTU-mapping actions
raise an error when SortMeRNA reports no assignments.

## Development

Run the test and lint suites with:

```shell
make test
make lint
```

Unit tests mock the external SortMeRNA boundary. End-to-end tests run only when
the `sortmerna` executable is available.

## Citation

If you use this plugin, cite:

Kopylova E, Noé L, Touzet H. SortMeRNA: fast and accurate filtering of
ribosomal RNAs in metatranscriptomic data. *Bioinformatics*. 2012;28(24):
3211–3217. <https://doi.org/10.1093/bioinformatics/bts611>.
