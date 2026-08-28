# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from dataclasses import dataclass
import gzip
import hashlib
from pathlib import Path
import tempfile
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from q2_types.feature_data import DNAFASTAFormat

from q2_sort_me_rna import __version__


SORTMERNA_DATABASE_RELEASE = "7.0.0"
_RELEASE_BASE_URL = (
    "https://github.com/sortmerna/sortmerna/releases/download/"
    f"v{SORTMERNA_DATABASE_RELEASE}"
)


@dataclass(frozen=True)
class _DatabaseAsset:
    filename: str
    sha256: str


_DATABASE_ASSETS = {
    "fast": _DatabaseAsset(
        filename="smr_v4.3_fast_db.fasta.gz",
        sha256=(
            "3ad47d6a9296e891f6165ec8152eb8e3"
            "8c7ed9f37a2b94e06f71b67ec694487b"
        ),
    ),
    "default": _DatabaseAsset(
        filename="smr_v4.3_default_db.fasta.gz",
        sha256=(
            "abd7e0ffdf4710800a954af595811049"
            "056705bac500999ca5d8f53b1c256752"
        ),
    ),
    "sensitive": _DatabaseAsset(
        filename="smr_v4.3_sensitive_db.fasta.gz",
        sha256=(
            "b9efe5fd8cf2b631c19f3b6ddf7b0e5"
            "c68351ed3409b9ce6f3c156052a1a830f"
        ),
    ),
    "sensitive-rfam-seeds": _DatabaseAsset(
        filename="smr_v4.3_sensitive_db_rfam_seeds.fasta.gz",
        sha256=(
            "29de2fd40280d6fb3832dcd69d61ab88"
            "aa8ea90065f0bf2a30eac3de23c80145"
        ),
    ),
}
DATABASE_CHOICES = tuple(_DATABASE_ASSETS)


def fetch_db(database: str = "default") -> DNAFASTAFormat:
    """Fetch and QIIME-normalize a SortMeRNA v7.0.0 rRNA database."""
    asset = _DATABASE_ASSETS[database]
    url = f"{_RELEASE_BASE_URL}/{asset.filename}"

    with tempfile.TemporaryDirectory(
        prefix="q2-sort-me-rna-database-"
    ) as temporary_directory:
        archive = Path(temporary_directory) / asset.filename
        observed_sha256 = _download(url, archive)
        if observed_sha256 != asset.sha256:
            raise RuntimeError(
                f"The downloaded SortMeRNA '{database}' database failed "
                "SHA-256 verification. "
                f"Expected {asset.sha256}, observed {observed_sha256}."
            )

        reference_sequences = DNAFASTAFormat()
        try:
            with gzip.open(archive, "rb") as input_file:
                with Path(str(reference_sequences)).open("wb") as output_file:
                    _normalize_rna_fasta(input_file, output_file)
        except (EOFError, gzip.BadGzipFile, OSError) as error:
            raise RuntimeError(
                f"The downloaded SortMeRNA '{database}' database could not "
                f"be decompressed: {error}"
            ) from error

    return reference_sequences


def _normalize_rna_fasta(input_file, output_file):
    for line in input_file:
        if not line.startswith(b">"):
            line = line.upper().replace(b"U", b"T")
        output_file.write(line)


def _download(url, destination):
    request = Request(
        url,
        headers={"User-Agent": f"q2-sort-me-rna/{__version__}"},
    )
    digest = hashlib.sha256()

    try:
        with urlopen(request, timeout=120) as response:
            with Path(destination).open("wb") as output_file:
                while chunk := response.read(1024 * 1024):
                    output_file.write(chunk)
                    digest.update(chunk)
    except HTTPError as error:
        raise RuntimeError(
            f"Unable to download the SortMeRNA database from {url}: "
            f"HTTP {error.code}."
        ) from error
    except (URLError, OSError) as error:
        reason = getattr(error, "reason", error)
        raise RuntimeError(
            f"Unable to download the SortMeRNA database from {url}: "
            f"{reason}."
        ) from error

    return digest.hexdigest()
