# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import hashlib
from io import BytesIO
from pathlib import Path
from unittest.mock import patch
from urllib.error import URLError

from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat

from q2_sort_me_rna._database import (
    _DATABASE_ASSETS,
    _DatabaseAsset,
    fetch_db,
)


class FetchDatabaseTests(TestPluginBase):
    package = "q2_sort_me_rna.tests"

    def setUp(self):
        super().setUp()
        self.reference = Path(
            self.get_data_path("seq/rna_references.fasta")
        ).read_bytes()
        self.archive = gzip.compress(self.reference)
        self.asset = _DatabaseAsset(
            filename="test-database.fasta.gz",
            sha256=hashlib.sha256(self.archive).hexdigest(),
        )

    def test_fetch_database(self):
        with patch.dict(_DATABASE_ASSETS, {"default": self.asset}):
            with patch(
                "q2_sort_me_rna._database.urlopen",
                return_value=BytesIO(self.archive),
            ) as mock_urlopen:
                observed = fetch_db()

        self.assertIsInstance(observed, DNAFASTAFormat)
        observed_lines = Path(str(observed)).read_bytes().splitlines()
        self.assertEqual(observed_lines[0], b">RNA_reference_U")
        self.assertEqual(observed_lines[1], b"ACGTACGTACGT")
        request = mock_urlopen.call_args.args[0]
        self.assertEqual(
            request.full_url,
            "https://github.com/sortmerna/sortmerna/releases/download/"
            "v7.0.0/test-database.fasta.gz",
        )

    def test_fetch_database_rejects_checksum_mismatch(self):
        asset = _DatabaseAsset(
            filename="test-database.fasta.gz",
            sha256="0" * 64,
        )
        with patch.dict(_DATABASE_ASSETS, {"default": asset}):
            with patch(
                "q2_sort_me_rna._database.urlopen",
                return_value=BytesIO(self.archive),
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "failed SHA-256 verification"
                ):
                    fetch_db()

    def test_fetch_database_reports_network_failure(self):
        with patch(
            "q2_sort_me_rna._database.urlopen",
            side_effect=URLError("network unavailable"),
        ):
            with self.assertRaisesRegex(
                RuntimeError, "network unavailable"
            ):
                fetch_db()

    def test_fetch_database_rejects_invalid_gzip(self):
        archive = b"not a gzip archive"
        asset = _DatabaseAsset(
            filename="test-database.fasta.gz",
            sha256=hashlib.sha256(archive).hexdigest(),
        )
        with patch.dict(_DATABASE_ASSETS, {"default": asset}):
            with patch(
                "q2_sort_me_rna._database.urlopen",
                return_value=BytesIO(archive),
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "could not be decompressed"
                ):
                    fetch_db()

    def test_fetch_database_is_registered(self):
        action = self.plugin.methods["fetch_db"]
        self.assertEqual(action.signature.parameters["database"].default,
                         "default")

        with patch.dict(_DATABASE_ASSETS, {"default": self.asset}):
            with patch(
                "q2_sort_me_rna._database.urlopen",
                return_value=BytesIO(self.archive),
            ):
                observed = action(database="default")

        self.assertEqual(
            str(observed.reference_db.type), "FeatureData[Sequence]"
        )
