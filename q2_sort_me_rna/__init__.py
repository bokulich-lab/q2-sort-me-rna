# flake8: noqa
# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

try:
    from ._version import __version__
except (ImportError, ModuleNotFoundError):
    try:
        from importlib.metadata import PackageNotFoundError, version

        __version__ = version("q2-sort-me-rna")
    except PackageNotFoundError:
        __version__ = "0.0.0+notfound"

from ._rna_sorter import sort_rna

__all__ = [
    "sort_rna",
]
