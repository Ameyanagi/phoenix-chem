"""
Reference management for thermodynamic data provenance.

This module provides structured literature citations for traceability.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from functools import lru_cache


@dataclass(frozen=True)
class Reference:
    """
    Structured literature citation.

    Attributes
    ----------
    author : str
        Author name(s)
    year : int
        Publication year
    title : str
        Publication title
    source_type : str
        Type: "book", "journal", "database", "software"
    doi : str | None
        Digital Object Identifier
    url : str | None
        URL for online resources
    """

    author: str
    year: int
    title: str
    source_type: str
    doi: str | None = None
    url: str | None = None

    def cite(self) -> str:
        """Return short citation string."""
        return f"{self.author} ({self.year})"

    def cite_full(self) -> str:
        """Return full citation string."""
        s = f"{self.author} ({self.year}). {self.title}"
        if self.doi:
            s += f" DOI: {self.doi}"
        elif self.url:
            s += f" {self.url}"
        return s


# Predefined standard references
BENSON_1976 = Reference(
    author="Benson, S.W.",
    year=1976,
    title="Thermochemical Kinetics, 2nd ed.",
    source_type="book",
    doi=None,
)

NIST_WEBBOOK = Reference(
    author="NIST",
    year=2024,
    title="NIST Chemistry WebBook, SRD 69",
    source_type="database",
    url="https://webbook.nist.gov/",
)

NIST_JANAF = Reference(
    author="Chase, M.W.",
    year=1998,
    title="NIST-JANAF Thermochemical Tables, 4th Ed.",
    source_type="book",
    doi="10.18434/T4D303",
)

PGRADD = Reference(
    author="Kang et al.",
    year=2021,
    title="pgradd: Python Group Additivity",
    source_type="journal",
    doi="10.1016/j.cpc.2021.108221",
)

CHEMICALS_LIB = Reference(
    author="Bell, C.",
    year=2024,
    title="chemicals Python library",
    source_type="software",
    url="https://github.com/CalebBell/chemicals",
)


class StandardReference(Enum):
    """Enum for quick access to standard references."""

    BENSON_1976 = BENSON_1976
    NIST_WEBBOOK = NIST_WEBBOOK
    NIST_JANAF = NIST_JANAF
    PGRADD = PGRADD
    CHEMICALS_LIB = CHEMICALS_LIB


@lru_cache(maxsize=100)
def get_reference(key: str) -> Reference:
    """
    Get reference by key name.

    Parameters
    ----------
    key : str
        Reference key (e.g., "BENSON_1976", "NIST_WEBBOOK")

    Returns
    -------
    Reference
        The reference object
    """
    try:
        return StandardReference[key.upper()].value
    except KeyError:
        raise ValueError(f"Unknown reference key: {key}")
