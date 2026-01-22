"""
Functional group detection for hazard alerts.

Uses SMARTS patterns to identify reactive functional groups that
are associated with chemical instability or explosive potential.

Reference: CHETAH documentation, NFPA 704 hazard classifications
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from rdkit import Chem

if TYPE_CHECKING:
    from phoenix.core.compound import Compound


@dataclass(frozen=True)
class FunctionalGroupAlert:
    """A detected reactive functional group."""

    name: str
    smarts: str
    count: int
    description: str


# SMARTS patterns for reactive functional groups
# Format: (name, SMARTS, description)
REACTIVE_GROUPS: list[tuple[str, str, str]] = [
    # Nitro compounds
    (
        "Nitro group",
        "[N+](=O)[O-]",
        "Nitro compounds can be shock-sensitive and explosive",
    ),
    (
        "Aromatic nitro",
        "c[N+](=O)[O-]",
        "Aromatic nitro compounds are potential explosives",
    ),
    # Peroxides
    (
        "Peroxide",
        "[OX2][OX2]",
        "Peroxides are sensitive to heat, shock, and friction",
    ),
    (
        "Organic peroxide",
        "[CX4][OX2][OX2][CX4]",
        "Organic peroxides can decompose explosively",
    ),
    (
        "Peroxy acid",
        "[CX3](=O)[OX2][OX2]",
        "Peroxy acids are oxidizers and can be explosive",
    ),
    (
        "Hydroperoxide",
        "[OX2][OX2H]",
        "Hydroperoxides are unstable and sensitive",
    ),
    # Azo and diazo compounds
    (
        "Azo group",
        "[NX2]=[NX2]",
        "Azo compounds can decompose with gas release",
    ),
    (
        "Diazo group",
        "[CX3]=[N+]=[N-]",
        "Diazo compounds are highly reactive",
    ),
    (
        "Diazonium",
        "[N+]#N",
        "Diazonium salts are shock-sensitive explosives",
    ),
    # Azides
    (
        "Azide group",
        "[N-]=[N+]=[N-]",
        "Azides are primary explosives",
    ),
    (
        "Organic azide",
        "[CX4][N-]=[N+]=[N-]",
        "Organic azides can detonate",
    ),
    # N-oxides and nitroso
    (
        "N-oxide",
        "[NX3+]([O-])",
        "N-oxides can be oxidizers",
    ),
    (
        "Nitroso group",
        "[NX2]=O",
        "Nitroso compounds can be unstable",
    ),
    # Nitrates and nitrites
    (
        "Nitrate ester",
        "[OX2][N+](=O)[O-]",
        "Nitrate esters are explosives (e.g., nitroglycerin)",
    ),
    (
        "Nitrite ester",
        "[OX2][NX2]=O",
        "Nitrite esters can be unstable",
    ),
    # Nitrogen-nitrogen bonds
    (
        "Hydrazine",
        "[NX3H2][NX3H2]",
        "Hydrazines are toxic and can be explosive",
    ),
    (
        "Hydrazide",
        "[CX3](=O)[NX3][NX3]",
        "Hydrazides can decompose violently",
    ),
    # Strained rings
    (
        "Epoxide",
        "C1OC1",
        "Epoxides are reactive and can polymerize violently",
    ),
    (
        "Aziridine",
        "C1NC1",
        "Aziridines are strained and reactive",
    ),
    # Acyl halides and anhydrides
    (
        "Acyl halide",
        "[CX3](=O)[F,Cl,Br,I]",
        "Acyl halides are highly reactive",
    ),
    (
        "Acid anhydride",
        "[CX3](=O)[OX2][CX3](=O)",
        "Acid anhydrides are reactive with water",
    ),
    # Metal-containing (catch for future)
    (
        "Acetylide",
        "[CX2]#[CX2-]",
        "Acetylides can be shock-sensitive",
    ),
]


def detect_functional_groups(compound: Compound) -> list[FunctionalGroupAlert]:
    """
    Detect reactive functional groups in a compound.

    Parameters
    ----------
    compound : Compound
        Compound to analyze

    Returns
    -------
    list[FunctionalGroupAlert]
        List of detected reactive functional groups
    """
    alerts: list[FunctionalGroupAlert] = []
    mol = compound.rdmol

    for name, smarts, description in REACTIVE_GROUPS:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue

        matches = mol.GetSubstructMatches(pattern)
        if matches:
            alerts.append(
                FunctionalGroupAlert(
                    name=name,
                    smarts=smarts,
                    count=len(matches),
                    description=description,
                )
            )

    return alerts


def get_alert_names(compound: Compound) -> list[str]:
    """
    Get names of all detected reactive functional groups.

    Parameters
    ----------
    compound : Compound
        Compound to analyze

    Returns
    -------
    list[str]
        Names of detected groups
    """
    alerts = detect_functional_groups(compound)
    return [alert.name for alert in alerts]
