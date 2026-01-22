"""
PHOENIX - Physicochemical Hazard Observation & Energetics Numerical Indexing eXpert

A Python library for reactive chemical hazard evaluation, successor to ASTM CHETAH.

Warning
-------
This library is for screening purposes only. Results must not be used as the
sole basis for safety decisions. Experimental validation is required before
handling energetic materials.

Basic Usage
-----------
>>> from phoenix import Compound
>>> compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")  # nitrobenzene
>>> print(f"Formula: {compound.formula}")
Formula: C6H5NO2
>>> print(f"OB%: {compound.oxygen_balance_percent:.1f}%")
OB%: -163.0%
>>> result = compound.evaluate_hazard()
>>> print(f"Hazard Class: {result.hazard_class}")
Hazard Class: MEDIUM

Batch Processing
----------------
>>> from phoenix import screen
>>> results = screen(["CCO", "c1ccccc1[N+](=O)[O-]"])
>>> print(results.dataframe[["formula", "hazard_class"]])

References
----------
- ASTM E659: Standard Test Method for CHETAH
- Benson, S.W. "Thermochemical Kinetics" (Wiley, 1976)
- Meyer et al., "Explosives" 6th Ed. (Wiley-VCH, 2007)
"""

from phoenix.batch import BatchResult, screen
from phoenix.core import Auto, Compound, Reaction, ReactionSpecies
from phoenix.thermo.models import ThermoProperty
from phoenix.exceptions import (
    BalanceError,
    DecompositionError,
    InvalidSmilesError,
    MissingGroupError,
    OverconstrainedError,
    PhoenixError,
    UnderconstrainedError,
    UnsupportedElementError,
    UnsupportedStructureError,
)
from phoenix.hazard import HazardResult

__version__ = "0.1.0"

__all__ = [
    # Core
    "Compound",
    "Reaction",
    "ReactionSpecies",
    "Auto",
    # Thermo
    "ThermoProperty",
    # Batch processing
    "screen",
    "BatchResult",
    # Results
    "HazardResult",
    # Exceptions
    "PhoenixError",
    "InvalidSmilesError",
    "UnsupportedElementError",
    "UnsupportedStructureError",
    "MissingGroupError",
    "DecompositionError",
    "BalanceError",
    "OverconstrainedError",
    "UnderconstrainedError",
]
