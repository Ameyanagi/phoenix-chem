"""Hazard evaluation, decomposition energetics, and classification."""

from phoenix.hazard.classification import HazardResult, evaluate_hazard
from phoenix.hazard.decomposition import (
    DecompositionComparison,
    DecompositionResult,
    calculate_max_decomposition,
)
from phoenix.hazard.functional_groups import detect_functional_groups
from phoenix.hazard.oxygen_balance import calculate_oxygen_balance

__all__ = [
    "HazardResult",
    "evaluate_hazard",
    "DecompositionComparison",
    "DecompositionResult",
    "calculate_max_decomposition",
    "detect_functional_groups",
    "calculate_oxygen_balance",
]
