"""
Hazard classification based on CHETAH criteria.

Evaluates compounds against ASTM E659 hazard criteria:
- Criterion 1: ΔHd < -300 cal/g (High instability)
- Criterion 2: ΔHd < -100 cal/g (Medium instability)
- Criterion 3: -200% < OB% < +100% (Oxidizer/fuel balance concern)
- Criterion 4: Functional group alerts (Known reactive moieties)

Reference: ASTM E659 - Standard Test Method for CHETAH
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

from phoenix.hazard.decomposition import DecompositionResult, calculate_max_decomposition
from phoenix.hazard.functional_groups import get_alert_names
from phoenix.hazard.oxygen_balance import calculate_oxygen_balance

if TYPE_CHECKING:
    from phoenix.core.compound import Compound

# CHETAH criteria thresholds
CRITERION_1_THRESHOLD = -300.0  # cal/g (High hazard)
CRITERION_2_THRESHOLD = -100.0  # cal/g (Medium hazard)
CRITERION_3_OB_LOW = -200.0  # % (Lower bound for OB concern)
CRITERION_3_OB_HIGH = 100.0  # % (Upper bound for OB concern)

HazardClass = Literal["HIGH", "MEDIUM", "LOW"]


@dataclass(frozen=True)
class HazardResult:
    """
    Complete hazard evaluation result.

    Attributes
    ----------
    smiles : str
        Canonical SMILES of the compound
    formula : str
        Molecular formula
    molecular_weight : float
        Molecular weight in g/mol
    delta_hf_kJ_mol : float
        Formation enthalpy in kJ/mol
    max_decomposition_kJ_mol : float
        Maximum heat of decomposition in kJ/mol
    max_decomposition_cal_g : float
        Maximum heat of decomposition in cal/g
    oxygen_balance_percent : float
        Oxygen balance in percent
    hazard_class : HazardClass
        Overall hazard classification (HIGH, MEDIUM, LOW)
    triggered_criteria : tuple[int, ...]
        List of triggered CHETAH criteria numbers (1-4)
    functional_group_alerts : tuple[str, ...]
        Names of detected reactive functional groups
    product_breakdown : dict[str, float]
        Decomposition product amounts
    gas_volume_L_g : float
        Gas generation at STP in L/g
    method : str
        Calculation method used
    """

    smiles: str
    formula: str
    molecular_weight: float
    delta_hf_kJ_mol: float
    max_decomposition_kJ_mol: float
    max_decomposition_cal_g: float
    oxygen_balance_percent: float
    hazard_class: HazardClass
    triggered_criteria: tuple[int, ...]
    functional_group_alerts: tuple[str, ...]
    product_breakdown: dict[str, float]
    gas_volume_L_g: float
    method: str = "CHETAH_analytical"


def evaluate_hazard(compound: Compound) -> HazardResult:
    """
    Perform complete hazard evaluation for a compound.

    Parameters
    ----------
    compound : Compound
        Compound to evaluate

    Returns
    -------
    HazardResult
        Complete hazard assessment
    """
    # Calculate decomposition
    decomp: DecompositionResult = calculate_max_decomposition(compound)

    # Calculate oxygen balance
    ob_percent = calculate_oxygen_balance(compound.composition, compound.molecular_weight)

    # Detect functional groups
    alerts = get_alert_names(compound)

    # Evaluate criteria
    triggered: list[int] = []

    # Criterion 1: ΔHd < -300 cal/g
    if decomp.delta_hd_cal_g < CRITERION_1_THRESHOLD:
        triggered.append(1)

    # Criterion 2: ΔHd < -100 cal/g (only if criterion 1 not triggered)
    elif decomp.delta_hd_cal_g < CRITERION_2_THRESHOLD:
        triggered.append(2)

    # Criterion 3: OB in reactive range
    if CRITERION_3_OB_LOW < ob_percent < CRITERION_3_OB_HIGH:
        triggered.append(3)

    # Criterion 4: Functional group alerts
    if alerts:
        triggered.append(4)

    # Determine hazard class
    hazard_class: HazardClass
    if 1 in triggered:
        hazard_class = "HIGH"
    elif 2 in triggered or (3 in triggered and 4 in triggered):
        hazard_class = "MEDIUM"
    elif triggered:
        hazard_class = "MEDIUM"  # Any criterion triggered = at least medium
    else:
        hazard_class = "LOW"

    return HazardResult(
        smiles=compound.canonical_smiles,
        formula=compound.formula,
        molecular_weight=compound.molecular_weight,
        delta_hf_kJ_mol=compound.delta_hf_kJ_mol,
        max_decomposition_kJ_mol=decomp.delta_hd_kJ_mol,
        max_decomposition_cal_g=decomp.delta_hd_cal_g,
        oxygen_balance_percent=ob_percent,
        hazard_class=hazard_class,
        triggered_criteria=tuple(sorted(triggered)),
        functional_group_alerts=tuple(alerts),
        product_breakdown=decomp.products,
        gas_volume_L_g=decomp.gas_volume_L_g,
        method="CHETAH_analytical",
    )
