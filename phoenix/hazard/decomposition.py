"""
Maximum heat of decomposition (Max ΔHd) calculation.

Implements both:
1. CHETAH analytical hierarchy - priority-based product selection
2. LP optimization - Linear Programming for mathematically optimal solution

Reference: ASTM E659 - Standard Test Method for Determining CHETAH

Thermodynamic Hierarchy (priority order):
1. F + H → HF (most exothermic H-X bond)
2. N → ½N₂ (always forms N₂)
3. P + O → ¼P₄O₁₀ (highly exothermic)
4. H + O → ½H₂O (water formation)
5. C + O → ½CO₂ (full oxidation if O available)
6. S + O → SO₂ (sulfur oxidation)
7. C + ½O → CO (partial oxidation if O limited)
8. Cl + H → HCl (after water)
9. Br + H → HBr (after HCl)
10. C → C(s) graphite (O-deficient)
11. H → ½H₂ (excess hydrogen)
12. S → S(s) (excess sulfur)
13. Excess halogens → X₂
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np
from scipy.optimize import linprog

from phoenix.exceptions import DecompositionError
from phoenix.thermo.data import get_formation_enthalpy

if TYPE_CHECKING:
    from phoenix.core.compound import Compound

# Ideal gas constant in L·atm/(mol·K)
R_GAS = 0.08206

# Default gas temperature (K)
DEFAULT_GAS_TEMP = 298.15

# Elements tracked for atom balance
ELEMENTS = ["C", "H", "N", "O", "S", "P", "F", "Cl", "Br"]

# Decomposition products with formation enthalpies and atom counts
# is_gas indicates if product is gaseous at standard conditions
DECOMPOSITION_PRODUCTS: dict[str, dict] = {
    "HF": {"delta_hf_kJ_mol": -273.30, "atoms": {"H": 1, "F": 1}, "is_gas": True},
    "N2": {"delta_hf_kJ_mol": 0.0, "atoms": {"N": 2}, "is_gas": True},
    "P4O10": {"delta_hf_kJ_mol": -2984.0, "atoms": {"P": 4, "O": 10}, "is_gas": False},
    "H2O": {"delta_hf_kJ_mol": -241.83, "atoms": {"H": 2, "O": 1}, "is_gas": True},
    "CO2": {"delta_hf_kJ_mol": -393.52, "atoms": {"C": 1, "O": 2}, "is_gas": True},
    "SO2": {"delta_hf_kJ_mol": -296.81, "atoms": {"S": 1, "O": 2}, "is_gas": True},
    "CO": {"delta_hf_kJ_mol": -110.53, "atoms": {"C": 1, "O": 1}, "is_gas": True},
    "HCl": {"delta_hf_kJ_mol": -92.31, "atoms": {"H": 1, "Cl": 1}, "is_gas": True},
    "HBr": {"delta_hf_kJ_mol": -36.29, "atoms": {"H": 1, "Br": 1}, "is_gas": True},
    "C": {"delta_hf_kJ_mol": 0.0, "atoms": {"C": 1}, "is_gas": False},  # graphite
    "H2": {"delta_hf_kJ_mol": 0.0, "atoms": {"H": 2}, "is_gas": True},
    "S": {"delta_hf_kJ_mol": 0.0, "atoms": {"S": 1}, "is_gas": False},  # rhombic
    "F2": {"delta_hf_kJ_mol": 0.0, "atoms": {"F": 2}, "is_gas": True},
    "Cl2": {"delta_hf_kJ_mol": 0.0, "atoms": {"Cl": 2}, "is_gas": True},
    "Br2": {"delta_hf_kJ_mol": 30.91, "atoms": {"Br": 2}, "is_gas": True},
    "O2": {"delta_hf_kJ_mol": 0.0, "atoms": {"O": 2}, "is_gas": True},
    "P": {"delta_hf_kJ_mol": 0.0, "atoms": {"P": 1}, "is_gas": False},  # white P
}

# Legacy dict for backward compatibility
PRODUCT_HF: dict[str, float] = {
    name: data["delta_hf_kJ_mol"] for name, data in DECOMPOSITION_PRODUCTS.items()
}


@dataclass(frozen=True)
class DecompositionResult:
    """
    Result of maximum heat of decomposition calculation.

    Attributes
    ----------
    delta_hd_kJ_mol : float
        Maximum heat of decomposition in kJ/mol (negative = exothermic)
    delta_hd_cal_g : float
        Maximum heat of decomposition in cal/g
    products : dict[str, float]
        Molar amounts of each decomposition product
    reactant_hf_kJ_mol : float
        Formation enthalpy of reactant in kJ/mol
    products_hf_kJ_mol : float
        Total formation enthalpy of products in kJ/mol
    gas_volume_L_g : float
        Gas volume per gram at gas_temperature_K
    gas_moles : float
        Total moles of gas per mole of compound
    gas_composition : dict[str, float]
        Gaseous products with mole fractions (sum to 1.0)
    gas_temperature_K : float
        Temperature used for gas volume calculation
    method : str
        Calculation method used ('hierarchy' or 'lp')
    """

    delta_hd_kJ_mol: float
    delta_hd_cal_g: float
    products: dict[str, float]
    reactant_hf_kJ_mol: float
    products_hf_kJ_mol: float
    gas_volume_L_g: float
    gas_moles: float
    gas_composition: dict[str, float]
    gas_temperature_K: float = DEFAULT_GAS_TEMP
    method: str = "hierarchy"


@dataclass(frozen=True)
class DecompositionComparison:
    """
    Result of comparing hierarchy and LP decomposition methods.

    Attributes
    ----------
    hierarchy_result : DecompositionResult
        Result from analytical hierarchy method
    lp_result : DecompositionResult
        Result from LP optimization method
    deviation_percent : float
        Percent deviation between methods: |ΔHd_lp - ΔHd_hier| / |ΔHd_hier| × 100
    """

    hierarchy_result: DecompositionResult
    lp_result: DecompositionResult
    deviation_percent: float

    @property
    def hierarchy_delta_hd(self) -> float:
        """ΔHd from hierarchy method in kJ/mol."""
        return self.hierarchy_result.delta_hd_kJ_mol

    @property
    def lp_delta_hd(self) -> float:
        """ΔHd from LP method in kJ/mol."""
        return self.lp_result.delta_hd_kJ_mol


def calculate_max_decomposition(
    compound: Compound,
    *,
    method: Literal["hierarchy", "lp", "both"] = "hierarchy",
    gas_temperature_K: float = DEFAULT_GAS_TEMP,
) -> DecompositionResult | DecompositionComparison:
    """
    Calculate maximum heat of decomposition.

    Parameters
    ----------
    compound : Compound
        Compound to analyze
    method : {'hierarchy', 'lp', 'both'}, default 'hierarchy'
        Calculation method:
        - 'hierarchy': CHETAH analytical priority rules (default)
        - 'lp': Linear Programming optimization
        - 'both': Run both and compare results
    gas_temperature_K : float, default 298.15
        Temperature for gas volume calculation using PV=nRT

    Returns
    -------
    DecompositionResult or DecompositionComparison
        Decomposition result with ΔHd, products, and gas generation data.
        Returns DecompositionComparison when method='both'.

    Examples
    --------
    >>> compound = Compound.from_smiles('CC')
    >>> result = compound.max_decomposition()  # hierarchy (default)
    >>> result_lp = compound.max_decomposition(method='lp')
    >>> result_both = compound.max_decomposition(method='both')
    >>> print(result_both.deviation_percent)
    """
    if method == "both":
        hier_result = _calculate_hierarchy(compound, gas_temperature_K)
        lp_result = _calculate_lp(compound, gas_temperature_K)

        # Calculate deviation
        if abs(hier_result.delta_hd_kJ_mol) > 1e-6:
            deviation = (
                abs(lp_result.delta_hd_kJ_mol - hier_result.delta_hd_kJ_mol)
                / abs(hier_result.delta_hd_kJ_mol)
                * 100
            )
        else:
            deviation = 0.0

        return DecompositionComparison(
            hierarchy_result=hier_result,
            lp_result=lp_result,
            deviation_percent=deviation,
        )
    elif method == "lp":
        return _calculate_lp(compound, gas_temperature_K)
    else:  # hierarchy (default)
        return _calculate_hierarchy(compound, gas_temperature_K)


def _calculate_hierarchy(
    compound: Compound, gas_temperature_K: float
) -> DecompositionResult:
    """Calculate decomposition using CHETAH analytical hierarchy."""
    comp = compound.composition.copy()
    mw = compound.molecular_weight
    reactant_hf = compound.delta_hf_kJ_mol

    products = _apply_hierarchy(comp)
    products_hf = _calculate_products_enthalpy(products)

    # ΔHd = ΔHf(products) - ΔHf(reactant)
    delta_hd_kJ_mol = products_hf - reactant_hf

    # Convert to cal/g: (kJ/mol × 1000 J/kJ) / (4.184 J/cal × MW g/mol)
    delta_hd_cal_g = (delta_hd_kJ_mol * 1000) / (4.184 * mw)

    # Calculate gas data
    gas_moles, gas_composition = _calculate_gas_composition(products)
    gas_volume = _calculate_gas_volume(gas_moles, mw, gas_temperature_K)

    return DecompositionResult(
        delta_hd_kJ_mol=delta_hd_kJ_mol,
        delta_hd_cal_g=delta_hd_cal_g,
        products=products,
        reactant_hf_kJ_mol=reactant_hf,
        products_hf_kJ_mol=products_hf,
        gas_volume_L_g=gas_volume,
        gas_moles=gas_moles,
        gas_composition=gas_composition,
        gas_temperature_K=gas_temperature_K,
        method="hierarchy",
    )


def _calculate_lp(compound: Compound, gas_temperature_K: float) -> DecompositionResult:
    """
    Calculate decomposition using Linear Programming optimization.

    Minimizes total product formation enthalpy subject to atom balance constraints.
    """
    comp = compound.composition
    mw = compound.molecular_weight
    reactant_hf = compound.delta_hf_kJ_mol

    # Build LP problem
    product_names = list(DECOMPOSITION_PRODUCTS.keys())
    n_products = len(product_names)

    # Cost vector: formation enthalpies (minimize = most exothermic)
    c = np.array([DECOMPOSITION_PRODUCTS[p]["delta_hf_kJ_mol"] for p in product_names])

    # Build atom balance matrix
    A_eq, b_eq = _build_atom_matrix(product_names, comp)

    # Bounds: x_i >= 0
    bounds = [(0, None) for _ in range(n_products)]

    # Solve LP
    result = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")

    if not result.success:
        warnings.warn(
            f"LP solver failed: {result.message}. Falling back to hierarchy method.",
            UserWarning,
        )
        return _calculate_hierarchy(compound, gas_temperature_K)

    # Parse results - filter negligible amounts
    products: dict[str, float] = {}
    for i, name in enumerate(product_names):
        amount = result.x[i]
        if amount > 1e-6:
            products[name] = float(amount)

    products_hf = _calculate_products_enthalpy(products)
    delta_hd_kJ_mol = products_hf - reactant_hf
    delta_hd_cal_g = (delta_hd_kJ_mol * 1000) / (4.184 * mw)

    # Calculate gas data
    gas_moles, gas_composition = _calculate_gas_composition(products)
    gas_volume = _calculate_gas_volume(gas_moles, mw, gas_temperature_K)

    return DecompositionResult(
        delta_hd_kJ_mol=delta_hd_kJ_mol,
        delta_hd_cal_g=delta_hd_cal_g,
        products=products,
        reactant_hf_kJ_mol=reactant_hf,
        products_hf_kJ_mol=products_hf,
        gas_volume_L_g=gas_volume,
        gas_moles=gas_moles,
        gas_composition=gas_composition,
        gas_temperature_K=gas_temperature_K,
        method="lp",
    )


def _build_atom_matrix(
    product_names: list[str], composition: dict[str, int]
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build atom balance constraint matrix for LP.

    Returns A_eq and b_eq for: A_eq @ x = b_eq

    Important: We include ALL elements that appear in any product,
    not just those in the compound. This ensures we can't create
    atoms that don't exist in the reactant.
    """
    n_elements = len(ELEMENTS)
    n_products = len(product_names)

    A_eq = np.zeros((n_elements, n_products))
    b_eq = np.zeros(n_elements)

    for i, elem in enumerate(ELEMENTS):
        b_eq[i] = composition.get(elem, 0)
        for j, product in enumerate(product_names):
            atoms = DECOMPOSITION_PRODUCTS[product]["atoms"]
            A_eq[i, j] = atoms.get(elem, 0)

    return A_eq, b_eq


def _apply_hierarchy(comp: dict[str, int]) -> dict[str, float]:
    """
    Apply thermodynamic hierarchy to determine product distribution.

    Modifies comp in place (consumes atoms) and returns product amounts.
    """
    products: dict[str, float] = {}

    # Get available atoms
    n_c = comp.get("C", 0)
    n_h = comp.get("H", 0)
    n_n = comp.get("N", 0)
    n_o = comp.get("O", 0)
    n_s = comp.get("S", 0)
    n_p = comp.get("P", 0)
    n_f = comp.get("F", 0)
    n_cl = comp.get("Cl", 0)
    n_br = comp.get("Br", 0)

    # Priority 1: F + H → HF
    hf_formed = min(n_f, n_h)
    if hf_formed > 0:
        products["HF"] = hf_formed
        n_f -= hf_formed
        n_h -= hf_formed

    # Priority 2: N → ½N₂
    if n_n > 0:
        products["N2"] = n_n / 2.0
        n_n = 0

    # Priority 3: P + 2.5O → ¼P₄O₁₀
    if n_p > 0 and n_o >= 2.5 * n_p:
        products["P4O10"] = n_p / 4.0
        n_o -= int(2.5 * n_p)
        n_p = 0
    elif n_p > 0 and n_o > 0:
        # Partial phosphorus oxidation - use available O
        p_oxidized = min(n_p, n_o / 2.5)
        if p_oxidized > 0:
            products["P4O10"] = p_oxidized / 4.0
            n_o -= int(2.5 * p_oxidized)
            n_p -= int(p_oxidized)

    # Priority 4: H + ½O → ½H₂O (2H + O → H₂O)
    h2o_formed = min(n_h // 2, n_o)
    if h2o_formed > 0:
        products["H2O"] = h2o_formed
        n_h -= 2 * h2o_formed
        n_o -= h2o_formed

    # Priority 5: C + O → ½CO₂ (C + 2O → CO₂)
    co2_formed = min(n_c, n_o // 2)
    if co2_formed > 0:
        products["CO2"] = co2_formed
        n_c -= co2_formed
        n_o -= 2 * co2_formed

    # Priority 6: S + O → SO₂ (S + 2O → SO₂)
    so2_formed = min(n_s, n_o // 2)
    if so2_formed > 0:
        products["SO2"] = so2_formed
        n_s -= so2_formed
        n_o -= 2 * so2_formed

    # Priority 7: C + ½O → CO
    co_formed = min(n_c, n_o)
    if co_formed > 0:
        products["CO"] = co_formed
        n_c -= co_formed
        n_o -= co_formed

    # Priority 8: Cl + H → HCl
    hcl_formed = min(n_cl, n_h)
    if hcl_formed > 0:
        products["HCl"] = hcl_formed
        n_cl -= hcl_formed
        n_h -= hcl_formed

    # Priority 9: Br + H → HBr
    hbr_formed = min(n_br, n_h)
    if hbr_formed > 0:
        products["HBr"] = hbr_formed
        n_br -= hbr_formed
        n_h -= hbr_formed

    # Priority 10: C → C(s) graphite
    if n_c > 0:
        products["C"] = n_c
        n_c = 0

    # Priority 11: H → ½H₂
    if n_h > 0:
        products["H2"] = n_h / 2.0
        n_h = 0

    # Priority 12: S → S(s)
    if n_s > 0:
        products["S"] = n_s
        n_s = 0

    # Priority 13: Excess halogens → X₂
    if n_f > 0:
        products["F2"] = n_f / 2.0
    if n_cl > 0:
        products["Cl2"] = n_cl / 2.0
    if n_br > 0:
        products["Br2"] = n_br / 2.0

    # Excess oxygen
    if n_o > 0:
        products["O2"] = n_o / 2.0

    # Excess phosphorus (shouldn't happen normally)
    if n_p > 0:
        products["P"] = n_p

    return products


def _calculate_products_enthalpy(products: dict[str, float]) -> float:
    """Calculate total formation enthalpy of products."""
    total = 0.0
    for formula, moles in products.items():
        hf = PRODUCT_HF.get(formula)
        if hf is None:
            # Try to get from data module
            try:
                hf_data = get_formation_enthalpy(formula)
                hf = hf_data.gas_kJ_mol or hf_data.solid_kJ_mol or 0.0
            except Exception:
                hf = 0.0
        total += hf * moles
    return total


def _calculate_gas_composition(
    products: dict[str, float],
) -> tuple[float, dict[str, float]]:
    """
    Calculate gas moles and composition from product distribution.

    Returns
    -------
    tuple[float, dict[str, float]]
        (gas_moles, gas_composition) where gas_composition has mole fractions
    """
    gas_products: dict[str, float] = {}
    total_gas_moles = 0.0

    for formula, moles in products.items():
        if DECOMPOSITION_PRODUCTS.get(formula, {}).get("is_gas", False):
            gas_products[formula] = moles
            total_gas_moles += moles

    # Calculate mole fractions
    if total_gas_moles > 0:
        gas_composition = {
            formula: moles / total_gas_moles for formula, moles in gas_products.items()
        }
    else:
        gas_composition = {}

    return total_gas_moles, gas_composition


def _calculate_gas_volume(
    gas_moles: float, mw: float, temperature_K: float = DEFAULT_GAS_TEMP
) -> float:
    """
    Calculate gas volume per gram using ideal gas law PV=nRT.

    Parameters
    ----------
    gas_moles : float
        Total moles of gas per mole of compound
    mw : float
        Molecular weight in g/mol
    temperature_K : float
        Temperature in Kelvin (default 298.15 K)

    Returns
    -------
    float
        Gas volume in L/g at specified temperature and 1 atm
    """
    # V = nRT/P, at P = 1 atm: V = n × R × T
    # R = 0.08206 L·atm/(mol·K)
    molar_volume = R_GAS * temperature_K  # L/mol at 1 atm

    # Volume per gram = (moles gas × molar_volume) / MW
    volume_L_g = (gas_moles * molar_volume) / mw

    return volume_L_g
