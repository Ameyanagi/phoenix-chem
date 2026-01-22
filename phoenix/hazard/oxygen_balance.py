"""
Oxygen balance (OB%) calculation.

Oxygen balance measures the degree to which a compound can oxidize completely
to CO₂, H₂O, SO₂, P₄O₁₀, etc. using only its own oxygen content.

Formula for C_a H_b N_c O_d S_e P_f Cl_g Br_h:
    OB% = -1600/MW × (2a + b/2 + 2e + 2.5f - g/2 - h/2 - d)

Reference: Meyer et al., "Explosives" 6th Ed. (Wiley-VCH, 2007)
"""

from __future__ import annotations


def calculate_oxygen_balance(
    composition: dict[str, int],
    molecular_weight: float,
) -> float:
    """
    Calculate oxygen balance (OB%) for a compound.

    Parameters
    ----------
    composition : dict[str, int]
        Elemental composition as {element: count}
    molecular_weight : float
        Molecular weight in g/mol

    Returns
    -------
    float
        Oxygen balance in percent. Positive = oxygen excess, negative = oxygen deficiency.

    Examples
    --------
    >>> # TNT: C7H5N3O6, MW = 227.13
    >>> calculate_oxygen_balance({"C": 7, "H": 5, "N": 3, "O": 6}, 227.13)
    -73.96...
    >>> # Nitroglycerin: C3H5N3O9, MW = 227.09
    >>> calculate_oxygen_balance({"C": 3, "H": 5, "N": 3, "O": 9}, 227.09)
    3.52...
    """
    # Get atom counts (default to 0 if not present)
    a = composition.get("C", 0)  # Carbon
    b = composition.get("H", 0)  # Hydrogen
    # c = composition.get("N", 0)  # Nitrogen (not used in OB formula)
    d = composition.get("O", 0)  # Oxygen
    e = composition.get("S", 0)  # Sulfur
    f = composition.get("P", 0)  # Phosphorus
    g = composition.get("Cl", 0)  # Chlorine
    h = composition.get("Br", 0)  # Bromine

    # Oxygen required for complete oxidation minus oxygen available
    # C -> CO2 needs 2O per C
    # H -> H2O needs 0.5O per H
    # S -> SO2 needs 2O per S
    # P -> P4O10 needs 2.5O per P
    # Cl -> HCl releases 0.5O equivalent (consumes H that could form H2O)
    # Br -> HBr releases 0.5O equivalent
    oxygen_needed = 2 * a + 0.5 * b + 2 * e + 2.5 * f - 0.5 * g - 0.5 * h
    oxygen_balance = d - oxygen_needed

    # Convert to percentage
    # OB% = (O available - O needed) × 1600 / MW
    # The factor 1600 = 100% × (16 g/mol for O atom)
    ob_percent = (oxygen_balance * 1600) / molecular_weight

    return ob_percent
