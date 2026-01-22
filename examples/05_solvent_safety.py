#!/usr/bin/env python3
"""
Solvent Safety Check Example
============================

This example demonstrates quick safety evaluation for solvents and
additives - a common workflow when selecting compounds for synthesis
or formulation.

You'll learn how to:
- Evaluate single compounds quickly
- Detect hazardous functional groups
- Interpret oxygen balance for mixtures
- Create a "traffic light" safety assessment

Target audience: Medicinal chemists, formulation scientists

Safety Note:
  Always verify results experimentally. This tool is for screening only.
  Peroxides and other hazardous compounds require special handling.
"""

from phoenix import Compound
from phoenix.hazard import detect_functional_groups

# =============================================================================
# Define Solvents to Evaluate
# =============================================================================

solvents = [
    ("Ethanol", "CCO"),
    ("Acetone", "CC(=O)C"),
    ("Tetrahydrofuran", "C1CCOC1"),
    ("Diethyl ether", "CCOCC"),
    ("Dimethyl sulfoxide", "CS(=O)C"),
    ("Acetonitrile", "CC#N"),
    # Potentially hazardous
    ("Nitromethane", "C[N+](=O)[O-]"),
]

print("=" * 70)
print("Solvent Safety Screening")
print("=" * 70)

# =============================================================================
# Safety Assessment Function
# =============================================================================

def assess_safety(compound: Compound) -> tuple[str, str]:
    """
    Simple traffic-light safety assessment.

    Returns (color, reason) where color is GREEN/YELLOW/RED.
    """
    hazard = compound.evaluate_hazard()

    # Check hazard class
    if hazard.hazard_class == "HIGH":
        return "RED", "High hazard classification"

    if hazard.hazard_class == "MEDIUM":
        return "YELLOW", "Medium hazard - use caution"

    # Check oxygen balance (near zero is concerning)
    ob = compound.oxygen_balance
    if -10 < ob < 10:
        return "YELLOW", f"Near-zero oxygen balance ({ob:.1f}%)"

    # Check functional groups
    groups = detect_functional_groups(compound)
    dangerous_group_names = {"nitro", "peroxide", "azide", "nitroso"}
    found_dangerous = [g.name for g in groups if g.name in dangerous_group_names]
    if found_dangerous:
        return "YELLOW", f"Contains {', '.join(found_dangerous)} groups"

    return "GREEN", "Standard handling"

# =============================================================================
# Evaluate Each Solvent
# =============================================================================

print(f"\n{'Solvent':<20} {'Formula':<12} {'MW':>8} {'OB%':>8} {'Status':<8} {'Notes'}")
print("-" * 85)

for name, smiles in solvents:
    try:
        compound = Compound.from_smiles(smiles)
        status, reason = assess_safety(compound)

        # Format status with indicator
        status_display = {
            "GREEN": "[OK]   ",
            "YELLOW": "[WARN] ",
            "RED": "[STOP] ",
        }[status]

        print(
            f"{name:<20} "
            f"{compound.formula:<12} "
            f"{compound.molecular_weight:>8.1f} "
            f"{compound.oxygen_balance:>8.1f} "
            f"{status_display}"
            f"{reason}"
        )

    except Exception as e:
        print(f"{name:<20} {'ERROR':<12} {'':<8} {'':<8} [ERR]   {e}")

# =============================================================================
# Detailed Analysis of Concerning Compounds
# =============================================================================

print("\n" + "=" * 70)
print("\n--- Detailed Analysis: Nitromethane ---")

nitromethane = Compound.from_smiles("C[N+](=O)[O-]")
print(f"\nFormula: {nitromethane.formula}")
print(f"Molecular weight: {nitromethane.molecular_weight:.2f} g/mol")

# Thermodynamics
print(f"\nThermodynamics:")
print(f"  ΔHf° = {nitromethane.enthalpy_of_formation.value:.1f} kJ/mol")
print(f"  S°   = {nitromethane.entropy.value:.1f} J/(mol·K)")

# Hazard indicators
print(f"\nHazard Indicators:")
print(f"  Oxygen Balance: {nitromethane.oxygen_balance:.1f}%")

# Decomposition
decomp = nitromethane.max_decomposition()
print(f"  Max ΔHd: {decomp.delta_hd_kJ_mol:.1f} kJ/mol")
print(f"  Gas generation: {decomp.gas_volume_L_g:.3f} L/g at 298 K")

# Functional groups
groups = detect_functional_groups(nitromethane)
group_names = [g.name for g in groups]
print(f"\nFunctional Groups: {', '.join(group_names) if group_names else 'None detected'}")

# Classification
hazard = nitromethane.evaluate_hazard()
print(f"\nHazard Classification: {hazard.hazard_class}")

print("\nSafety Notes for Nitromethane:")
print("  - Used as racing fuel and solvent")
print("  - Can detonate under confinement with heat/shock")
print("  - Requires special storage and handling")

# =============================================================================
# Key Takeaways
# =============================================================================

print("\n" + "=" * 70)
print("Key Takeaways:")
print("  1. Quick evaluation: Compound.from_smiles() + evaluate_hazard()")
print("  2. detect_functional_groups() identifies structural hazards")
print("  3. Oxygen balance near 0% indicates self-oxidizing potential")
print("  4. Traffic light system: GREEN/YELLOW/RED for quick decisions")
print("  5. Always verify with experimental data before scale-up")
print("=" * 70)
