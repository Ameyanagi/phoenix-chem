#!/usr/bin/env python3
"""
PHOENIX Quick Start Example
===========================

This example introduces the basic PHOENIX API in under 5 minutes.
You'll learn how to:
- Create compounds from SMILES strings
- Access molecular properties (formula, MW, composition)
- Calculate thermodynamic properties
- Evaluate hazard classification

Target audience: New users getting started with PHOENIX
"""

from phoenix import Compound
from phoenix.thermo import get_formation_enthalpy, check_data_sources

# =============================================================================
# Creating Compounds
# =============================================================================

# Create a compound from a SMILES string
# Ethanol is a common, safe compound - good for testing
ethanol = Compound.from_smiles("CCO")

print("=" * 60)
print("PHOENIX Quick Start")
print("=" * 60)

# Basic molecular properties
print(f"\n--- Ethanol (CCO) ---")
print(f"Formula:          {ethanol.formula}")
print(f"Molecular weight: {ethanol.molecular_weight:.2f} g/mol")
print(f"Composition:      {ethanol.composition}")

# =============================================================================
# Thermodynamic Properties
# =============================================================================

# Enthalpy of formation (estimated via Benson Group Additivity)
hf = ethanol.enthalpy_of_formation
print(f"\nThermodynamics (298.15 K):")
print(f"  ΔHf° = {hf.value:.1f} {hf.unit}")

# Entropy
s = ethanol.entropy
print(f"  S°   = {s.value:.1f} {s.unit}")

# =============================================================================
# Hazard Indicators
# =============================================================================

# Oxygen balance - key indicator for explosive potential
# Negative = oxygen-deficient (needs external O2 for complete combustion)
# Positive = oxygen-rich (self-oxidizing)
# Near zero = most dangerous for explosives
print(f"\nHazard Indicators:")
print(f"  Oxygen Balance: {ethanol.oxygen_balance:.1f}%")

# =============================================================================
# Full Hazard Evaluation
# =============================================================================

# evaluate_hazard() runs the complete CHETAH-style assessment
result = ethanol.evaluate_hazard()
print(f"\nHazard Classification:")
print(f"  Class: {result.hazard_class}")
print(f"  Max ΔHd: {result.max_decomposition_kJ_mol:.1f} kJ/mol")

# =============================================================================
# Now try a more energetic compound: Nitrobenzene
# =============================================================================

print("\n" + "=" * 60)
nitrobenzene = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")

print(f"\n--- Nitrobenzene ---")
print(f"Formula:          {nitrobenzene.formula}")
print(f"Molecular weight: {nitrobenzene.molecular_weight:.2f} g/mol")

# Note the more negative oxygen balance - very oxygen-deficient
print(f"\nHazard Indicators:")
print(f"  Oxygen Balance: {nitrobenzene.oxygen_balance:.1f}%")

# Nitrobenzene should show higher hazard than ethanol
result_nb = nitrobenzene.evaluate_hazard()
print(f"\nHazard Classification:")
print(f"  Class: {result_nb.hazard_class}")
print(f"  Max ΔHd: {result_nb.max_decomposition_kJ_mol:.1f} kJ/mol")

# =============================================================================
# NIST Reference Data Comparison
# =============================================================================

print("\n" + "=" * 60)
print("\n--- NIST Reference Data ---")

# Check available data sources
sources = check_data_sources()
print(f"\nData sources available:")
for source, available in sources.items():
    status = "✓" if available else "✗"
    print(f"  {status} {source}")

# Compare Benson GA estimates vs NIST reference for common products
print(f"\nReference data for decomposition products:")
print(f"{'Formula':<10} {'ΔHf° (kJ/mol)':>15} {'Source':<30}")
print("-" * 55)

for formula in ["CO2", "H2O", "CO", "N2", "CH4"]:
    hf = get_formation_enthalpy(formula)
    value = hf.gas_kJ_mol
    if value is not None:
        print(f"{formula:<10} {value:>15.2f} {hf.source:<30}")

print("\nNote: PHOENIX uses NIST reference data for product enthalpies")
print("      and Benson Group Additivity for compound estimation.")

# =============================================================================
# Key Takeaways
# =============================================================================

print("\n" + "=" * 60)
print("Key Takeaways:")
print("  1. Compound.from_smiles() creates compounds from SMILES")
print("  2. .oxygen_balance gives a quick hazard indicator")
print("  3. .evaluate_hazard() runs full CHETAH assessment")
print("  4. More negative ΔHd = more energy release on decomposition")
print("  5. get_formation_enthalpy() provides NIST reference data")
print("=" * 60)
