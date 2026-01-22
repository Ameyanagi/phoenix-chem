#!/usr/bin/env python3
"""
Decomposition Pathway Analysis Example
======================================

This example analyzes what happens when energetic compounds decompose:
- What products form (N2, CO2, CO, H2O, C(s), etc.)
- How much energy is released
- How much gas is generated (important for vent sizing)

You'll learn how to:
- Use max_decomposition() with hierarchy and LP methods
- Interpret product distributions
- Calculate gas generation at different temperatures
- Compare analytical vs optimized decomposition

Target audience: Process safety engineers, explosion hazard analysts

Chemical Background:
  The CHETAH hierarchy determines products in priority order:
  1. HF (strongest H-X bond)
  2. N2 (always forms)
  3. P4O10, H2O, CO2, SO2, CO, HCl, HBr
  4. C(s), H2, S(s), excess halogens, O2

  Oxygen balance determines whether you get:
  - CO2 (enough oxygen)
  - CO (limited oxygen)
  - C(s) graphite (very oxygen-deficient)
"""

from phoenix import Compound

# =============================================================================
# Example 1: TNT (Oxygen-Deficient)
# =============================================================================

print("=" * 70)
print("Decomposition Pathway Analysis")
print("=" * 70)

# TNT: C7H5N3O6 - severely oxygen-deficient (OB = -74%)
tnt = Compound.from_smiles("Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]")

print(f"\n--- TNT (C7H5N3O6) ---")
print(f"Oxygen Balance: {tnt.oxygen_balance:.1f}%")
print("(Negative = oxygen-deficient, needs external O2 for complete combustion)")

# Calculate decomposition using analytical hierarchy
result = tnt.max_decomposition(method="hierarchy")

print(f"\nDecomposition Products (CHETAH Hierarchy):")
for product, moles in sorted(result.products.items(), key=lambda x: -x[1]):
    if moles > 0.01:
        print(f"  {product:>6}: {moles:.2f} mol")

print(f"\nEnergy Release:")
print(f"  ΔHd = {result.delta_hd_kJ_mol:.1f} kJ/mol")
print(f"  ΔHd = {result.delta_hd_cal_g:.1f} cal/g")

# Note: TNT is oxygen-deficient, so it produces C(s) and limited CO2 instead of all CO2
print(f"\nInterpretation:")
print(f"  - TNT is severely oxygen-deficient (OB = {tnt.oxygen_balance:.0f}%)")
print(f"  - Not enough O to convert all C to CO2")
print(f"  - Most carbon remains as solid C (graphite/soot)")

# =============================================================================
# Example 2: Nitroglycerin (Oxygen-Balanced)
# =============================================================================

print("\n" + "=" * 70)

# Nitroglycerin: C3H5N3O9 - nearly oxygen-balanced (OB = +3.5%)
ng = Compound.from_smiles("O=[N+]([O-])OCC(O[N+](=O)[O-])CO[N+](=O)[O-]")

print(f"\n--- Nitroglycerin (C3H5N3O9) ---")
print(f"Oxygen Balance: {ng.oxygen_balance:.1f}%")
print("(Near zero = self-oxidizing, most energetic)")

result_ng = ng.max_decomposition(method="hierarchy")

print(f"\nDecomposition Products (CHETAH Hierarchy):")
for product, moles in sorted(result_ng.products.items(), key=lambda x: -x[1]):
    if moles > 0.01:
        print(f"  {product:>6}: {moles:.2f} mol")

print(f"\nEnergy Release:")
print(f"  ΔHd = {result_ng.delta_hd_kJ_mol:.1f} kJ/mol")
print(f"  ΔHd = {result_ng.delta_hd_cal_g:.1f} cal/g")

# Nitroglycerin is oxygen-balanced, so it produces mostly CO2 and H2O
print(f"\nInterpretation:")
print(f"  - Nitroglycerin is nearly oxygen-balanced (OB = {ng.oxygen_balance:.0f}%)")
print(f"  - Enough O to oxidize all C to CO2 and H to H2O")
print(f"  - No solid carbon - all products are gaseous")

# =============================================================================
# Compare Hierarchy vs LP Methods
# =============================================================================

print("\n" + "=" * 70)
print("\n--- Method Comparison: Hierarchy vs LP ---")

# Use method="both" to compare
comparison = tnt.max_decomposition(method="both")

print(f"\nTNT Decomposition:")
print(f"  Hierarchy ΔHd: {comparison.hierarchy_delta_hd:.1f} kJ/mol")
print(f"  LP Optimized:  {comparison.lp_delta_hd:.1f} kJ/mol")
print(f"  Deviation:     {comparison.deviation_percent:.2f}%")

print(f"\nNote: LP optimization finds the mathematically optimal product")
print(f"distribution. For standard compounds, it should match the hierarchy.")

# =============================================================================
# Gas Generation Analysis
# =============================================================================

print("\n" + "=" * 70)
print("\n--- Gas Generation (for Vent Sizing) ---")

# Gas volume depends on temperature (PV = nRT)
# Compare at different temperatures

print(f"\nTNT Gas Generation:")
for temp in [298.15, 500, 1000]:
    result_t = tnt.max_decomposition(gas_temperature_K=temp)
    print(f"  T = {temp:6.1f} K: {result_t.gas_volume_L_g:.3f} L/g")
    print(f"                   {result_t.gas_moles:.2f} mol gas/mol TNT")

print(f"\nGas Composition (mole fractions):")
for gas, fraction in sorted(result.gas_composition.items(), key=lambda x: -x[1]):
    print(f"  {gas:>6}: {fraction*100:.1f}%")

# =============================================================================
# Key Takeaways
# =============================================================================

print("\n" + "=" * 70)
print("Key Takeaways:")
print("  1. Oxygen balance determines product distribution:")
print("     - OB >> 0: Excess O2 released")
print("     - OB ~ 0:  CO2 + H2O (most energetic)")
print("     - OB << 0: C(s) + limited CO2 (soot formation)")
print("  2. Gas volume scales with temperature (PV = nRT)")
print("  3. LP and hierarchy methods should agree for typical compounds")
print("  4. Gas generation is critical for vent sizing calculations")
print("=" * 70)
