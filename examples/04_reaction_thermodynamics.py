#!/usr/bin/env python3
"""
Reaction Thermodynamics Example
===============================

This example demonstrates temperature-dependent reaction thermodynamics
for process safety evaluation.

You'll learn how to:
- Create balanced reactions from formulas
- Calculate enthalpy of reaction (ΔH_rxn)
- Use temperature-dependent thermodynamics
- Evaluate Gibbs energy for spontaneity
- Assess heat removal requirements

Target audience: Process chemists, reaction engineers

Chemical Background:
  ΔG = ΔH - TΔS
  - ΔG < 0: Spontaneous (thermodynamically favorable)
  - ΔG > 0: Non-spontaneous
  - ΔH < 0: Exothermic (releases heat)
  - ΔH > 0: Endothermic (absorbs heat)
"""

from phoenix import Compound, Reaction, Auto
from phoenix.thermo import get_formation_enthalpy

# =============================================================================
# Example 1: Combustion of Methane
# =============================================================================

print("=" * 70)
print("Reaction Thermodynamics Analysis")
print("=" * 70)

print("\n--- Example 1: Methane Combustion ---")
print("CH4 + 2 O2 → CO2 + 2 H2O")

# Create reaction with auto-balanced stoichiometry
# Auto means "figure out the coefficient"
rxn = Reaction.from_smiles(
    reactants=[("C", 1), ("O=O", Auto)],  # CH4 = C (with H added), O2 = O=O
    products=[("O=C=O", Auto), ("O", Auto)],  # CO2, H2O
)

print(f"\nBalanced equation: {rxn}")
print(f"Enthalpy of reaction: ΔH_rxn = {rxn.delta_h:.1f} kJ/mol")

# Negative ΔH means exothermic
if rxn.delta_h < 0:
    print("  → Exothermic (releases heat)")
else:
    print("  → Endothermic (absorbs heat)")

# Show NIST reference data used for products
print(f"\n  NIST Reference Data (products):")
for formula in ["CO2", "H2O"]:
    hf = get_formation_enthalpy(formula)
    gas_val = hf.gas_kJ_mol
    print(f"    {formula}: ΔHf° = {gas_val:.2f} kJ/mol ({hf.source})")

# =============================================================================
# Example 2: Temperature-Dependent Properties
# =============================================================================

print("\n" + "=" * 70)
print("\n--- Example 2: Temperature Dependence ---")

# Analyze ethanol properties across temperature range
ethanol = Compound.from_smiles("CCO")

print(f"\nEthanol thermodynamic properties vs temperature:")
print(f"{'T (K)':<10} {'H (kJ/mol)':<15} {'S (J/mol·K)':<15} {'G (kJ/mol)':<15}")
print("-" * 55)

for T in [298.15, 350, 400, 450, 500]:
    state = ethanol.thermo_at(T=T)
    print(f"{T:<10.1f} {state.H.value:<15.1f} {state.S.value:<15.1f} {state.G.value:<15.1f}")

print("\nNote: G = H - TS (Gibbs free energy)")
print("      G becomes more negative at higher T due to entropy term")

# =============================================================================
# Example 3: Reaction Enthalpy from Components
# =============================================================================

print("\n" + "=" * 70)
print("\n--- Example 3: Hydrogenation Reaction ---")
print("Benzene + 3 H2 → Cyclohexane")

# Create compounds
benzene = Compound.from_smiles("c1ccccc1")
hydrogen = Compound.from_smiles("[H][H]")
cyclohexane = Compound.from_smiles("C1CCCCC1")

# Calculate reaction enthalpy manually
# ΔH_rxn = Σ ΔHf(products) - Σ ΔHf(reactants)
hf_benzene = benzene.enthalpy_of_formation.value
hf_h2 = hydrogen.enthalpy_of_formation.value
hf_cyclohexane = cyclohexane.enthalpy_of_formation.value

delta_h_rxn = hf_cyclohexane - (hf_benzene + 3 * hf_h2)

print(f"\nFormation enthalpies (kJ/mol):")
print(f"  Benzene:     ΔHf° = {hf_benzene:.1f}")
print(f"  H2:          ΔHf° = {hf_h2:.1f}")
print(f"  Cyclohexane: ΔHf° = {hf_cyclohexane:.1f}")

print(f"\nReaction enthalpy:")
print(f"  ΔH_rxn = {hf_cyclohexane:.1f} - ({hf_benzene:.1f} + 3×{hf_h2:.1f})")
print(f"  ΔH_rxn = {delta_h_rxn:.1f} kJ/mol")

if delta_h_rxn < 0:
    print("  → Exothermic hydrogenation")
    print(f"  → Heat removal needed: {abs(delta_h_rxn):.1f} kJ per mol benzene")

# =============================================================================
# Example 4: Process Safety Assessment
# =============================================================================

print("\n" + "=" * 70)
print("\n--- Example 4: Process Safety Assessment ---")

# For a batch reactor processing 100 mol of benzene
moles_benzene = 100
heat_released = abs(delta_h_rxn) * moles_benzene

print(f"\nBatch hydrogenation of {moles_benzene} mol benzene:")
print(f"  Total heat released: {heat_released:.0f} kJ")
print(f"  = {heat_released/4.184:.0f} kcal")
print(f"  = {heat_released/3600:.2f} kWh")

# Adiabatic temperature rise estimate (rough)
# Assuming Cp ~ 150 J/(mol·K) for organic liquids
cp_estimate = 150  # J/(mol·K)
delta_T_adiabatic = (delta_h_rxn * 1000) / cp_estimate

print(f"\nAdiabatic temperature rise (estimate):")
print(f"  ΔT_ad ≈ {abs(delta_T_adiabatic):.0f} K")
print(f"  → Significant cooling required to prevent runaway!")

# =============================================================================
# Example 5: NIST Reference Data Table
# =============================================================================

print("\n" + "=" * 70)
print("\n--- Example 5: NIST Reference Data ---")
print("Formation enthalpies for common reaction products (298.15 K):\n")

print(f"{'Formula':<8} {'ΔHf°(g) kJ/mol':>15} {'ΔHf°(l) kJ/mol':>15} {'Source':<25}")
print("-" * 65)

reference_compounds = ["CO2", "CO", "H2O", "H2", "N2", "O2", "CH4", "NH3"]
for formula in reference_compounds:
    hf = get_formation_enthalpy(formula)
    gas_str = f"{hf.gas_kJ_mol:.2f}" if hf.gas_kJ_mol is not None else "N/A"
    liq_str = f"{hf.liquid_kJ_mol:.2f}" if hf.liquid_kJ_mol is not None else "N/A"
    print(f"{formula:<8} {gas_str:>15} {liq_str:>15} {hf.source:<25}")

print("\nNote: These values are from NIST-JANAF Thermochemical Tables")
print("      or the chemicals library (CalebBell/ChEDL).")

# =============================================================================
# Key Takeaways
# =============================================================================

print("\n" + "=" * 70)
print("Key Takeaways:")
print("  1. Reaction.from_smiles() auto-balances stoichiometry")
print("  2. compound.thermo_at(T=...) gives T-dependent properties")
print("  3. ΔH_rxn = Σ ΔHf(products) - Σ ΔHf(reactants)")
print("  4. Negative ΔH = exothermic = heat removal needed")
print("  5. Adiabatic ΔT helps assess runaway risk")
print("  6. get_formation_enthalpy() provides NIST reference data")
print("=" * 70)
