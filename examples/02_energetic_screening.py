#!/usr/bin/env python3
"""
Energetic Material Screening Example
=====================================

This example demonstrates batch hazard screening - a common workflow
for safety engineers evaluating multiple candidate compounds.

You'll learn how to:
- Screen multiple compounds at once using screen()
- Analyze results with pandas DataFrame
- Filter compounds by hazard class
- Handle invalid SMILES gracefully
- Export results to CSV

Target audience: Safety engineers, process chemists evaluating libraries
"""

import pandas as pd

from phoenix import screen

# =============================================================================
# Define Compounds to Screen
# =============================================================================

# A mix of well-known energetic materials and some common compounds
# SMILES from PubChem/ChemSpider
compounds = {
    # High explosives
    "TNT": "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]",
    "RDX": "O=[N+]([O-])N1CN([N+](=O)[O-])CN([N+](=O)[O-])C1",
    "PETN": "O=[N+]([O-])OCC(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]",
    "Nitroglycerin": "O=[N+]([O-])OCC(O[N+](=O)[O-])CO[N+](=O)[O-]",
    # Mildly energetic
    "Nitrobenzene": "c1ccc([N+](=O)[O-])cc1",
    "AIBN": "CC(C)(C#N)N=NC(C)(C)C#N",  # Azo initiator
    # Common (low hazard)
    "Ethanol": "CCO",
    "Acetone": "CC(=O)C",
    "Toluene": "Cc1ccccc1",
    # Invalid SMILES (to demonstrate error handling)
    "Invalid": "not_a_smiles",
}

print("=" * 70)
print("Energetic Material Batch Screening")
print("=" * 70)

# =============================================================================
# Run Batch Screening
# =============================================================================

# screen() processes multiple compounds and returns a BatchResult
smiles_list = list(compounds.values())
names = list(compounds.keys())

print(f"\nScreening {len(smiles_list)} compounds...")
results = screen(smiles_list)

# =============================================================================
# Analyze Results
# =============================================================================

# BatchResult contains a DataFrame with all results
df = results.dataframe

# Add compound names for readability
df.insert(0, "name", names)

print(f"\n--- Screening Summary ---")
print(f"Total compounds:  {len(smiles_list)}")
print(f"Successfully processed: {results.successful}")
print(f"Failed: {results.failed}")

# =============================================================================
# Display Results Table
# =============================================================================

print(f"\n--- Results ---")
print(f"{'Name':<15} {'Formula':<12} {'OB%':>8} {'ΔHd kJ/mol':>12} {'Hazard':>10}")
print("-" * 60)

for idx, row in df.iterrows():
    name = names[idx][:14]
    formula = str(row.get("formula", "?"))[:11] if pd.notna(row.get("formula")) else "N/A"
    ob = row.get("ob_percent", float("nan"))
    dhd = row.get("delta_hd_kJ_mol", float("nan"))
    hazard = str(row.get("hazard_class", "N/A")) if pd.notna(row.get("hazard_class")) else "ERROR"

    # Format with handling for NaN values
    ob_str = f"{ob:8.1f}" if pd.notna(ob) else "     N/A"
    dhd_str = f"{dhd:12.1f}" if pd.notna(dhd) else "         N/A"

    print(f"{name:<15} {formula:<12} {ob_str} {dhd_str} {hazard:>10}")

# =============================================================================
# Filter by Hazard Class
# =============================================================================

print(f"\n--- High Hazard Compounds ---")
high_hazard = df[df["hazard_class"] == "HIGH"]
if len(high_hazard) > 0:
    for idx, row in high_hazard.iterrows():
        print(f"  - {names[idx]}: ΔHd = {row['delta_hd_kJ_mol']:.1f} kJ/mol")
else:
    print("  None found")

print(f"\n--- Medium Hazard Compounds ---")
medium_hazard = df[df["hazard_class"] == "MEDIUM"]
if len(medium_hazard) > 0:
    for idx, row in medium_hazard.iterrows():
        print(f"  - {names[idx]}: ΔHd = {row['delta_hd_kJ_mol']:.1f} kJ/mol")
else:
    print("  None found")

# =============================================================================
# Export to CSV
# =============================================================================

# Uncomment to save results:
# df.to_csv("hazard_screening_results.csv", index=False)
# print("\nResults saved to hazard_screening_results.csv")

# =============================================================================
# Key Takeaways
# =============================================================================

print("\n" + "=" * 70)
print("Key Takeaways:")
print("  1. screen() processes multiple SMILES and returns BatchResult")
print("  2. results.dataframe gives pandas DataFrame for analysis")
print("  3. Invalid SMILES are handled gracefully (NaN values)")
print("  4. Oxygen balance near 0% = most dangerous (self-oxidizing)")
print("  5. More negative ΔHd = more energy released on decomposition")
print("=" * 70)
