"""Regression tests: small molecules with no Benson thermochem groups must
fall back to the chemicals library instead of silently reporting S = 0.

See: Sabatier reaction (CO2 + 4H2 -> CH4 + 2H2O) previously returned
delta_S = +377 J/(mol K) instead of the literature -172 J/(mol K).
"""

from __future__ import annotations

import pytest

from phoenix import Auto, Compound, Reaction

SMALL_MOLECULES = {
    "O=C=O": 213.8,  # CO2
    "[H][H]": 130.7,  # H2
    "C": 186.3,  # CH4
}


class TestEntropyFallback:
    @pytest.mark.parametrize(("smiles", "expected_s"), SMALL_MOLECULES.items())
    def test_small_molecule_entropy_nonzero(self, smiles: str, expected_s: float):
        compound = Compound.from_smiles(smiles)
        s = compound.thermo_at(T=298.15).S
        assert s.value == pytest.approx(expected_s, abs=1.0)
        assert s.source == "chemicals (Lookup)"

    def test_sabatier_reaction_thermodynamics(self):
        rxn = Reaction.from_smiles(
            reactants=[("O=C=O", 1), ("[H][H]", Auto)],
            products=[("C", Auto), ("O", Auto)],
        ).balance()
        assert rxn.delta_h == pytest.approx(-165, abs=1)
        assert rxn.delta_s == pytest.approx(-172, abs=3)
        assert rxn.delta_g == pytest.approx(-113, abs=2)

    def test_small_molecule_cp_nonzero(self):
        for smiles in SMALL_MOLECULES:
            compound = Compound.from_smiles(smiles)
            assert compound.thermo_at(T=298.15).Cp.value > 0
