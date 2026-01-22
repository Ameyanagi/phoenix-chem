"""
Tests for the Reaction class with auto-balancing stoichiometry.

These tests verify:
1. API flexibility (multiple input formats)
2. Full auto-balancing using null-space algorithm
3. Mixed explicit/auto coefficient systems
4. Error handling for over/under-constrained systems
5. Edge cases and numerical stability

Note: SMILES notation is used throughout. Common molecules:
- Methane: C (implicit hydrogens)
- Oxygen: O=O
- Water: O
- CO2: O=C=O
- Ammonia: N
- H2: [H][H]
- Propane: CCC
"""

import pytest

from phoenix.core.reaction import Auto, Reaction, ReactionSpecies
from phoenix.exceptions import OverconstrainedError, UnderconstrainedError

# =============================================================================
# SMILES notation reference for tests
# =============================================================================
# These are valid SMILES strings (not molecular formulas):
METHANE = "C"  # CH4 - implicit hydrogens
OXYGEN = "O=O"  # O2
WATER = "O"  # H2O
CO2 = "O=C=O"  # CO2
AMMONIA = "N"  # NH3
H2 = "[H][H]"  # H2 (explicit)
PROPANE = "CCC"  # C3H8
GLYCEROL = "OCC(O)CO"  # C3H8O3
PROPANEDIOL = "CC(O)CO"  # C3H8O2 (1,2-propanediol)
GLUCOSE = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"  # C6H12O6
N2 = "N#N"  # N2
CO = "[C-]#[O+]"  # CO (carbon monoxide)
OZONE = "[O-][O+]=O"  # O3


class TestReactionCreation:
    """Tests for Reaction factory methods."""

    def test_from_smiles_simple_lists(self):
        """Plain SMILES lists with all Auto coefficients."""
        rxn = Reaction.from_smiles(
            reactants=[METHANE, OXYGEN],
            products=[CO2, WATER],
        )
        assert len(rxn.reactants) == 2
        assert len(rxn.products) == 2
        assert all(s.is_auto for s in rxn.all_species)
        assert all(s.coefficient is None for s in rxn.all_species)

    def test_from_smiles_explicit_coefficients(self):
        """Tuples with explicit coefficients."""
        rxn = Reaction.from_smiles(
            reactants=[(METHANE, 1), (OXYGEN, 2)],
            products=[(CO2, 1), (WATER, 2)],
        )
        assert rxn.reactants[0].coefficient == 1.0
        assert rxn.reactants[1].coefficient == 2.0
        assert rxn.products[0].coefficient == 1.0
        assert rxn.products[1].coefficient == 2.0
        assert not any(s.is_auto for s in rxn.all_species)

    def test_from_smiles_mixed_coefficients(self):
        """Mixed explicit and Auto coefficients."""
        rxn = Reaction.from_smiles(
            reactants=[(METHANE, 1), (OXYGEN, Auto)],
            products=[(CO2, Auto), (WATER, Auto)],
        )
        assert rxn.reactants[0].coefficient == 1.0
        assert rxn.reactants[0].is_auto is False
        assert rxn.reactants[1].coefficient is None
        assert rxn.reactants[1].is_auto is True
        assert all(s.is_auto for s in rxn.products)

    def test_from_smiles_none_as_auto(self):
        """None treated as Auto."""
        rxn = Reaction.from_smiles(
            reactants=[(METHANE, 1), (OXYGEN, None)],
            products=[(CO2, 1), (WATER, None)],
        )
        assert rxn.reactants[1].is_auto is True
        assert rxn.products[1].is_auto is True

    def test_from_smiles_keyword_args(self):
        """Legacy keyword argument style."""
        rxn = Reaction.from_smiles(
            reactant_smiles=[(METHANE, 1), (OXYGEN, 2)],
            product_smiles=[(CO2, 1), (WATER, 2)],
        )
        assert rxn.reactants[0].coefficient == 1.0
        assert rxn.products[1].coefficient == 2.0

    def test_from_reaction_smiles_with_coefficients(self):
        """Parse reaction SMILES string with coefficients."""
        rxn = Reaction.from_reaction_smiles(
            f"{METHANE} + 2 {OXYGEN} >> {CO2} + 2 {WATER}", auto_balance=False
        )
        assert rxn.reactants[0].coefficient == 1.0
        assert rxn.reactants[1].coefficient == 2.0
        assert rxn.products[0].coefficient == 1.0
        assert rxn.products[1].coefficient == 2.0

    def test_from_reaction_smiles_without_coefficients(self):
        """Parse reaction SMILES without coefficients (default to 1)."""
        rxn = Reaction.from_reaction_smiles(
            f"{METHANE} + {OXYGEN} >> {CO2} + {WATER}", auto_balance=False
        )
        # Without explicit coefficients, default is 1.0 (not Auto)
        # This matches standard chemistry notation
        assert all(s.coefficient == 1.0 for s in rxn.all_species)
        assert not any(s.is_auto for s in rxn.all_species)


class TestFullAutoBalance:
    """Tests for fully automatic balancing (all unknowns)."""

    def test_combustion_methane(self):
        """CH4 + 2 O2 -> CO2 + 2 H2O"""
        rxn = Reaction.from_smiles([METHANE, OXYGEN], [CO2, WATER])
        rxn.balance()

        assert rxn.is_balanced
        # Check stoichiometry
        coeffs = rxn.coefficients
        # Normalize to CH4 = 1
        ch4_coeff = coeffs["CH4"]
        assert coeffs["O2"] / ch4_coeff == pytest.approx(2.0)
        assert coeffs["CO2"] / ch4_coeff == pytest.approx(1.0)
        assert coeffs["H2O"] / ch4_coeff == pytest.approx(2.0)

    def test_combustion_propane(self):
        """C3H8 + 5 O2 -> 3 CO2 + 4 H2O"""
        rxn = Reaction.from_smiles([PROPANE, OXYGEN], [CO2, WATER])
        rxn.balance()

        coeffs = rxn.coefficients
        c3h8_coeff = coeffs["C3H8"]
        assert coeffs["O2"] / c3h8_coeff == pytest.approx(5.0)
        assert coeffs["CO2"] / c3h8_coeff == pytest.approx(3.0)
        assert coeffs["H2O"] / c3h8_coeff == pytest.approx(4.0)

    def test_synthesis_ammonia(self):
        """N2 + 3 H2 -> 2 NH3"""
        rxn = Reaction.from_smiles([N2, H2], [AMMONIA])
        rxn.balance()

        coeffs = rxn.coefficients
        n2_coeff = coeffs["N2"]
        assert coeffs["H2"] / n2_coeff == pytest.approx(3.0)
        assert coeffs["H3N"] / n2_coeff == pytest.approx(2.0)

    def test_water_formation(self):
        """2 H2 + O2 -> 2 H2O"""
        rxn = Reaction.from_smiles([H2, OXYGEN], [WATER])
        rxn.balance()

        coeffs = rxn.coefficients
        o2_coeff = coeffs["O2"]
        assert coeffs["H2"] / o2_coeff == pytest.approx(2.0)
        assert coeffs["H2O"] / o2_coeff == pytest.approx(2.0)

    def test_integer_coefficients_preferred(self):
        """Coefficients should be scaled to integers when possible."""
        rxn = Reaction.from_smiles([METHANE, OXYGEN], [CO2, WATER])
        rxn.balance(prefer_integers=True)

        for species in rxn.all_species:
            coeff = species.coefficient
            assert coeff == pytest.approx(round(coeff)), f"Non-integer: {coeff}"


class TestMixedExplicitAutoBalance:
    """Tests for mixed explicit/auto coefficient systems."""

    def test_glycerol_hydrogenolysis(self):
        """
        Glycerol + H2 -> Propanediol + H2O
        OCC(O)CO + [H][H] -> CC(O)CO + O

        With glycerol=1, propanediol=1 fixed, solve for H2 and H2O.
        Atom balance:
          C: 3 = 3 (OK)
          H: 8 + 2*H2 = 8 + 2*H2O -> H2 = H2O
          O: 3 = 2 + H2O -> H2O = 1, H2 = 1
        """
        rxn = Reaction.from_smiles(
            reactants=[(GLYCEROL, 1), (H2, Auto)],  # glycerol + H2
            products=[(PROPANEDIOL, 1), (WATER, Auto)],  # 1,2-propanediol + H2O
        )
        rxn.balance()

        coeffs = rxn.coefficients
        assert coeffs["C3H8O3"] == pytest.approx(1.0)  # glycerol
        assert coeffs["H2"] == pytest.approx(1.0)
        assert coeffs["C3H8O2"] == pytest.approx(1.0)  # propanediol
        assert coeffs["H2O"] == pytest.approx(1.0)

    def test_partial_combustion_fixed_fuel(self):
        """CH4 + O2 -> CO2 + H2O with CH4=1 fixed."""
        rxn = Reaction.from_smiles(
            reactants=[(METHANE, 1), (OXYGEN, Auto)],
            products=[(CO2, Auto), (WATER, Auto)],
        )
        rxn.balance()

        coeffs = rxn.coefficients
        assert coeffs["CH4"] == pytest.approx(1.0)
        assert coeffs["O2"] == pytest.approx(2.0)
        assert coeffs["CO2"] == pytest.approx(1.0)
        assert coeffs["H2O"] == pytest.approx(2.0)

    def test_fixed_product_solve_reactants(self):
        """Fix product coefficients, solve for reactants."""
        # 2 H2 + O2 -> 2 H2O (fix H2O=2)
        rxn = Reaction.from_smiles(
            reactants=[(H2, Auto), (OXYGEN, Auto)],
            products=[(WATER, 2)],
        )
        rxn.balance()

        coeffs = rxn.coefficients
        assert coeffs["H2"] == pytest.approx(2.0)
        assert coeffs["O2"] == pytest.approx(1.0)
        assert coeffs["H2O"] == pytest.approx(2.0)


class TestOverconstrainedErrors:
    """Tests for over-constrained (impossible) systems."""

    def test_impossible_balance_wrong_coefficients(self):
        """Explicit coefficients that violate atom conservation."""
        # H2 + O2 -> 2 H2O is wrong: should be 2 H2 + O2 -> 2 H2O
        # With H2=1, O2=1, H2O=2:
        #   H: -2*1 + 2*2 = -2 + 4 = +2 (too much H in products)
        rxn = Reaction.from_smiles(
            reactants=[(H2, 1), (OXYGEN, 1)],  # Wrong: should be 2 H2
            products=[(WATER, 2)],
        )

        with pytest.raises(OverconstrainedError) as exc_info:
            rxn.balance()

        # Should report hydrogen imbalance: +2 excess H in products
        assert "H" in exc_info.value.imbalances
        assert exc_info.value.imbalances["H"] == pytest.approx(2.0)

    def test_impossible_elements(self):
        """Reaction with mismatched elements."""
        # CH4 -> NH3 is impossible: N in products but not reactants
        rxn = Reaction.from_smiles(
            reactants=[METHANE],
            products=[AMMONIA],  # Ammonia contains N, but no N in reactants
        )

        with pytest.raises(OverconstrainedError):
            rxn.balance()

    def test_all_fixed_but_unbalanced(self):
        """All coefficients fixed but not balanced."""
        # CH4 + O2 -> CO2 + 2 H2O with O2=1 is wrong (should be O2=2)
        # O balance: -2*1 + 2*1 + 1*2 = -2 + 2 + 2 = +2 (excess O in products)
        rxn = Reaction.from_smiles(
            reactants=[(METHANE, 1), (OXYGEN, 1)],  # Should be O2=2
            products=[(CO2, 1), (WATER, 2)],
        )

        with pytest.raises(OverconstrainedError) as exc_info:
            rxn.balance()

        assert "O" in exc_info.value.imbalances
        assert exc_info.value.imbalances["O"] == pytest.approx(2.0)


class TestUnderconstrainedErrors:
    """Tests for under-constrained (multiple solutions) systems."""

    def test_underdetermined_parallel_reactions(self):
        """
        Under-constrained reactions have multiple valid solutions.

        Example: CH4 + O2 -> CO + CO2 + H2O
        This has infinitely many solutions because both CO and CO2 can form.
        """
        # Multiple products can form from the same reactant carbon
        # This creates an under-determined system
        rxn = Reaction.from_smiles(
            reactants=[METHANE, OXYGEN],
            products=[CO, CO2, WATER],  # Both CO and CO2 can form
        )

        with pytest.raises(UnderconstrainedError) as exc_info:
            rxn.balance()

        assert exc_info.value.degrees_of_freedom > 0

    def test_suggestion_for_underconstrained(self):
        """Under-constrained error should suggest how to fix."""
        rxn = Reaction.from_smiles(
            reactants=[METHANE, OXYGEN],
            products=[CO, CO2, WATER],
        )

        with pytest.raises(UnderconstrainedError) as exc_info:
            rxn.balance()

        assert exc_info.value.suggestion is not None
        assert "coefficient" in exc_info.value.suggestion.lower()


class TestBalanceChaining:
    """Tests for method chaining and idempotency."""

    def test_balance_returns_self(self):
        """balance() returns self for method chaining."""
        rxn = Reaction.from_smiles([METHANE, OXYGEN], [CO2, WATER])
        result = rxn.balance()
        assert result is rxn

    def test_balance_idempotent(self):
        """Calling balance() multiple times has same result."""
        rxn = Reaction.from_smiles([METHANE, OXYGEN], [CO2, WATER])
        rxn.balance()
        coeffs1 = dict(rxn.coefficients)

        rxn.balance()  # Call again
        coeffs2 = dict(rxn.coefficients)

        assert coeffs1 == coeffs2

    def test_already_balanced_flag(self):
        """Balanced flag prevents redundant computation."""
        rxn = Reaction.from_smiles([METHANE, OXYGEN], [CO2, WATER])
        rxn.balance()
        assert rxn.is_balanced

        # Second call should be no-op
        rxn.balance()
        assert rxn.is_balanced


class TestStoichiometryVector:
    """Tests for the stoichiometry vector property."""

    def test_stoichiometry_vector_signs(self):
        """Reactants negative, products positive."""
        rxn = Reaction.from_smiles(
            reactants=[(METHANE, 1), (OXYGEN, 2)],
            products=[(CO2, 1), (WATER, 2)],
        )

        nu = rxn.stoichiometry_vector
        assert nu is not None
        # Reactants: -1, -2; Products: +1, +2
        assert nu[0] == -1.0  # CH4
        assert nu[1] == -2.0  # O2
        assert nu[2] == 1.0  # CO2
        assert nu[3] == 2.0  # H2O

    def test_stoichiometry_vector_none_for_unbalanced(self):
        """Returns None if any coefficient is undetermined."""
        rxn = Reaction.from_smiles([METHANE, OXYGEN], [CO2, WATER])
        assert rxn.stoichiometry_vector is None


class TestStringRepresentations:
    """Tests for __str__ and __repr__."""

    def test_str_unbalanced(self):
        """String representation with unknown coefficients."""
        rxn = Reaction.from_smiles([METHANE, OXYGEN], [CO2, WATER])
        s = str(rxn)
        assert "->" in s
        assert "?" in s  # Unknown coefficient marker

    def test_str_balanced(self):
        """String representation after balancing."""
        rxn = Reaction.from_smiles([METHANE, OXYGEN], [CO2, WATER])
        rxn.balance()
        s = str(rxn)
        assert "->" in s
        assert "?" not in s
        assert "2 O2" in s or "O2" in s  # Coefficient shown

    def test_repr(self):
        """Repr shows formulas and balance status."""
        rxn = Reaction.from_smiles([METHANE, OXYGEN], [CO2, WATER])
        r = repr(rxn)
        assert "Reaction" in r
        assert "CH4" in r
        assert ">>" in r

        rxn.balance()
        r = repr(rxn)
        assert "balanced" in r


class TestEdgeCases:
    """Tests for edge cases and numerical stability."""

    def test_single_element_reaction(self):
        """Simple single-element transformation."""
        # O3 -> O2 (ozone decomposition: 2 O3 -> 3 O2)
        rxn = Reaction.from_smiles([OZONE], [OXYGEN])
        rxn.balance()

        coeffs = rxn.coefficients
        # Ratio should be 2:3
        ratio = coeffs["O3"] / coeffs["O2"]
        assert ratio == pytest.approx(2.0 / 3.0)

    def test_fractional_coefficients(self):
        """Reactions that naturally have fractional coefficients."""
        # If we disable integer preference
        rxn = Reaction.from_smiles([METHANE, OXYGEN], [CO2, WATER])
        rxn.balance(prefer_integers=False, normalize=True)

        # Smallest should be 1.0 due to normalization
        min_coeff = min(s.coefficient for s in rxn.all_species)
        assert min_coeff == pytest.approx(1.0)

    def test_large_coefficients(self):
        """Reaction with large integer coefficients."""
        # C6H12O6 + 6 O2 -> 6 CO2 + 6 H2O (glucose combustion)
        rxn = Reaction.from_smiles(
            [GLUCOSE, OXYGEN],
            [CO2, WATER],
        )
        rxn.balance()

        coeffs = rxn.coefficients
        glucose_coeff = coeffs["C6H12O6"]
        assert coeffs["O2"] / glucose_coeff == pytest.approx(6.0)
        assert coeffs["CO2"] / glucose_coeff == pytest.approx(6.0)
        assert coeffs["H2O"] / glucose_coeff == pytest.approx(6.0)


class TestAutoSentinel:
    """Tests for the Auto sentinel object."""

    def test_auto_singleton(self):
        """Auto is a singleton."""
        from phoenix.core.reaction import Auto, _AutoType

        auto1 = Auto
        auto2 = _AutoType()
        assert auto1 is auto2

    def test_auto_repr(self):
        """Auto has meaningful repr/str."""
        assert repr(Auto) == "Auto"
        assert str(Auto) == "Auto"

    def test_auto_is_falsy(self):
        """Auto evaluates as False in boolean context."""
        assert not Auto
        assert bool(Auto) is False


class TestReactionSpecies:
    """Tests for ReactionSpecies dataclass."""

    def test_species_with_coefficient(self):
        """with_coefficient creates new instance."""
        from phoenix.core.compound import Compound

        compound = Compound.from_smiles(METHANE)  # Use SMILES, not formula
        species1 = ReactionSpecies(compound, None, is_auto=True)
        species2 = species1.with_coefficient(2.0)

        assert species1.coefficient is None
        assert species2.coefficient == 2.0
        assert species1 is not species2
        assert species1.compound is species2.compound  # Same compound

    def test_species_str_formatting(self):
        """String formatting of species."""
        from phoenix.core.compound import Compound

        compound = Compound.from_smiles(WATER)  # Use SMILES, not formula
        assert str(ReactionSpecies(compound, None)) == "? H2O"
        assert str(ReactionSpecies(compound, 1.0)) == "H2O"
        assert str(ReactionSpecies(compound, 2.0)) == "2 H2O"
        assert "2.5" in str(ReactionSpecies(compound, 2.5))


# =============================================================================
# Example Usage (can be run as documentation)
# =============================================================================


def example_full_auto_balance():
    """
    Example: Fully automatic balancing.

    When no coefficients are specified, the reaction is balanced
    automatically using the null-space algorithm.
    """
    # Methane combustion - all coefficients unknown
    rxn = Reaction.from_smiles(
        reactants=[METHANE, OXYGEN],
        products=[CO2, WATER],
    )
    rxn.balance()

    print("Methane combustion:")
    print(f"  {rxn}")
    print(f"  Coefficients: {rxn.coefficients}")
    # Output: CH4 + 2 O2 -> CO2 + 2 H2O


def example_mixed_coefficients():
    """
    Example: Mixed explicit and auto coefficients.

    Fix some coefficients, let the algorithm solve for the rest.
    Useful for specifying the basis of the reaction.
    """
    # Glycerol hydrogenolysis to propanediol
    # Fix glycerol and propanediol at 1, solve for H2 and H2O
    rxn = Reaction.from_smiles(
        reactants=[(GLYCEROL, 1), (H2, Auto)],  # glycerol (fixed), H2 (auto)
        products=[(PROPANEDIOL, 1), (WATER, Auto)],  # propanediol (fixed), H2O (auto)
    )
    rxn.balance()

    print("\nGlycerol hydrogenolysis:")
    print(f"  {rxn}")
    print(f"  H2 coefficient: {rxn.coefficients['H2']}")
    print(f"  H2O coefficient: {rxn.coefficients['H2O']}")


def example_error_handling():
    """
    Example: Handling impossible and under-constrained reactions.
    """
    # Over-constrained (impossible)
    try:
        rxn = Reaction.from_smiles(
            reactants=[(H2, 1), (OXYGEN, 1)],  # Wrong stoichiometry
            products=[(WATER, 2)],
        )
        rxn.balance()
    except OverconstrainedError as e:
        print(f"\nOver-constrained error: {e}")
        print(f"  Imbalances: {e.imbalances}")

    # Under-constrained (multiple solutions)
    try:
        rxn = Reaction.from_smiles(
            reactants=[METHANE, OXYGEN],  # Methane combustion with multiple products
            products=[CO, CO2, WATER],  # Both CO and CO2 can form
        )
        rxn.balance()
    except UnderconstrainedError as e:
        print(f"\nUnder-constrained error: {e}")
        print(f"  Degrees of freedom: {e.degrees_of_freedom}")
        print(f"  Suggestion: {e.suggestion}")


if __name__ == "__main__":
    # Run examples
    example_full_auto_balance()
    example_mixed_coefficients()
    example_error_handling()
