"""
Reaction class for chemical reactions with auto-balancing stoichiometry.

This module provides the core Reaction class for PHOENIX, supporting:
- Explicit stoichiometric coefficients (化学量論係数)
- Automatic coefficient detection via atom conservation (null-space algorithm)
- Mixed explicit/auto coefficients for partially constrained systems
- Comprehensive error handling for impossible balances

Theory
------
For a chemical reaction with n species and m elements, atom conservation gives:

    sum_i (nu_i * a_ij) = 0   for each element j

where nu_i is the stoichiometric coefficient (positive for products, negative
for reactants) and a_ij is the count of element j in species i.

This forms a homogeneous linear system: A * nu = 0, where A is the m x n
composition matrix. The solution lies in the null space of A.

For auto-balancing:
1. If all coefficients are unknown: find null space basis, choose simplest integer solution
2. If some are known: solve the constrained system A_unknown * x = -A_known * x_known

Examples
--------
>>> # Full auto-balance (all coefficients determined)
>>> rxn = Reaction.from_smiles(
...     reactants=["CH4", "O2"],
...     products=["CO2", "H2O"]
... )
>>> rxn.balance()
>>> print(rxn)  # CH4 + 2 O2 -> CO2 + 2 H2O

>>> # Mixed explicit/auto coefficients
>>> from phoenix import Auto
>>> rxn = Reaction.from_smiles(
...     reactants=[("OCC(O)CO", 1), ("[H][H]", Auto)],  # glycerol + H2
...     products=[("CC(O)CO", 1), ("O", Auto)]          # propanediol + H2O
... )
>>> rxn.balance()
>>> print(rxn.coefficients)  # {glycerol: 1, H2: 1, propanediol: 1, H2O: 1}
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from math import gcd
from typing import TYPE_CHECKING, overload

import numpy as np
from scipy import linalg as scipy_linalg

from phoenix.core.compound import Compound
from phoenix.exceptions import (
    BalanceError,
    OverconstrainedError,
    UnderconstrainedError,
)

if TYPE_CHECKING:
    from collections.abc import Sequence


# =============================================================================
# Sentinel for Auto-Balancing
# =============================================================================


class _AutoType:
    """
    Sentinel value indicating a coefficient should be auto-calculated.

    Use `Auto` (the singleton instance) to mark coefficients for auto-balancing:

        rxn = Reaction.from_smiles(
            reactants=[("CH4", 1), ("O2", Auto)],
            products=[("CO2", Auto), ("H2O", Auto)]
        )

    The balance() method will solve for coefficients marked with Auto.
    """

    _instance: _AutoType | None = None

    def __new__(cls) -> _AutoType:
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __repr__(self) -> str:
        return "Auto"

    def __str__(self) -> str:
        return "Auto"

    def __bool__(self) -> bool:
        # Auto is falsy so `if coeff:` returns False for Auto
        return False


# Singleton instance
Auto = _AutoType()

# Type alias for coefficient specification
CoeffSpec = float | int | _AutoType | None  # None is treated as Auto


# =============================================================================
# Data Classes for Reaction Components
# =============================================================================


@dataclass(frozen=True, slots=True)
class ReactionSpecies:
    """
    A species (compound) in a reaction with its stoichiometric coefficient.

    Attributes
    ----------
    compound : Compound
        The chemical compound
    coefficient : float | None
        Stoichiometric coefficient (None if not yet determined)
    is_auto : bool
        Whether coefficient was marked for auto-calculation
    """

    compound: Compound
    coefficient: float | None
    is_auto: bool = False

    @property
    def formula(self) -> str:
        """Molecular formula of the compound."""
        return self.compound.formula

    @property
    def composition(self) -> dict[str, int]:
        """Elemental composition of the compound."""
        return self.compound.composition

    def with_coefficient(self, coeff: float) -> ReactionSpecies:
        """Return a new ReactionSpecies with the specified coefficient."""
        return ReactionSpecies(self.compound, coeff, is_auto=self.is_auto)

    def __str__(self) -> str:
        if self.coefficient is None:
            return f"? {self.formula}"
        # Use approximate comparison for floating point
        coeff = float(self.coefficient)
        if abs(coeff - 1.0) < 1e-9:
            return self.formula
        elif abs(coeff - round(coeff)) < 1e-9:
            return f"{int(round(coeff))} {self.formula}"
        else:
            return f"{coeff:.4g} {self.formula}"


# =============================================================================
# Main Reaction Class
# =============================================================================


class Reaction:
    """
    Chemical reaction with stoichiometric coefficients and thermodynamic evaluation.

    The Reaction class supports:
    - Explicit coefficients for all species
    - Automatic coefficient calculation using atom conservation
    - Mixed explicit/auto coefficients for partially constrained systems

    Coefficient Types
    -----------------
    - Numeric (int, float): Explicit coefficient value
    - `Auto` or `None`: Coefficient to be auto-calculated during balance()

    Parameters
    ----------
    reactants : list[ReactionSpecies]
        List of reactant species with coefficients
    products : list[ReactionSpecies]
        List of product species with coefficients

    Examples
    --------
    >>> # All explicit coefficients
    >>> rxn = Reaction.from_smiles(
    ...     reactants=[("CH4", 1), ("O2", 2)],
    ...     products=[("CO2", 1), ("H2O", 2)]
    ... )

    >>> # Full auto-balance
    >>> rxn = Reaction.from_smiles(
    ...     reactants=["CH4", "O2"],  # No coefficients -> all Auto
    ...     products=["CO2", "H2O"]
    ... )
    >>> rxn.balance()

    >>> # Mixed explicit/auto
    >>> rxn = Reaction.from_smiles(
    ...     reactants=[("OCC(O)CO", 1), ("[H][H]", Auto)],
    ...     products=[("CC(O)CO", 1), ("O", Auto)]
    ... )
    >>> rxn.balance()
    """

    def __init__(
        self,
        reactants: list[ReactionSpecies],
        products: list[ReactionSpecies],
    ):
        self._reactants = list(reactants)
        self._products = list(products)
        self._balanced = False
        self._balance_info: dict | None = None

    # =========================================================================
    # Factory Methods
    # =========================================================================

    @overload
    @classmethod
    def from_smiles(
        cls,
        reactants: Sequence[str],
        products: Sequence[str],
    ) -> Reaction: ...

    @overload
    @classmethod
    def from_smiles(
        cls,
        reactants: Sequence[tuple[str, CoeffSpec]],
        products: Sequence[tuple[str, CoeffSpec]],
    ) -> Reaction: ...

    @overload
    @classmethod
    def from_smiles(
        cls,
        *,
        reactant_smiles: Sequence[str] | Sequence[tuple[str, CoeffSpec]],
        product_smiles: Sequence[str] | Sequence[tuple[str, CoeffSpec]],
    ) -> Reaction: ...

    @classmethod
    def from_smiles(
        cls,
        reactants: Sequence[str] | Sequence[tuple[str, CoeffSpec]] | None = None,
        products: Sequence[str] | Sequence[tuple[str, CoeffSpec]] | None = None,
        *,
        reactant_smiles: Sequence[str] | Sequence[tuple[str, CoeffSpec]] | None = None,
        product_smiles: Sequence[str] | Sequence[tuple[str, CoeffSpec]] | None = None,
    ) -> Reaction:
        """
        Create a Reaction from SMILES strings with optional coefficients.

        Supports multiple input formats for flexibility:

        1. Simple SMILES lists (all coefficients Auto):
           >>> Reaction.from_smiles(["CH4", "O2"], ["CO2", "H2O"])

        2. Tuples with explicit coefficients:
           >>> Reaction.from_smiles(
           ...     reactants=[("CH4", 1), ("O2", 2)],
           ...     products=[("CO2", 1), ("H2O", 2)]
           ... )

        3. Mixed explicit/Auto coefficients:
           >>> Reaction.from_smiles(
           ...     reactants=[("OCC(O)CO", 1), ("[H][H]", Auto)],
           ...     products=[("CC(O)CO", 1), ("O", Auto)]
           ... )

        4. Legacy keyword argument style:
           >>> Reaction.from_smiles(
           ...     reactant_smiles=[("CH4", 1), ("O2", 2)],
           ...     product_smiles=[("CO2", 1), ("H2O", 2)]
           ... )

        Parameters
        ----------
        reactants : Sequence[str] | Sequence[tuple[str, CoeffSpec]]
            Reactant SMILES, optionally with coefficients
        products : Sequence[str] | Sequence[tuple[str, CoeffSpec]]
            Product SMILES, optionally with coefficients
        reactant_smiles : Sequence (keyword, legacy)
            Alternative way to specify reactants
        product_smiles : Sequence (keyword, legacy)
            Alternative way to specify products

        Returns
        -------
        Reaction
            New Reaction instance

        Raises
        ------
        ValueError
            If neither positional nor keyword arguments provided
        InvalidSmilesError
            If any SMILES string is invalid
        """
        # Handle both positional and keyword arguments for backwards compatibility
        r_input = reactants if reactants is not None else reactant_smiles
        p_input = products if products is not None else product_smiles

        if r_input is None or p_input is None:
            raise ValueError(
                "Must provide both reactants and products. "
                "Use positional args or reactant_smiles/product_smiles keywords."
            )

        def parse_species_list(
            items: Sequence[str] | Sequence[tuple[str, CoeffSpec]],
        ) -> list[ReactionSpecies]:
            """Parse a list of SMILES or (SMILES, coeff) tuples."""
            species_list = []
            for item in items:
                if isinstance(item, str):
                    # Plain SMILES string -> Auto coefficient
                    compound = Compound.from_smiles(item)
                    species_list.append(ReactionSpecies(compound, None, is_auto=True))
                elif isinstance(item, tuple) and len(item) == 2:
                    smiles, coeff = item
                    compound = Compound.from_smiles(smiles)

                    if coeff is None or coeff is Auto:
                        species_list.append(ReactionSpecies(compound, None, is_auto=True))
                    elif isinstance(coeff, (int, float)):
                        species_list.append(
                            ReactionSpecies(compound, float(coeff), is_auto=False)
                        )
                    else:
                        raise ValueError(
                            f"Invalid coefficient type: {type(coeff)}. "
                            f"Use int, float, Auto, or None."
                        )
                else:
                    raise ValueError(
                        f"Invalid species format: {item}. "
                        f"Use SMILES string or (SMILES, coefficient) tuple."
                    )
            return species_list

        reactant_species = parse_species_list(r_input)
        product_species = parse_species_list(p_input)

        return cls(reactant_species, product_species)

    @classmethod
    def from_reaction_smiles(cls, reaction_smiles: str, auto_balance: bool = True) -> Reaction:
        """
        Create a Reaction from a reaction SMILES string.

        Format: 'coeff SMILES + coeff SMILES >> coeff SMILES + coeff SMILES'
        Coefficients are optional; if omitted, marked as Auto.

        Parameters
        ----------
        reaction_smiles : str
            Reaction SMILES with '>>' separator
        auto_balance : bool
            If True, automatically balance after parsing (default: True)

        Returns
        -------
        Reaction
            Parsed reaction

        Examples
        --------
        >>> rxn = Reaction.from_reaction_smiles("CH4 + 2 O2 >> CO2 + 2 H2O")
        >>> rxn = Reaction.from_reaction_smiles("OCC(O)CO + [H][H] >> CC(O)CO + O")
        """
        if ">>" not in reaction_smiles:
            raise ValueError("Reaction SMILES must contain '>>' separator")

        left, right = reaction_smiles.split(">>", 1)

        def parse_side(side: str) -> list[tuple[str, CoeffSpec]]:
            """Parse one side of the reaction."""
            species = []
            parts = re.split(r"\s*\+\s*", side.strip())
            for part in parts:
                part = part.strip()
                if not part:
                    continue

                # Match optional coefficient followed by SMILES
                # Supports: "2 H2O", "1.5 CH4", "H2O", "2H2O" (no space)
                match = re.match(r"^(\d+\.?\d*)?\s*(.+)$", part)
                if not match:
                    raise ValueError(f"Invalid compound format: {part}")

                coeff_str, smiles = match.groups()
                # If no coefficient specified, default to 1.0 (not Auto)
                # This matches standard chemistry notation
                if coeff_str:
                    coeff: CoeffSpec = float(coeff_str)
                else:
                    coeff = 1.0
                species.append((smiles, coeff))
            return species

        reactant_specs = parse_side(left)
        product_specs = parse_side(right)

        rxn = cls.from_smiles(reactant_specs, product_specs)

        if auto_balance:
            rxn.balance()

        return rxn

    # =========================================================================
    # Properties
    # =========================================================================

    @property
    def reactants(self) -> list[ReactionSpecies]:
        """List of reactant species."""
        return list(self._reactants)

    @property
    def products(self) -> list[ReactionSpecies]:
        """List of product species."""
        return list(self._products)

    @property
    def all_species(self) -> list[ReactionSpecies]:
        """All species in order: reactants then products."""
        return self._reactants + self._products

    @property
    def is_balanced(self) -> bool:
        """Whether the reaction has been balanced."""
        return self._balanced

    @property
    def elements(self) -> set[str]:
        """Set of all elements in the reaction."""
        elements: set[str] = set()
        for species in self.all_species:
            elements.update(species.composition.keys())
        return elements

    @property
    def coefficients(self) -> dict[str, float | None]:
        """
        Dictionary mapping formula to coefficient.

        Returns coefficients in order: reactants (positive), products (positive).
        Note: For thermodynamic calculations, reactant coefficients are negated.
        """
        return {s.formula: s.coefficient for s in self.all_species}

    @property
    def stoichiometry_vector(self) -> np.ndarray | None:
        """
        Stoichiometry vector nu where nu_i < 0 for reactants, > 0 for products.

        Returns None if any coefficient is undetermined.
        """
        coeffs = []
        for species in self._reactants:
            if species.coefficient is None:
                return None
            coeffs.append(-species.coefficient)  # Negative for reactants
        for species in self._products:
            if species.coefficient is None:
                return None
            coeffs.append(species.coefficient)  # Positive for products
        return np.array(coeffs)

    # =========================================================================
    # Balancing Algorithm
    # =========================================================================

    def balance(self, *, normalize: bool = True, prefer_integers: bool = True) -> Reaction:
        """
        Balance the reaction by solving for unknown coefficients.

        Uses the null-space algorithm for atom conservation:
        - For fully unconstrained: finds basis of null(A) and selects simplest solution
        - For partially constrained: solves the linear system with known coefficients

        Parameters
        ----------
        normalize : bool
            If True, normalize so smallest coefficient is 1 (default: True)
        prefer_integers : bool
            If True, scale to get integer coefficients when possible (default: True)

        Returns
        -------
        Reaction
            Self (for method chaining)

        Raises
        ------
        OverconstrainedError
            If the constraints are inconsistent (no solution exists)
        UnderconstrainedError
            If multiple solutions exist (need more constraints)

        Notes
        -----
        The balancing problem is formulated as:

        For m elements and n species with composition matrix A (m x n):
            A @ nu = 0  (atom conservation)

        where nu is the stoichiometry vector (negative for reactants).

        If some coefficients are known, we partition: A @ [nu_known; nu_unknown] = 0
        which gives: A_unknown @ nu_unknown = -A_known @ nu_known

        Examples
        --------
        >>> rxn = Reaction.from_smiles(["CH4", "O2"], ["CO2", "H2O"])
        >>> rxn.balance()
        >>> print(rxn)
        CH4 + 2 O2 -> CO2 + 2 H2O

        >>> rxn = Reaction.from_smiles(
        ...     reactants=[("OCC(O)CO", 1), ("[H][H]", Auto)],
        ...     products=[("CC(O)CO", 1), ("O", Auto)]
        ... )
        >>> rxn.balance()
        >>> print(rxn.coefficients)
        {'C3H8O3': 1.0, 'H2': 1.0, 'C3H8O2': 1.0, 'H2O': 1.0}
        """
        if self._balanced:
            return self

        # Collect indices of unknown coefficients
        n_reactants = len(self._reactants)
        n_species = len(self.all_species)
        species_list = self.all_species

        unknown_indices = []
        known_indices = []
        known_coeffs = []

        for i, species in enumerate(species_list):
            if species.is_auto or species.coefficient is None:
                unknown_indices.append(i)
            else:
                known_indices.append(i)
                known_coeffs.append(species.coefficient)

        # Get sorted elements for consistent ordering
        elements = sorted(self.elements)
        n_elements = len(elements)
        n_unknowns = len(unknown_indices)

        # Build the composition matrix A where A[j, i] = count of element j in species i
        # Sign convention: reactants negative, products positive
        def get_signed_composition(species_idx: int, elem: str) -> float:
            species = species_list[species_idx]
            count = species.composition.get(elem, 0)
            if species_idx < n_reactants:
                return -count  # Reactant (consumed)
            else:
                return count  # Product (produced)

        if not unknown_indices:
            # All coefficients known - just verify balance
            self._verify_balance(elements)
            self._balanced = True
            return self

        if not known_indices:
            # All coefficients unknown - use null space method
            self._balance_full_null_space(
                elements, species_list, n_reactants, normalize, prefer_integers
            )
        else:
            # Mixed known/unknown - solve constrained system
            self._balance_mixed_system(
                elements,
                species_list,
                n_reactants,
                unknown_indices,
                known_indices,
                known_coeffs,
                normalize,
                prefer_integers,
            )

        self._balanced = True
        return self

    def _balance_full_null_space(
        self,
        elements: list[str],
        species_list: list[ReactionSpecies],
        n_reactants: int,
        normalize: bool,
        prefer_integers: bool,
    ) -> None:
        """
        Balance when all coefficients are unknown using null space.

        The stoichiometry vector must lie in null(A) where A is the composition matrix.
        """
        n_elements = len(elements)
        n_species = len(species_list)

        # Build composition matrix
        A = np.zeros((n_elements, n_species))
        for j, elem in enumerate(elements):
            for i, species in enumerate(species_list):
                count = species.composition.get(elem, 0)
                # Reactants negative, products positive
                sign = -1 if i < n_reactants else 1
                A[j, i] = sign * count

        # Compute null space using SVD
        # A @ nu = 0 => nu is in null(A)
        null_space = self._compute_null_space(A)
        null_dim = null_space.shape[1]

        if null_dim == 0:
            # No solution exists
            raise OverconstrainedError(
                "No valid stoichiometry exists for this reaction. "
                "Check that reactants can actually form products."
            )

        if null_dim > 1:
            # Multiple independent solutions
            raise UnderconstrainedError(
                degrees_of_freedom=null_dim,
                suggestion=(
                    f"Fix at least {null_dim - 1} coefficient(s) to get a unique solution. "
                    f"Consider specifying coefficients for some species."
                ),
            )

        # Single solution (null_dim == 1)
        nu = null_space[:, 0]

        # The null space vector satisfies A @ nu = 0, where:
        # A[j,i] = sign_i * count_ij, with sign = -1 for reactants, +1 for products
        #
        # For CH4 + 2 O2 -> CO2 + 2 H2O:
        # C balance: -1*nu[CH4] + 1*nu[CO2] = 0
        # H balance: -4*nu[CH4] + 2*nu[H2O] = 0
        # O balance: -2*nu[O2] + 2*nu[CO2] + 1*nu[H2O] = 0
        #
        # Solution: nu = [1, 2, 1, 2] (all positive!)
        # These ARE the stoichiometric coefficients directly.
        #
        # The null space solution might be scaled or flipped. We want all positive.
        # Check if any element is negative - if so, flip the entire vector.

        # Count positive and negative entries
        n_positive = np.sum(nu > 1e-10)
        n_negative = np.sum(nu < -1e-10)

        # If more negative than positive, flip
        if n_negative > n_positive:
            nu = -nu

        # The coefficients are directly the absolute values of nu
        # Since all should be positive after proper orientation
        coefficients = list(np.abs(nu))

        # Process coefficients (normalize, integerize)
        coefficients = self._process_coefficients(
            coefficients, normalize, prefer_integers
        )

        # Validate: all coefficients should be positive
        for i, coeff in enumerate(coefficients):
            if coeff <= 1e-10:
                species = species_list[i]
                side = "reactant" if i < n_reactants else "product"
                raise OverconstrainedError(
                    f"Invalid solution: {side} {species.formula} has non-positive "
                    f"coefficient {coeff:.4f}. The reaction may be written incorrectly."
                )

        # Update species with solved coefficients
        self._update_coefficients(coefficients)

    def _balance_mixed_system(
        self,
        elements: list[str],
        species_list: list[ReactionSpecies],
        n_reactants: int,
        unknown_indices: list[int],
        known_indices: list[int],
        known_coeffs: list[float],
        normalize: bool,
        prefer_integers: bool,
    ) -> None:
        """
        Balance with some coefficients known using constrained linear solve.

        The atom conservation equation is: A @ coeff = 0
        where A[j,i] = -count[j,i] for reactants, +count[j,i] for products,
        and coeff[i] is the (always positive) stoichiometric coefficient.

        Partitioning into known and unknown:
            A_k @ coeff_k + A_u @ coeff_u = 0
            A_u @ coeff_u = -A_k @ coeff_k

        We solve for coeff_u directly (no sign flipping needed).
        """
        n_elements = len(elements)
        n_unknowns = len(unknown_indices)
        n_known = len(known_indices)

        # Build composition matrices for known and unknown species
        # A[j,i] = sign * count, where sign = -1 for reactants, +1 for products
        A_unknown = np.zeros((n_elements, n_unknowns))
        A_known = np.zeros((n_elements, n_known))

        for j, elem in enumerate(elements):
            for ui, species_idx in enumerate(unknown_indices):
                species = species_list[species_idx]
                count = species.composition.get(elem, 0)
                sign = -1 if species_idx < n_reactants else 1
                A_unknown[j, ui] = sign * count

            for ki, species_idx in enumerate(known_indices):
                species = species_list[species_idx]
                count = species.composition.get(elem, 0)
                sign = -1 if species_idx < n_reactants else 1
                A_known[j, ki] = sign * count

        # Known coefficients (always positive, as given by user)
        coeff_known = np.array(known_coeffs)

        # Right-hand side: b = -A_known @ coeff_known
        b = -A_known @ coeff_known

        # Solve A_unknown @ coeff_unknown = b using least squares
        coeff_unknown, residuals, rank, s = np.linalg.lstsq(A_unknown, b, rcond=None)

        # Check solution quality
        residual_norm = np.linalg.norm(A_unknown @ coeff_unknown - b)

        if residual_norm > 1e-8:
            # System is inconsistent - compute imbalances
            imbalances = {}
            # Build full coefficient vector
            full_coeff = np.zeros(len(species_list))
            for i, species_idx in enumerate(known_indices):
                full_coeff[species_idx] = known_coeffs[i]
            for i, species_idx in enumerate(unknown_indices):
                full_coeff[species_idx] = coeff_unknown[i]

            A_full = self._build_composition_matrix(elements, species_list, n_reactants)
            imbalance_vec = A_full @ full_coeff
            for j, elem in enumerate(elements):
                if abs(imbalance_vec[j]) > 1e-8:
                    imbalances[elem] = float(imbalance_vec[j])

            raise OverconstrainedError(imbalances=imbalances)

        # Check for under-constrained system
        null_space = self._compute_null_space(A_unknown)
        if null_space.shape[1] > 0:
            raise UnderconstrainedError(
                degrees_of_freedom=null_space.shape[1],
                suggestion=(
                    f"Specify {null_space.shape[1]} more coefficient(s) to get a unique solution."
                ),
            )

        # The solved coefficients should be positive directly
        unknown_coeffs_list = list(coeff_unknown)

        # Validate: all unknowns should be positive
        for i, coeff in enumerate(unknown_coeffs_list):
            if coeff <= 1e-10:
                species_idx = unknown_indices[i]
                species = species_list[species_idx]
                side = "reactant" if species_idx < n_reactants else "product"
                raise OverconstrainedError(
                    f"Invalid solution: {side} {species.formula} has non-positive "
                    f"coefficient {coeff:.6f}. Check constraint consistency."
                )

        # Reconstruct full coefficient list
        full_coeffs = [0.0] * len(species_list)
        for i, species_idx in enumerate(known_indices):
            full_coeffs[species_idx] = known_coeffs[i]
        for i, species_idx in enumerate(unknown_indices):
            full_coeffs[species_idx] = unknown_coeffs_list[i]

        # Process coefficients
        if normalize or prefer_integers:
            # Only process the auto-calculated ones
            min_unknown = min(unknown_coeffs_list)
            if normalize and min_unknown > 0:
                scale = 1.0 / min_unknown
                for i, species_idx in enumerate(unknown_indices):
                    full_coeffs[species_idx] *= scale

            if prefer_integers:
                # Try to find integer scaling for unknowns only
                unknown_vals = [full_coeffs[idx] for idx in unknown_indices]
                int_scale = self._find_integer_scale(unknown_vals)
                if int_scale is not None:
                    for idx in unknown_indices:
                        full_coeffs[idx] *= int_scale

        # Update species
        self._update_coefficients(full_coeffs)

    def _compute_null_space(self, A: np.ndarray, tol: float = 1e-10) -> np.ndarray:
        """
        Compute null space of matrix A using SVD.

        Returns matrix where columns span null(A).
        """
        U, S, Vh = scipy_linalg.svd(A, full_matrices=True)
        n = A.shape[1]
        rank = np.sum(S > tol)
        null_space = Vh[rank:, :].T
        return null_space

    def _build_composition_matrix(
        self, elements: list[str], species_list: list[ReactionSpecies], n_reactants: int
    ) -> np.ndarray:
        """Build the full composition matrix with sign convention."""
        n_elements = len(elements)
        n_species = len(species_list)
        A = np.zeros((n_elements, n_species))
        for j, elem in enumerate(elements):
            for i, species in enumerate(species_list):
                count = species.composition.get(elem, 0)
                sign = -1 if i < n_reactants else 1
                A[j, i] = sign * count
        return A

    def _reconstruct_nu(
        self,
        species_list: list[ReactionSpecies],
        n_reactants: int,
        unknown_indices: list[int],
        known_indices: list[int],
        nu_unknown: np.ndarray,
        known_coeffs: list[float],
    ) -> np.ndarray:
        """Reconstruct full stoichiometry vector from known and unknown parts."""
        n_species = len(species_list)
        nu = np.zeros(n_species)

        for i, species_idx in enumerate(known_indices):
            coeff = known_coeffs[i]
            sign = -1 if species_idx < n_reactants else 1
            nu[species_idx] = sign * coeff

        for i, species_idx in enumerate(unknown_indices):
            nu[species_idx] = nu_unknown[i]

        return nu

    def _process_coefficients(
        self, coefficients: list[float], normalize: bool, prefer_integers: bool
    ) -> list[float]:
        """Normalize and/or convert coefficients to integers."""
        coeffs = list(coefficients)

        if normalize:
            min_coeff = min(c for c in coeffs if c > 0)
            coeffs = [c / min_coeff for c in coeffs]

        if prefer_integers:
            int_scale = self._find_integer_scale(coeffs)
            if int_scale is not None:
                coeffs = [c * int_scale for c in coeffs]

        return coeffs

    def _find_integer_scale(
        self, coefficients: list[float], max_denominator: int = 100
    ) -> float | None:
        """
        Find scaling factor to convert coefficients to integers.

        Uses continued fraction approximation to find rational representations.
        """
        # Try to express each coefficient as a fraction
        fractions = []
        for c in coefficients:
            if c <= 0:
                continue
            try:
                frac = Fraction(c).limit_denominator(max_denominator)
                # Check if approximation is close enough
                if abs(float(frac) - c) > 1e-6:
                    return None
                fractions.append(frac)
            except (ValueError, ZeroDivisionError):
                return None

        if not fractions:
            return None

        # LCM of denominators
        def lcm(a: int, b: int) -> int:
            return abs(a * b) // gcd(a, b)

        denominators = [f.denominator for f in fractions]
        common_denom = reduce(lcm, denominators, 1)

        # Check if scaled values are close to integers
        for c in coefficients:
            if c <= 0:
                continue
            scaled = c * common_denom
            if abs(scaled - round(scaled)) > 1e-6:
                return None

        return float(common_denom)

    def _update_coefficients(self, coefficients: list[float]) -> None:
        """Update internal species lists with new coefficients."""
        n_reactants = len(self._reactants)

        new_reactants = []
        for i, species in enumerate(self._reactants):
            new_reactants.append(species.with_coefficient(coefficients[i]))
        self._reactants = new_reactants

        new_products = []
        for i, species in enumerate(self._products):
            new_products.append(species.with_coefficient(coefficients[n_reactants + i]))
        self._products = new_products

    def _verify_balance(self, elements: list[str]) -> None:
        """Verify that explicit coefficients satisfy atom conservation."""
        imbalances = {}

        for elem in elements:
            balance = 0.0
            for species in self._reactants:
                if species.coefficient is None:
                    raise ValueError(f"Cannot verify: {species.formula} has no coefficient")
                balance -= species.coefficient * species.composition.get(elem, 0)
            for species in self._products:
                if species.coefficient is None:
                    raise ValueError(f"Cannot verify: {species.formula} has no coefficient")
                balance += species.coefficient * species.composition.get(elem, 0)

            if abs(balance) > 1e-8:
                imbalances[elem] = balance

        if imbalances:
            raise OverconstrainedError(imbalances=imbalances)

    # =========================================================================
    # Thermodynamic Properties
    # =========================================================================

    @property
    def enthalpy(self) -> float:
        """
        Standard enthalpy of reaction (Delta H_r) in kJ/mol.

        Auto-balances if needed.
        """
        if not self._balanced:
            self.balance()

        h_products = sum(
            s.coefficient * float(s.compound.enthalpy_of_formation)
            for s in self._products
            if s.coefficient is not None
        )
        h_reactants = sum(
            s.coefficient * float(s.compound.enthalpy_of_formation)
            for s in self._reactants
            if s.coefficient is not None
        )
        return h_products - h_reactants

    @property
    def delta_h(self) -> float:
        """Alias for enthalpy."""
        return self.enthalpy

    @property
    def entropy(self) -> float:
        """
        Standard entropy of reaction (Delta S_r) in J/(mol*K).

        Auto-balances if needed.
        """
        if not self._balanced:
            self.balance()

        s_products = sum(
            s.coefficient * float(s.compound.entropy)
            for s in self._products
            if s.coefficient is not None
        )
        s_reactants = sum(
            s.coefficient * float(s.compound.entropy)
            for s in self._reactants
            if s.coefficient is not None
        )
        return s_products - s_reactants

    @property
    def delta_s(self) -> float:
        """Alias for entropy."""
        return self.entropy

    @property
    def gibbs_free_energy(self) -> float:
        """
        Standard Gibbs free energy of reaction (Delta G_r) in kJ/mol at 298.15 K.

        Calculated as: Delta G = Delta H - T * Delta S
        """
        if not self._balanced:
            self.balance()
        T = 298.15  # K
        # Convert entropy from J/(mol*K) to kJ/(mol*K)
        return self.enthalpy - T * (self.entropy / 1000.0)

    @property
    def delta_g(self) -> float:
        """Alias for gibbs_free_energy."""
        return self.gibbs_free_energy

    # =========================================================================
    # String Representations
    # =========================================================================

    def __str__(self) -> str:
        """Human-readable reaction equation."""
        r_strs = []
        for species in self._reactants:
            r_strs.append(str(species))
        p_strs = []
        for species in self._products:
            p_strs.append(str(species))

        return f"{' + '.join(r_strs)} -> {' + '.join(p_strs)}"

    def __repr__(self) -> str:
        r_formulas = [s.formula for s in self._reactants]
        p_formulas = [s.formula for s in self._products]
        balanced_str = ", balanced" if self._balanced else ""
        return f"Reaction({' + '.join(r_formulas)} >> {' + '.join(p_formulas)}{balanced_str})"

    def to_equation(self, *, use_names: bool = False) -> str:
        """
        Format as a chemical equation string.

        Parameters
        ----------
        use_names : bool
            If True, use compound names instead of formulas (not yet implemented)

        Returns
        -------
        str
            Formatted equation like "CH4 + 2 O2 -> CO2 + 2 H2O"
        """
        return str(self)
