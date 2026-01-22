"""
Compound class for molecular representation and property access.

The Compound class is the primary entry point for PHOENIX, providing
SMILES parsing, validation, and lazy-evaluated property access.
"""

from __future__ import annotations

import warnings
from functools import cached_property
from typing import TYPE_CHECKING

from rdkit import Chem
from rdkit.Chem import Descriptors, inchi, rdMolDescriptors

from phoenix.exceptions import (
    InvalidSmilesError,
    UnsupportedElementError,
    UnsupportedStructureError,
)

from phoenix.thermo.models import ThermoProperty, ThermoPropertyAccessor, ThermoState

if TYPE_CHECKING:
    from phoenix.hazard.classification import HazardResult

# Supported elements for MVP (C, H, N, O, S, P, F, Cl, Br)
SUPPORTED_ELEMENTS = frozenset({"C", "H", "N", "O", "S", "P", "F", "Cl", "Br"})

# Large molecule warning threshold
LARGE_MOLECULE_THRESHOLD = 100


class Compound:
    """
    A chemical compound parsed from SMILES with lazy property evaluation.

    The Compound class provides:
    - SMILES parsing and canonicalization via RDKit
    - Input validation (element support, charge, radical checks)
    - Elemental composition extraction
    - Molecular weight and formula calculation
    - Thermodynamic property estimation (via phoenix.thermo)
    - Hazard evaluation (via phoenix.hazard)

    Parameters
    ----------
    rdmol : rdkit.Chem.Mol
        RDKit molecule object (use from_smiles() for construction)

    Examples
    --------
    >>> compound = Compound.from_smiles("CCO")
    >>> compound.formula
    'C2H6O'
    >>> compound.molecular_weight
    46.07
    >>> compound.composition
    {'C': 2, 'H': 6, 'O': 1}
    """

    def __init__(self, rdmol: Chem.Mol, original_smiles: str | None = None, phase: str = "g"):
        """
        Initialize Compound from RDKit Mol object.

        Use Compound.from_smiles() for the standard entry point.
        """
        self._rdmol = rdmol
        self._smiles = original_smiles
        self._phase = phase
        self._warnings: list[str] = []

    @classmethod
    def from_smiles(cls, smiles: str, phase: str = "g") -> Compound:
        """
        Create a Compound from a SMILES string.

        Parameters
        ----------
        smiles : str
            SMILES string representing the molecule
        phase : str
            Phase of the compound ('g' for gas, 'l' for liquid, 's' for solid).
            Default is 'g'.

        Returns
        -------
        Compound
            Validated compound object
        """
        # Parse SMILES
        rdmol = Chem.MolFromSmiles(smiles)
        if rdmol is None:
            raise InvalidSmilesError(smiles)

        # Add hydrogens for proper atom counting
        rdmol = Chem.AddHs(rdmol)

        # Validate the molecule
        compound = cls(rdmol, original_smiles=smiles, phase=phase)
        compound._validate()

        return compound

    def _validate(self) -> None:
        """Validate the molecule for PHOENIX compatibility."""
        # Check for unsupported elements
        unsupported = []
        for atom in self._rdmol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol not in SUPPORTED_ELEMENTS:
                unsupported.append(symbol)

        if unsupported:
            # Check specifically for iodine
            if "I" in unsupported:
                raise UnsupportedElementError(
                    list(set(unsupported)),
                    message="Iodine (I) is not supported in MVP. "
                    f"Unsupported elements: {', '.join(set(unsupported))}",
                )
            raise UnsupportedElementError(list(set(unsupported)))

        # Check for charged species
        total_charge = Chem.GetFormalCharge(self._rdmol)
        if total_charge != 0:
            raise UnsupportedStructureError(
                f"Charged species not supported (formal charge: {total_charge})",
                smiles=self._smiles,
            )

        # Check for radicals
        for atom in self._rdmol.GetAtoms():
            if atom.GetNumRadicalElectrons() > 0:
                raise UnsupportedStructureError(
                    "Radical species not supported",
                    smiles=self._smiles,
                )

        # Warn for large molecules
        num_atoms = self._rdmol.GetNumAtoms()
        if num_atoms > LARGE_MOLECULE_THRESHOLD:
            warning_msg = (
                f"Large molecule ({num_atoms} atoms). "
                f"Benson GA accuracy may be reduced for molecules >{LARGE_MOLECULE_THRESHOLD} atoms."
            )
            self._warnings.append(warning_msg)
            warnings.warn(warning_msg, UserWarning, stacklevel=3)

    @property
    def rdmol(self) -> Chem.Mol:
        """Access the underlying RDKit Mol object."""
        return self._rdmol

    @property
    def original_smiles(self) -> str | None:
        """Original SMILES string used to create this compound."""
        return self._smiles

    @cached_property
    def canonical_smiles(self) -> str:
        """Canonical SMILES representation."""
        # Remove Hs for canonical SMILES
        mol_no_h = Chem.RemoveHs(self._rdmol)
        return Chem.MolToSmiles(mol_no_h, canonical=True)

    @cached_property
    def inchikey(self) -> str | None:
        """InChIKey for compound lookup."""
        try:
            mol_no_h = Chem.RemoveHs(self._rdmol)
            return inchi.MolToInchiKey(mol_no_h)
        except Exception:
            return None

    @cached_property
    def formula(self) -> str:
        """Molecular formula in Hill notation."""
        return rdMolDescriptors.CalcMolFormula(self._rdmol)

    @cached_property
    def molecular_weight(self) -> float:
        """Molecular weight in g/mol."""
        return Descriptors.MolWt(self._rdmol)

    @cached_property
    def composition(self) -> dict[str, int]:
        """
        Elemental composition as atom counts.

        Returns
        -------
        dict[str, int]
            Mapping of element symbol to atom count
        """
        counts: dict[str, int] = {}
        for atom in self._rdmol.GetAtoms():
            symbol = atom.GetSymbol()
            counts[symbol] = counts.get(symbol, 0) + 1
        return counts

    @cached_property
    def num_atoms(self) -> int:
        """Total number of atoms in the molecule."""
        return self._rdmol.GetNumAtoms()

    @property
    def warnings(self) -> list[str]:
        """List of warnings generated during compound creation."""
        return list(self._warnings)

    @property
    def phase(self) -> str:
        """Phase of the compound ('g', 'l', or 's')."""
        return self._phase

    # =========================================================================
    # Thermodynamic Properties (delegated to phoenix.thermo)
    # =========================================================================

    @property
    def enthalpy_of_formation(self) -> ThermoPropertyAccessor:
        """
        Enthalpy of formation accessor.

        Can be used as property (298.15 K) or method with temperature:

        Examples
        --------
        >>> compound = Compound.from_smiles('CCO')

        >>> # As property (298.15 K) - backward compatible
        >>> hf = compound.enthalpy_of_formation
        >>> print(hf.value)  # kJ/mol

        >>> # As method with temperature
        >>> hf_500 = compound.enthalpy_of_formation(T=500)
        >>> print(hf_500.value)  # kJ/mol at 500 K

        >>> # Vectorized with NumPy
        >>> import numpy as np
        >>> temps = np.linspace(300, 1000, 100)
        >>> values = compound.enthalpy_of_formation(T=temps)  # ndarray
        """
        from phoenix.thermo.benson import estimate_delta_hf

        return ThermoPropertyAccessor(
            compound=self,
            calc_func=lambda c, T: estimate_delta_hf(c, T=T),
            name="enthalpy_of_formation",
        )

    @property
    def delta_hf_kJ_mol(self) -> float:
        """ΔHf° in kJ/mol at 298.15 K (backwards compatible float access)."""
        return float(self.enthalpy_of_formation)

    @property
    def entropy(self) -> ThermoPropertyAccessor:
        """
        Standard entropy accessor.

        Can be used as property (298.15 K) or method with temperature:

        Examples
        --------
        >>> compound = Compound.from_smiles('CCO')

        >>> # As property (298.15 K)
        >>> s = compound.entropy
        >>> print(s.value)  # J/(mol·K)

        >>> # As method with temperature
        >>> s_500 = compound.entropy(T=500)
        >>> print(s_500.value)  # J/(mol·K) at 500 K
        """
        from phoenix.thermo.benson import estimate_entropy

        return ThermoPropertyAccessor(
            compound=self,
            calc_func=lambda c, T: estimate_entropy(c, T=T),
            name="entropy",
        )

    @property
    def entropy_J_mol_K(self) -> float:
        """S° in J/(mol·K) at 298.15 K (backwards compatible float access)."""
        return float(self.entropy)

    def heat_capacity(self, temperature_K: float = 298.15) -> ThermoProperty:
        """
        Heat capacity Cp at specified temperature.

        Parameters
        ----------
        temperature_K : float
            Temperature in Kelvin (default: 298.15 K)

        Returns
        -------
        ThermoProperty
            Heat capacity with value in J/(mol·K), uncertainty, and source
        """
        from phoenix.thermo.benson import estimate_heat_capacity

        return estimate_heat_capacity(self, temperature_K)

    def thermo_at(self, *, T: float) -> ThermoState:
        """
        Get thermodynamic state at specified temperature.

        Returns an immutable ThermoState object with all properties
        (H, S, Cp, G) calculated at the same temperature.

        Parameters
        ----------
        T : float
            Temperature in Kelvin (keyword-only)

        Returns
        -------
        ThermoState
            Immutable state with H, S, Cp, G properties

        Examples
        --------
        >>> compound = Compound.from_smiles('CCO')
        >>> state = compound.thermo_at(T=500)
        >>> print(f"H = {state.H.value:.2f} kJ/mol")
        >>> print(f"S = {state.S.value:.2f} J/(mol·K)")
        >>> print(f"Cp = {state.Cp.value:.2f} J/(mol·K)")
        >>> print(f"G = {state.G.value:.2f} kJ/mol")
        """
        return ThermoState(_compound=self, temperature=T)

    # =========================================================================
    # Hazard Evaluation (delegated to phoenix.hazard)
    # =========================================================================

    @cached_property
    def oxygen_balance(self) -> float:
        """
        Oxygen balance (OB%) for the compound.

        Positive values indicate oxygen excess, negative values oxygen deficiency.
        """
        from phoenix.hazard.oxygen_balance import calculate_oxygen_balance

        return calculate_oxygen_balance(self.composition, self.molecular_weight)

    @property
    def oxygen_balance_percent(self) -> float:
        """Alias for oxygen_balance (backwards compatibility)."""
        return self.oxygen_balance

    def max_decomposition(
        self,
        *,
        method: str = "hierarchy",
        gas_temperature_K: float = 298.15,
    ) -> DecompositionResult:
        """
        Calculate maximum heat of decomposition.

        Parameters
        ----------
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
        >>> result = compound.max_decomposition()  # hierarchy
        >>> result_lp = compound.max_decomposition(method='lp')
        >>> result_both = compound.max_decomposition(method='both')
        >>> print(result.gas_volume_L_g)  # L/g at 298.15 K
        >>> print(result.gas_moles)       # moles gas per mol compound
        >>> print(result.gas_composition) # {"H2": 1.0}
        """
        from phoenix.hazard.decomposition import calculate_max_decomposition

        return calculate_max_decomposition(
            self, method=method, gas_temperature_K=gas_temperature_K
        )

    def evaluate_hazard(self) -> HazardResult:
        """
        Perform full hazard evaluation.

        Returns
        -------
        HazardResult
            Complete hazard assessment with classification and criteria
        """
        from phoenix.hazard.classification import evaluate_hazard

        return evaluate_hazard(self)

    # =========================================================================
    # String Representations
    # =========================================================================

    def __repr__(self) -> str:
        return f"Compound('{self.canonical_smiles}')"

    def __str__(self) -> str:
        return f"{self.formula} (MW={self.molecular_weight:.2f})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Compound):
            return NotImplemented
        return self.canonical_smiles == other.canonical_smiles

    def __hash__(self) -> int:
        return hash(self.canonical_smiles)


# Import decomposition types for type hints
from phoenix.hazard.decomposition import (  # noqa: E402, F401
    DecompositionComparison,
    DecompositionResult,
)
