"""
Thermodynamic property models with CHETAH-style breakdown support.

This module provides data structures for thermodynamic properties with
full provenance tracking, group contribution breakdown, and reference
comparison (Benson GA estimate vs experimental database values).

Temperature-Dependent API
-------------------------
The module provides two patterns for temperature-dependent calculations:

1. Direct method call with T parameter:
   >>> compound.enthalpy_of_formation(T=500)

2. State object for grouped access:
   >>> state = compound.thermo_at(T=500)
   >>> state.H, state.S, state.Cp, state.G
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING, Any, Callable

import numpy as np

if TYPE_CHECKING:
    from phoenix.core.compound import Compound
    from phoenix.thermo.references import Reference


# Temperature validation constants
TEMP_MIN_WARN = 200.0  # K - warn below this
TEMP_MAX_WARN = 6000.0  # K - warn above this
TEMP_DEFAULT = 298.15  # K - standard state


def _validate_temperature(T: float) -> None:
    """Validate temperature and issue warnings for extreme values."""
    if T < TEMP_MIN_WARN:
        warnings.warn(
            f"Temperature {T} K is below {TEMP_MIN_WARN} K. "
            "Benson GA extrapolation may be unreliable.",
            UserWarning,
            stacklevel=3,
        )
    elif T > TEMP_MAX_WARN:
        warnings.warn(
            f"Temperature {T} K is above {TEMP_MAX_WARN} K. "
            "NASA polynomial extrapolation may be unreliable.",
            UserWarning,
            stacklevel=3,
        )


@dataclass(frozen=True)
class GroupContribution:
    """
    Single Benson group contribution to a thermodynamic property.

    Attributes
    ----------
    group_name : str
        Group notation (e.g., "C-(H)3(C)", "Cb-H")
    count : int
        Number of occurrences in molecule
    contribution : float
        Contribution per group (in property units)
    property_type : str
        Property type: "Hf", "S", "Cp"
    source : str
        Data source (e.g., "Benson 1976", "pgradd")
    """

    group_name: str
    count: int
    contribution: float
    property_type: str
    source: str = "pgradd"

    @property
    def total(self) -> float:
        """Total contribution (count * contribution)."""
        return self.count * self.contribution


@dataclass(frozen=True)
class ThermoValue:
    """
    A thermodynamic value with full provenance.

    Used for reference values from databases.
    """

    value: float
    unit: str
    uncertainty: float | None = None
    method: str = ""
    references: tuple[Reference, ...] = ()

    def __float__(self) -> float:
        return self.value


@dataclass
class ThermoProperty:
    """
    Enhanced thermodynamic property with CHETAH-style breakdown support.

    Provides:
    - Value with uncertainty
    - Group contribution breakdown (when available)
    - Literature references
    - Comparison with database reference values

    Backward compatible: all enhanced fields are optional.
    """

    # Original fields (backward compatible)
    value: float
    unit: str
    uncertainty: float | None = None
    source: str | None = None
    phase: str = "g"

    # Enhanced fields
    breakdown: tuple[GroupContribution, ...] = ()
    references: tuple[Reference, ...] = ()
    reference_value: ThermoValue | None = None
    estimation_method: str = ""
    temperature_K: float = 298.15

    def __float__(self) -> float:
        return self.value

    def __repr__(self) -> str:
        s = f"{self.value:.2f} {self.unit}"
        if self.uncertainty is not None:
            s += f" ± {self.uncertainty:.2f}"
        if self.source:
            s += f" ({self.source})"
        return s

    def has_breakdown(self) -> bool:
        """Check if detailed breakdown is available."""
        return len(self.breakdown) > 0

    def has_reference(self) -> bool:
        """Check if reference value is available for comparison."""
        return self.reference_value is not None

    @property
    def deviation(self) -> float | None:
        """Deviation from reference value (estimate - reference)."""
        if self.reference_value is None:
            return None
        return self.value - self.reference_value.value

    @property
    def deviation_percent(self) -> float | None:
        """Percent deviation from reference value."""
        if self.reference_value is None or abs(self.reference_value.value) < 1e-10:
            return None
        return 100.0 * (self.value - self.reference_value.value) / abs(
            self.reference_value.value
        )

    def format_breakdown(self, property_name: str = "ENTHALPY OF FORMATION") -> str:
        """
        Return CHETAH-style breakdown table.

        Parameters
        ----------
        property_name : str
            Name to show in header

        Returns
        -------
        str
            Formatted breakdown table

        Example Output
        --------------
        ENTHALPY OF FORMATION (GAS) CALCULATION
        ========================================
        Group                    Count    Contribution    Total
        -------------------------------------------------------
        C-(H)3(C)                   2       -10.00       -20.00
        C-(H)2(C)2                  1        -5.00        -5.00
        -------------------------------------------------------
        ESTIMATED VALUE (Benson GA):                     -25.00 kJ/mol

        REFERENCE VALUE (NIST WebBook):                  -26.50 kJ/mol
        DEVIATION:                                        +1.50 kJ/mol

        References:
        - Benson, S.W. (1976)
        - NIST Chemistry WebBook, SRD 69
        """
        lines = []

        # Header
        phase_str = {"g": "GAS", "l": "LIQUID", "s": "SOLID"}.get(self.phase, self.phase)
        lines.append(f"{property_name} ({phase_str}) CALCULATION")
        lines.append("=" * 60)

        if self.breakdown:
            # Column headers
            lines.append(
                f"{'Group':<28} {'Count':>5} {'Contribution':>14} {'Total':>10}"
            )
            lines.append("-" * 60)

            # Group contributions
            running_total = 0.0
            for group in self.breakdown:
                total = group.total
                running_total += total
                lines.append(
                    f"{group.group_name:<28} {group.count:>5} "
                    f"{group.contribution:>+14.2f} {total:>+10.2f}"
                )

            lines.append("-" * 60)

            # Estimated value
            method = self.estimation_method or "Benson GA"
            lines.append(f"ESTIMATED VALUE ({method}):{self.value:>+28.2f} {self.unit}")
        else:
            # No breakdown available
            lines.append(f"VALUE: {self.value:>+.2f} {self.unit}")
            if self.source:
                lines.append(f"SOURCE: {self.source}")

        # Reference value comparison
        if self.reference_value is not None:
            lines.append("")
            ref_method = self.reference_value.method or "Database"
            lines.append(
                f"REFERENCE VALUE ({ref_method}):"
                f"{self.reference_value.value:>+23.2f} {self.unit}"
            )
            if self.deviation is not None:
                lines.append(f"DEVIATION:{self.deviation:>+49.2f} {self.unit}")

        # References section
        if self.references:
            lines.append("")
            lines.append("References:")
            for ref in self.references:
                lines.append(f"  - {ref.cite()}")

        return "\n".join(lines)


@dataclass(frozen=True)
class ThermoState:
    """
    Immutable thermodynamic state at a specific temperature.

    Provides grouped access to all temperature-dependent properties
    at a single temperature. Use `Compound.thermo_at(T=...)` to create.

    Attributes
    ----------
    temperature : float
        Temperature in Kelvin

    Properties
    ----------
    H : ThermoProperty
        Enthalpy of formation [kJ/mol]
    S : ThermoProperty
        Standard entropy [J/(mol·K)]
    Cp : ThermoProperty
        Heat capacity at constant pressure [J/(mol·K)]
    G : ThermoProperty
        Gibbs free energy of formation [kJ/mol]

    Examples
    --------
    >>> from phoenix import Compound
    >>> ethanol = Compound.from_smiles('CCO')
    >>> state = ethanol.thermo_at(T=500)
    >>> print(f"H = {state.H.value:.2f} kJ/mol")
    >>> print(f"S = {state.S.value:.2f} J/(mol·K)")
    >>> print(f"G = {state.G.value:.2f} kJ/mol")
    """

    _compound: Any  # Compound, but avoid circular import
    temperature: float

    def __post_init__(self) -> None:
        """Validate temperature on creation."""
        _validate_temperature(self.temperature)

    @cached_property
    def H(self) -> ThermoProperty:
        """Enthalpy of formation at this temperature [kJ/mol]."""
        from phoenix.thermo.benson import estimate_delta_hf

        return estimate_delta_hf(self._compound, T=self.temperature)

    @cached_property
    def S(self) -> ThermoProperty:
        """Standard entropy at this temperature [J/(mol·K)]."""
        from phoenix.thermo.benson import estimate_entropy

        return estimate_entropy(self._compound, T=self.temperature)

    @cached_property
    def Cp(self) -> ThermoProperty:
        """Heat capacity at this temperature [J/(mol·K)]."""
        from phoenix.thermo.benson import estimate_heat_capacity

        return estimate_heat_capacity(self._compound, self.temperature)

    @cached_property
    def G(self) -> ThermoProperty:
        """
        Gibbs free energy of formation at this temperature [kJ/mol].

        Calculated as G = H - T*S, with appropriate unit conversion.
        """
        H_kJ = self.H.value
        S_J = self.S.value
        T = self.temperature

        # G = H - T*S, with S in J/(mol·K) -> kJ/(mol·K)
        G_kJ = H_kJ - T * (S_J / 1000.0)

        return ThermoProperty(
            value=G_kJ,
            unit="kJ/mol",
            uncertainty=None,  # Propagated uncertainty would require correlation
            source=f"G = H - T*S at {T:.1f} K",
            phase=self._compound.phase,
            temperature_K=T,
        )

    @property
    def enthalpy(self) -> ThermoProperty:
        """Alias for H (enthalpy of formation)."""
        return self.H

    @property
    def entropy(self) -> ThermoProperty:
        """Alias for S (standard entropy)."""
        return self.S

    @property
    def heat_capacity(self) -> ThermoProperty:
        """Alias for Cp (heat capacity)."""
        return self.Cp

    @property
    def gibbs_energy(self) -> ThermoProperty:
        """Alias for G (Gibbs free energy)."""
        return self.G


class ThermoPropertyAccessor:
    """
    Enables dual property/method access for temperature-dependent properties.

    When accessed as attribute: returns value at 298.15 K (backward compatible)
    When called as method: returns value at specified temperature

    Examples
    --------
    >>> compound = Compound.from_smiles('CCO')

    >>> # As property (298.15 K)
    >>> hf = compound.enthalpy_of_formation
    >>> print(hf)  # Shows value at 298.15 K

    >>> # As method with temperature
    >>> hf_500 = compound.enthalpy_of_formation(T=500)
    >>> print(hf_500)  # Shows value at 500 K

    >>> # Vectorized with NumPy
    >>> import numpy as np
    >>> temps = np.linspace(300, 1000, 100)
    >>> values = compound.enthalpy_of_formation(T=temps)  # Returns ndarray
    """

    __slots__ = ("_compound", "_calc_func", "_name", "_default_T")

    def __init__(
        self,
        compound: Compound,
        calc_func: Callable[[Compound, float], ThermoProperty],
        name: str,
    ):
        """
        Initialize accessor.

        Parameters
        ----------
        compound : Compound
            The compound to calculate properties for
        calc_func : callable
            Function (compound, T) -> ThermoProperty
        name : str
            Property name for display
        """
        self._compound = compound
        self._calc_func = calc_func
        self._name = name
        self._default_T = TEMP_DEFAULT

    def _get_at_default(self) -> ThermoProperty:
        """Get property at default temperature (cached via compound)."""
        return self._calc_func(self._compound, self._default_T)

    def __repr__(self) -> str:
        """Return string representation at default temperature."""
        return repr(self._get_at_default())

    def __str__(self) -> str:
        """Return string at default temperature."""
        return str(self._get_at_default())

    def __float__(self) -> float:
        """Return float value at default temperature."""
        return float(self._get_at_default().value)

    def __call__(self, *, T: float | np.ndarray = TEMP_DEFAULT) -> ThermoProperty | np.ndarray:
        """
        Calculate property at specified temperature(s).

        Parameters
        ----------
        T : float or numpy.ndarray
            Temperature in Kelvin (keyword-only)

        Returns
        -------
        ThermoProperty or numpy.ndarray
            Single ThermoProperty for scalar T, array of values for array T
        """
        if isinstance(T, np.ndarray):
            # Vectorized calculation
            return np.array([self._calc_func(self._compound, t).value for t in T])

        _validate_temperature(T)
        return self._calc_func(self._compound, T)

    # Forward common attributes to the default-temperature value
    @property
    def value(self) -> float:
        """Property value at default temperature."""
        return self._get_at_default().value

    @property
    def unit(self) -> str:
        """Property unit."""
        return self._get_at_default().unit

    @property
    def uncertainty(self) -> float | None:
        """Property uncertainty at default temperature."""
        return self._get_at_default().uncertainty

    @property
    def source(self) -> str | None:
        """Property source."""
        return self._get_at_default().source

    @property
    def temperature_K(self) -> float:
        """Temperature at which default value is calculated."""
        return self._default_T
