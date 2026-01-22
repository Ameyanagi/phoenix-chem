"""
Thermodynamic data access layer for PHOENIX.

Uses the `chemicals` library (CalebBell/ChEDL) as the primary data source
for formation enthalpies and other thermodynamic properties.

References:
    - chemicals library: https://github.com/CalebBell/chemicals
    - NIST-JANAF Thermochemical Tables, 4th Ed. (Chase, 1998)
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Literal

# Lazy import to avoid startup overhead
_chemicals_available: bool | None = None


def _check_chemicals() -> bool:
    """Check if chemicals library is available."""
    global _chemicals_available
    if _chemicals_available is None:
        try:
            import chemicals  # noqa: F401

            _chemicals_available = True
        except ImportError:
            _chemicals_available = False
    return _chemicals_available


# =============================================================================
# FALLBACK DATA (when chemicals library unavailable)
# =============================================================================

# Minimal fallback data for core decomposition products
# Source: NIST-JANAF Thermochemical Tables, 4th Ed. (Chase, 1998)
# Units: kJ/mol at 298.15 K
_FALLBACK_HF_GAS: dict[str, float] = {
    "CO2": -393.52,
    "CO": -110.53,
    "H2O": -241.83,
    "H2": 0.0,
    "N2": 0.0,
    "O2": 0.0,
    "NH3": -45.90,
    "H2S": -20.60,
    "SO2": -296.81,
    "SO3": -395.77,
    "HF": -273.30,
    "HCl": -92.31,
    "HBr": -36.29,
    "HI": 26.50,
    "HCN": 135.14,
    "NO": 91.27,
    "NO2": 33.10,
    "N2O": 81.60,
    "CH4": -74.87,
    "CF4": -933.20,
    "F2": 0.0,
    "Cl2": 0.0,
    "Br2": 30.91,
    "I2": 62.42,
}

_FALLBACK_HF_SOLID: dict[str, float] = {
    "C": 0.0,  # graphite
    "S": 0.0,  # rhombic
    "P": 0.0,  # white
    "P4O10": -2984.0,
    "P4O6": -1640.1,
    "I2": 0.0,
}

_FALLBACK_HF_LIQUID: dict[str, float] = {
    "H2O": -285.83,
    "Br2": 0.0,  # reference state
    "H3PO4": -1271.7,
}

_FALLBACK_MW: dict[str, float] = {
    "CO2": 44.01,
    "CO": 28.01,
    "H2O": 18.015,
    "H2": 2.016,
    "N2": 28.014,
    "O2": 32.0,
    "NH3": 17.031,
    "H2S": 34.08,
    "SO2": 64.07,
    "SO3": 80.07,
    "HF": 20.01,
    "HCl": 36.46,
    "HBr": 80.91,
    "HI": 127.91,
    "HCN": 27.03,
    "NO": 30.01,
    "NO2": 46.01,
    "N2O": 44.01,
    "CH4": 16.04,
    "CF4": 88.00,
    "C": 12.011,
    "S": 32.065,
    "P": 30.974,
    "F2": 38.00,
    "Cl2": 70.90,
    "Br2": 159.81,
    "I2": 253.81,
    "P4O10": 283.89,
    "P4O6": 219.89,
    "H3PO4": 98.00,
}

# CAS numbers for common products
_CAS_NUMBERS: dict[str, str] = {
    "CO2": "124-38-9",
    "CO": "630-08-0",
    "H2O": "7732-18-5",
    "H2": "1333-74-0",
    "N2": "7727-37-9",
    "O2": "7782-44-7",
    "NH3": "7664-41-7",
    "H2S": "7783-06-4",
    "SO2": "7446-09-5",
    "SO3": "7446-11-9",
    "HF": "7664-39-3",
    "HCl": "7647-01-0",
    "HBr": "10035-10-6",
    "HI": "10034-85-2",
    "HCN": "74-90-8",
    "NO": "10102-43-9",
    "NO2": "10102-44-0",
    "N2O": "10024-97-2",
    "CH4": "74-82-8",
    "CF4": "75-73-0",
    "C": "7782-42-5",
    "S": "7704-34-9",
    "P": "7723-14-0",
    "F2": "7782-41-4",
    "Cl2": "7782-50-5",
    "Br2": "7726-95-6",
    "I2": "7553-56-2",
    "P4O10": "1314-56-3",
    "P4O6": "12440-00-5",
    "H3PO4": "7664-38-2",
}


# =============================================================================
# DATA CLASSES
# =============================================================================


@dataclass(frozen=True)
class FormationEnthalpy:
    """Formation enthalpy data for a species."""

    formula: str
    gas_kJ_mol: float | None = None
    liquid_kJ_mol: float | None = None
    solid_kJ_mol: float | None = None
    source: str = "unknown"

    def get(self, phase: Literal["gas", "liquid", "solid"] = "gas") -> float | None:
        """Get formation enthalpy for specified phase."""
        if phase == "gas":
            return self.gas_kJ_mol
        elif phase == "liquid":
            return self.liquid_kJ_mol
        elif phase == "solid":
            return self.solid_kJ_mol
        return None


# =============================================================================
# DATA ACCESS FUNCTIONS
# =============================================================================


@lru_cache(maxsize=256)
def get_formation_enthalpy(
    identifier: str,
    use_cas: bool = True,
) -> FormationEnthalpy:
    """
    Get formation enthalpy for a species.

    Uses the `chemicals` library if available, falls back to built-in data.

    Parameters
    ----------
    identifier : str
        Chemical formula (e.g., "CO2") or CAS number (e.g., "124-38-9")
    use_cas : bool
        If True, try to look up by CAS number first (more reliable)

    Returns
    -------
    FormationEnthalpy
        Formation enthalpy data for all available phases

    Examples
    --------
    >>> hf = get_formation_enthalpy("CO2")
    >>> hf.gas_kJ_mol
    -393.52
    >>> hf = get_formation_enthalpy("H2O")
    >>> hf.liquid_kJ_mol
    -285.83
    """
    formula = identifier
    cas = _CAS_NUMBERS.get(identifier)

    if _check_chemicals():
        return _get_from_chemicals(formula, cas, use_cas)
    else:
        return _get_from_fallback(formula)


def _get_from_chemicals(
    formula: str,
    cas: str | None,
    use_cas: bool,
) -> FormationEnthalpy:
    """Fetch formation enthalpy from chemicals library."""
    from chemicals import Hfg, Hfl, Hfs

    lookup_id = cas if (use_cas and cas) else formula

    hf_gas = None
    hf_liquid = None
    hf_solid = None

    # Try gas phase
    try:
        val = Hfg(lookup_id)
        if val is not None:
            hf_gas = val / 1000.0  # J/mol -> kJ/mol
    except Exception:
        pass

    # Try liquid phase
    try:
        val = Hfl(lookup_id)
        if val is not None:
            hf_liquid = val / 1000.0
    except Exception:
        pass

    # Try solid phase
    try:
        val = Hfs(lookup_id)
        if val is not None:
            hf_solid = val / 1000.0
    except Exception:
        pass

    # If chemicals library returned nothing, use fallback
    if hf_gas is None and hf_liquid is None and hf_solid is None:
        return _get_from_fallback(formula)

    return FormationEnthalpy(
        formula=formula,
        gas_kJ_mol=hf_gas,
        liquid_kJ_mol=hf_liquid,
        solid_kJ_mol=hf_solid,
        source="chemicals (CalebBell/ChEDL)",
    )


def _get_from_fallback(formula: str) -> FormationEnthalpy:
    """Get formation enthalpy from fallback data."""
    return FormationEnthalpy(
        formula=formula,
        gas_kJ_mol=_FALLBACK_HF_GAS.get(formula),
        liquid_kJ_mol=_FALLBACK_HF_LIQUID.get(formula),
        solid_kJ_mol=_FALLBACK_HF_SOLID.get(formula),
        source="NIST-JANAF (fallback)",
    )


@lru_cache(maxsize=256)
def get_molecular_weight(identifier: str) -> float | None:
    """
    Get molecular weight for a species.

    Parameters
    ----------
    identifier : str
        Chemical formula or CAS number

    Returns
    -------
    float | None
        Molecular weight in g/mol, or None if not found
    """
    cas = _CAS_NUMBERS.get(identifier)

    if _check_chemicals():
        try:
            from chemicals import MW

            lookup_id = cas if cas else identifier
            return MW(lookup_id)
        except Exception:
            pass

    return _FALLBACK_MW.get(identifier)


@lru_cache(maxsize=256)
def get_entropy(identifier: str, phase: str = "gas") -> float | None:
    """
    Get standard entropy for a species.

    Parameters
    ----------
    identifier : str
        Chemical formula or CAS number
    phase : str
        Phase: "gas", "liquid", or "solid"

    Returns
    -------
    float | None
        Standard entropy in J/(mol·K), or None if not found
    """
    cas = _CAS_NUMBERS.get(identifier)

    if _check_chemicals():
        try:
            from chemicals import S0g, S0l, S0s

            lookup_id = cas if cas else identifier

            if phase == "gas":
                return S0g(lookup_id)
            elif phase == "liquid":
                return S0l(lookup_id)
            elif phase == "solid":
                return S0s(lookup_id)
        except Exception:
            pass

    return None


def get_heat_capacity(
    identifier: str,
    temperature_K: float = 298.15,
    phase: str = "gas",
) -> float | None:
    """
    Get heat capacity at specified temperature.

    Parameters
    ----------
    identifier : str
        Chemical formula or CAS number
    temperature_K : float
        Temperature in Kelvin
    phase : str
        Phase: "gas", "liquid", or "solid"

    Returns
    -------
    float | None
        Heat capacity Cp in J/(mol·K), or None if not found
    """
    cas = _CAS_NUMBERS.get(identifier)

    if _check_chemicals():
        try:
            from chemicals import HeatCapacityGas, HeatCapacityLiquid, HeatCapacitySolid

            lookup_id = cas if cas else identifier

            if phase == "gas":
                cp_obj = HeatCapacityGas(CASRN=lookup_id)
                return cp_obj.T_dependent_property(temperature_K)
            elif phase == "liquid":
                cp_obj = HeatCapacityLiquid(CASRN=lookup_id)
                return cp_obj.T_dependent_property(temperature_K)
            elif phase == "solid":
                cp_obj = HeatCapacitySolid(CASRN=lookup_id)
                return cp_obj.T_dependent_property(temperature_K)
        except Exception:
            pass

    return None


# =============================================================================
# DECOMPOSITION PRODUCTS INTERFACE
# =============================================================================


class DecompositionProducts:
    """
    Interface for accessing thermodynamic data of decomposition products.

    This class provides convenient access to formation enthalpies for
    the standard set of decomposition products used in max ΔHd calculations.
    """

    # Standard product set for analytical decomposition
    PRODUCTS = [
        "CO2",
        "CO",
        "H2O",
        "H2",
        "N2",
        "O2",
        "HF",
        "HCl",
        "HBr",
        "HI",
        "SO2",
        "H2S",
        "NH3",
        "C",
        "S",
        "P",
        "F2",
        "Cl2",
        "Br2",
        "I2",
        "P4O10",
        "CF4",
    ]

    @classmethod
    def get_hf(cls, formula: str, phase: str = "gas") -> float:
        """
        Get formation enthalpy for a decomposition product.

        Parameters
        ----------
        formula : str
            Product formula (e.g., "CO2", "H2O", "C")
        phase : str
            Phase: "gas", "liquid", or "solid"

        Returns
        -------
        float
            Formation enthalpy in kJ/mol

        Raises
        ------
        ValueError
            If product not found or no data for specified phase
        """
        hf = get_formation_enthalpy(formula)
        value = hf.get(phase)

        if value is None:
            # Try alternative phases
            for alt_phase in ["gas", "solid", "liquid"]:
                if alt_phase != phase:
                    value = hf.get(alt_phase)
                    if value is not None:
                        break

        if value is None:
            raise ValueError(f"No formation enthalpy data for {formula} ({phase})")

        return value

    @classmethod
    def get_all_hf(cls, phase: str = "gas") -> dict[str, float]:
        """
        Get formation enthalpies for all standard products.

        Parameters
        ----------
        phase : str
            Preferred phase (falls back to other phases if unavailable)

        Returns
        -------
        dict[str, float]
            Mapping of formula to ΔHf° in kJ/mol
        """
        result = {}
        for formula in cls.PRODUCTS:
            try:
                result[formula] = cls.get_hf(formula, phase)
            except ValueError:
                pass
        return result

    @classmethod
    def is_available(cls) -> bool:
        """Check if chemicals library is available for full data access."""
        return _check_chemicals()


# =============================================================================
# MODULE INITIALIZATION
# =============================================================================


def check_data_sources() -> dict[str, bool]:
    """
    Check availability of data sources.

    Returns
    -------
    dict[str, bool]
        Mapping of source name to availability status
    """
    sources = {
        "chemicals": _check_chemicals(),
        "fallback": True,  # Always available
    }

    # Check janaf if needed for validation
    try:
        import janaf  # noqa: F401

        sources["janaf"] = True
    except ImportError:
        sources["janaf"] = False

    return sources
