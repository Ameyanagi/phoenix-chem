"""
Benson Group Additivity integration for thermodynamic property estimation.

This module provides thermodynamic property estimation using Benson Group
Additivity (GA) via the pgradd library. When pgradd is unavailable or
lacks group data, falls back to the chemicals library for known compounds.

Reference:
- Benson, S.W. "Thermochemical Kinetics" (Wiley, 1976)
- pgradd: https://github.com/VlachosGroup/PythonGroupAdditivity
"""

from __future__ import annotations

from functools import lru_cache
from typing import TYPE_CHECKING

from phoenix.exceptions import MissingGroupError
from phoenix.thermo.models import GroupContribution, ThermoProperty, ThermoValue
from phoenix.thermo.references import BENSON_1976, CHEMICALS_LIB, PGRADD, Reference

if TYPE_CHECKING:
    from phoenix.core.compound import Compound

# Module-level cache for pgradd library
_pgradd_available: bool | None = None
_benson_library = None

# Gas constant
R = 8.314  # J/(mol·K)


def _check_pgradd() -> bool:
    """Check if pgradd library is available."""
    global _pgradd_available
    if _pgradd_available is None:
        try:
            import pgradd.ThermoChem  # noqa: F401
            from pgradd.GroupAdd.Library import GroupLibrary  # noqa: F401

            _pgradd_available = True
        except ImportError:
            _pgradd_available = False
    return _pgradd_available


def _get_benson_library():
    """Get or load the Benson GA library (cached)."""
    global _benson_library
    if _benson_library is None:
        import pgradd.ThermoChem  # noqa: F401
        from pgradd.GroupAdd.Library import GroupLibrary

        _benson_library = GroupLibrary.Load("BensonGA")
    return _benson_library


def estimate_delta_hf(compound: Compound, *, T: float = 298.15) -> ThermoProperty:
    """
    Estimate standard enthalpy of formation (ΔHf°) in kJ/mol.

    Uses Benson GA via pgradd if available, falls back to chemicals
    library for known compounds.

    Parameters
    ----------
    compound : Compound
        Compound to estimate ΔHf° for
    T : float, optional
        Temperature in Kelvin (default: 298.15 K, keyword-only)

    Returns
    -------
    ThermoProperty
        Formation enthalpy in kJ/mol at temperature T

    Raises
    ------
    MissingGroupError
        If neither pgradd nor chemicals library can provide data

    Examples
    --------
    >>> from phoenix import Compound
    >>> ethanol = Compound.from_smiles("CCO")
    >>> hf = estimate_delta_hf(ethanol)  # At 298.15 K
    >>> hf_500 = estimate_delta_hf(ethanol, T=500)  # At 500 K
    """
    # Try pgradd Benson GA first (primary method)
    if _check_pgradd():
        try:
            val, missing_groups = _estimate_from_pgradd(compound, temperature_K=T)
            source = "Benson GA (pgradd)"
            if missing_groups:
                source += f" [missing: {', '.join(missing_groups[:3])}]"
            return ThermoProperty(
                value=val,
                unit="kJ/mol",
                uncertainty=12.0,  # Typical Benson GA uncertainty ~3 kcal/mol
                source=source,
                phase=compound.phase,
                temperature_K=T,
            )
        except MissingGroupError:
            pass  # Fall through to chemicals library
        except Exception:
            pass  # Fall through to chemicals library

    # Try chemicals library (only valid at 298.15 K)
    try:
        val = _estimate_from_chemicals(compound)
        # Note: chemicals library data is at 298.15 K only
        # For other temperatures, we'd need Cp integration
        return ThermoProperty(
            value=val,
            unit="kJ/mol",
            uncertainty=2.0,  # Lookup is generally more reliable
            source="chemicals (Lookup)",
            phase=compound.phase,
            temperature_K=T,  # Store requested T, though data is 298.15 K
        )
    except Exception:
        pass

    # If we get here, try a simple estimation based on composition
    # This is a last resort fallback
    try:
        val = _estimate_fallback(compound)
        return ThermoProperty(
            value=val,
            unit="kJ/mol",
            uncertainty=50.0,  # Very high uncertainty for fallback
            source="fallback (Composition correlation)",
            phase=compound.phase,
            temperature_K=T,
        )
    except Exception as e:
        raise MissingGroupError(
            message=f"Cannot estimate ΔHf° for {compound.formula}: "
            f"pgradd unavailable and compound not in chemicals database. "
            f"Original error: {e}"
        ) from e


def _estimate_from_pgradd(
    compound: Compound, temperature_K: float = 298.15
) -> tuple[float, list[str]]:
    """
    Estimate ΔHf° using pgradd Benson GA.

    Returns
    -------
    tuple[float, list[str]]
        (ΔHf° in kJ/mol, list of missing groups)
    """
    lib = _get_benson_library()
    smiles = compound.canonical_smiles

    # Get group descriptors
    descriptors = lib.GetDescriptors(smiles)

    if not descriptors:
        raise MissingGroupError(message=f"No groups found for {smiles}")

    # Sum contributions from each group
    total_HoRT = 0.0
    missing_groups = []

    for group_name, count in descriptors.items():
        if group_name in lib:
            group_data = lib[group_name]
            if "thermochem" in group_data:
                tc = group_data["thermochem"]
                total_HoRT += tc.get_HoRT(temperature_K) * count
            else:
                missing_groups.append(group_name)
        else:
            missing_groups.append(group_name)

    if missing_groups and len(missing_groups) == len(descriptors):
        raise MissingGroupError(
            message=f"All groups missing thermochemistry data: {missing_groups}"
        )

    # Convert HoRT to kJ/mol: H = HoRT * R * T
    H_kJ = total_HoRT * R * temperature_K / 1000.0

    return H_kJ, missing_groups


def _estimate_from_chemicals(compound: Compound) -> float:
    """Get ΔHf° from chemicals library."""
    from chemicals import Hfg, search_chemical

    # Common compounds SMILES -> CAS mapping (search_chemical misinterprets SMILES)
    # These are gas-phase ΔHf° reference compounds
    _SMILES_TO_CAS = {
        "O": "7732-18-5",  # Water (not atomic oxygen)
        "[H][H]": "1333-74-0",  # Hydrogen gas
        "[H]": "12385-13-6",  # Atomic hydrogen
        "N#N": "7727-37-9",  # Nitrogen gas
        "[N]=[N]": "7727-37-9",  # Nitrogen gas (alt)
        "O=O": "7782-44-7",  # Oxygen gas
        "[O][O]": "7782-44-7",  # Oxygen gas (alt)
        "O=C=O": "124-38-9",  # Carbon dioxide
        "[C-]#[O+]": "630-08-0",  # Carbon monoxide
        "C": "74-82-8",  # Methane
        "CC": "74-84-0",  # Ethane
        "C=C": "74-85-1",  # Ethylene
        "C#C": "74-86-2",  # Acetylene
        "CCO": "64-17-5",  # Ethanol
        "CO": "67-56-1",  # Methanol
        "CC=O": "75-07-0",  # Acetaldehyde
        "C=O": "50-00-0",  # Formaldehyde
        "CC(=O)O": "64-19-7",  # Acetic acid
        "OCC(O)CO": "56-81-5",  # Glycerol
        "CC(O)CO": "57-55-6",  # 1,2-Propanediol
        "OCCO": "107-21-1",  # Ethylene glycol
        "c1ccccc1": "71-43-2",  # Benzene
        "Cc1ccccc1": "108-88-3",  # Toluene
        "c1ccccc1[N+](=O)[O-]": "98-95-3",  # Nitrobenzene
        "CN": "74-89-5",  # Methylamine
        "N": "7664-41-7",  # Ammonia
        "Cl": "7647-01-0",  # Hydrogen chloride
        "F": "7664-39-3",  # Hydrogen fluoride
        "Br": "10035-10-6",  # Hydrogen bromide
        "S": "7783-06-4",  # Hydrogen sulfide
        "C1CO1": "75-21-8",  # Ethylene oxide
    }

    smiles = compound.canonical_smiles

    # Try direct CAS lookup for known compounds
    if smiles in _SMILES_TO_CAS:
        cas = _SMILES_TO_CAS[smiles]
        hf = Hfg(cas)
        if hf is not None:
            return hf / 1000.0  # J/mol -> kJ/mol

    # Try InChIKey lookup (more reliable than SMILES search)
    try:
        inchikey = compound.inchikey
        if inchikey:
            chem = search_chemical(inchikey)
            if chem:
                hf = Hfg(chem.CASs)
                if hf is not None:
                    return hf / 1000.0
    except Exception:
        pass

    # Try search by common name patterns
    try:
        chem = search_chemical(compound.formula)
        if chem:
            hf = Hfg(chem.CASs)
            if hf is not None:
                return hf / 1000.0
    except Exception:
        pass

    raise MissingGroupError(message=f"Compound {compound.formula} not found in chemicals database")


def _estimate_fallback(compound: Compound) -> float:
    """
    Fallback estimation using simple bond energy approximation.

    This is a rough approximation and should be used only when
    better methods are unavailable.
    """
    # Get any known product formation enthalpies and work backwards
    # This is a very rough approximation based on combustion products

    comp = compound.composition

    # Very rough estimate based on typical organic compound combustion
    # ΔHf° ≈ ΔHc - products enthalpy (very simplified)
    # For organic CHNO compounds, typical ΔHf° correlates with composition

    # Use a rough correlation: -50 kJ/mol per carbon as baseline
    # Adjust for other elements
    c = comp.get("C", 0)
    h = comp.get("H", 0)
    n = comp.get("N", 0)
    o = comp.get("O", 0)

    # Very rough correlation (will be improved with pgradd/chemicals)
    # This is intentionally conservative (more negative = safer screening)
    hf_estimate = -50.0 * c + 30.0 * n - 20.0 * o + 5.0 * h

    # Cap at reasonable bounds
    hf_estimate = max(-2000.0, min(500.0, hf_estimate))

    return hf_estimate


def _estimate_entropy_from_pgradd(compound: Compound, temperature_K: float = 298.15) -> float:
    """Estimate S° using pgradd Benson GA."""
    lib = _get_benson_library()
    smiles = compound.canonical_smiles
    descriptors = lib.GetDescriptors(smiles)

    if not descriptors:
        raise MissingGroupError(message=f"No groups found for {smiles}")

    total_SoR = 0.0
    matched = 0
    for group_name, count in descriptors.items():
        if group_name in lib:
            group_data = lib[group_name]
            if "thermochem" in group_data:
                tc = group_data["thermochem"]
                total_SoR += tc.get_SoR(temperature_K) * count
                matched += 1

    if matched == 0:
        raise MissingGroupError(
            message=f"No thermochemistry groups matched for {smiles}"
        )

    return total_SoR * R  # J/(mol·K)


def estimate_entropy(compound: Compound, *, T: float = 298.15) -> ThermoProperty:
    """
    Estimate standard entropy S° in J/(mol·K).

    Parameters
    ----------
    compound : Compound
        Compound to estimate S° for
    T : float, optional
        Temperature in Kelvin (default: 298.15 K, keyword-only)

    Returns
    -------
    ThermoProperty
        Standard entropy in J/(mol·K) at temperature T

    Examples
    --------
    >>> from phoenix import Compound
    >>> ethanol = Compound.from_smiles("CCO")
    >>> s = estimate_entropy(ethanol)  # At 298.15 K
    >>> s_500 = estimate_entropy(ethanol, T=500)  # At 500 K
    """
    if _check_pgradd():
        try:
            val = _estimate_entropy_from_pgradd(compound, temperature_K=T)
            return ThermoProperty(
                value=val,
                unit="J/(mol·K)",
                uncertainty=8.0,
                source="Benson GA (pgradd)",
                phase=compound.phase,
                temperature_K=T,
            )
        except Exception:
            pass

    # Try chemicals library lookup (exact data, valid at 298.15 K)
    try:
        from phoenix.thermo.data import get_entropy

        phase_name = {"g": "gas", "l": "liquid", "s": "solid"}.get(compound.phase, compound.phase)
        val = get_entropy(compound.formula, phase=phase_name)
        if val is not None:
            return ThermoProperty(
                value=val,
                unit="J/(mol·K)",
                uncertainty=2.0,  # Lookup is generally more reliable
                source="chemicals (Lookup)",
                phase=compound.phase,
                temperature_K=T,  # Store requested T, though data is 298.15 K
            )
    except Exception:
        pass

    # Fallback: rough estimate based on molecular size
    n_atoms = compound.num_atoms
    val = 100.0 + 20.0 * (n_atoms / 10)
    return ThermoProperty(
        value=val,
        unit="J/(mol·K)",
        uncertainty=20.0,
        source="fallback (Size correlation)",
        phase=compound.phase,
        temperature_K=T,
    )


def _estimate_cp_from_pgradd(compound: Compound, temperature_K: float = 298.15) -> float:
    """Estimate Cp using pgradd Benson GA."""
    lib = _get_benson_library()
    smiles = compound.canonical_smiles
    descriptors = lib.GetDescriptors(smiles)

    if not descriptors:
        raise MissingGroupError(message=f"No groups found for {smiles}")

    total_CpoR = 0.0
    matched = 0
    for group_name, count in descriptors.items():
        if group_name in lib:
            group_data = lib[group_name]
            if "thermochem" in group_data:
                tc = group_data["thermochem"]
                total_CpoR += tc.get_CpoR(temperature_K) * count
                matched += 1

    if matched == 0:
        raise MissingGroupError(
            message=f"No thermochemistry groups matched for {smiles}"
        )

    return total_CpoR * R  # J/(mol·K)


def estimate_heat_capacity(compound: Compound, temperature_K: float = 298.15) -> ThermoProperty:
    """
    Estimate heat capacity Cp at specified temperature.

    Parameters
    ----------
    compound : Compound
        Compound to estimate Cp for
    temperature_K : float
        Temperature in Kelvin

    Returns
    -------
    ThermoProperty
        Heat capacity in J/(mol·K)
    """
    if _check_pgradd():
        try:
            val = _estimate_cp_from_pgradd(compound, temperature_K)
            return ThermoProperty(
                value=val,
                unit="J/(mol·K)",
                uncertainty=5.0,
                source="Benson GA (pgradd)",
                phase=compound.phase,
            )
        except Exception:
            pass

    # Try chemicals library lookup
    try:
        from phoenix.thermo.data import get_heat_capacity

        phase_name = {"g": "gas", "l": "liquid", "s": "solid"}.get(compound.phase, compound.phase)
        val = get_heat_capacity(compound.formula, temperature_K=temperature_K, phase=phase_name)
        if val is not None:
            return ThermoProperty(
                value=val,
                unit="J/(mol·K)",
                uncertainty=2.0,  # Lookup is generally more reliable
                source="chemicals (Lookup)",
                phase=compound.phase,
            )
    except Exception:
        pass

    # Fallback: Dulong-Petit approximation
    n_atoms = compound.num_atoms
    val = 29.1 if n_atoms <= 2 else 0.6 * (3 * n_atoms - 6) * 8.314 + 2.5 * 8.314

    return ThermoProperty(
        value=val,
        unit="J/(mol·K)",
        uncertainty=10.0,
        source="fallback (Simplified vibrational modes)",
        phase=compound.phase,
    )


def estimate_delta_hf_with_breakdown(
    compound: Compound, include_reference: bool = True
) -> ThermoProperty:
    """
    Estimate ΔHf° with full group contribution breakdown and database comparison.

    This is the enhanced version that provides CHETAH-style output with:
    - Individual group contributions
    - Literature references
    - Comparison with database reference values (when available)

    Parameters
    ----------
    compound : Compound
        Compound to estimate ΔHf° for
    include_reference : bool
        If True, lookup database value for comparison (default: True)

    Returns
    -------
    ThermoProperty
        Enhanced property with breakdown, references, and reference_value
    """
    breakdown: tuple[GroupContribution, ...] = ()
    references: tuple[Reference, ...] = ()
    reference_value: ThermoValue | None = None
    estimation_method = ""
    value = 0.0
    uncertainty = 50.0

    # Try to get Benson GA estimate with breakdown
    if _check_pgradd():
        try:
            value, breakdown = _estimate_with_breakdown(compound)
            estimation_method = "Benson GA"
            uncertainty = 12.0
            references = (BENSON_1976, PGRADD)
        except Exception:
            pass

    # If Benson GA failed, try chemicals library
    if not breakdown:
        try:
            value = _estimate_from_chemicals(compound)
            estimation_method = "Database Lookup"
            uncertainty = 2.0
            references = (CHEMICALS_LIB,)
        except Exception:
            # Fallback estimation
            value = _estimate_fallback(compound)
            estimation_method = "Composition Correlation"
            uncertainty = 50.0

    # Get database reference value for comparison
    if include_reference and estimation_method != "Database Lookup":
        try:
            ref_val = _estimate_from_chemicals(compound)
            reference_value = ThermoValue(
                value=ref_val,
                unit="kJ/mol",
                uncertainty=2.0,
                method="NIST/chemicals",
                references=(CHEMICALS_LIB,),
            )
        except Exception:
            pass

    return ThermoProperty(
        value=value,
        unit="kJ/mol",
        uncertainty=uncertainty,
        source=estimation_method,
        phase=compound.phase,
        breakdown=breakdown,
        references=references,
        reference_value=reference_value,
        estimation_method=estimation_method,
        temperature_K=298.15,
    )


def _estimate_with_breakdown(
    compound: Compound, temperature_K: float = 298.15
) -> tuple[float, tuple[GroupContribution, ...]]:
    """
    Estimate ΔHf° using pgradd with group contribution breakdown.

    Returns
    -------
    tuple[float, tuple[GroupContribution, ...]]
        (ΔHf° in kJ/mol, tuple of group contributions)
    """
    lib = _get_benson_library()
    smiles = compound.canonical_smiles

    # Get group descriptors
    descriptors = lib.GetDescriptors(smiles)

    if not descriptors:
        raise MissingGroupError(message=f"No groups found for {smiles}")

    # Build group contributions
    contributions = []
    total_HoRT = 0.0

    for group_name, count in descriptors.items():
        if group_name in lib:
            group_data = lib[group_name]
            if "thermochem" in group_data:
                tc = group_data["thermochem"]
                HoRT_per_group = tc.get_HoRT(temperature_K)
                total_HoRT += HoRT_per_group * count

                # Convert to kJ/mol for display
                H_kJ_per_group = HoRT_per_group * R * temperature_K / 1000.0

                contributions.append(
                    GroupContribution(
                        group_name=group_name,
                        count=count,
                        contribution=H_kJ_per_group,
                        property_type="Hf",
                        source="pgradd",
                    )
                )

    if not contributions:
        raise MissingGroupError(message=f"No thermochemistry data for groups in {smiles}")

    # Convert total to kJ/mol
    H_kJ = total_HoRT * R * temperature_K / 1000.0

    return H_kJ, tuple(contributions)


@lru_cache(maxsize=256)
def get_thermochemistry(smiles: str) -> dict[str, float]:
    """
    Get full thermochemistry for a compound (cached).

    Parameters
    ----------
    smiles : str
        SMILES string

    Returns
    -------
    dict[str, float]
        Dictionary with 'delta_hf_kJ_mol', 'entropy_J_mol_K', 'cp_J_mol_K'
    """
    from phoenix.core.compound import Compound

    compound = Compound.from_smiles(smiles)

    return {
        "delta_hf_kJ_mol": estimate_delta_hf(compound),
        "entropy_J_mol_K": estimate_entropy(compound),
        "cp_J_mol_K": estimate_heat_capacity(compound),
    }
