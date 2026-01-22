"""Unit tests for temperature-dependent thermodynamic API."""

import warnings

import numpy as np
import pytest

from phoenix import Compound
from phoenix.thermo.models import (
    TEMP_DEFAULT,
    TEMP_MAX_WARN,
    TEMP_MIN_WARN,
    ThermoProperty,
    ThermoPropertyAccessor,
    ThermoState,
)


class TestThermoState:
    """Tests for ThermoState frozen dataclass."""

    def test_thermo_state_creation(self):
        """Test creating a ThermoState."""
        compound = Compound.from_smiles("CC")  # Ethane
        state = ThermoState(_compound=compound, temperature=500.0)

        assert state.temperature == 500.0

    def test_thermo_state_is_frozen(self):
        """Test that ThermoState is immutable."""
        compound = Compound.from_smiles("CC")
        state = ThermoState(_compound=compound, temperature=500.0)

        with pytest.raises(AttributeError):
            state.temperature = 600.0  # type: ignore

    def test_thermo_state_H(self):
        """Test enthalpy access via ThermoState."""
        compound = Compound.from_smiles("CC")
        state = compound.thermo_at(T=500)

        H = state.H
        assert isinstance(H, ThermoProperty)
        assert H.unit == "kJ/mol"
        assert H.temperature_K == 500.0

    def test_thermo_state_S(self):
        """Test entropy access via ThermoState."""
        compound = Compound.from_smiles("CC")
        state = compound.thermo_at(T=500)

        S = state.S
        assert isinstance(S, ThermoProperty)
        assert S.unit == "J/(mol·K)"
        assert S.temperature_K == 500.0

    def test_thermo_state_Cp(self):
        """Test heat capacity access via ThermoState."""
        compound = Compound.from_smiles("CC")
        state = compound.thermo_at(T=500)

        Cp = state.Cp
        assert isinstance(Cp, ThermoProperty)
        assert Cp.unit == "J/(mol·K)"

    def test_thermo_state_G(self):
        """Test Gibbs free energy calculation."""
        compound = Compound.from_smiles("CC")
        state = compound.thermo_at(T=500)

        G = state.G
        assert isinstance(G, ThermoProperty)
        assert G.unit == "kJ/mol"

        # Verify G = H - T*S (with unit conversion)
        H_kJ = state.H.value
        S_J = state.S.value
        T = state.temperature
        expected_G = H_kJ - T * (S_J / 1000.0)

        assert G.value == pytest.approx(expected_G, rel=1e-6)

    def test_thermo_state_aliases(self):
        """Test property aliases (enthalpy, entropy, etc.)."""
        compound = Compound.from_smiles("CC")
        state = compound.thermo_at(T=500)

        assert state.enthalpy == state.H
        assert state.entropy == state.S
        assert state.heat_capacity == state.Cp
        assert state.gibbs_energy == state.G

    def test_thermo_state_caching(self):
        """Test that ThermoState caches property calculations."""
        compound = Compound.from_smiles("CC")
        state = compound.thermo_at(T=500)

        # Access H multiple times - should be cached
        H1 = state.H
        H2 = state.H

        # Same object due to caching
        assert H1 is H2


class TestThermoPropertyAccessor:
    """Tests for ThermoPropertyAccessor dual behavior."""

    def test_accessor_as_property(self):
        """Test accessing as property returns value at 298.15 K."""
        compound = Compound.from_smiles("CC")
        hf = compound.enthalpy_of_formation

        assert isinstance(hf, ThermoPropertyAccessor)
        assert hf.temperature_K == TEMP_DEFAULT

    def test_accessor_as_method(self):
        """Test calling as method with T parameter."""
        compound = Compound.from_smiles("CC")
        hf_500 = compound.enthalpy_of_formation(T=500)

        assert isinstance(hf_500, ThermoProperty)
        assert hf_500.temperature_K == 500.0

    def test_accessor_float_conversion(self):
        """Test float() returns value at default temperature."""
        compound = Compound.from_smiles("CC")
        hf = compound.enthalpy_of_formation

        val = float(hf)
        assert isinstance(val, float)

    def test_accessor_repr(self):
        """Test string representation."""
        compound = Compound.from_smiles("CC")
        hf = compound.enthalpy_of_formation

        s = repr(hf)
        assert "kJ/mol" in s

    def test_accessor_value_property(self):
        """Test .value property access."""
        compound = Compound.from_smiles("CC")
        hf = compound.enthalpy_of_formation

        assert isinstance(hf.value, float)

    def test_accessor_unit_property(self):
        """Test .unit property access."""
        compound = Compound.from_smiles("CC")
        hf = compound.enthalpy_of_formation

        assert hf.unit == "kJ/mol"

    def test_accessor_keyword_only(self):
        """Test that T must be keyword-only."""
        compound = Compound.from_smiles("CC")

        # This should raise TypeError (positional not allowed)
        with pytest.raises(TypeError):
            compound.enthalpy_of_formation(500)  # type: ignore


class TestBackwardCompatibility:
    """Tests for backward compatibility with existing API."""

    def test_enthalpy_of_formation_property(self):
        """Test that enthalpy_of_formation still works as property."""
        compound = Compound.from_smiles("CCO")
        hf = compound.enthalpy_of_formation

        # Should be usable as before
        assert hf.value < 0  # Ethanol has negative ΔHf°
        assert hf.unit == "kJ/mol"

    def test_delta_hf_kJ_mol(self):
        """Test backward compatible float property."""
        compound = Compound.from_smiles("CCO")
        val = compound.delta_hf_kJ_mol

        assert isinstance(val, float)
        assert val < 0

    def test_entropy_property(self):
        """Test that entropy still works as property."""
        compound = Compound.from_smiles("CCO")
        s = compound.entropy

        assert s.value > 0
        assert s.unit == "J/(mol·K)"

    def test_entropy_J_mol_K(self):
        """Test backward compatible float property."""
        compound = Compound.from_smiles("CCO")
        val = compound.entropy_J_mol_K

        assert isinstance(val, float)
        assert val > 0

    def test_heat_capacity_method(self):
        """Test that heat_capacity() still works."""
        compound = Compound.from_smiles("CCO")
        cp = compound.heat_capacity()

        assert isinstance(cp, ThermoProperty)
        assert cp.value > 0
        assert cp.unit == "J/(mol·K)"


class TestTemperatureValidation:
    """Tests for temperature validation and warnings."""

    def test_low_temperature_warning(self):
        """Test warning for T < 200 K."""
        compound = Compound.from_smiles("CC")

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            compound.enthalpy_of_formation(T=150)

            assert len(w) == 1
            assert "below" in str(w[0].message).lower()
            assert "200" in str(w[0].message)

    def test_high_temperature_warning(self):
        """Test warning for T > 6000 K."""
        compound = Compound.from_smiles("CC")

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            compound.enthalpy_of_formation(T=7000)

            assert len(w) == 1
            assert "above" in str(w[0].message).lower()
            assert "6000" in str(w[0].message)

    def test_normal_temperature_no_warning(self):
        """Test no warning for normal temperature range."""
        compound = Compound.from_smiles("CC")

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            compound.enthalpy_of_formation(T=500)

            # No warnings for 500 K
            assert len(w) == 0


class TestNumpyVectorization:
    """Tests for NumPy array support."""

    def test_enthalpy_array_input(self):
        """Test enthalpy calculation with NumPy array."""
        compound = Compound.from_smiles("CC")
        temps = np.linspace(300, 1000, 10)

        values = compound.enthalpy_of_formation(T=temps)

        assert isinstance(values, np.ndarray)
        assert len(values) == 10
        assert all(isinstance(v, (float, np.floating)) for v in values)

    def test_entropy_array_input(self):
        """Test entropy calculation with NumPy array."""
        compound = Compound.from_smiles("CC")
        temps = np.array([300, 400, 500, 600])

        values = compound.entropy(T=temps)

        assert isinstance(values, np.ndarray)
        assert len(values) == 4

    def test_array_values_vary_with_temperature(self):
        """Test that values actually change with temperature."""
        compound = Compound.from_smiles("CC")
        temps = np.array([300, 500, 800])

        enthalpies = compound.enthalpy_of_formation(T=temps)

        # Values should be different at different temperatures
        # (Benson GA with proper thermochem should show T dependence)
        assert not np.allclose(enthalpies[0], enthalpies[2])


class TestCompoundThermoAt:
    """Tests for Compound.thermo_at() method."""

    def test_thermo_at_returns_state(self):
        """Test that thermo_at returns ThermoState."""
        compound = Compound.from_smiles("CCO")
        state = compound.thermo_at(T=500)

        assert isinstance(state, ThermoState)
        assert state.temperature == 500.0

    def test_thermo_at_keyword_only(self):
        """Test that T is keyword-only."""
        compound = Compound.from_smiles("CCO")

        with pytest.raises(TypeError):
            compound.thermo_at(500)  # type: ignore

    def test_thermo_at_consistency(self):
        """Test that thermo_at and direct accessor give same results."""
        compound = Compound.from_smiles("CCO")

        state = compound.thermo_at(T=500)
        direct_H = compound.enthalpy_of_formation(T=500)

        assert state.H.value == pytest.approx(direct_H.value, rel=1e-6)


class TestGibbsEnergyCalculation:
    """Tests for Gibbs free energy G = H - TS."""

    def test_gibbs_energy_formula(self):
        """Test that G = H - T*S with correct units."""
        compound = Compound.from_smiles("CCO")
        state = compound.thermo_at(T=400)

        H = state.H.value  # kJ/mol
        S = state.S.value  # J/(mol·K)
        T = state.temperature  # K
        G = state.G.value  # kJ/mol

        # G = H - T*S, with S converted from J to kJ
        expected = H - T * (S / 1000.0)

        assert G == pytest.approx(expected, rel=1e-6)

    def test_gibbs_energy_temperature_dependence(self):
        """Test that G changes with temperature."""
        compound = Compound.from_smiles("CCO")

        state_300 = compound.thermo_at(T=300)
        state_600 = compound.thermo_at(T=600)

        # G should be different at different temperatures
        assert state_300.G.value != pytest.approx(state_600.G.value, rel=0.01)
