"""Unit tests for thermodynamic models with CHETAH-style breakdown."""

import pytest

from phoenix.thermo.models import GroupContribution, ThermoProperty, ThermoValue
from phoenix.thermo.references import (
    BENSON_1976,
    NIST_WEBBOOK,
    Reference,
    StandardReference,
    get_reference,
)


class TestReference:
    """Tests for Reference dataclass."""

    def test_reference_creation(self):
        """Test creating a reference."""
        ref = Reference(
            author="Test Author",
            year=2024,
            title="Test Title",
            source_type="journal",
            doi="10.1234/test",
        )
        assert ref.author == "Test Author"
        assert ref.year == 2024
        assert ref.doi == "10.1234/test"

    def test_reference_cite(self):
        """Test short citation format."""
        assert BENSON_1976.cite() == "Benson, S.W. (1976)"
        assert NIST_WEBBOOK.cite() == "NIST (2024)"

    def test_reference_cite_full(self):
        """Test full citation format."""
        full = BENSON_1976.cite_full()
        assert "Benson" in full
        assert "1976" in full
        assert "Thermochemical Kinetics" in full

    def test_standard_reference_enum(self):
        """Test StandardReference enum access."""
        assert StandardReference.BENSON_1976.value == BENSON_1976
        assert StandardReference.NIST_WEBBOOK.value == NIST_WEBBOOK

    def test_get_reference(self):
        """Test get_reference helper."""
        ref = get_reference("BENSON_1976")
        assert ref == BENSON_1976

        ref_lower = get_reference("nist_webbook")
        assert ref_lower == NIST_WEBBOOK

    def test_get_reference_unknown(self):
        """Test get_reference with unknown key."""
        with pytest.raises(ValueError, match="Unknown reference"):
            get_reference("UNKNOWN_REF")


class TestGroupContribution:
    """Tests for GroupContribution dataclass."""

    def test_group_contribution_creation(self):
        """Test creating a group contribution."""
        group = GroupContribution(
            group_name="C-(H)3(C)",
            count=2,
            contribution=-10.5,
            property_type="Hf",
            source="pgradd",
        )
        assert group.group_name == "C-(H)3(C)"
        assert group.count == 2
        assert group.contribution == -10.5

    def test_group_contribution_total(self):
        """Test total property calculation."""
        group = GroupContribution(
            group_name="C-(H)3(C)",
            count=3,
            contribution=-10.0,
            property_type="Hf",
        )
        assert group.total == -30.0

    def test_group_contribution_frozen(self):
        """Test that GroupContribution is immutable."""
        group = GroupContribution(
            group_name="Test",
            count=1,
            contribution=5.0,
            property_type="Hf",
        )
        with pytest.raises(AttributeError):
            group.count = 2


class TestThermoValue:
    """Tests for ThermoValue dataclass."""

    def test_thermo_value_creation(self):
        """Test creating a ThermoValue."""
        val = ThermoValue(
            value=-234.5,
            unit="kJ/mol",
            uncertainty=2.0,
            method="NIST",
            references=(NIST_WEBBOOK,),
        )
        assert val.value == -234.5
        assert val.unit == "kJ/mol"
        assert val.method == "NIST"

    def test_thermo_value_float(self):
        """Test float conversion."""
        val = ThermoValue(value=-100.5, unit="kJ/mol")
        assert float(val) == -100.5


class TestThermoProperty:
    """Tests for enhanced ThermoProperty."""

    def test_thermo_property_basic(self):
        """Test basic ThermoProperty creation (backward compatible)."""
        prop = ThermoProperty(
            value=-234.0,
            unit="kJ/mol",
            uncertainty=5.0,
            source="test",
            phase="g",
        )
        assert prop.value == -234.0
        assert float(prop) == -234.0

    def test_thermo_property_repr(self):
        """Test string representation."""
        prop = ThermoProperty(
            value=-234.0,
            unit="kJ/mol",
            uncertainty=5.0,
            source="test",
        )
        s = repr(prop)
        assert "-234.00" in s
        assert "kJ/mol" in s
        assert "± 5.00" in s
        assert "test" in s

    def test_has_breakdown_false(self):
        """Test has_breakdown when no breakdown."""
        prop = ThermoProperty(value=-100.0, unit="kJ/mol")
        assert not prop.has_breakdown()

    def test_has_breakdown_true(self):
        """Test has_breakdown when breakdown present."""
        group = GroupContribution("Test", 1, -50.0, "Hf")
        prop = ThermoProperty(
            value=-100.0,
            unit="kJ/mol",
            breakdown=(group,),
        )
        assert prop.has_breakdown()

    def test_has_reference(self):
        """Test has_reference method."""
        prop = ThermoProperty(value=-100.0, unit="kJ/mol")
        assert not prop.has_reference()

        ref_val = ThermoValue(value=-95.0, unit="kJ/mol")
        prop_with_ref = ThermoProperty(
            value=-100.0,
            unit="kJ/mol",
            reference_value=ref_val,
        )
        assert prop_with_ref.has_reference()

    def test_deviation(self):
        """Test deviation calculation."""
        ref_val = ThermoValue(value=-95.0, unit="kJ/mol")
        prop = ThermoProperty(
            value=-100.0,
            unit="kJ/mol",
            reference_value=ref_val,
        )
        assert prop.deviation == -5.0

    def test_deviation_none(self):
        """Test deviation when no reference."""
        prop = ThermoProperty(value=-100.0, unit="kJ/mol")
        assert prop.deviation is None

    def test_deviation_percent(self):
        """Test percent deviation calculation."""
        ref_val = ThermoValue(value=-100.0, unit="kJ/mol")
        prop = ThermoProperty(
            value=-105.0,
            unit="kJ/mol",
            reference_value=ref_val,
        )
        assert prop.deviation_percent == pytest.approx(-5.0)

    def test_format_breakdown_no_groups(self):
        """Test format_breakdown without group contributions."""
        prop = ThermoProperty(
            value=-234.0,
            unit="kJ/mol",
            source="Test Source",
            phase="g",
        )
        output = prop.format_breakdown()
        assert "ENTHALPY OF FORMATION (GAS)" in output
        assert "VALUE:" in output
        assert "-234.00" in output

    def test_format_breakdown_with_groups(self):
        """Test format_breakdown with group contributions."""
        groups = (
            GroupContribution("C-(H)3(C)", 2, -10.0, "Hf"),
            GroupContribution("C-(H)2(O)", 1, -5.0, "Hf"),
        )
        prop = ThermoProperty(
            value=-25.0,
            unit="kJ/mol",
            breakdown=groups,
            estimation_method="Benson GA",
            references=(BENSON_1976,),
        )
        output = prop.format_breakdown()

        assert "ENTHALPY OF FORMATION" in output
        assert "C-(H)3(C)" in output
        assert "C-(H)2(O)" in output
        assert "ESTIMATED VALUE" in output
        assert "Benson GA" in output
        assert "References:" in output
        assert "Benson" in output

    def test_format_breakdown_with_reference(self):
        """Test format_breakdown with reference value comparison."""
        ref_val = ThermoValue(
            value=-26.0,
            unit="kJ/mol",
            method="NIST",
        )
        prop = ThermoProperty(
            value=-25.0,
            unit="kJ/mol",
            reference_value=ref_val,
        )
        output = prop.format_breakdown()

        assert "REFERENCE VALUE" in output
        assert "DEVIATION" in output
        assert "+1.00" in output  # deviation = -25 - (-26) = +1
