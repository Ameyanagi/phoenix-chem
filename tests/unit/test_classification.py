"""Unit tests for hazard classification."""

import pytest

from phoenix import Compound, HazardResult
from phoenix.hazard.classification import (
    CRITERION_1_THRESHOLD,
    CRITERION_2_THRESHOLD,
    evaluate_hazard,
)
from phoenix.hazard.functional_groups import detect_functional_groups, get_alert_names


class TestHazardClassification:
    """Tests for CHETAH hazard classification."""

    def test_hazard_result_structure(self):
        """Test HazardResult has all required fields."""
        compound = Compound.from_smiles("CCO")
        result = evaluate_hazard(compound)

        assert isinstance(result, HazardResult)
        assert hasattr(result, "smiles")
        assert hasattr(result, "formula")
        assert hasattr(result, "molecular_weight")
        assert hasattr(result, "delta_hf_kJ_mol")
        assert hasattr(result, "max_decomposition_kJ_mol")
        assert hasattr(result, "max_decomposition_cal_g")
        assert hasattr(result, "oxygen_balance_percent")
        assert hasattr(result, "hazard_class")
        assert hasattr(result, "triggered_criteria")
        assert hasattr(result, "functional_group_alerts")
        assert hasattr(result, "product_breakdown")
        assert hasattr(result, "gas_volume_L_g")

    def test_hazard_class_values(self):
        """Test hazard class is one of valid values."""
        compound = Compound.from_smiles("CCO")
        result = evaluate_hazard(compound)
        assert result.hazard_class in ("HIGH", "MEDIUM", "LOW")

    def test_low_hazard_compound(self):
        """Test low hazard classification for simple alcohol."""
        compound = Compound.from_smiles("CCO")  # ethanol
        result = evaluate_hazard(compound)
        # Ethanol should be low hazard
        assert result.hazard_class == "LOW" or result.hazard_class == "MEDIUM"

    def test_nitro_compound_triggers_alert(self):
        """Test nitro compounds trigger functional group alert."""
        compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
        result = evaluate_hazard(compound)
        assert 4 in result.triggered_criteria  # Criterion 4 = functional groups
        assert any("Nitro" in alert for alert in result.functional_group_alerts)

    def test_criterion_3_oxygen_balance_range(self):
        """Test criterion 3 triggers for OB in reactive range."""
        # Nitroglycerin has OB ≈ +3.5%
        compound = Compound.from_smiles("C(C(CO[N+](=O)[O-])O[N+](=O)[O-])O[N+](=O)[O-]")
        result = evaluate_hazard(compound)
        # OB between -200% and +100% triggers criterion 3
        assert -200 < result.oxygen_balance_percent < 100
        assert 3 in result.triggered_criteria

    def test_triggered_criteria_is_tuple(self):
        """Test triggered_criteria is immutable tuple."""
        compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
        result = evaluate_hazard(compound)
        assert isinstance(result.triggered_criteria, tuple)

    def test_functional_group_alerts_is_tuple(self):
        """Test functional_group_alerts is immutable tuple."""
        compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
        result = evaluate_hazard(compound)
        assert isinstance(result.functional_group_alerts, tuple)


class TestFunctionalGroupDetection:
    """Tests for functional group alert detection."""

    def test_nitro_group_detection(self):
        """Test nitro group is detected."""
        compound = Compound.from_smiles("c1ccccc1[N+](=O)[O-]")
        alerts = get_alert_names(compound)
        assert any("Nitro" in alert for alert in alerts)

    def test_peroxide_detection(self):
        """Test peroxide is detected."""
        # Benzoyl peroxide
        compound = Compound.from_smiles("c1ccc(cc1)C(=O)OOC(=O)c2ccccc2")
        alerts = get_alert_names(compound)
        assert any("eroxide" in alert.lower() for alert in alerts)

    def test_azo_group_detection(self):
        """Test azo group is detected."""
        # Azobenzene
        compound = Compound.from_smiles("c1ccc(cc1)N=Nc2ccccc2")
        alerts = get_alert_names(compound)
        assert any("Azo" in alert for alert in alerts)

    def test_epoxide_detection(self):
        """Test epoxide is detected."""
        # Ethylene oxide
        compound = Compound.from_smiles("C1CO1")
        alerts = get_alert_names(compound)
        assert any("Epoxide" in alert for alert in alerts)

    def test_no_alerts_for_simple_compound(self):
        """Test no alerts for non-reactive compound."""
        compound = Compound.from_smiles("CCO")  # ethanol
        alerts = get_alert_names(compound)
        # Ethanol has no reactive groups
        assert len(alerts) == 0

    def test_multiple_alerts(self):
        """Test multiple alerts can be detected."""
        # Compound with nitro and another group
        compound = Compound.from_smiles("c1ccc(cc1[N+](=O)[O-])N=Nc2ccccc2")
        alerts = get_alert_names(compound)
        # Should detect both nitro and azo
        assert len(alerts) >= 2

    def test_alert_count(self):
        """Test alert count matches detected groups."""
        compound = Compound.from_smiles("Cc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]")  # TNT
        alerts = detect_functional_groups(compound)
        # TNT has 3 nitro groups
        nitro_alerts = [a for a in alerts if "Nitro" in a.name]
        assert len(nitro_alerts) > 0
        # At least one of them should have count >= 3
        total_nitro = sum(a.count for a in nitro_alerts)
        assert total_nitro >= 3


class TestHazardCriteria:
    """Tests for specific CHETAH criteria."""

    def test_criterion_thresholds(self):
        """Test criterion threshold values are correct."""
        assert CRITERION_1_THRESHOLD == -300.0  # cal/g
        assert CRITERION_2_THRESHOLD == -100.0  # cal/g

    def test_result_includes_ob_percent(self):
        """Test result includes oxygen balance."""
        compound = Compound.from_smiles("CCO")
        result = evaluate_hazard(compound)
        assert result.oxygen_balance_percent is not None
        assert isinstance(result.oxygen_balance_percent, float)

    def test_result_includes_product_breakdown(self):
        """Test result includes decomposition products."""
        compound = Compound.from_smiles("CCO")
        result = evaluate_hazard(compound)
        assert result.product_breakdown is not None
        assert isinstance(result.product_breakdown, dict)
        assert len(result.product_breakdown) > 0


class TestHazardResultImmutability:
    """Tests for HazardResult immutability."""

    def test_result_is_frozen(self):
        """Test HazardResult cannot be modified."""
        compound = Compound.from_smiles("CCO")
        result = evaluate_hazard(compound)

        with pytest.raises(AttributeError):
            result.hazard_class = "HIGH"  # type: ignore

    def test_triggered_criteria_cannot_be_modified(self):
        """Test triggered_criteria tuple is immutable."""
        compound = Compound.from_smiles("CCO")
        result = evaluate_hazard(compound)

        with pytest.raises(TypeError):
            result.triggered_criteria[0] = 1  # type: ignore
