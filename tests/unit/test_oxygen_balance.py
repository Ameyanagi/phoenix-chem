"""Unit tests for oxygen balance calculation."""

from phoenix.hazard.oxygen_balance import calculate_oxygen_balance


class TestOxygenBalance:
    """Tests for oxygen balance calculation."""

    def test_tnt_oxygen_balance(self):
        """Test OB% for TNT (C7H5N3O6, MW=227.13)."""
        comp = {"C": 7, "H": 5, "N": 3, "O": 6}
        mw = 227.13
        ob = calculate_oxygen_balance(comp, mw)
        # TNT is oxygen-deficient, OB ≈ -74%
        assert -76 < ob < -72

    def test_nitroglycerin_oxygen_balance(self):
        """Test OB% for nitroglycerin (C3H5N3O9, MW=227.09)."""
        comp = {"C": 3, "H": 5, "N": 3, "O": 9}
        mw = 227.09
        ob = calculate_oxygen_balance(comp, mw)
        # Nitroglycerin is slightly oxygen-positive, OB ≈ +3.5%
        assert 2 < ob < 5

    def test_ammonium_nitrate_oxygen_balance(self):
        """Test OB% for ammonium nitrate (H4N2O3, MW=80.04)."""
        comp = {"H": 4, "N": 2, "O": 3}
        mw = 80.04
        ob = calculate_oxygen_balance(comp, mw)
        # NH4NO3 is oxygen-positive, OB ≈ +20%
        assert 18 < ob < 22

    def test_ethanol_oxygen_balance(self):
        """Test OB% for ethanol (C2H6O, MW=46.07)."""
        comp = {"C": 2, "H": 6, "O": 1}
        mw = 46.07
        ob = calculate_oxygen_balance(comp, mw)
        # Ethanol is very oxygen-deficient, OB ≈ -209%
        assert -215 < ob < -205

    def test_methane_oxygen_balance(self):
        """Test OB% for methane (CH4, MW=16.04)."""
        comp = {"C": 1, "H": 4}
        mw = 16.04
        ob = calculate_oxygen_balance(comp, mw)
        # Methane needs 4 O for complete combustion, OB ≈ -400%
        assert -410 < ob < -390

    def test_carbon_dioxide_oxygen_balance(self):
        """Test OB% for CO2 (MW=44.01)."""
        comp = {"C": 1, "O": 2}
        mw = 44.01
        ob = calculate_oxygen_balance(comp, mw)
        # CO2 is fully oxidized, OB = 0% (no H to form water)
        # Actually: OB = -1600/44.01 * (2*1 - 2) = 0
        assert abs(ob) < 1

    def test_water_oxygen_balance(self):
        """Test OB% for water (H2O, MW=18.015)."""
        comp = {"H": 2, "O": 1}
        mw = 18.015
        ob = calculate_oxygen_balance(comp, mw)
        # Water: OB = -1600/18.015 * (0.5*2 - 1) = 0
        assert abs(ob) < 1

    def test_sulfur_containing_compound(self):
        """Test OB% with sulfur (needs 2O per S for SO2)."""
        # H2S: needs 3O total (1 for H2O, 2 for SO2)
        comp = {"H": 2, "S": 1}
        mw = 34.08
        ob = calculate_oxygen_balance(comp, mw)
        # OB = -1600/34.08 * (0.5*2 + 2*1 - 0) = -1600/34.08 * 3 = -140.8%
        assert -145 < ob < -135

    def test_phosphorus_containing_compound(self):
        """Test OB% with phosphorus (needs 2.5O per P for P4O10)."""
        # PH3: needs 4O total (1.5 for H2O, 2.5 for P4O10 equivalent)
        comp = {"H": 3, "P": 1}
        mw = 34.0
        ob = calculate_oxygen_balance(comp, mw)
        # OB = -1600/34.0 * (0.5*3 + 2.5*1) = -1600/34.0 * 4 = -188%
        assert -195 < ob < -180

    def test_chlorine_releases_oxygen(self):
        """Test that chlorine reduces O requirement (forms HCl)."""
        # CH3Cl vs CH4: Cl replaces H, reducing H2O formation
        ch4_ob = calculate_oxygen_balance({"C": 1, "H": 4}, 16.04)
        ch3cl_ob = calculate_oxygen_balance({"C": 1, "H": 3, "Cl": 1}, 50.49)
        # CH3Cl should be less negative (needs less O)
        assert ch3cl_ob > ch4_ob

    def test_bromine_releases_oxygen(self):
        """Test that bromine reduces O requirement (forms HBr)."""
        ch4_ob = calculate_oxygen_balance({"C": 1, "H": 4}, 16.04)
        ch3br_ob = calculate_oxygen_balance({"C": 1, "H": 3, "Br": 1}, 94.94)
        # CH3Br should be less negative
        assert ch3br_ob > ch4_ob

    def test_exact_stoichiometry(self):
        """Test exact stoichiometric calculation."""
        # TNT: C7H5N3O6
        # O needed = 2*7 + 0.5*5 = 14 + 2.5 = 16.5
        # O available = 6
        # OB = (6 - 16.5) * 1600 / 227.13 = -73.96%
        comp = {"C": 7, "H": 5, "N": 3, "O": 6}
        mw = 227.13
        ob = calculate_oxygen_balance(comp, mw)
        expected = (6 - 16.5) * 1600 / 227.13
        assert abs(ob - expected) < 0.1
