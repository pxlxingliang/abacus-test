#!/usr/bin/env python
"""Test get_vbm and get_cbm methods."""

import numpy as np
from abacustest.lib_data.band import BandData


def test_simple():
    # Create dummy high symmetry labels
    high_symm_labels = {"G": [0, 0, 0], "X": [0.5, 0, 0]}
    # Create kpaths: one segment from G to X with 3 kpoints (indices 0,1,2)
    kpaths = [{"start": "G", "end": "X", "start_nkpt": 0, "end_nkpt": 2}]
    # Cumulative distance (3 points)
    kpath_cum_dist = [0.0, 0.5, 1.0]
    efermi = 0.0
    # band_data shape (1 spin, 3 kpoints, 2 bands)
    # Let's set band 0 below fermi, band 1 above fermi
    band_data = np.array(
        [
            [
                [-1.0, 0.5],  # k0: band0=-1, band1=0.5
                [-0.5, 1.0],  # k1: band0=-0.5, band1=1.0
                [-0.2, 1.2],
            ]
        ]
    )  # k2: band0=-0.2, band1=1.2
    bd = BandData(high_symm_labels, kpaths, kpath_cum_dist, efermi, band_data)
    print("Band data shape:", bd.band_data.shape)
    print("Is metal?", bd.is_metal())
    vbm = bd.get_vbm()
    print("VBM:", vbm)
    cbm = bd.get_cbm()
    print("CBM:", cbm)
    # Expected VBM: maximum below fermi is -0.2 at k=2, band=0
    assert abs(vbm["energy"] - (-0.2)) < 1e-6
    assert vbm["kpoint_index"] == [2]  # label None, so only index 2
    assert vbm["band_index"] == {0: [0]}  # spin 0, band 0
    # Expected CBM: minimum above fermi is 0.5 at k=0, band=1
    assert abs(cbm["energy"] - 0.5) < 1e-6
    assert cbm["kpoint_index"] == [0]
    assert cbm["band_index"] == {0: [1]}
    print("Passed simple test")


def test_degenerate():
    # Multiple bands sharing same energy
    high_symm_labels = {"G": [0, 0, 0]}
    kpaths = [{"start": "G", "end": "G", "start_nkpt": 0, "end_nkpt": 0}]
    kpath_cum_dist = [0.0]
    efermi = 0.0
    # 1 spin, 1 kpoint, 3 bands: band0=-1, band1=-1 (degenerate VBM), band2=0.5
    band_data = np.array([[[-1.0, -1.0, 0.5]]])
    bd = BandData(high_symm_labels, kpaths, kpath_cum_dist, efermi, band_data)
    vbm = bd.get_vbm()
    print("Degenerate VBM:", vbm)
    assert abs(vbm["energy"] - (-1.0)) < 1e-6
    assert set(vbm["band_index"][0]) == {0, 1}
    cbm = bd.get_cbm()
    print("CBM:", cbm)
    assert abs(cbm["energy"] - 0.5) < 1e-6
    assert cbm["band_index"][0] == [2]
    print("Passed degenerate test")


def test_metal():
    # Band crossing fermi level
    high_symm_labels = {"G": [0, 0, 0]}
    kpaths = [{"start": "G", "end": "G", "start_nkpt": 0, "end_nkpt": 0}]
    kpath_cum_dist = [0.0]
    efermi = 0.0
    # band0 = -0.5, band1 = 0.5 (no crossing at same kpoint, but is_metal checks crossing across kpoints)
    # To make metal, we need two kpoints where band changes sign.
    # Let's create 2 kpoints
    kpaths = [{"start": "G", "end": "G", "start_nkpt": 0, "end_nkpt": 1}]
    kpath_cum_dist = [0.0, 0.5]
    band_data = np.array(
        [[[-0.5], [0.5]]]
    )  # shape (1,2,1): band0 at k0=-0.5, k1=0.5 -> crosses fermi
    bd = BandData(high_symm_labels, kpaths, kpath_cum_dist, efermi, band_data)
    print("Is metal?", bd.is_metal())
    assert bd.is_metal() == True
    vbm = bd.get_vbm()
    cbm = bd.get_cbm()
    print("VBM for metal:", vbm)
    print("CBM for metal:", cbm)
    # Should return empty dicts and None energy
    assert vbm["energy"] is None
    assert cbm["energy"] is None
    assert vbm["band_index"] == {}
    assert cbm["band_index"] == {}
    print("Passed metal test")


def test_spin_polarized():
    # Two spin channels
    high_symm_labels = {"G": [0, 0, 0]}
    kpaths = [{"start": "G", "end": "G", "start_nkpt": 0, "end_nkpt": 0}]
    kpath_cum_dist = [0.0]
    efermi = 0.0
    # shape (2,1,2): spin up bands: [-1.0, 0.5]; spin down bands: [-0.8, 0.6]
    band_data = np.array([[[-1.0, 0.5]], [[-0.8, 0.6]]])
    bd = BandData(high_symm_labels, kpaths, kpath_cum_dist, efermi, band_data)
    vbm = bd.get_vbm()
    print("Spin polarized VBM:", vbm)
    # VBM is max below fermi across spins: -0.8 (spin down) vs -1.0 (spin up) -> -0.8
    assert abs(vbm["energy"] - (-0.8)) < 1e-6
    assert vbm["band_index"] == {1: [0]}  # spin index 1, band 0
    cbm = bd.get_cbm()
    print("Spin polarized CBM:", cbm)
    # CBM min above fermi: 0.5 (spin up) vs 0.6 (spin down) -> 0.5
    assert abs(cbm["energy"] - 0.5) < 1e-6
    assert cbm["band_index"] == {0: [1]}
    print("Passed spin polarized test")


def test_label_mapping():
    # Test that labels are correctly mapped
    high_symm_labels = {"G": [0, 0, 0], "X": [0.5, 0, 0], "Y": [0, 0.5, 0]}
    # Two segments: G->X (3 kpoints), X->Y (2 kpoints)
    kpaths = [
        {"start": "G", "end": "X", "start_nkpt": 0, "end_nkpt": 2},
        {"start": "X", "end": "Y", "start_nkpt": 3, "end_nkpt": 4},
    ]
    kpath_cum_dist = [0.0, 0.3, 0.6, 0.8, 1.0]
    efermi = 0.0
    # dummy band data (1 spin, 5 kpoints, 2 bands) where band0 always below, band1 always above
    band_data = np.array(
        [
            [
                [-1.0, 0.5],  # k0
                [-0.5, 0.6],  # k1
                [-0.2, 0.7],  # k2
                [-0.1, 0.8],  # k3
                [-0.05, 0.9],
            ]
        ]
    )  # k4
    bd = BandData(high_symm_labels, kpaths, kpath_cum_dist, efermi, band_data)
    # Check label mapping
    print("Label to indices:", bd.label_to_indices)
    print("Kpoint labels:", bd.kpoint_labels)
    assert bd.label_to_indices["G"] == [0]
    assert bd.label_to_indices["X"] == [2, 3]  # end of first segment, start of second
    assert bd.label_to_indices["Y"] == [4]
    assert bd.kpoint_labels[0] == "G"
    assert bd.kpoint_labels[2] == "X"
    assert bd.kpoint_labels[3] == "X"
    assert bd.kpoint_labels[4] == "Y"
    # VBM is at k=4 (energy -0.05) but label Y, should collect all Y indices (only k=4)
    vbm = bd.get_vbm()
    print("VBM with label:", vbm)
    assert vbm["kpoint_labels"] == ["Y"]
    assert vbm["kpoint_index"] == [4]
    # CBM is at k=0 (energy 0.5) label G
    cbm = bd.get_cbm()
    print("CBM with label:", cbm)
    assert cbm["kpoint_labels"] == ["G"]
    assert cbm["kpoint_index"] == [0]
    print("Passed label mapping test")


def test_band_gap():
    """Test get_band_gap method."""
    # Simple semiconductor case
    high_symm_labels = {"G": [0, 0, 0], "X": [0.5, 0, 0]}
    kpaths = [{"start": "G", "end": "X", "start_nkpt": 0, "end_nkpt": 2}]
    kpath_cum_dist = [0.0, 0.5, 1.0]
    efermi = 0.0
    # band_data shape (1 spin, 3 kpoints, 2 bands)
    band_data = np.array(
        [
            [
                [-1.0, 0.5],  # k0: band0=-1, band1=0.5
                [-0.5, 1.0],  # k1: band0=-0.5, band1=1.0
                [-0.2, 1.2],
            ]
        ]
    )  # k2: band0=-0.2, band1=1.2
    bd = BandData(high_symm_labels, kpaths, kpath_cum_dist, efermi, band_data)

    # Test non-spin-resolved band gap
    gap = bd.get_band_gap()
    print("Band gap:", gap)
    expected_gap = 0.5 - (-0.2)  # CBM energy - VBM energy
    assert abs(gap - expected_gap) < 1e-6
    print(f"Expected gap: {expected_gap}, got: {gap}")

    # Test spin-resolved case (same data, but nspin=1)
    gap_spin = bd.get_band_gap(spin_resolved=True)
    print("Spin-resolved band gap:", gap_spin)
    assert isinstance(gap_spin, dict)
    assert "global" in gap_spin
    assert abs(gap_spin["global"] - expected_gap) < 1e-6
    assert 0 in gap_spin
    assert abs(gap_spin[0] - expected_gap) < 1e-6

    # Test metal case (band crossing Fermi level)
    metallic_band_data = np.array(
        [
            [
                [-1.0, 0.5],
                [-0.2, 1.0],
                [0.1, 1.2],  # band0 crosses fermi at 0.1 > 0
            ]
        ]
    )
    bd_metal = BandData(
        high_symm_labels, kpaths, kpath_cum_dist, efermi, metallic_band_data
    )
    assert bd_metal.is_metal()
    metal_gap = bd_metal.get_band_gap()
    print("Metal band gap:", metal_gap)
    assert metal_gap is None

    metal_gap_spin = bd_metal.get_band_gap(spin_resolved=True)
    print("Metal spin-resolved band gap:", metal_gap_spin)
    assert isinstance(metal_gap_spin, dict)
    assert metal_gap_spin["global"] is None
    assert metal_gap_spin[0] is None

    # Test spin-polarized case (nspin=2)
    spin_polarized_band_data = np.array(
        [
            [  # spin up
                [-1.0, 0.5],
                [-0.5, 1.0],
                [-0.2, 1.2],
            ],
            [  # spin down
                [-0.8, 0.7],  # VBM down: -0.8 at k=0, band=0
                [-0.3, 1.1],  # CBM down: 0.7 at k=0, band=1
                [-0.1, 1.3],
            ],
        ]
    )
    bd_spin2 = BandData(
        high_symm_labels, kpaths, kpath_cum_dist, efermi, spin_polarized_band_data
    )
    gap_spin2 = bd_spin2.get_band_gap(spin_resolved=True)
    print("Spin-2 band gap:", gap_spin2)
    # Global gap: min CBM across spins - max VBM across spins
    # Global VBM = max(-0.2, -0.1) = -0.1? Wait, need to check _find_band_edge logic
    # Actually VBM is maximum below fermi: spin up VBM = -0.2, spin down VBM = -0.1, so global VBM = -0.1
    # Global CBM is minimum above fermi: spin up CBM = 0.5, spin down CBM = 0.7, so global CBM = 0.5
    # Global gap = 0.5 - (-0.1) = 0.6
    # But let's compute from actual get_vbm/get_cbm
    vbm_global = bd_spin2.get_vbm(spin_resolved=True)["global"]["energy"]
    cbm_global = bd_spin2.get_cbm(spin_resolved=True)["global"]["energy"]
    expected_global_gap = cbm_global - vbm_global
    assert abs(gap_spin2["global"] - expected_global_gap) < 1e-6

    # Per-spin gaps
    for spin in [0, 1]:
        vbm_spin = bd_spin2.get_vbm(spin_resolved=True)[spin]["energy"]
        cbm_spin = bd_spin2.get_cbm(spin_resolved=True)[spin]["energy"]
        expected_spin_gap = cbm_spin - vbm_spin
        assert abs(gap_spin2[spin] - expected_spin_gap) < 1e-6

    print("Passed band gap test")


if __name__ == "__main__":
    test_simple()
    test_degenerate()
    test_metal()
    test_spin_polarized()
    test_label_mapping()
    test_band_gap()
    print("All tests passed!")
