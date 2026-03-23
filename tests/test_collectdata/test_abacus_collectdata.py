import unittest
from pathlib import Path

from abacustest.lib_collectdata.collectdata import RESULT


class TestAbacusCollectdata(unittest.TestCase):
    """Test Abacus result collection from SCF calculation"""

    def test_abacus_scf(self):
        """Test all extractable metrics from abacus-scf example"""
        ab = RESULT(path=Path(__file__).parent / "abacus-scf", fmt="abacus")

        input_dict = ab["INPUT"]
        self.assertIsInstance(input_dict, dict)
        self.assertEqual(input_dict["calculation"], "scf")
        self.assertEqual(input_dict["esolver_type"], "ksdft")
        self.assertEqual(input_dict["kspacing"], "0.1 0.1 0.1")
        self.assertAlmostEqual(ab["absolute_mag"], 0.000300325, places=3)
        self.assertEqual(ab["absolute_mags"], [0.000300325])
        self.assertEqual(len(ab["atomlabel_list"]), 27)
        self.assertEqual(ab["atomlabel_list"][0], "C")
        self.assertEqual(ab["atomlabel_list"][len(ab["atomlabel_list"]) // 2], "Sb")
        self.assertEqual(ab["atomlabel_list"][-1], "Y")
        self.assertEqual(len(ab["band"]), 2)
        self.assertEqual(len(ab["band"][0]), 26)
        self.assertEqual(len(ab["band"][0][0]), 163)
        band_flat = [item for sublist in ab["band"] for subsublist in sublist for item in subsublist]
        self.assertEqual(len(band_flat), 8476)
        self.assertAlmostEqual(band_flat[0], -84.8346, places=3)
        self.assertAlmostEqual(band_flat[len(band_flat) // 2], -84.8346, places=3)
        self.assertAlmostEqual(band_flat[-1], 4.92066, places=3)
        self.assertEqual(ab["band_gap"], 0)
        # band_plot: skipped
        self.assertEqual(len(ab["band_weight"]), 2)
        self.assertEqual(len(ab["band_weight"][0]), 26)
        self.assertEqual(len(ab["band_weight"][0][0]), 163)
        band_weight_flat = [item for sublist in ab["band_weight"] for subsublist in sublist for item in subsublist]
        self.assertEqual(len(band_weight_flat), 8476)
        self.assertAlmostEqual(band_weight_flat[0], 0.02, places=3)
        self.assertAlmostEqual(band_weight_flat[len(band_weight_flat) // 2], 0.02, places=3)
        self.assertAlmostEqual(band_weight_flat[-1], 0.0, places=3)
        self.assertEqual(len(ab["cell"]), 3)
        self.assertEqual(len(ab["cell"][0]), 3)
        cell_flat = [item for sublist in ab["cell"] for item in sublist]
        self.assertEqual(len(cell_flat), 9)
        self.assertAlmostEqual(cell_flat[0], 8.3004, places=3)
        self.assertAlmostEqual(cell_flat[len(cell_flat) // 2], 7.18836, places=3)
        self.assertAlmostEqual(cell_flat[-1], 29.2889, places=3)
        self.assertEqual(len(ab["cell_init"]), 3)
        self.assertEqual(len(ab["cell_init"][0]), 3)
        cell_init_flat = [item for sublist in ab["cell_init"] for item in sublist]
        self.assertEqual(len(cell_init_flat), 9)
        self.assertAlmostEqual(cell_init_flat[0], 8.3004, places=3)
        self.assertAlmostEqual(cell_init_flat[len(cell_init_flat) // 2], 7.18836, places=3)
        self.assertAlmostEqual(cell_init_flat[-1], 29.2889, places=3)
        self.assertEqual(len(ab["cells"]), 1)
        self.assertEqual(len(ab["cells"][0]), 3)
        self.assertEqual(len(ab["cells"][0][0]), 3)
        cells_flat = [item for sublist in ab["cells"] for subsublist in sublist for item in subsublist]
        self.assertEqual(len(cells_flat), 9)
        self.assertAlmostEqual(cells_flat[0], 8.3004, places=3)
        self.assertAlmostEqual(cells_flat[len(cells_flat) // 2], 7.18836, places=3)
        self.assertAlmostEqual(cells_flat[-1], 29.2889, places=3)
        self.assertEqual(ab["converge"], True)
        self.assertEqual(len(ab["coordinate"]), 27)
        self.assertEqual(len(ab["coordinate"][0]), 3)
        coordinate_flat = [item for sublist in ab["coordinate"] for item in sublist]
        self.assertEqual(len(coordinate_flat), 81)
        self.assertAlmostEqual(coordinate_flat[0], 3.9536412166, places=3)
        self.assertAlmostEqual(coordinate_flat[len(coordinate_flat) // 2], 3.9520368576, places=3)
        self.assertAlmostEqual(coordinate_flat[-1], 13.598405838, places=3)
        self.assertEqual(len(ab["coordinate_init"]), 27)
        self.assertEqual(len(ab["coordinate_init"][0]), 3)
        coordinate_init_flat = [item for sublist in ab["coordinate_init"] for item in sublist]
        self.assertEqual(len(coordinate_init_flat), 81)
        self.assertAlmostEqual(coordinate_init_flat[0], 3.9536412166, places=3)
        self.assertAlmostEqual(coordinate_init_flat[len(coordinate_init_flat) // 2], 3.9520368576, places=3)
        self.assertAlmostEqual(coordinate_init_flat[-1], 13.598405838, places=3)
        self.assertEqual(len(ab["coordinates"]), 1)
        self.assertEqual(len(ab["coordinates"][0]), 27)
        self.assertEqual(len(ab["coordinates"][0][0]), 3)
        coordinates_flat = [item for sublist in ab["coordinates"] for subsublist in sublist for item in subsublist]
        self.assertEqual(len(coordinates_flat), 81)
        self.assertAlmostEqual(coordinates_flat[0], 3.9536412166, places=3)
        self.assertAlmostEqual(coordinates_flat[len(coordinates_flat) // 2], 3.9520368576, places=3)
        self.assertAlmostEqual(coordinates_flat[-1], 13.598405838, places=3)
        self.assertEqual(len(ab["denergy"]), 56)
        self.assertAlmostEqual(ab["denergy"][0], 0.0, places=3)
        self.assertAlmostEqual(ab["denergy"][len(ab["denergy"]) // 2], 3.73493806e-06, places=3)
        self.assertAlmostEqual(ab["denergy"][-1], -9.48614623e-08, places=3)
        self.assertAlmostEqual(ab["denergy_last"], -9.48614623e-08, places=3)
        self.assertEqual(len(ab["drho"]), 56)
        self.assertAlmostEqual(ab["drho"][0], 0.149730032735, places=3)
        self.assertAlmostEqual(ab["drho"][len(ab["drho"]) // 2], 9.38740190449e-06, places=3)
        self.assertAlmostEqual(ab["drho"][-1], 8.73780825487e-08, places=3)
        self.assertAlmostEqual(ab["drho_last"], 8.73780825487e-08, places=3)
        self.assertEqual(ab["ds_lambda_rms"], [])
        self.assertEqual(ab["ds_lambda_step"], [])
        self.assertAlmostEqual(ab["efermi"], 1.7590656467, places=3)
        self.assertEqual(len(ab["element"]), 27)
        self.assertEqual(ab["element"][0], "C")
        self.assertEqual(ab["element"][len(ab["element"]) // 2], "Sb")
        self.assertEqual(ab["element"][-1], "Y")
        self.assertEqual(len(ab["element_list"]), 27)
        self.assertEqual(ab["element_list"][0], "C")
        self.assertEqual(ab["element_list"][len(ab["element_list"]) // 2], "Sb")
        self.assertEqual(ab["element_list"][-1], "Y")
        self.assertEqual(ab["energies"], [-28364.401228])
        self.assertAlmostEqual(ab["energy"], -28364.40122753045, places=3)
        self.assertAlmostEqual(ab["energy_ks"], -28364.4012275304, places=3)
        self.assertAlmostEqual(ab["energy_per_atom"], -1050.5333787974241, places=3)
        self.assertEqual(ab["fft_grid"], [100.0, 100.0, 360.0])
        self.assertEqual(len(ab["force"]), 81)
        self.assertAlmostEqual(ab["force"][0], -0.0725480006, places=3)
        self.assertAlmostEqual(ab["force"][len(ab["force"]) // 2], 0.0020492624, places=3)
        self.assertAlmostEqual(ab["force"][-1], -0.0046602777, places=3)
        self.assertEqual(len(ab["forces"]), 1)
        self.assertEqual(len(ab["forces"][0]), 81)
        forces_flat = [item for sublist in ab["forces"] for item in sublist]
        self.assertEqual(len(forces_flat), 81)
        self.assertAlmostEqual(forces_flat[0], -0.0725480006, places=3)
        self.assertAlmostEqual(forces_flat[len(forces_flat) // 2], 0.0020492624, places=3)
        self.assertAlmostEqual(forces_flat[-1], -0.0046602777, places=3)
        self.assertEqual(ab["ibzk"], 26)
        self.assertEqual(ab["kpt"], [5, 5, 2])
        self.assertEqual(len(ab["label"]), 27)
        self.assertEqual(ab["label"][0], "C")
        self.assertEqual(ab["label"][len(ab["label"]) // 2], "Sb")
        self.assertEqual(ab["label"][-1], "Y")
        self.assertEqual(ab["largest_gradient"], [1.542298352])
        self.assertEqual(ab["largest_gradient_stress"], [52.8786932479])
        self.assertEqual(len(ab["lattice_constant"]), 6)
        self.assertAlmostEqual(ab["lattice_constant"][0], 8.3004, places=3)
        self.assertAlmostEqual(ab["lattice_constant"][len(ab["lattice_constant"]) // 2], 90.0, places=3)
        self.assertAlmostEqual(ab["lattice_constant"][-1], 60.00000945136993, places=3)
        self.assertEqual(len(ab["lattice_constants"]), 1)
        self.assertEqual(len(ab["lattice_constants"][0]), 6)
        lattice_constants_flat = [item for sublist in ab["lattice_constants"] for item in sublist]
        self.assertEqual(len(lattice_constants_flat), 6)
        self.assertAlmostEqual(lattice_constants_flat[0], 8.3004, places=3)
        self.assertAlmostEqual(lattice_constants_flat[len(lattice_constants_flat) // 2], 90.0, places=3)
        self.assertAlmostEqual(lattice_constants_flat[-1], 60.00000945136993, places=3)
        self.assertEqual(ab["natom"], 27)
        self.assertEqual(ab["nbands"], 163)
        self.assertEqual(ab["nbase"], 547)
        self.assertEqual(ab["ncore"], 16)
        self.assertAlmostEqual(ab["nelec"], 244.0, places=3)
        self.assertEqual(ab["nelec_dict"], {'C': 4.0, 'Fe': 16.0, 'H': 1.0, 'O': 6.0, 'Sb': 15.0, 'Y': 11.0})
        self.assertEqual(ab["nkstot"], 50)
        self.assertEqual(ab["noccu_band"], 122)
        self.assertEqual(ab["normal_end"], True)
        self.assertAlmostEqual(ab["pressure"], -48.548971477133335, places=3)
        self.assertEqual(ab["pressures"], [-48.548971477133335])
        self.assertEqual(ab["relax_steps"], 1)
        self.assertEqual(ab["scf_steps"], 56)
        self.assertAlmostEqual(ab["scf_time"], 866.9699999999999, places=3)
        self.assertEqual(len(ab["scf_time_each_step"]), 56)
        self.assertAlmostEqual(ab["scf_time_each_step"][0], 16.28, places=3)
        self.assertAlmostEqual(ab["scf_time_each_step"][len(ab["scf_time_each_step"]) // 2], 15.28, places=3)
        self.assertAlmostEqual(ab["scf_time_each_step"][-1], 14.91, places=3)
        self.assertAlmostEqual(ab["step1_time"], 16.28, places=3)
        self.assertEqual(len(ab["stress"]), 9)
        self.assertAlmostEqual(ab["stress"][0], -52.8027809021, places=3)
        self.assertAlmostEqual(ab["stress"][len(ab["stress"]) // 2], -52.8786932479, places=3)
        self.assertAlmostEqual(ab["stress"][-1], -39.9654402814, places=3)
        self.assertAlmostEqual(ab["stress_time"], 46.2901955291803, places=3)
        self.assertEqual(len(ab["stresses"]), 1)
        self.assertEqual(len(ab["stresses"][0]), 9)
        stresses_flat = [item for sublist in ab["stresses"] for item in sublist]
        self.assertEqual(len(stresses_flat), 9)
        self.assertAlmostEqual(stresses_flat[0], -52.8027809021, places=3)
        self.assertAlmostEqual(stresses_flat[len(stresses_flat) // 2], -52.8786932479, places=3)
        self.assertAlmostEqual(stresses_flat[-1], -39.9654402814, places=3)
        self.assertAlmostEqual(ab["total_mag"], -0.000264531, places=3)
        self.assertEqual(ab["total_mags"], [-0.000264531])
        self.assertAlmostEqual(ab["total_time"], 932.927, places=3)
        self.assertEqual(ab["version"], "v3.10.1(f71921fe8 (Fri Nov 21 13:49:34 2025 +0800))")
        self.assertEqual(len(ab["virial"]), 9)
        self.assertAlmostEqual(ab["virial"][0], -57.594128079741175, places=3)
        self.assertAlmostEqual(ab["virial"][len(ab["virial"]) // 2], -57.67692874463316, places=3)
        self.assertAlmostEqual(ab["virial"][-1], -43.59192161863308, places=3)
        self.assertEqual(len(ab["virials"]), 1)
        self.assertEqual(len(ab["virials"][0]), 9)
        virials_flat = [item for sublist in ab["virials"] for item in sublist]
        self.assertEqual(len(virials_flat), 9)
        self.assertAlmostEqual(virials_flat[0], -57.594128079741175, places=3)
        self.assertAlmostEqual(virials_flat[len(virials_flat) // 2], -57.67692874463316, places=3)
        self.assertAlmostEqual(virials_flat[-1], -43.59192161863308, places=3)
        self.assertAlmostEqual(ab["volume"], 1747.56, places=3)

        # Metrics that returned None for this example:
        # atom_elec: None
        # atom_elec_mul: None
        # atom_elecs: None
        # atom_elecs_mul: None
        # atom_mag: None
        # atom_mag_mul: None
        # atom_mags: None
        # atom_mags_mul: None
        # atom_orb_elec: None
        # atom_orb_elec_mul: None
        # atom_orb_elecs: None
        # atom_orb_elecs_mul: None
        # atom_orb_mag: None
        # atom_orb_mags: None
        # denergy_womix: None
        # denergy_womix_last: None
        # dos: None
        # ds_mag: None
        # ds_mag_force: None
        # ds_mag_forces: None
        # ds_mags: None
        # ds_time: None
        # e_bandgap: None
        # force_time: None
        # k_coord: None
        # mem_psipw: None
        # mem_vkb: None
        # omp_num: None
        # pdos: None
        # point_group: None
        # point_group_in_space_group: None
        # relax_converge: None


if __name__ == '__main__':
    unittest.main()