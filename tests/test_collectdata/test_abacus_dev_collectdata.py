import unittest
from pathlib import Path

from abacustest.lib_collectdata.collectdata import RESULT


class TestAbacusCollectdata(unittest.TestCase):
    """Test Abacus result collection from SCF calculation"""

    def test_abacus_scf(self):
        """Test all extractable metrics from abacus-scf example"""
        ab = RESULT(path=Path(__file__).parent / "abacus-scf-dev/", fmt="abacus")

        input_dict = ab["INPUT"]
        self.assertIsInstance(input_dict, dict)
        self.assertEqual(input_dict["calculation"], "scf")
        self.assertEqual(input_dict["esolver_type"], "ksdft")
        self.assertEqual(input_dict["kspacing"], "0.14 0.14 0.14")
        self.assertAlmostEqual(ab["absolute_mag"], 2.49978e-06, places=3)
        self.assertEqual(len(ab["absolute_mags"]), 1)
        self.assertAlmostEqual(ab["absolute_mags"][0], 2.49978e-06, places=3)
        self.assertEqual(len(ab["atom_mag"]), 27)
        self.assertEqual(len(ab["atom_mag"][0]), 1)
        atom_mag_flat = [item for sublist in ab["atom_mag"] for item in sublist]
        self.assertEqual(len(atom_mag_flat), 27)
        self.assertAlmostEqual(atom_mag_flat[0], 6e-10, places=3)
        self.assertAlmostEqual(atom_mag_flat[len(atom_mag_flat) // 2], 5.4e-09, places=3)
        self.assertAlmostEqual(atom_mag_flat[-1], 8.43e-08, places=3)
        self.assertEqual(len(ab["atom_mags"]), 1)
        self.assertEqual(len(ab["atom_mags"][0]), 27)
        self.assertEqual(len(ab["atom_mags"][0][0]), 1)
        atom_mags_flat = [item for sublist in ab["atom_mags"] for subsublist in sublist for item in subsublist]
        self.assertEqual(len(atom_mags_flat), 27)
        self.assertAlmostEqual(atom_mags_flat[0], 6e-10, places=3)
        self.assertAlmostEqual(atom_mags_flat[len(atom_mags_flat) // 2], 5.4e-09, places=3)
        self.assertAlmostEqual(atom_mags_flat[-1], 8.43e-08, places=3)
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
        self.assertAlmostEqual(coordinate_flat[0], 3.953641216622, places=3)
        self.assertAlmostEqual(coordinate_flat[len(coordinate_flat) // 2], 3.952036857609, places=3)
        self.assertAlmostEqual(coordinate_flat[-1], 13.598405837985, places=3)
        self.assertEqual(len(ab["coordinate_init"]), 27)
        self.assertEqual(len(ab["coordinate_init"][0]), 3)
        coordinate_init_flat = [item for sublist in ab["coordinate_init"] for item in sublist]
        self.assertEqual(len(coordinate_init_flat), 81)
        self.assertAlmostEqual(coordinate_init_flat[0], 3.953641216622, places=3)
        self.assertAlmostEqual(coordinate_init_flat[len(coordinate_init_flat) // 2], 3.952036857609, places=3)
        self.assertAlmostEqual(coordinate_init_flat[-1], 13.598405837985, places=3)
        self.assertEqual(len(ab["coordinates"]), 1)
        self.assertEqual(len(ab["coordinates"][0]), 27)
        self.assertEqual(len(ab["coordinates"][0][0]), 3)
        coordinates_flat = [item for sublist in ab["coordinates"] for subsublist in sublist for item in subsublist]
        self.assertEqual(len(coordinates_flat), 81)
        self.assertAlmostEqual(coordinates_flat[0], 3.953641216622, places=3)
        self.assertAlmostEqual(coordinates_flat[len(coordinates_flat) // 2], 3.952036857609, places=3)
        self.assertAlmostEqual(coordinates_flat[-1], 13.598405837985, places=3)
        self.assertEqual(len(ab["denergy"]), 47)
        self.assertAlmostEqual(ab["denergy"][0], 0.0, places=3)
        self.assertAlmostEqual(ab["denergy"][len(ab["denergy"]) // 2], -0.000736085913, places=3)
        self.assertAlmostEqual(ab["denergy"][-1], -1.42422124e-07, places=3)
        self.assertAlmostEqual(ab["denergy_last"], -1.42422124e-07, places=3)
        self.assertEqual(len(ab["drho"]), 47)
        self.assertAlmostEqual(ab["drho"][0], 0.14941, places=3)
        self.assertAlmostEqual(ab["drho"][len(ab["drho"]) // 2], 0.00041724, places=3)
        self.assertAlmostEqual(ab["drho"][-1], 8.2814e-08, places=3)
        self.assertAlmostEqual(ab["drho_last"], 8.2814e-08, places=3)
        self.assertEqual(ab["ds_lambda_rms"], [])
        self.assertEqual(ab["ds_lambda_step"], [])
        self.assertEqual(len(ab["ds_mag"]), 27)
        self.assertEqual(len(ab["ds_mag"][0]), 1)
        ds_mag_flat = [item for sublist in ab["ds_mag"] for item in sublist]
        self.assertEqual(len(ds_mag_flat), 27)
        self.assertAlmostEqual(ds_mag_flat[0], 6e-10, places=3)
        self.assertAlmostEqual(ds_mag_flat[len(ds_mag_flat) // 2], 5.4e-09, places=3)
        self.assertAlmostEqual(ds_mag_flat[-1], 8.43e-08, places=3)
        self.assertEqual(len(ab["ds_mags"]), 1)
        self.assertEqual(len(ab["ds_mags"][0]), 27)
        self.assertEqual(len(ab["ds_mags"][0][0]), 1)
        ds_mags_flat = [item for sublist in ab["ds_mags"] for subsublist in sublist for item in subsublist]
        self.assertEqual(len(ds_mags_flat), 27)
        self.assertAlmostEqual(ds_mags_flat[0], 6e-10, places=3)
        self.assertAlmostEqual(ds_mags_flat[len(ds_mags_flat) // 2], 5.4e-09, places=3)
        self.assertAlmostEqual(ds_mags_flat[-1], 8.43e-08, places=3)
        self.assertAlmostEqual(ab["efermi"], 1.6247208577, places=3)
        self.assertEqual(len(ab["energies"]), 1)
        self.assertAlmostEqual(ab["energies"][0], -28316.628295020262, places=3)
        self.assertAlmostEqual(ab["energy"], -28316.628295020262, places=3)
        self.assertAlmostEqual(ab["energy_ks"], -28316.6282950203, places=3)
        self.assertAlmostEqual(ab["energy_per_atom"], -1048.7640109266763, places=3)
        self.assertEqual(ab["ibzk"], 10)
        self.assertEqual(ab["kpt"], [5, 5, 2])
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
        self.assertEqual(ab["ncore"], 16)
        self.assertAlmostEqual(ab["nelec"], 244.0, places=3)
        self.assertEqual(ab["nkstot"], 16)
        self.assertEqual(ab["normal_end"], True)
        self.assertEqual(ab["omp_num"], 1)
        self.assertEqual(ab["scf_steps"], 47)
        self.assertAlmostEqual(ab["scf_time"], 631.93, places=3)
        self.assertEqual(len(ab["scf_time_each_step"]), 47)
        self.assertAlmostEqual(ab["scf_time_each_step"][0], 14.13, places=3)
        self.assertAlmostEqual(ab["scf_time_each_step"][len(ab["scf_time_each_step"]) // 2], 13.32, places=3)
        self.assertAlmostEqual(ab["scf_time_each_step"][-1], 15.75, places=3)
        self.assertAlmostEqual(ab["step1_time"], 14.13, places=3)
        self.assertAlmostEqual(ab["total_mag"], -5.88021e-07, places=3)
        self.assertEqual(len(ab["total_mags"]), 1)
        self.assertAlmostEqual(ab["total_mags"][0], -5.88021e-07, places=3)
        self.assertAlmostEqual(ab["total_time"], 645.799, places=3)
        self.assertEqual(ab["version"], "v3.9.0.24(14746d392 (Mon Feb 2 14:11:41 2026 +0800))")
        self.assertAlmostEqual(ab["volume"], 1747.558574478872, places=3)

        # Metrics that returned None for this example:
        # atom_elec: None
        # atom_elec_mul: None
        # atom_elecs: None
        # atom_elecs_mul: None
        # atom_mag_mul: None
        # atom_mags_mul: None
        # atom_orb_elec: None
        # atom_orb_elec_mul: None
        # atom_orb_elecs: None
        # atom_orb_elecs_mul: None
        # atom_orb_mag: None
        # atom_orb_mags: None
        # atomlabel_list: None
        # band: None
        # band_gap: None
        # band_plot: None
        # band_weight: None
        # denergy_womix: None
        # denergy_womix_last: None
        # dos: None
        # ds_mag_force: None
        # ds_mag_forces: None
        # ds_time: None
        # e_bandgap: None
        # element: None
        # element_list: None
        # fft_grid: None
        # force: None
        # force_time: None
        # forces: None
        # k_coord: None
        # label: None
        # largest_gradient: None
        # largest_gradient_stress: None
        # mem_psipw: None
        # mem_vkb: None
        # nbase: None
        # nelec_dict: None
        # noccu_band: None
        # pdos: None
        # point_group: None
        # point_group_in_space_group: None
        # pressure: None
        # pressures: None
        # relax_converge: None
        # relax_steps: None
        # stress: None
        # stress_time: None
        # stresses: None
        # virial: None
        # virials: None


if __name__ == '__main__':
    unittest.main()