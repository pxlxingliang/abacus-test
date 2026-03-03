import os
from collections import defaultdict

from typing import Dict, List, Tuple, Optional, Union, Any
from pathlib import Path

import numpy as np

from abacustest import ReadInput, ReadKpt, RESULT, AbacusSTRU
from abacustest.lib_prepare.comm import Cartesian2Direct, real2rec


class BandData:
    """Class for managing band data from ABACUS band calculations."""

    def __init__(
        self,
        high_symm_labels: Dict[str, List[float]],
        kpaths: List[Dict[str, Any]],
        kpath_cum_dist: Union[np.ndarray, List[float]],
        efermi: Optional[float] = None,
        band_data: Optional[Union[np.ndarray, List[np.ndarray]]] = None,
    ):
        """
        Initialize BandData object.
        Args:
            high_symm_labels (Dict[str, List[float]]): High symmetry labels and coordinates. For example: {"GAMMA": [0.0, 0.0, 0.0], "X": [0.5, 0.0, 0.0], "L": [0.5, 0.5, 0.0]}.
            kpaths (List[Dict[str, Union[float, str]]]): K-path information. For example: [{"start": "GAMMA", "end": "X", "start_nkpt": 0, "end_nkpt": 30}, {"start": "X", "end": "U", "start_nkpt": 30, "end_nkpt": 60}].
                For each kpath, band data are stored in the range [start_nkpt, end_nkpt] of `band_data` for each spin channel.
            kpath_cum_dist (Union[np.ndarray, List[float]]): K-path cumulative distance. For example: [0.0, 0.1, 0.2, 0.3, 0.8, 0.85, 0.9]. Discontinuity in k-path distance will be automatically removed.
            efermi (float): Fermi energy used during the initialization of the BandData object. Will subtract this value from the band data. If None, the band data will not be modified.
            band_data (np.ndarray): Energy level data for each k-point. For nspin=1 or 4, it is an np.ndarray shaped (1, nkpts, nbands). For nspin=2, it is shaped (2, nkpts, nbands).
        """
        self.high_symm_labels = high_symm_labels
        self.kpaths = kpaths
        self.kpath_lengths = self._remove_jump_dist(kpath_cum_dist, kpaths)

        self.efermi = efermi

        if (
            isinstance(band_data, np.ndarray)
            and len(band_data.shape) == 3
            and band_data.shape[0] in (1, 2)
        ):
            if efermi is not None:
                self.band_data = band_data - efermi
            else:
                print("efermi is None, band data will not be shifted")
                self.band_data = band_data
        else:
            raise ValueError(
                "band_data must be an np.ndarray with shape (1, nkpts, nbands) or (2, nkpts, nbands)"
            )

        self.nspin, self.nkpts, self.nbands = self.band_data.shape
        self._build_label_mapping()
        self._build_kpath_segments()

    def _get_spin_label(self, spin_index):
        """Get human-readable label for spin index.

        Args:
            spin_index (int): Spin index (0, 1, ...).

        Returns:
            str: Label for the spin channel. For nspin=2: "spin_up" for index 0,
                 "spin_down" for index 1. For other nspin values, returns str(spin_index).
        """
        if self.nspin == 2:
            return "spin_up" if spin_index == 0 else "spin_down"
        return str(spin_index)

    def _build_kpath_segments(self):
        """
        Build mapping from kpoint index to segment for coordinate interpolation.
        """
        self.kpoint_to_segment = {}
        for kpath in self.kpaths:
            start_idx, end_idx = kpath["start_nkpt"], kpath["end_nkpt"]
            for i in range(start_idx, end_idx + 1):
                self.kpoint_to_segment[i] = kpath

    def _get_kpoint_coord(self, idx):
        """
        Get direct coordinate of kpoint at index idx.

        Args:
            idx: kpoint index

        Returns:
            List[float]: [kx, ky, kz] in direct coordinates
        """
        if idx not in self.kpoint_to_segment:
            # If index not in any segment (shouldn't happen), return None
            return None

        segment = self.kpoint_to_segment[idx]
        start_label, end_label = segment["start"], segment["start"]
        start_coord = self.high_symm_labels.get(start_label)
        end_coord = self.high_symm_labels.get(end_label)

        if start_coord is None or end_coord is None:
            return None

        start_idx, end_idx = segment["start_nkpt"], segment["end_nkpt"]

        if start_idx == end_idx:
            return start_coord  # Single point segment

        # Linear interpolation
        t = (idx - start_idx) / (end_idx - start_idx)
        coord = []
        for i in range(3):
            coord.append(start_coord[i] + t * (end_coord[i] - start_coord[i]))
        return coord

    def _build_label_mapping(self):
        """
        Build mapping from label to kpoint indices and list of labels per kpoint.
        """
        self.label_to_indices = defaultdict(list)
        self.kpoint_labels = [None] * self.nkpts
        for kpath in self.kpaths:
            start_label, start_idx = kpath["start"], kpath["start_nkpt"]
            end_label, end_idx = kpath["end"], kpath["end_nkpt"]
            self.kpoint_labels[start_idx] = start_label
            self.kpoint_labels[end_idx] = end_label
            self.label_to_indices[start_label].append(start_idx)
            self.label_to_indices[end_label].append(end_idx)

    @staticmethod
    def _remove_jump_dist(
        kpath_cum_dist: Union[np.ndarray, List[float]],
        kpaths: List[Dict[str, Any]],
    ) -> np.ndarray:
        """
        Remove the jump in the k-path distance.
        Args:
            kpath_cum_dist (Union[np.ndarray, List[float]]): K-path cumulative distance. For example: [0.0, 0.1, 0.2, 0.3, 0.8, 0.85, 0.9].
            kpaths (List[Dict[str, Union[int, str]]]): K-path information. For example: [{"start": "GAMMA", "end": "X", "start_nkpt": 0, "end_nkpt": 30}, {"start": "X", "end": "U", "start_nkpt": 30, "end_nkpt": 60}]
        Returns:
            np.ndarray: K-path cumulative distance without jump
        """
        # Convert to numpy array if needed
        if not isinstance(kpath_cum_dist, np.ndarray):
            kpath_cum_dist = np.array(kpath_cum_dist)

        for ipath in range(len(kpaths) - 1):
            # Check if there's a jump in cumulative distance between segments
            if kpaths[ipath]["end"] != kpaths[ipath + 1]["start"]:
                start_nkpt_next = kpaths[ipath + 1]["start_nkpt"]
                end_nkpt_curr = kpaths[ipath]["end_nkpt"]
                jump_dist = (
                    kpath_cum_dist[start_nkpt_next] - kpath_cum_dist[end_nkpt_curr]
                )
                kpath_cum_dist[start_nkpt_next:] -= jump_dist

        return kpath_cum_dist

    @staticmethod
    def ReadFromAbacusJob(
        abacusjob_dir: str,
        efermi: Optional[float] = None,
        high_symm_labels: Optional[List[str]] = None,
    ):
        """
        Read band data from output directory of a finished ABACUS band calculation.

        Args:
            abacusjob_dir (str): Path to the output directory of a finished ABACUS band calculation.
            efermi (float): Fermi energy used during the initialization of the BandData object. Will subtract this value from the band data.
            high_symm_labels (List[str]): List of high symmetry labels. If None, the labels will be read from the KPT file, using comments after "#" at each line.

        Returns:
            BandData: A BandData object containing the band data and high symmetry labels.
        """
        input_params = ReadInput(os.path.join(abacusjob_dir, "INPUT"))
        suffix = input_params.get("suffix", "ABACUS")
        nspin = input_params.get("nspin", 1)
        stru_file = os.path.join(abacusjob_dir, input_params.get("stru_file", "STRU"))
        stru = AbacusSTRU.read(stru_file)
        rec_cell = real2rec(stru.cell)

        kpt_result = ReadKpt(abacusjob_dir)
        if kpt_result is None:
            raise ValueError(f"Failed to read KPT file from {abacusjob_dir}")
        kpt_data, model = kpt_result
        if model not in ["line", "line_cartesian"]:
            raise ValueError(
                f"KPT file must be in 'line' or 'line_cartesian' mode for band calculation, got '{model}'"
            )

        # Process kpt_data to extract kpaths and high symmetry points
        high_symm_kpts = {}
        insert_kpt_nums = []
        high_symm_labels_list = []

        # Ensure kpt_data is a list
        if not isinstance(kpt_data, list):
            raise TypeError(f"Expected kpt_data to be a list, got {type(kpt_data)}")

        for i, kpt in enumerate(kpt_data):
            # kpt format: [kx, ky, kz, npoints, label] or [kx, ky, kz, npoints]
            if not isinstance(kpt, (list, tuple)):
                raise TypeError(
                    f"Expected kpt to be a list or tuple, got {type(kpt)} at index {i}"
                )

            if len(kpt) < 4:
                raise ValueError(
                    f"Must have at least 4 elements in line-mode KPT file, got {len(kpt)} in line {i}"
                )

            kx, ky, kz, npoints = kpt[:4]
            if model == "line":
                kpt_coord_direct = [float(kx), float(ky), float(kz)]
            else:
                kpt_coord_cart = [float(kx), float(ky), float(kz)]
                kpt_coord_direct = Cartesian2Direct([kpt_coord_cart], rec_cell)[0]

            # Get label
            if len(kpt) >= 5:
                # Use label is provided in KPT file
                label = str(kpt[4]).strip()
                if label.startswith("#"):
                    label = label[1:]
            elif high_symm_labels is not None and i < len(high_symm_labels):
                label = high_symm_labels[i]
            else:
                label = f"K{i}"  # Generate default label

            insert_kpt_nums.append(int(npoints))
            high_symm_labels_list.append(label)
            high_symm_kpts[label] = kpt_coord_direct

        # Build kpaths
        start_nkpt, kpaths = 0, []
        for i in range(len(high_symm_labels_list)):
            # Loop over all labels and jump over the discontinuous points in the k-path
            if insert_kpt_nums[i] > 1:
                end_nkpt = start_nkpt + insert_kpt_nums[i] - 1

                # If the next segment has only 1 k-point, include it in this segment.
                # Used to treat discountinuous k-paths and end points
                if i < len(high_symm_labels_list) - 2 and insert_kpt_nums[i + 1] == 1:
                    end_nkpt += 1

                kpaths.append(
                    {
                        "start": high_symm_labels_list[i],
                        "end": high_symm_labels_list[i + 1],
                        "start_nkpt": start_nkpt,
                        "end_nkpt": end_nkpt,
                    }
                )
                start_nkpt = end_nkpt + 1

        abacusresult = RESULT(fmt="abacus", path=abacusjob_dir)
        if abacusresult is None:
            raise ValueError("Failed to read ABACUS results")
        if efermi is None:
            try:
                efermi = abacusresult["efermi"]
            except (KeyError, TypeError):
                raise ValueError(
                    "Fermi energy (efermi) not found in ABACUS results and not provided"
                )

        # Read band data for nspin=1 or 4
        band_file = os.path.join(abacusjob_dir, f"OUT.{suffix}/BANDS_1.dat")
        with open(band_file, "r") as f:
            original_band_data = np.loadtxt(f)
            kpath_cum_dist, band = original_band_data[:, 1], original_band_data[:, 2:]
            band_data = np.expand_dims(
                band, axis=0
            )  # Transform to ndarray shaped (1, nkpt, nband)
        # Read band data for nspin=2
        if nspin == 2:
            band_file_dw = os.path.join(abacusjob_dir, f"OUT.{suffix}/BANDS_2.dat")
            with open(band_file_dw, "r") as f:
                original_band_data = np.loadtxt(f)
                band_dw = original_band_data[:, 2:]
                band_data = np.stack(
                    (band, band_dw), axis=0
                )  # Transform to ndarray shaped (2, nkpt, nband)

        return BandData(high_symm_kpts, kpaths, kpath_cum_dist, efermi, band_data)

    def is_metal(self, efermi: Optional[float] = None) -> bool:
        """
        Check if the band structure is metallic.
        Returns:
            bool: True if the band structure is metallic, False otherwise.
        """
        if efermi is not None:
            band_data = self.band_data - efermi
        elif self.efermi is not None:
            band_data = self.band_data  # Already shifted
        else:
            raise RuntimeError(
                "Fermi energy is not given or not present in the band data"
            )

        nspin_channel, nkpt, nband = band_data.shape
        for ispin_channel in range(nspin_channel):
            for ikpt in range(nkpt - 1):
                for iband in range(nband):
                    # Check whether the band crosses fermi energy. Also work if band is discontinous
                    if (
                        band_data[ispin_channel, ikpt, iband]
                        * band_data[ispin_channel, ikpt + 1, iband]
                        < 0
                    ):
                        return True

        return False

    def get_band_branch_data(self, branch: List[str], extra_end_point: bool = False):
        """
        Get band data of the selected branch.
        Args:
            branch: A list containing 2 high-symmetry labels to define the band. e.g.: ['G', 'H']
            extra_end_point: Whether to include extra band data of the first kpoint of the next branch.
                Used in getting band data for plotting bands to avoid discountinuites at continous high-symmetry points.
        Returns:
            branch_band_data: band data of the assigned band branch. If the band is reversed, the band data will be reversed.
            shifted_kpath_coord (np.ndarray): A 1D np.ndarray containing distance of each kpoint to starting point of the branch.
        """
        if isinstance(branch, List):
            assert len(branch) == 2
            for i in branch:
                assert isinstance(i, str)

        start, end = branch
        for ikpath, kpath in enumerate(self.kpaths):
            if (
                kpath["start"] == start and kpath["end"] == end
            ):  # Find the matched branch
                start_nkpt, end_nkpt = kpath["start_nkpt"], kpath["end_nkpt"]
                if extra_end_point:
                    if ikpath < len(self.kpaths) - 1:  # Not the last band branch
                        if (
                            kpath["end"] == self.kpaths[ikpath + 1]["start"]
                        ):  # Continus with the next branch
                            end_nkpt += 1
                        else:
                            start_nkpt, end_nkpt = (
                                kpath["start_nkpt"],
                                kpath["end_nkpt"],
                            )

                branch_band_data = self.band_data[:, start_nkpt : end_nkpt + 1, :]
                kpath_coord = self.kpath_lengths[start_nkpt : end_nkpt + 1]
                shifted_kpath_coord = kpath_coord - min(kpath_coord)

                return branch_band_data, shifted_kpath_coord
            elif kpath["start"] == end and kpath["end"] == start:
                # Revert the reqested branch to get data, and reverse the band data
                reverse_branch_band_data, shifted_kpath_coord = (
                    self.get_band_branch_data([end, start], extra_end_point)
                )
                if reverse_branch_band_data is None or shifted_kpath_coord is None:
                    return None, None
                return reverse_branch_band_data[:, ::-1, :], shifted_kpath_coord

        return None, None  # If no matched branch found

    def _process_indices(self, spin_indices, energy=None):
        """Convert list of (spin, k, b) tuples to result dict.
        If energy is provided, use it; otherwise compute from first index."""
        if not spin_indices:
            return {
                "band_index": {},
                "kpoint_index": [],
                "kpoint_labels": [],
                "kpoint_coord": [],
                "energy": None,
            }

        # Group by spin and band index
        band_index = defaultdict(set)
        kpoint_indices = set()
        for spin, k, b in spin_indices:
            band_index[spin].add(b)
            kpoint_indices.add(k)

        # Convert sets to sorted lists
        band_index_dict = {spin: sorted(bands) for spin, bands in band_index.items()}

        kpoint_indices = sorted(kpoint_indices)
        # Get labels for each kpoint
        kpoint_labels = [self.kpoint_labels[idx] for idx in kpoint_indices]

        kpoint_coord = [self._get_kpoint_coord(idx) for idx in kpoint_indices]

        # Use provided energy if available, otherwise compute from first index
        if energy is None:
            if spin_indices:
                # Use first index to get energy
                spin_idx, k_idx, b_idx = spin_indices[0]
                energy = self.band_data[spin_idx, k_idx, b_idx]
            else:
                energy = None

        return {
            "band_index": band_index_dict,
            "kpoint_index": kpoint_indices,
            "kpoint_labels": kpoint_labels,
            "kpoint_coord": kpoint_coord,
            "energy": energy,
        }

    def _find_band_edge(self, below_fermi=True, tol=1e-4):
        """
        Find band edge (VBM if below_fermi=True, CBM if below_fermi=False).
        Returns:
            edge_energy_global (float): global edge energy (-inf/+inf if not found)
            edge_indices_global (list): list of (spin, k, b) tuples for global edge
            edge_energy_per_spin (list): per-spin edge energies
            edge_indices_per_spin (list): per-spin lists of (spin, k, b) tuples
        """
        if below_fermi:
            # VBM: look for bands below Fermi level, maximize energy
            energy_init = -float("inf")
            compare = lambda val, best: val > best
            condition = lambda val: val < -tol
        else:
            # CBM: look for bands above Fermi level, minimize energy
            energy_init = float("inf")
            compare = lambda val, best: val < best
            condition = lambda val: val >= -tol

        edge_energy_global = energy_init
        edge_indices_global = []
        edge_energy_per_spin = [energy_init] * self.nspin
        edge_indices_per_spin = [[] for _ in range(self.nspin)]

        for spin in range(self.nspin):
            for k in range(self.nkpts):
                for b in range(self.nbands):
                    val = self.band_data[spin, k, b]
                    if condition(val):
                        # Update global edge
                        if abs(val - edge_energy_global) < tol:
                            edge_indices_global.append((spin, k, b))
                        elif compare(val, edge_energy_global):
                            edge_energy_global = val
                            edge_indices_global = [(spin, k, b)]

                        # Update per-spin edge
                        if abs(val - edge_energy_per_spin[spin]) < tol:
                            edge_indices_per_spin[spin].append((spin, k, b))
                        elif compare(val, edge_energy_per_spin[spin]):
                            edge_energy_per_spin[spin] = val
                            edge_indices_per_spin[spin] = [(spin, k, b)]

        return (
            edge_energy_global,
            edge_indices_global,
            edge_energy_per_spin,
            edge_indices_per_spin,
        )

    def _get_edge(self, below_fermi, tol=1e-4, spin_resolved=False):
        """
        Internal helper to get VBM or CBM.
        """
        null_result = {
            "band_index": {},
            "kpoint_index": [],
            "kpoint_labels": [],
            "kpoint_coord": [],
            "energy": None,
        }
        if self.is_metal():
            if spin_resolved:
                return {
                    self._get_spin_label(spin): null_result
                    for spin in range(self.nspin)
                }
            else:
                return null_result

        (
            edge_energy_global,
            edge_indices_global,
            edge_energy_per_spin,
            edge_indices_per_spin,
        ) = self._find_band_edge(below_fermi=below_fermi, tol=tol)

        # Check if no edge found
        if below_fermi:
            no_edge = edge_energy_global == -float("inf")
        else:
            no_edge = edge_energy_global == float("inf")
        if no_edge:
            if spin_resolved:
                return {
                    self._get_spin_label(spin): null_result
                    for spin in range(self.nspin)
                }
            else:
                return null_result

        # Process results
        if spin_resolved:
            result = {}
            for spin in range(self.nspin):
                label = self._get_spin_label(spin)
                result[label] = self._process_indices(
                    edge_indices_per_spin[spin], energy=edge_energy_per_spin[spin]
                )
            result["global"] = self._process_indices(
                edge_indices_global, energy=edge_energy_global
            )
            return result
        else:
            return self._process_indices(edge_indices_global, energy=edge_energy_global)

    def get_vbm(self, tol=1e-4, spin_resolved=False):
        """
        Get data about the valence band maximum (VBM).

        Args:
            tol (float): Tolerance for energy comparison (eV).
            spin_resolved (bool): If True, return results for each spin channel separately.

        Returns:
            If spin_resolved=False:
                dict with keys "band_index", "kpoint_index", "kpoint_labels", "kpoint_coord", "energy":
                    - "band_index": dict mapping spin index (int) to list of band indices.
                    - "kpoint_index": list of indices of kpoints where VBM occurs.
                    - "kpoint_labels": list of labels (str or None) for each kpoint in kpoint_index.
                    - "kpoint_coord": list of coordinates [kx, ky, kz] for each kpoint in kpoint_index.
                    - "energy": energy of VBM relative to Fermi level (eV).
                      If Fermi level is not set (efermi=None), returns absolute energy.

            If spin_resolved=True:
                dict with keys for each spin label (e.g., "spin_up", "spin_down" for nspin=2, otherwise str(spin_index)) and "global":
                    - spin_label (str): dict with same keys as above for that spin channel.
                    - "global": dict with same keys as above for the global VBM (across all spins).
        """
        return self._get_edge(below_fermi=True, tol=tol, spin_resolved=spin_resolved)

    def get_cbm(self, tol=1e-4, spin_resolved=False):
        """
        Get data about the conduction band minimum (CBM).

        Args:
            tol (float): Tolerance for energy comparison (eV).
            spin_resolved (bool): If True, return results for each spin channel separately.

        Returns:
            If spin_resolved=False:
                dict with keys "band_index", "kpoint_index", "kpoint_labels", "kpoint_coord", "energy":
                    - "band_index": dict mapping spin index (int) to list of band indices.
                    - "kpoint_index": list of indices of kpoints where CBM occurs.
                    - "kpoint_labels": list of labels (str or None) for each kpoint in kpoint_index.
                    - "kpoint_coord": list of coordinates [kx, ky, kz] for each kpoint in kpoint_index.
                    - "energy": energy of CBM relative to Fermi level (eV).
                      If Fermi level is not set (efermi=None), returns absolute energy.

            If spin_resolved=True:
                dict with keys for each spin label (e.g., "spin_up", "spin_down" for nspin=2, otherwise str(spin_index)) and "global":
                    - spin_label (str): dict with same keys as above for that spin channel.
                    - "global": dict with same keys as above for the global CBM (across all spins).
        """
        return self._get_edge(below_fermi=False, tol=tol, spin_resolved=spin_resolved)

    def _calculate_gap_from_results(self, vbm_dict, cbm_dict):
        """Calculate band gap from VBM and CBM result dictionaries.

        Args:
            vbm_dict: Dictionary from get_vbm()
            cbm_dict: Dictionary from get_cbm()

        Returns:
            float or None: Band gap energy (CBM - VBM) or None if not available.
        """
        if not (isinstance(vbm_dict, dict) and isinstance(cbm_dict, dict)):
            return None

        vbm_energy = vbm_dict.get("energy")
        cbm_energy = cbm_dict.get("energy")
        if vbm_energy is not None and cbm_energy is not None:
            return cbm_energy - vbm_energy
        return None

    def get_band_gap(self, tol=1e-4, spin_resolved=False):
        """
        Get the band gap energy.

        Args:
            tol (float): Tolerance for energy comparison (eV).
            spin_resolved (bool): If True, return results for each spin channel separately.

        Returns:
            If spin_resolved=False:
                float or None: Band gap energy (CBM energy - VBM energy) in eV.
                Returns None if the system is metallic (no band gap) or if VBM/CBM not found.

            If spin_resolved=True:
                dict with keys for each spin label (e.g., "spin_up", "spin_down" for nspin=2, otherwise str(spin_index)) and "global":
                    - spin_label (str): Band gap for that spin channel (CBM_energy[spin] - VBM_energy[spin]).
                    - "global": Global band gap (global CBM energy - global VBM energy).
                Returns None values for metallic systems or when VBM/CBM not found.
        """
        # Get VBM and CBM data
        vbm_result = self.get_vbm(tol=tol, spin_resolved=spin_resolved)
        cbm_result = self.get_cbm(tol=tol, spin_resolved=spin_resolved)

        if not spin_resolved:
            # Non-spin-resolved case: both results are dictionaries with "energy" key
            return self._calculate_gap_from_results(vbm_result, cbm_result)

        # Spin-resolved case: results are dictionaries with spin labels and "global" key
        result = {}

        # Handle global band gap
        if "global" in vbm_result and "global" in cbm_result:
            result["global"] = self._calculate_gap_from_results(
                vbm_result["global"], cbm_result["global"]
            )
        else:
            result["global"] = None

        # Handle per-spin band gaps
        for spin in range(self.nspin):
            label = self._get_spin_label(spin)
            if label in vbm_result and label in cbm_result:
                result[label] = self._calculate_gap_from_results(
                    vbm_result[label], cbm_result[label]
                )
            else:
                result[label] = None

        return result

    def plot_band(
        self,
        emin: float = -10,
        emax: float = 10,
        fig_name: str = "band.png",
        plot_kpaths: Optional[List[List[str]]] = None,
    ):
        """
        Plot the band structure.
        Args:
            emin (float): Minimum energy in the band plot.
            emax (float): Maximum energy in the band plot.
            fig_name (str): Name of the figure to save.
            plot_kpaths: List of k-point labels to plot. If provided, should be in the form like [["G", "H"], ["H", "K"]].
                If not provided, will use kpath stored in the object.
        """
        if plot_kpaths is None:
            plot_kpaths = []
            for kpath in self.kpaths:
                plot_kpaths.append([kpath["start"], kpath["end"]])

        # Prepare data for plot_multiple_bands
        branches, shifted_kpath_coord, band_data = [], [], []

        for plot_kpath in plot_kpaths:
            branch_band_data, shifted_coord = self.get_band_branch_data(
                plot_kpath, extra_end_point=True
            )
            if branch_band_data is None:
                raise ValueError("Provided kpath not found in the band structure")

            branches.append(plot_kpath)
            shifted_kpath_coord.append(shifted_coord)
            band_data.append(branch_band_data)

        plot_band_datas = [
            {
                "branches": branches,
                "shifted_kpath_coord": shifted_kpath_coord,
                "band_data": band_data,
            }
        ]

        # Call plot_multiple_bands
        plot_multiple_bands(
            plot_band_datas=plot_band_datas,
            mat_names=None,
            fig_name=fig_name,
            emin=emin,
            emax=emax,
        )

    def get_effective_mass(
        self,
        edge_type="CBM",
        direction_labels=None,
        num_fit_points=5,
        bilateral_fit=False,
        plot_fit=False,
    ):
        """
        Calculate effective mass for CBM or VBM along a specified direction.

        Args:
            edge_type (str): "CBM" or "VBM".
            direction_labels (list): List of two high-symmetry labels, e.g., ['G', 'X'].
                The first label should be the high-symmetry point where the CBM/VBM is located.
            num_fit_points (int): Number of k-points to use for quadratic fitting.
            bilateral_fit (bool): If True, use data from both sides of the extremum.
                Requires the reverse direction branch to exist.
            plot_fit (bool): If True, plot the fitting curve and data points.

        Note:
            This method does not support nspin=2 (spin-polarized calculations).

        Returns:
            list: List of effective mass results for each degenerate band at the edge.
                  Each element is a dict containing:
                - effective_mass (float): effective mass in units of electron mass.
                  For CBM (electron), effective mass is positive.
                  For VBM (hole), effective mass is negative (positive if taking absolute value).
                - curvature (float): second derivative d²E/dk² in eV/Å².
                  Positive for CBM (upward curvature), negative for VBM (downward curvature).
                - fit_coeffs (list): quadratic fit coefficients [a, b, c] for E = a*(k - k0)^2 + b.
                - band_index (int): band index used.
        """
        import numpy as np
        import matplotlib.pyplot as plt

        if direction_labels is None or len(direction_labels) != 2:
            raise ValueError(
                "direction_labels must be a list of two high-symmetry labels"
            )

        start_label, end_label = direction_labels

        if self.nspin == 2:
            raise ValueError("Effective mass calculation not supported for nspin=2")

        # Get CBM/VBM data
        if edge_type.upper() == "CBM":
            edge_result = self.get_cbm(spin_resolved=False)
        elif edge_type.upper() == "VBM":
            edge_result = self.get_vbm(spin_resolved=False)
        else:
            raise ValueError("edge_type must be 'CBM' or 'VBM'")

        if edge_result is None or edge_result["energy"] is None:
            raise ValueError(f"No {edge_type} found in band data")

        if not edge_result["kpoint_index"]:
            raise ValueError(f"No k-point indices for {edge_type}")

        # Get band indices (for nspin=1, band_index dict maps spin 0 to list of bands)
        # There may be multiple degenerate bands at CBM/VBM
        band_indices = list(edge_result["band_index"].values())[0]
        if not band_indices:
            raise ValueError(f"No band indices found for {edge_type}")

        # Get branch data for forward direction
        forward_branch = [start_label, end_label]
        forward_data, forward_coord = self.get_band_branch_data(
            forward_branch, extra_end_point=False
        )
        if forward_data is None:
            raise ValueError(f"Branch {forward_branch} not found in band structure")
        assert forward_coord is not None
        assert self.nspin == 1

        # Find which kpath segment corresponds to the forward branch
        forward_segment = None
        for kpath in self.kpaths:
            if kpath["start"] == start_label and kpath["end"] == end_label:
                forward_segment = kpath
                break
        if forward_segment is None:
            # Check reverse direction
            for kpath in self.kpaths:
                if kpath["start"] == end_label and kpath["end"] == start_label:
                    forward_segment = kpath
                    break
        if forward_segment is None:
            raise ValueError(f"No segment found for branch {forward_branch}")

        # Select a k-point index that lies on this segment
        global_k_idx = None
        for k_idx in edge_result["kpoint_index"]:
            if forward_segment["start_nkpt"] <= k_idx <= forward_segment["end_nkpt"]:
                global_k_idx = k_idx
                break
        if global_k_idx is None:
            # None of the extremum k-points are on the requested branch
            raise ValueError(f"No {edge_type} k-point found on branch {forward_branch}")

        # Determine local index within branch
        # forward_coord is shifted_kpath_coord for the branch, starting at 0 at start_label
        # If extremum is at start_label, local index is 0.
        # If extremum is at end_label (reverse direction), local index is last.
        if forward_segment["start"] == start_label:
            local_idx = global_k_idx - forward_segment["start_nkpt"]
        else:
            # segment is reversed relative to forward_branch
            local_idx = forward_segment["end_nkpt"] - global_k_idx

        # Ensure local_idx is within forward_data shape
        nkpt_branch = forward_data.shape[1]
        if local_idx < 0 or local_idx >= nkpt_branch:
            raise ValueError("Extremum k-point index mismatch with branch data")

        # Collect data for bilateral fitting
        reverse_data = None
        reverse_coord = None
        rev_local_idx = None
        if bilateral_fit:
            # Try to get reverse branch data
            reverse_branch = [end_label, start_label]
            reverse_data, reverse_coord = self.get_band_branch_data(
                reverse_branch, extra_end_point=False
            )
            if reverse_data is None:
                # Reverse branch not found, fallback to unilateral
                bilateral_fit = False
                print("Warning: Reverse branch not found, using unilateral fitting")
            else:
                # Determine local index in reverse branch (extremum should be at start_label of reverse branch)
                # Find segment for reverse branch
                rev_segment = None
                for kpath in self.kpaths:
                    if kpath["start"] == end_label and kpath["end"] == start_label:
                        rev_segment = kpath
                        break
                if rev_segment is None:
                    bilateral_fit = False
                    print(
                        "Warning: Cannot locate reverse segment, using unilateral fitting"
                    )
                else:
                    # Extremum should be at start of reverse branch (which is end_label)
                    rev_local_idx = global_k_idx - rev_segment["start_nkpt"]
                    # Ensure rev_local_idx is within reverse_data shape
                    nkpt_rev = reverse_data.shape[1]
                    if rev_local_idx < 0 or rev_local_idx >= nkpt_rev:
                        bilateral_fit = False
                        print(
                            "Warning: Extremum index mismatch in reverse branch, using unilateral fitting"
                        )

        # Calculate effective mass for each degenerate band
        # Forward k distances relative to extremum (in reduced coordinates, i.e., k/(2π))
        # This is the same for all bands at the same k-point
        forward_k_reduced = forward_coord - forward_coord[local_idx]

        results = []
        for band_idx in band_indices:
            # Forward branch energies for the band (nspin=1, use spin index 0)
            forward_energy = forward_data[0, :, band_idx]

            if bilateral_fit and reverse_data is not None:
                assert (
                    reverse_data is not None
                    and reverse_coord is not None
                    and rev_local_idx is not None
                )
                reverse_energy = reverse_data[0, :, band_idx]
                # Reverse k distances: negative relative to extremum (reduced coordinates)
                reverse_k_reduced = -(reverse_coord - reverse_coord[rev_local_idx])
                # Combine forward and reverse points
                combined_k_reduced = np.concatenate(
                    [reverse_k_reduced, forward_k_reduced]
                )
                combined_energy = np.concatenate([reverse_energy, forward_energy])
                # Select num_fit_points closest to extremum
                sorted_indices = np.argsort(np.abs(combined_k_reduced))
                selected_indices = sorted_indices[:num_fit_points]
                k_fit_reduced = combined_k_reduced[selected_indices]
                energy_fit = combined_energy[selected_indices]
            else:
                # Unilateral: select points only in the direction from start_label to end_label
                # Determine which side of the extremum corresponds to the forward direction
                # Positive forward_k_reduced means moving from start_label toward end_label
                # Negative forward_k_reduced means moving toward start_label
                # For unilateral fitting, we want points in the forward direction only.
                if local_idx == 0:
                    # extremum at start, all points are forward direction
                    valid_mask = forward_k_reduced >= 0
                elif local_idx == nkpt_branch - 1:
                    # extremum at end, forward direction is actually backward along the branch
                    # but the direction from start to end is forward, so points toward start_label are negative forward_k_reduced
                    # we need points moving from end toward start? Actually unilateral fitting along the specified direction
                    # from start_label to end_label, but extremum at end, so there are no points beyond end.
                    # All points are in the forward direction (from start to end), i.e., forward_k_reduced <= 0.
                    valid_mask = forward_k_reduced <= 0
                else:
                    # extremum in middle, take points toward end_label (forward_k_reduced >= 0)
                    valid_mask = forward_k_reduced >= 0

                k_valid = forward_k_reduced[valid_mask]
                energy_valid = forward_energy[valid_mask]
                if len(k_valid) == 0:
                    # fallback to closest points regardless of direction
                    distances = np.abs(forward_k_reduced)
                    sorted_indices = np.argsort(distances)
                    selected_indices = sorted_indices[:num_fit_points]
                    k_fit_reduced = forward_k_reduced[selected_indices]
                    energy_fit = forward_energy[selected_indices]
                else:
                    # select num_fit_points closest to extremum among valid points
                    distances = np.abs(k_valid)
                    sorted_indices = np.argsort(distances)
                    selected_indices = sorted_indices[
                        : min(num_fit_points, len(k_valid))
                    ]
                    k_fit_reduced = k_valid[selected_indices]
                    energy_fit = energy_valid[selected_indices]

            # Convert reduced k coordinates to physical units (Å⁻¹): k_phys = k_reduced * 2π
            k_fit_phys = k_fit_reduced * 2 * np.pi

            # Quadratic fitting in physical units: E = a * k_phys^2 + b * k_phys + c
            coeffs = np.polyfit(k_fit_phys, energy_fit, 2)
            a, b, c = coeffs[0], coeffs[1], coeffs[2]
            # Second derivative d²E/dk² = 2*a (in eV/Å²)
            curvature = a  # eV/Å²

            from scipy.constants import hbar, e, m_e
            convfactor = hbar**2 / e * 1e20
            effective_mass = convfactor / curvature / m_e if curvature != 0 else np.inf

            # If plot_fit, generate plot
            if plot_fit:
                plt.figure()
                # Plot fitted curve
                k_fine = np.linspace(np.min(k_fit_phys), np.max(k_fit_phys), 100)
                energy_fine = np.polyval(coeffs, k_fine)
                plt.plot(k_fine, energy_fine, "r-", label="Quadratic fit")
                plt.scatter(k_fit_phys, energy_fit, label="Data points")
                plt.xlabel("k (Å⁻¹)")
                plt.ylabel("E (eV)")
                plt.title(f"{edge_type} effective mass fit (band {band_idx})")
                plt.tight_layout()
                plt.legend()
                plt.savefig(f"effective_mass_fit_{edge_type}_band{band_idx}.png", dpi=300)
                plt.close()

            results.append(
                {
                    "effective_mass": effective_mass,
                    "curvature": curvature,
                    "fit_coeffs": coeffs.tolist(),
                    "band_index": band_idx,
                }
            )

        return results


def plot_multiple_bands(
    plot_band_datas: List[Dict[str, Any]],
    mat_names: Optional[List[str]] = None,
    fig_name: str = "band.png",
    emin: float = -10,
    emax: float = 10,
):
    """
    Plot band using the same kpath of different materials into a single figure.
    Args:
        plot_band_datas: A list of dictionary, where each dictionary contains band data for a material and contains the following keys:
            - branches: List[List[str]]: List of branches. e.g.: [['G', 'H'], ['H', 'P'], ['N', 'P']]
            - shifted_kpath_coord: List[List[float]]: Coordinates of kpoints in each branch start from 0
            - band_data: List[np.ndarray]: band data for each branch. Each element is an array of shape (nspin, nkpt_per_branch, nbands)
        mat_names: name for each material, used in legends to distinguish band curves of different materials
        fig_name: name of the plotted figure.
        emin: minimum energy for y-axis limit.
        emax: maximum energy for y-axis limit.
    """
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    import numpy as np

    def get_colors(n):
        """Get colors used in the bandplot"""
        base = ["#d62728", "#1f77b4", "#ff7f0e", "#2ca02c"]
        if n <= len(base):
            return base[:n]
        # Add extra color from HSV if not enough
        return base + [
            mcolors.hsv_to_rgb((i / (n - 2) + 0.1, 0.85, 0.95)) for i in range(n - 2)
        ]

    # Check whether the path of all materials are all same
    for i in range(len(plot_band_datas)):
        for j in range(i + 1, len(plot_band_datas)):
            branches_i, branches_j = (
                plot_band_datas[i]["branches"],
                plot_band_datas[j]["branches"],
            )
            assert len(branches_i) == len(branches_j)
            for ib in range(len(branches_i)):
                assert branches_i[ib][0] == branches_j[ib][0]
                assert branches_i[ib][1] == branches_j[ib][1]

    show_labels = mat_names is not None
    if show_labels:
        assert len(plot_band_datas) == len(mat_names)

    # Determine if single material with spin polarization.
    # Use to provide different color of spin up and spin down bands in ABACUS-agent-tools
    single_material_spin2, spin_colors = False, None
    if len(plot_band_datas) == 1:
        # Check nspin of first branch of first material
        first_data = plot_band_datas[0]
        if len(first_data["band_data"]) > 0:
            first_branch_data = first_data["band_data"][
                0
            ]  # shape (nspin, nkpt, nbands)
            if first_branch_data.shape[0] == 2:  # nspin == 2
                single_material_spin2 = True
                spin_colors = get_colors(2)  # Different colors for spin up and down

    # Assign colors to materials
    colors = get_colors(len(plot_band_datas))

    # Determine high symmetry points positions from the first material
    first_data = plot_band_datas[0]
    high_symm_labels, high_symm_poses = [], []
    start_x_pos = 0.0
    for ipath, branch in enumerate(first_data["branches"]):
        branch_length = np.max(first_data["shifted_kpath_coord"][ipath])
        start_label, end_label = branch[0], branch[1]
        if ipath == 0:
            high_symm_labels.append(start_label)
            high_symm_poses.append(start_x_pos)
        # Check continuity with previous branch
        if len(high_symm_labels) > 0 and high_symm_labels[-1] != start_label:
            # discontinuous: merge start label with previous label
            high_symm_labels[-1] += f"|{start_label}"
            high_symm_labels.append(end_label)
            high_symm_poses.append(start_x_pos + branch_length)
        else:
            # continuous: add end label only (start already exists)
            high_symm_labels.append(end_label)
            high_symm_poses.append(start_x_pos + branch_length)
        start_x_pos += branch_length

    # Plot band for each material
    for imat, plot_band_data in enumerate(plot_band_datas):
        # Determine colors for this material
        if single_material_spin2 and imat == 0 and spin_colors is not None:
            # Set different color for spin up and down band if only 1 material, used in ABACUS-agent-tools
            spin_up_color = spin_colors[0]
            spin_down_color = spin_colors[1]
        else:
            spin_up_color = colors[imat]
            spin_down_color = colors[imat]  # same color for spin down (if any)

        start_x_pos = 0.0
        for ipath, branch in enumerate(plot_band_data["branches"]):
            shifted_coord = plot_band_data["shifted_kpath_coord"][ipath]
            branch_data = plot_band_data["band_data"][
                ipath
            ]  # shape (nspin, nkpt, nbands)
            nspin = branch_data.shape[0]  # Merge treatment for nspin=1 and 4
            nbands = branch_data.shape[2]
            for iband in range(nbands):
                # spin up (or nspin=1)
                plt.plot(
                    shifted_coord + start_x_pos,
                    branch_data[0, :, iband],
                    "-",
                    color=spin_up_color,
                    linewidth=1.0,
                    label=mat_names[imat]
                    if (show_labels and ipath == 0 and iband == 0)
                    else None,
                )
                if nspin == 2:  # spin down
                    spin_down_label = None
                    if (
                        show_labels
                        and single_material_spin2
                        and imat == 0
                        and ipath == 0
                        and iband == 0
                    ):
                        spin_down_label = f"{mat_names[imat]} (down)"  # Provide labe for only the first segment of band if nspin=2
                    plt.plot(
                        shifted_coord + start_x_pos,
                        branch_data[1, :, iband],
                        "--",
                        color=spin_down_color,
                        linewidth=1.0,
                        label=spin_down_label,
                    )
            start_x_pos += np.max(shifted_coord)

    # Set plot limits and labels
    plt.xlim(0, high_symm_poses[-1])
    plt.ylim(emin, emax)
    plt.ylabel(r"$E-E_\text{F}$/eV")
    plt.xticks(high_symm_poses, high_symm_labels)
    for pos in high_symm_poses:
        plt.axvline(pos, color="black", linestyle="-", lw=0.5, alpha=0.5)
    plt.axhline(0, color="black", linestyle="--", lw=0.5)
    # Add legend if multiple materials or single material with spin polarization, and labels are requested
    if show_labels and (len(plot_band_datas) > 1 or single_material_spin2):
        plt.legend()
    plt.title("Band structure")
    plt.tight_layout()
    plt.savefig(fig_name, dpi=300)
    plt.close()
