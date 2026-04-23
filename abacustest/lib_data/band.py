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
        start_label, end_label = segment["start"], segment["end"]
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
                    label = label[1:].strip()
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
                if i < len(high_symm_labels_list) - 1 and insert_kpt_nums[i + 1] == 1:
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

        # Group by spin, kpoint, and band index
        # Structure: {spin: {k: [bands]}}
        spin_k_bands = {}
        kpoint_indices = set()

        for spin, k, b in spin_indices:
            if spin not in spin_k_bands.keys():
                spin_k_bands[spin] = dict({})
            if k not in spin_k_bands[spin].keys():
                spin_k_bands[spin][k] = []

            spin_k_bands[spin][k].append(b)
            kpoint_indices.add(k)

        # Convert to three-level list structure
        band_index = []
        for spin in sorted(spin_k_bands.keys()):
            spin_bands = []
            for k in sorted(spin_k_bands[spin].keys()):
                # Append the list of band indices for this (spin, k)
                spin_bands.append(sorted(spin_k_bands[spin][k]))
            band_index.append(spin_bands)

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
            "band_index": band_index,
            "kpoint_index": kpoint_indices,
            "kpoint_labels": kpoint_labels,
            "kpoint_coord": kpoint_coord,
            "energy": energy,
        }

    def _find_band_edge(self, below_fermi=True, tol=1e-6):
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

    def _get_edge(self, below_fermi, tol=1e-6, spin_resolved=False):
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

    def get_vbm(self, tol=1e-6, spin_resolved=False):
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

    def get_cbm(self, tol=1e-6, spin_resolved=False):
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

    def get_band_gap(self, tol=1e-6, spin_resolved=False):
        """
        Get the band gap energy.

        Args:
            tol (float): Tolerance for energy comparison (eV).
            spin_resolved (bool): If True, return results for each spin channel separately.

        Returns:
            If spin_resolved=False:
                float: Band gap energy (CBM energy - VBM energy) in eV.
                Returns 0.0 if the system is metallic (no band gap) or if VBM/CBM not found.

            If spin_resolved=True:
                dict with keys for each spin label (e.g., "spin_up", "spin_down" for nspin=2, otherwise str(spin_index)) and "global":
                    - spin_label (str): Band gap for that spin channel (CBM_energy[spin] - VBM_energy[spin]).
                    - "global": Global band gap (global CBM energy - global VBM energy).
                Returns 0.0 for metallic systems or when VBM/CBM not found.
        """
        # Get VBM and CBM data
        vbm_result = self.get_vbm(tol=tol, spin_resolved=spin_resolved)
        cbm_result = self.get_cbm(tol=tol, spin_resolved=spin_resolved)

        def safe_calculate_gap(vbm, cbm):
            gap = self._calculate_gap_from_results(vbm, cbm)
            return gap if gap is not None else 0.0

        if not spin_resolved:
            # Non-spin-resolved case: both results are dictionaries with "energy" key
            return safe_calculate_gap(vbm_result, cbm_result)

        # Spin-resolved case: results are dictionaries with spin labels and "global" key
        result = {}

        # Handle global band gap
        if "global" in vbm_result and "global" in cbm_result:
            result["global"] = safe_calculate_gap(
                vbm_result["global"], cbm_result["global"]
            )
        else:
            result["global"] = 0.0

        # Handle per-spin band gaps
        for spin in range(self.nspin):
            label = self._get_spin_label(spin)
            if label in vbm_result and label in cbm_result:
                result[label] = self._calculate_gap_from_results(
                    vbm_result[label], cbm_result[label]
                )
            else:
                result[label] = 0.0

        return result

    def write_to_file(
        self,
        filename: str = "band.dat",
        spin_up_suffix: str = "_up",
        spin_down_suffix: str = "_down",
    ):
        """
        Write band data to file(s).

        Args:
            filename (str): Base filename for output. Default is "band.dat".
            spin_up_suffix (str): Suffix for spin-up file when nspin=2. Default is "_up".
            spin_down_suffix (str): Suffix for spin-down file when nspin=2. Default is "_down".

        Output format:
            First column: k-path cumulative distance
            Subsequent columns: Band energies for each band
        """
        import os
        from pathlib import Path

        filepath = Path(filename)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        if self.nspin == 1 or self.nspin == 4:
            self._write_single_file(filepath, 0)
        elif self.nspin == 2:
            base_name = filepath.stem
            ext = filepath.suffix

            spin_up_filepath = filepath.parent / f"{base_name}{spin_up_suffix}{ext}"
            self._write_single_file(spin_up_filepath, 0)

            spin_down_filepath = filepath.parent / f"{base_name}{spin_down_suffix}{ext}"
            self._write_single_file(spin_down_filepath, 1)

    def _write_single_file(self, filepath: Path, spin_index: int):
        """
        Helper method to write band data for a single spin channel.

        Args:
            filepath (Path): Path to output file.
            spin_index (int): Index of spin channel to write.
        """
        import numpy as np

        with open(filepath, "w") as f:
            f.write("#The first column is accumulated path length, and the rest columns are band data, and has minused Fermi energy.\n")

        kpath_col = self.kpath_lengths.reshape(-1, 1)
        bands = self.band_data[spin_index, :, :]
        data_to_write = np.hstack([kpath_col, bands])

        with open(filepath, "a") as f:
            np.savetxt(f, data_to_write, fmt="%12.8f")

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
        direction_labels,
        band_index,
        num_fit_points=5,
        plot_fit=False,
    ):
        """
        Calculate effective mass for a specific band along a specified direction.

        Args:
            direction_labels (list): List of two high-symmetry labels, e.g., ['G', 'X'].
                The first label should be the high-symmetry point where the band extremum is located.
            band_index (int): Band index to calculate effective mass for (must be at CBM or VBM).
            num_fit_points (int): Number of k-points to use for quadratic fitting.
            plot_fit (bool): If True, plot the fitting curve and data points.

        Note:
            This method does not support nspin=2 (spin-polarized calculations).

        Returns:
            dict: Effective mass result for the specified band containing:
                - effective_mass (float): effective mass in units of electron mass.
                  Positive for electron-like bands (upward curvature), negative for hole-like bands (downward curvature).
                - fit_coeffs (list): quadratic fit coefficients [a, b, c] for E = a*x**2 + b*x + c.
                - a_std: Standard deviation of a.
        """
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy.constants import hbar, e, m_e
        from scipy.optimize import curve_fit

        def parabola_fit(x, a, b, c):
            return a * x**2 + b * x + c

        if direction_labels is None or len(direction_labels) != 2:
            raise ValueError("direction_labels must be a list of two high-symmetry labels")

        if self.nspin == 2:
            raise ValueError("Effective mass calculation not supported for nspin=2")

        # Get branch data
        band_data, band_seg_coord = self.get_band_branch_data(direction_labels, extra_end_point=False)
        band_seg_data = band_data[0, :, band_index]  # For nspin=1 currently
        if band_seg_data is None:
            raise ValueError(f"Branch {direction_labels} not found in band structure")
        assert band_seg_coord is not None

        if len(band_seg_data) < num_fit_points:
            raise ValueError(f"Not enough valid points for band {band_index}")

        energy_fit, k_coord_fit = band_seg_data[:num_fit_points], band_seg_coord[:num_fit_points]
        k_coord_fit *= 2 * np.pi

        popt, pcov = curve_fit(parabola_fit, k_coord_fit, energy_fit)  # Quadratic fitting
        a, b, c = popt
        p_err = np.sqrt(np.diag(pcov))  # The std deviation of the fit coefficients is the sqrt of the diagonal of the covariance matrix
        a_std = p_err[0]
        curvature = a * 2  # Second derivative d^2E/dk^2 = 2*a (in eV * Angstrom^2)

        convfactor = hbar**2 / e * 1e20  # 1e20 is the conversion factor from Angstrom^-2 to m^-2
        effective_mass = convfactor / curvature / m_e if curvature != 0 else np.inf

        # If plot_fit, generate plot
        if plot_fit:
            plt.figure()
            k_fine = np.linspace(np.min(k_coord_fit), np.max(k_coord_fit), 100)
            energy_fine = np.polyval(popt, k_fine)
            formula_text = f"$E - E_F = {a:.2f}k^2 {b:+.2f}k {c:+.2f}$\n$\\sigma_a = {a_std:.2f}, m^* = {effective_mass:.4f} m_e$"
            plt.plot(k_fine, energy_fine, "r-", label=formula_text)
            plt.plot(k_coord_fit, energy_fit, "kD")
            plt.xlabel(r"k ($\AA^{-1}$)")
            plt.ylabel(r"E - E$_{F}$ (eV)")
            plt.title(f"Effective mass fit: {direction_labels[0]} to {direction_labels[1]} (band {band_index+1})")  # band index starts from 1
            plt.tight_layout()
            plt.legend()
            plt.savefig(f"effective_mass_fit_band_{direction_labels[0]}_to_{direction_labels[1]}_{band_index+1}.png", dpi=300)
            plt.close()

        return {
            "effective_mass": effective_mass,
            "fit_coeffs": [a, b, c],
            "a_std": a_std,
        }

    def get_effective_mass_at_edge(
        self,
        direction_labels,
        edge_type="CBM",
        num_fit_points=5,
        plot_fit=False,
    ):
        """
        Calculate effective mass for all degenerate bands at CBM or VBM along a specified direction.

        Args:
            direction_labels (list): List of two high-symmetry labels, e.g., ['G', 'X'].
                The first label should be the high-symmetry point where the band extremum is located.
            edge_type (str): "CBM" or "VBM".
            num_fit_points (int): Number of k-points to use for quadratic fitting.
            plot_fit (bool): If True, plot the fitting curve and data points for each band.

        Note:
            This method does not support nspin=2 (spin-polarized calculations).

        Returns:
            list: List of effective mass results for each degenerate band at the edge.
                  Each element is a dict containing:
                - effective_mass (float): effective mass in units of electron mass.
                  Positive for electron-like bands (upward curvature), negative for hole-like bands (downward curvature).
                - curvature (float): second derivative d²E/dk² in eV/Å².
                  Positive for upward curvature, negative for downward curvature.
                - fit_coeffs (list): quadratic fit coefficients [a, b, c] for E = a*(k - k0)^2 + b.
                - band_index (int): band index used.
        """
        # Get CBM/VBM data
        edge_type_upper = edge_type.upper()
        if edge_type_upper == "CBM":
            edge_result = self.get_cbm(spin_resolved=False)
        elif edge_type_upper == "VBM":
            edge_result = self.get_vbm(spin_resolved=False)
        else:
            raise ValueError("edge_type must be 'CBM' or 'VBM'")

        # Get band indices (for nspin=1, band_index dict maps spin 0 to list of bands)
        # There may be multiple degenerate bands at CBM/VBM
        band_indices = list(edge_result["band_index"].values())[0]

        # Calculate effective mass for each band at the edge
        results = []
        for band_idx in band_indices:
            result = self.get_effective_mass(
                direction_labels=direction_labels,
                band_index=band_idx,
                num_fit_points=num_fit_points,
                plot_fit=plot_fit,
            )
            results.append(result)

        return results

    def get_high_symmetry_points_info(self):
        """
        Get high symmetry points information including k-point number, label and cumulative distance.
        
        Returns:
            List[Dict]: List of dictionaries with keys 'kpoint_num', 'label', 'cumulative_distance'
        """
        high_symm_info = []
        seen_points = set()

        # Process the first kpath start point
        first_kpath = self.kpaths[0]
        start_label = first_kpath["start"]
        start_idx = first_kpath["start_nkpt"]
        cum_dist = self.kpath_lengths[start_idx]
        kpoint_num = start_idx + 1  # Convert to 1-based index
        high_symm_info.append(
            {
                "kpoint_num": kpoint_num,
                "label": start_label,
                "cumulative_distance": float(cum_dist),
            }
        )
        seen_points.add((start_label, start_idx))

        # Process all kpaths
        for i, kpath in enumerate(self.kpaths):
            end_label = kpath["end"]
            end_idx = kpath["end_nkpt"]

            # Add end label of current kpath (if not already added)
            if (end_label, end_idx) not in seen_points:
                cum_dist = self.kpath_lengths[end_idx]
                kpoint_num = end_idx + 1  # Convert to 1-based index
                high_symm_info.append(
                    {
                        "kpoint_num": kpoint_num,
                        "label": end_label,
                        "cumulative_distance": float(cum_dist),
                    }
                )
                seen_points.add((end_label, end_idx))

            # Check if this is a discontinuous point (different from next start)
            if i < len(self.kpaths) - 1:
                next_kpath = self.kpaths[i + 1]
                next_start_label = next_kpath["start"]
                next_start_idx = next_kpath["start_nkpt"]

                if end_label != next_start_label:
                    # Discontinuous point: add next start label separately
                    if (next_start_label, next_start_idx) not in seen_points:
                        cum_dist_next = self.kpath_lengths[next_start_idx]
                        kpoint_num_next = next_start_idx + 1  # Convert to 1-based index
                        high_symm_info.append(
                            {
                                "kpoint_num": kpoint_num_next,
                                "label": next_start_label,
                                "cumulative_distance": float(cum_dist_next),
                            }
                        )
                        seen_points.add((next_start_label, next_start_idx))
                # If continuous, next_start_label is same as end_label, already added
            # Last kpath end point already added above

        return high_symm_info

    def print_high_symmetry_labels(self):
        """
        Print high symmetry labels and coordinates to console.
        """
        print("-" * 60 + f"\n{'label':<10} {'Fractional coordinates in reciprocal space':<20}\n" + "-" * 60)
        for label, pos in self.high_symm_labels.items():
            print(f"{label:<10} ({pos[0]:12.9f}, {pos[1]:12.9f}, {pos[2]:12.9f})")
        print("-" * 60)

    def print_kpath_info(self):
        """
        Print high symmetry points information to console.
        """
        high_symm_info = self.get_high_symmetry_points_info()
        
        print("\nBand path information:")
        print("-" * 60)
        print(f"{'K-point No.':<12} {'Label':<15} {'Cumulative distance':<20}")
        print("-" * 60)
        
        for info in high_symm_info:
            print(f"{info['kpoint_num']:<12} {info['label']:<15} {info['cumulative_distance']:<20.8f}")
        
        print("-" * 60)

    def write_kpath_info(self, filename="KPATH.txt"):
        """
        Write high symmetry points information to a file.
        
        Args:
            filename (str): Output filename
        """
        high_symm_info = self.get_high_symmetry_points_info()
        
        filepath = Path(filename)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        with open(filepath, "w") as f:
            f.write("High symmetry points information\n")
            f.write("-" * 60 + "\n")
            f.write(f"{'K-point No.':<12} {'Label':<15} {'Cumulative distance':<20}\n")
            f.write("-" * 60 + "\n")
            for info in high_symm_info:
                f.write(f"{info['kpoint_num']:<12} {info['label']:<15} {info['cumulative_distance']:<20.8f}\n")
        
        return str(filepath)


class ProjBandData(BandData):
    """
    Projected band data class.
    """

    def __init__(
        self,
        high_symm_labels: Dict[str, List[float]],
        kpaths: List[Dict[str, Any]],
        kpath_cum_dist: Union[np.ndarray, List[float]],
        efermi: Optional[float] = None,
        band_data: Optional[Union[np.ndarray, List[np.ndarray]]] = None,
        proj_band_data: Optional[List[Dict[str, Any]]] = None,
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
            proj_band_data (List[Dict[str, Any]]): Projected band data for each k-path. Should be a list of norb dictionaries. Each dict contains:
                - "data": np.ndarray shaped (nspin, nkpts, nbands) for each spin channel.
                - "atom_index: index of the atom starts from 1.
                - "l": orbital quantum number.
                - "m": magnetic quantum number.
                - "z": zeta of the orbital for the same l and m.
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
                self.band_data = band_data
        else:
            raise ValueError("band_data must be an np.ndarray with shape (1, nkpts, nbands) or (2, nkpts, nbands)")

        self.nspin, self.nkpts, self.nbands = self.band_data.shape

        for pbd in proj_band_data:
            nspin, nkpts, nbands = pbd["data"].shape
            assert nspin == self.nspin
            assert nkpts == self.nkpts
            assert nbands == self.nbands

        self.proj_band_data = proj_band_data
        self._build_label_mapping()
        self._build_kpath_segments()

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
        import xml.etree.ElementTree as ET
        from io import StringIO

        band_data = BandData.ReadFromAbacusJob(abacusjob_dir, efermi, high_symm_labels)

        input_params = ReadInput(os.path.join(abacusjob_dir, "INPUT"))
        abacusjob_outdir = os.path.join(abacusjob_dir, f"OUT.{input_params.get('suffix', 'ABACUS')}")
        proj_band_file_up = os.path.join(abacusjob_outdir, "PBANDS_1")
        tree_up = ET.parse(proj_band_file_up)
        if input_params.get("nspin", 1) == 2:
            proj_band_file_dn = os.path.join(abacusjob_outdir, "PBANDS_2")
            tree_dn = ET.parse(proj_band_file_dn)

        all_orbital_projband_data = []

        # Read projected band data for nspin=1 case
        root_up = tree_up.getroot()
        root_up_orbs = root_up.findall("orbital")
        if input_params.get("nspin", 1) == 2:
            root_dn = tree_dn.getroot()
            root_dn_orbs = root_dn.findall("orbital")

        for iorb, orb in enumerate(root_up_orbs):
            raw_pband_data_up = np.loadtxt(StringIO(orb.find("data").text))
            if input_params.get("nspin", 1) in [1, 4]:
                raw_pband_data = raw_pband_data_up[np.newaxis, :, :]
            if input_params.get("nspin", 1) == 2:
                raw_pband_data_dn = np.loadtxt(StringIO(root_dn_orbs[iorb].find("data").text))
                raw_pband_data = np.array([raw_pband_data_up, raw_pband_data_dn])

            # Ignore Band Data in PBANDS_* for simplicity
            orbital_info = {
                "index": int(orb.get("index")),
                "atom_index": int(orb.get("atom_index")),
                "species": orb.get("species"),
                "l": int(orb.get("l")),
                "m": int(orb.get("m")),
                "z": int(orb.get("z")),
                "data": raw_pband_data,
            }

            all_orbital_projband_data.append(orbital_info)

        return ProjBandData(
            high_symm_labels,
            band_data.kpaths,
            band_data.kpath_lengths,
            efermi,
            band_data.band_data,
            all_orbital_projband_data,
        )

    # Constants for orbital mapping (same as in comm_dos.py)
    l_map = ["s", "p", "d", "f", "g"]
    orbital_names = {
        (0, 0): "s",
        (1, 0): r"$p_z$",
        (1, 1): r"$p_x$",
        (1, 2): r"$p_y$",
        (2, 0): r"$d_{z^2}$",
        (2, 1): r"$d_{xz}$",
        (2, 2): r"$d_{yz}$",
        (2, 3): r"$d_{x^2-y^2}$",
        (2, 4): r"$d_{xy}$",
        (3, 0): r"$f_{z^3}$",
        (3, 1): r"$f_{xz^2}$",
        (3, 2): r"$f_{yz^2}$",
        (3, 3): r"$f_{zx^2-zy^2}$",
        (3, 4): r"$f_{xyz}$",
        (3, 5): r"$f_{x^3-3xy^2}$",
        (3, 6): r"$f_{3yx^2-y^3}$",
        (4, 0): r"$g_1$",
        (4, 1): r"$g_2$",
        (4, 2): r"$g_3$",
        (4, 3): r"$g_4$",
        (4, 4): r"$g_5$",
        (4, 5): r"$g_6$",
        (4, 6): r"$g_7$",
        (4, 7): r"$g_8$",
        (4, 8): r"$g_9$",
    }

    @staticmethod
    def sum_proj_data(proj_datas: List[Dict]) -> np.ndarray:
        """
        Sum the projected band data from a list of extracted projected band data.

        Parameters:
        -----------
        proj_datas : list of dict
            List of projected band data dictionaries

        Returns:
        --------
        np.ndarray
            Summed projected band data with shape (nspin, nkpts, nbands)
        """
        if len(proj_datas) == 0:
            raise ValueError("No projected band data provided")

        proj_sum = np.zeros_like(proj_datas[0]["data"])
        for proj_data in proj_datas:
            proj_sum += np.array(proj_data["data"])
        return proj_sum

    @classmethod
    def _get_orbital_name(cls, l: int, m: int) -> str:
        """Get orbital name from l and m quantum numbers."""
        return cls.orbital_names.get((l, m), f"l{l}_m{m}")

    @classmethod
    def _parse_orbital_name(cls, name: str) -> Tuple[int, int]:
        """Parse orbital name to (l, m) tuple.

        Args:
            name: Orbital name (e.g., 's', 'px', 'py', 'pz', 'dz2', etc.)

        Returns:
            (l, m) tuple

        Raises:
            ValueError: If name cannot be parsed
        """
        # Simple mapping for common orbital names
        name_to_lm = {
            "s": (0, 0),
            "px": (1, 1),
            "py": (1, 2),
            "pz": (1, 0),
            "dz2": (2, 0),
            "dxz": (2, 1),
            "dyz": (2, 2),
            "dx2-y2": (2, 3),
            "dxy": (2, 4),
            # Add more mappings as needed
        }

        # Try exact match
        if name in name_to_lm:
            return name_to_lm[name]

        # Try case-insensitive match
        name_lower = name.lower()
        for key, value in name_to_lm.items():
            if key.lower() == name_lower:
                return value

        # Try to match against orbital_names values (LaTeX format)
        for (l, m), orbital_name in cls.orbital_names.items():
            # Remove LaTeX formatting for comparison
            clean_name = (
                orbital_name.replace("$", "")
                .replace("{", "")
                .replace("}", "")
                .replace("^", "")
                .replace("_", "")
            )
            if clean_name == name or clean_name.replace(" ", "") == name.replace(
                " ", ""
            ):
                return (l, m)

        raise ValueError(f"Cannot parse orbital name: {name}")

    def get_proj_by_orbital(
        self,
        species: str,
        atom_index: int,
        l: int,
        m: Optional[int] = None,
        z: Optional[int] = None,
    ) -> Dict:
        """
        Get projected band data for a specific orbital.

        Parameters:
        -----------
        species : str
            Species name (e.g. 'Fe', 'O')
        atom_index : int
            Index of the atom (1-based as in ABACUS)
        l : int
            Angular momentum quantum number (0=s, 1=p, 2=d, 3=f, ...)
        m : int, optional
            Magnetic quantum number
        z : int, optional
            Orbital split index

        Returns:
        --------
        dict or None
            Dictionary containing orbital information and projected band data
        """
        for orbital in self.proj_band_data:
            if (
                orbital["species"] == species
                and orbital["atom_index"] == atom_index
                and orbital["l"] == l
            ):
                # If m and z are specified, match those as well
                if m is not None and orbital["m"] != m:
                    continue
                if z is not None and orbital["z"] != z:
                    continue

                return orbital

        return None

    def get_proj_by_atom(
        self, atom_index: int, sum_only: bool = True
    ) -> Union[np.ndarray, List[Dict]]:
        """
        Get all projected band data for a specific atom.

        Parameters:
        -----------
        atom_index : int
            Index of the atom (1-based as in ABACUS)
        sum_only : bool, optional
            If True, sum the projected band data for each orbital, otherwise return the original data. Default is True.

        Returns:
        --------
        np.ndarray or list of dict
            If sum_only is True, returns summed projected band data array with shape (nspin, nkpts, nbands).
            Otherwise, returns list of orbitals belonging to the specified atom.
        """
        proj_datas = [
            orb for orb in self.proj_band_data if orb["atom_index"] == atom_index
        ]
        if sum_only:
            return self.sum_proj_data(proj_datas)
        else:
            return proj_datas

    def get_proj_by_atom_shell(
        self, atom_index: int, l: Union[int, str], sum_only: bool = True
    ) -> Union[np.ndarray, List[Dict]]:
        """
        Get projected band data for a specific atom shell.

        Parameters:
        -----------
        atom_index : int
            Index of the atom (1-based as in ABACUS)
        l : int or str
            Angular momentum quantum number (0 or 's', 1 or 'p', 2 or 'd', 3 or 'f', ...)
        sum_only : bool, optional
            If True, sum the projected band data for each orbital, otherwise return the original data. Default is True.

        Returns:
        --------
        np.ndarray or list of dict
            If sum_only is True, returns summed projected band data array with shape (nspin, nkpts, nbands).
            Otherwise, returns list of orbitals belonging to the specified atom and shell.
        """
        if l in self.l_map:
            l = self.l_map.index(l)
        proj_datas = [
            orb
            for orb in self.proj_band_data
            if (orb["atom_index"] == atom_index and orb["l"] == l)
        ]
        if sum_only:
            return self.sum_proj_data(proj_datas)
        else:
            return proj_datas

    def get_proj_by_atom_orbital(
        self,
        atom_index: int,
        l: Optional[Union[int, str]],
        m: Union[int, str],
        sum_only: bool = True,
    ) -> Union[np.ndarray, List[Dict]]:
        """
        Get projected band data for a specific atom orbital.

        Parameters:
        -----------
        atom_index : int
            Index of the atom (1-based as in ABACUS)
        l : int or str, optional
            Angular momentum quantum number (0 or 's', 1 or 'p', 2 or 'd', 3 or 'f', ...).
            If m is explicitly specified by atomic orbital name, l is ignored.
        m : int or str
            Magnetic quantum number used by PBANDS file in ABACUS.
            If m is an integer, it should be between 0 and 2l.
            If m is a string, it should be the name of the atomic orbital.
        sum_only : bool, optional
            If True, sum all the projected band data for the specified atom orbital, otherwise return all the individual projected band data.

        Returns:
        --------
        np.ndarray or list of dict
            If sum_only is True, returns summed projected band data array with shape (nspin, nkpts, nbands).
            Otherwise, returns list of orbitals belonging to the specified atom and orbital.
        """
        if isinstance(m, str):
            # Parse orbital name to get l and m
            try:
                l_parsed, m_parsed = self._parse_orbital_name(m)
                proj_datas = [
                    orb
                    for orb in self.proj_band_data
                    if (
                        orb["atom_index"] == atom_index
                        and orb["l"] == l_parsed
                        and orb["m"] == m_parsed
                    )
                ]
            except ValueError:
                # If cannot parse, try to match against orbital names in data (if available)
                # Note: proj_band_data doesn't have 'atomic_orbital_name' key, so we fall back to empty result
                proj_datas = []
        elif isinstance(m, int):
            if l in self.l_map:
                l = self.l_map.index(l)
            proj_datas = [
                orb
                for orb in self.proj_band_data
                if (orb["atom_index"] == atom_index and orb["l"] == l and orb["m"] == m)
            ]
        else:
            raise ValueError("Invalid m value")

        if sum_only:
            return self.sum_proj_data(proj_datas)
        else:
            return proj_datas

    def get_proj_by_species(
        self, species: str, sum_only: bool = True
    ) -> Union[np.ndarray, List[Dict]]:
        """
        Get projected band data for a specific species.

        Parameters:
        -----------
        species : str
            Species name (e.g. 'Fe', 'O')
        sum_only : bool, optional
            If True, sum the projected band data for each orbital, otherwise return the original data. Default is True.

        Returns:
        --------
        np.ndarray or list of dict
            If sum_only is True, returns summed projected band data array with shape (nspin, nkpts, nbands).
            Otherwise, returns list of orbitals belonging to the specified species.
        """
        proj_datas = [orb for orb in self.proj_band_data if orb["species"] == species]
        if sum_only:
            return self.sum_proj_data(proj_datas)
        else:
            return proj_datas

    def get_proj_by_species_shell(
        self, species: str, l: Union[int, str], sum_only: bool = True
    ) -> Union[np.ndarray, List[Dict]]:
        """
        Get projected band data for shell of a specific species.

        Parameters:
        -----------
        species : str
            Species name (e.g. 'Fe', 'O')
        l : int or str
            Angular momentum quantum number (0 or 's', 1 or 'p', 2 or 'd', 3 or 'f', ...)
        sum_only : bool, optional
            If True, sum the projected band data for each orbital, otherwise return the original data. Default is True.

        Returns:
        --------
        np.ndarray or list of dict
            If sum_only is True, returns summed projected band data array with shape (nspin, nkpts, nbands).
            Otherwise, returns list of orbitals belonging to the specified species and shell.
        """
        if l in self.l_map:
            l = self.l_map.index(l)
        proj_datas = [
            orb
            for orb in self.proj_band_data
            if (orb["species"] == species and orb["l"] == l)
        ]

        if sum_only:
            return self.sum_proj_data(proj_datas)
        else:
            return proj_datas

    def get_proj_by_species_orbital(
        self,
        species: str,
        l: Optional[Union[int, str]],
        m: Union[int, str],
        sum_only: bool = True,
    ) -> Union[np.ndarray, List[Dict]]:
        """
        Get projected band data for a specific species orbital.

        Parameters:
        -----------
        species : str
            Species name (e.g. 'Fe', 'O')
        l : int or str, optional
            Angular momentum quantum number (0 or 's', 1 or 'p', 2 or 'd', 3 or 'f', ...).
            If m is explicitly specified by atomic orbital name, l is ignored.
        m : int or str
            Magnetic quantum number used by PBANDS file in ABACUS.
            If m is an integer, it should be between 0 and 2l.
            If m is a string, it should be the name of the atomic orbital (e.g. 's', 'px', 'py', 'pz', 'dz2').
        sum_only : bool, optional
            If True, sum all the projected band data for the specified species orbital, otherwise return all the individual projected band data.

        Returns:
        --------
        np.ndarray or list of dict
            If sum_only is True, returns summed projected band data array with shape (nspin, nkpts, nbands).
            Otherwise, returns list of orbitals belonging to the specified species and orbital.
        """
        if isinstance(m, str):
            # Parse orbital name to get l and m
            try:
                l_parsed, m_parsed = self._parse_orbital_name(m)
                proj_datas = [
                    orb
                    for orb in self.proj_band_data
                    if (
                        orb["species"] == species
                        and orb["l"] == l_parsed
                        and orb["m"] == m_parsed
                    )
                ]
            except ValueError:
                # If cannot parse, try to match against orbital names in data (if available)
                # Note: proj_band_data doesn't have 'atomic_orbital_name' key, so we fall back to empty result
                proj_datas = []
        elif isinstance(m, int):
            if l in self.l_map:
                l = self.l_map.index(l)
            proj_datas = [
                orb
                for orb in self.proj_band_data
                if (orb["species"] == species and orb["l"] == l and orb["m"] == m)
            ]
        else:
            raise ValueError("Invalid m value")

        if sum_only:
            return self.sum_proj_data(proj_datas)
        else:
            return proj_datas

    def get_species(self) -> List[str]:
        """Get all species in the projected band data."""
        return list(set([orb["species"] for orb in self.proj_band_data]))

    def get_species_shell(self, species: str) -> List[int]:
        """Get all shells for a specific species in the projected band data."""
        return list(
            set([orb["l"] for orb in self.proj_band_data if orb["species"] == species])
        )

    def get_species_shell_orbital(self, species: str, l: int) -> List[int]:
        """Get all orbitals for a specific species and shell in the projected band data."""
        return list(
            set(
                [
                    orb["m"]
                    for orb in self.proj_band_data
                    if (orb["species"] == species and orb["l"] == l)
                ]
            )
        )

    def get_atom_species(self, atom_index: int) -> str:
        """Get the species for a specific atom in the projected band data."""
        species = list(
            set(
                [
                    orb["species"]
                    for orb in self.proj_band_data
                    if orb["atom_index"] == atom_index
                ]
            )
        )
        assert len(species) == 1
        return species[0]

    def get_atom_shell(self, atom_index: int) -> List[int]:
        """Get all shells for a specific atom in the projected band data."""
        return list(
            set(
                [
                    orb["l"]
                    for orb in self.proj_band_data
                    if orb["atom_index"] == atom_index
                ]
            )
        )

    def get_atom_shell_orbital(self, atom_index: int, l: int) -> List[int]:
        """Get all orbitals for a specific atom and shell in the projected band data."""
        return list(
            set(
                [
                    orb["m"]
                    for orb in self.proj_band_data
                    if (orb["atom_index"] == atom_index and orb["l"] == l)
                ]
            )
        )

    # ------------------------------------------------------------------
    # Plotting helpers
    # ------------------------------------------------------------------

    def _get_plot_kpaths(
        self, plot_kpaths: Optional[List[List[str]]]
    ) -> List[List[str]]:
        """Return plot_kpaths, defaulting to all stored kpaths."""
        if plot_kpaths is None:
            return [[kp["start"], kp["end"]] for kp in self.kpaths]
        return plot_kpaths

    def _prepare_proj_plot_data(
        self,
        proj_arrays: List[np.ndarray],
        plot_kpaths: List[List[str]],
    ):
        """
        Slice full-length projection arrays (nspin, nkpts, nbands) into
        per-segment arrays aligned with the k-path segments used for plotting.

        Args:
            proj_arrays: List of arrays shaped (nspin, nkpts, nbands).
            plot_kpaths: List of [start_label, end_label] for each segment.

        Returns:
            branches: same as plot_kpaths (pass-through).
            shifted_kpath_coords: List[np.ndarray] k-distances per segment.
            band_data_segs: List[np.ndarray] energy data per segment.
            proj_weights_list: List[List[np.ndarray]] weights per projection per segment.
        """
        branches = []
        shifted_kpath_coords = []
        band_data_segs = []
        # n_proj lists, each will be filled with one entry per segment
        proj_weights_list: List[List[np.ndarray]] = [[] for _ in proj_arrays]

        for kpath_labels in plot_kpaths:
            branch_band, k_coord = self.get_band_branch_data(
                kpath_labels, extra_end_point=True
            )
            if branch_band is None:
                raise ValueError(
                    f"k-path segment {kpath_labels} not found in the band structure."
                )
            # Determine kpt slice indices for this segment
            start_label, end_label = kpath_labels[0], kpath_labels[1]
            seg_info = None
            for kp in self.kpaths:
                if kp["start"] == start_label and kp["end"] == end_label:
                    seg_info = kp
                    break
                if kp["start"] == end_label and kp["end"] == start_label:
                    # reversed branch – mirror handled by get_band_branch_data
                    seg_info = {
                        "start_nkpt": kp["start_nkpt"],
                        "end_nkpt": kp["end_nkpt"],
                        "_reversed": True,
                    }
                    break

            if seg_info is None:
                raise ValueError(
                    f"k-path segment {kpath_labels} not found when slicing projections."
                )

            s_nkpt = seg_info["start_nkpt"]
            # extra_end_point may add +1 k-point – match branch_band shape
            nkpt_seg = branch_band.shape[1]
            e_nkpt = s_nkpt + nkpt_seg  # exclusive

            branches.append(kpath_labels)
            shifted_kpath_coords.append(k_coord)
            band_data_segs.append(branch_band)

            for iproj, proj_arr in enumerate(proj_arrays):
                seg_proj = proj_arr[:, s_nkpt:e_nkpt, :]  # (nspin, nkpt_seg, nbands)
                if seg_info.get("_reversed"):
                    seg_proj = seg_proj[:, ::-1, :]
                proj_weights_list[iproj].append(seg_proj)

        return branches, shifted_kpath_coords, band_data_segs, proj_weights_list

    @staticmethod
    def _get_proj_colors(n: int) -> List:
        """Return n distinct colors (reuses BandData color palette logic)."""
        import matplotlib.colors as mcolors

        base = [
            "#1f77b4",
            "#d62728",
            "#ff7f0e",
            "#2ca02c",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]
        if n <= len(base):
            return base[:n]
        extra = [
            mcolors.hsv_to_rgb(((i / (n - len(base))) * 0.7 + 0.05, 0.85, 0.90))
            for i in range(n - len(base))
        ]
        return base + extra

    def _call_plot(
        self,
        proj_arrays: List[np.ndarray],
        proj_labels: List[str],
        plot_kpaths: Optional[List[List[str]]],
        emin: float,
        emax: float,
        fig_name: str,
        spin_mode: str,
        bg_alpha: float,
    ) -> None:
        """Common driver: prepare data and call plot_proj_bands_subplots."""
        kpaths = self._get_plot_kpaths(plot_kpaths)
        branches, k_coords, band_segs, proj_weights_list = self._prepare_proj_plot_data(
            proj_arrays, kpaths
        )
        colors = self._get_proj_colors(len(proj_arrays))
        plot_proj_bands_subplots(
            branches=branches,
            shifted_kpath_coord=k_coords,
            band_data=band_segs,
            proj_weights_list=proj_weights_list,
            proj_labels=proj_labels,
            proj_colors=colors,
            emin=emin,
            emax=emax,
            fig_name=fig_name,
            spin_mode=spin_mode,
            bg_alpha=bg_alpha,
        )

    # ------------------------------------------------------------------
    # Public plotting methods
    # ------------------------------------------------------------------

    def plot_proj_band_species(
        self,
        emin: float = -10,
        emax: float = 10,
        fig_name: str = "proj_band_species.png",
        plot_kpaths: Optional[List[List[str]]] = None,
        spin_mode: str = "combined",
        bg_alpha: float = 0.15,
    ) -> None:
        """
        Plot projected band structure with one subplot per species.

        Each subplot shows all orbitals of that species summed together,
        coloured with alpha proportional to the projection weight.

        Args:
            emin: Lower energy bound (eV, relative to Fermi level).
            emax: Upper energy bound (eV, relative to Fermi level).
            fig_name: Output figure filename.
            plot_kpaths: k-path segments to plot (default: all stored segments).
            spin_mode: Spin handling mode – "combined", "two_panel", "up_only", "down_only".
            bg_alpha: Alpha of the grey background band lines.
        """
        species_list = sorted(self.get_species())
        arrays = [self.get_proj_by_species(s, sum_only=True) for s in species_list]
        self._call_plot(
            arrays, species_list, plot_kpaths, emin, emax, fig_name, spin_mode, bg_alpha
        )

    def plot_proj_band_species_shell(
        self,
        emin: float = -10,
        emax: float = 10,
        fig_name: str = "proj_band_species_shell.png",
        plot_kpaths: Optional[List[List[str]]] = None,
        spin_mode: str = "combined",
        bg_alpha: float = 0.15,
    ) -> None:
        """
        Plot projected band structure with one subplot per (species, shell) pair.

        Args:
            emin: Lower energy bound (eV).
            emax: Upper energy bound (eV).
            fig_name: Output figure filename.
            plot_kpaths: k-path segments to plot (default: all stored segments).
            spin_mode: "combined", "two_panel", "up_only", or "down_only".
            bg_alpha: Alpha of the grey background band lines.
        """
        arrays, labels = [], []
        for s in sorted(self.get_species()):
            for l in sorted(self.get_species_shell(s)):
                shell_name = self.l_map[l] if l < len(self.l_map) else f"l{l}"
                arrays.append(self.get_proj_by_species_shell(s, l, sum_only=True))
                labels.append(f"{s}-{shell_name}")
        self._call_plot(
            arrays, labels, plot_kpaths, emin, emax, fig_name, spin_mode, bg_alpha
        )

    def plot_proj_band_species_orbital(
        self,
        emin: float = -10,
        emax: float = 10,
        fig_name: str = "proj_band_species_orbital.png",
        plot_kpaths: Optional[List[List[str]]] = None,
        spin_mode: str = "combined",
        bg_alpha: float = 0.15,
    ) -> None:
        """
        Plot projected band structure with one subplot per (species, l, m) triplet.

        Each magnetic quantum-number channel is shown separately.

        Args:
            emin: Lower energy bound (eV).
            emax: Upper energy bound (eV).
            fig_name: Output figure filename.
            plot_kpaths: k-path segments to plot (default: all stored segments).
            spin_mode: "combined", "two_panel", "up_only", or "down_only".
            bg_alpha: Alpha of the grey background band lines.
        """
        arrays, labels = [], []
        for s in sorted(self.get_species()):
            for l in sorted(self.get_species_shell(s)):
                for m in sorted(self.get_species_shell_orbital(s, l)):
                    orb_name = self.orbital_names.get((l, m), f"l{l}m{m}")
                    # Strip LaTeX for the title
                    orb_name_clean = (
                        orb_name.replace("$", "")
                        .replace("{", "")
                        .replace("}", "")
                        .replace("^", "")
                        .replace("_", "")
                    )
                    arrays.append(
                        self.get_proj_by_species_orbital(s, l, m, sum_only=True)
                    )
                    labels.append(f"{s}-{orb_name_clean}")
        self._call_plot(
            arrays, labels, plot_kpaths, emin, emax, fig_name, spin_mode, bg_alpha
        )

    def plot_proj_band_atoms(
        self,
        atom_indices: List[int],
        emin: float = -10,
        emax: float = 10,
        fig_name: str = "proj_band_atoms.png",
        plot_kpaths: Optional[List[List[str]]] = None,
        spin_mode: str = "combined",
        bg_alpha: float = 0.15,
    ) -> None:
        """
        Plot projected band structure with one subplot per specified atom.

        Each subplot shows all orbitals of that atom summed together.

        Args:
            atom_indices: List of atom indices (1-based) to plot.
            emin: Lower energy bound (eV).
            emax: Upper energy bound (eV).
            fig_name: Output figure filename.
            plot_kpaths: k-path segments to plot (default: all stored segments).
            spin_mode: "combined", "two_panel", "up_only", or "down_only".
            bg_alpha: Alpha of the grey background band lines.
        """
        arrays, labels = [], []
        for idx in atom_indices:
            species = self.get_atom_species(idx)
            arrays.append(self.get_proj_by_atom(idx, sum_only=True))
            labels.append(f"{species}{idx}")
        self._call_plot(
            arrays, labels, plot_kpaths, emin, emax, fig_name, spin_mode, bg_alpha
        )


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


def plot_proj_bands_subplots(
    branches: List[List[str]],
    shifted_kpath_coord: List[np.ndarray],
    band_data: List[np.ndarray],
    proj_weights_list: List[List[np.ndarray]],
    proj_labels: List[str],
    proj_colors: List,
    emin: float = -10,
    emax: float = 10,
    fig_name: str = "proj_band.png",
    spin_mode: str = "combined",
    bg_alpha: float = 0.15,
):
    """
    Plot projected band structure with alpha-weighted color overlays, one subplot per projection.

    Args:
        branches: List of [start_label, end_label] for each k-path segment.
        shifted_kpath_coord: List of 1-D arrays with cumulative k-distances for each segment (starting from 0).
        band_data: List of arrays shaped (nspin, nkpt_per_seg, nbands) for each segment.
        proj_weights_list: Outer list over projections; inner list over k-path segments.
            Each element is an array shaped (nspin, nkpt_per_seg, nbands) with values in [0, 1].
        proj_labels: Label for each projection (used as subplot title).
        proj_colors: Color for each projection.
        emin: Lower energy bound for the plot.
        emax: Upper energy bound for the plot.
        fig_name: Output file name.
        spin_mode: How to handle spin channels.
            - "combined": spin-up solid line + spin-down dashed line in the same panel.
            - "two_panel": separate left/right panels sharing y-axis.
            - "up_only": only spin-up.
            - "down_only": only spin-down.
        bg_alpha: Alpha value for the background grey band lines.
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.collections import LineCollection

    nspin = band_data[0].shape[0]  # actual number of spin channels in data

    # Determine which spin channels to draw and how many panels per projection
    if nspin == 1 or spin_mode in ("up_only",):
        spin_indices = [0]
        two_panel = False
    elif spin_mode == "down_only":
        spin_indices = [1]
        two_panel = False
    elif spin_mode == "two_panel":
        spin_indices = [0, 1]
        two_panel = True
    else:  # "combined"
        spin_indices = list(range(nspin))
        two_panel = False

    # --- Build x-axis info from branches (same logic as plot_multiple_bands) ---
    high_symm_labels_plot, high_symm_poses = [], []
    start_x_pos = 0.0
    for ipath, branch in enumerate(branches):
        branch_length = float(np.max(shifted_kpath_coord[ipath]))
        start_label, end_label = branch[0], branch[1]
        if ipath == 0:
            high_symm_labels_plot.append(start_label)
            high_symm_poses.append(start_x_pos)
        if high_symm_labels_plot[-1] != start_label:
            high_symm_labels_plot[-1] += f"|{start_label}"
            high_symm_labels_plot.append(end_label)
            high_symm_poses.append(start_x_pos + branch_length)
        else:
            high_symm_labels_plot.append(end_label)
            high_symm_poses.append(start_x_pos + branch_length)
        start_x_pos += branch_length
    total_k_length = high_symm_poses[-1]

    n_proj = len(proj_weights_list)

    # Determine subplot grid: panels_per_proj = 2 if two_panel else 1
    panels_per_proj = 2 if two_panel else 1
    n_total_cols_needed = n_proj * panels_per_proj  # used only for two_panel layout

    # For two_panel: ncols = 2*n_proj; rows = 1 (each projection fills one pair)
    # For non-two_panel: use PDOS-style nrow/ncol logic
    START_MULTI_COL = 4
    TWO_COL_MAX = 6

    if two_panel:
        # Each projection occupies a pair of side-by-side axes sharing y
        nrow, ncol = 1, n_total_cols_needed
        if n_proj > 3:  # stack rows if many projections
            nrow = n_proj
            ncol = 2
    else:
        if n_proj >= START_MULTI_COL:
            if n_proj > TWO_COL_MAX:
                ncol = max(3, int(np.sqrt(n_proj)))
            else:
                ncol = 2
            nrow = int(np.ceil(n_proj / ncol))
        else:
            nrow, ncol = n_proj, 1

    # --- Helper: draw background and projection bands into an Axes ---
    def _draw_single_panel(ax, spin_idx, proj_weight_segs, proj_color, linestyle="-"):
        """Draw background grey bands + colored projection bands for one spin into ax."""
        x_offset = 0.0
        for iseg, seg_coord in enumerate(shifted_kpath_coord):
            seg_energy = band_data[iseg][spin_idx]  # (nkpt, nbands)
            seg_weight = np.clip(
                proj_weight_segs[iseg][spin_idx], 0.0, 1.0
            )  # (nkpt, nbands)
            x = seg_coord + x_offset
            nbands = seg_energy.shape[1]

            for iband in range(nbands):
                y = seg_energy[:, iband]
                w = seg_weight[:, iband]

                # Background grey line
                ax.plot(
                    x,
                    y,
                    color="grey",
                    alpha=bg_alpha,
                    linewidth=0.8,
                    linestyle=linestyle,
                    zorder=1,
                )

                # Colored projection line via LineCollection (alpha per segment)
                # Build segments: each segment connects two adjacent k-points
                points = np.array([x, y]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                # alpha of each segment = mean weight of its two endpoints, clipped
                seg_alphas = np.clip(0.5 * (w[:-1] + w[1:]), 0.0, 1.0)

                # Build RGBA array for the collection
                r, g, b, _ = mcolors.to_rgba(proj_color)
                rgba = np.zeros((len(segments), 4))
                rgba[:, 0] = r
                rgba[:, 1] = g
                rgba[:, 2] = b
                rgba[:, 3] = seg_alphas

                lc = LineCollection(
                    segments, colors=rgba, linewidth=1.2, linestyle=linestyle, zorder=2
                )
                ax.add_collection(lc)

            x_offset += float(np.max(seg_coord))

    # --- Create figure and axes ---
    fig_width = 6 * ncol
    fig_height = 4 * nrow

    if two_panel:
        # Each row: a pair of axes sharing y (spin_up | spin_down)
        fig, axes_grid = plt.subplots(
            nrow,
            ncol,
            figsize=(fig_width, fig_height),
            sharey="row",
        )
        if nrow == 1:
            # axes_grid shape: (ncol,) – reshape to (1, ncol)
            axes_grid = axes_grid[np.newaxis, :]
        # axes_grid[iproj_row, 0] = spin_up panel; [iproj_row, 1] = spin_down panel
        # when n_proj > 3, nrow=n_proj, ncol=2
        flat_axes_up = [axes_grid[i, 0] for i in range(nrow)]
        flat_axes_down = [axes_grid[i, 1] for i in range(nrow)]
    else:
        fig, axes_flat = plt.subplots(
            nrow,
            ncol,
            figsize=(fig_width, fig_height),
            squeeze=False,
        )
        axes_list = axes_flat.flatten().tolist()

    spin_labels = {0: "↑", 1: "↓"}

    for iproj, (proj_weight_segs, label, color) in enumerate(
        zip(proj_weights_list, proj_labels, proj_colors)
    ):
        if two_panel:
            ax_up = flat_axes_up[iproj]
            ax_down = flat_axes_down[iproj]

            _draw_single_panel(ax_up, 0, proj_weight_segs, color, linestyle="-")
            _draw_single_panel(ax_down, 1, proj_weight_segs, color, linestyle="-")

            for ax, sidx in [(ax_up, 0), (ax_down, 1)]:
                _configure_ax(
                    ax,
                    high_symm_poses,
                    high_symm_labels_plot,
                    total_k_length,
                    emin,
                    emax,
                    title=f"{label} ({spin_labels[sidx]})",
                    show_ylabel=(sidx == 0),
                )
        else:
            ax = axes_list[iproj]
            if spin_mode == "combined" and nspin == 2:
                _draw_single_panel(ax, 0, proj_weight_segs, color, linestyle="-")
                _draw_single_panel(ax, 1, proj_weight_segs, color, linestyle="--")
                # Legend handles for spin up/down
                from matplotlib.lines import Line2D

                handles = [
                    Line2D([0], [0], color=color, lw=1.2, linestyle="-", label="↑"),
                    Line2D([0], [0], color=color, lw=1.2, linestyle="--", label="↓"),
                ]
                ax.legend(handles=handles, fontsize=8, loc="upper right")
            else:
                sidx = spin_indices[0]
                _draw_single_panel(ax, sidx, proj_weight_segs, color, linestyle="-")

            _configure_ax(
                ax,
                high_symm_poses,
                high_symm_labels_plot,
                total_k_length,
                emin,
                emax,
                title=label,
                show_ylabel=True,
            )

    # Hide unused axes
    if not two_panel:
        for idx in range(n_proj, len(axes_list)):
            axes_list[idx].set_visible(False)

    plt.tight_layout()
    plt.savefig(fig_name, dpi=300)
    plt.close()


def _configure_ax(
    ax,
    high_symm_poses,
    high_symm_labels,
    total_k_length,
    emin,
    emax,
    title,
    show_ylabel,
):
    """Apply common axis formatting for a projected band subplot."""
    ax.set_xlim(0, total_k_length)
    ax.set_ylim(emin, emax)
    ax.set_xticks(high_symm_poses)
    ax.set_xticklabels(high_symm_labels)
    for pos in high_symm_poses:
        ax.axvline(pos, color="black", linestyle="-", lw=0.5, alpha=0.5)
    ax.axhline(0, color="black", linestyle="--", lw=0.5)
    if show_ylabel:
        ax.set_ylabel(r"$E-E_\text{F}$/eV")
    ax.set_title(title)
