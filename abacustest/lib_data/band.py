import os

from typing import Dict, List, Tuple, Optional, Union, Any
from pathlib import Path

import numpy as np

from abacustest import ReadInput, ReadKpt, RESULT


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
        if nspin not in (1, 2):
            raise NotImplementedError("Band plot for nspin=4 is not supported yet")

        kpt_result = ReadKpt(abacusjob_dir)
        if kpt_result is None:
            raise ValueError(f"Failed to read KPT file from {abacusjob_dir}")
        kpt_data, model = kpt_result
        if model not in ["line", "line_cartesian"]:
            raise ValueError(f"KPT file must be in 'line' or 'line_cartesian' mode for band calculation, got '{model}'")

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
                raise TypeError(f"Expected kpt to be a list or tuple, got {type(kpt)} at index {i}")

            if len(kpt) < 4:
                raise ValueError(f"Must have at least 4 elements in line-mode KPT file, got {len(kpt)} in line {i}")

            kx, ky, kz, npoints = kpt[:4]
            if model == "line":
                kpt_coord_direct = [float(kx), float(ky), float(kz)]
            else:
                # TODO: transform to cartesian coordinate in recipord space
                # For now, treat as direct coordinates
                kpt_coord_direct = [float(kx), float(ky), float(kz)]

            # Get label
            if len(kpt) >= 5:
                # Use label is provided in KPT file
                label = str(kpt[4]).strip()
                if label.startswith("#"):
                    label = label[1:]
            elif high_symm_labels is not None and i < len(high_symm_labels):
                # Use provided labels
                label = high_symm_labels[i]
            else:
                # Generate default label
                label = f"K{i}"

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
                raise ValueError("Fermi energy (efermi) not found in ABACUS results and not provided")

        # Read band data for nspin=1 or 4
        band_file = os.path.join(abacusjob_dir, f"OUT.{suffix}/BANDS_1.dat")
        with open(band_file, "r") as f:
            original_band_data = np.loadtxt(f)
            kpath_cum_dist, band = original_band_data[:, 1], original_band_data[:, 2:]
            band_data = np.expand_dims(band, axis=0) # Transform to ndarray shaped (1, nkpt, nband)
        # Read band data for nspin=2
        if nspin == 2:
            band_file_dw = os.path.join(abacusjob_dir, f"OUT.{suffix}/BANDS_2.dat")
            with open(band_file_dw, "r") as f:
                original_band_data = np.loadtxt(f)
                band_dw = original_band_data[:, 2:]
                band_data = np.stack((band, band_dw), axis=0) # Transform to ndarray shaped (2, nkpt, nband)

        return BandData(high_symm_kpts, kpaths, kpath_cum_dist, efermi, band_data)

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
            if kpath['start'] == start and kpath['end'] == end: # Find the matched branch
                start_nkpt, end_nkpt = kpath['start_nkpt'], kpath['end_nkpt']
                if extra_end_point:
                    if ikpath < len(self.kpaths) - 1: # Not the last band branch
                        if kpath['end'] == self.kpaths[ikpath+1]['start']: # Continus with the next branch
                            end_nkpt += 1
                        else:
                            start_nkpt, end_nkpt = kpath['start_nkpt'], kpath['end_nkpt']

                branch_band_data = self.band_data[:, start_nkpt:end_nkpt+1, :]
                kpath_coord = self.kpath_lengths[start_nkpt:end_nkpt+1]
                shifted_kpath_coord = kpath_coord - min(kpath_coord)
                return branch_band_data, shifted_kpath_coord
            elif kpath['start'] == end and kpath['end'] == start:
                # Revert the reqested branch to get data, and reverse the band data
                reverse_branch_band_data, shifted_kpath_coord = self.get_band_branch_data([end, start], extra_end_point)
                if reverse_branch_band_data is None or shifted_kpath_coord is None:
                    return None, None
                return reverse_branch_band_data[:, ::-1, :], shifted_kpath_coord

        return None, None  # If no matched branch found

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
            branch_band_data, shifted_coord = self.get_band_branch_data(plot_kpath, extra_end_point=True)
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
            first_branch_data = first_data["band_data"][0]  # shape (nspin, nkpt, nbands)
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
                    label=mat_names[imat] if (show_labels and ipath == 0 and iband == 0) else None,
                )
                if nspin == 2: # spin down
                    spin_down_label = None
                    if show_labels and single_material_spin2 and imat == 0 and ipath == 0 and iband == 0:
                        spin_down_label = f"{mat_names[imat]} (down)"
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
