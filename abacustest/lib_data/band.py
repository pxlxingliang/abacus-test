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
        kpaths: List[Dict[str, Union[int, str]]],
        kpath_cum_dist: Union[np.ndarray, List[float]],
        efermi: Optional[float] = None,
        band_data: Optional[Union[np.ndarray, List[np.ndarray]]] = None,
    ):
        """
        Initialize BandData object.
        Args:
            high_symm_labels (Dict[str, List[float]]): High symmetry labels and coordinates. For example: {"GAMMA": [0.0, 0.0, 0.0], "X": [0.5, 0.0, 0.0], "L": [0.5, 0.5, 0.0]}.
            kpaths (List[Dict[str, Union[float, str]]]): K-path information. For example: [{"start": "GAMMA", "end": "X", "start_nkpt": 0, "end_nkpt": 30}, {"start": "X", "end": "U", "start_nkpt": 30, "end_nkpt": 60}].
                For each kpath, band data are stored in the range [start_nkpt, end_nkpt) of `band_data` for each spin channel.
            kpath_cum_dist (Union[np.ndarray, List[float]]): K-path cumulative distance. For example: [0.0, 0.1, 0.2, 0.3, 0.8, 0.85, 0.9]. Discontinuity in k-path distance will be automatically removed.
            efermi (float): Fermi energy used during the initialization of the BandData object. Will subtract this value from the band data. If None, the band data will not be modified.
            band_data (np.ndarray): Energy level data for each k-point. For nspin=1 or 4, it is an np.ndarray shaped (1, nkpts, nbands). For nspin=2, it is shaped (2, nkpts, nbands).
        """
        self.high_symm_labels = high_symm_labels
        self.kpaths = kpaths
        self.kpath_lengths = self._remove_jump_dist(kpath_cum_dist, kpaths)

        self.efermi = efermi

        if isinstance(band_data, np.ndarray) and len(band_data.shape) == 3 and band_data.shape[0] in (1, 2):
            if efermi is not None:
                self.band_data = band_data - efermi
            else:
                print("efermi is None, band data will not be shifted")
                self.band_data = band_data
        else:
            raise ValueError("band_data must be an np.ndarray with shape (1, nkpts, nbands) or (2, nkpts, nbands)")

    @staticmethod
    def _remove_jump_dist(
        kpath_cum_dist: np.ndarray,
        kpaths: List[Dict[str, Union[int, str]]],
    ) -> np.ndarray:
        """
        Remove the jump in the k-path distance.
        Args:
            kpath_cum_dist (Union[np.ndarray, List[float]]): K-path cumulative distance. For example: [0.0, 0.1, 0.2, 0.3, 0.8, 0.85, 0.9].
            kpaths (List[Dict[str, Union[int, str]]]): K-path information. For example: [{"start": "GAMMA", "end": "X", "start_nkpt": 0, "end_nkpt": 30}, {"start": "X", "end": "U", "start_nkpt": 30, "end_nkpt": 60}]
        Returns:
            np.ndarray: K-path cumulative distance without jump
        """
        for ipath in range(len(kpaths) - 1):
            # Check if there's a jump in cumulative distance between segments
            if kpaths[ipath]["end"] != kpaths[ipath + 1]["start"]:
                jump_dist = kpath_cum_dist[kpaths[ipath+1]["start_nkpt"]] - kpath_cum_dist[kpaths[ipath]["end_nkpt"] - 1]
                print(jump_dist, kpaths[ipath]["end_nkpt"]-1, kpaths[ipath+1]["start_nkpt"])
                kpath_cum_dist[kpaths[ipath+1]["start_nkpt"]:] -= jump_dist
        
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

        kpt_data, model = ReadKpt(abacusjob_dir)
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
                #TODO: transform to cartesian coordinate in recipord space
                pass

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
                end_nkpt = start_nkpt + insert_kpt_nums[i]

                # If the next segment has only 1 k-point, include it in this segment.
                # Used to treat discountinuous k-paths and end points
                if i < len(high_symm_labels_list) - 2 and insert_kpt_nums[i+1] == 1:
                    end_nkpt += 1

                kpaths.append(
                    {
                        "start": high_symm_labels_list[i],
                        "end": high_symm_labels_list[i+1],
                        "start_nkpt": start_nkpt,
                        "end_nkpt": end_nkpt,
                    }
                )
                start_nkpt = end_nkpt

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

    def plot_band(
        self, emin: float = -10, emax: float = 10, fig_name: str = "band.png",
        plot_kpaths: Optional[List[str]] = None,
    ):
        """
        Plot the band structure.
        Args:
            emin (float): Minimum energy in the band plot.
            emax (float): Maximum energy in the band plot.
            fig_name (str): Name of the figure to save.
            plot_kpath: List of k-point labels to plot. If provided, should be in the form like [["G", "H"], ["H", "K"]]. 
                If not provided, will use kpath stored in the object.
        """
        import matplotlib.pyplot as plt

        high_symm_poses, high_symm_labels = [], []

        if plot_kpaths is None:
            for ikpath, kpath in enumerate(self.kpaths):
                # Append high-symmetry point labels and positions to list
                start_label, end_label = self.kpaths[ikpath]["start"], self.kpaths[ikpath]["end"]
                start_nkpt, end_nkpt = self.kpaths[ikpath]["start_nkpt"], self.kpaths[ikpath]["end_nkpt"]
                if len(high_symm_labels) > 0 and high_symm_labels[-1] != start_label:
                    high_symm_labels[-1] += f"|{start_label}"
                    high_symm_labels.append(end_label)
                    high_symm_poses.append(self.kpath_lengths[end_nkpt])
                else:
                    high_symm_labels.extend([start_label, end_label])
                    high_symm_poses.extend([self.kpath_lengths[start_nkpt], self.kpath_lengths[end_nkpt]])
                
                if ikpath < len(self.kpaths) - 1 and self.kpaths[ikpath]["end"] == self.kpaths[ikpath+1]["start"]:
                    # Include first k-point of the next kpath for continous point in the plot to avoid discontinuties in the plot
                    start_nkpt, end_nkpt = kpath["start_nkpt"], kpath["end_nkpt"] + 1
                else:
                    start_nkpt, end_nkpt = kpath["start_nkpt"], kpath["end_nkpt"]
                
                for iband in range(self.band_data.shape[2]):
                    # nspin = 1, 4 and spin up channel of nspin = 2
                    plt.plot(
                        self.kpath_lengths[start_nkpt:end_nkpt],
                        self.band_data[0, start_nkpt:end_nkpt, iband],
                        "r-",
                        linewidth=1.0,
                    )
                    if self.band_data.shape[0] == 2:
                        plt.plot(
                            self.kpath_lengths[start_nkpt:end_nkpt],
                            self.band_data[1, start_nkpt:end_nkpt, iband],
                            "b--",
                            linewidth=1.0,
                        )
            
            plt.xlim(self.kpath_lengths[0], self.kpath_lengths[-1])

        else:
            start_x_pos = 0
            for plot_kpath in plot_kpaths:
                start_label, end_label = plot_kpath[0], plot_kpath[1]
                
                plot_band_data = None

                for kpath in self.kpaths:
                    if start_label == kpath["start"] and end_label == kpath["end"]:
                        plot_band_data = self.band_data[:, kpath["start_nkpt"]:kpath["end_nkpt"], :]
                        break
                    elif start_label == kpath["end"] and end_label == kpath["start"]:
                        plot_band_data = self.band_data[:, kpath["end_nkpt"]:kpath["start_nkpt"]:-1, :]
                        break
                
                if plot_band_data is None:
                    raise ValueError("Provided kpath not found in the band structure")
                    
                band_x_min = min(self.kpath_lengths[kpath["start_nkpt"]:kpath["end_nkpt"]])
                band_x_max = max(self.kpath_lengths[kpath["start_nkpt"]:kpath["end_nkpt"]])
                for iband in range(plot_band_data.shape[2]):
                    plt.plot(
                        self.kpath_lengths[kpath["start_nkpt"]:kpath["end_nkpt"]] - band_x_min + start_x_pos,
                        plot_band_data[0, :, iband],
                        "r-",
                        linewidth=1.0,
                    )
                    if plot_band_data.shape[0] == 2:
                        plt.plot(
                            self.kpath_lengths[kpath["start_nkpt"]:kpath["end_nkpt"]] - band_x_min + start_x_pos,
                            plot_band_data[1, :, iband],
                            "b--",
                            linewidth=1.0,
                        )
                
                # Append high-symmetry point labels and positions to list
                if len(high_symm_labels) > 0 and high_symm_labels[-1] != start_label:
                    high_symm_labels[-1] += f"|{start_label}"
                    high_symm_labels.append(end_label)
                    high_symm_poses.append(start_x_pos + band_x_max - band_x_min)
                else:
                    high_symm_labels.extend([start_label, end_label])
                    high_symm_poses.extend([start_x_pos, start_x_pos + band_x_max - band_x_min])

                start_x_pos += band_x_max - band_x_min
            
            plt.xlim(0, start_x_pos)
                    
        plt.ylim(emin, emax)
        plt.ylabel(r"$E-E_\text{F}$/eV")
        plt.xticks(high_symm_poses, high_symm_labels)
        plt.grid()
        plt.title(f"Band structure")
        plt.tight_layout()
        plt.savefig(fig_name, dpi=300)
        plt.close()

