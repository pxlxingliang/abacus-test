"""Module for managing DOS and PDOS data from ABACUS calculations."""

import os
from io import StringIO
import numpy as np
from typing import Dict, List, Tuple, Optional, Union
import xml.etree.ElementTree as ET

from abacustest.lib_prepare.abacus import ReadInput
from abacustest.lib_collectdata.collectdata import RESULT


l_map = ["s", "p", "d", "f", "g"]


class PDOSData:
    """Class for managing DOS and PDOS data from ABACUS calculations."""

    def __init__(self, path: str = "./"):
        """
        Initialize DosObj to manage DOS and PDOS data from ABACUS output files.
        Parameters:
        -----------
        path : str
            Path to the directory containing ABACUS output files
        suffix : str, optional
            Suffix of the ABACUS output directory. If None, will auto-detect.
        """
        self._path = path
        self._suffix = self._get_suffix()
        self._nspin = self._get_nspin()
        self._fermi_energy = self._get_fermi_energy()
        self._dos: Optional[np.ndarray] = None  # Energy and total DOS [E, TDOS]
        self.projected_dos: List[Dict] = []  # List of projected DOS data
        self._energy: Optional[np.ndarray] = None  # Energy grid
        self.atom_orb_names: Dict[
            Tuple[int, int, int, int], str
        ] = {}  # Maps (atom_idx, l, m, z) to atomic orbital name

        # Try to read from different DOS file formats
        self._read_total_dos()
        self._read_orbital_info()  # Read orbital information from Orbital file
        self._read_projected_dos()  # Read projected DOS data

    def _get_suffix(self) -> str:
        input_file = os.path.join(self._path, "INPUT")
        input_params = ReadInput(input_file)
        return input_params.get("suffix", "ABACUS")

    def _get_nspin(self) -> int:
        input_file = os.path.join(self._path, "INPUT")
        input_params = ReadInput(input_file)
        return int(input_params.get("nspin", 1))

    def _get_fermi_energy(self) -> float:
        """Get the Fermi energy from the ABACUS output files."""
        results = RESULT(fmt="abacus", path=self._path)
        efermi = results["efermi"]
        if efermi is None:
            raise ValueError("Fermi energy not found in ABACUS output files")
        return float(efermi)

    def _read_total_dos(self):
        """Read total DOS from ABACUS output files, with support for spin-polarized calculations."""

        dos1_file = os.path.join(self._path, f"OUT.{self._suffix}", "DOS1_smearing.dat")
        dos1_data = np.loadtxt(dos1_file)
        self._energy = dos1_data[:, 0] - self._fermi_energy
        dos_data = [dos1_data[:, 1]]
        if self._nspin == 2:
            dos2_file = os.path.join(
                self._path, f"OUT.{self._suffix}", "DOS2_smearing.dat"
            )
            dos2_data = np.loadtxt(dos2_file)
            dos_data.append(dos2_data[:, 1])
        self._dos = np.array(dos_data)

    def _read_orbital_info(self):
        """Read orbital information from the Orbital file to get symmetry information."""
        # Look for the Orbital file in the output directory
        orbital_file = os.path.join(self._path, f"OUT.{self._suffix}", "Orbital")

        if os.path.isfile(orbital_file):
            try:
                with open(orbital_file, "r") as f:
                    lines = f.readlines()

                # Parse the orbital file
                for line in lines:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue

                    parts = line.split()
                    if len(parts) >= 6:
                        try:
                            atom_idx = (
                                int(parts[0]) + 1
                            )  # Use atom index starting from 1
                            species = parts[1]
                            l = int(parts[2])
                            m = int(parts[3])
                            z = int(parts[4])
                            atom_orb_name = parts[5]

                            # Store the mapping from (atom_idx, l, m, z) to symmetry name
                            self.atom_orb_names[(atom_idx, l, m, z)] = atom_orb_name
                        except (ValueError, IndexError) as e:
                            # Skip lines that don't match the expected format
                            continue

            except Exception as e:
                print(f"Warning: Failed to read orbital file {orbital_file}: {e}")
        else:
            print(f"Info: Orbital file not found at {orbital_file}")

    def _read_projected_dos(self):
        """Read projected DOS from ABACUS output files (XML format)."""
        pdos_file = os.path.join(self._path, f"OUT.{self._suffix}", "PDOS")

        if os.path.isfile(pdos_file):
            try:
                tree = ET.parse(pdos_file)
                root = tree.getroot()

                nspin_elem = root.find("nspin")
                if nspin_elem is None or nspin_elem.text is None:
                    raise ValueError("nspin element not found or empty in PDOS file")
                self.nspin = int(nspin_elem.text)

                energy_elem = root.find("energy_values")
                if energy_elem is None or energy_elem.text is None:
                    raise ValueError(
                        "energy_values element not found or empty in PDOS file"
                    )
                energy = [float(i) for i in energy_elem.text.split()]
                self._energy = np.array(energy) - self._fermi_energy

                all_orbitals = []
                for iorb in root.findall("orbital"):
                    data_elem = iorb.find("data")
                    if data_elem is None or data_elem.text is None:
                        print("Warning: Skipping orbital with missing data element")
                        continue
                    pdos_data = np.loadtxt(StringIO(data_elem.text))

                    # Get the orbital index, atom index, l, m, and z values
                    index_attr = iorb.get("index")
                    atom_index_attr = iorb.get("atom_index")
                    l_attr = iorb.get("l")
                    m_attr = iorb.get("m")
                    z_attr = iorb.get("z")
                    species = iorb.get("species")

                    if (
                        index_attr is None
                        or atom_index_attr is None
                        or l_attr is None
                        or m_attr is None
                        or z_attr is None
                        or species is None
                    ):
                        print(
                            f"Warning: Skipping orbital with missing attributes: {iorb.attrib}"
                        )
                        continue

                    index = int(index_attr)
                    atom_index = int(atom_index_attr)
                    l = int(l_attr)
                    m = int(m_attr)
                    z = int(z_attr)

                    # Look up the atomic orbital name from the orbital file
                    orbital_name_key = (atom_index, l, m, z)
                    atomic_orbital_name = self.atom_orb_names[orbital_name_key]

                    orbital_info = {
                        "index": index,
                        "atom_index": atom_index,
                        "species": species,
                        "l": l,
                        "m": m,
                        "z": z,
                        "atomic_orbital_name": atomic_orbital_name,  # Store the atomic orbital name
                        "data": pdos_data,
                    }

                    all_orbitals.append(orbital_info)

                self.projected_dos = all_orbitals
                print(f"Info: Successfully read PDOS from {pdos_file}")
            except Exception as e:
                print(f"Error: Failed to read PDOS from {pdos_file}: {e}")
                import traceback

                traceback.print_exc()
        else:
            print(f"Info: PDOS file not found at {pdos_file}")

    def get_pdos_by_orbital(
        self,
        species: str,
        atom_index: int,
        l: int,
        m: Optional[int] = None,
        z: Optional[int] = None,
    ) -> Optional[Dict]:
        """
        Get PDOS data for a specific orbital.

        Parameters:
        -----------
        species : str
            Species name (e.g. 'Fe', 'O')
        atom_index : int
            Index of the atom (0-based)
        l : int
            Angular momentum quantum number (0=s, 1=p, 2=d, 3=f, ...)
        m : int, optional
            Magnetic quantum number
        z : int, optional
            Orbital split index

        Returns:
        --------
        dict or None
            Dictionary containing orbital information and PDOS data
        """
        for orbital in self.projected_dos:
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

    def get_pdos_by_atom(
        self, atom_index: int, sum_only: bool = True
    ) -> Union[np.ndarray, List[Dict]]:
        """
        Get all PDOS data for a specific atom.

        Parameters:
        -----------
        atom_index : int
            Index of the atom (0-based)
        sum_only : bool, optional
            If True, sum the PDOS data for each orbital, otherwise return the original data. Default is True.

        Returns:
        --------
        list of dict
            List of orbitals belonging to the specified atom
        """
        pdos_datas = [
            orb for orb in self.projected_dos if orb["atom_index"] == atom_index
        ]
        if sum_only:
            return sum_pdos_data(pdos_datas)
        else:
            return pdos_datas

    def get_pdos_by_atom_shell(
        self, atom_index: int, l: Union[int, str], sum_only: bool = True
    ) -> Union[np.ndarray, List[Dict]]:
        """
        Get PDOS data for a specific atom shell.

        Parameters:
        -----------
        atom_index : int
            Index of the atom (0-based)
        l : int or str
            Angular momentum quantum number (0 or s, 1 or p, 2 or d, 3 or f, ...)
        sum_only : bool, optional
            If True, sum the PDOS data for each orbital, otherwise return the original data. Default is True.

        Returns:
        --------
        list of dict
            Dictionary containing orbital information and PDOS data
        """
        if l in l_map:
            l = l_map.index(l)
        pdos_datas = [
            orb
            for orb in self.projected_dos
            if (orb["atom_index"] == atom_index and orb["l"] == l)
        ]
        if sum_only:
            return sum_pdos_data(pdos_datas)
        else:
            return pdos_datas

    def get_pdos_by_atom_orbital(
        self,
        atom_index: int,
        l: Optional[Union[int, str]],
        m: Union[int, str],
        sum_only: bool = True,
    ) -> Union[np.ndarray, List[Dict]]:
        """
        Get PDOS data for a specific atom orbital.

        Parameters:
        -----------
        atom_index : int
            Index of the atom (0-based)
        l : int or str
            Angular momentum quantum number (0 or s, 1 or p, 2 or d, 3 or f, ...). If m is explicitly specified by atomic orbital name, l is ignored.
        m : int or str
            Magnetic quantum number used by PDOS file in ABACUS.
            If m is an integer, it should be between 0 an 2l.
            If m is a string, it should be the name of the atomic orbital.
        sum_only : bool, optional
            If True, sum all the PDOS data for the specified atom orbital, otherwise return all the individual PDOS data.

        Returns:
        --------
        dict or None
            Dictionary containing orbital information and PDOS data
        """
        if isinstance(m, str):
            pdos_datas = [
                orb
                for orb in self.projected_dos
                if (orb["atom_index"] == atom_index and orb["atomic_orbital_name"] == m)
            ]
        elif isinstance(m, int):
            if l in l_map:
                l = l_map.index(l)
            pdos_datas = [
                orb
                for orb in self.projected_dos
                if (orb["atom_index"] == atom_index and orb["l"] == l and orb["m"] == m)
            ]
        else:
            raise ValueError("Invalid m value")

        if sum_only:
            return sum_pdos_data(pdos_datas)
        else:
            return pdos_datas

    def get_pdos_by_species_shell(
        self, species: str, l: Union[int, str], sum_only: bool = True
    ) -> Union[np.ndarray, List[Dict]]:
        """
        Get PDOS data for shell of a specific species.

        Parameters:
        -----------
        species : str
            Species name (e.g. 'Fe', 'O')
        l : int or str
            Angular momentum quantum number (0 or 's', 1 or 'p', 2 or 'd', 3 or 'f', ...)
        sum_only : bool, optional
            If True, sum the PDOS data for each orbital, otherwise return the original data. Default is True.

        Returns:
        --------
        list of dict or np.ndarray
            If sum_only is True, returns summed PDOS data array with shape (nspin, len(energy))
            Otherwise, returns list of orbitals belonging to the specified species and shell
        """
        if l in l_map:
            l = l_map.index(l)
        pdos_datas = [
            orb
            for orb in self.projected_dos
            if (orb["species"] == species and orb["l"] == l)
        ]

        if sum_only:
            return sum_pdos_data(pdos_datas)
        else:
            return pdos_datas

    def get_pdos_by_species_orbital(
        self,
        species: str,
        l: Optional[Union[int, str]],
        m: Union[int, str],
        sum_only: bool = True,
    ) -> Union[np.ndarray, List[Dict]]:
        """
        Get PDOS data for a specific species orbital.

        Parameters:
        -----------
        species : str
            Species name (e.g. 'Fe', 'O')
        l : int or str, optional
            Angular momentum quantum number (0 or 's', 1 or 'p', 2 or 'd', 3 or 'f', ...).
            If m is explicitly specified by atomic orbital name, l is ignored.
        m : int or str
            Magnetic quantum number used by PDOS file in ABACUS.
            If m is an integer, it should be between 0 and 2l.
            If m is a string, it should be the name of the atomic orbital (e.g. 's', 'px', 'py', 'pz', 'dz2').
        sum_only : bool, optional
            If True, sum all the PDOS data for the specified species orbital, otherwise return all the individual PDOS data.

        Returns:
        --------
        list of dict or np.ndarray
            If sum_only is True, returns summed PDOS data array with shape (nspin, len(energy))
            Otherwise, returns list of orbitals belonging to the specified species and orbital
        """
        if isinstance(m, str):
            pdos_datas = [
                orb
                for orb in self.projected_dos
                if (orb["species"] == species and orb["atomic_orbital_name"] == m)
            ]
        elif isinstance(m, int):
            if l in l_map:
                l = l_map.index(l)
            pdos_datas = [
                orb
                for orb in self.projected_dos
                if (orb["species"] == species and orb["l"] == l and orb["m"] == m)
            ]
        else:
            raise ValueError("Invalid m value")

        if sum_only:
            return sum_pdos_data(pdos_datas)
        else:
            return pdos_datas

    @property
    def dos(self) -> np.ndarray:
        """Get the total DOS data."""
        if self._dos is None:
            raise ValueError("DOS data not loaded. Call _read_total_dos() first.")
        return self._dos

    @property
    def energy(self) -> np.ndarray:
        """Get the energy grid."""
        if self._energy is None:
            raise ValueError(
                "Energy data not loaded. Call _read_total_dos() or _read_projected_dos() first."
            )
        return self._energy


def sum_pdos_data(pdos_datas: List[Dict]) -> np.ndarray:
    """
    Sum the PDOS data from a list of extracted PDOS data.

    Parameters:
    -----------
    pdos_datas : list of dict
        List of PDOS data dictionaries. Each dictionary must contain a 'data' key
        with numpy array containing PDOS values.

    Returns:
    --------
    np.ndarray
        Summed PDOS array with shape (nspin, len(energy)). The shape depends on
        the spin configuration (nspin=1 or nspin=2).

    Raises:
    -------
    ValueError
        If pdos_datas is empty or if the 'data' key is missing from any dictionary.
    """
    if not pdos_datas:
        raise ValueError("pdos_datas list is empty")

    first_data = pdos_datas[0].get("data")
    if first_data is None:
        raise ValueError("First PDOS data dictionary missing 'data' key")

    pdos_sum = np.zeros_like(first_data)
    for pdos_data in pdos_datas:
        data = pdos_data.get("data")
        if data is None:
            raise ValueError(f"PDOS data dictionary missing 'data' key: {pdos_data}")
        pdos_sum += np.array(data)
    return pdos_sum
