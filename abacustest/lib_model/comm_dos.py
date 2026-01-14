"""Module for managing DOS and PDOS data from ABACUS calculations."""

import os
from io import StringIO
import numpy as np
from typing import Dict, List, Tuple, Optional, Union
import xml.etree.ElementTree as ET

from abacustest.lib_prepare.abacus import ReadInput
from abacustest.lib_collectdata.collectdata import RESULT


l_map = ['s', 'p', 'd', 'f', 'g']

class PDOSData:
    """Class for managing DOS and PDOS data from ABACUS calculations."""

    def __init__(self, nspin: int,
                 energy: np.ndarray, 
                 dosdata: np.ndarray, 
                 pdosdata: List[Dict[str, Tuple[int, np.ndarray]]], 
                 orbital_names: Dict[Tuple[int, int, int, str], str], efermi: float=0.0):
        """
        Initialize DosObj to manage DOS and PDOS data from ABACUS output files.
        Parameters:
        -----------
        nspin: int
            Number of spin channels. Should be among [1, 2, 4].
        energy (np.ndarray):
            Original energy data outputed in DOS and PDOS calculation.
        dosdata (np.ndarray):
            Total DOS data. If nspin is 1 or 4, dosdata is a 1D array. If nspin is 2, dosdata is a 2D array with shape (nenergy, nspin).
        pdosdata (List[Dict[str, Tuple[int, np.ndarray]]]):
            PDOS data for each orbital. Each element is a dictionary containing orbital information and PDOS data, with keys "atom_index", "l", "m", "z" and "data". The format of "data" is same with dosdata.
        orbital_names (Dict[Tuple[int, int, int, str], str]):
            Dictionary mapping orbital indices to their names. Each key is a tuple containing atom index, l, m, and z values, and each value is the name of the atomic orbital.
        """
        assert len(energy) == dosdata.shape[0] # Length of energy and DOS data must match
        if len(dosdata.shape) == 1:
            assert nspin in [1, 4]
        elif len(dosdata.shape) == 2:
            assert nspin == 2
        else:
            raise ValueError("Invalid shape of DOS data")
        
        self.nspin = nspin
        self.energy = np.array(energy) - efermi
        self.dos = np.array(dosdata)

        for orbital_pdos in pdosdata:
            assert len(energy) == orbital_pdos['data'].shape[0] # Length of energy and PDOS data of one orbital must match
            if len(orbital_pdos['data'].shape) == 1:
                assert nspin in [1, 4]
            elif len(orbital_pdos['data'].shape) == 2:
                assert nspin == 2
            self.projected_dos = pdosdata

        for orbital_pdos in self.projected_dos:
            atom_idx = orbital_pdos['atom_index']
            l = orbital_pdos['l']
            m = orbital_pdos['m']
            z = orbital_pdos['z']
            key = (atom_idx, l, m, z)
            if key in orbital_names:
                orbital_pdos['atomic_orbital_name'] = orbital_names[key]
            else:
                orbital_pdos['atomic_orbital_name'] = None

    def read_total_dos(abacusjob_dir: str) -> Tuple[np.ndarray, np.ndarray]:
        """Read total DOS from ABACUS output files, with support for spin-polarized calculations."""
        results = RESULT(fmt="abacus", path=abacusjob_dir)
        return results['dos']['energy'], results['dos']['data']

    @staticmethod
    def read_orbital_info(orbital_file):
        """Read orbital information from the Orbital file to get symmetry information."""
        atom_orb_names = {}
        if os.path.isfile(orbital_file):
            try:
                with open(orbital_file, 'r') as f:
                    lines = f.readlines()
                
                # Parse the orbital file
                for line in lines:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 6:
                        try:
                            atom_idx = int(parts[0]) + 1 # Use atom index starting from 1
                            species = parts[1]
                            l = int(parts[2])
                            m = int(parts[3])
                            z = int(parts[4])
                            atom_orb_name = parts[5]
                            
                            # Store the mapping from (atom_idx, l, m, z) to symmetry name
                            atom_orb_names[(atom_idx, l, m, z)] = atom_orb_name
                        except (ValueError, IndexError) as e:
                            # Skip lines that don't match the expected format
                            continue
            
                return atom_orb_names
    
            except Exception as e:
                print(f"Failed to read orbital file {orbital_file}: {e}")
        else:
            print(f"Orbital file not found at {orbital_file}")
    
    def read_projected_dos(abacusjob_dir):
        """Read projected DOS from ABACUS output files (XML format)."""
        results = RESULT(fmt="abacus", path=abacusjob_dir)
        return results['pdos']['energy'], results['pdos']['orbitals']
    
    @staticmethod
    def ReadFromAbacusJob(abacus_job: str) -> 'PDOSData':
        """
        Create a PDOSData object from an ABACUS job object.
        
        Parameters:
        -----------
        abacus_job : AbacusJob
            An instance of the AbacusJob class representing an ABACUS calculation.
        
        Returns:
        --------
        PDOSData
            An instance of the PDOSData class with data read from the ABACUS job output.
        """
        input_params = ReadInput(os.path.join(abacus_job, "INPUT"))
        suffix = input_params.get('suffix', 'ABACUS')
        nspin = input_params.get('nspin', 1)
        energy, dos = PDOSData.read_total_dos(abacus_job)
        results = RESULT(fmt="abacus", path=abacus_job)
        efermi = results['efermi']
        energy, all_orbitals = PDOSData.read_projected_dos(abacus_job)

        atom_orb_names = PDOSData.read_orbital_info(os.path.join(abacus_job, f"OUT.{suffix}", "Orbital"))

        return PDOSData(nspin=nspin, energy=energy, dosdata=dos, pdosdata=all_orbitals, orbital_names=atom_orb_names, efermi=efermi)

    def get_pdos_by_orbital(self, species: str, atom_index: int, l: int, 
                            m: Optional[int] = None, z: Optional[int] = None) -> Dict:
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
            if (orbital['species'] == species and 
                orbital['atom_index'] == atom_index and 
                orbital['l'] == l):
                
                # If m and z are specified, match those as well
                if m is not None and orbital['m'] != m:
                    continue
                if z is not None and orbital['z'] != z:
                    continue
                
                return orbital
        
        return None
    
    def get_pdos_by_atom(self, atom_index: int, sum_only: bool=True) -> List[Dict]:
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
        pdos_datas = [orb for orb in self.projected_dos if orb['atom_index'] == atom_index]
        if sum_only:
            return PDOSData.sum_pdos_data(pdos_datas)
        else:
            return pdos_datas

    def get_pdos_by_atom_shell(self, atom_index: int, l: Tuple[int, str], sum_only: bool=True) -> List[Dict]:
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
        pdos_datas = [orb for orb in self.projected_dos if (orb['atom_index'] == atom_index and orb['l'] == l)]
        if sum_only:
            return PDOSData.sum_pdos_data(pdos_datas)
        else:
            return pdos_datas
    
    def get_pdos_by_atom_orbital(self, atom_index: int, l: Optional[Tuple[int, str]], m: Tuple[int, str], sum_only: bool=True) -> List[Dict]:
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
            pdos_datas = [orb for orb in self.projected_dos if (orb['atom_index'] == atom_index and orb['atomic_orbital_name'] == m)]
        elif isinstance(m, int):
            if l in l_map:
                l = l_map.index(l)
            pdos_datas = [orb for orb in self.projected_dos if (orb['atom_index'] == atom_index and orb['l'] == l and orb['m'] == m)]
        else:
            raise ValueError("Invalid m value")
        
        if sum_only:
            return PDOSData.sum_pdos_data(pdos_datas)
        else:
            return pdos_datas

    def get_pdos_by_species_shell(self, species: str, l: Union[int, str], sum_only: bool=True) -> Union[np.ndarray, List[Dict]]:
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
        pdos_datas = [orb for orb in self.projected_dos if (orb['species'] == species and orb['l'] == l)]

        if sum_only:
            return PDOSData.sum_pdos_data(pdos_datas)
        else:
            return pdos_datas

    def get_pdos_by_species_shell(self, species: str, l: Union[int, str], sum_only: bool=True) -> Union[np.ndarray, List[Dict]]:
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
        pdos_datas = [orb for orb in self.projected_dos if (orb['species'] == species and orb['l'] == l)]
    
        if sum_only:
            return PDOSData.sum_pdos_data(pdos_datas)
        else:
            return pdos_datas
    
    def get_pdos_by_species_orbital(self, species: str, l: Optional[Union[int, str]], m: Union[int, str], sum_only: bool=True) -> Union[np.ndarray, List[Dict]]:
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
            pdos_datas = [orb for orb in self.projected_dos if (orb['species'] == species and orb['atomic_orbital_name'] == m)]
        elif isinstance(m, int):
            if l in l_map:
                l = l_map.index(l)
            pdos_datas = [orb for orb in self.projected_dos if (orb['species'] == species and orb['l'] == l and orb['m'] == m)]
        else:
            raise ValueError("Invalid m value")
        
        if sum_only:
            return PDOSData.sum_pdos_data(pdos_datas)
        else:
            return pdos_datas

    @staticmethod
    def sum_pdos_data(pdos_datas: List[Dict]) -> np.ndarray:
        """
        Sum the PDOS data from a list of extracted PDOS data.
    
        Parameters:
        -----------
        pdos_datas : list of dict
            List of PDOS data dictionaries
    
        Returns:
        --------
            Summed PDOS will be an array with shape (nspin, len(energy))
        """
        if len(pdos_datas) == 0:
            raise ValueError("No PDOS data provided")
    
        pdos_sum = np.zeros_like(pdos_datas[0]['data'])
        for pdos_data in pdos_datas:
            pdos_sum += np.array(pdos_data['data'])
        return pdos_sum
