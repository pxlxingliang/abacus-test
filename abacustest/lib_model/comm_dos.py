"""Module for managing DOS and PDOS data from ABACUS calculations."""

import numpy as np
from typing import Dict, List, Tuple, Optional, Union

from abacustest.lib_collectdata.collectdata import RESULT


l_map = ['s', 'p', 'd', 'f', 'g']

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

color_list = [
    '#1f77b4',  # 蓝色
    '#ff7f0e',  # 橙色
    '#2ca02c',  # 绿色
    '#d62728',  # 红色
    '#9467bd',  # 紫色
    '#8c564b',  # 棕色
    '#e377c2',  # 粉色
    '#7f7f7f',  # 灰色
    '#bcbd22',  # 黄绿色
    '#17becf',  # 青色
    '#aec7e8',  # 浅蓝色
    '#ffbb78',  # 浅橙色
    '#98df8a',  # 浅绿色
    '#ff9896',  # 浅红色
]

class DOSData:
    """Class for managing DOS data from ABACUS calculations."""
    def __init__(self,
                 energy: np.ndarray = None,
                 dosdata: np.ndarray = None,
                 efermi: float=None):
        """
        Args:
            energy: Original energy values for the DOS.
            dosdata: Total DOS data. If nspin=1 or 4, should be a 1D array with shape (nenergy,), or nspin=2 with shape (nenergy, 2).
            efermi: Fermi energy.
        """
        if efermi is not None:
            self.energy = energy - efermi
        else:
            self.energy = energy

        self.dosdata = dosdata
    
    @staticmethod
    def ReadFromAbacusJob(abacusjob_dir: str, efermi: Optional[float] = None) -> 'DOSData':
        """
        Read DOS outputed by ABACUS from directory of finished ABACUS calculation.
        """
        results = RESULT(fmt="abacus", path=abacusjob_dir)
        energy = np.array(results['dos']['energy'])
        dosdata = np.array(results['dos']['data'])
        if efermi is None:
            efermi = results['efermi']

        return DOSData(energy=energy, dosdata=dosdata, efermi=efermi)
    
    @staticmethod
    def plot_dos(dosdata: np.ndarray, energy: np.ndarray, emin: float, emax: float, title: str, fname: str):
        """
        Plot DOS using the given energy and dosdata.
        """
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 6))

        is_spin_polarized = False
        if len(dosdata.shape) == 1: #nspin = 1 or 4
            plt.plot(energy, dosdata, label='Total DOS', linestyle='-', linewidth=1)
            plt.ylim(bottom=0)
        elif len(dosdata.shape) == 2: #nspin = 2
            is_spin_polarized = True
            plt.plot(energy, dosdata[:, 0], label='Spin up', linestyle='-', linewidth=1)
            plt.plot(energy, -dosdata[:, 1], label='Spin down', linestyle='--', linewidth=1)
        plt.xlim(emin, emax)
        if not is_spin_polarized:
            plt.ylim(bottom=0)
        plt.xlabel('Energy (eV)', fontsize=12)
        plt.ylabel(r"States ($eV^{-1}$)", fontsize=12)
        plt.axvline(x=0, color="k", linestyle=":", alpha=0.5)
        plt.title(title, fontsize=12)
        plt.legend(loc='best')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(fname, dpi=300)
        plt.close()
    
    @staticmethod
    def write_dos(dosdata: np.ndarray, energy: np.ndarray, fname: str):
        """
        Write DOS data to a file.
        """
        with open(fname, "w") as dosfile:
            if len(dosdata.shape) == 1:
                dosfile.write(f"{'Energy (eV)':>14s}{'Total DOS':>14s}\n")
                for i in range(len(energy)):
                    dosfile.write(f"{energy[i]:>14.6f}{dosdata[i]:>14.6f}\n")
            else:
                dosfile.write(f"{'Energy (eV)':>14s}{'Spin up':>14s}{'Spin down':>14s}\n")
                for i in range(len(energy)):
                    dosfile.write(f"{energy[i]:>14.6f}{dosdata[i, 0]:>14.6f}{dosdata[i, 1]:>14.6f}\n")


class PDOSData:
    """Class for managing PDOS data from ABACUS calculations."""

    def __init__(self,
                 energy: np.ndarray = None, 
                 pdosdata: List[Dict[str, Tuple[int, np.ndarray]]] = None,
                 efermi: float=None):
        """
        Initialize DosObj to manage DOS and PDOS data from ABACUS output files.
        Parameters:
        -----------
        energy (np.ndarray):
            Original energy data outputed in DOS and PDOS calculation.
        pdosdata (List[Dict[str, Tuple[int, np.ndarray]]]):
            PDOS data for each orbital. Each element is a dictionary containing orbital information and PDOS data, with keys "atom_index", "l", "m", "z" and "data". The format of "data" is same with dosdata.
        efermi (float):
            Fermi energy. Original energy data is subtracted by this value. If None, will use the original energy data.
        """
        if efermi is not None:
            self.energy = np.array(energy) - efermi

        for orbital_pdos in pdosdata:
            assert len(energy) == orbital_pdos['data'].shape[0] # Length of energy and PDOS data of one orbital must match
            
        self.projected_dos = pdosdata

    def read_projected_dos(abacusjob_dir):
        """Read projected DOS from ABACUS output files (XML format)."""
        results = RESULT(fmt="abacus", path=abacusjob_dir)
        energy = np.array(results['pdos']['energy'])
        for i in range(len(results['pdos']['orbitals'])):
            results['pdos']['orbitals'][i]['data'] = np.array(results['pdos']['orbitals'][i]['data'])
        return energy, results['pdos']['orbitals']
    
    @staticmethod
    def ReadFromAbacusJob(abacus_job: str, efermi: Optional[float]=None) -> 'PDOSData':
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
        results = RESULT(fmt="abacus", path=abacus_job)
        energy, pdosdata = PDOSData.read_projected_dos(abacus_job)
        if efermi is None:
            efermi = results['efermi']

        return PDOSData(energy=energy, pdosdata=pdosdata, efermi=efermi)

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
    def get_pdos_by_species(self, species: str, sum_only: bool=True) -> Union[np.ndarray, List[Dict]]:
        """
        Get PDOS data for a specific species.
    
        Parameters:
        -----------
        species : str
            Species name (e.g. 'Fe', 'O')
        sum_only : bool, optional
            If True, sum the PDOS data for each orbital, otherwise return the original data. Default is True.
    
        Returns:
        --------
        list of dict or np.ndarray
            If sum_only is True, returns summed PDOS data array with shape (nspin, len(energy)).
            Otherwise, returns list of orbitals belonging to the specified species
        """
        pdos_datas = [orb for orb in self.projected_dos if orb['species'] == species]
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
    
    @staticmethod
    def plot_pdos(pdosdatas: List[List[np.ndarray]],
                  labels: List[List[str]],
                  titles: List[str],
                  energy: np.ndarray,
                  energy_min: float,
                  energy_max: float,
                  pdos_fig_name: str):
        """
        Plot PDOS using given PDOS data, titles and labels.

        Args:
            - pdosdata (List[List[np.ndarray]]): PDOS datas used in the plot. The length of PDOS data is the number of subplots in the plot, and
              the length of each element of the list is number of PDOS data plotted in each subplot. If nspin=2, the spin up and spin down PDOS
              will be plotted in a subplot separately.
            - labels (List[List[str]]): Labels for each PDOS data in each subplot. The length of labels is the number of subplots in the plot, and
              the length of each element of the list is number of PDOS data plotted in each subplot. If nspin=2, up and down arrows will be added automatically.
            - titles (List[str]): Titles for each subplot. The length of titles is the number of subplots in the plot.
            - energy (np.ndarray): Energy values for the PDOS data.
            - energy_min (float): Minimum energy value for the plot.
            - energy_max (float): Maximum energy value for the plot.
            - pdos_fig_name (str): Name of the PDOS figure.
        """
        import matplotlib.pyplot as plt

        num_subplots = len(pdosdatas)
        fig, axes = plt.subplots(num_subplots, 1, figsize=(8, 4*num_subplots))

        for idx, pdosdata in enumerate(pdosdatas):
            assert len(pdosdata) == len(labels[idx]) # Check if number of PDOS data matches number of labels

            if num_subplots > 1:
                ax = axes[idx]
            else:
                ax = axes
            
            if_spin_polarized = False
            for i, (data, label) in enumerate(zip(pdosdata, labels[idx])):
                if len(data.shape) == 1: # nspin = 1 or 4
                    ax.plot(energy, data, label=label, linestyle='-', color=color_list[i], linewidth=1.0)
                elif len(data.shape) == 2: # nspin = 2
                    if_spin_polarized = True
                    ax.plot(energy, data[:,0], label=f"{label} "+r"$\uparrow$", linestyle='-', color=color_list[i], linewidth=1.0)
                    ax.plot(energy, -data[:,1], label=f"{label} "+r"$\downarrow$", linestyle='--', color=color_list[i], linewidth=1.0)
            
            ax.set_xlim(energy_min, energy_max)
            ax.relim(visible_only=True)
            ax.autoscale_view()
            ax.set_xlabel("Energy (eV)", fontsize=12)
            ax.set_ylabel(r"States ($eV^{-1}$)", fontsize=12)
            ax.axvline(x=0, color="k", linestyle=":", alpha=0.5)
            ax.legend(loc='best', fontsize=8, ncol=if_spin_polarized+1)
            ax.grid(alpha=0.3)
            ax.set_title(titles[idx])
        
        plt.tight_layout()
        plt.savefig(pdos_fig_name, dpi=300)
        plt.close()
    
    @staticmethod
    def write_pdos(pdosdatas: List[np.ndarray], energy: np.ndarray, labels: List[str], filename: str):
        """
        Write processed PDOS data to a file.

        Args:
            - pdosdata (List[np.ndarray]): Processed PDOS data. The shape of each element is same with self.projected_dos[i]['data'].
            - energy (np.ndarray): Energy values for the PDOS data.
            - labels (List[str]): Labels for each PDOS data. If nspin=2. up and down PDOS will be written in different columns.
            - filename (str): Name of the file to write the data to.
        """
        assert len(pdosdatas) == len(labels)

        with open(filename, "w") as writefile:
            header_line = f"{'Energy (eV)':>14s}"
            for pdosdata, label in zip(pdosdatas, labels):
                if len(pdosdata.shape) == 1:
                    header_line += f"{label:>20s}"
                elif len(pdosdata.shape) == 2:
                    header_line += f"{label+'_up':>20s}{label+'_dn':>20s}"
            writefile.write(header_line + "\n")

            for i in range(len(energy)):
                line = f"{energy[i]:>14.6f}"
                for pdosdata in pdosdatas:
                    if len(pdosdata.shape) == 1:
                        line += f"{pdosdata[i]:>20.6f}"
                    elif len(pdosdata.shape) == 2:
                        line += f"{pdosdata[i,0]:>20.6f}{pdosdata[i,1]:>20.6f}"
                writefile.write(line + "\n")

