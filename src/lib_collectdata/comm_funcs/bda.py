from abc import ABC, abstractmethod
from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.io.vasp.outputs import Oszicar, Outcar, Vasprun
from pymatgen.core.periodic_table import Element
from pathlib import Path
import pandas as pd
from typing import List, Dict, Literal, Union

class BasicProperty(ABC):
    def __init__(self, magnetization: tuple = None, 
                 energy:float = None, 
                 volume: float=None, 
                 structrue: Poscar = None):
        self._magnetization = magnetization
        self.energy = energy
        self.volume = volume
        self.structure = structrue

    @abstractmethod
    def bond_length(
            self, 
        ) -> pd.DataFrame:
        r"""Calculated bond length of all transition metal atoms to their nearest neighbors 

        Returns
        -------
        - dataframe : (`pandas.DataFrame`) The bond length of all transition metal atoms 
                to their nearest neighbors
        """
        pass

    @abstractmethod
    def magnetic_moment(
            self, 
            element: Literal['TM', 'O', 'ALL']='ALL'
        ) -> pd.DataFrame:
        r"""Grasp the magnetic moment of elements interested
        Parameters
        ----------
        - element : (`str`) interested element, now support 'TM', 'O'(oxygen) and 'ALL', stands for exporting
                magnetic moment of 'Transition Metal', 'Oxygen' and 'ALL element'

        Returns
        -------
        - mag : (`pandas.DataFrame`) magnetization of interested element
        """
        pass

    @abstractmethod
    def inter_layer_distance(
        self,
        ) -> dict:
        r"""Calculate inter layer distance of both TM-layer and alkali-layer in layered structure

        Returns
        -------
        - argument : (`Dict[float]`) Dictionary of averaged inter layer distance
                    e.g. {"TM-layer": 2.5, "Alkali-layer": 2.1}
        """
        pass

    def __cut(self, sortedList: List[float], num: float) -> List[List[float]]:
        r"""Cut a sorted list into List[List] by the distance of its elements
            e.g. [1,1,1,3,3,7] -> [[1,1,1], [3,3], [7]] if num <= 1
            e.g. [1,1,1,3,3,7] -> [[1,1,1,3,3],[7]] if 2 < num <= 4

        Parameters
        ----------
        - sortedList : (`List`) A sorted list
        - num : (`List`) A float decide how the sorted list devided

        Returns
        -------
        - argument : (`List[List[float]]`)
        """
        ptr = 0
        list = []
        for i in range(0, len(sortedList) - 1):
            if (sortedList[i+1] - sortedList[i]) > num:
                list.append(sortedList[ptr:i+1])
                ptr = i + 1
        list.append(sortedList[ptr::])
        return list
    
    def __get_num_index(self) -> dict:
        r"""Get the number of atoms, index in structure file of all elements in structure

        Returns
        -------
        -dict : (`dict`) Contains element symbol, start index and number of atoms in structure file
                e.g. {"Fe": {"start_index": 0, "number": 20}, 
                      "O": {"start_index": 20, "number":10}, 
                      "F": {"start_index": 30, "number": 10}}
        """
        num_index = {}
        site_symbols = self.structure.site_symbols
        natoms = self.structure.natoms
        index_cnt = 1
        for _id, ele in enumerate(site_symbols):
            num_index[ele] = {"start_index": index_cnt - 1, "number": natoms[_id]}
            index_cnt += natoms[_id]
        return num_index

    def __is_transition(self, ele: str) -> bool:
        TM = ('Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pb', 'Ag', 'Cd')
        return (ele in TM)

    def bond_length(
            self 
        ) -> pd.DataFrame:
        r"""Calculated bond length of all transition metal atoms to their nearest neighbors 

        Returns
        -------
        - dataframe : (`pandas.DataFrame`) The bond length of all transition metal atoms 
                to their nearest neighbors
        """
        structure = self.structure.structure
        num_index = self.__get_num_index()
        TM_bonds = []
        for ele in num_index.keys():
            if self.__is_transition(ele):
                start_index = num_index[ele]["start_index"]
                number = num_index[ele]["number"]
                TM_sites = structure.sites[start_index:start_index + number]
                neighbor = structure.get_all_neighbors(r=2.5, sites=TM_sites)
                for ii in range(len(TM_sites)):
                    temp_dict = {"Total_index": start_index + ii + 1, "Atomic_index": ii + 1, "Element": ele}
                    total_length = 0
                    for jj in range(len(neighbor[ii])):
                        length = TM_sites[ii].distance(neighbor[ii][jj])
                        temp_dict["bond_length_{}".format(jj + 1)] = length
                        total_length += length
                    temp_dict["bond_length_ave"] = total_length/(len(neighbor[ii]))
                    TM_bonds.append(temp_dict)
        #return pd.DataFrame.from_dict(TM_bonds)
        return TM_bonds
          
    def magnetic_moment(
            self, 
            element: Literal['TM', 'O', 'ALL']='ALL'
        ) -> list:
        r"""Grasp the magnetic moment of elements interested
        Parameters
        ----------
        - element : (`str`) interested element, now support 'TM', 'O'(oxygen) and 'ALL', stands for exporting
                magnetic moment of 'Transition Metal', 'Oxygen' and 'ALL element'

        Returns
        -------
        - mag : (`pandas.DataFrame`) magnetization of interested element
        """
        if element.upper() not in ('TM', 'O', 'ALL'):
            raise ValueError("Invalid element")
        magnetization = self._magnetization
        num_index = self.__get_num_index()
        for ele in num_index.keys():
            ele_num = num_index[ele].get("number")
            for num in range(ele_num):
                real_index = num_index[ele]["start_index"] + num
                magnetization[real_index]["Total_index"] = real_index + 1
                magnetization[real_index]["Element"] = ele
                magnetization[real_index]["Atomic_index"] = num + 1
        if element == 'ALL':
            return list(magnetization)
            #return pd.DataFrame.from_dict(list(magnetization))
        elif element == 'TM':
            mag = ()
            for ele in num_index.keys():
                if self.__is_transition(ele):
                    start_index = num_index[ele].get("start_index")
                    number = num_index[ele].get("number")
                    mag += magnetization[start_index:start_index + number]
            #return pd.DataFrame.from_dict(list(mag))
            return list(mag)
        else:
            start_index = num_index['O'].get("start_index")
            number = num_index['O'].get("number")
            mag = magnetization[start_index:start_index + number]
            #return pd.DataFrame.from_dict(list(mag))
            return list(mag)

    def inter_layer_distance(self) -> dict:
        from statistics import mean
        r"""Calculated inter-layer distance of layered structure like LiCoO2
        Parameters
        ----------

        Returns
        -------
        - dict : (`dict`) Dictionary contains inter-layer distance of TM-layers and alkali-layers
        """
        cart_coord = self.structure.structure.cart_coords
        num_index = self.__get_num_index()
        start_index = num_index['O'].get("start_index")
        number = num_index['O'].get("number")
        oxy_coord = cart_coord[start_index:start_index+number]
        oxy_coord = sorted(oxy_coord[:,-1])
        cluster_oxy_coord = self.__cut(oxy_coord, 1)
        height_position = []
        for ii in cluster_oxy_coord:
            height_position.append(mean(ii))
        distance = []
        for i in range(1, len(height_position)):
            distance.append(height_position[i] - height_position[i - 1])
        return {"d_TM": min(distance[0], distance[1]), "d_alkali": max(distance[0], distance[1])}

class BasicPropertyVasp(BasicProperty):
    def __init__(self, out_file, struc_file):  
        from pymatgen.io.vasp.outputs import Outcar
        from pymatgen.io.vasp.inputs import Poscar
        self.energy = Outcar(out_file).final_energy
        self.volume = Poscar.from_file(struc_file).structure.volume
        self._magnetization = Outcar(out_file).magnetization
        self.structure = Poscar.from_file(struc_file)
    
if __name__ == "__main__":
    output = BasicPropertyVasp('OUTCAR', 'CONTCAR')
    energy = output.energy
    volume = output.volume
    mag_moment = output.magnetic_moment('TM')
    bond_length = output.bond_length()
    print(energy)
    print(volume)
    print(mag_moment)
    print(bond_length)
