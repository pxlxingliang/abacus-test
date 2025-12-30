from pydantic import BaseModel, Field
from typing import List, Tuple, Union, Dict, Any, Optional, Literal
import numpy as np
import copy

from abacustest.constant import MASS_DICT, A2BOHR, BOHR2A, ABACUS_STRU_KEY_WORD
from abacustest.lib_prepare.comm import Cartesian2Direct, Direct2Cartesian

import traceback
import os
import sys

class AbacusAtomType(BaseModel):
    """ABACUS Atom Type class, which defines the common properties for a type of atoms.

    Args:
        label (str): Atom label.
        element (str): Element symbol.
        mass (float): Atomic mass.
        pp (str, optional): Pseudopotential file name.
        orb (str, optional): Orbital file name.
        paw (str, optional): PAW file name.
        type_mag (float, optional): Default magnetic moment for the atom type. Default is 0.0.
        natom (int, optional): Number of atoms of this type. Default is 0.
    """
    label: str
    element: Optional[str] = None
    mass: Optional[float] = None
    pp: Optional[str] = None
    orb: Optional[str] = None
    paw: Optional[str] = None
    type_mag: float = 0.0
    natom: int = 0

    # define equivalence method
    def __eq__(self, other):
        if not isinstance(other, AbacusAtomType):
            return NotImplemented
        mass1 = round(self.mass, 6) if self.mass is not None else None
        mass2 = round(other.mass, 6) if other.mass is not None else None
        mag1 = round(self.type_mag, 6)
        mag2 = round(other.type_mag, 6)
        return (self.label == other.label and
                self.element == other.element and
                mass1 == mass2 and
                self.pp == other.pp and
                self.orb == other.orb and
                self.paw == other.paw and
                mag1 == mag2)

    def _default_list(self):
        out = [self.label]
        out.append(self.element if self.element is not None else "")
        out.append(self.mass if self.mass is not None else 0.0)
        out.append(self.pp if self.pp is not None else "")
        out.append(self.orb if self.orb is not None else "")
        out.append(self.paw if self.paw is not None else "")
        out.append(self.type_mag)
        return out

    def __lt__(self, other):
        return self._default_list() < other._default_list()

class AbacusATOM(BaseModel):
    """ABACUS ATOM class.

    Args:
        label (str): Atom label.
        coord (Tuple[float, float, float] or np.ndarray): Atom coordinate in cartesian type with unit Angstrom.
        element (str, optional): Element symbol. If not provided, it will be inferred from the label.
        mass (float, optional): Atomic mass. If not provided, it will be set according to the element.
        pp (str, optional): Pseudopotential file name.
        orb (str, optional): Orbital file name.
        paw (str, optional): PAW file name.
        type_mag (float, optional): Default magnetic moment for the atom type. Default is 0.0.
        move (Tuple[bool, bool, bool], optional): Movement constraints along x, y, z directions. Default is (True, True, True).
        mag (float or Tuple[float, float, float], optional): Magnetic moment for the atom. Can be a scalar or a 3D vector.
        angle1 (float, optional): For non-collinear magnetism, the magnetic moment direction with respect to z-axis.
        angle2 (float, optional): For non-collinear magnetism, the magnetic moment direction projection on xy-plane with respect to x-axis.
        velocity (Tuple[float, float, float], optional): Velocity of the atom.
        constrain (bool or Tuple[bool, bool, bool], optional): If do delta_spin constrain on this atom.
        lambda_ (float or Tuple[float, float, float], optional): Lambda parameter for the atom. Can be a scalar or a 3D tuple.

    NOTE: 
        1. type_mag is referred as the value defined in STRU file for each atom type, such as:
            O   # atom label
            0.0  # magnetic moment for this atom type
            2    # atom number
            0 0 0  # coord of atom 1
            0.5 0.5 0.5  # coord of atom 2
          If mag is defined for the atom, then type_mag will be ignored.
        2. mag/angle1/angle2 are the values defined for each atom individually, such as below:
            0 0 0 mag 2 angle1 90 angle2 45  # for atom 1
            0.5 0.5 0.5 mag 0 0 2  # for atom 2
        3. for non-collinear case, ABACUS supports to set mag and angle1/angle2 at the same time, and the real magnetic moment should calculated as:
            mag_z = |mag| * cos(angle1)
            mag_x = |mag| * sin(angle1) * cos(angle2)
            mag_y = |mag| * sin(angle1) * sin(angle2)
            When only one of angle1/angle2 is provided, the another angle will be set the default value 0.0.
        4. When read/write STRU file, the type_mag/mag/angle1/angle2 is exactly the same as in the file. None value means not set in the file.
        5. Use atommag to get the real magnetic moment vector (mag_x, mag_y, mag_z) for noncolinear case, and mag for collinear case. If mag is None will return type_mag.
        6. Use set_atommag() method to set the mag/angle1/angle2 values for the atom.
    """
    label: str
    coord: Union[Tuple[float, float, float]]
    # optional
    element: Optional[str] = None
    mass: Optional[float] = None
    pp: Optional[str] = None
    orb: Optional[str] = None
    paw: Optional[str] = None
    type_mag: Optional[float] = 0.0
    move: Optional[Tuple[bool, bool, bool]] = (True, True, True)
    velocity: Optional[Tuple[float, float, float]] = None
    constrain: Optional[Union[bool, Tuple[bool, bool, bool]]] = None
    lambda_: Optional[Union[float, Tuple[float, float, float]]] = None

    mag: Optional[Union[float, Tuple[float, float, float]]] = Field(default=None, frozen=True)  # magnitude of magnetic moment
    angle1: Optional[float] = Field(default=None, frozen=True)  # angle between magnetic moment and z-axis
    angle2: Optional[float] = Field(default=None, frozen=True)  # angle between projection of magnetic moment on xy-plane and x-axis

    def __init__(self, **data):
        super().__init__(**data)

        if self.element is None:
            self.element = self._infer_element_from_label()
        if self.mass is None:
            self.mass = MASS_DICT.get(self.element, 0.0)

    def _infer_element_from_label(self) -> Optional[str]:
        """Infer the element symbol from the atom label."""
        if len(self.label) >= 2 and self.label[:2].capitalize() in MASS_DICT:
            return self.label[:2].capitalize()
        elif self.label[0].capitalize() in MASS_DICT:
            return self.label[0].capitalize()
        return None

    def __str__(self):
        c = f"Label: {self.label}, Element: {self.element}, Mass: {self.mass:.4f}\n"
        c += f"PP: {self.pp}\nORB: {self.orb}\n"
        c += f"Type mag: {self.type_mag}, mag: {self.mag}, angle1: {self.angle1}, angle2: {self.angle2}\n"
        c += f"Coord: {self.coord[0]:%12.7f} {self.coord[1]:%12.7f} {self.coord[2]:%12.7f}\n"
        c += f"Move: {self.move}, Velocity: {self.velocity}\n"
        c += f"Constrain: {self.constrain}, Lambda: {self.lambda_}\n"
        c += f"Real atomic mag/magnitude: {self.atommag} / {self.atommag_magnitude}\n"
        return c

    @property
    def atomtype(self) -> AbacusAtomType:
        """Get the AbacusAtomType object for this atom."""
        return AbacusAtomType(
            label=self.label,
            element=self.element,
            mass=self.mass,
            pp=self.pp,
            orb=self.orb,
            paw=self.paw,
            type_mag=self.type_mag,
            natom=1
        )

    @property
    def noncolinear(self) -> bool:
        """Check if the atom is non-collinear magnetic. Only when angle1 or angle2 is set, or mag is a 3D vector, the atom is non-collinear magnetic."""
        
        if self.angle1 is not None or self.angle2 is not None:
            return True
        if isinstance(self.mag, (list, tuple)) and len(self.mag) == 3:
            return True
        return False

    @property
    def atommag_magnitude(self) -> float:
        """Get the magnitude of the magnetic moment.

        Returns:
            float: Magnitude of the magnetic moment.
        """
        import numpy as np
        if self.mag is None:
            return abs(self.type_mag)
        else:
            return np.linalg.norm(self.mag)

    @property
    def atommag(self) -> Union[float, Tuple[float, float, float]]:
        """Get the atomic magnetic moment: one float for colinear case, or a tuple of three floats for non-colinear case.
        """
        import numpy as np
        if self.noncolinear:
            mag_norm = self.atommag_moment
            angle1 = angle2 = 0.0
            if self.angle1 is not None:
                angle1 = self.angle1    
            if self.angle2 is not None:
                angle2 = self.angle2
            mag_z = mag_norm * np.cos(np.radians(angle1))
            mag_x = mag_norm * np.sin(np.radians(angle1)) * np.cos(np.radians(angle2))
            mag_y = mag_norm * np.sin(np.radians(angle1)) * np.sin(np.radians(angle2))
            return (mag_x, mag_y, mag_z)
        else:
            if self.mag is None:
                return self.type_mag
            else:
                return self.mag

    def set_atommag(self, mag: Union[float,Tuple[float,float,float]],
                    angle1: float = None,
                    angle2: float = None):
        """Set magnetic moment for the atom.

        Args:
            mag (float or Tuple[float,float,float]): Magnetic moment value.
            angle1 (float, optional): The angle between the magnetic moment and the z-axis in degrees.
            angle2 (float, optional): The angle between the projection of the magnetic moment on the xy-plane and the x-axis in degrees.
        
        If angle1 or angle2 is None, then this value will not write to STRU file.
        """
        if isinstance(mag, (list, tuple)):
            assert len(mag) == 3, f"Magnetic moment tuple must have three components, got {mag}."

        self.mag = mag
        self.angle1 = angle1
        self.angle2 = angle2

    @staticmethod
    def sort(atom_list: List["AbacusATOM"],
             keep_first_order=True,
             only_label=False) -> List["AbacusATOM"]:
        """Rearrange atom list according to atomtype.

        Args:
            atom_list (List[AbacusATOM]): List of AbacusATOM objects.
            keep_first_order (bool): If True, keep the order of first appearance of atom types. Default is True.
            only_label (bool): If True, only sort according to label, ignoring other properties (element/mass/pp/orb/paw/type_mag).

        Returns:
            List[AbacusATOM]: Rearranged list of AbacusATOM objects.
            List[int]: Indices of the rearranged atoms in the original list.
        """
        key = "label" if only_label else "atomtype"
        unique_atom = []
        for i in atom_list:
            attr = getattr(i, key)
            if attr not in unique_atom:
                unique_atom.append(attr)

        if not keep_first_order:
            unique_atom.sort()
        new_atom_list = []
        index = []
        for atomtype in unique_atom:
            for i, atom in enumerate(atom_list):
                attr = getattr(atom, key)
                if attr == atomtype:
                    new_atom_list.append(atom)
                    index.append(i)
        return new_atom_list, index

    @staticmethod
    def find_uniq_atomtypes(atom_list: List['AbacusATOM']
                            ) -> List[AbacusAtomType]:
        """Return a list of unique AbacusAtomType objects for atom list.
        NOTE: It only merge the consecutive same atom types, which means there may be duplicate types in the returned list.
        If you want to get the true unique atom types, please use sort() method first.
        """
        unique_types = []
        if not atom_list:
            return unique_types
        unique_types.append(atom_list[0].atomtype)
        natom = 1
        for i in range(1, len(atom_list)):
            if atom_list[i].atomtype != atom_list[i-1].atomtype:
                unique_types[-1].natom = natom
                unique_types.append(atom_list[i].atomtype)
                natom = 1
            else:
                natom += 1
        unique_types[-1].natom = natom
        assert sum([atype.natom for atype in unique_types]) == len(atom_list), "Sum of unique atom types' natom does not equal to total number of atoms."
        return unique_types

    @staticmethod
    def same_type(atomlist: List['AbacusATOM']) -> bool:
        """Check if all atoms in the list are of the same type by comparing their label/element/mass/pp/orb/paw/type_mag.

        Args:
            atom_list (List[AbacusATOM]): List of AbacusATOM objects.
        Returns:
            bool: True if all atoms are of the same type, False otherwise.
        """
        if len(atomlist) in [0, 1]:
            return True
        
        for i in range(1, len(atomlist)):
            if atomlist[i].atomtype != atomlist[0].atomtype:
                return False
        return True

class AbacusSTRU:
    """ABACUS STRU class

    """

    def __init__(self,
                 cell: Union[List[Tuple[float,float,float]], np.ndarray],
                 atoms: List[AbacusATOM],
                 dpks: Optional[str] = None,
                 metadata: Dict[str, Any] = {}, 
                 ):
        self._cell = cell
        self._atoms = atoms
        self.dpks = dpks
        if isinstance(self.cell, np.ndarray):
            assert self.cell.ndim == 2 and self.cell.shape == (3,3), "Cell numpy array must be two-dimensional with shape (3,3)."
            self.cell = [tuple(self.cell[i].tolist()) for i in range(3)]
        self.metadata = {
            "lattice_constant": 1.0,
            "atom_type": "cartesian", # or "direct"
        }
        if metadata:
            self.metadata.update(metadata)

    def __str__(self):
        chem_labels = "".join([self.labels_uniq[i]+str(self.labels.count(self.labels_uniq[i])) for i in range(len(self.labels_uniq))])
        c = f"ABACUS STRU object:{chem_labels}\n"
        c += f"NAtoms: {self.natoms}\n"
        c += "Cell vectors (Angstrom):\n"
        for vec in self.cell:
            c += f"  {vec[0]:%12.7f} {vec[1]:%12.7f} {vec[2]:%12.7f}\n"
        return c

    # write len function
    def __len__(self):
        return len(self._atoms)

    def sort(self,
             keep_first_order=True) -> List[int]:
        """Classify atoms according to their atom types.
        The new atom type order is the order of first appearance in the original atom list.

        Returns:
            List[int]: Indices of the rearranged atoms in the original list.
        """
        new_atom_list, index = AbacusATOM.sort(self._atoms, keep_first_order=keep_first_order)
        self._atoms = new_atom_list
        return index
    
    @property
    def natoms(self):
        # Number of atoms in the structure
        return len(self._atoms)
    
    @property
    def cell(self):
        # Return the cell vectors
        return self._cell
    
    @property
    def atoms(self):
        # Return the list of AbacusATOM objects
        return self._atoms

    @property
    def atomtypes(self):
        """Return a list of AbacusAtomType objects for unique atom types in the structure.
        """
        return [atom.atomtype for atom in self._atoms]

    @property
    def uniq_atomtypes(self):
        """Return a list of unique AbacusAtomType objects for the structure.
        NOTE: It only merge the consecutive same atom types, which means there may be duplicate types in the returned list.
        If you want to get the true unique atom types, please use sort() method first.
        """
        return AbacusATOM.find_uniq_atomtypes(self._atoms)

    @property
    def labels(self):
        # Return the list of labels of all atoms
        return [atom.label for atom in self._atoms]

    @property
    def elements(self):
        # Return the list of elements of all atoms
        return [atom.element for atom in self._atoms]
    
    @property
    def masses(self):
        # Return the list of masses of all atoms
        return [atom.mass for atom in self._atoms]
    
    @property
    def pps(self):
        # Return the list of pseudopotential file names of all atoms
        return [atom.pp for atom in self._atoms]
    
    @property
    def orbs(self):
        # Return the list of orbital file names of all atoms
        return [atom.orb for atom in self._atoms]
    
    @property
    def paws(self):
        # Return the list of PAW file names of all atoms
        return [atom.paw for atom in self._atoms]
    
    @property
    def coords(self):
        # Return the list of coordinates of all atoms
        return [atom.coord for atom in self._atoms]
    
    @property
    def coords_direct(self):
        return Cartesian2Direct(self.coords, self.cell)
    
    @property
    def moves(self):
        # Return the list of move flags of all atoms
        return [atom.move for atom in self._atoms]
    
    @property
    def atom_mags(self) -> List[Union[float,Tuple[float,float,float]]]:
        # return magnetic moments of all atoms, will transfer to same type if needed
        # e.g., if one atom is noncolinear, all will be noncolinear
        # if atom.mag is None, value will be atom.type_mag
        noncolinear = any(atom.noncolinear for atom in self._atoms)
        mags = []
        for atom in self._atoms:
            if noncolinear and not atom.noncolinear:
                mags.append( (0.0, 0.0, atom.type_mag) )
            else:
                mags.append(atom.atommag)
        return mags

    @cell.setter
    def cell(self, value: Union[List[Tuple[float,float,float]], np.ndarray]):
        if isinstance(value, np.ndarray):
            assert value.ndim == 2 and value.shape == (3,3), "Cell numpy array must be two-dimensional with shape (3,3)."
            value = [tuple(value[i].tolist()) for i in range(3)]
        else:
            assert len(value) == 3 and all(len(vec) == 3 for vec in value), "Cell must be a list of three 3D vectors."
        self._cell = value
    
    @coords.setter
    def coords(self, value: List[Tuple[float,float,float]]):
        assert len(value) == self.natoms, "Number of coordinates must match number of atoms."
        for i in range(self.natoms):
            self._atoms[i].coord = value[i]
    
    @coords_direct.setter
    def coords_direct(self, value: List[Tuple[float,float,float]]):
        assert len(value) == self.natoms, "Number of coordinates must match number of atoms."
        cart_coords = Direct2Cartesian(value, self.cell)
        for i in range(self.natoms):
            self._atoms[i].coord = cart_coords[i]
    
    def set_pp(self, pp_dict: Dict[str, str], key_type:Literal["element","label"]="element"):
        """Set pseudopotential file names for atoms based on a provided dictionary.

        Args:
            pp_dict (Dict[str, str]): Dictionary mapping element symbols to pseudopotential file names.
        """
        for atom in self._atoms:
            key = atom.element if key_type == "element" else atom.label
            if key in pp_dict:
                atom.pp = pp_dict[key]
            else:
                raise KeyError(f"Pseudopotential for {key_type}: {key} not found in provided dictionary.")

    def set_orb(self, orb_dict: Dict[str, str], key_type: Literal["element","label"]="element"):
        """Set orbital file names for atoms based on a provided dictionary.

        Args:
            orb_dict (Dict[str, str]): Dictionary mapping element symbols to orbital file names.
        """
        for atom in self._atoms:
            key = atom.element if key_type == "element" else atom.label
            if key in orb_dict:
                atom.orb = orb_dict[key]
            else:
                raise KeyError(f"Orbital for {key_type}: {key} not found in provided dictionary.")
    
    def set_paw(self, paw_dict: Dict[str, str], key_type: Literal["element","label"]="element"):
        """Set PAW file names for atoms based on a provided dictionary.

        Args:
            paw_dict (Dict[str, str]): Dictionary mapping element symbols to PAW file names.
        """
        for atom in self._atoms:
            key = atom.element if key_type == "element" else atom.label
            if key in paw_dict:
                atom.paw = paw_dict[key]
            else:
                raise KeyError(f"PAW for {key_type}: {key} not found in provided dictionary.")
    
    def set_coords(self, coords: List[Tuple[float,float,float]], direct: bool=False):
        """
        Set coordinates of all atoms using given coordinates.
        Args:
            coords (List[Tuple[float,float,float]]): List of coordinates for each atom. If direct is False, the unit of coordinates is Angstrom.
            direct (bool): Whether the coordinates are in direct or cartesian coordinates. Default is False.
        """
        if direct:
            coords = Direct2Cartesian(coords, self.cell)
        for i in range(len(self._atoms)):
            self._atoms[i].coord = coords[i]

    @staticmethod
    def read(filename: str, fmt: Literal["stru", "abacus/stru", "poscar","vasp", "cif"]="stru") -> "AbacusSTRU":
        """Read structure from a file in the specified format.

        Args:
            filename (str): Input file name.
            fmt (str): Format of the input file. Options are "stru", "poscar", "cif". Default is "stru".

        Returns:
            AbacusSTRU: An AbacusSTRU object representing the structure.
        """
        fmt = fmt.lower()
        if fmt in ["stru", "abacus/stru"]:
            stru_data = read_stru_file(stru=filename)
            cell = (np.array(stru_data["cell"]) * stru_data['lattice_constant'] * BOHR2A).tolist()
            if stru_data["cartesian"]:
                coords = (np.array(stru_data["coord"]) * stru_data['lattice_constant'] * BOHR2A).tolist()
            else:
                coords = Direct2Cartesian(stru_data["coord"], cell)
            atom_list = []
            label_tot = get_total_property(stru_data, "label")
            pp_tot = get_total_property(stru_data, "pp")
            orb_tot = get_total_property(stru_data, "orb")
            paw_tot = get_total_property(stru_data, "paw")
            type_mag_tot = get_total_property(stru_data, "magmom")

            for i in range(len(coords)):
                atom = AbacusATOM(
                    label=label_tot[i],
                    coord=tuple(coords[i]),
                    element=None,
                    mass=None,
                    pp=None if len(stru_data['pp']) == 0 else pp_tot[i],
                    orb=None if len(stru_data['orb']) == 0 else orb_tot[i],
                    paw=None if len(stru_data['paw']) == 0 else paw_tot[i],
                    type_mag=type_mag_tot[i],
                    move=stru_data["move"][i],
                    mag=stru_data["magmom_atom"][i],
                    angle1=stru_data["angle1"][i],
                    angle2=stru_data["angle2"][i],
                    velocity=stru_data["velocity"][i],
                    constrain=stru_data["constrain"][i],
                    lambda_=stru_data["lambda_"][i],
                )
                atom_list.append(atom)
            dpks = stru_data.get("dpks", None)
            metadata = {
                "lattice_constant": stru_data.get("lattice_constant", 1.0),
                "atom_type": "cartesian" if stru_data.get("cartesian", True) else "direct",
            }
            return AbacusSTRU(cell=cell, atoms=atom_list, dpks=dpks, metadata=metadata)
        elif fmt in ["poscar", "vasp", "cif"]:
            # use ase to read poscar/vasp/cif file
            from ase.io import read as ase_read
            if fmt == "poscar":
                fmt = "vasp"
            atom_type = "direct"
            if fmt == "vasp":
                with open(filename, 'r') as f: lines = f.readlines()
                if lines[5].strip().lower().startswith("c"):
                    atom_type = "cartesian"
    
            ase_stru = ase_read(filename, format=fmt)
            return AbacusSTRU.from_ase(ase_stru, meta_data={
                "lattice_constant": A2BOHR,
                "atom_type": atom_type,
            })
        else:
            raise ValueError(f"Unsupported format: {fmt}")

    def write(self, filename: str,
              fmt: Literal["stru", "poscar","vasp", "cif"]="stru",
              empty2x: bool = False,
              direct: Optional[bool]=None):
        """Write the structure to a file in the specified format.
        Args:
            filename (str): Output file name.
            fmt (str): Format of the output file. Options are "stru", "poscar", "cif". Default is "stru".
            empty2x (bool): If True, convert 'empty' atoms to 'X' element before writing. Default is False.
            direct (bool, optional): If True, write atomic positions in direct coordinates. If False, write in cartesian coordinates. If None, use the value from metadata. Default is None.
        """
        fmt = fmt.lower()
        if direct is None:
            direct = (self.metadata.get("atom_type","cartesian").lower() == "direct")
        else:
            direct = direct
        
        elements = copy.deepcopy(self.elements)
        if empty2x:
            for i in range(len(self.labels)):
                if "empty" in self.labels[i]:
                    elements[i] = 'X'


        if fmt in ["stru", "abacus/stru"]:
            atom_list = copy.deepcopy(self._atoms)
            unique_types = AbacusATOM.find_uniq_atomtypes(atom_list)
            lc = self.metadata.get("lattice_constant", 1.0)
            cell = np.array(self.cell) * A2BOHR / lc
            coord = np.array([atom.coord for atom in atom_list]) * A2BOHR / lc
            if direct:
                coord = Cartesian2Direct(coord.tolist(), self.cell)
            else:
                coord = coord.tolist()
            cell = cell.tolist()

            write_stru_file(cell=cell, coord=coord, 
                            label=[ut.label for ut in unique_types],
                            atom_number=[ut.natom for ut in unique_types],
                            struf=filename,
                            direct=direct,
                            pp =[ut.pp for ut in unique_types],
                            mass = [ut.mass for ut in unique_types],
                            orb = [ut.orb for ut in unique_types],
                            paw = [ut.paw for ut in unique_types],
                            magmom_global=[ut.type_mag for ut in unique_types],
                            lattice_constant=lc,
                            move = [atom.move for atom in  atom_list],
                            magmom = [atom.mag for atom in atom_list],
                            velocity = [atom.velocity for atom in atom_list],
                            angle1 = [atom.angle1 for atom in atom_list],
                            angle2 = [atom.angle2 for atom in atom_list],
                            constrain = [atom.constrain for atom in atom_list],
                            lambda_ = [atom.lambda_ for atom in atom_list],
                            dpks = self.dpks)
        elif fmt in  ["poscar", "vasp"]:
            write_poscar(cell = self.cell,
                         coord=self.coords if not direct else self.coords_direct,
                         label=elements, poscar=filename, direct=direct, move=self.moves)
        elif fmt == "cif":
            ase_stru = self.to(fmt="ase", empty2x=empty2x)
            ase_stru.write(filename, format="cif")
        else:
            raise ValueError(f"Unsupported format: {fmt}")
    
    @staticmethod
    def from_ase(ase_stru,
                 meta_data: Dict[str, Any] = {}) -> "AbacusSTRU":
        """Create an AbacusSTRU object from an ASE Atoms object.
        Args:
            ase_stru: ASE Atoms object.
            meta_data (Dict[str, Any], optional): Additional metadata for the AbacusSTRU object. Default is an empty dictionary.
        Returns:
            AbacusSTRU: An AbacusSTRU object representing the structure.
        """
        from ase import Atoms
        assert isinstance(ase_stru, Atoms), "Input structure must be an ASE Atoms object."
        cell = ase_stru.get_cell().tolist()
        atom_list = []
        mags = ase_stru.get_magnetic_moments() if ase_stru.has_magnetic_moments() else [None]*len(ase_stru)
        for i in range(len(ase_stru)):
            atom = AbacusATOM(
                label=ase_stru[i].symbol,
                coord=tuple(ase_stru[i].position.tolist()),
                element=ase_stru[i].symbol,
                mass=ase_stru[i].mass,
                pp=ase_stru.info.get("pp", {}).get(ase_stru[i].symbol, None),
                orb=ase_stru.info.get("orb", {}).get(ase_stru[i].symbol, None),
                paw=ase_stru.info.get("paw", {}).get(ase_stru[i].symbol, None),
                type_mag= 0.0 if mags[i] is None else mags[i],
                move=(True, True, True),
                mag= mags[i],
            )
            atom_list.append(atom)

        return AbacusSTRU(cell=cell, atom_list=atom_list, meta_data=meta_data)

    @staticmethod
    def from_pymatgen(pymatgen_stru,
        meta_data: Dict[str, Any] = {}) -> "AbacusSTRU":
        """Create an AbacusSTRU object from a Pymatgen Structure object.
        Args:
            pymatgen_stru: Pymatgen Structure object.
            meta_data (Dict[str, Any], optional): Additional metadata for the AbacusSTRU object. Default is an empty dictionary.
        Returns:
            AbacusSTRU: An AbacusSTRU object representing the structure.
        """
        from pymatgen.core import Structure
        assert isinstance(pymatgen_stru, Structure), "Input structure must be a Pymatgen Structure object."
        cell = pymatgen_stru.lattice.matrix.tolist()
        atom_list = []
        for site in pymatgen_stru.sites:
            atom = AbacusATOM(
                label=site.specie.symbol,
                coord=tuple(site.coords.tolist()),
                element=site.specie.symbol,
                mass=site.specie.atomic_mass,
                type_mag=0.0,
                move=(True, True, True),
            )
            atom_list.append(atom)

        return AbacusSTRU(cell=cell, atom_list=atom_list, meta_data=meta_data)
    
    @staticmethod
    def from_dpdata(dpdata_stru,
                    index : int = 0,
                    meta_data: Dict[str, Any] = {}) -> "AbacusSTRU":
        """Create an AbacusSTRU object from a DPData System object.
        Args:
            dpdata_stru: DPData System object.
            index (int): Index of the configuration in the DPData System to convert. Default is 0.
            meta_data (Dict[str, Any], optional): Additional metadata for the AbacusSTRU object. Default is an empty dictionary.
        Returns:
            AbacusSTRU: An AbacusSTRU object representing the structure.
        """
        from dpdata import System
        assert isinstance(dpdata_stru, System), "Input structure must be a DPData System object."
        data = dpdata_stru.data
        coords = data['coords'][index]
        atom_list = []
        for i in range(dpdata_stru.get_natoms()):
            atom_list.append(AbacusATOM(
                label=data['atom_names'][data["atom_types"][i]],
                coord=tuple(coords[i].tolist()),
                element=data['atom_names'][data["atom_types"][i]],
                mass=data["masses"][data["atom_types"][i]],
                type_mag=0.0,
                move=(True, True, True) if "move" not in data else tuple(data["move"][index][i]),
                mag = None if "spins" not in data else data["spins"][index][i]
            ))
        cell = data['cells'][index].tolist()
        return AbacusSTRU(cell=cell, atom_list=atom_list, meta_data=meta_data)
        
    @staticmethod
    def from_phonopy(phonopy_stru,
                     meta_data: Dict[str, Any] = {}) -> "AbacusSTRU":
        """Create an AbacusSTRU object from a Phonopy Atoms object.
        Args:
            phonopy_stru: Phonopy Atoms object.
            meta_data (Dict[str, Any], optional): Additional metadata for the AbacusSTRU object
        Returns:
            AbacusSTRU: An AbacusSTRU object representing the structure.
        """
        from phonopy.structure.atoms import PhonopyAtoms
        assert isinstance(phonopy_stru, PhonopyAtoms), "Input structure must be a Phonopy Atoms object."
        cell = phonopy_stru.cell.tolist()
        atom_list = []
        for i in range(len(phonopy_stru)):
            atom = AbacusATOM(
                label=phonopy_stru.get_chemical_symbols()[i],
                coord=tuple(phonopy_stru.get_positions()[i].tolist()),
                element=phonopy_stru.get_chemical_symbols()[i],
                mass=phonopy_stru.get_masses()[i],
                type_mag=0.0,
                move=(True, True, True),
            )
            atom_list.append(atom)
        return AbacusSTRU(cell=cell, atom_list=atom_list, meta_data=meta_data)
        

    def to(self, fmt: Literal["ase", "pymatgen", "dpdata", "phonopy"],
           empty2x: bool = False):
        """Convert the structure to another format.

        Args:
            fmt (str): Target format. Options are:
             - "ase": ASE Atoms object.
             - "pymatgen": Pymatgen Structure object.
             - "dpdata": DPData System object.
             - "phonopy": Phonopy Atoms object.
            empty2x (bool): If True, convert 'empty' atoms to 'X' element. Default is False.

        Returns:
            Converted structure object.
        """
        elements = self.elements
        if empty2x:
            labels = self.labels
            for i in range(len(labels)):
                if 'empty' in labels[i]:
                    elements[i] = 'X'

        if fmt == "ase":
            from ase import Atoms
            return Atoms(symbols=elements, positions=self.coords, cell=self.cell, pbc=True, magmoms=self.atom_mags,
                          info={"pp": self.pps,"orb": self.orbs,"paw": self.paws,"dpks": self.dpks})
        elif fmt == "pymatgen":
            from pymatgen.core import Structure
            return Structure(lattice=self.cell, species=elements, coords=self.coords, coords_are_cartesian=True,labels=self.labels,)
        elif fmt == "dpdata":
            raise NotImplementedError("DPData format conversion is not implemented yet.")
        elif fmt == "phonopy":
            from phonopy.structure.atoms import PhonopyAtoms
            return PhonopyAtoms(symbols=elements, positions=self.coords, cell=self.cell, masses=self.masses,
                            magnetic_moments=self.atom_mags)
        else:
            raise ValueError(f"Unsupported format: {fmt}")

    def fix_atom_by_index(self, indices: List[int],
                          move: Optional[Tuple[bool, bool, bool]] = (False, False, False),
                          only: Optional[bool]=False):
        """
        Fix atoms by index.
        
        Args:
            indices (List[int]): List of indices of atoms to fix. Starts from 0.
            move (Tuple[bool, bool, bool]): Tuple indicating whether the atom is allowed to move in each direction. Default is (False, False, False), which means all 3 directions are fixed.
            only (bool): If True, override the move settings for unselected atoms and allow all unselected atoms to move. Default is False.
        """
        for i in range(self.natoms):
            if i in indices:
                self._atoms[i].move = move
            elif only:
                self._atoms[i].move = (True, True, True)

    def fix_atom_by_coord(self,
                          min: float,
                          max: float,
                          cartesian: Optional[bool]=True,
                          direction: Optional[Literal[0, 1, 2]]=2, 
                          move: Optional[Tuple[bool, bool, bool]]=(False, False, False),
                          only: Optional[bool]=False):
        """
        Fix atoms by coordinates (cartesian or direct).
        
        Args:
            min (float): Minimum height of atoms to fix.
            max (float): Maximum height of atoms to fix.
            cartesian (bool): Whether to use cartesian coordinates. Default is True.
            direction (int): Direction of atoms to fix. Can be 0, 1 or 2, means 'x', 'y' and 'z' respectively. Default is 2.
            move: Tuple[bool, bool, bool]: Tuple indicating whether the selected atoms are allowed to move in each direction. Default is (False, False, False), which means all 3 directions are fixed.
            only (bool): If True, override the move settings for unselected atoms and allow all unselected atoms to move. Default is False.
        """
        for i in range(self.natoms):
            if cartesian:
                if self.coords_angs[i][direction] >= min and self.coords_angs[i][direction] <= max:
                    self._atoms[i].move = move
                elif only:
                    self._atoms[i].move = (True, True, True)
            else:
                if self.coords_direct[i][direction] >= min and self.coords_direct[i][direction] <= max:
                    self._atoms[i].move = move
                elif only:
                    self._atoms[i].move = (True, True, True)



def parse_stru_position(pos_line):
    '''
  The content in atom position block:
  - `m` or NO key word: three numbers, which take value in 0 or 1, control how the atom move in geometry relaxation calculations. In example below, the numbers `0 0 0` following the coordinates of the first atom means this atom are *not allowed* to move in all three directions, and the numbers `1 1 1` following the coordinates of the second atom means this atom *can* move in all three directions.
  - `v` or `vel` or `velocity`: set the three components of initial velocity of atoms in geometry relaxation calculations(e. g. `v 1.0 1.0 1.0`).
  - `mag` or `magmom` : set the start magnetization for each atom. In colinear case only one number should be given. In non-colinear case one have two choice:either set one number for the norm of magnetization here and specify two polar angle later(e. g. see below), or set three number for the xyz commponent of magnetization here (e. g. `mag 0.0 0.0 1.0`). Note that if this parameter is set, the initial magnetic moment setting in the second line will be overrided.
    - `angle1`: in non-colinear case, specify the angle between c-axis and real spin, in angle measure instead of radian measure
    - `angle2`: in non-colinear case, specify angle between a-axis and real spin in projection in ab-plane , in angle measure instead of radian measure

      e.g.:

      ```
      Fe
      1.0
      2
      0.0 0.0 0.0 m 0 0 0 mag 1.0 angle1 90 angle2 0 cs 0 0 0
      0.5 0.5 0.5 m 1 1 1 mag 1.0 angle1 90 angle2 180
      ```
    '''
    sline = pos_line.split()
    pos = [float(i) for i in sline[:3]]
    move = None
    velocity = None
    magmom = None
    angle1 = None
    angle2 = None
    constrain = None
    lambda1 = None
    if len(sline) > 3:
        mag_list = []
        velocity_list = []
        move_list = []
        angle1_list = []
        angle2_list = []
        constrain_list = []
        lambda_list = []
        label = "move"
        for i in range(3,len(sline)):
            # firstly read the label
            if sline[i] == "m":
                label = "move"
                move_list = []
            elif sline[i] in ["v","vel","velocity"]:
                label = "velocity"
                velocity_list = []
            elif sline[i] in ["mag","magmom"]:
                label = "magmom"
                mag_list = []
            elif sline[i] == "angle1":
                label = "angle1"
                angle1_list = []
            elif sline[i] == "angle2":
                label = "angle2"
                angle2_list = []
            elif sline[i] in ["constrain","sc"]:
                label = "constrain"
                constrain_list = []
            elif sline[i] in ["lambda"]:
                label = "lambda"
                lambda_list = []
            
            # the read the value to the list    
            elif label == "move":
                move_list.append(int(sline[i]))
            elif label == "velocity":
                velocity_list.append(float(sline[i]))
            elif label == "magmom":
                mag_list.append(float(sline[i]))
            elif label == "angle1":
                angle1_list.append(float(sline[i]))
            elif label == "angle2":
                angle2_list.append(float(sline[i]))
            elif label == "constrain":
                constrain_list.append(bool(int(sline[i])))
            elif label == "lambda":
                lambda_list.append(float(sline[i]))
        if len(move_list) == 3:
            move = move_list
        if len(velocity_list) == 3:
            velocity = velocity_list
        if len(mag_list) in [1,3]:
            magmom = mag_list if len(mag_list) == 3 else mag_list[0]
        if len(angle1_list) == 1:
            angle1 = angle1_list[0]
        if len(angle2_list) == 1:
            angle2 = angle2_list[0]
        if len(constrain_list) == 3:
            constrain = constrain_list
        elif len(constrain_list) == 1:
            constrain = constrain_list[0]
        if len(lambda_list) == 3:
            lambda1 = lambda_list
        elif len(lambda_list) == 1:
            lambda1 = lambda_list[0]
            
            
    return pos,move,velocity,magmom,angle1,angle2,constrain,lambda1

def read_stru_file(stru:str = "STRU"):
    '''Read ABACUS STRU file and return a dictionary with structure information.
    Args:
        stru (str): Path to the STRU file. Default is "STRU".
    Returns:
        dict: A dictionary containing structure information with keys:
            - label: list of labels for each atom type
            - atom_number: list of atom numbers for each atom type
            - cell: 3x3 list of cell vectors
            - coords: list of coordinates for each atom
            - pp: list of pseudopotential files for each atom type
            - orb: list of orbital files for each atom type
            - paw: list of PAW files for each atom type
            - lattice_constant: lattice constant value
            - move: list of move flags for each atom
            - magmom_global: list of global magnetic moments for each atom type
            - magmom: list of magnetic moments for each atom
            - velocity: list of velocities for each atom
            - angle1: list of angle1 values for each atom
            - angle2: list of angle2 values for each atom
            - constrain: list of constraints for each atom
            - lambda1: list of lambda values for each atom
    
    NOTE:
        1. Do not support bravais lattice now.
        2. the value for cell/coords/lattice_constant are the exact values in STRU file
        3. Only direct/cartessian coordinate type is supported.
    
    '''
    def get_block(keyname):
        block = []
        for i,line in enumerate(lines):
            if line.strip() == "": continue
            elif line.split('#')[0].strip() == keyname:
                for ij in range(i+1,len(lines)):
                    if lines[ij].strip() == "" or \
                        lines[ij].strip()[0] in ["#"] or\
                        ("//" in lines[ij] and lines[ij].strip()[:2] in ["//"]): continue
                    elif lines[ij].strip() in ABACUS_STRU_KEY_WORD:
                        return block
                    else:
                        block.append(lines[ij].split("#")[0].split("//")[0].strip())
                return block
        return None
    
    if not os.path.isfile(stru):
        return None
    with open(stru) as f1: lines = f1.readlines()  
    atomic_species = get_block("ATOMIC_SPECIES")
    numerical_orbital = get_block("NUMERICAL_ORBITAL")
    lattice_constant = get_block("LATTICE_CONSTANT")
    lattice_vector = get_block("LATTICE_VECTORS")
    atom_positions = get_block("ATOMIC_POSITIONS")
    dpks = get_block("NUMERICAL_DESCRIPTOR")
    pawf = get_block("PAW_FILES")
    lattice_constant = 1.0 if lattice_constant == None else float(lattice_constant[0].split()[0]) 
    dpks = None if dpks == None else dpks[0].strip()
    
    #read species
    pp = []
    labels = []
    mass = []
    for line in atomic_species:
        sline = line.split()
        labels.append(sline[0])
        mass.append(float(sline[1]))
        if len(sline) > 2: 
            pp.append(sline[2])
    if len(pp) == 0:
        pp = None
        
    #read orbital
    if numerical_orbital == None:
        orb = None
    else:
        orb = []
        for line in numerical_orbital:
            orb.append(line.split()[0])
    
    # read paw files
    if pawf == None:
        paw = None
    else:
        paw = []
        for line in pawf:
            paw.append(line.split("#")[0].strip())
    
    #read cell
    cell = []
    try:
        for line in lattice_vector:
            cell.append([float(i) for i in line.split()[:3]])
    except:
        traceback.print_exc()
        print("WARNING: LATTICE_VECTORS is incorrect !!!!!!")
    #read coordinate and coordinate type and atom number of each type
    atom_number = []
    coords = []
    magmom_global = [] # the initial magmom of each type
    magmom = [] # the initial magmom of each atom
    move = []
    velocity = []
    angle1 = []
    angle2 = []
    constrain = [] # the constrain of each atom
    lambda1 = [] # the lambda for delta spin
    coord_type = atom_positions[0].split("#")[0].strip().lower()
    if coord_type.startswith("dire"):
        cartesian = False
    elif coord_type.startswith("cart"):
        cartesian = True
    else:
        print("Not support coordinate type %s now." % atom_positions[0].strip())
        sys.exit(1)
    i = 1
    real_label = []
    real_pp = []
    real_orb = []
    real_paw = []
    while i < len(atom_positions):
        label = atom_positions[i].strip()
        if label not in labels:
            print("label '%s' is not matched that in ATOMIC_SPECIES" % label)
            sys.exit(1)
        an = int(atom_positions[i+2].split()[0])
        if an == 0:
            i += 3
            continue
        
        real_label.append(label)
        label_idx = labels.index(label)
        if pp:
            real_pp.append(pp[label_idx])
        if orb:
            real_orb.append(orb[label_idx])
        if paw:
            real_paw.append(paw[label_idx])
            
        magmom_global.append(float(atom_positions[i+1].split()[0]))
            
        atom_number.append(an)
            
        i += 3
        for j in range(atom_number[-1]):
            pos,imove,ivelocity,imag,iangle1,iangle2,iconstrain,ilambda1 = parse_stru_position(atom_positions[i+j])
            coords.append(pos)
            move.append(imove)
            velocity.append(ivelocity)
            magmom.append(imag)
            angle1.append(iangle1)
            angle2.append(iangle2)
            constrain.append(iconstrain)
            lambda1.append(ilambda1)
            
        i += atom_number[-1]

    return {
        "label": real_label,   # list of labels for each atom type
        "atom_number": atom_number,  # list of atom numbers for each atom type
        "cell": cell,              # 3x3 list of cell vectors
        "coord": coords,       # list of coordinates for each atom
        "pp": real_pp,          # list of pseudopotential files for each atom type
        "orb": real_orb,        # list of orbital files for each atom type
        "paw": real_paw,        # list of PAW files for each atom type
        "lattice_constant": lattice_constant,  # lattice constant value
        "move": move,              # list of move flags for each atom
        "magmom": magmom_global, # list of global magnetic moments for each atom type
        "magmom_atom": magmom,       # list of magnetic moments for each atom
        "velocity": velocity,    # list of velocities for each atom
        "angle1": angle1,       # list of angle1 values for each atom
        "angle2": angle2,       # list of angle2 values for each atom
        "constrain": constrain, # list of constrain flags for each atom
        "lambda_": lambda1,     # list of lambda1 values for each atom
        "dpks": dpks,           # dpks value
        "cartesian": cartesian  # boolean indicating if coordinates are cartesian
    }



def write_stru_file(
          cell:List[Tuple[float,float,float]],
          coord: List[Tuple[float,float,float]],
          label:List[str],
          atom_number:List[int],
          struf: str ="STRU", 
          direct:bool=False,
          pp:Optional[List[str]] = None,
          mass:Optional[List[float]] = None,
          orb:Optional[List[str]] = None,
          paw:Optional[List[str]] = None,
          magmom_global:Optional[List[float]] = None,
          lattice_constant:float=1.0,
          move:Optional[List[Optional[Tuple[bool,bool,bool]]]] = None,
          magmom:Optional[List[float]] = None,
          velocity:Optional[List[Optional[Tuple[float,float,float]]]] = None,
          angle1:Optional[List[Optional[float]]] = None,
          angle2:Optional[List[Optional[float]]] = None,
          constrain:Optional[List[Optional[Union[bool,Tuple[bool,bool,bool]]]]] = None,
          lambda_:Optional[List[Optional[Union[float,Tuple[float,float,float]]]]] = None,
          dpks:Optional[str] = None,
          ):
    '''Write to ABACUS STRU file. 
    1. The value will be directly write to STRU file, no unit conversion.
       So the unit of cell/coord/lattice_constant should be bohr, and the coord should be in direct/cartesian type according to `direct` parameter.
       The real cell should be cell * lattice_constant.
    2. The order of pp/mass/orb/paw/magmom_global should be the same as label list.
    3. The coord should be in the order of label and atom_number.
    4. The order of coord/move/magmom/velocity/angle1/angle2/constrain/lambda_ should be the same as coord list.
    '''
    # check parameters
    assert len(label) == len(atom_number), "label and atom_number length mismatch"
    natoms = sum(atom_number)
    assert len(coord) == natoms, "coord length mismatch with atom_number"
    for key, lst in [("pp",pp),("mass",mass),("magmom_global",magmom_global),
                     ("orb",orb),("paw",paw)]:
        if lst is not None:
            assert len(lst) == len(label), f"{key} length mismatch with label"
    for key, lst in [("move",move),("magmom",magmom),("velocity",velocity),
                     ("angle1",angle1),("angle2",angle2),("constrain",constrain),("lambda_",lambda_)]:
        if lst is not None:
            assert len(lst) == natoms, f"{key} length mismatch with natoms"
    
    cc = ""
    #write species
    cc += "ATOMIC_SPECIES\n"
    for i, ilabel in enumerate(label):
        cc += f"{ilabel} "
        if mass and mass[i] is not None:
            cc += f"{mass[i]} "
        else:
            cc += "1.0 "
        if pp and pp[i] is not None:
            cc += f"{pp[i]}\n"
        else:
            cc += "\n"
    
    if orb and None not in orb:
        #write orb
        cc += "\nNUMERICAL_ORBITAL\n"
        for i in orb:
            cc += i + "\n"
    elif all([i is None for i in orb]):
        pass
    elif orb and None in orb:
        print("WARNING: orb list contains None, skip writing NUMERICAL_ORBITAL block")
    
    if paw and None not in paw:
        # write pawfile
        cc += "\nPAW_FILES\n"
        for i in paw:
            cc += i + "\n"
    elif all([i is None for i in paw]):
        pass
    elif paw and None in paw:
        print("WARNING: paw list contains None, skip writing PAW_FILES block")
    
    # write LATTICE_CONSTANT
    cc += "\nLATTICE_CONSTANT\n%f\n" % lattice_constant
    #write LATTICE_VECTORS
    cc += "\nLATTICE_VECTORS\n"
    for i in cell:
        cc += "%17.11f %17.11f %17.11f\n" % tuple(i)
    #write ATOMIC_POSITIONS
    cc += "\nATOMIC_POSITIONS\n"
    if not direct:
        cc += "Cartesian\n"
    else:
        cc += "Direct\n"
    icoord = 0
    for i,ilabel in enumerate(label):
        cc += "\n%s\n" % ilabel
        if magmom_global and magmom_global[i] is not None:
            cc += "%f\n" % magmom_global[i]
        else:
            cc += "0.0\n"
        cc += "%d\n" % atom_number[i]
        for j in range(atom_number[i]):
            cc += "%17.11f %17.11f %17.11f " % tuple(coord[icoord + j])
            if move and move[icoord + j] and len(move[icoord + j]) == 3:
                cc += "%d %d %d " % tuple(move[icoord + j])
            if magmom and magmom[icoord + j] is not None:
                if isinstance(magmom[icoord + j],list):
                    if len(magmom[icoord + j]) == 3:
                        cc += "mag %12.8f %12.8f %12.8f " % tuple(magmom[icoord + j])
                    elif len(magmom[icoord + j]) == 1:
                        cc += "mag %12.8f " % magmom[icoord + j][0]
                elif magmom[icoord + j] != None:
                    cc += "mag %12.8f " % magmom[icoord + j]
            if velocity and velocity[icoord + j] and len(velocity[icoord + j]) == 3:
                cc += "v %f %f %f " % tuple(velocity[icoord + j])
            if angle1 and angle1[icoord + j] != None:
                    cc += "angle1 %f " % angle1[icoord + j]
            if angle2 and angle2[icoord + j] != None:
                    cc += "angle2 %f " % angle2[icoord + j]
            if constrain and constrain[icoord + j]:
                if isinstance(constrain[icoord + j],list) and len(constrain[icoord + j]) == 3:
                    cc += "sc " + " ".join(["1" if ic else "0" for ic in constrain[icoord + j]]) + " " 
                elif isinstance(constrain[icoord + j],list) and len(constrain[icoord + j]) == 1:
                    cc += "sc " + ("1" if constrain[icoord + j][0] else "0") + " "
                elif not isinstance(constrain[icoord + j],list):
                    cc += "sc " + ("1" if constrain[icoord + j] else "0") + " "
                else:
                    print("ERROR: the constrain is not a list or a bool value, skip it")
                    print("\t\tconstrain:",constrain[icoord + j])
            if lambda_ and lambda_[icoord + j]:
                if isinstance(lambda_[icoord + j],list):
                    cc += "lambda " + " ".join([str(ic) for ic in lambda_[icoord + j]]) + " "
                else:
                    cc += "lambda " + str(lambda_[icoord + j]) + " "
            cc += "\n"
        icoord += atom_number[i]

    if dpks:
        cc += "\nNUMERICAL_DESCRIPTOR\n"
        cc += dpks    

    os.makedirs(os.path.dirname(struf), exist_ok=True)
    with open(struf,"w") as f1:
        f1.write(cc)
    
    return cc

def write_poscar(
          cell:List[Tuple[float,float,float]],
          coord: List[Tuple[float,float,float]],
          label:List[str],
          poscar: str ="POSCAR", 
          direct:bool=True,
          move:Optional[List[Optional[Tuple[bool,bool,bool]]]] = None,
          ):
    '''Write to VASP POSCAR file. 
    
    Args:
        cell (List[Tuple[float,float,float]]): 3x3 list of cell vectors in Angstrom.
        coord (List[Tuple[float,float,float]]): List of coordinates for each atom in Angstrom (if direct=False) or in direct coordinates (if direct=True).
        label (List[str]): List of labels for each atom.
        poscar (str): Path to the POSCAR file. Default is "POSCAR".
        direct (bool): If True, write atomic positions in direct coordinates. If False, write in cartesian coordinates. Default is True.
        move (Optional[List[Optional[Tuple[bool,bool,bool]]]]): List of move flags for each atom. Default is None.
    '''
    # check parameters
    assert len(label) == len(coord), "label and coord length mismatch"
    if move:
        assert len(move) == len(coord), "move length mismatch with coord"
    
    # find unique labels and their counts
    unique_labels = [label[0]]
    atom_number = []
    natoms = 1
    for i in range(1,len(label)):
        if label[i] != label[i-1]:
            atom_number.append(natoms)
            unique_labels.append(label[i])
            natoms = 1
        else:
            natoms += 1
    atom_number.append(natoms) 

    cc = ""
    #write header
    cc += "poscar\n"
    cc += "1.0\n"
    #write cell
    for i in cell:
        cc += "%17.11f %17.11f %17.11f\n" % tuple(i)
    #write atom types and numbers
    cc += " ".join(unique_labels) + "\n"
    cc += " ".join([str(num) for num in atom_number]) + "\n"
    #write coordinates
    if direct:
        cc += "Direct\n"
    else:
        cc += "Cartesian\n"
    for i in range(len(coord)):
        cc += "%17.11f %17.11f %17.11f" % tuple(coord[i])
        if move and move[i] and len(move[i]) == 3:
            cc += " " + " ".join(["T" if mv else "F" for mv in move[i]])
        cc += "\n"
    os.makedirs(os.path.dirname(poscar), exist_ok=True)
    with open(poscar,"w") as f1:
        f1.write(cc)
    return cc

def get_total_property(stru_data: Dict[str, Any],
                       prop: Literal['label', 'magmom', 'pp', 'orb', 'paw']):
    """
    Get selected property in stru_data (read by read_stru_file) for each atom stored for each type .
    Args:
        stru_data (Dict[str, Any]): The structure data.
        prop (Literal['label', 'magmom', 'pp', 'orb', 'paw']): The property name.
    Returns:
        float: The total property.
    """
    result = []
    if len(stru_data[prop]) == 0:
        return []
    else:
        assert len(stru_data[prop]) == len(stru_data['atom_number'])
        for item, count in zip(stru_data[prop], stru_data['atom_number']):
            result.extend([item] * count)

        return result
