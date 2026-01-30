import numpy as np
from typing import List, Union, Optional, Tuple, Literal
import os
from abacustest.constant import BOHR2A, RY2EV, PERIOD_DICT_NUMBER
from abacustest import AbacusStru
from abacustest.lib_model.comm import chg2pot


class Grid:
    """The class for charge density and potential data
    
    Attributes:
        data (np.ndarray): 3D numpy array of the grid data (charge density or potential), shape (Nx, Ny, Nz)
        cell (np.ndarray): 3x3 numpy array representing the lattice vectors
        atom_positions (np.ndarray): Nx3 numpy array of atomic positions in Cartesian coordinates
        atom_types (List[int]): List of atomic types (atomic numbers)
        atom_charges (List[float]): the charge of each atom
        origin (np.ndarray): 1D numpy array of length 3 representing the origin of the grid in Cartesian coordinates
    """
    def __init__(self, 
                 data: np.ndarray,
                 cell: np.ndarray,
                 atom_positions: Optional[np.ndarray] = None,
                 atom_types: Optional[np.ndarray] = None,
                 atom_charges: Optional[np.ndarray] = None,
                 origin: Optional[np.ndarray] = np.zeros(3),
                 
                 ):
        assert data.ndim == 3, "Data should be a 3D numpy array"
        assert cell.shape == (3,3), "Cell should be a 3x3 numpy array"
        
        if atom_positions is not None:
            assert atom_positions.ndim == 2 and atom_positions.shape[1] == 3, "Atom positions should be a Nx3 numpy array"

            if atom_types is not None:
                assert len(atom_types) == atom_positions.shape[0], "Atom types should be a list of length N"
            else:
                atom_types = [0] * atom_positions.shape[0]  # Default type if not provided
                
            if atom_charges is not None:
                assert len(atom_charges) == atom_positions.shape[0], "Atom charges should be a list of length N"
            else:
                atom_charges = [0.0] * atom_positions.shape[0]  # Default charge if not provided
        else:
            atom_positions = np.zeros(3)
            atom_types = [0]
            atom_charges = [0]
                
        assert origin.shape == (3,), "Origin should be a 1D numpy array of length 3"
        
        self._data_dict = {
            'data': data,
            'cell': cell,
            'atom_positions': atom_positions,
            'atom_types': atom_types,
            'atom_charges': atom_charges,
            'origin': origin
        }
    
    @property
    def data(self):
        return self._data_dict['data']
    
    @data.setter
    def data(self, value: np.ndarray):
        assert value.shape == self._data_dict['data'].shape, "New data must have the same shape as the original data"
        self._data_dict['data'] = value

    @property
    def data_coord(self):
        """Get the coordinates of the grid points in Cartesian coordinates.
        
        Returns:
            np.ndarray: A (Nx*Ny*Nz, 3) array of grid point coordinates.
        """
        nx, ny, nz = self.data.shape
        x = np.linspace(0, 1, nx, endpoint=False)
        y = np.linspace(0, 1, ny, endpoint=False)
        z = np.linspace(0, 1, nz, endpoint=False)
        xv, yv, zv = np.meshgrid(x, y, z, indexing='ij')
        fractional_coords = np.vstack([xv.ravel(), yv.ravel(), zv.ravel()]).T
        cartesian_coords = fractional_coords @ self.cell
        # reshape to (Nx, Ny, Nz, 3)
        cartesian_coords = cartesian_coords.reshape((nx, ny, nz, 3))
        return cartesian_coords
    
    @property
    def cell(self):
        return self._data_dict['cell']
    
    @property
    def atom_positions(self):
        return self._data_dict['atom_positions']
    
    @property
    def atom_types(self):
        return self._data_dict['atom_types']
    
    @property
    def atom_charges(self):
        return self._data_dict['atom_charges']
    
    @property
    def origin(self):
        return self._data_dict['origin']
    
    @property
    def volume(self):
        return np.abs(np.linalg.det(self.cell))
    
    def _save_cube(self, filename: str,
                  data_factor: float = 1.0,
                  box_factor: float = 1.0):
        """Save the grid data to a cube file.
        
        Args:
            filename (str): Path to the output cube file.
            data_factor (float): Factor to multiply the data values (e.g., to convert units).
            box_factor (float): Factor to multiply the cell vectors (e.g., to convert units
        """
        with open(filename, 'w') as f:
            f.write("CUBE FILE\n")
            f.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
            f.write(f"{len(self.atom_types):5d} {self.origin[0] * box_factor:12.6f} {self.origin[1] * box_factor:12.6f} {self.origin[2] * box_factor:12.6f}\n")
            for i in range(3):
                f.write(f"{self.data.shape[i]:5d} {self.cell[i,0] * box_factor/self.data.shape[i]:12.6f} {self.cell[i,1] * box_factor/self.data.shape[i]:12.6f} {self.cell[i,2] * box_factor/self.data.shape[i]:12.6f}\n")
            for i in range(len(self.atom_types)):
                f.write(f"{self.atom_types[i]:5d} {self.atom_charges[i]:12.6f} {self.atom_positions[i,0] * box_factor:12.6f} {self.atom_positions[i,1] * box_factor:12.6f} {self.atom_positions[i,2] * box_factor:12.6f}\n")
            flat_data = self.data.flatten() * data_factor
            for i in range(0, len(flat_data), 6):
                line_data = flat_data[i:i+6] 
                f.write(" ".join(f"{x:17.11e}" for x in line_data) + "\n")
    
    def save_cube(self, filename: str):
        """Save the grid data to a cube file in its original units."""
        self._save_cube(filename, data_factor=1.0, box_factor=1.0)
    
    @staticmethod
    def from_cube(cube_file: str, 
                  data_factor: float = 1.0,
                  box_factor: float = 1.0,
                  ):
        """Create a Cube object from a cube file.
        
        Args:
            cube_file (str): Path to the cube file.
            data_factor (float): Factor to multiply the data values (e.g., to convert units).
            box_factor (float): Factor to multiply the cell vectors (e.g., to convert units).
        """
        assert os.path.isfile(cube_file), f"{cube_file} does not exist"
        
        with open(cube_file, 'r') as f:
            lines = f.readlines()

        natom, x_origin, y_origin, z_origin = map(float, lines[2].split())
        origin =[x_origin, y_origin, z_origin]

        grid_size = []
        cell = []
        for i in range(3):
            sl = lines[3 + i].split()
            grid_size.append(int(sl[0]))
            cell.append([float(x) * int(sl[0]) for x in sl[1:4]])

        atomz, chg, coords = [], [], []
        for line in lines[6:6+int(natom)]:
            atomz.append(int(line.split()[0]))
            chg.append(float(line.split()[1]))
            coords.append(list(map(float, line.split()[2:])))

        data = []
        for line in lines[6+int(natom):]:
            data.extend(list(map(float, line.split())))
        data = np.array(data) 
        data = data.reshape(grid_size)
        
        return Grid(data * data_factor, 
                    np.array(cell) * box_factor, 
                    np.array(coords) * box_factor, 
                    np.array(atomz), 
                    np.array(chg), 
                    np.array(origin) * box_factor)
    
    
    def save_npz(self, filename: str):
        """Save the data to a .npz file"""
        np.savez_compressed(filename, **self._data_dict)
    
    @classmethod
    def from_npz(cls, filename: str):
        """Load the grid data from a .npz file"""
        loaded = np.load(filename, allow_pickle=True)
        return cls(loaded['data'], loaded['cell'], loaded['atom_positions'], loaded['atom_types'], loaded['atom_charges'], loaded['origin'])
    
    @classmethod
    def from_stru(cls, stru_file: str, grid_size: Tuple[int, int, int]):
        """Create a Grid object from a STRU file.
        
        Args:
            stru_file (str): Path to the STRU file.
            grid_size (tuple): Grid size along each cell direction (nx, ny, nz).
        """
        

        stru = AbacusStru.ReadStru(stru_file)
        if stru is None:
            raise ValueError(f"Failed to read STRU file {stru_file}")
        
        cell = stru.get_cell(bohr=False)  # in Angstrom
        coord = stru.get_coord(bohr=False, direct=False)  # in Angstrom
        element = stru.get_element(number=True, total=True)
        
        return cls(np.zeros(grid_size), np.array(cell), np.array(coord), np.array(element), np.zeros(len(element)), np.zeros(3))


    def supercell(self, sc: Tuple[int, int, int]):
        """Create a supercell of the grid data.
        
        Args:
            sc (tuple): Supercell size in each dimension (sx, sy, sz).
        
        Returns:
            Grid: A new Grid object representing the supercell.
        """
        assert len(sc) == 3 and all(isinstance(x, int) and x > 0 for x in sc), "Supercell size should be a tuple of three positive integers"
        
        new_data = np.tile(self.data, sc)
        new_cell = self.cell * np.array(sc)[:, None]
        
        if self.atom_positions is not None and len(self.atom_positions) > 0:
            new_atom_positions = []
            new_atom_types = []
            new_atom_charges = []
            for i in range(sc[0]):
                for j in range(sc[1]):
                    for k in range(sc[2]):
                        shift = i * self.cell[0] + j * self.cell[1] + k * self.cell[2]
                        new_atom_positions.append(self.atom_positions + shift)
                        new_atom_types.extend(self.atom_types)
                        new_atom_charges.extend(self.atom_charges)
            new_atom_positions = np.vstack(new_atom_positions)
        else:
            new_atom_positions = np.zeros(3)
            new_atom_types = np.array([0])
            new_atom_charges = np.array([0.0])
        
        new_origin = self.origin
        
        return Grid(new_data, new_cell, new_atom_positions, new_atom_types, new_atom_charges, new_origin)
    
    def profile1d(self, axis: Literal['a', 'b', 'c'] = 'c', average: bool=False, cartesian: bool=False):
        """Integrate the 3D cube data to 2D plane.
        Args:
            axis (str): the axis to be integrated. 'a' means integrate bc plane, 'b' means ac plane, 'c' means ab plane.
            average (bool): whether to take the average of the integrated values. If False, the integrated values will be summed.
            cartesian (bool): whether to return the coordinates of profile in cartesian coordinates. If False, the profile will be returned in direct coordinates.
        
        Returns:
            tuple: A tuple containing the integrated values and the coordinates of the profile.
        """
        import numpy as np

        func = np.mean if average else np.sum
        if axis == "a":
            val = func(self.data, axis=2) # integrate along c
            val = func(val, axis=1) # integrate along b (axis 1 in integrated val)
        elif axis == "b":
            val = func(self.data, axis=0) # integrate along a
            val = func(val, axis=1) # integrate along c (axis 1 in integrated val)
        elif axis == "c":
            val = func(self.data, axis=0) # integrate along a
            val = func(val, axis=0) # integrate along b (axis 0 in integrated val)
        else:
            raise ValueError(f"Invalid axis: {axis} is not a, b or c")
        
        ngrid = self.data.shape[0] if axis == "a" else self.data.shape[1] if axis == "b" else self.data.shape[2]
        if cartesian:
            vec_length = np.linalg.norm(self.cell[0]) if axis == "a" else np.linalg.norm(self.cell[1]) if axis == "b" else np.linalg.norm(self.cell[2])
            coord = np.linspace(0, vec_length, ngrid)
        else:
            coord = np.linspace(0, 1, ngrid)

        return val, coord

class Charge(Grid):
    """Subclass for charge density data.
    
    By using from_cube method, the unit of charge density will be converted to e/Ang^3,
    and the unit of cell and positions will be converted to Angstrom.
    
    For example:
    >>> chg = Charge.from_cube("charge.cube")
    >>> print(chg.data.shape)
    (51, 51, 51)
    >>> print(chg.cell)
    [[ 5.00000000e+01  0.00000000e+00  0.00000000e+00]
     [ 0.00000000e+00  5.00000000e+01  0.00000000e+00]
     [ 0.00000000e+00  0.00000000e+00  5.00000000e+01]]
    >>> supercell = chg.supercell((2, 2, 2)) # create a supercell of size (2, 2, 2)
    >>> print(supercell.data.shape)
    (102, 102, 102)
    >>> supercell.save_cube("supercell.cube") # save the supercell to a cube file

    """
    def __init__(self, 
                 data: np.ndarray,
                 cell: np.ndarray,
                 atom_positions: Optional[np.ndarray] = None,
                 atom_types: Optional[List[str]] = None,
                 atom_charges: Optional[List[float]] = None,
                 origin: np.ndarray = np.zeros(3)
                 ):
        super().__init__(data, cell, atom_positions, atom_types, atom_charges, origin)

    @staticmethod        
    def from_cube(cube_file: str, format: str = "abacus"):
        """Load charge density from a cube file"""
        
        if format == "abacus":
            data_factor = 1 / BOHR2A**3  # ABACUS charge density is in e/Bohr^3, convert to e/Ang^3
            box_factor = BOHR2A  # ABACUS cell is in Bohr, convert to Angstrom
        else:
            raise ValueError(f"Unsupported format {format}. Supported formats are 'abacus'.")
            
        cube = Grid.from_cube(cube_file, data_factor, box_factor)
        return Charge(cube.data, cube.cell, cube.atom_positions, cube.atom_types, cube.atom_charges, cube.origin)
    
    def save_cube(self, filename, format: str = "abacus"):
        """Save the charge density to a cube file"""
        if format == "abacus":
            data_factor = BOHR2A**3  # Convert back to e/Bohr^3
            box_factor = 1 / BOHR2A  # Convert back to Bohr
            return self._save_cube(filename, data_factor, box_factor)
        else:   
            raise ValueError(f"Unsupported format {format}. Supported formats is 'abacus'.")
    
    def to_pot(self):
        """Solve the Poisson equation to get the electrostatic potential from the charge density.
        Returns:
            Potential: The electrostatic potential object.
        """
        pot = chg2pot(self.data, self.cell)
        return Potential(pot, self.cell, self.atom_positions, self.atom_types, self.atom_charges, self.origin)  
    
    def supercell(self, sc: Tuple[int, int, int]):
        """Create a supercell of the charge density data.
        
        Args:
            sc (tuple): Supercell size in each dimension (sx, sy, sz).
        
        Returns:
            Charge: A new Charge object representing the supercell.
        """
        grid = super().supercell(sc)
        return Charge(grid.data, grid.cell, grid.atom_positions, grid.atom_types, grid.atom_charges, grid.origin)

class Potential(Grid):
    """Subclass for electrostatic potential data.
    
    By using from_cube method, the unit of potential will be converted to eV,
    and the unit of cell and positions will be converted to Angstrom.
    
    For example:
    >>> pot = Potential.from_cube("potential.cube")
    >>> print(pot.data.shape)
    (51, 51, 51)
    >>> supercell = pot.supercell((2, 2, 2))
    >>> print(supercell.data.shape)
    (102, 102, 102)
    >>> supercell.save_cube("supercell.cube")
    """

    def __init__(
        self,
        data: np.ndarray,
        cell: np.ndarray,
        atom_positions: Optional[np.ndarray] = None,
        atom_types: Optional[List[str]] = None,
        atom_charges: Optional[List[float]] = None,
        origin: np.ndarray = np.zeros(3),
    ):
        super().__init__(data, cell, atom_positions, atom_types, atom_charges, origin)

    @staticmethod
    def from_cube(cube_file: str, format: str = "abacus"):
        """Load potential from a cube file"""

        if format == "abacus":
            data_factor = (
                -1 * RY2EV
            )  # ABACUS potential is in Ry, convert to eV. In abacus the electron is positive, so we need a negative sign here.
            box_factor = BOHR2A  # ABACUS cell is in Bohr, convert to Angstrom
        else:
            raise ValueError(f"Unsupported format {format}. Supported formats is 'abacus'.")
        
        cube = Grid.from_cube(cube_file, data_factor, box_factor)
        return Potential(
            cube.data,
            cube.cell,
            cube.atom_positions,
            cube.atom_types,
            cube.atom_charges,
            cube.origin,
        )

    @staticmethod
    def from_locpot(locpot_file: str):
        """Read the local potential from a LOCPOT file written by VASP"""
        from ase.io import read
        import uuid

        with open(locpot_file, "r") as f:
            lines = [line.strip() for line in f if line.strip() != "" or line == "\n"]

        # Get total number of atoms
        atom_nums = [int(x) for x in lines[6].split()]
        total_atom_nums = sum(atom_nums)

        # Dump POSCAR file in the head of LOCPOT
        end_line_idx = total_atom_nums + 7
        if lines[7].lower == "selective dynamics":
            end_line_idx += 1

        poscar_dump = "_dumped_POSCAR"
        while os.path.exists(poscar_dump):
            poscar_dump = "_dumped_POSCAR_" + str(uuid.uuid4())[:8]

        with open(poscar_dump, "w") as f:
            for line_idx in range(end_line_idx + 1):
                f.write(lines[line_idx] + "\n")

        pos = read(poscar_dump, format="vasp")
        os.unlink(poscar_dump)

        # Read size of gird data
        line_idx = end_line_idx + 1
        grid = tuple(int(x) for x in lines[line_idx].split())
        nx, ny, nz = grid
        total_grid_points = nx * ny * nz
        line_idx += 1

        def read_grid_data(start_idx: int, n_points: int) -> np.ndarray:
            values = []
            idx = start_idx
            while len(values) < n_points:
                if idx >= len(lines):
                    raise RuntimeError("No sufficient data in LOCPOT file")
                values.extend([float(x) for x in lines[idx].split()])
                idx += 1
            if len(values) != n_points:
                raise RuntimeError( f"expected {n_points} data points, read {len(values)} data points")
            # reshape data - VASP stores data in z, y, x order (fastest to slowest: x, y, z)
            return (np.array(values).reshape((nz, ny, nx)).transpose(2, 1, 0))  # Convert to (nx, ny, nz)

        data_first = read_grid_data(line_idx, total_grid_points)
        line_idx += int(np.ceil(total_grid_points / 5))

        # check if there are any remaining lines
        remaining_lines = len(lines) - line_idx
        min_spin_lines = (1 + int(np.ceil(total_atom_nums / 5)) + int(np.ceil(total_grid_points / 5)))
        is_spin_polarized = remaining_lines >= min_spin_lines

        if is_spin_polarized:
            # skip atom count lines of "1" (5 per line)
            spin_marker_lines = int(np.ceil(total_atom_nums / 5))
            line_idx += spin_marker_lines

            grid2 = tuple(int(x) for x in lines[line_idx].split())
            if grid2 != grid:
                raise Warning(f"Mesh size of second data is {grid2}, not same with first data ({grid})")
            line_idx += 1

            data_second = read_grid_data(line_idx, total_grid_points)
            data = data_first + data_second
        else:
            data = data_first

        # Convert to potential of electrons
        data *= -1

        # Create atom_charges (default to 0) and origin (default to [0, 0, 0])
        atom_charges = np.array([0.0] * total_atom_nums)
        origin = np.zeros(3)

        return Potential(data, np.array(pos.get_cell()), pos.get_positions(), pos.get_atomic_numbers(), atom_charges, origin)

    def save_cube(self, filename, format: str = "abacus"):
        """Save the potential to a cube file"""
        if format == "abacus":
            data_factor = -1 / RY2EV  # Convert back to Ry
            box_factor = 1 / BOHR2A  # Convert back to Bohr
            return self._save_cube(filename, data_factor, box_factor)
        else:   
            raise ValueError(f"Unsupported format {format}. Supported formats is 'abacus'.")
        
    def supercell(self, sc: Tuple[int, int, int]):
        """Create a supercell of the potential data.
        
        Args:
            sc (tuple): Supercell size in each dimension (sx, sy, sz).
        
        Returns:
            Potential: A new Potential object representing the supercell.
        """
        grid = super().supercell(sc)
        return Potential(grid.data, grid.cell, grid.atom_positions, grid.atom_types, grid.atom_charges, grid.origin)

        
        
        