from typing import List, Dict, Union
import os,sys,traceback, re
import numpy as np
from pathlib import Path
from .. import constant

class AbacusStru:
    def __init__(self,
                 label:List[str],
                 cell:List[List[float]] ,
                 coord:List[List[float]] ,
                 pp:List[str] = None,
                 atom_number:List[int] = None,
                 orb:List[str] = None,
                 paw:List[str] = None,
                 mass:List[float] = None,
                 element:List[str] = None,
                 lattice_constant:float = 1,
                 move:List[List[int]] =None,
                 magmom: List[float] = None,
                 magmom_atom: List[Union[float,List[float]]] = None,
                 velocity: List[List[float]] = None,
                 angle1: List[float] = None,
                 angle2: List[float] = None,
                 dpks:str = None,
                 cartesian: bool = False):    
        """ABACUS STRU class, the unit is Bohr

        Parameters
        ----------
        label : List[str]
            the label name of each type; can has natom number or n atom type number.
            if has natom number, then the atom_number will be ignored, and the pp/orb is the order of the first label.
            if has n atom type number, then the atom_number will be the number of each type, and the pp/orb is the order of the label, and the coord is the order of the label.
        atom_number : List[int]
            the atom number of each type
        pp : List[str]
            pseudopotential file of each type
        orb : List[str], optional
            orbital file of each type, by default None
        paw: List[str], optional
            paw file of each type, by default None. Should has pp or paw.
        mass : List[float], optional
            mass of each type, by default None
        element : List[str], optional
            element name of each type, by default None
        cell : List[List[float]]
            cell of a, b, c
        coord : List[List[float]]
            coordinate of each atom
        lattice_constant : float, optional
            the lattice constance, by default 1
        move: List[List[int]], optional
            if do move of each atom, by default None
        magmom : float, optional
            if the number of magmom equal to type number, then set it to each type, by default None
            else if the number of magmom equal to coord number, then set it to each atom, by default None
        magmom_atom: List[Union[float,List[float]]], optional
            if the magmom_atom is not None, then set the magmom of each atom, by default None
        velocity: List[List[float]], optional
            the velocity of each atom, by default None
        angle1: List[float], optional
            the angle1 of each atom, by default None
        angle2: float, optional 
            the angle2 of each atom, by default None
        dpks : str, optional
            the deepks descriptor file name, by default None
        cartesian : bool, optional
            if the coordinate is cartesian type, default is direct type, by default False
        """        
        #check if label number is equal to coord number
        if len(label) == len(coord):
            self._label = []
            for i in label:
                if i not in self._label:
                    self._label.append(i)
            self._atom_number = [label.count(i) for i in self._label]
            # should sort coord to the order of label
            coords = [[] for i in self._label]
            for i,ilabel in enumerate(label):
                coords[self._label.index(ilabel)].append(coord[i])
            self._coord = []
            for i in coords:
                self._coord += i
        else:
            assert(atom_number != None), "ERROR: label number is not equal to coord number, but atom_number is not defined "
            assert(len(label) == len(atom_number)), "ERROR: label number is not equal to atom_number number"
            self._label = label
            self._atom_number = atom_number
            self._coord = coord

        total_atom = np.array(self._atom_number).sum()
        assert(total_atom == len(self._coord)), f"ERROR: the total atom number is not equal to coord number, {total_atom} != {len(coord)}"
        self._cell = cell
        
        # check pp orb paw
        self._pp = pp if pp else None
        self._orb = orb if orb else None
        self._paw = paw if paw else None
        
        if (not self._pp):
            print("WARNING: pp is not defined!!!")
            self._pp = ["" for i in range(len(self._label))]
            
        
        self._lattice_constant = lattice_constant
        self._dpks = dpks if dpks else None
        self._cartesian = cartesian
        self._move = move
        self._magmom = magmom if magmom else [0]*len(label)
        self._magmom_atom = magmom_atom
        self._velocity = velocity
        self._angle1 = angle1
        self._angle2 = angle2
        
        if element != None:
            self._element = element
        else:
            #print("WARNING: element is not defined, will use the first one/two letter of label as element name")
            self._element = []
            for i in self._label:
                ele = i[0]
                if len(i) > 1 and i[1].islower():
                    ele += i[1]
                self._element.append(ele)
        
        if mass != None:
            self._mass = mass
        else:
            self._mass = [constant.MASS_DICT.get(i,1.0) for i in self._element]
        
        self._check()
    
    def _check(self):
        '''check the structure'''
        n_label = len(self._label)
        n_atom = len(self._coord)
        
        assert(n_label == len(self._atom_number)), "ERROR: the length of label is not equal to the length of atom_number"
        assert(n_atom == sum(self._atom_number)), "ERROR: the length of coord is not equal to sum of atom_number"
        
        if self._pp != None:
            assert(n_label == len(self._pp)), "ERROR: the length of label is not equal to pp"
        if self._orb != None:
            assert(n_label == len(self._orb)), "ERROR: the length of label is not equal to orb"
        if self._paw != None:
            assert(n_label == len(self._paw)), "ERROR: the length of label is not equal to paw"
        if self._move:
            assert(n_atom == len(self._move)), "ERROR: the length of coord is not equal to move"
        if self._magmom:
            assert(n_label == len(self._magmom)), "ERROR: the length of label is not equal to magmom"
        if self._magmom_atom:
            assert(n_atom == len(self._magmom_atom)), "ERROR: the length of coord is not equal to magmom_atom"
        if self._velocity:
            assert(n_atom == len(self._velocity)), "ERROR: the length of coord is not equal to velocity"
        if self._angle1:
            assert(n_atom == len(self._angle1)), "ERROR: the length of coord is not equal to angle1"
        if self._angle2:
            assert(n_atom == len(self._angle2)), "ERROR: the length of coord is not equal to angle2"
        return True    
    
    def get_pp(self):
        return self._pp

    def get_orb(self):
        return self._orb  
    
    def get_paw(self):
        return self._paw  

    def get_label(self,total=True):
        '''return the label name of each atom'''
        if total:
            label = []
            for idx,i in enumerate(self._atom_number):
                label += [self._label[idx]] * i
            return label
        else:
            return self._label
    
    def get_dpks(self):
        return self._dpks
    
    def get_mass(self):
        return self._mass

    def get_atommag(self):
        # return the magmom of each atom
        magmom = []
        for i,ilabel in enumerate(self._label):
            magmom += [self._magmom[i]] * self._atom_number[i]
        if self._magmom_atom:
            for ii,i in enumerate(self._magmom_atom):
                if i != None:
                    magmom[ii] = i
        return magmom
    
    def get_mag(self):
        # return the magmom of each type
        return self._magmom
    
    def get_element(self,number=True,total=True):
        '''return the element name of each atom'''
        if total:
            element = []
            for idx,i in enumerate(self._atom_number):
                element += [self._element[idx]] * i
        else:
            element = self._element
            
        if not number:
            return element
        else:
            return [constant.PERIOD_DICT_NUMBER[i] for i in element]
    
    def get_cell(self,bohr = False):
        "return the cell matrix, in unit of Angstrom"
        transfer_unit = 1 if bohr else constant.BOHR2A
        cell = np.array(self._cell) * self._lattice_constant * transfer_unit
        return cell.tolist()
    
    def get_coord(self,bohr = False, direct=False):
        '''return the coordinate matrix, in cartesian or direct type, in unit of Angstrom or Bohr.
        If direct is True, then return the direct type coordinate, and bohr will be ignored.'''
        transfer_unit = 1 if bohr else constant.BOHR2A
        if self._cartesian:
            if not direct:
                coord = np.array(self._coord) * transfer_unit * self._lattice_constant
                return coord.tolist()
            else:
                coord = np.array(self._coord)
                return coord.dot(np.linalg.inv(np.array(self._cell))).tolist()
        else:
            if not direct:
                coord = np.array(self._coord) * transfer_unit * self._lattice_constant
                return coord.dot(np.array(self._cell)).tolist()
            else:
                return self._coord
    
    def get_stru(self):
        return {
            "cell": self._cell,
            "coord": self._coord,
            "lat": self._lattice_constant,
            "cartesian": self._cartesian,
        }
    
    def get_kline(self,orig_cell=False,
                  new_stru_file=None,
                  kpt_file=None,
                  point_number=20,
                  with_time_reversal=True,
                  recipe="hpkot", 
                  threshold=1e-7, 
                  symprec=1e-5, 
                  angle_tolerance=-1.0):
        '''
        Use seekpath to generate the k-line: https://github.com/giovannipizzi/seekpath
        orig_cell: if True, then use the input cell to generate the k-line, else use primitive cell
        return a new_stru, point_coords, path, kpath

            "atom_label": label, # a list of the label of each atom
            "point_coords": {"Gamma": [0,0,0], "X": [0.5,0,0], ...}, # a dict of the coordinate of each high symmetry point
            "path": [("Gamma","X"),("X","M"),...], # a list of the k-line, each element is a tuple of two high symmetry point
        
        If orig_cell = False, then the cell and coord is the primitive cell, else the input cell.
        
        Point_number: the number of k-point between two high symmetry point, 
        And if the high symmetry of previous and next point is not same, then the number will be 1.
        
        with_time_reversal/recipe/threshold/symprec/angle_tolerance: the parameter of seekpath
        '''
        import seekpath
        cell = self.get_cell(bohr=True)
        coord = self.get_coord(bohr=True,direct=True)
        labels = self.get_label()
        label = []
        for i in labels: 
            if i not in label:
                label.append(i)
        number = [label.index(i) for i in labels]
        stru = (cell,coord,number)
        
        while symprec <= 1e-2:
            try:
                if orig_cell:
                    kpath = seekpath.get_path_orig_cell(stru,with_time_reversal=with_time_reversal,recipe=recipe,threshold=threshold,symprec=symprec,angle_tolerance=angle_tolerance)
                else:
                    kpath = seekpath.get_path(stru,with_time_reversal=with_time_reversal,recipe=recipe,threshold=threshold,symprec=symprec,angle_tolerance=angle_tolerance)
                break
            except:
                traceback.print_exc()
                print("WARNING: get_path failed, increase symprec to %e" % (symprec*10))
                symprec *= 10
        
        if symprec > 1e-2:
            print("ERROR: get_path failed, please check the structure")
            sys.exit(1)
            
        if orig_cell:
            new_cell = cell
            new_coord = coord
            new_label = labels
        else:
            new_cell = kpath["primitive_lattice"]
            # resort the new_coord and new_label to the order of label
            new_label = []
            new_coord = []
            for i in range(len(label)):
                for j in range(len(kpath["primitive_types"])):
                    if kpath["primitive_types"][j] == i:
                        new_label.append(label[i])
                        new_coord.append(kpath["primitive_positions"][j])
        
        point_coords = kpath["point_coords"]
        path = kpath["path"]
        
        # transfer GAMMA to G
        if "G" not in point_coords:
            if "GAMMA" in point_coords:
                point_coords["G"] = point_coords["GAMMA"]
                del point_coords["GAMMA"]
                
                for idx,ipath in enumerate(path):
                    if ipath[0] == "GAMMA":
                        path[idx] = ("G",ipath[1])
                    if ipath[1] == "GAMMA":
                        path[idx] = (ipath[0],"G")
        
        # generate the k-line
        if kpt_file != None:
            kpt = [point_coords[path[0][0]] + [point_number,path[0][0]]]
            for idx,ipath in enumerate(path[:-1]):
                if ipath[1] == path[idx+1][0]:
                    kpt.append(point_coords[ipath[1]] + [point_number,ipath[1]])
                else:
                    kpt.append(point_coords[ipath[1]] + [1,ipath[1]])
                    kpt.append(point_coords[path[idx+1][0]] + [point_number,path[idx+1][0]])
            kpt.append(point_coords[path[-1][1]] + [1,path[-1][1]])        
            WriteKpt(kpt,kpt_file,model="line")
        
        # write the new stru file
        if new_stru_file != None and not orig_cell:
            new_stru = AbacusStru(label=new_label,
                                  cell=new_cell,
                                  coord=new_coord,
                                  pp=self._pp,
                                  orb=self._orb,
                                  paw = self._paw,
                                  lattice_constant=1.0,
                                  magmom=self._magmom,
                                  cartesian=False)
            new_stru.write(new_stru_file)
        else:
            new_stru = self
        
        return new_stru,point_coords,path,kpath   

    def get_kline_ase(self, point_number=10, kpt_file=None):
        from ase import Atoms
        from ase.cell import Cell
        from ase.geometry.dimensionality import (analyze_dimensionality,isolate_components)
        cell = self.get_cell(bohr=False)
        coord =[tuple(i) for i in self.get_coord(bohr=False,direct=True)]
        labels = self.get_label()
        atoms = Atoms(symbols=labels, scaled_positions=coord, cell=cell, pbc=True)
        intervals = analyze_dimensionality(atoms,method='RDA')
        kpath = None
        if intervals[0].dimtype == "3D":
            kpath = atoms.cell.bandpath(npoints=100)
        elif intervals[0].dimtype == "2D":
            lat_length0 = atoms.cell.lengths()
            lat_length1 = []
            result = isolate_components(atoms,kcutoff=1.5)
            for dim,components in result.items():
                for atoms in components:
                    lat_length1 = atoms.cell.lengths()
            for i in range(len(lat_length0)):
                lat_length0[i] = round(lat_length0[i],5)
                lat_length1[i] = round(lat_length1[i],5)
            pbc = []
            for i in range(len(lat_length0)):
                if lat_length0[i] in lat_length1:
                    pbc.append(1)
                else:
                    pbc.append(0)
            kpath = lat2d_pbc.bandpath(npoints=100)

        if kpath:
            pattern = '([A-Za-z]\\d?|,)'  # pattern to match points of form like 'A' or 'B1'
            kpt_new = []
            # Split path by high symmetry point names allowing for points like 'B1'
            path_elements = re.findall(pattern, kpath.path)
            #print(path_elements)
            insert_details = []
            for i, point in enumerate(path_elements):
                if point == ',':
                    continue
                elif point not in kpath.special_points:
                    raise ValueError(f"The point {point} is not defined in your special_points dictionary.")
                else:
                    #print(path_elements[i])
                    #print(kpath.special_points[point])
                    kpt_new.append(kpath.special_points[point])
                    if i == len(path_elements) - 1 or path_elements[i+1] == ',':
                        insert_details.append(str(1)+" # "+path_elements[i])
                    else:
                        insert_details.append(str(point_number)+" # "+path_elements[i])

            lines = ['K_POINTS\n', str(len(kpt_new))+'\n', 'Line\n']
            lines.extend(f"{k[0]:.12f} {k[1]:.12f} {k[2]:.12f} {mult}\n" for k, mult in zip(kpt_new, insert_details))
            with open(kpt_file, 'w') as f:
                f.write(''.join(lines))
        else:
            raise ValueError(f"kpath not found")

    def set_pp(self,pplist):
        if pplist and len(pplist) != len(self._label):
            print("ERROR: the length of pplist is not equal to label number")
            sys.exit(1)
        self._pp = pplist if pplist else None
    
    def set_orb(self,orblist):
        if orblist and len(orblist) != len(self._label):
            print("ERROR: the length of orblist is not equal to label number")
            sys.exit(1)
        self._orb = orblist if orblist else None
    
    def set_paw(self,pawlist):
        if pawlist and len(pawlist) != len(self._label):
            print("ERROR: the length of pawlist is not equal to label number")
            sys.exit(1)
        self._paw = pawlist if pawlist else None

    def set_dpks(self,descriptor):
        self._dpks = descriptor if descriptor else None
    
    def set_mass(self,mass):
        if len(mass) != len(self._label):
            print("ERROR: the length of mass is not equal to label number")
            sys.exit(1)
        self._mass = mass
        
    def set_element(self,element):
        self._element = element
    
    def set_coord(self,coord,direct=False,bohr=True):
        '''
        set the coordinate of each atom
        
        if direct is True, then the coord is direct type, and will modify self._cartesian to False
        else the coord is cartesian type, and will modify self._cartesian to True, and will set lattice_constant to 1.0, and modify cell *= lattice_constant
        
        if bohr is False, then will transfer the coord to Bohr unit
        '''
        if direct:
            self._cartesian = False
            self._coord = coord
        else:
            self._cartesian = True
            self._cell = np.array(self._cell) * self._lattice_constant
            self._lattice_constant = 1.0
            if bohr:
                self._coord = coord
            else:
                self._coord = np.array(coord) / constant.BOHR2A
                self._coord = self._coord.tolist()
        self._check()
    
    def set_cell(self,cell,bohr=True,change_coord=True):
        '''
        set the lattice of a, b, c
        
        Will set the lattice_constant to 1.0
        
        if bohr is False, then will transfer the cell to Bohr unit
        
        if change_coord is True, then will set the coord to direct type
        else will firstly transfer the coord to cartesian type.
        '''
        if change_coord:
            if self._cartesian:
                coord = self.get_coord(bohr=bohr,direct=True)
                self._coord = coord
                self._cartesian = False
                
        self._lattice_constant = 1.0
        if bohr:
            self._cell = cell
        else:
            self._cell = np.array(cell) / constant.BOHR2A
            self._cell = self._cell.tolist()
        self._check()        
        
     
    def write(self,struf="STRU"):
        cc = ""
        #write species
        cc += "ATOMIC_SPECIES\n"
        for i,ilabel in enumerate(self._label):
            if self._pp :
                cc += "%s %f %s\n" % (ilabel,self._mass[i],self._pp[i])
            else:
                cc += "%s %f\n" % (ilabel,self._mass[i])
        
        #write orb
        if self._orb:
            cc += "\nNUMERICAL_ORBITAL\n"
            for i in self._orb:
                cc += i + "\n"
        
        # write pawfile
        if self._paw:
            cc += "\nPAW_FILES\n"
            for i in self._paw:
                cc += i + "\n"
        
        #write  LATTICE_CONSTANT
        cc += "\nLATTICE_CONSTANT\n%f\n" % self._lattice_constant

        #write LATTICE_VECTORS
        cc += "\nLATTICE_VECTORS\n"
        for i in self._cell:
            cc += "%17.11f %17.11f %17.11f\n" % tuple(i)
        
        #write ATOMIC_POSITIONS
        cc += "\nATOMIC_POSITIONS\n"
        if self._cartesian:
            cc += "Cartesian\n"
        else:
            cc += "Direct\n"
        icoord = 0
        for i,ilabel in enumerate(self._label):
            cc += "\n%s\n%f\n%d\n" % (ilabel,self._magmom[i],self._atom_number[i])
            for j in range(self._atom_number[i]):
                cc += "%17.11f %17.11f %17.11f " % tuple(self._coord[icoord + j])
                if self._move:
                    cc += "%d %d %d " % tuple(self._move[icoord + j])
                if self._magmom_atom:
                    if isinstance(self._magmom_atom[icoord + j],list):
                        if len(self._magmom_atom[icoord + j]) == 3:
                            cc += "mag %f %f %f " % tuple(self._magmom_atom[icoord + j])
                        elif len(self._magmom_atom[icoord + j]) == 1:
                            cc += "mag %f " % self._magmom_atom[icoord + j][0]
                    elif self._magmom_atom[icoord + j] != None:
                        cc += "mag %f " % self._magmom_atom[icoord + j]
                if self._velocity:
                    if self._velocity[icoord + j]:
                        cc += "v %f %f %f " % tuple(self._velocity[icoord + j])
                if self._angle1:
                    if self._angle1[icoord + j]:
                        cc += "angle1 %f " % self._angle1[icoord + j]
                if self._angle2:
                    if self._angle2[icoord + j]:
                        cc += "angle2 %f " % self._angle2[icoord + j]
                cc += "\n"
            icoord += self._atom_number[i]
        
        #write dpks
        if self._dpks:
            cc += "\nNUMERICAL_DESCRIPTOR\n"
            cc += self._dpks    

        Path(struf).write_text(cc)  
    
    def write2poscar(self,poscar="POSCAR"):
        '''
        write to POSCAR file for VASP
        '''
        cc = "STRUCTURE translated by abacustest\n"
        cc += f"{self._lattice_constant * constant.BOHR2A}\n"
        for i in self._cell:
            cc += "%17.11f %17.11f %17.11f\n" % tuple(i)
        cc += " ".join(self._label)+"\n"
        cc += " ".join([str(i) for i in self._atom_number])+"\n"
        if self._cartesian:
            cc += "Cartesian\n"
        else:
            cc += "Direct\n"
        for i in self._coord:
            cc += "%17.11f %17.11f %17.11f\n" % tuple(i)
            
        if poscar != None:
            Path(poscar).write_text(cc)
        return cc
    
    @staticmethod
    def parse_stru_pos(pos_line):
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
      0.0 0.0 0.0 m 0 0 0 mag 1.0 angle1 90 angle2 0
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
        if len(sline) > 3:
            mag_list = []
            velocity_list = []
            move_list = []
            angle1_list = []
            angle2_list = []
            label = "move"
            for i in range(3,len(sline)):
                if sline[i] == "m":
                    label = "move"
                elif sline[i] in ["v","vel","velocity"]:
                    label = "velocity"
                elif sline[i] in ["mag","magmom"]:
                    label = "magmom"
                elif sline[i] == "angle1":
                    label = "angle1"
                elif sline[i] == "angle2":
                    label = "angle2"
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
                    
            if len(move_list) == 3:
                move = move_list
            if len(velocity_list) == 3:
                velocity = velocity_list
            if len(mag_list) in [1,3]:
                magmom = mag_list
            if len(angle1_list) == 1:
                angle1 = angle1_list[0]
            if len(angle2_list) == 1:
                angle2 = angle2_list[0]
                
        return pos,move,velocity,magmom,angle1,angle2

    @staticmethod
    def ReadStru(stru:str = "STRU"):
        '''
        read the label, pp, orb, cell, coord, deepks-descriptor
        

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
                        elif lines[ij].strip() in constant.ABACUS_STRU_KEY_WORD:
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
        coord_type = atom_positions[0].split("#")[0].strip().lower()
        if coord_type.startswith("dire"):
            cartesian = False
        elif coord_type.startswith("cart"):
            cartesian = True
        else:
            print("Not support coordinate type %s now." % atom_positions[0].strip())
            sys.exit(1)
        i = 1
        while i < len(atom_positions):
            label = atom_positions[i].strip()
            if label not in labels:
                print("label '%s' is not matched that in ATOMIC_SPECIES" % label)
                sys.exit(1)
            magmom_global.append(float(atom_positions[i+1].split()[0]))
            atom_number.append(int(atom_positions[i+2].split()[0]))
            i += 3
            for j in range(atom_number[-1]):
                pos,imove,ivelocity,imag,iangle1,iangle2 = AbacusStru.parse_stru_pos(atom_positions[i+j])
                coords.append(pos)
                move.append(imove)
                velocity.append(ivelocity)
                magmom.append(imag)
                angle1.append(iangle1)
                angle2.append(iangle2)
            i += atom_number[-1]
        
        return AbacusStru(label=labels,
                          atom_number=atom_number,
                          cell=cell,
                          coord=coords,
                          pp=pp,
                          orb=orb,
                          paw = paw,
                          lattice_constant=lattice_constant,
                          move=move,
                          magmom=magmom_global,
                          magmom_atom=magmom,
                          velocity=velocity,
                          angle1=angle1,
                          angle2=angle2,
                          dpks=dpks,
                          cartesian=cartesian)
        

def WriteKpt(kpoint_list:List = [1,1,1,0,0,0],file_name:str = "KPT", model="gamma"):
    '''
    Docs for KPT file: https://abacus.deepmodeling.com/en/latest/advanced/input_files/kpt.html
    ABACUS KPT support three models:
    - gamma/mp: is the Monkhorst-Pack method, such as:
    
        K_POINTS //keyword for start
        0 //total number of k-point, `0' means generate automatically
        Gamma //which kind of Monkhorst-Pack method, `Gamma' or `MP'
        2 2 2 0 0 0 //first three number: subdivisions along recpri. vectors
                    //last three number: shift of the mesh
        
    - direct/cartessian: set the k-point explicitly, such as:
        
            K_POINTS
            1
            Direct
            0.0 0.0 0.0 1.0  // the last number is the weight of this k-point
            
    - line: set the k-point along a line, such as:
            
                K_POINTS
                2  // number of high symmetry k-points along the
                Line
                0.0 0.0 0.0 10  // Gamma the last number is number of k-points between this and next k-point
                0.5 0.0 0.0 1 
    
    For gamma/mp model, the k-point_list should be a list of 6 values, such as:
    [2,2,2,0,0,0]
    
    For explicitly model, the k-point_list should be a list of list of 3 or 4 values, such as:
    [[0.0,0.0,0.0,1.0],[0.5,0.0,0.0,1.0]]
    
    For line model, the k-point_list should be a list of list of 4 or 5 values (the last value is a string of comment), such as:
    [[0.0,0.0,0.0,10],[0.5,0.0,0.0,1]] or
    [[0.0,0.0,0.0,10 "#Gamma"],[0.5,0.0,0.0,1,"//"],[0.5,0.5,0.0,1,"//"]]         
    '''
    if model.lower() in ["gamma","mp"]:
        with open(file_name,'w') as f1:
            f1.write("K_POINTS\n0\nGamma\n")
            f1.write(" ".join([str(i) for i in kpoint_list]))
    elif model.lower() in ["direct","cartessian"]:
        # normalize the weight
        kpt = []
        if len(kpoint_list[0]) == 3:
            kpt = [i+[1.0/len(kpoint_list)] for i in kpoint_list]
        elif len(kpoint_list[0]) == 4:
            total_weight = sum([i[3] for i in kpoint_list])
            kpt = [i[:3]+[i[3]/total_weight] for i in kpoint_list]
        else:
            print(f"ERROR: model is {model}, the kpoint_list is not correct!!!\n{kpoint_list}")
            sys.exit(1)
            
        with open(file_name,'w') as f1:
            f1.write(f"K_POINTS\n{len(kpoint_list)}\n{model.capitalize()}\n")
            for i in kpt:
                f1.write("%17.11f %17.11f %17.11f %17.11f\n" % tuple(i))
    elif model.lower() in ["line"]:
        with open(file_name,'w') as f1:
            f1.write(f"K_POINTS\n{len(kpoint_list)}\nLine\n")
            for i in kpoint_list:
                if len(i) == 4:
                    f1.write("%17.11f %17.11f %17.11f %4d\n" % tuple(i))
                elif len(i) == 5:
                    if not (i[-1].startswith("#") or i[-1].startswith("//")):
                        i[-1] = "#"+i[-1]
                    f1.write("%17.11f %17.11f %17.11f %4d %s\n" % tuple(i))
    else:
        print(f"ERROR: model is {model}, not support now!!!")
        sys.exit(1)


def ReadInput(INPUTf: str = None, input_lines: str = None) -> Dict[str,any]:
    #read the INPUT file
    input_context = {}
    if INPUTf != None:
        if not os.path.isfile(INPUTf):
            print("Can not find the file '%s'" % INPUTf)
            return input_context
        with open(INPUTf) as f1: input_lines = f1.readlines()
    if input_lines == None:
        print(INPUTf)
        print("Please provide the INPUT file name of INPUT lines")
        return input_context
    
    def str2intfloat(ii):
        try:
            return int(ii)
        except:
            pass
        try:
            return float(ii)
        except:
            return ii
    
    readinput = False
    for i,iline in enumerate(input_lines):
        if iline.strip()[:16] == 'INPUT_PARAMETERS':
            readinput = True
        elif iline.strip() == '' or iline.strip()[0] in ['#']:
            continue
        elif readinput:
            sline = re.split('[ \t]',iline.split("#")[0].strip(),maxsplit=1)
            if len(sline) == 2:
                input_context[sline[0].lower().strip()] = str2intfloat(sline[1].strip())
    return input_context

def WriteInput(input_context:Dict[str,any],
               INPUTf: str = "INPUT"):
    out = "INPUT_PARAMETERS\n"
    for k,v in input_context.items():
        if v != None:
            out += "%s\t%s\n" % (str(k),str(v))
        else:
            out += "#%s\t \n" % (str(k))
    with open(INPUTf,'w') as f1: f1.write(out)        