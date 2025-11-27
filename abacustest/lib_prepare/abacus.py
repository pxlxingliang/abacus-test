from typing import List, Dict, Union, Tuple, Any
import os,sys,traceback, re,copy, json, shutil
import numpy as np
from pathlib import Path
from .. import constant
from . import comm
from abacustest.lib_prepare.stru import read_stru_file, write_stru_file, write_poscar

def gen_stru(stru_files, stru_type, pp_path, orb_path, tpath = ".", copy_pp_orb=False):
    """
    Generate the structure files for ABACUS.
    
    Args:
        stru_files (list): List of structure files.
        stru_type (str): Type of the structure files.
        pp_path (str): Path to the pseudopotential files.
        orb_path (str): Path to the orbital files.
    
    Returns:
        jobs (dict): Dictionary, key is the job path, value is:
            {
                "element": element,
                "pp": pp,
                "orb": orb,
                "recommand_ecutwfc":List[float|None] # the recommended ecutwfc for each element
            }
    """
    # Translate structure files to ABACUS STRU format
    stru_paths = comm.translate_strus(stru_files, stru_type, output_path = tpath)
    if stru_paths is None:
        print("Error: tanslate structure files failed.")
        return None
    print("Structures are translated to ABACUS STRU and saved in",stru_paths)
    
    # Find the pseudopotential files and orbital files
    pp_paths = comm.collect_pp(pp_path) # the return file name is pp_path/pp
    orb_paths = comm.collect_pp(orb_path)
    recommand_ecutwfc = None if not pp_path or not os.path.isfile(os.path.join(pp_path, "ecutwfc.json")) else json.load(open(os.path.join(pp_path, "ecutwfc.json"), "r"))
    jobs = {}
    for ipath in stru_paths:
        istru = os.path.join(ipath, "STRU")
        stru = AbacusStru.ReadStru(istru)
        if stru is None:
            print("Error: read structure failed.")
            continue
        element = stru.get_element(number=False, total=False)
        pp = []
        orb = []
        element_no_pp = [ie for ie in element if ie not in pp_paths]
        if len(element_no_pp) > 0:
            print(f"Error: some elements of {ipath} have no pseudopotentials.\n  {', '.join(element_no_pp)}")
        for ie in element:
            if ie in pp_paths:
                pp.append(os.path.basename(pp_paths[ie]))
                t_file = os.path.join(ipath, os.path.basename(pp_paths[ie]))
                if os.path.isfile(t_file):
                    os.remove(t_file)
                if copy_pp_orb:
                    shutil.copy(os.path.abspath(pp_paths[ie]), t_file)
                else:
                    os.symlink(os.path.abspath(pp_paths[ie]), t_file)
            else:
                pp.append(None)
            if ie in orb_paths:
                orb.append(os.path.basename(orb_paths[ie]))
                t_file = os.path.join(ipath, os.path.basename(orb_paths[ie]))
                if os.path.isfile(t_file):
                    os.remove(t_file)
                if copy_pp_orb:
                    shutil.copy(os.path.abspath(orb_paths[ie]), t_file)
                else:
                    os.symlink(os.path.abspath(orb_paths[ie]), t_file)
            else:
                orb.append(None)
        if None in pp:
            e = " ".join([ie for ie,ipp in zip(element,pp) if ipp is None])
            print(f"WARNING: some elements ({e}) have no pseudopotentials, will NOT set pseudopotential in STRU file.")
        else:
            stru.set_pp(pp)
            
        if None in orb:
            e = " ".join([ie for ie,iorb in zip(element,orb) if iorb is None])
            print(f"WARNING: some elements ({e}) have no orbitals, will NOT set orbital in STRU file.")
        else:
            stru.set_orb(orb)
        stru.write(os.path.join(ipath, "STRU"))
        jobs[ipath] = {
            "element": element,
            "pp": pp,
            "orb": orb,
            "recommand_ecutwfc": [recommand_ecutwfc.get(ie,None) for ie in element] if recommand_ecutwfc else None
        } 
    return jobs

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
                 constrain: List[Union[bool,List[bool]]] = None,
                 lambda_: List[Union[float,List[float]]] = None,
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
            if the magmom of each type, by default None
        magmom_atom: List[Union[float,List[float]]], optional
            if the magmom_atom is not None, then set the magmom of each atom, by default None
            each element should be one float or three float (for non-colinear case), or None (not set)
        velocity: List[List[float]], optional
            the velocity of each atom, by default None
        angle1: List[float], optional
            the angle1 of each atom, by default None
        angle2: float, optional 
            the angle2 of each atom, by default None
        constrain: List[List[bool],bool], optional
            if constrain the magnetization of each atom, by default None
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
            sort_idx = [[] for i in self._label]
            for i,ilabel in enumerate(label):
                coords[self._label.index(ilabel)].append(coord[i])
                sort_idx[self._label.index(ilabel)].append(i)
            self._coord = []
            for i in coords:
                self._coord += i
            sort_idx = [j for i in sort_idx for j in i]  # the index of the original coord after sorting
        else:
            assert(atom_number != None), "ERROR: label number is not equal to coord number, but atom_number is not defined "
            assert(len(label) == len(atom_number)), "ERROR: label number is not equal to atom_number number"
            self._label = label
            self._atom_number = atom_number
            self._coord = coord
            sort_idx = list(range(len(coord)))

        total_atom = np.array(self._atom_number).sum()
        assert(total_atom == len(self._coord)), f"ERROR: the total atom number is not equal to coord number, {total_atom} != {len(coord)}"
        self._cell = cell
        
        # check pp orb paw
        self._pp = self._clean_pporb(pp, label)
        self._orb = self._clean_pporb(orb, label)
        self._paw = self._clean_pporb(paw, label)
        
        if (not self._pp):
            #print("WARNING: pp is not defined!!!")
            self._pp = ["" for i in range(len(self._label))]
            
        
        self._lattice_constant = lattice_constant
        self._dpks = dpks if dpks else None
        self._cartesian = cartesian
        self._move = None if move is None else [move[i] for i in sort_idx]
        self._magmom = magmom if magmom else [0]*len(self._label)
        self._magmom_atom = None if magmom_atom is None else [magmom_atom[i] for i in sort_idx]
        self._velocity = None if velocity is None else [velocity[i] for i in sort_idx]
        self._angle1 = None if angle1 is None else [angle1[i] for i in sort_idx]
        self._angle2 = None if angle2 is None else [angle2[i] for i in sort_idx]
        self._constrain = None if constrain is None else [constrain[i] for i in sort_idx]
        self._lambda = None if lambda_ is None else [lambda_[i] for i in sort_idx]
        
        if element != None:
            self._element = element if len(element) == len(self._label) else [element[label.index(i)] for i in self._label]
        else:
            #print("WARNING: element is not defined, will use the first one/two letter of label as element name")
            self._element = []
            for i in self._label:
                ele = i[0]
                if len(i) > 1 and i[1].islower():
                    ele += i[1]
                self._element.append(ele)
        
        if mass != None:
            self._mass = mass if len(mass) == len(self._label) else [mass[label.index(i)] for i in self._label]
        else:
            self._mass = [constant.MASS_DICT.get(i,1.0) for i in self._element]
        
        self._check()
    
    def _clean_pporb(self,pporb, label_input):
        if pporb:
            if len(pporb) == len(self._coord):
                new_pporb = []
                for i in range(len(self._label)):
                    new_pporb.append(pporb[label_input.index(self._label[i])])
                return new_pporb
            else:
                return pporb
        else:
            return None
    
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
        if self._constrain:
            assert(n_atom == len(self._constrain)), "ERROR: the length of coord is not equal to constrain"
        if self._lambda:
            assert(n_atom == len(self._lambda)), "ERROR: the length of coord is not equal to lambda"
        return True    
    
    def get_natoms(self):
        '''return the total number of atoms'''
        return sum(self._atom_number)
    
    def get_pp(self, total=False):
        if total:
            pp = []
            for idx,i in enumerate(self._atom_number):
                pp += [self._pp[idx]] * i
            return pp
        else:
            return self._pp

    def get_orb(self, total=False):
        if total:
            orb = []
            for idx,i in enumerate(self._atom_number):
                orb += [self._orb[idx]] * i
            return orb
        else:
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
    
    def get_mass(self,total=False):
        if total:
            mass = []
            for idx,i in enumerate(self._atom_number):
                mass += [self._mass[idx]] * i
            return mass
        else:
            return self._mass

    def globalidx2labelidx(self, idx):
        """Get the label and index of the atom in that label from the global index
        """
        if idx < 0:
            raise ValueError(f"The idx {idx} is less than 0!")
        elif idx >= self.get_natoms():
            raise ValueError(f"The idx {idx} is larger than the total atom numbers!")
        pre_atom_num = 0
        for i in range(len(self._atom_number)):
            kind_i_num = self._atom_number[i]
            if idx < pre_atom_num + kind_i_num:
                return self._label[i], idx - pre_atom_num
            pre_atom_num += kind_i_num
    
    def labelidx2globalidx(self, label, idx):
        # Get the global index of the atom from the label and index of that label
        assert (label in self._label), f"{label} is not valid, should be one of {self._label}"
        label_idx = self._label.index(label)
        assert (idx <= self._atom_number[label_idx]), f"{idx} is larger than the atom number of {label}"
        if label_idx == 0:
            return idx
        else:
            return sum(self._atom_number[:label_idx]) + idx
    
    def get_atommag(self,norm = False):
        # return the magmom of each atom
        # if non-colinear, then return a list of three float for that atom.
        # like: [0,0,0,[1,1,1],[1,1,1],0,0,0], which means the atom 1-3,6-8 has mag 0, and atom 4-5 has mag [1,1,1] 
        # norm: if return the norm of the magmom of each atom
        # for non-collinear case, if set the mag and angle1/angle2 at the same time, then firstly calculate the norm of magmom, then trasfer to magx, magy, magz based on angle1 and angle2
        magmom = []
        for i,ilabel in enumerate(self._label):
            magmom += [self._magmom[i]] * self._atom_number[i]
        noncollinear = False
        if self._magmom_atom:
            for ii,i in enumerate(self._magmom_atom):
                if i != None:
                    imag = None
                    if isinstance(i,(list,tuple)) and len(i) == 3: # for non-colinear case
                        imag = list(i)
                        noncollinear = True
                    elif isinstance(i,(list,tuple)) and len(i) == 1 and isinstance(i[0],(int,float)):
                        imag = float(i[0])
                    elif isinstance(i,(int,float)):
                        imag = float(i)
                    
                    # if has set the total mag and angle1 and angle2, then transfer to magx,magy,magz 
                    if imag is not None:
                        angle1 = None if self._angle1 is None else self._angle1[ii]
                        angle2 = None if self._angle2 is None else self._angle2[ii]
                        if angle1 is not None or angle2 is not None:
                            imag = self.angle_to_mag(np.linalg.norm(imag),angle1,angle2)   
                        magmom[ii] = imag
                        
        if noncollinear:
            magmom = [[0,0,i] if not isinstance(i,list) else i for i in magmom]
            
        if norm:
            magmom = [np.linalg.norm(i) for i in magmom]

        return magmom
    
    def set_atommag(self,new_maglist):
        '''set the magmom of each atom
        The new_maglist should be a list of float or list of three float for each atom, or None
        Will set angle1 and angle2 to None
        '''
        if new_maglist is None:
            self._magmom_atom = None
            self.angle1 = None
            self.angle2 = None
            return
        
        if len(new_maglist) != len(self._coord):
            print("ERROR: the length of new_maglist is not equal to coord number")
            sys.exit(1)
        self._magmom_atom = new_maglist
        self.angle1 = None
        self.angle2 = None
    
    def set_atommag_from_mulliken(self, mullikenf, 
                                  step=-1,
                                  check_label=True):
        '''set the magmom of each atom from mulliken population analysis file
        step: the step index to read, default is -1, which means the last step
        check_label: if check the label of mulliken file and stru file
        Will set angle1 and angle2 to None
        '''
        from abacustest.lib_collectdata.comm import get_mulliken
        try:
            atom_mag,atom_elec,atom_label = get_mulliken(mullikenf)
        except Exception as e:
            traceback.print_exc()
            print(f"ERROR: read mulliken file {mullikenf} failed.")
            raise e
        
        if len(atom_mag) == 0:
            print(f"ERROR: no mulliken magmom data in file {mullikenf}.")
            raise RuntimeError(f"no mulliken magmom data in file {mullikenf}.")

        atom_mag = atom_mag[step]
        atom_label = atom_label[step]
        
        if len(atom_mag) != len(self._coord):
            print("ERROR: the atom number in mulliken file is not equal to coord number")
            sys.exit(1)
        
        if check_label and atom_label != self.get_label(total=True):
            print("ERROR: the atom label in mulliken file is not equal to coord label")
            sys.exit(1)
        
        self._magmom_atom = atom_mag
        self.angle1 = None
        self.angle2 = None
    
    def set_angle1(self,new_angle1):
        if new_angle1 is None:
            self._angle1 = None
            return
        if len(new_angle1) != len(self._coord):
            print("ERROR: the length of new_angle1 is not equal to coord number")
            sys.exit(1)
        self._angle1 = new_angle1
        
    def set_angle2(self,new_angle2):
        if new_angle2 is None:
            self._angle2 = None
            return
        if len(new_angle2) != len(self._coord):
            print("ERROR: the length of new_angle2 is not equal to coord number")
            sys.exit(1)
        self._angle2 = new_angle2
    
    def set_constrain(self,new_constrain):
        '''set the constrain of each atom
        The new_constrain should be a list of bool or list of three bool for each atom, or None
        '''
        if len(new_constrain) != len(self._coord):
            print("ERROR: the length of new_constrain is not equal to coord number")
            sys.exit(1)
        self._constrain = new_constrain
    
    def get_constrain(self):
        # return a list of N bool values, each value is True or False or a list of three True or False
        if self._constrain:
            c = []
            for i in self._constrain:
                if i == None:
                    c.append(False)
                else:
                    if isinstance(i,list):
                        c.append([bool(j) for j in i])
                    else:
                        c.append(bool(i))
            return c
        else:
            return [False] * len(self._coord)
    
    def get_isconstrain (self):
        # return a list of N bool values, each value is True or False
        # for noncollinear case, if one component is constrained, then return True
        if self._constrain:
            c = []
            for i in self._constrain:
                if i == None:
                    c.append(False)
                elif isinstance(i,list):
                    if True in [bool(j) for j in i]:
                        c.append(True)
                    else:
                        c.append(False)
                elif isinstance(i,(bool,int)):
                    c.append(bool(i))
                else:
                    c.append(False in i)
            return c
        else:
            return [False] * len(self._coord)
    
    def get_lambda(self):
        # return a list
        if self._lambda:
            l = []
            for i in self._lambda:
                if i == None:
                    l.append(0)
                else:
                    if isinstance(i,list):
                        l.append([float(j) for j in i])
                    else:
                        l.append(float(i))
            return l
        else:
            return [0] * len(self._coord)
    
    def get_angle1(self):
        return self._angle1
    
    def get_angle2(self):
        return self._angle2

    @staticmethod
    def mag_to_angle(magx,magy,magz):
        '''return the angle1 and angle2 of the magnetization'''
        mag = np.array([magx,magy,magz])
        mag /= np.linalg.norm(mag)
        angle1 = np.arccos(mag[2]) * 180 / np.pi
        angle2 = np.arctan2(mag[1],mag[0]) * 180 / np.pi
        return angle1,angle2
    
    @staticmethod
    def angle_to_mag(totmag,angle1,angle2):
        '''return the magnetization of the magnetization'''
        if angle1 is None:
            angle1 = 0
        if angle2 is None:
            angle2 = 0
        angle1 = angle1 * np.pi / 180
        angle2 = angle2 * np.pi / 180
        mag = np.array([np.sin(angle1)*np.cos(angle2),np.sin(angle1)*np.sin(angle2),np.cos(angle1)])
        mag *= totmag
        return mag.tolist()
        
    
    def get_mag(self):
        # return the magmom of each type
        return self._magmom
    
    def get_move(self):
        # return the constraint of each atom
        # 0: not move, 1: move
        if self._move == None:
            return [1,1,1] * len(self._coord)
        else:
            m = []
            for i in self._move:
                if i == None:
                    m.append([1,1,1])
                else:
                    m.append(i)
            return m
    
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
    
    def get_cell_param(self):
        # return the box parameter: a,b,c,alpha,beta,gamma, unit is Angstrom and degree
        cell = np.array(self.get_cell(bohr=False))
        a = np.linalg.norm(cell[0])
        b = np.linalg.norm(cell[1])
        c = np.linalg.norm(cell[2])
        alpha = np.arccos(np.dot(cell[1],cell[2])/(b*c)) * 180 / np.pi
        beta = np.arccos(np.dot(cell[0],cell[2])/(a*c)) * 180 / np.pi
        gamma = np.arccos(np.dot(cell[0],cell[1])/(a*b)) * 180 / np.pi
        return a,b,c,alpha,beta,gamma

    def get_volume(self):
        # return the volume of the cell, in unit of Angstrom^3
        cell = np.array(self.get_cell(bohr=False))
        volume = np.abs(np.linalg.det(cell))
        return volume
    
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
    
    def perturb_stru(self,pert_number,cell_pert_frac=None, atom_pert_dist=None,mag_rotate_angle=None,mag_tilt_angle=None,mag_norm_dist=None):
        '''
        perturb the structure
        cell_pert_frac: the perturb fraction of the cell. Will generate a random perturb matrix.
                        the diagonal part is between -cell_pert_frac and cell_pert_frac, and off-diagonal part is between -cell_pert_frac/2 and cell_pert_frac/2
        atom_pert_dis: the perturb distance of each atom. Will generate a random perturb vector for each atom.
        mag_rotate_angle: the rotation angle of the magnetization of all atom. 
        mag_tilt_angle: the tilt angle of the magnetization of all atom.
        mag_norm_dist: the perturb distance of the norm of the magnetization of all atom.
        '''
        def transfer_range(number):
            if number is None:
                return None
            elif isinstance(number, (float,int)):
                return [0, abs(number)]
            elif isinstance(number, (list,tuple)):
                return [min(number), max(number)]
            else:
                print(f"ERROR: the type of {number} {type(number)} is not supported")
                return None
        
        if pert_number <= 0:
            return [self]
        
        cell_pert_frac = transfer_range(cell_pert_frac)
        atom_pert_dist = transfer_range(atom_pert_dist)
        mag_rotate_angle = transfer_range(mag_rotate_angle)
        mag_tilt_angle = transfer_range(mag_tilt_angle)
        mag_norm_dist = transfer_range(mag_norm_dist)
        
        new_stru = [copy.deepcopy(self) for i in range(pert_number)]
        print(f"perturb {pert_number} new structures")
        print(f"cell_pert_frac: {cell_pert_frac}")
        print(f"atom_pert_dist: {atom_pert_dist}")
        print(f"mag_rotate_angle: {mag_rotate_angle}")
        print(f"mag_tilt_angle: {mag_tilt_angle}")
        print(f"mag_norm_dist: {mag_norm_dist}")
        
        if cell_pert_frac is not None:
            for i in range(pert_number):
                icell = new_stru[i].get_cell(bohr=False)
                new_cell,_ = comm.perturb_cell(icell,cell_pert_frac)
                new_stru[i].set_cell(new_cell,bohr=False,change_coord=True)
        
        if atom_pert_dist is not None:
            for i in range(pert_number):
                new_coord = comm.perturb_coord(new_stru[i].get_coord(bohr=False,direct=False),atom_pert_dist)
                new_stru[i].set_coord(new_coord,direct=False,bohr=False)
        
        # for perturbation of magnetization, only performed on the atom who's magmom is constrained
        atom_mag = self.get_atommag()
        is_sc = self.get_isconstrain()
        if True not in is_sc:
            return new_stru
        
        noncollinear = False
        for i in atom_mag:
            if isinstance(i,list) and len(i) == 3:
                noncollinear = True
                break
    
        if noncollinear and mag_rotate_angle is not None:
            for i in range(pert_number):
                # perturb all magmom with a same random angle
                atom_mag = new_stru[i].get_atommag()
                new_mag = comm.pert_vector(atom_mag,mag_rotate_angle)
                new_mag = [new_mag[j] if is_sc[j] else atom_mag[j] for j in range(len(new_mag))]
                new_stru[i].set_atommag(new_mag)
                new_stru[i].set_angle1(None)
                new_stru[i].set_angle2(None)

        if noncollinear and mag_tilt_angle is not None:
            for i in range(pert_number):
                atom_mag = new_stru[i].get_atommag()
                new_mag = []
                for j in range(len(atom_mag)):
                    if is_sc[j]:
                        new_mag.append(comm.pert_vector([atom_mag[j]],mag_tilt_angle)[0])
                    else:
                        new_mag.append(atom_mag[j])
                new_stru[i].set_atommag(new_mag)
                new_stru[i].set_angle1(None)
                new_stru[i].set_angle2(None)
        
        if mag_norm_dist is not None:
            for i in range(pert_number):
                atom_mag = new_stru[i].get_atommag() 
                new_mag = []
                for j in range(len(atom_mag)):
                    if is_sc[j]:
                        mag_norm_random = np.random.uniform(mag_norm_dist[0],mag_norm_dist[1]) * np.random.choice([-1,1])
                        mag_norm = np.linalg.norm(atom_mag[j]) + mag_norm_random
                        if mag_norm != mag_norm_random:
                            new_magj = np.array(atom_mag[j]) / np.linalg.norm(atom_mag[j]) * mag_norm
                        else: # if the original mag is 0, then generate a random mag
                            tmp_mag = np.random.rand(3) - 0.5
                            new_magj = tmp_mag / np.linalg.norm(tmp_mag) * mag_norm
                        new_mag.append(new_magj.tolist())
                    else:
                        new_mag.append(atom_mag[j])
                        
                new_stru[i].set_atommag(new_mag)
                new_stru[i].set_angle1(None)
                new_stru[i].set_angle2(None)
        return new_stru
         
        
    
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
            print("pplist:",pplist)
            print("label:",self._label)
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
    
    def set_coord(self,coord,direct=False,bohr=True, keep_lattice_constant=True):
        '''
        set the coordinate of each atom
        
        if direct is True, then the coord is direct type, and will modify self._cartesian to False
        else the coord is cartesian type, and will modify self._cartesian to True, and will set lattice_constant to 1.0, and modify cell *= lattice_constant
        
        if bohr is False, then will transfer the coord to Bohr unit
        
        if keep_lattice_constant is True, then will keep the lattice constant, else will set the lattice constant to 1.0
        '''
        unit_coef = 1 if bohr else constant.A2BOHR # need save the unit with bohr
        if direct:
            self._cartesian = False
            self._coord = coord
        else:
            self._cartesian = True
            if keep_lattice_constant:
                new_coord = np.array(coord) / self._lattice_constant * unit_coef
                self._coord = new_coord.tolist()
            else:
                self._cell = (np.array(self._cell) * self._lattice_constant).tolist()
                self._coord = (np.array(coord) * unit_coef).tolist()
                self._lattice_constant = 1.0  
        self._check()
    
    def set_cell(self,cell,bohr=True,change_coord=True,keep_lattice_constant=True):
        '''
        set the lattice of a, b, c
        
        if bohr is False, then will transfer the cell to Bohr unit
        
        if change_coord is True, then will transfer the coord based on the new cell and keep the relative position of each atom
        if change_coord is False, then will not modify the coord, which means the coord is not changed whatever is direct or cartesian type
        '''
        unit_coef = 1 if bohr else constant.A2BOHR # need save the unit with bohr
        
        cell = np.array(cell) * unit_coef # now cell is in Bohr unit
        if keep_lattice_constant:
            if change_coord and self._cartesian: 
                coord = self.get_coord(bohr=bohr,direct=True)
                new_coord = np.array(coord).dot(cell) / self._lattice_constant
                self._coord = new_coord.tolist()

            self._cell = (cell / self._lattice_constant ).tolist()
        else:
            # will change the lattice constant to 1.0
            if self._cartesian:
                if change_coord:
                    coord = self.get_coord(bohr=bohr,direct=True)
                    new_coord = np.array(coord).dot(cell)
                    self._coord = new_coord.tolist()
                else:
                    self._coord = (np.array(self._coord) * self._lattice_constant).tolist()
                    
            self._cell = cell.tolist()
            self._lattice_constant = 1.0

        self._check()        
    
    def split_list(self, alist, indices):
        '''
        split a list to several sublists based on the indices
        '''
        if not indices:
            return [alist]
        if len(alist) != sum(indices):
            print("ERROR: the sum of indices is not equal to the length of alist")
            print("indices:",indices)
            print("alist",alist)
            sys.exit(1)
        new_list = []
        start_idx = 0
        for i in indices:
            new_list.append(alist[start_idx:start_idx+i])
            start_idx += i
        return new_list
    
    def supercell(self,nabc = [1,1,1]):
        '''
        generate a supercell
        
        na,nb,nc: the number of the supercell in a,b,c direction
        
        return a new AbacusStru object of the supercell.
        The atom label order is same as the original structure, and the atom order of each label is repeatedly added based on the original order.
        '''
        # be careful, the real coord should be coord * lattice_constant
        na, nb, nc = nabc
        new_stru = copy.deepcopy(self)
        cell = np.array(new_stru._cell)
        coord = np.array(new_stru._coord)
        new_stru._cell = cell * np.array([[na],[nb],[nc]])
        
        if not new_stru._cartesian:
            # if the coord is direct type, then transfer the original coord based on new cell
            coord = coord / np.array([na,nb,nc])

        # split the coord to coords by each atom type
        coords = self.split_list(coord.tolist(),new_stru._atom_number)
        key_p = ["_move","_magmom_atom","_velocity","_angle1","_angle2","_constrain","_lambda"]
        kv = {i:None if getattr(new_stru,i) is None else self.split_list(getattr(new_stru,i), new_stru._atom_number) for i in key_p}

        for ia in range(na):
            for ib in range(nb):
                for ic in range(nc):
                    if ia == 0 and ib == 0 and ic == 0:
                        continue
                    
                    start_idx = 0
                    for idx, iatomn in enumerate(new_stru._atom_number):
                        if new_stru._cartesian:
                            add_coord = coord[start_idx:start_idx+iatomn] + np.array([ia,ib,ic]).dot(cell) 
                        else:
                            add_coord = coord[start_idx:start_idx+iatomn] + np.array([ia/na,ib/nb,ic/nc])
                        add_coord = add_coord.tolist()
                        coords[idx] += add_coord

                        for attri in key_p:
                            if kv[attri] is not None:
                                kv[attri][idx] += getattr(self,attri)[start_idx:start_idx+iatomn]                               
                        start_idx += iatomn

        new_stru._coord = [j for i in coords for j in i]
        for attri in key_p:
            if kv[attri] is not None:
                setattr(new_stru,attri,[j for i in kv[attri] for j in i])
        new_stru._atom_number = [len(i) for i in coords]
        new_stru._check()
        return new_stru 
    
    def delete_atom(self, 
                    label:str, 
                    idx: Union[int, List[int], None] = None):
        """Delete the atom of the element in del_element.
        
        label: str, the atomic label to be deleted, e.g. "H"
        idx: int or list of int, the index of the atom to be deleted, if None, then delete all the atom of the element
        """
        if label not in self._label:
            print(f"ERROR: the label {label} is not in the structure")
            print("label:",self._label)
            raise ValueError(f"the label {label} is not in the structure")
        if isinstance(idx,int):
            idx = [idx]
        elif idx is None:
            idx = list(range(self._atom_number[self._label.index(label)]))
        elif isinstance(idx,list):
            pass
        else:
            print(f"ERROR: the type of idx {type(idx)} is not supported")
            raise ValueError(f"the type of idx {type(idx)} is not supported")
        
        if max(idx) >= self._atom_number[self._label.index(label)] or min(idx) < 0:
            print(f"ERROR: the index in idx is out of range")
            print("idx:",idx)
            print("atom number of",label,":",self._atom_number[self._label.index(label)])
            raise ValueError(f"the index in idx is out of range")

        new_stru = copy.deepcopy(self)
        global_idx = [sum(new_stru._atom_number[:new_stru._label.index(label)]) + i for i in idx]
        new_stru._coord = [i for j, i in enumerate(new_stru._coord) if j not in global_idx]
        new_stru._atom_number[new_stru._label.index(label)] -= len(idx)
        if new_stru._atom_number[new_stru._label.index(label)] == 0:
            del_idx = new_stru._label.index(label)
            new_stru._label.pop(del_idx)
            new_stru._mass.pop(del_idx)
            new_stru._magmom.pop(del_idx)
            if new_stru._pp:
                new_stru._pp.pop(del_idx)
            if new_stru._orb:
                new_stru._orb.pop(del_idx)
            if new_stru._paw:
                new_stru._paw.pop(del_idx)
            new_stru._atom_number.pop(del_idx)
        key_p = ["_move","_magmom_atom","_velocity","_angle1","_angle2","_constrain","_lambda"]
        for attri in key_p:
            if getattr(new_stru,attri) is not None:
                setattr(new_stru,attri,[i for j, i in enumerate(getattr(new_stru,attri)) if j not in global_idx])
        new_stru._check()
        return new_stru
    
    def set_empty_atom(self, label:str, idx: Union[int, List[int], None] = None):
        """Set the specified atom to be empty by adding 'empty' to the label. Such as: H -> H_empty.
        
        This is useful for BSSE calculation, where some atoms are set to be empty.
        
        Args:
            label: str, the atomic label to be set to empty, e.g. "H"
            idx: int or list of int, the index of the atom to be set to empty, if None, then set all the atom of the element to empty
        """    
        if label not in self._label:
            print(f"ERROR: the label {label} is not in the structure")
            print("label:",self._label)
            raise ValueError(f"the label {label} is not in the structure")
        
        label_idx = self._label.index(label)
            
        if idx is None:
            new_label = label + "_empty"
            self._label[label_idx] = new_label
            return
        
        if isinstance(idx,int):
            idx = [idx]
        
        if len(idx) == self._atom_number[label_idx]:
            new_label = label + "_empty"
            self._label[label_idx] = new_label
            return
        
        atom_idx = [sum(self._atom_number[:label_idx]) + i for i in idx]
        new_label = label + "_empty"
        
        self._label.append(new_label)
        self._element.append(self._element[label_idx])
        self._mass.append(self._mass[label_idx])
        self._magmom.append(self._magmom[label_idx])
        if self._pp:
            self._pp.append(self._pp[label_idx])
        if self._orb:
            self._orb.append(self._orb[label_idx])
        if self._paw:
            self._paw.append(self._paw[label_idx])
        self._atom_number.append(len(idx))
        self._atom_number[label_idx] -= len(idx)

        key_p = ["_coord","_move","_magmom_atom","_velocity","_angle1","_angle2","_constrain","_lambda"]
        for attri in key_p:
            if getattr(self,attri) is not None:
                old_value = getattr(self,attri)
                new_value = []
                for j in range(len(old_value)):
                    if j not in atom_idx:
                        new_value.append(old_value[j])
                for j in atom_idx:
                    new_value.append(old_value[j])
                setattr(self,attri,new_value)
        self._check()
                
            
        
     
    def write(self,struf="STRU", direct=None):
        '''
        write to STRU file
        direct: if True, then write the coord in direct type, if False, then write in cartesian type
                if None, then write in the original type
        '''
        if direct is None:
            if self._cartesian:
                direct = False
            else:
                direct = True
        
        if direct:
            coord = self.get_coord(direct=True)
        else:
            coord = np.array(self.get_coord(direct=False, bohr = True)) / self._lattice_constant

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

        if not direct:
            cc += "Cartesian\n"
        else:
            cc += "Direct\n"
        icoord = 0
        for i,ilabel in enumerate(self._label):
            cc += "\n%s\n%f\n%d\n" % (ilabel,self._magmom[i],self._atom_number[i])
            for j in range(self._atom_number[i]):
                cc += "%17.11f %17.11f %17.11f " % tuple(coord[icoord + j])
                if self._move and self._move[icoord + j] and len(self._move[icoord + j]) == 3:
                    cc += "%d %d %d " % tuple(self._move[icoord + j])
                if self._magmom_atom:
                    if isinstance(self._magmom_atom[icoord + j],list):
                        if len(self._magmom_atom[icoord + j]) == 3:
                            cc += "mag %12.8f %12.8f %12.8f " % tuple(self._magmom_atom[icoord + j])
                        elif len(self._magmom_atom[icoord + j]) == 1:
                            cc += "mag %12.8f " % self._magmom_atom[icoord + j][0]
                    elif self._magmom_atom[icoord + j] != None:
                        cc += "mag %12.8f " % self._magmom_atom[icoord + j]
                if self._velocity:
                    if self._velocity[icoord + j]:
                        cc += "v %f %f %f " % tuple(self._velocity[icoord + j])
                if self._angle1 and self._angle1[icoord + j] != None:
                        cc += "angle1 %f " % self._angle1[icoord + j]
                if self._angle2 and self._angle2[icoord + j] != None:
                        cc += "angle2 %f " % self._angle2[icoord + j]
                if self._constrain and self._constrain[icoord + j]:
                    if isinstance(self._constrain[icoord + j],list) and len(self._constrain[icoord + j]) == 3:
                        cc += "sc " + " ".join(["1" if ic else "0" for ic in self._constrain[icoord + j]]) + " " 
                    elif isinstance(self._constrain[icoord + j],list) and len(self._constrain[icoord + j]) == 1:
                        cc += "sc " + "1" if self._constrain[icoord + j][0] else "0" + " "
                    elif not isinstance(self._constrain[icoord + j],list):
                        cc += "sc " + "1" if self._constrain[icoord + j] else "0" + " "
                    else:
                        print("ERROR: the constrain is not a list or a bool value, skip it")
                        print("\t\tconstrain:",self._constrain[icoord + j])
                if self._lambda and self._lambda[icoord + j]:
                        cc += "lambda " + " ".join([str(ic) for ic in self._lambda[icoord + j]]) + " "
                cc += "\n"
            icoord += self._atom_number[i]
        
        #write dpks
        if self._dpks:
            cc += "\nNUMERICAL_DESCRIPTOR\n"
            cc += self._dpks    

        Path(struf).parent.mkdir(parents=True, exist_ok=True)
        Path(struf).write_text(cc) 
    
    def write2poscar(self,poscar="POSCAR",empty2x=False):
        '''
        write to POSCAR file for VASP
        Arguments:
            poscar: the POSCAR file name
            empty2x: if True, then transfer the label with 'empty' to 'X' in POSCAR file
        '''
        cc = "STRUCTURE translated by abacustest\n"
        cc += f"{self._lattice_constant * constant.BOHR2A}\n"
        for i in self._cell:
            cc += "%17.11f %17.11f %17.11f\n" % tuple(i)
        
        if empty2x:
            labels = [i if "empty" not in i else "X" for i in self.get_label(total=True)]
            new_label = []
            atom_number = []
            for i in labels:
                if i not in new_label:
                    new_label.append(i)
            for i in range(len(new_label)):
                atom_number.append(labels.count(new_label[i]))
            cc += " ".join(new_label)+"\n"
            cc += " ".join([str(i) for i in atom_number])+"\n"    
        else:
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
    
    def write2cif(self,cif="STRU.cif", empty2x=False):
        '''
        write to CIF file
        '''
        from ase.io import write

        stru_cif = self.to_ase(empty2x=empty2x)
        write(cif, stru_cif, format='cif')
    
    def to_ase(self, empty2x=False):
        '''
        convert the AbacusStru to ase.Atoms
        
        Returns:
            atoms: the ase.Atoms object
        '''
        from ase import Atoms
        
        cell = self.get_cell(bohr=False)
        coord = self.get_coord(bohr=False,direct=False)
        elements = self.get_element(total=True, number=False)
        if empty2x:
            labels = self.get_label(total=True)
            for i in range(len(labels)):
                if 'empty' in labels[i]:
                    elements[i] = 'X'

        magmom = self.get_atommag()
        return Atoms(symbols=elements, positions=coord, cell=cell, pbc=True, magmoms=magmom,
                      info={"pp": self.get_pp(),
                            "orb": self.get_orb(),
                            "paw": self.get_paw(),
                            "dpks": self.get_dpks()
                      })
    
    def to_pymatgen(self):
        '''
        convert the AbacusStru to pymatgen.Structure
        
        Returns:
            structure: the pymatgen.Structure object
        '''
        from pymatgen.core import Structure
        
        cell = self.get_cell(bohr=False)
        coord = self.get_coord(bohr=False,direct=False)
        elements = self.get_element(number=False,total=True)
        labels = self.get_label(total=True)
        return Structure(lattice=cell, species=elements, coords=coord, coords_are_cartesian=True,labels=labels,)

    def to_phonopy(self):
        '''
        convert the AbacusStru to phonopy.Phonopy
        
        Returns:
            phonopy: the phonopy.Phonopy object
        '''
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms
        
        cell = self.get_cell(bohr=False)
        coord = self.get_coord(bohr=False,direct=False)
        elements = self.get_element(number=False,total=True)
        mass = self.get_mass(total=True)
        atom_mag = self.get_atommag()
        all_zero = True
        for i in atom_mag:
            if i != 0 and i != [0,0,0]:
                all_zero = False
                break
        if all_zero:
            atom_mag = None
        
        return PhonopyAtoms(symbols=elements, positions=coord, cell=cell, masses=mass,
                            magnetic_moments=atom_mag)

    
    @staticmethod
    def stru2ase(stru):
        '''
        convert the AbacusStru to ase.Atoms
        
        Arguments:
            stru: the STRU file
        
        Returns:
            atoms: the ase.Atoms object
        '''
        
        stru = AbacusStru.ReadStru(stru)
        if stru is not None:
            return stru.to_ase()
        else:
            raise ValueError("The STRU file is not valid or cannot be read.")

    @staticmethod
    def stru2pymatgen(stru):
        '''
        convert the AbacusStru to pymatgen.Structure
        
        Arguments:
            stru: the STRU file
            
        Returns:
            structure: the pymatgen.Structure object
        '''
        
        stru = AbacusStru.ReadStru(stru)
        if stru is not None:
            return stru.to_pymatgen()
        else:
            raise ValueError("The STRU file is not valid or cannot be read.")
    
    @staticmethod
    def stru2phonopy(stru):
        '''
        convert the AbacusStru to phonopy.Phonopy
        
        Arguments:
            stru: the STRU file
            
        Returns:
            phonopy: the phonopy.Phonopy object
        '''
        
        stru = AbacusStru.ReadStru(stru)
        if stru is not None:
            return stru.to_phonopy()
        else:
            raise ValueError("The STRU file is not valid or cannot be read.")

    @staticmethod
    def ReadStru(stru:str = "STRU"):
        '''Read the structure from STRU file, and return an AbacusStru object.
        '''
        stru_values = read_stru_file(stru)

        return AbacusStru(**stru_values)
    
    @staticmethod
    def FromDpdata(input_stru: Union[str, Path],
                   input_fmt: str):
        '''Read the structure from dpdata supportted file, such as POSCAR, etc.
        
        Args:
            input_stru (str or Path): the file name of the structure file
            input_fmt (str): the format of the structure file, such as "poscar" etc.
        '''
        import dpdata
        import uuid
        
        stru = dpdata.System(input_stru, fmt=input_fmt)
        tmp_file_name = "tmp." + str(uuid.uuid4())[:8]
        stru.to("abacus/stru", tmp_file_name)
        
        new_stru = AbacusStru.ReadStru(tmp_file_name)
        if os.path.isfile(tmp_file_name):
            os.remove(tmp_file_name)
        if new_stru is None:
            raise ValueError(f"Can not read the structure from {input_stru} with format {input_fmt}")
        return new_stru


def ReadKpt(kptpath):
    '''
    kptpath should be a file name of KPT file or a path of ABACUS inputs.
    
    return kpt,model
    - kpt is a list of k-point + shift, such as [1,1,1,0,0,0] 
    - model is the model of k-point, such as "mp","gamma","direct","line"
    '''
    if os.path.isdir(kptpath):
        # try to file input
        if os.path.isfile(os.path.join(kptpath,"INPUT")):
            input_param = ReadInput(os.path.join(kptpath,"INPUT"))
            kptf = os.path.join(kptpath,input_param.get("kpoint_file","KPT"))
            struf = os.path.join(kptpath,input_param.get("stru_file","STRU"))
            kspacing = input_param.get("kspacing",None)
            if input_param.get("basis_type","").lower() == "lcao" and comm.IsTrue(input_param.get("gamma_only",False)):
                print("Have set gamma_only in INPUT file, will use 1 1 1 for KPOINT.")
                return [1,1,1,0,0,0],"gamma"
            elif kspacing != None and kspacing != 0:
                
                if not os.path.isfile(struf):
                    print("  Can not find the STRU file, and try to read KPOINT from KPT file")
                    if os.path.isfile(kptf):
                        return ReadKpt(kptf)
                    else:
                        print("  Can not find the KPT file:",kptf)
                        sys.exit(1)
                else:
                    stru = AbacusStru.ReadStru(struf)
                    if stru == None:
                        print("  Can not read the STRU file:",struf)
                        sys.exit(1)
                    cell = stru.get_cell(bohr=True)
                    kpt = comm.kspacing2kpt(kspacing,cell) + [0,0,0]
                    print(f"Transfer kspacing: {kspacing} to K points {kpt[:3]}.")
                    return kpt,"gamma"
            elif os.path.isfile(kptf):
                return ReadKpt(kptf)
            else:
                print("ERROR: Can not find the KPT file:",kptf)
                sys.exit(1)
        elif os.path.isfile(os.path.join(kptpath,"KPT")):
            return ReadKpt(os.path.join(kptpath,"KPT"))
        else:
            print("ERROR: Can not find the INPUT/KPT file in the path:",kptpath)
            sys.exit(1)
    elif os.path.isfile(kptpath):
        with open(kptpath) as f1: lines = [i for i in f1.readlines() if i.split("#")[0].strip() != ""]
        model = lines[2].split()[0].lower()
        if model.startswith("m"):
            model = "mp"
        elif model.startswith("g"):
            model = "gamma"
        elif model.startswith("d"):
            model = "direct"
        elif model.startswith("c"):
            model = "cartesian"
        elif model.lower() == "line":
            model = "line"
        elif model.lower() == "line_cartesian":
            model = "line_cartesian"
        else:
            print(f"ERROR: the model of KPT file is not support now!!!\n{model}")
            sys.exit(1)
            
        if model in ["mp","gamma"]:
            kpt = [int(i) for i in lines[3].split()[:3]] + [float(i) for i in lines[3].split()[3:6]]
            return kpt,model
        elif model in ["direct", "cartesian"]:
            nk = int(lines[1].split()[0])
            kpt = []
            for i in range(nk):
                kpt.append([float(ii) for ii in lines[2+i].split()[:4]])
            return kpt,model
        elif model in ["line", "line_cartesian"]:
            kpt = []
            for line in lines[3:]:
                if line.strip() == "":
                    break
                sline = line.split()
                ik = [float(i) for i in sline[:3]] + [int(sline[3])]
                if len(sline) == 4:
                    kpt.append(ik)
                elif len(sline) > 4:
                    kpt.append(ik + [" ".join(sline[4:])])
            return kpt,model
    else:
        print(f"ERROR: {kptpath} is not a file or path!!!")
        sys.exit(1)

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
        model_ = "Gamma" if model.lower() == "gamma" else "MP"
        with open(file_name,'w') as f1:
            f1.write(F"K_POINTS\n0\n{model_}\n")
            if len(kpoint_list) == 3:
                kpoint_list += [0,0,0]
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
    elif model.lower() in ["line", "line_cartesian"]:
        with open(file_name,'w') as f1:
            f1.write(f"K_POINTS\n{len(kpoint_list)}\n")
            if model.lower() == "line":
                f1.write("Line\n")
            else:
                f1.write("Line_Cartesian\n")
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
    """Read the INPUT file and return a dictionary of input parameters.
    
    Args:
        INPUTf (str): the file name of the INPUT file.
        input_lines (str): the lines of the INPUT file.
    
    Returns:
        Dict[str, any]: a dictionary of input parameters.
        
    If both `INPUTf` and `input_lines` are provided, `input_lines` will be ignored.
    The value of each parameter will be converted to int or float as appropriate.
    """
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
        iis = ii.split()
        if len(iis) > 1:
            try:
                return [int(i) for i in iis]
            except:
                pass
            try:
                return [float(i) for i in iis]
            except:
                return iis
   
        try:
            return int(ii)
        except:
            pass
        try:
            return float(ii)
        except:
            return ii
    
    for i,iline in enumerate(input_lines):
        if iline.strip() == '' or iline.strip()[0] in ['#']:
            continue
        else:
            sline = re.split('[ \t]',iline.split("#")[0].strip(),maxsplit=1)
            if len(sline) == 2:
                k = sline[0].lower().strip()
                v = str2intfloat(sline[1].strip())
                input_context[k] = v
    return input_context

def WriteInput(input_context:Dict[str,any],
               INPUTf: str = "INPUT"):
    """Write the input parameters to the INPUT file.
    Args:
        input_context (Dict[str, any]): a dictionary of input parameters.
        INPUTf (str): the file name of the INPUT file to write.
    Returns:
        None: the function will write the input parameters to the INPUT file.
        
    If the value of a parameter is None, it will be written as a comment line.
    If the value is a list or tuple, it will be joined with spaces.

    """
    out = "INPUT_PARAMETERS\n"
    for k,v in input_context.items():
        if v != None:
            if isinstance(v, (list,tuple)):
                v = ' '.join([str(i) for i in v])
            out += f"{k}     {v}\n"    
        else:
            out += f"#{k}  \n"
    with open(INPUTf,'w') as f1: f1.write(out)        