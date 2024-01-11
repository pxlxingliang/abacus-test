from typing import List, Dict
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
                 magmom: List[float] = None,
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
        magmom : float, optional
            magmom setting of each type, by default None
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
        
        if (not self._pp) and (not self._paw):
            print("ERROR: Please define the pseudopotential or paw file")
            sys.exit(1)
        if self._pp:
            assert(len(self._label) == len(self._pp)), "ERROR: label number is not equal to pp number"
        if self._paw:
            assert(len(self._label) == len(self._paw)), "ERROR: label number is not equal to paw number"
        if self._orb:
            assert(len(self._label) == len(self._orb)), "ERROR: label number is not equal to orb number"
            
        
        self._lattice_constant = lattice_constant
        self._dpks = dpks if dpks else None
        self._cartesian = cartesian
        self._magmom = magmom if magmom else [0]*len(label)
        
        if element != None:
            self._element = element
        else:
            print("WARNING: element is not defined, will use the first one/two letter of label as element name")
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

    def get_mag(self):
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
                  angle_tolerance=-1.0,):
        '''
        Use seekpath to generate the k-line: https://github.com/giovannipizzi/seekpath
        orig_cell: if True, then use the input cell to generate the k-line, else use primitive cell
        return a dict:
        {
            "cell": cell, # unit: Bohr
            "coord": coord, # direct type
            "atom_label": label, # a list of the label of each atom
            "point_coords": {"Gamma": [0,0,0], "X": [0.5,0,0], ...}, # a dict of the coordinate of each high symmetry point
            "path": [("Gamma","X"),("X","M"),...], # a list of the k-line, each element is a tuple of two high symmetry point
        }
        
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
        
        if orig_cell:
            kpath = seekpath.get_path_orig_cell(stru,with_time_reversal=with_time_reversal,recipe=recipe,threshold=threshold,symprec=symprec,angle_tolerance=angle_tolerance)
            new_cell = cell
            new_coord = coord
            new_label = labels
        else:
            kpath = seekpath.get_path(stru,with_time_reversal=with_time_reversal,recipe=recipe,threshold=threshold,symprec=symprec,angle_tolerance=angle_tolerance)
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
                cc += "%17.11f %17.11f %17.11f 1 1 1\n" % tuple(self._coord[icoord + j])
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
    def ReadStru(stru:str = "STRU"):
        "read the label, pp, orb, cell, coord, deepks-descriptor"
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
            print("WARNING: LATTICE_VECTORS is not correct !!!!!!")

        #read coordinate and coordinate type and atom number of each type
        atom_number = []
        coords = []
        magmom = []
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
            magmom.append(float(atom_positions[i+1].split()[0]))
            atom_number.append(int(atom_positions[i+2].split()[0]))
            i += 3
            for j in range(atom_number[-1]):
                coords.append([float(k) for k in atom_positions[i+j].split()[:3]])
            i += atom_number[-1]
        
        return AbacusStru(label=labels,
                          atom_number=atom_number,
                          cell=cell,
                          coord=coords,
                          pp=pp,
                          orb=orb,
                          paw = paw,
                          lattice_constant=lattice_constant,
                          magmom=magmom,
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