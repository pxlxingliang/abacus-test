import os,sys,argparse,glob,json,traceback,re
import numpy as np
from typing import Dict, List, Optional, Union
import copy
from pathlib import Path
import dpdata

MASS_DICT = {
    "H": 1.0079,
    "He": 4.0026,
    "Li": 6.941,
    "Be": 9.0122,
    "B": 10.811,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.9984,
    "Ne": 20.1797,
    "Na": 22.9897,
    "Mg": 24.305,
    "Al": 26.9815,
    "Si": 28.0855,
    "P": 30.9738,
    "S": 32.065,
    "Cl": 35.453,
    "K": 39.0983,
    "Ar": 39.948,
    "Ca": 40.078,
    "Sc": 44.9559,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938,
    "Fe": 55.845,
    "Ni": 58.6934,
    "Co": 58.9332,
    "Cu": 63.546,
    "Zn": 65.39,
    "Ga": 69.723,
    "Ge": 72.64,
    "As": 74.9216,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.8,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.9059,
    "Zr": 91.224,
    "Nb": 92.9064,
    "Mo": 95.94,
    "Tc": 98,
    "Ru": 101.07,
    "Rh": 102.9055,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.411,
    "In": 114.818,
    "Sn": 118.71,
    "Sb": 121.76,
    "I": 126.9045,
    "Te": 127.6,
    "Xe": 131.293,
    "Cs": 132.9055,
    "Ba": 137.327,
    "La": 138.9055,
    "Ce": 140.116,
    "Pr": 140.9077,
    "Nd": 144.24,
    "Pm": 145,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.9253,
    "Dy": 162.5,
    "Ho": 164.9303,
    "Er": 167.259,
    "Tm": 168.9342,
    "Yb": 173.04,
    "Lu": 174.967,
    "Hf": 178.49,
    "Ta": 180.9479,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.078,
    "Au": 196.9665,
    "Hg": 200.59,
    "Tl": 204.3833,
    "Pb": 207.2,
    "Bi": 208.9804,
    "Po": 209,
    "At": 210,
    "Rn": 222,
    "Fr": 223,
    "Ra": 226,
    "Ac": 227,
    "Pa": 231.0359,
    "Th": 232.0381,
    "Np": 237,
    "U": 238.0289,
    "Am": 243,
    "Pu": 244,
    "Cm": 247,
    "Bk": 247,
    "Cf": 251,
    "Es": 252,
    "Fm": 257,
    "Md": 258,
    "No": 259,
    "Rf": 261,
    "Lr": 262,
    "Db": 262,
    "Bh": 264,
    "Sg": 266,
    "Mt": 268,
    "Rg": 272,
    "Hs": 277,
    "H": 1.0079,
    "He": 4.0026,
    "Li": 6.941,
    "Be": 9.0122,
    "B": 10.811,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.9984,
    "Ne": 20.1797,
    "Na": 22.9897,
    "Mg": 24.305,
    "Al": 26.9815,
    "Si": 28.0855,
    "P": 30.9738,
    "S": 32.065,
    "Cl": 35.453,
    "K": 39.0983,
    "Ar": 39.948,
    "Ca": 40.078,
    "Sc": 44.9559,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938,
    "Fe": 55.845,
    "Ni": 58.6934,
    "Co": 58.9332,
    "Cu": 63.546,
    "Zn": 65.39,
    "Ga": 69.723,
    "Ge": 72.64,
    "As": 74.9216,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.8,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.9059,
    "Zr": 91.224,
    "Nb": 92.9064,
    "Mo": 95.94,
    "Tc": 98,
    "Ru": 101.07,
    "Rh": 102.9055,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.411,
    "In": 114.818,
    "Sn": 118.71,
    "Sb": 121.76,
    "I": 126.9045,
    "Te": 127.6,
    "Xe": 131.293,
    "Cs": 132.9055,
    "Ba": 137.327,
    "La": 138.9055,
    "Ce": 140.116,
    "Pr": 140.9077,
    "Nd": 144.24,
    "Pm": 145,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.9253,
    "Dy": 162.5,
    "Ho": 164.9303,
    "Er": 167.259,
    "Tm": 168.9342,
    "Yb": 173.04,
    "Lu": 174.967,
    "Hf": 178.49,
    "Ta": 180.9479,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.078,
    "Au": 196.9665,
    "Hg": 200.59,
    "Tl": 204.3833,
    "Pb": 207.2,
    "Bi": 208.9804,
    "Po": 209,
    "At": 210,
    "Rn": 222,
    "Fr": 223,
    "Ra": 226,
    "Ac": 227,
    "Pa": 231.0359,
    "Th": 232.0381,
    "Np": 237,
    "U": 238.0289,
    "Am": 243,
    "Pu": 244,
    "Cm": 247,
    "Bk": 247,
    "Cf": 251,
    "Es": 252,
    "Fm": 257,
    "Md": 258,
    "No": 259,
    "Rf": 261,
    "Lr": 262,
    "Db": 262,
    "Bh": 264,
    "Sg": 266,
    "Mt": 268,
    "Rg": 272,
    "Hs": 277,
}

ABACUS_STRU_KEY_WORD = [
    "ATOMIC_SPECIES",
    "NUMERICAL_ORBITAL",
    "LATTICE_CONSTANT",
    "LATTICE_VECTORS",
    "ATOMIC_POSITIONS",
    "NUMERICAL_DESCRIPTOR",
]

BOHR2ANGSTROM = 0.5291875321107901

class AbacusStru:
    def __init__(self,
                 label:List[str],
                 atom_number:List[int],
                 cell:List[List[float]] ,
                 coord:List[List[float]] ,
                 pp:List[str],
                 orb:List[str] = None,
                 lattice_constance:float = 1,
                 magmom: List[float] = None,
                 dpks:str = None,
                 cartessian: bool = False):    
        """ABACUS STRU class, the unit is Bohr

        Parameters
        ----------
        label : List[str]
            the label name of each type
        atom_number : List[int]
            the atom number of each type
        pp : List[str]
            pseudopotential file of each type
        orb : List[str], optional
            orbital file of each type, by default None
        cell : List[List[float]]
            cell of a, b, c
        coord : List[List[float]]
            coordinate of each atom
        lattice_constance : float, optional
            the lattice constance, by default 1
        magmom : float, optional
            magmom setting of each type, by default None
        dpks : str, optional
            the deepks descriptor file name, by default None
        cartessian : bool, optional
            if the coordinate is cartessian type, default is direct type, by default False
        """        
        #check if label number is equal to pp number and orb number
        assert(len(label) == len(pp))
        if orb != None:
            assert(len(label) == len(orb))
        total_atom = np.array(atom_number).sum()
        assert(total_atom == len(coord))
        
        self._label = label
        self._atom_number = atom_number
        self._cell = cell
        self._coord = coord
        self._pp = pp
        self._orb = orb
        self._lattice_constance = lattice_constance
        self._dpks = dpks
        self._cartessian = cartessian
        self._magmom = magmom
    
    def get_pp(self):
        return self._pp

    def get_orb(self):
        return self._orb    

    def get_label(self):
        return self._label
    
    def get_cell(self,bohr = False):
        "return the cell matrix, in unit of Angstrom"
        cell = np.array(self._cell) * self._lattice_constance * 
    
    def set_pp(self,pplist):
        self.pp = pplist
    
    def set_orb(self,orblist):
        self.orb = orblist

class PrepareAbacus:
    def __init__(self,
                 save_path: str = Path("."),
                 input_template: str = None,
                 kpt_template: str = None,
                 stru_template: str = None,
                 mix_input: Dict[str,any] = {},
                 mix_kpt: List = List[Union(list,int)],
                 mix_stru: List = List[str],
                 pp_dict: Dict[str,str]= {},
                 orb_dict: Dict[str,str]= {},
                 extra_files = List[str]):
        """To prepare the inputs of abacus

        Parameters
        ----------
        save_path : str, optional
            the path to store inputs, by default Path(".")
        input_template : str, optional
            file name of INPUT template, by default None
        kpt_template : str, optional
            file name of KPT template, by default None
        stru_template : str, optional
            file name of STRU template, by default None
        mix_input : Dict[str,any], optional
            mix setting of input parameters, by default {}. 
            Such as:{"ecutwfc":[50, 60, 70], 
                     "mixing_beta":[0.3, 0.4, 0.5]}, 
            then will prepare 3*3 = 9 inputs.
        mix_kpt : List, optional
            mixing setting of k point list, by default []. If this parameter
            is defined, then kpt_template will be invalide.
            Such as: [[1,1,1,0,0,0],
                      [2,2,2,0,0,0]], 
            then will prepare 2 KPT files.
            The kpoint list can be a list of 6 values (3 int number for K point
            and 3 float values for the shift value) like: [1,1,1,0,0,0], 
            or 3 values like: [1,1,1] which means [1,1,1,0,0,0], or one int number
            like: 3 which means [3,3,3,0,0,0].
        mix_stru : List, optional
            a list of STRU files, by default []. If this parameter is defined
            then stru_template will be invalide.
            Such as: ["a/STRU","b/STRU"]
        pp_dict : List, optional
            a dictionary specify the pseudopotential files, by default {}
        orb_dict : List, optional
            a dictionary specify the orbital files, by default {}
        extra_files : List, optional
            a list of somes extra files that will be putted in each input, by default []   
        """        
        self.save_path = save_path
        self.input_template = input_template
        self.kpt_template = kpt_template
        self.stru_template = stru_template
        self.mix_input = mix_input
        self.pp_dict = pp_dict
        self.orb_dict = orb_dict
        self.extra_files = extra_files
        self.mix_kpt = mix_kpt
        self.mix_stru = mix_stru
        
        self.input_list = self.Construct_input_list()
        self.kpt_list = self.Construct_kpt_list()
        self.stru_list = self.Construct_stru_list()
    
    def Construct_input_list(self):
        all_inputs = []
        if self.input_template == None:
            input_constant = {}
        else:
            input_constant = PrepareAbacus.ReadInput(self.input_template)
        
        list_param = {}
        for k,v in self.mix_input.items():
            #if value is list type, then we need prepare INPUT for each value
            #input_constant stores the parameters whose value is constant
            #list_param stores the parameters that has several values.
            if isinstance(v,(int,float,str)):
                input_constant[k] = v
            elif isinstance(v,list):
                list_param[k] = v
            else:
                print("NOTICE: type of '%s' is" % str(v),type(v))
                input_constant[k] = v
        
        all_inputs.append(input_constant)
        for k,v in list_param.items():
            for i,iv in enumerate(v):
                if i == 0:
                    #if v only one value, add in all_inputs directly
                    for iinput in all_inputs:
                        iinput[k] = iv   
                else:
                    #if v has more than two values, deepcopy all_inputs and modify v, and then add to all_inputs
                    if i == 1:  
                        tmp_inputs = copy.deepcopy(all_inputs)
                    for iinput in tmp_inputs:
                        iinput[k] = iv
                    all_inputs.append(copy.deepcopy(tmp_inputs))
        
        return all_inputs

    def Construct_kpt_list(self):
        """all_kpt is a list of list (6 element) or filename of KPT"""
        all_kpt = []
        if len(self.mix_kpt) == 0:
            if self.kpt_template != None:
                if os.path.isfile(self.kpt_template):
                    all_kpt.append(os.path.abspath(self.kpt_template))
                else:
                    print("Has define the KPT file '%s', but file is not exist!!!" % self.kpt_template)
            return all_kpt
        else:
            for ikpt in self.mix_kpt:
                if isinstance(ikpt,int):
                    all_kpt.append([ikpt,ikpt,ikpt,0,0,0])
                elif isinstance(ikpt,list):
                    if len(ikpt) == 3:
                        all_kpt.append(ikpt+[0,0,0])
                    elif len(ikpt) == 6:
                        all_kpt.append(ikpt)
                    else:
                        print("mix_kpt should be a int or list of 3/6 element, but not ",ikpt)
                else:
                    print("element of mix_kpt should be int or list of 3/6 elements, but not ",ikpt)
        return all_kpt
    
    def Construct_stru_list(self):
        all_stru = []
        if len(self.mix_stru) == 0:
            if self.stru_template != None:
                if os.path.isfile(self.stru_template):
                    all_stru.append(os.path.abspath(self.stru_template))
                else:
                    print("Has define the STRU file '%s', but file is not exist!!!" % self.stru_template)
            return all_stru
        else:
            for istru in self.mix_stru:
                if os.path.isfile(istru):
                    all_stru.append(os.path.abspath(istru))
                else:
                    print("Structure file '%s' is not exist" % istru)
        return all_stru
    
    def prepare(self):
        if not self.stru_list:
            print("No stru files, skip!!!")
            return False
        
        if not self.kpt_list:
            print("WARNING: not set KPT")
            kpt_list = [None]
        else:
            kpt_list = self.kpt_list
        
        if not self.input_list:
            input_list = [None]
        else:
            input_list = self.input_list 
                   
        ipath = 0
        for istru in self.stru_list: 
            stru_data = PrepareAbacus.ReadStru(istru)
            for ikpt in kpt_list:
                for iinput in input_list:
                    #create folder
                    save_path = os.path.join(self.save_path,ipath.zfill(5))
                    if not os.path.isdir(save_path):
                        os.makedirs(save_path)
                    
                    #create INPUT   
                    if iinput != None: 
                        PrepareAbacus.WriteInput(iinput,os.path.join(save_path,"INPUT"))
                    
                    #create KPT
                    if ikpt != None:
                        kptf = os.path.join(save_path,"KPT")
                        if isinstance(ikpt,str):
                            if os.path.isfile(kptf):
                                os.unlink(kptf)
                            os.symlink(os.path.abspath(ikpt),kptf)
                        elif isinstance(ikpt,list):
                            PrepareAbacus.WriteKpt(ikpt,kptf)
                    
                    #create STRU
                    struf = os.path.join(os.path.join(save_path,"STRU"))
                    if os.path.isfile(struf):
                        os.unlink(struf)
                    os.symlink(os.path.abspath(istru),struf)
                    
                    #create pp/orb
                    for ielement in element:
                        if ielement in self.pp_dict:
                            ppf = os.path.join(save_path, os.path.split(self.pp_dict[ielement])[1])
                            if os.path.isfile(ppf):
                                os.unlink(ppf)
                            os.symlink(os.path.abspath(self.pp_dict[ielement]),ppf)
                        else:
                            print("WARING: %s: element '%s' is found in STRU, but not specified the pseudopotential file" % (save_path, ielement))
                        
                        if ielement in self.orb_dict:
                            orbf = os.path.join(save_path, os.path.split(self.orb_dict[ielement])[1])
                            if os.path.isfile(orbf):
                                os.unlink(orbf)
                            os.symlink(os.path.abspath(self.orb_dict[ielement]),orbf)
                        else:
                            print("WARING: %s: element '%s' is found in STRU, but not specified the pseudopotential file" % (save_path, ielement))
                        
                    for ipp in self.pp_dict:
                        

    @staticmethod
    def ReadStru(stru:str = "STRU") -> AbacusStru:
        "read the label, pp, orb, cell, coord, deepks-descriptor"
        if not os.path.isfile(stru):
            return None

        with open(stru) as f1: lines = f1.readlines()
        atomic_species = []
        numerical_orbital = []
        lattice_constant = 1.0
        lattice_vector = []
        atom_positions = []
        
        def get_block(keyname):
            block = []
            for i,line in enumerate(lines):
                if line.strip() == "": continue
                elif line.strip() == keyname:
                    for ij in range(i+1,len(lines)):
                        if lines[ij].strip() == "": continue
                        elif lines[ij].strip() in ABACUS_STRU_KEY_WORD:
                            return block
                        else:
                            block.append(lines[ij])
                    break
            return None
        
        atomic_species = get_block("ATOMIC_SPECIES")
        numerical_orbital = get_block("NUMERICAL_ORBITAL")
        lattice_constant = get_block("LATTICE_CONSTANT")
        lattice_vector = get_block("LATTICE_VECTORS")
        atom_positions = get_block("ATOMIC_POSITIONS")
        dpks = get_block("NUMERICAL_DESCRIPTOR")
        lattice_constant = 1.0 if lattice_constant == None else float(lattice_constant[0].split()[0]) 
        dpks = None if dpks == None else dpks[0].strip()
        
        #read species
        pp = []
        labels = []
        mass = []
        for line in atomic_species:
            sline = line.split()
            labels.append(sline[0])
            pp.append(sline[2])
            mass.append(float(sline[1]))
            
        #read orbital
        if numerical_orbital == None:
            orb = None
        else:
            orb = []
            for line in numerical_orbital:
                orb.append(line.split()[0])

        #read cell
        cell = []
        for line in lattice_vector:
            cell.append([float(i) for i in line.split()[:3]])

        #read coordinate and coordinate type and atom number of each type
        atom_number = []
        coords = []
        magmom = []
        coord_type = atom_positions[0].strip().lower()
        if coord_type == "direct":
            cartessian = False
        elif coord_type == "cartessian":
            cartessian = True
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
                          lattice_constance=lattice_constant,
                          magmom=magmom,
                          dpks=dpks,
                          cartessian=cartessian)    
            

    @staticmethod
    def WriteKpt(kpoint_list:List = [1,1,1,0,0,0],file_name:str = "KPT"):
        with open(file_name,'w') as f1:
            f1.write("K_POINTS\n0\nGamma\n")
            f1.write(" ".join([str(i) for i in kpoint_list]))
    
    @staticmethod
    def ReadInput(INPUTf: str = None, input_lines: str = None) -> Dict[str,any]:
        #read the INPUT file
        input_context = {}
        if INPUTf != None:
            if not os.path.isfile(INPUTf):
                print("Can not find the file '%s'" % INPUTf)
                return input_context
            with open(INPUTf) as f1: input_lines = f1.readlines()
        if input_lines == None:
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
        
        for i,iline in enumerate(input_lines):
            if iline.strip() == 'INPUT_PARAMETERS':
                readinput = True
            elif iline.strip() == '' or iline.strip()[0] in ['#']:
                continue
            elif readinput:
                sline = re.split('[ \t]',iline.split("#")[0],maxsplit=1)
                if len(sline) == 2:
                    input_context[sline[0].lower().strip()] = str2intfloat(sline[1].strip())
        return input_context

    @staticmethod
    def WriteInput(input_context:Dict[str,any],
                   INPUTf: str = "INPUT"):
        out = "INPUT_PARAMETERS\n"
        for k,v in input_context.items():
            out += "%s\t%s\n" % (str(k),str(v))
        with open(INPUTf,'w') as f1: f1.write(out)            


def PrepareInput(param):
    param_setting_file = param.param
    save_folder = param.save
    
    if not os.path.isfile(param_setting_file):
        print("Can not find file %s!!!" % param_setting_file)
        sys.exit(1)

    param_setting = json.load(open(param_setting_file))
    
    

def PrepareArgs(parser):  
    parser.description = "This script is used to prepare the INPUTS OF ABACUS JOBS"
    parser.add_argument('-p', '--param', type=str, help='the parameter file, should be .json type',required=True)
    parser.add_argument('-s', '--save', type=str,  default="abacustest",help='where to store the inputs, default is abacustest ')
    return parser

def main():
    parser = argparse.ArgumentParser()
    param = PrepareArgs(parser)
    Prepare(param)
    
if __name__ == "__main__":
    main()
