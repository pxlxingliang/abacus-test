import os,sys,argparse,glob,json,traceback,re,shutil
import numpy as np
from typing import Dict, List, Optional, Union
import copy
from pathlib import Path
from abacustest import constant

class AbacusStru:
    def __init__(self,
                 label:List[str],
                 atom_number:List[int],
                 cell:List[List[float]] ,
                 coord:List[List[float]] ,
                 pp:List[str],
                 orb:List[str] = None,
                 mass:List[float] = None,
                 element:List[str] = None,
                 lattice_constant:float = 1,
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
        self._orb = orb if orb else []
        self._lattice_constant = lattice_constant
        self._dpks = dpks if dpks else None
        self._cartessian = cartessian
        self._magmom = magmom if magmom else []
        
        if element != None:
            self._element = element
        else:
            self._element = []
            for i in label:
                while i[-1].isdigit(): i = i[:-1]
                self._element.append(i)
        
        if mass != None:
            self._mass = mass
        else:
            self._mass = [constant.MASS_DICT.get(i,1.0) for i in self._element]
    
    def get_pp(self):
        return self._pp

    def get_orb(self):
        return self._orb    

    def get_label(self):
        return self._label
    
    def get_dpks(self):
        return self._dpks
    
    def get_mass(self):
        return self._mass
    
    def get_element(self):
        return self._element
    
    def get_cell(self,bohr = False):
        "return the cell matrix, in unit of Angstrom"
        transfer_unit = 1 if bohr else constant.BOHR2A
        cell = np.array(self._cell) * self._lattice_constant * transfer_unit
        return cell.tolist()
    
    def set_pp(self,pplist):
        self._pp = pplist
    
    def set_orb(self,orblist):
        self._orb = orblist

    def set_dpks(self,descriptor):
        self._dpks = descriptor
    
    def set_mass(self,mass):
        self._mass = mass
        
    def set_element(self,element):
        self._element = element
        
    def write(self,struf="STRU"):
        cc = ""
        #write species
        cc += "ATOMIC_SPECIES\n"
        for i,ilabel in enumerate(self._label):
            cc += "%s %f %s\n" % (ilabel,self._mass[i],self._pp[i])
        
        #write orb
        if self._orb:
            cc += "\nNUMERICAL_ORBITAL\n"
            for i in self._orb:
                cc += i + "\n"
        
        #write  LATTICE_CONSTANT
        cc += "\nLATTICE_CONSTANT\n%f\n" % self._lattice_constant

        #write LATTICE_VECTORS
        cc += "\nLATTICE_VECTORS\n"
        for i in self._cell:
            cc += "%17.11f %17.11f %17.11f\n" % tuple(i)
        
        #write ATOMIC_POSITIONS
        cc += "\nATOMIC_POSITIONS\n"
        if self._cartessian:
            cc += "Cartessian\n"
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

    @staticmethod
    def ReadStru(stru:str = "STRU"):
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
                elif line.split('#')[0].strip() == keyname:
                    for ij in range(i+1,len(lines)):
                        if lines[ij].strip() == "": continue
                        elif lines[ij].strip() in constant.ABACUS_STRU_KEY_WORD:
                            return block
                        else:
                            block.append(lines[ij])
                    return block
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
        coord_type = atom_positions[0].split("#")[0].strip().lower()
        if coord_type == "direct":
            cartessian = False
        elif coord_type == "cartesian":
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
                          lattice_constant=lattice_constant,
                          magmom=magmom,
                          dpks=dpks,
                          cartessian=cartessian)
    
        
class PrepareAbacus:
    def __init__(self,
                 save_path: str = Path("."),
                 example_template: str = None,
                 input_template: str = None,
                 kpt_template: str = None,
                 stru_template: str = None,
                 mix_input: Dict[str,any] = {},
                 mix_kpt: List[Union[int,List[int]]] = [],
                 mix_stru: List[str]=[],
                 pp_dict: Dict[str,str]= {},
                 orb_dict: Dict[str,str]= {},
                 pp_path:str = None,
                 orb_path:str = None,
                 dpks_descriptor:str=None,
                 extra_files: List[str] = []):
        """To prepare the inputs of abacus

        Parameters
        ----------
        save_path : str, optional
            the path to store inputs, by default Path(".")

        example_template : str, optional
            folder name of an abacus inputs template, by default None.
            If specify input/kpt/stru_template, the related file in exmaple_template will
            be ingored. The pp/orb defined in example_template/STRU will also added to 
            pp/orb_dict.

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

        pp_dict : Dict, optional
            a dictionary specify the pseudopotential files, by default {}

        orb_dict : Dict, optional
            a dictionary specify the orbital files, by default {}

        pp_path : str, optional
            the path of your pseudopotential lib, by default None. 
            The file name shuold be started with the element name and followed by "_".
            Such as: Ag_ONCV_PBE-1.0.upf.

        orb_path : str, optional
            the path of your orbital lib, by default None. 
            The file name shuold be started with the element name and followed by "_".
            Such as: Ag_gga_7au_100Ry_4s2p2d1f.orb.

        dpks_descriptor : str, optional
            file name of deepks descriptor, by default None

        extra_files : List, optional
            a list of somes extra files that will be putted in each input, by default []   
        """        
        self.save_path = save_path
        self.example_template = example_template
        self.input_template = input_template
        self.kpt_template = kpt_template
        self.stru_template = stru_template
        self.mix_input = mix_input
        self.pp_dict = pp_dict
        self.orb_dict = orb_dict
        self.extra_files = extra_files if extra_files else []
        self.mix_kpt = mix_kpt
        self.mix_stru = mix_stru
        self.dpks_descriptor = dpks_descriptor if dpks_descriptor else None
        
        self.pp_path = None if not pp_path else pp_path.strip()
        self.orb_path = None if not orb_path else orb_path.strip()
        self.CollectPP()
        self.CollectOrb()

        self.input_list,self.input_mix_param = self.Construct_input_list()
        self.kpt_list = self.Construct_kpt_list()
        self.stru_list = self.Construct_stru_list()
    
    def CollectPP(self):
        if not self.pp_path:
            return
        if os.path.isdir(self.pp_path):
            allfiles = os.listdir(self.pp_path)
            for ifile in allfiles:
                if not os.path.isfile(os.path.join(self.pp_path,ifile)): continue
                element_name = ifile.split("_")[0]
                if element_name not in self.pp_dict:
                    self.pp_dict[element_name] = os.path.join(self.pp_path,ifile)
        else:
            print("Not find pp dir: %s" % self.pp_path)

    def CollectOrb(self):
        if not self.orb_path:
            return
        if os.path.isdir(self.orb_path):
            allfiles = os.listdir(self.orb_path)
            for ifile in allfiles:
                if not os.path.isfile(os.path.join(self.orb_path,ifile)): continue
                element_name = ifile.split("_")[0]
                if element_name not in self.orb_dict:
                    self.orb_dict[element_name] = os.path.join(self.orb_path,ifile)
        else:
            print("Not find orb dir: %s" % self.orb_path)

    def Construct_input_list(self):
        all_inputs = []
        inputf = None
        if self.example_template != None:
            if os.path.isfile(os.path.join(self.example_template,"INPUT")):
                inputf = os.path.join(self.example_template,"INPUT")
        if self.input_template != None:
            if os.path.isfile(self.input_template):
                inputf = self.input_template 
            else:
                print("WARNING: File %s is not found" % self.input_template)

        print("INPUT template file: %s" % str(inputf))
        if inputf != None:
            input_constant = PrepareAbacus.ReadInput(inputf)
        else:
            print("Please define the INPUT template!!!")
            return None,None 
        
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
                print("WARNING: type of '%s' is" % str(v),type(v),"will not add to INPUT")
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
                    all_inputs += copy.deepcopy(tmp_inputs)

        return all_inputs,[i for i in list_param.keys()]

    def Construct_kpt_list(self):
        """
        all_kpt is a list of list (6 element) or filename of KPT.
        If not specify mix_kpt and has specify the KPT template, will retrun a list
        of KPT filename, or will return a null list if no KPT template.
        """
        all_kpt = []

        kptf = None
        if self.example_template != None:
            if os.path.isfile(os.path.join(self.example_template,"KPT")):
                kptf = os.path.join(self.example_template,"KPT")
        if self.kpt_template != None:
            if os.path.isfile(self.kpt_template):
                kptf = self.kpt_template 
            else:
                print("WARNING: File %s is not found" % self.kpt_template)

        print("KPT template file: %s" % str(kptf))

        if len(self.mix_kpt) == 0:
            if kptf != None:
                all_kpt.append(kptf)
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
                elif isinstance(ikpt,str):
                    for iikpt in glob.glob(ikpt):
                        all_kpt.append(iikpt)
                else:
                    print("element of mix_kpt should be int or list of 3/6 elements, but not ",ikpt)
        return all_kpt
    
    def Construct_stru_list(self):
        all_stru = []
        struf = None
        if self.example_template != None:
            if os.path.isfile(os.path.join(self.example_template,"STRU")):
                struf = os.path.join(self.example_template,"STRU")
        if self.stru_template != None:
            if os.path.isfile(self.stru_template):
                struf = self.stru_template 
            else:
                print("WARNING: File %s is not found" % self.stru_template)

        print("STRU template file: %s" % str(struf))

        if len(self.mix_stru) == 0:
            if struf != None:
                all_stru.append(struf)
            else:
                print("Please set the stru_template!")
            return all_stru
        else:
            for stru in self.mix_stru:
                allstrus = glob.glob(stru)
                for istru in allstrus:
                    all_stru.append(os.path.abspath(istru))
                if len(allstrus) == 0:
                    print("Structure file '%s' is not exist" % stru)
        return all_stru
    
    def prepare(self):
        "Return a dict, and the key is the path of inputs,"
        "and the value is list of structure, kpt, and input setting information."
        "The input setting information is a dict of param name and value"

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
                   
        ipath = -1
        param_setting = {}
        has_create_savepath = False
        for istru in self.stru_list:  #iteration of STRU
            cwd = os.getcwd()
            stru_data = AbacusStru.ReadStru(istru)
            stru_path = os.path.split(istru)[0]
            labels = stru_data.get_label()
            linkstru = True
            skipstru = False
            allfiles = []  #files that will be linked
            pp_list = []   #pp file name
            orb_list = []  #orb file name
            dpks = None    #dpks file name
            for i,ilabel in enumerate(labels):
                #check pp file 
                if ilabel not in self.pp_dict:
                    #print("label '%s' is found in '%s', but not defined in pp_dict." % (ilabel,istru))
                    os.chdir(stru_path)
                    pp_in_stru = stru_data.get_pp()[i]
                    if os.path.isfile(pp_in_stru):
                        print("label '%s': link the pseudopotential file '%s' defined in %s" % (ilabel,pp_in_stru,istru))
                        pp_list.append(os.path.split(pp_in_stru)[1])
                        allfiles.append(os.path.abspath(pp_in_stru))
                    else:
                        print("label '%s': the pseudopotential file '%s' defined in %s is not found, skip this structure" % (ilabel,pp_in_stru,istru))
                        skipstru = True
                        break
                    os.chdir(cwd)
                else:
                    pp_list.append(os.path.split(self.pp_dict[ilabel])[1])  #only store the file name to pp_list
                    allfiles.append(self.pp_dict[ilabel]) #store the whole pp file to allfiles                    
                if not stru_data.get_pp() or pp_list[-1] != stru_data.get_pp()[i]:
                    linkstru = False

                #check orbital file    
                if ilabel not in self.orb_dict:
                    if stru_data.get_orb():
                        #print("label '%s' is found in '%s', but not defined in orb_dict." % (ilabel,istru))
                        os.chdir(stru_path)
                        orb_in_stru = stru_data.get_orb()[i]
                        if os.path.isfile(orb_in_stru):
                            print("label '%s': link the orb file '%s' defined in %s" % (ilabel,orb_in_stru,istru))
                            orb_list.append(os.path.split(orb_in_stru)[1])
                            allfiles.append(os.path.abspath(orb_in_stru))
                            if orb_list[-1] != stru_data.get_orb()[i]:  
                                linkstru = False
                        else:
                            print("label '%s': the orbital file '%s' defined in %s is not found." % (ilabel,orb_in_stru,istru))
                        os.chdir(cwd)
                else:
                    orb_list.append(os.path.split(self.orb_dict[ilabel])[1]) 
                    allfiles.append(self.orb_dict[ilabel])
                    if not stru_data.get_orb() or orb_list[-1] != stru_data.get_orb()[i]:  
                        linkstru = False

            if skipstru:
                continue
            
            #check dpks 
            if self.dpks_descriptor:
                if os.path.isfile(self.dpks_descriptor):
                    dpks = os.path.split(self.dpks_descriptor)[1]
                    allfiles.append(self.dpks_descriptor)
                else:
                    print("Error: Can not find file %s, skip the prepare of dpks_descriptor" % self.dpks_descriptor)
                    dpks = None
            else:
                os.chdir(stru_path)
                dpks = stru_data.get_dpks()
                if dpks:
                    if os.path.isfile(dpks):
                        allfiles.append(os.path.abspath(dpks))
                        dpks = os.path.split(dpks)[1]
                    else:
                        print("Error: deepks descriptor is defined in %s/STRU, but can not find the file, skip the prepare" % stru_path)
                        dpks = None
                else:
                    dpks = None
                os.chdir(cwd)
            if dpks != stru_data.get_dpks():
                linkstru = False

            #stru_data will be writen in new folder
            stru_data.set_pp(pp_list)
            stru_data.set_orb(orb_list)
            stru_data.set_dpks(dpks)
            
            for ikpt in kpt_list:  #iteration of KPT 
                for iinput in input_list: #iteration of INPUT
                    ipath += 1
                    #create folder
                    if not has_create_savepath:
                        if os.path.isdir(self.save_path):
                            from myflow.comm import GetBakFile
                            bk = GetBakFile(self.save_path)
                            shutil.move(self.save_path,bk)
                        has_create_savepath = True
                    save_path = os.path.join(self.save_path,str(ipath).zfill(5))

                    #store the param setting and path
                    param_setting[save_path] = [istru,ikpt,{}]
                    for input_param in self.input_mix_param:
                        param_setting[save_path][-1][input_param] = iinput.get(input_param)

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
                    if linkstru:
                        os.symlink(os.path.abspath(istru),struf)
                    else:
                        stru_data.write(struf)
                    
                    #link other files
                    for ifile in allfiles:
                        ifile = os.path.abspath(ifile)
                        filename = os.path.split(ifile)[1]
                        target_file = os.path.join(save_path,filename)
                        if os.path.isfile(target_file):
                            os.unlink(target_file)
                        os.symlink(ifile,target_file)         
        
        return param_setting

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
            if v != None:
                out += "%s\t%s\n" % (str(k),str(v))
            else:
                out += "#%s\t \n" % (str(k))

        with open(INPUTf,'w') as f1: f1.write(out)            


def DoPrepare(param_setting: Dict[str, any], save_folder: str) -> List[Dict[str, dict]]:  
    """
    param_setting is a dictionary like:
    {
        "example_template":["example_path"]
        "input_template":"INPUT",
        "kpt_template":"KPT",
        "stru_template":"STRU",
        "mix_input":{
            "ecutwfc":[50,60,70],
            "kspacing":[0.1,0.12,0.13]
        },
        "mix_kpt":[],
        "mix_stru":[],
        "pp_dict":{},
        "orb_dict":{},
        "pp_path": str,
        "orb_path": str,
        "dpks_descriptor":"",
        "extra_files":[],
        "mix_input_comment":"Do mixing for several input parameters. The diffrent values of one parameter should put in a list.",
        "mix_kpt_comment":"If need the mixing of several kpt setting. The element should be an int (like 2 means [2,2,2,0,0,0]), or a list of 3 or 6 elements (like [2,2,2] or [2,2,2,0,0,0]).",
        "mix_stru_commnet":"If need the mixing of several stru files. Like: ["a/stru","b/stru"],
    },

    save_folder specifies the destination of the new created examples.

    Return a list of dict, which is related to each example_template element. The key of the dict is the newcreated example path, and the value
    is the STRU/KPT/INPUT settings. 
    """
    example_template = param_setting.get("example_template",None)
    if isinstance(example_template,str):
        example_template = glob.glob(example_template)
    elif isinstance(example_template,list):
        tmp = copy.deepcopy(example_template)
        example_template = []
        for i in tmp:
            example_template += glob.glob(i)
    else:
        example_template = [example_template]
    
    all_path_setting = []
    for iexample in example_template:
        if len(example_template) > 1:
            save_path = os.path.join(save_folder,iexample)
            print("\n%s" % iexample)
        else:
            save_path = save_folder
        prepareabacus = PrepareAbacus(save_path=save_path,
                                  example_template=iexample,
                                  input_template=param_setting.get("input_template",None),
                                  kpt_template=param_setting.get("kpt_template",None),
                                  stru_template=param_setting.get("stru_template",None),
                                  mix_input=param_setting.get("mix_input",{}),
                                  mix_kpt=param_setting.get("mix_kpt",[]),
                                  mix_stru=param_setting.get("mix_stru",[]),
                                  pp_dict=param_setting.get("pp_dict",{}),
                                  orb_dict=param_setting.get("orb_dict",{}),
                                  pp_path=param_setting.get("pp_path",None),
                                  orb_path=param_setting.get("orb_path",None),
                                  dpks_descriptor=param_setting.get("dpks_descriptor",None),
                                  extra_files=param_setting.get("extra_files",[])
                                  )
        all_path_setting.append(prepareabacus.prepare())

    return all_path_setting

def PrepareInput(param):
    param_setting_file = param.param
    save_folder = param.save
    
    if not os.path.isfile(param_setting_file):
        print("Can not find file %s!!!" % param_setting_file)
        sys.exit(1)

    param_setting = json.load(open(param_setting_file))
    if "prepare" in param_setting:
        param_setting = param_setting["prepare"]
    
    all_path_setting = DoPrepare(param_setting,save_folder)
    for path_setting in all_path_setting:
        if path_setting:
            for k,v in path_setting.items():
                print("%s:%s" % (k,str(v)))

def PrepareArgs(parser):  
    parser.description = "This script is used to prepare the INPUTS OF ABACUS JOB"
    parser.add_argument('-p', '--param', type=str, help='the parameter file, should be .json type',required=True)
    parser.add_argument('-s', '--save', type=str,  default="abacustest",help='where to store the inputs, default is abacustest ')
    return parser

def main():
    parser = argparse.ArgumentParser()
    param = PrepareArgs(parser).parse_args()
    PrepareInput(param)
    
if __name__ == "__main__":
    main()
