from typing import Dict, List
import numpy as np
import re
import sys
from .. import constant


class QeInputs:
    def __init__(self,
                 param:Dict[str,any],
                 label:List[str],
                 atom_number:List[int],
                 cell:List[List[float]] ,
                 coord:List[List[float]] ,
                 pp:List[str],
                 kpt:List[any] = None,
                 cell_type:str = "alat",
                 coord_type: str = "alat",
                 kpt_type = "automatic",
                 mass:List[float] = None,
                 element:List[str] = None,  
    ):
        """QE INPUT class, the unit is Bohr
        
        Parameters
        ----------
        param : Dict[str,any]
            the input parameters of QE.
            key is the card name: control, system, electrons, ions, cell
            value is a dictionary of the card parameters.
        label : List[str]
            the label name of each type; 
        atom_number : List[int]
            the atom number of each type
        pp : List[str]
            pseudopotential file of each type
        coord : List[List[float]]
            coordinate of each atom
        cell : List[List[float]]
            cell of a, b, c 
        kpt : List[any], optional
            k point list, by default None, which means [1,1,1,0,0,0] and kpt_type is automatic.
        cell_type : str, optional
            the type of cell, by default "alat"
        coord_type : str, optional
            the type of coordinate, by default "alat"
        kpt_type : str, optional    
            the type of k point, by default "automatic"
        mass : List[float], optional
            mass of each type, by default None
        element : List[str], optional
            element name of each type, by default None
        """
 
        self.label = label
        self.atom_number = atom_number
        self.cell = cell
        self.coord = coord
        self.pp = pp
        self.kpt = kpt
        self.cell_type = cell_type
        self.coord_type = coord_type
        self.kpt_type = kpt_type
        
        self.param = {}
        comm_cards = ["control","system","electrons","ions","cell"]
        for ik,iv in param.items():
            if ik.lower() in comm_cards:
                self.param[ik.lower()] = iv if iv != None else {}
            else:
                self.param[ik] = iv if iv != None else []
        for ikey in comm_cards:
            if ikey not in self.param:
                self.param[ikey] = {}
                
        if element == None:
            print("WARNING: element is not defined, will use the first one/two letter of label as element name")
            self.element = []
            for i in label:
                while i[-1].isdigit(): i = i[:-1]
                self.element.append(i)
        else:
            self.element = element
        
        if mass == None:
            self.mass = [constant.MASS_DICT.get(i,1.0) for i in self.element]
        else:
            self.mass = mass
        
        if kpt == None:
            self.kpt = [1,1,1,0,0,0]
            self.kpt_type = "automatic"

    def get_param(self):
        return self.param
    
    def get_label(self):
        # return the label name of each atom
        label = []
        for idx,i in enumerate(self.atom_number):
            label += [self.label[idx]] * i
        return label
    
    def get_pp(self):
        return self.pp
    
    def get_mass(self):
        return self.mass
    
    def get_kpt(self):
        return self.kpt
    
    def get_kpt_type(self):
        return self.kpt_type
    
    def get_cell(self,bohr = False):
        if self.cell_type == "bohr":
            transfer_unit = 1
        elif self.cell_type == "angstrom":
            transfer_unit = constant.A2BOHR
        elif self.cell_type == "alat":
            transfer_unit = float(self.param["system"]["celldm(1)"])
        else:
            print("ERROR: cell_type is not defined")
            sys.exit(1)
        if not bohr:
            transfer_unit *= constant.BOHR2A
            
        cell = np.array(self.cell) * transfer_unit
        return cell.tolist()
    
    def get_coord(self,bohr = False, direct=False):
        '''return the coordinate matrix, in cartesian or direct type, in unit of Angstrom or Bohr.
        If direct is True, then return the direct type coordinate, and bohr will be ignored.'''
        
        if direct and self.coord_type == "crystal":
            return self.coord
        
        if not direct:
            if self.coord_type == "crystal":
                return np.array(self.coord).dot(np.array(self.get_cell(bohr=bohr))).tolist()
            else:
                if self.coord_type == "bohr":
                    transfer_unit = 1
                elif self.coord_type == "angstrom":
                    transfer_unit = constant.A2BOHR
                elif self.coord_type == "alat":
                    transfer_unit = float(self.param["system"]["celldm(1)"])
                if not bohr:
                    transfer_unit *= constant.BOHR2A
                return (np.array(self.coord) * transfer_unit).tolist()
        else:
            if self.coord_type == "crystal":
                return self.coord
            else:
                coord_bohr = self.get_coord(bohr=True,direct=False)
                cell_bohr = self.get_cell(bohr=True)
                return np.array(coord_bohr).dot(np.linalg.inv(np.array(cell_bohr))).tolist()
    
    def _str_value(self,value):
        # if value is float or int then return str(value)
        # if value is string and not ".true.", ".false." then return "'%s'" % value
        
        if isinstance(value,str):
            try:
                float(value)
                return value
            except:
                if value.lower() in [".true.",".false."]:
                    return value
                else:
                    return "'%s'" % value
        else:
            return str(value)
    
    def write(self,filename="input"):
        with open(filename,"w") as f:
            for param_key in ["control","system","electrons","ions","cell"]:
                f.write("&%s\n" % param_key.upper())
                for key,value in self.param.get(param_key,{}).items():
                    f.write("%17s = %s\n" % (key,self._str_value(value)))
                f.write("/\n\n")
            
            f.write("\nCELL_PARAMETERS (%s)\n" % self.cell_type)
            for i in self.cell:
                f.write("%17.11f %17.11f %17.11f\n" % tuple(i))
            
            f.write("\nATOMIC_SPECIES\n")
            for i,ilabel in enumerate(self.label):
                f.write("%s %f %s\n" % (ilabel,self.mass[i],self.pp[i]))
            
            f.write("\nATOMIC_POSITIONS (%s)\n" % self.coord_type)
            coord_idx = 0
            for i,ilabel in enumerate(self.label):
                for j in range(self.atom_number[i]):
                    f.write("%5s %17.11f %17.11f %17.11f\n" % tuple([ilabel] + self.coord[coord_idx]))
                    coord_idx += 1
                    
            f.write("\nK_POINTS %s\n" % self.kpt_type)
            if self.kpt_type == "automatic":
                kpoints = self.kpt[:3]
                offset = [0 if ioff == 0 else 1 for ioff in self.kpt[3:]]
                f.write("%d %d %d %d %d %d\n" % tuple(kpoints + offset))   
            elif self.kpt_type == "gamma":
                f.write("0\n")
            else:
                print("ERROR: kpt_type is not defined")
                sys.exit(1)
            
            # write special card
            for key,value in self.param.items():
                if key.lower() not in ["control","system","electrons","ions","cell"]:
                    if isinstance(value,dict):
                        f.write("\n&%s\n" % key)
                        for ik,iv in value.items():
                            f.write("%17s = %s\n" % (ik,str(iv)))
                        f.write("/\n\n")
                    elif isinstance(value,list):
                        f.write("\n%s\n" % key)
                        for iv in value:
                            f.write("%s\n" % str(iv))
                        f.write("\n")
                    else:
                        print("ERROR: special card '%s' is not defined" % key)
                        sys.exit(1)
            
        
    @staticmethod
    def ReadInput(filename="input"):
        def get_card_type(line):
            line = re.split("[#!]",line)[0]
            if "(" in line:
                return line.split("(")[1].split(")")[0].strip().lower()
            else:
                return line.split()[1].lower()
            
        with open(filename) as f:
            lines = f.readlines()
        param = {}
        label = []
        atom_number = []
        cell = []
        coord = []
        pp = []
        kpt = None
        cell_type = "alat"
        coord_type = "alat"
        kpt_type = "automatic"
        mass = []
        element = []
        for i, line in enumerate(lines):
            if line.strip() == "" or line.strip()[0] in ["!", "#"]: continue
            elif line.strip()[0] == "&":
                block = line.split()[0].strip()[1:].lower()
                if block in ["control","system","electrons","ions","cell"]:
                    param[block] = {}
                    for ij in range(i+1,len(lines)):
                        if lines[ij].strip() == "" or lines[ij].strip()[0] in ["!","#"]: continue
                        elif lines[ij].strip() == "/": break
                        else:
                            for iline in re.split("[#!]",lines[ij])[0].split(","):
                                if iline.strip() == "":continue
                                sline = iline.split("=")
                                param[block][sline[0].strip()] = sline[1].strip()
            elif line.split("#")[0].split()[0] == "CELL_PARAMETERS":
                cell_type = get_card_type(line)
                for ij in range(i+1,len(lines)):
                    if lines[ij].strip() == "" or lines[ij].strip()[0] in ["!","#"]: continue
                    else:
                        cell.append([float(i) for i in lines[ij].split()[:3]])
                        if len(cell) == 3: break
            elif line.split("#")[0].split()[0] == "ATOMIC_SPECIES":
                for ij in range(i+1,len(lines)):
                    if lines[ij].strip() == "" or lines[ij].strip()[0] in ["!","#"]: continue
                    else:
                        sline = lines[ij].split("#")[0].split()
                        label.append(sline[0].strip())
                        atom_number.append(0)
                        mass.append(float(sline[1].strip()))
                        pp.append(sline[2].strip())
                        if len(label) == int(param["system"]["ntyp"]): break
            elif line.split("#")[0].split()[0] == "ATOMIC_POSITIONS":
                coord_type = get_card_type(line)
                label_coord = {}
                for ij in range(i+1,len(lines)):
                    if lines[ij].strip() == "" or lines[ij].strip()[0] in ["!","#"]: continue
                    elif lines[ij].strip() in ["crystal","alat","bohr","angstrom"]: continue
                    else:
                        sline = lines[ij].split("#")[0].split()
                        ilabel = sline[0].strip()
                        if ilabel not in label:
                            print("ERROR: label '%s' is not matched that in ATOMIC_SPECIES" % ilabel)
                            sys.exit(1)
                        atom_number[label.index(ilabel)] += 1
                        if ilabel not in label_coord:
                            label_coord[ilabel] = []
                        label_coord[ilabel].append([float(i) for i in sline[1:4]])
                        if sum(atom_number) == int(param["system"]["nat"]): break
                for ilabel in label:
                    coord += label_coord[ilabel]
            elif line.split("#")[0].split()[0] == "K_POINTS":
                kpt_type = get_card_type(line)
                if kpt_type == "automatic":
                    kpt = [int(i) if idx < 3 else float(i) for idx,i in enumerate(lines[lines.index(line)+1].split()[:6])]
                elif kpt_type == "crystal":
                    kpt = [int(i) for i in lines[lines.index(line)+1].split()]
                elif kpt_type == "gamma":
                    kpt = [0]
                else:
                    print("ERROR: kpt_type is not defined")
                    sys.exit(1)
        return QeInputs(param=param,
                        label=label,
                        atom_number=atom_number,
                        cell=cell,
                        coord=coord,
                        pp=pp,
                        kpt=kpt,
                        cell_type=cell_type,
                        coord_type=coord_type,
                        kpt_type=kpt_type,
                        mass=mass,
                        element=element)
