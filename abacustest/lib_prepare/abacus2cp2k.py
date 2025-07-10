from .abacus import AbacusStru, ReadInput, ReadKpt
import os
from . import comm,cp2k
import copy

def update_dict(dict1, dict2):
    '''
    Update dict1 with dict2, if the key in dict1 is also in dict2, then update the value of dict1 with the value of dict2.
    If the value of one key in dict1 and dict2 are both dict, then merge the two dicts.
    
    Need check the case for section&value, such as KIND&element, which should be treated as a section of KIND.
    For one key not in dict1:
        If the key is type section&value, then
            If key is in repeat_section, then we add the section to dict1.
            else, check if key is in dict1, if yes, del the key in dict1, and add the new section to dict1.
        else
            add the key to dict1.
    
    For repeat_section,
    '''
    if (not isinstance(dict1,dict)) or (not isinstance(dict2,dict)):
        return dict2
    repeat_section = ["KIND"]
    for k,v in dict2.items():
        if k not in dict1:
            if "&" in k:
                real_k = k.split("&")[0]
                if real_k in repeat_section and v != None:
                    dict1[k] = v
                elif real_k in dict1:
                    if v != None:
                        dict1[k] = update_dict(dict1[real_k],v)
                    del dict1[real_k]
                else:
                    dict1_keys = [i.split("&")[0] for i in dict1.keys()]
                    if real_k not in dict1_keys:
                        if v != None:
                            dict1[k] = v
                    else:
                        idx = dict1_keys.index(real_k)
                        del dict1[list(dict1.keys())[idx]]
                        if v!= None:
                            dict1[k] = v
            else: 
                dict1_keys = [i.split("&")[0] for i in dict1.keys()]
                if k not in dict1_keys:
                    if v != None:  
                        dict1[k] = v
                else:
                    idx = dict1_keys.index(k)
                    if v != None:
                        dict1[k] = update_dict(dict1[list(dict1.keys())[idx]],v)
                    del dict1[list(dict1.keys())[idx]]
        else:
            if isinstance(v,dict) and isinstance(dict1[k],dict):
                dict1[k] = update_dict(dict1[k],v)
            elif v == None:
                del dict1[k]
            else:
                dict1[k] = v
    return dict1

def param2cp2k(param:dict, save_path:str=None):
    '''
    write the param dict to cp2k input file.
    
    The param dict is a nested dict, the key is the section name, the value is the parameters in this section.
    For each key and value:
        1. If the value is a dict, it is a section, and add "&" before the section name, and "&END" after the section.
            Special for the "&KIND" section, the key should be "KIND&element", and the value should be a dict, which contains the element, basis_set, potential,....
        2. If the value is a list, it is a list of values.
        3. If the value is a string or int or float, it is a parameter.
    '''
    def write_section(iparam, layer=0):
        cc = ""
        for k,v in iparam.items():
            if isinstance(v,dict):
                if "&" in k:
                    ks = k.split("&")
                    cc += "  "*layer + f"&{ks[0]} {ks[1]}\n"
                    cc += write_section(v,layer+1)
                    cc += "  "*layer + f"&END {ks[0]}\n"
                else:
                    cc += "  "*layer + f"&{k}\n"
                    cc += write_section(v,layer+1)
                    cc += "  "*layer + f"&END {k}\n"
                if layer < 2:
                    cc += "\n"
            elif isinstance(v,list):
                cc += "  "*layer + f"{k} "+ " ".join([str(i) for i in v]) + "\n"
            elif isinstance(v,(str,int,float)):
                cc += "  "*layer + f"{k} {v}\n"
            else:
                print(f"Warning: the value of {k} is {v}, which type ({type(v)}) is not supported now, just transfer {v} to a string.")
                cc += "  "*layer + f"{k} {v}\n"
        return cc
    cc = write_section(param)
    if save_path:
        print(f"Write CP2K input to {save_path}")
        with open(save_path,"w") as f:
            f.write(cc)
    return cc

def kpt2cp2k(kpt,kpt_model):
    if kpt_model in ["gamma","mp"]:
        return {"FORCE_EVAL": {"DFT": {"KPOINTS": {"SCHEME": "MONKHORST-PACK " + " ".join([str(i) for i in kpt[:3]])}}}}
    elif kpt_model == "direct":
        return {"FORCE_EVAL": {"DFT": {"KPOINTS": {"KPOINT": kpt}}}}
    else:
        print("Warning: do not support KPT model {model} now! Do not set the kpoints.")
        return {}
    

def stru2cp2k(stru:AbacusStru):
    '''
    Convert abacus structure to cp2k input.
    
    Args:
        stru: the abacus structure
    '''
    cell = stru.get_cell(bohr=False)
    coord = stru.get_coord(bohr=False,direct=False)
    label = stru.get_label(total=True)
    element = stru.get_element(number=False,total=True)
    mag = stru.get_mag()
    
    tmp = {
        "FORCE_EVAL": {
            "SUBSYS": {
                "CELL": {
                    "A": cell[0],
                    "B": cell[1],
                    "C": cell[2],
                },
                "COORD": {
                    label[i] : coord[i] for i in range(len(label))
                },
            },
            "DFT": {
                "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
                "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
            }
        }
    }
    
    label_uniq = []
    for idx,ilabel in enumerate(label):
        if ilabel not in label_uniq:
            label_uniq.append(ilabel)
            tmp["FORCE_EVAL"]["SUBSYS"][f"KIND&{ilabel}"] = {
                "ELEMENT": element[idx],
                "MAGNETIZATION": mag[idx],
                "BASIS_SET": "DZVP-MOLOPT-GTH" if element[idx] not in cp2k.BASIS else cp2k.BASIS[element[idx]],
                "POTENTIAL": "GTH-BLYP" if  element[idx] not in cp2k.POTENTIAL else cp2k.POTENTIAL[element[idx]],
            }
    return tmp

def calculation2cp2k(calculation, force, stress):
    if calculation in ["scf","md","relax","cell-relax"]:
        tmp = {"GLOBAL": {"RUN_TYPE": {"scf": "ENERGY", "md": "MD", "relax": "GEO_OPT","cell-relax":"CELL_OPT"}[calculation]}}
    else:
        print(f"Warning: the value of calculation is {calculation}, which set GLOBAL/RUN_TYPE to ENERGY")
        tmp = {"GLOBAL": {"RUN_TYPE": "ENERGY"}}
    
    if calculation == "scf" and force:
        tmp["GLOBAL"]["RUN_TYPE"] = "ENERGY_FORCE"
    
    if stress or calculation == "cell-relax":
        tmp["FORCE_EVAL"] = {"STRESS_TENSOR": "ANALYTICAL"}
    
    print_section = {}
    if calculation in ["relax", "cell-realx"] or force:
        print_section["FORCES&ON"] = {}
    if calculation in ["cell-relax"] or stress:
        print_section["STRESS_TENSOR&ON"] = {}
    if print_section:
        if "FORCE_EVAL" not in tmp:
            tmp["FORCE_EVAL"] = {}
        tmp["FORCE_EVAL"]["PRINT"] = print_section
    return tmp

def esolver2cp2k(esolver,dft_functional):
    tmp = {}
    if esolver == "ksdft":
        tmp = {"FORCE_EVAL": {"METHOD": "QS"}}
    else:
        print(f"Warning: the value of esolver is {esolver}, which is not supported now. Set FORCE_EVAL/METHOD to QS.")
        tmp = {"FORCE_EVAL": {"METHOD": "QS"}}
    
    dft_func = dft_functional.upper()
    tmp["FORCE_EVAL"]["DFT"] = {"XC": {f"XC_FUNCTIONAL&{dft_func}": {}}}
    
    return tmp

def kssolver2cp2k(kssolver):
    if kssolver == None:
        return {}
    tmp = {}
    if kssolver.lower() in ["genelpa","scalapack_gvx","cg", "dav"]:
        tmp = {"FORCE_EVAL": {"DFT": {"SCF": {"DIAGONALIZATION":{"ALGORITHM": "STANDARD"}}}}}
    else:
        print(f"Warning: the value of ks_solver is {kssolver}, which is not supported now. Set FORCE_EVAL/DFT/SCF/ALGORITHM to STANDARD.")
        tmp = {"FORCE_EVAL": {"DFT": {"SCF": {"DIAGONALIZATION":{"ALGORITHM": "STANDARD"}}}}}
    return tmp

def mix2cp2k(mixing_type,mixing_beta, mixing_ndim):
    if mixing_type == None:
        return {}
    
    tmp = {"FORCE_EVAL": {"DFT": {"SCF": {"MIXING": {}}}}}
    if mixing_type == "plain":
        tmp["FORCE_EVAL"]["DFT"]["SCF"]["MIXING"]["METHOD"] = "DIRECT_P_MIXING"
    elif mixing_type == "pulay":
        tmp["FORCE_EVAL"]["DFT"]["SCF"]["MIXING"]["METHOD"] = "PULAY_MIXING"
    elif mixing_type == "broyden":
        tmp["FORCE_EVAL"]["DFT"]["SCF"]["MIXING"]["METHOD"] = "BROYDEN_MIXING"
    else:
        print(f"Warning: the value of mixing_type is {mixing_type}, which is not supported now. Set FORCE_EVAL/DFT/SCF/MIXING/METHOD to DIRECT_P_MIXING.")
    
    if mixing_beta != None:
        tmp["FORCE_EVAL"]["DFT"]["SCF"]["MIXING"]["ALPHA"] = mixing_beta
    
    if mixing_ndim != None:
        tmp["FORCE_EVAL"]["DFT"]["SCF"]["MIXING"]["NBUFFER"] = mixing_ndim
    
    return tmp

def smearing2cp2k(smearing_method, smearing_sigma):
    if smearing_method == None:
        return {}
    
    if smearing_method == "fd":
        tmp = {"FORCE_EVAL": {"DFT": {"SCF": {"SMEAR&T": {"METHOD": "FERMI_DIRAC"}}}}}
    else:
        print(f"Warning: the value of smearing_method is {smearing_method}, which is not supported now. Will not set smear.")
        return {}
    
    return tmp

def scf2cp2k(scf_thr, scf_nmax, basis): 
    '''
    As the scf_thr in ABACUS is for difference of charge, while in CP2K is for energy, so based on a general sense,
    if basis is PW, then set the EPS_SCF to scf_thr*1e3,
    else, basis is lcao, then set the EPS_SCF to scf_thr*1e2
    '''
    if scf_thr == None:
        return {}
    
    if basis == None or basis == "pw":
        tmp = {"FORCE_EVAL": {"DFT": {"SCF": {"EPS_SCF": scf_thr*1e3}}}}   
    elif basis in ["lcao","lcao_in_pw"]:
        tmp = {"FORCE_EVAL": {"DFT": {"SCF": {"EPS_SCF": scf_thr*1e2}}}}
    else: 
        print(f"Warning: the value of basis is {basis}, which is not supported now. Will not set EPS_SCF.")
        return {}
    if scf_nmax != None:
        tmp["FORCE_EVAL"]["DFT"]["SCF"]["MAX_SCF"] = scf_nmax
    return tmp

def relaxthr2cp2k(calculation,force_thr, force_thr_ev,stress_thr):
    # MAX_FORCE bohr^-1*hartree
    if calculation not in ["relax","cell-relax"]:
        return {}
    
    if calculation == "relax":
        sub_section = "GEO_OPT"
    else:
        sub_section = "CELL_OPT"
        
    force = None
    if force_thr_ev != None:
        force = force_thr_ev * 0.03889 / 2
    if force_thr:
        force = force_thr / 2
    
    thr_set = {}
    if force != None:
        thr_set["MAX_FORCE"] = force
        thr_set["RMS_FORCE"] = force*2/3
    if stress_thr != None:
        thr_set["PRESSURE_TOLERANCE"] = stress_thr*1000
    
    if thr_set:
        return {"MOTION": {sub_section: thr_set}}
    else:
        return {}
    
def plusu2cp2k(plusu=None):
    return {}

def vdw2cp2k(vdw=None):
    if vdw == None:
        return {}   
            
def spin2cp2k(nspin):
    if nspin == None or nspin == 1:
        return {}
    elif nspin ==2:
        return {"FORCE_EVAL": {"DFT": {"UKS": "T"}}}
    else:
        print(f"Warning: the value of nspin is {nspin}, which is not supported now. Will not set UKS.")
        return {}

def input2cp2k(input_param:dict):
    
    # one by one parameter dict:
    oneone = {
        "suffix": lambda x: {"GLOBAL": {"PROJECT": x}},
    }
    
    tmp = {"GLOBAL": {"PRINT_LEVEL": "MEDIUM"}}
    
    # transfer oneone
    allkeys = list(input_param.keys())
    for k in allkeys:
        k = k.lower()
        if k in oneone:
            itmp = oneone[k](input_param[k])
            if itmp is not None:
                update_dict(tmp,itmp)
            del input_param[k]
    
    # transfer calculation, force, stress
    calculation = input_param.pop("calculation","scf")
    update_dict(tmp,calculation2cp2k(calculation,
                                        input_param.pop("cal_force",False),
                                        input_param.pop("cal_stress",False)))
    update_dict(tmp,scf2cp2k(input_param.pop("scf_thr",None), input_param.pop("scf_nmax",None), input_param.pop("basis_type",None)))
    update_dict(tmp,smearing2cp2k(input_param.pop("smearing_method",None), input_param.pop("smearing_sigma",None)))
    update_dict(tmp,mix2cp2k(input_param.pop("mixing_type",None), input_param.pop("mixing_beta",None), input_param.pop("mixing_ndim",None)))
    update_dict(tmp,kssolver2cp2k(input_param.pop("ks_solver",None)))
    update_dict(tmp,esolver2cp2k(input_param.pop("esolver","ksdft"), input_param.pop("dft_functional","pbe")))
    update_dict(tmp,plusu2cp2k(input_param.pop("plusu",None)))
    update_dict(tmp,vdw2cp2k(input_param.pop("vdw",None)))
    update_dict(tmp,spin2cp2k(input_param.pop("nspin",None)))
    update_dict(tmp,relaxthr2cp2k(calculation,input_param.pop("force_thr",None), input_param.pop("force_thr_ev",None), input_param.pop("stress_thr",None)))
    
    not_supportted_keys = list(input_param.keys())
    if len(not_supportted_keys) > 0:
        print(f"Warning: below keys are not supported now, and will ignore them:\n    ", ", ".join(not_supportted_keys))
    
    return tmp
    
def Abacus2Cp2k(abacus_path:str, save_path:str=None, cp2k_setting={}):
    '''
    Convert abacus input to cp2k input.
    
    1. read the abacus input file, and STRU/KPT
    2. generate the cp2k input file
    
    Args:
        abacus_path: the path of the abacus input file
        save_path: the path to save the cp2k input file, default is the same as the abacus input file
        cp2k_setting: the extra setting for cp2k input file
    '''
    cp2k_setting = copy.deepcopy(cp2k_setting)
    tmp = {}
    abacus_input = ReadInput(os.path.join(abacus_path,"INPUT"))
    abacus_stru = AbacusStru.ReadStru(os.path.join(abacus_path,"STRU"))
    abacus_kpt, kpoint_model = ReadKpt(os.path.join(abacus_path,"KPT"))
    
    if abacus_input == None or abacus_stru == None:
        print("Error: the abacus input file or stru file is not found.")
        return
    
    update_dict(tmp,stru2cp2k(abacus_stru))
    update_dict(tmp,kpt2cp2k(abacus_kpt,kpoint_model))
    update_dict(tmp,input2cp2k(abacus_input))
    update_dict(tmp,cp2k_setting)
    
    if save_path == None:
        save_path = os.path.join(abacus_path,"cp2k.inp")
    
    param2cp2k(tmp,save_path)
    
    
    