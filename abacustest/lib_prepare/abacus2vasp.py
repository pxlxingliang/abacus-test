import os,sys
from . import abacus as MyAbacus
from abacustest import constant
from . import comm

def ParamAbacus2Vasp(abacus_input,incar= None):
    '''
    Generate the vasp input from abacus input.
    
    As the parameters are not totally related, so we just convert the parameters we need.
    - calculation
        - SCF
        - RELAX
        - CELL_RELAX
        - NSCF
    - smearing
    - symmetry
    - kspacing
    
    Below cases will not be considered:
    fixed_axes/
    
    for scf_thr:
    pw: decrease 1e-2
        < 1e-8 -> 1e-6
    lcao: decrease 1e-1
        < 1e-7 -> 1e-6
    '''
    vasp_input = {
        "LREAL": "Auto",
        "PREC":  "Normal",
        "ALGO":  "Normal" 
    }
    # job type
    calculation = abacus_input.get("calculation","scf")
    force = abacus_input.get("force",False)
    stress = abacus_input.get("stress",False)
    
    basis_type = abacus_input.get("basis_type","pw")
    if basis_type == "pw":
        vasp_input["EDIFF"] = abacus_input.get("scf_thr",1e-9)/1E-2
    else:
        vasp_input["EDIFF"] = abacus_input.get("scf_thr",1e-7) / 1E-1
    
    if calculation == "scf":
        vasp_input["IBRION"] = -1
        if stress:
            vasp_input["ISIF"] = 2
        else:
            vasp_input["ISIF"] = 0
    elif calculation in ["relax","cell-relax"]:
        vasp_input["IBRION"] = 2
        if calculation == "relax":
            vasp_input["ISIF"] = 2
        else:
            fixed_axes = abacus_input.get("fixed_axes","None")
            if fixed_axes == "None":
                vasp_input["ISIF"] = 3
            elif fixed_axes == "volume":
                vasp_input["ISIF"] = 4
            elif fixed_axes == "shape":
                vasp_input["ISIF"] = 7
            else:
                vasp_input["ISIF"] = 3  # FOR other cases, just set to 3
    elif calculation == "md":
        vasp_input["IBRION"] = 0
    elif calculation == "nscf":
        vasp_input["IBRION"] = -1
    else:
        raise ValueError(f"calculation {calculation} is not supported.")
    
    if "scf_nmax" in abacus_input:
        vasp_input["NELM"] = abacus_input.get("scf_nmax",100)
        
    if "smearing_method" in abacus_input and abacus_input["smearing_method"] != "fixed":
        if abacus_input["smearing_method"] == "gaussian":
            vasp_input["ISMEAR"] = 0
        elif abacus_input["smearing_method"] == "mp":
            vasp_input["ISMEAR"] = 2
        elif abacus_input["smearing_method"] == "fd":
            vasp_input["ISMEAR"] = -1
        
        vasp_input["SIGMA"] = abacus_input.get("smearing_sigma",0.015)*constant.RY2EV

    if "symmetry" in abacus_input:
        if abacus_input["symmetry"] == 0:
            vasp_input["ISYM"] = 0
        elif abacus_input["symmetry"] == 1:
            vasp_input["ISYM"] = 1
        elif abacus_input["symmetry"] == -1:
            vasp_input["ISYM"] = -1
        else:
            raise ValueError(f"symmetry {abacus_input['symmetry']} is not supported.")
    
    if calculation in ["relax","cell-relax"]:
        if "force_thr_ev" in abacus_input: 
            vasp_input["EDIFFG"] = abacus_input["force_thr_ev"]*-1
        if "force_thr" in abacus_input:
            vasp_input["EDIFFG"] = abacus_input["force_thr"]*-1*constant.RY2EV
           
    if "kspacing" in abacus_input:
        kspacing = abacus_input["kspacing"]
        if isinstance(kspacing,list):
            kspacing = sum(kspacing)/3
        elif isinstance(kspacing,str):
            kspacing = kspacing.split()
            if len(kspacing) == 1:
                kspacing = float(kspacing[0])
            elif len(kspacing) == 3:
                kspacing = sum([float(i) for i in kspacing])/3
            else:
                raise ValueError(f"kspacing {kspacing} is not supported.")
        vasp_input["KSPACING"] = kspacing / constant.BOHR2A
        
    if calculation in ["relax","cell-relax"] and "relax_nmax" in abacus_input:
        vasp_input["NSW"] = abacus_input["relax_nmax"]
    
    if calculation == "md" and "md_nstep" in abacus_input:
        vasp_input["NSW"] = abacus_input["md_nstep"]   
    
    if "nspin" in abacus_input:
        vasp_input["ISPIN"] = abacus_input["nspin"]
    
    if "dft_plus_u" in abacus_input:
        if comm.IsTrue(abacus_input["dft_plus_u"]):
            vasp_input["LDAU"] = ".TRUE."
            vasp_input["LDAUTYPE"] = 2
            if "orbital_corr" in abacus_input:
                vasp_input["LDAUL"] = abacus_input["orbital_corr"]
            if "hubbard_u" in abacus_input:
                vasp_input["LDAUU"] = abacus_input["hubbard_u"]
                vasp_input["LDAUJ"] = " ".join(["0" for i in abacus_input["hubbard_u"].split()])
    
    if "nupdown" in abacus_input:
        vasp_input["NUPDOWN"] = abacus_input["nupdown"]
    
    # lspinorb and nonclinear
    if "lspinorb" in abacus_input and comm.IsTrue(abacus_input["lspinorb"]):
        vasp_input["LSORBIT"] = ".TRUE."
    if "noncolin" in abacus_input and comm.IsTrue(abacus_input["noncolin"]):
        vasp_input["LNONCOLLINEAR"] = ".TRUE."

    if incar != None:
        with open(incar,"w") as f:
            for key,value in vasp_input.items():
                f.write(f"{key} = {value}\n")
    return vasp_input
    

def KptAbacus2Vasp(abacus_kpt,vasp_kpt="KPOINTS"):
    """
    Convert the kpt from abacus to vasp.
    
    Parameters
    ----------
    abacus_kpt : str
        the kpt file of abacus.
    vasp_kpt : str
        the file to save the vasp kpt.
        
    Only support the automatic kpt.
    """
    if not os.path.isfile(abacus_kpt):
        print(f"ERROR: Not find kpt file {abacus_kpt}!!!")
        sys.exit(1)
    with open(abacus_kpt) as f: aba_kpt = f.readlines()
    
    model = aba_kpt[2].strip().lower()
    
    if model in ["gamma","mp"]:
        kpt = [int(i) for i in aba_kpt[3].split()[:3]]
        shift_ = [float(i) for i in aba_kpt[3].split()[3:6]]

        # write the kpt to vasp_kpt
        cc = "K_POINTS\n0\nGamma\n"
        cc += " ".join([str(i) for i in kpt])+"\n"
        cc += " ".join([str(i) for i in shift_])
        if vasp_kpt != None:
            with open(vasp_kpt,"w") as f:
                f.write(cc)
        return cc       
    else:
        print(f"ERROR: Not support kpt model {model}!!!")
        return ""

def gen_potcar(potcar,label,save_path="POTCAR"):
    """
    Generate the potcar from potcar.
    
    Parameters
    ----------
    potcar : str|dict
        the path of potcar.
        If is a str, it should be the path of PBE-PAW POTCARS, and will use the potcar recommended by vaspwiki for each element.
            https://www.vasp.at/wiki/index.php/Available_PAW_potentials 
        If is a dict, it should be the dict of potcar, and the key is the label of atoms, and the value is the path of potcar.
    label : list
        the label of atoms.
    save_path : str
        the path to save the potcar.
    """
    if potcar == None:
        print("ERROR: Not find potcar!!!")
        return None
    print("potcar/label:",potcar,label)
    all_pots = []
    if isinstance(potcar,str):
        if not os.path.isdir(potcar):
            print(f"ERROR: Not find potcar folder {potcar}!!!")
            return None
        from .vasp import RECOMMEND_PAW
        for i in label:
            if i not in RECOMMEND_PAW:
                print(f"ERROR: Not find potcar for {i}!!!")
                return None
            ipot = os.path.join(potcar,RECOMMEND_PAW[i],"POTCAR")
            if not os.path.isfile(ipot):
                print(f"ERROR: Not find potcar for {i}: {ipot}!!!")
                return None
            all_pots.append(ipot)
    elif isinstance(potcar,dict):
        for i in label:
            if i not in potcar:
                print(f"ERROR: Not find potcar for {i}!!!")
                return None
            ipot = potcar[i]
            if not os.path.isfile(ipot):
                print(f"ERROR: Not find potcar for {i}: {ipot}!!!")
                return None
            all_pots.append(ipot)
    else:
        print("ERROR: potcar must be str or dict!!!")
        return None

    emax = []
    potcar_cc = ""
    for ipot in all_pots:
        with open(ipot) as g:
            lines = g.readlines()
            potcar_cc += "".join(lines)
            for iline in lines:
                if iline.startswith("   ENMAX  ="):
                    emax.append(float(iline.split(";")[0].split("=")[1]))
                    #   ENMAX  =  250.000; ENMIN  =  200.000 eV
                    break
    if save_path != None:
        with open(save_path,"w") as f:
            f.write(potcar_cc)
    emax = max(emax)
    return emax    
    
def Abacus2Vasp(abacus_path:str, save_path:str=None, potcar=None, vasp_setting={}):
    """
    Convert the abacus input to vasp input.
    
    Need convert four files:
    INPUT -> INCAR
    STRU -> POSCAR
    KPT -> KPOINTS
    POTCAR          # the prepare of this file based on potcar
    
    Parameters
    ----------
    abacus_path : str
        the path of abacus input.
    save_path : str
        the path to save vasp input.
    potcar : str|dict
        the path of potcar.
    emax_coef : float
        the coefficient to calculate the ENCUT. 
            ENCUT = emax * emax_coef (default is 1.5)
            emax is the max ENMAX in POTCAR.
    """
    abacus_input = MyAbacus.ReadInput(os.path.join(abacus_path,"INPUT"))
    stru = MyAbacus.AbacusStru.ReadStru(os.path.join(abacus_path,abacus_input.get("stru_file","STRU")))
    
    # 1. write KPT
    if abacus_input.get("basis_type","lcao") and abacus_input.get("gamma_only",False):
        kpoints = "K_POINTS\n0\nGamma\n1 1 1\n0 0 0"
        if save_path != None:
            with open(os.path.join(save_path,"KPOINTS"),"w") as f:
                f.write(kpoints)
    elif "kspacing" not in abacus_input:
        kpoints = KptAbacus2Vasp(os.path.join(abacus_path,abacus_input.get("kpt_file","KPT")),None if save_path==None else os.path.join(save_path,"KPOINTS"))
    else:
        kpoints = None

    # 2. write POSCAR
    poscar = stru.write2poscar(None if save_path==None else os.path.join(save_path,"POSCAR"))
    
    # 3. gen POTCAR
    emax = gen_potcar(potcar,stru.get_element(number=False,total=False), None if save_path==None else os.path.join(save_path,"POTCAR"))
    if emax == None:
        print("ERROR: gen POTCAR failed!!!")
    
    # 4. WRITE INCAR
    vasp_input = ParamAbacus2Vasp(abacus_input)
    
    # 4.1 special setting in vasp_setting
    emax_coef = vasp_setting.pop("emax_coef",1.5)
    
    # 4.2 set ENCUT
    if emax != None:
        vasp_input["ENCUT"] = emax * emax_coef
        
    # 4.3 if nspin != 1, then set the initial magmom
    if abacus_input.get("nspin",1) != 1:
        mag = stru.get_mag()
        labels = stru.get_label()
        label = []
        for i in labels:
            if i not in label:
                label.append(i)
        atom_number = [labels.count(i) for i in label]
        vasp_input["MAGMOM"] = " ".join([f"{atom_number[i]}*{mag[i]}" for i in range(len(mag))])
        
    # 4.4 update the vasp_input by vasp_setting    
    if vasp_setting:
        vasp_input.update(vasp_setting)
    if save_path != None:
        with open(os.path.join(save_path,"INCAR"),"w") as f:
            for key,value in vasp_input.items():
                f.write(f"{key} = {value}\n")
            
    return {
        "incar": vasp_input,  # a dict of incar
        "kpoints": kpoints,     # string of kpoints
        "poscar": poscar        # string of poscar
    }
