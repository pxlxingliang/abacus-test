import os,sys
from . import abacus as MyAbacus
from abacustest import constant
from . import comm
import copy

def ParamAbacus2Vasp(abacus_input):
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
    vasp_input = {}
    vasp_input["SYSTEM"] = abacus_input.pop("suffix","VASP")
    vasp_input["LREAL"] = "Auto"
    vasp_input["PREC"] = "Normal"
    vasp_input["ALGO"] = "Normal"
    
    # job type
    calculation = abacus_input.pop("calculation","scf")
    force = abacus_input.pop("force",False)
    stress = abacus_input.pop("stress",False)

    if "ecutwfc" in abacus_input:
        vasp_input["ENCUT"] = abacus_input.pop("ecutwfc") * constant.RY2EV
    
    basis_type = abacus_input.pop("basis_type","pw")
    if basis_type == "pw":
        vasp_input["EDIFF"] = abacus_input.pop("scf_thr",1e-9)/1E-2
    else:
        vasp_input["EDIFF"] = abacus_input.pop("scf_thr",1e-7) / 1E-1
    
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
            fixed_axes = abacus_input.pop("fixed_axes","None")
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
        vasp_input["NELM"] = abacus_input.pop("scf_nmax",100)
    
    smearing_method = abacus_input.pop("smearing_method",None)    
    if smearing_method != None:
        if smearing_method.startswith("gau"):
            vasp_input["ISMEAR"] = 0
        elif smearing_method == "mp":
            vasp_input["ISMEAR"] = 1
        elif smearing_method == "mp2":
            vasp_input["ISMEAR"] = 2
        elif smearing_method == "fd":
            vasp_input["ISMEAR"] = -1
        
        vasp_input["SIGMA"] = abacus_input.pop("smearing_sigma",0.015)*constant.RY2EV

    if "symmetry" in abacus_input:
        symmetry = abacus_input.pop("symmetry")
        if symmetry == 0:
            vasp_input["ISYM"] = 0
        elif symmetry == 1:
            vasp_input["ISYM"] = 1
        elif symmetry == -1:
            vasp_input["ISYM"] = -1
        else:
            raise ValueError(f"symmetry {abacus_input['symmetry']} is not supported.")
    
    if calculation in ["relax","cell-relax"]:
        if "force_thr_ev" in abacus_input: 
            vasp_input["EDIFFG"] = abacus_input.pop("force_thr_ev")*-1
        if "force_thr" in abacus_input:
            vasp_input["EDIFFG"] = abacus_input.pop("force_thr")*-1*constant.RY2EV
           
    if "kspacing" in abacus_input:
        kspacing = abacus_input.pop("kspacing")
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
        vasp_input["NSW"] = abacus_input.pop("relax_nmax")
    
    if calculation == "md" and "md_nstep" in abacus_input:
        vasp_input["NSW"] = abacus_input.pop("md_nstep")   
    
    if "nspin" in abacus_input:
        vasp_input["ISPIN"] = abacus_input.pop("nspin")
        if vasp_input["ISPIN"] == 4: vasp_input["ISPIN"] = 2
        if vasp_input["ISPIN"] in [2,4]:
            vasp_input["LORBIT"] = 11 # OUTPUT THE MAGNETIC MOMENT of each atom
    
    if "dft_plus_u" in abacus_input:
        if comm.IsTrue(abacus_input.pop("dft_plus_u")):
            vasp_input["LDAU"] = ".TRUE."
            vasp_input["LDAUTYPE"] = 2
            if "orbital_corr" in abacus_input:
                vasp_input["LDAUL"] = abacus_input.pop("orbital_corr")
                if "2" in vasp_input["LDAUL"]:
                    vasp_input["LMAXMIX"] = 4
                if "3" in vasp_input["LDAUL"]:
                    vasp_input["LMAXMIX"] = 6
            if "hubbard_u" in abacus_input:
                vasp_input["LDAUU"] = abacus_input.pop("hubbard_u")
                vasp_input["LDAUJ"] = " ".join(["0" for i in vasp_input["LDAUU"].split()])
        else:
            vasp_input["LDAU"] = ".FALSE."
            abacus_input.pop("orbital_corr",None)
            abacus_input.pop("hubbard_u",None)
    
    if "nupdown" in abacus_input:
        vasp_input["NUPDOWN"] = abacus_input.pop("nupdown")
    
    # lspinorb and nonclinear
    if "lspinorb" in abacus_input:
        if comm.IsTrue(abacus_input.pop("lspinorb")):
            vasp_input["LSORBIT"] = ".TRUE."
            vasp_input["ISYM"] = -1
        else:
            vasp_input["LSORBIT"] = ".FALSE."
    if "noncolin" in abacus_input:
        if comm.IsTrue(abacus_input.pop("noncolin")):
            vasp_input["LNONCOLLINEAR"] = ".TRUE."
        else:
            vasp_input["LNONCOLLINEAR"] = ".FALSE."
            
    # delta spin
    if "sc_mag_switch" in abacus_input:
        if comm.IsTrue(abacus_input.pop("sc_mag_switch")):
            vasp_input["SCTYPE"] = 1
            if "nsc" in abacus_input:
                vasp_input["NSC"] = abacus_input.pop("nsc")
            if "nsc_min" in abacus_input:
                vasp_input["NSCMIN"] = abacus_input.pop("nsc_min")
            if "sc_thr" in abacus_input:
                vasp_input["SCDIFF"] = abacus_input.pop("sc_thr")
            if "alpha_trial" in abacus_input:
                vasp_input["INISC"] = abacus_input.pop("alpha_trial")
            if "sccut" in abacus_input:
                vasp_input["SCCUT"] = abacus_input.pop("sccut")
        else:
            vasp_input["SCTYPE"] = 0
            abacus_input.pop("nsc",None)
            abacus_input.pop("nsc_min",None)
            abacus_input.pop("sc_thr",None)
            abacus_input.pop("alpha_trial",None)
            abacus_input.pop("sccut",None)
    
    if len(abacus_input) > 0:
        print("WARNING: The following parameters are not converted to VASP:")
        for ip in abacus_input.keys():
            if ip in ["kspacing"]: continue
            print("         %s" % ip)
    
    return vasp_input

def GenIncarFromStru(stru,vasp_setting):
    # if nspin != 1, then set the initial magmom
    if vasp_setting.get("ISPIN",1) != 1:
        mag = stru.get_atommag()  # mag is a list of list, eahc list is the magmom of one atom, for noncolinear, it should be [magx,magy,magz]
        labels = stru.get_label(total=True)
        label = []
        for i in labels:
            if i not in label:
                label.append(i)
        magmom = []
        if vasp_setting.get("LNONCOLLINEAR",False):
            # should create the magmom, a lisf of 3*N, N is the number of atoms
            for imag in mag:
                if isinstance(imag,float):
                    magmom += [0,0,imag]
                elif isinstance(imag,list):
                    if len(imag) == 1:
                        magmom += [0,0,imag[0]]
                    elif len(imag) == 3:
                        magmom += imag
                    else:
                        print("ERROR: The magmom should be a list of 1 or 3!!!")
                        sys.exit(1)
                else:
                    print("ERROR: The magmom should be a list of float or list!!!")
                    sys.exit(1)
        else:
            for imag in mag:
                if isinstance(imag,float):
                    magmom.append(imag)
                elif isinstance(imag,list):
                    if len(imag) == 1:
                        magmom.append(imag[0])
                    elif len(imag) == 3:
                        magmom.append((imag[0]**2 + imag[1]**2 + imag[2]**2)**0.5)
                    else:
                        print("ERROR: The magmom should be a list of 1 or 3!!!")
                        sys.exit(1)
                else:
                    print("ERROR: The magmom should be a list of float!!!")
                    sys.exit(1)
        
        vasp_setting["MAGMOM"] = VaspList2String(magmom)
        
        if comm.IsTrue(vasp_setting.get("SCTYPE")):
            # if sc_mag_switch is true, then set the initial magmom
            vasp_setting["M_CONSTR"] = vasp_setting["MAGMOM"]
            constrain = stru.get_constrain()
            lambda_ = stru.get_lambda()
            c = []
            l = []
            if not vasp_setting.get("LNONCOLLINEAR",False):
                for ic in constrain:
                    if isinstance(ic,list):
                        if False in ic:
                            c.append(0)
                    else:
                        c.append(1 if ic else 0)
                for il in lambda_:
                    if isinstance(il,list):
                        l += il
                    else:
                        l.append(il)
            else:
                for ic in constrain:
                    if isinstance(ic,bool):
                        c += [1 if ic else 0]*3
                    else:
                        c += [1 if i else 0 for i in ic]
                for il in lambda_:
                    if isinstance(il,(float,int)):
                        l += [il]*3
                    elif isinstance(il,list):
                        l += il
                    else:
                        print("ERROR: The lambda should be a list of float or list!!!")
                        sys.exit(1)

            vasp_setting["CONSTRL"] = VaspList2String(c)
            vasp_setting["LAMBDA"] = VaspList2String(l)

def VaspList2String(list1):
    # convert the list1 to string
    # need to merge the same value like: [1, 1, 1, 0, 0, 0, 0] to "3*1 4*0"
    #list1 = [str(i) for i in list1]
    if len(list1) == 0:
        return ""
    
    c = ""
    v = list1[0]
    n = 1
    for i in list1[1:]:
        if i == v or (isinstance(i,(float,int)) and isinstance(v,(float,int)) and abs(i-v) < 1E-6):
            n += 1
        else:
            if n == 1:
                c += f"{v} "
            else:
                c += f"{n}*{v} "
            v = i
            n = 1
    if n == 1:
        c += f"{v} "
    else:
        c += f"{n}*{v} "
    return c
    
    

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
    vasp_setting = copy.deepcopy(vasp_setting)
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
    if "LWAVE" not in vasp_input: vasp_input["LWAVE"] = ".FALSE."
    if "LCHARG" not in vasp_input: vasp_input["LCHARG"] = ".FALSE."
    
    # 4.1 special setting in vasp_setting
    emax_coef = vasp_setting.pop("emax_coef",None)
    
    # 4.2 set ENCUT
    if emax != None and emax_coef != None:
        vasp_input["ENCUT"] = emax * emax_coef
    
    # 4.3 set MAGMOM and deltaspin    
    GenIncarFromStru(stru,vasp_input)
        
    # 4.4 update the vasp_input by vasp_setting    
    if vasp_setting:
        for ik,iv in vasp_setting.items():
            vasp_input[ik.upper()] = iv
    if save_path != None:
        with open(os.path.join(save_path,"INCAR"),"w") as f:
            for key,value in vasp_input.items():
                f.write(f"{key} = {value}\n")
            
    return {
        "incar": vasp_input,  # a dict of incar
        "kpoints": kpoints,     # string of kpoints
        "poscar": poscar        # string of poscar
    }
