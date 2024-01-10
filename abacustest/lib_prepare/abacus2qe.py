from typing import Dict
import os
import numpy as np
from . import abacus as MyAbacus
from . import qe as MyQe
from . import comm


def ParamAbacus2Qe(param:Dict[str,any]):
    '''transfer the abacus input to qe input
    
    Parameters
    ----------
    param : Dict[str,any]
        the abacus input parameters
        
    Returns
    -------
    dict of qe input parameters
    
    
    '''
    print("Transfer abacus input to qe input")
    qp = {"control":{},
          "system":{},
          "electrons":{},
          "ions":{},
          "cell":{}}
    
    if param.pop("basis_type","pw") != "pw":
        print("WARNING: The basis_type is not pw, is %s" % param.get("basis_type","pw"))
    
    calculation = param.pop("calculation","scf").lower()
    if calculation == "scf":
        qp["control"]["calculation"] = "scf"
    elif calculation == "relax":
        qp["control"]["calculation"] = "relax"
    elif calculation == "cell-relax":
        qp["control"]["calculation"] = "vc-relax"
    else:
        print("ERROR: calculation %s is not supported now!" % calculation)
        return None
    
    if qp["control"]["calculation"] in ["relax","cell-relax"]:
        qp["control"]["forc_conv_thr"] = 0.001
        if param.get("force_thr_ev"):
            qp["control"]["forc_conv_thr"] = float(param.get("force_thr_ev")) / 25.7112
            param.pop("force_thr_ev")
        if param.get("force_thr"):
            qp["control"]["forc_conv_thr"] = param.get("force_thr")
            param.pop("force_thr")
    if qp["control"]["calculation"] == "vc-relax":
        qp["calculation"]["press_conv_thr"] = param.pop("press_conv_thr",0.5)
    
    # set common parameters
    qp["control"]["pseudo_dir"] = param.pop("pseudo_dir",".")
    qp["control"]["tstress"] = ".true." if bool(int(param.pop("cal_stress",0))) else ".false."
    qp["control"]["tprnfor"] = ".true." if bool(int(param.pop("cal_force",0))) else ".false."
    qp["system"]["nosym"] = ".false." if bool(int(param.pop("symmetry",0))) else ".true."
    # these parameters will use the same value.
    for para_aba,para_qe in [["scf_nmax",("electrons","electron_maxstep")],
                             ["ecutwfc", ("system","ecutwfc")],
                             ["scf_thr", ("electrons","conv_thr")],
                             ["nspin", ("system","nspin")],
                             ]:
        if para_aba in param:
            qp[para_qe[0]][para_qe[1]] = param.pop(para_aba)
    
    # set smearing  
    smearing = param.pop("smearing_method","gauss").lower()
    if smearing == "fixed":
        qp["system"]["occupations"] = "fixed"
    elif smearing.startswith("gauss"):
        qp["system"]["occupations"] = "smearing"
        qp["system"]["smearing"] = "gauss"
        qp["system"]["degauss"] = param.pop("smearing_sigma","0.015")
    else:
        print("WARNING: smearing_method %s is not supported now! Will not set occupations in QE." % smearing)
    
    # set ks_solver
    kssolver = param.pop("ks_solver","cg").lower()
    if kssolver == "cg":
        qp["electrons"]["diagonalization"] = "cg"
    elif kssolver == "dav":
        qp["electrons"]["diagonalization"] = "david"
    else:
        print("WARNING: ks_solver %s is not supported now, will not set diagonalization is QE." % kssolver)
    
    # set mixing_type
    mixing = param.pop("mixing_type","broyden").lower()
    if mixing == "broyden":
        qp["electrons"]["mixing_mode"] = "plain"
        qp["electrons"]["mixing_beta"] = param.pop("mixing_beta","0.7")
    else:
        print("WARNING: mixing_type %s is not supported now, will not set mixing_mode in QE." % mixing)
        
    # these paramters will be ignored    
    for ip in ["suffix","basis_type","gamma_only","kpt_file"]:
        param.pop(ip,None)
    
    # print out warning for the left parameters
    left_param = []
    for ip in param:
        if ip.lower().startswith("out_"):
            continue
        left_param.append(ip)
    if len(left_param) > 0:    
        print("WARNING: The following parameters are not supported in QE, will be ignored:")
        for ip in left_param:
            print("         %s" % ip)
    return qp 

def Abacus2Qe(path: str= ".", save_path: str = None):
    '''transfer the abacus input to qe input
    
    Parameters
    ----------
    path : str, optional
        the path of abacus input, by default ".".
    save_path : str, optional
        the path of qe input, by default "input" in path.
        
    Returns
    -------
    QeInputs/None
    
    
    '''
    print("Transfer abacus input to qe input in %s" % path)
    if not os.path.isfile(os.path.join(path,"INPUT")):
        print("ERROR: Not find INPUT file in %s" % path)
        return None
    
    abacus_input = MyAbacus.ReadInput(os.path.join(path,"INPUT"))
    stru_file = os.path.join(path,abacus_input.get("stru_file","STRU"))
    if not os.path.isfile(stru_file):
        print("ERROR: Not find STRU file in %s" % path)
        return None
    
    abacus_stru = MyAbacus.AbacusStru.ReadStru(stru_file)
    if abacus_stru == None:
        print("ERROR: Read STRU file failed in %s" % path)
        return None
    
    param = ParamAbacus2Qe(abacus_input)
    if not param:
        print("ERROR: Transfer abacus input to qe input failed in %s" % path)
        return None
    
    # analysis the abacus_stru
    labels = abacus_stru.get_label()
    label = []
    for i in labels:
        if i not in label:
            label.append(i)
    atom_number = [labels.count(i) for i in label]
    stru = abacus_stru.get_stru()
    cell = stru["cell"]
    coord = stru["coord"]
    cartesian = stru["cartesian"]
    
    mag = abacus_stru.get_mag()
    pp = abacus_stru.get_pp()
    
    # set structure related param
    param["system"]["nat"] = len(coord)
    param["system"]["ntyp"] = len(label)
    param["system"]["ibrav"] = 0
    param["system"]["celldm(1)"] = stru["lat"]
    if param["system"].get("nspin") == 2:
        for i in range(len(mag)):
            if mag[i] != 0:
                param["system"]["starting_magnetization(%d)" % (i+1)] = mag[i]
    
    # analysis the K_POINTS
    # if set basis_type to lcao and set gamma_only, then set kpt_type to gamma
    if abacus_input.get("basis_type","pw") == "lcao" and int(abacus_input.get("gamma_only",0)):
        kpt = [1,1,1,0,0,0]
    elif abacus_input.get("kspacing",0):
        kspacing = abacus_input.get("kspacing",0)
        if isinstance(kspacing,str):
            try:
                kspacing = [float(i) for i in kspacing.split()]
                assert(len(kspacing) in [1,3])
            except:
                print("ERROR: kspacing should be a float or three floats")
                return None
            if len(kspacing) == 1:
                kspacing = [kspacing[0]]*3
        elif isinstance(kspacing,float):
            kspacing = [kspacing]*3
        else:
            print("ERROR: kspacing should be a float or three floats")
            return None
        kpt = comm.kspacing2kpt(kspacing,[[i*stru["lat"] for i in j] for j in cell])
        kpt.extend([0,0,0])
    else:
        kptf = abacus_input.get("kpt_file","KPT")
        if not os.path.isfile(os.path.join(path,kptf)):
            print("ERROR: Not find KPT file in %s" % path)
            return None
        with open(os.path.join(path,kptf)) as f1: lines = f1.readlines()
        if lines[2].strip().lower() == "gamma":
            iline = lines[3].split()
            kpt = [int(i) for i in iline[:3]] + [float(i) for i in iline[3:6]]
        else:
            print("ERROR: Not support KPT type: %s now!!!" % lines[2].strip())
            return None

    qe = MyQe.QeInputs(param=param,
                    label=label,
                    atom_number=atom_number,
                    cell=cell,
                    coord=coord,
                    pp=pp,
                    kpt = kpt,
                    kpt_type ="automatic",
                    mass=abacus_stru.get_mass(),
                    cell_type = "alat",
                    
                    coord_type = "alat" if cartesian else "crystal")
    
    save_path = save_path if save_path != None else os.path.join(path,"input")
    qe.write(save_path)
    return qe