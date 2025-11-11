from typing import Dict
import os,copy
import numpy as np
from abacustest.constant import PERIOD_DICT_NUMBER

from abacustest.lib_collectdata.abacus.abacus import Abacus
from . import abacus as MyAbacus
from . import qe as MyQe
from . import comm

def set_mixing_type(param,qp):
    mixing = param.pop("mixing_type","").lower()
    if not mixing:
        pass
    elif mixing == "broyden":
        qp["electrons"]["mixing_mode"] = "plain"
        if "mixing_beta" in param:
            qp["electrons"]["mixing_beta"] = param.pop("mixing_beta")
    else:
        print("WARNING: mixing_type %s is not supported now, will not set mixing_mode in QE." % mixing)

def set_smearing(param,qp):
    smearing = param.pop("smearing_method","")
    if not smearing:
        pass
    elif smearing.lower() == "fixed":
        qp["system"]["occupations"] = "fixed"
    elif smearing.lower().startswith("gauss") or smearing.lower() in ["mp", "fd", "mv"]:
        qp["system"]["occupations"] = "smearing"
        qp["system"]["smearing"] = "gauss" if smearing.lower() == "gauss" else smearing.lower()
        if "smearing_sigma" in param:
            qp["system"]["degauss"] = param.pop("smearing_sigma")
        else:
            print("WARNING: smearing_sigma is not set in abacus input, will set degauss to 0.015 in QE.")
            qp["system"]["degauss"] = 0.015 
    else:
        print("WARNING: smearing_method %s is not supported now! Will not set occupations in QE." % smearing)

def set_ks_solver(param,qp):
    kssolver = param.pop("ks_solver","")
    if not kssolver:
        pass
    elif kssolver.lower() == "cg":
        qp["electrons"]["diagonalization"] = "cg"
    elif kssolver.lower() in ["dav","dav_subspace"]:
        qp["electrons"]["diagonalization"] = "david"
    else:
        print("WARNING: ks_solver %s is not supported now, will not set diagonalization is QE." % kssolver)

def set_relax_param(param,qp,calculation):
    # set relax method
    if calculation in ["relax", "cell-relax"] and "relax_method" in param:
        relax_method = param.pop("relax_method")
        if relax_method == "bfgs":
            qp["ions"]["ion_dynamics"] = "bfgs"
            if calculation == "cell-relax":
                qp["cell"]["cell_dynamics"] = "bfgs"
            if "relax_bfgs_w1" in param:
                qp["ions"]["w_1"] = param.pop("relax_bfgs_w1")
            if "relax_bfgs_w2" in param:
                qp["ions"]["w_2"] = param.pop("relax_bfgs_w2")
            if "relax_bfgs_rmax" in param:
                qp["ions"]["trust_radius_max"] = param.pop("relax_bfgs_rmax")
            if "relax_bfgs_rmin" in param:
                qp["ions"]["trust_radius_min"] = param.pop("relax_bfgs_rmin")
            if "relax_bfgs_init" in param:
                qp["ions"]["trust_radius_ini"] = param.pop("relax_bfgs_init")
        else:
            print("WARNING: relax_method %s is not supported now, will not set ion_dynamics in QE." % relax_method)

def set_spin_soc(param,qp):
    if "noncolin" in param:
        noncolin = bool(param.pop("noncolin",1))
        if noncolin:
            qp["system"]["noncolin"] = ".true."
            qp["system"]["nspin"] = 4
        else:
            qp["system"]["noncolin"] = ".false."   
    
    if "lspinorb" in param:
        if param.pop("lspinorb"):
            qp["system"]["lspinorb"] = ".true."
            qp["system"]["noinv"] = ".true."
            qp["system"]["noncolin"] = ".true."
            qp["system"]["nspin"] = 4
        else:
            qp["system"]["lspinorb"] = ".false."   

def set_fixed_axes(param,qp):
    if "fixed_axes" in param:
        fixed_axes = param.pop("fixed_axes")
        if fixed_axes == "None":
            qp["cell"]["cell_dofree"] = "all"
        elif fixed_axes == "volume":
            qp["cell"]["cell_dofree"] = "shape"
        elif fixed_axes == "shape":
            qp["cell"]["cell_dofree"] = "volume"
        elif fixed_axes in ["a","b","c"]:
            qp["cell"]["cell_dofree"] = f"fix{fixed_axes}"
        elif fixed_axes in ["ab","ac","bc"]:
            qp["cell"]["cell_dofree"] = {"a","b","c"}.difference(set(fixed_axes)).pop()
        else:
            print("WARNING: fixed_axes %s is not supported now, will not set cell_dofree in QE." % fixed_axes)

def set_cdft(param, qp):
    if "sc_mag_switch" in param:
        sc_mag_switch = param.pop("sc_mag_switch")
        sc_direction_only = param.pop("sc_direction_only", False)
        if sc_mag_switch:
            if sc_direction_only:
                qp["system"]["constrained_magnetization"] = "atomic direction"
            else:
                qp["system"]["constrained_magnetization"] = "atomic"

def ParamAbacus2Qe(input_param:Dict[str,any],version=7.0,qe_param={}):
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
    param = copy.deepcopy(input_param)
    qp = {"control":{},
          "system":{},
          "electrons":{},
          "ions":{},
          "cell":{}}
    basis_type = param.pop("basis_type","pw") 
    if basis_type != "pw":
        print("WARNING: The basis_type is not pw, is %s" % basis_type)
    
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
    
    if "relax_nmax" in param and calculation in ["relax","cell-relax","md"]:
        qp["control"]["nstep"] = param.pop("relax_nmax")

    if qp["control"]["calculation"] in ["relax","vc-relax"]:
        qp["control"]["forc_conv_thr"] = 0.001
        if param.get("force_thr_ev"):
            qp["control"]["forc_conv_thr"] = float(param.get("force_thr_ev")) / 25.7112
            param.pop("force_thr_ev")
        if param.get("force_thr"):
            qp["control"]["forc_conv_thr"] = param.get("force_thr")
            param.pop("force_thr")
    if qp["control"]["calculation"] == "vc-relax":
        qp["cell"]["press_conv_thr"] = param.pop("stress_thr",0.5)
    
    # set common parameters
    qp["control"]["pseudo_dir"] = param.pop("pseudo_dir",".")
    qp["control"]["tstress"] = ".true." if bool(int(param.pop("cal_stress",0))) else ".false."
    qp["control"]["tprnfor"] = ".true." if bool(int(param.pop("cal_force",0))) else ".false."
    isymm = param.pop("symmetry",None)
    if isymm is not None:
        isymm = int(isymm)
        if isymm > 0:
            qp["system"]["nosym"] = ".false."
        else:
            qp["system"]["nosym"] = ".true."
    # these parameters will use the same value.
    for para_aba,para_qe in [["scf_nmax",("electrons","electron_maxstep")],
                             ["ecutwfc", ("system","ecutwfc")],
                             ["scf_thr", ("electrons","conv_thr")],
                             ["nspin", ("system","nspin")],
                             ["nx",("system","nr1")],
                             ["ny",("system","nr2")],
                             ["nz",("system","nr3")],
                             ["dft_functional",("system","input_dft")],
                             ["ecutrho", ("system","ecutrho")],
                             ]:
        if para_aba in param:
            qp[para_qe[0]][para_qe[1]] = param.pop(para_aba)
            
    # set potential 
    if "chg_extrap" in param:
        aba2qe_dict = {"atomic": "atomic", "first-order": "first_order", "second-order": "second_order"}
        qp["ions"]["pot_extrapolation"] = aba2qe_dict.get(param.pop("chg_extrap"))

    set_smearing(param,qp)
    set_ks_solver(param,qp)
    set_mixing_type(param,qp)
    set_relax_param(param,qp,calculation)
    set_spin_soc(param,qp)
    set_fixed_axes(param,qp)
    set_cdft(param,qp)
    
    # these paramters will be ignored    
    for ip in ["suffix","basis_type","gamma_only","kpt_file","dft_plus_u","orbital_corr","hubbard_u"]:
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
            if ip in ["kspacing"]: continue
            print("         %s" % ip)
    
    if qe_param:
        for ik,iv in qe_param.items():
            if ik.lower() in ["system","control","electrons","ions","cell"]:
                for ikk,ivv in iv.items():
                    qp[ik.lower()][ikk] = ivv
            else:
                qp[ik] = iv
                
    return qp 

def transfer_mag(imag,nspin=2):
    if nspin == 2:
        if imag == None:
            return 0
        elif isinstance(imag,float):
            return imag
        elif isinstance(imag,list) and len(imag) == 1:
            return imag[0]
        elif isinstance(imag,list) and len(imag) == 3:
            return (imag[0]**2 + imag[1]**2 + imag[2]**2)**0.5
        else:
            print("ERROR: The magnetization is not correct")
            return 0
    elif nspin == 4:
        # return (total_mag, angle1, angle2)
        if imag == None:
            return 0,0,0
        elif isinstance(imag,float):
            return imag,0,0
        elif isinstance(imag,list) and len(imag) == 1:
            return imag[0],0,0
        elif isinstance(imag,list) and len(imag) == 3:
            tot = (imag[0]**2 + imag[1]**2 + imag[2]**2)**0.5
            if abs(tot) < 1e-8:
                return 0,0,0
            else:
                angle1,angle2 = MyAbacus.AbacusStru.mag_to_angle(imag[0],imag[1],imag[2])
                return tot,angle1,angle2
        else:
            print("ERROR: The magnetization is not correct")
            return 0,0,0
    else:
        return 0

def set_stru(abacus_stru,param):
    # containing mag/U
    labels = abacus_stru.get_label(total=True)
    label =abacus_stru.get_label(total=False)
    mass = abacus_stru.get_mass()
    stru = abacus_stru.get_stru()
    coord = stru["coord"]
    mag = abacus_stru.get_atommag() # the mag of each atom
    pp = abacus_stru.get_pp() # the pp of each atom type
    
    nspin = param["system"].get("nspin",1)
    if nspin in [2,4]:
        '''
        1. find all (label,atom_mag) pairs and store in label_mag,
        2. rename the same label in label_mag, and rearrange the coord
        3. 
        '''
        label_mag = []
        for i in range(len(mag)):
            imag = transfer_mag(mag[i],nspin)
            label_mag.append((labels[i],imag))
        
        # need rearrange the coord/label/pp
        new_coord = []
        new_label = []
        new_pp = []
        new_mass = []
        org_atomtype_idx = [] # save the original index of atomtype for plus U
        
        label_mag_types = [] # element is (label,atom_mag)
        label_types = [] # element is label of each new atom type
        
        for il,ilabel in enumerate(label_mag):
            if ilabel not in label_mag_types:
                label_mag_types.append(ilabel)
                label_types.append(ilabel[0])
                new_coord.append([])
            new_coord[label_mag_types.index(ilabel)].append(coord[il])
        
        for i_new_label in label_types:
            new_pp.append(pp[label.index(i_new_label)])
            new_mass.append(mass[label.index(i_new_label)])
            org_atomtype_idx.append(label.index(i_new_label))

        new_label = copy.deepcopy(label_types)
        for i in range(len(label_types)):
            if label_types.count(label_types[i]) > 1:
                new_label[i] = label_types[i] + str(label_types[:i+1].count(label_types[i]))
        
        label = new_label
        pp = new_pp
        mass = new_mass
        atomtype_mag = []
        coord = []
        atom_number = []
        for i,icoord in enumerate(new_coord):
            atom_number.append(len(icoord))
            coord += icoord
            atomtype_mag.append(label_mag_types[i][1])
    else:
        org_atomtype_idx = list(range(len(label)))
        atomtype_mag = abacus_stru.get_mag()
        atom_number = [labels.count(i) for i in label]
    
    # set structure related param
    param["system"]["nat"] = len(coord)
    param["system"]["ntyp"] = len(label)
    param["system"]["ibrav"] = 0
    param["system"]["celldm(1)"] = stru["lat"]
    
    # set init mag
    if nspin == 2:
        for i in range(len(atomtype_mag)):
            if abs(atomtype_mag[i]) > 1e-8:
                param["system"]["starting_magnetization(%d)" % (i+1)] = atomtype_mag[i]
    elif nspin == 4:
        for i in range(len(atomtype_mag)):
            if abs(atomtype_mag[i][0]) > 1e-8:
                param["system"]["starting_magnetization(%d)" % (i+1)] = atomtype_mag[i][0]
                param["system"]["angle1(%d)" % (i+1)] = atomtype_mag[i][1]
                param["system"]["angle2(%d)" % (i+1)] = atomtype_mag[i][2]
    
    new_stru_dict = {
        "cell":stru["cell"],
        "coord":coord,
        "la":stru["lat"],
        "cartesian":stru["cartesian"],
    }
    
    return label,atom_number,pp,mass,new_stru_dict,org_atomtype_idx

def set_plusU(abacus_input,param,label,org_atomtype_idx,version):
    # set plus U
    # label should be the new label list (considering the labels are re-arranged because of the magnetization)
    # org_atomtype_idx is the original index of atomtype for each new label
    
    def find_e(ilabel):
        if len(ilabel) == 1:
            return ilabel[0]
        else:
            if ilabel[:2] in PERIOD_DICT_NUMBER:
                return ilabel[:2]
            elif ilabel[:1] in PERIOD_DICT_NUMBER:
                return ilabel[:1]
            else:
                print("ERROR: Can not find the element of %s" % ilabel)
                return None
            
    if "dft_plus_u" in abacus_input:
        if comm.IsTrue(abacus_input.pop("dft_plus_u")) and "hubbard_u" in abacus_input:
            us = abacus_input.pop("hubbard_u").split()
            if version <= 7.0:
                param["system"]["lda_plus_u"] = ".TRUE."
                for i in range(len(org_atomtype_idx)):
                    param["system"][" Hubbard_U(%d)" % (i+1)] = us[org_atomtype_idx[i]]
            else:
                # U Fe-2d 3.0
                orbital_corr = abacus_input.pop("orbital_corr","").split()
                if orbital_corr:
                    pdf = {"1":"p","2":"d","3":"f"}
                    param["HUBBARD (ortho-atomic)"] = []
                    for i in range(len(org_atomtype_idx)):
                        org_idx = org_atomtype_idx[i]
                        if orbital_corr[org_idx] in ["1","2","3"]:
                            il = label[i]
                            iu = us[org_idx]
                            ie = find_e(il)
                            if not ie:
                                print("ERROR: Can not find the element of %s" % il)
                                continue
                            uorbital = comm.get_period(ie)-int(orbital_corr[org_idx])+1
                            orb = pdf[orbital_corr[org_idx]]
                            param["HUBBARD (ortho-atomic)"].append("U %s-%d%s %s" % (il,uorbital,orb,iu))


def set_kpt(abacus_input,cell,path):
    # stru is a dict from abacus_stru.get_stru()    
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
        kpt = comm.kspacing2kpt(kspacing,[[i for i in j] for j in cell])
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
    return kpt

def Abacus2Qe(path: str= ".", save_path: str = None, qe_param:Dict[str,any] = {}):
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
    
    qe_param_input = copy.deepcopy(qe_param)
    version = qe_param_input.pop("version",7.0)
    param = ParamAbacus2Qe(abacus_input,version,qe_param_input)
    if not param:
        print("ERROR: Transfer abacus input to qe input failed in %s" % path)
        return None
    
    kpt = set_kpt(abacus_input,abacus_stru.get_cell(bohr=True),path)
    if not kpt:
        print("ERROR: Transfer KPT failed in %s" % path)
        return None
    
    # set structure, as well as mag, and may rearrange the label because of the mag
    # the lattice constant lat0 is also setted in this function
    label,atom_number,pp,mass,new_stru_dict,org_atomtype_idx = set_stru(abacus_stru,param)
    
    # set plus U
    set_plusU(abacus_input,param,label,org_atomtype_idx,version)  

    qe = MyQe.QeInputs(param=param,
                    label=label,
                    atom_number=atom_number,
                    cell=new_stru_dict["cell"],
                    coord=new_stru_dict["coord"],
                    pp=pp,
                    kpt = kpt,
                    kpt_type ="automatic",
                    mass=mass,
                    cell_type = "alat",
                    coord_type = "alat" if new_stru_dict["cartesian"] else "crystal")
    
    save_path = save_path if save_path != None else os.path.join(path,"input")
    qe.write(save_path)
    return qe
