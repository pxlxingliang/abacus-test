import numpy as np
from typing import List
import os, glob, re, json
import traceback
from abacustest.constant import MASS_DICT

def Direct2Cartesian(coord:List[List[float]],cell:List[List[float]]):
    return np.array(coord).dot(np.array(cell)).tolist()

def Cartesian2Direct(coord:List[List[float]],cell:List[List[float]]):
    return np.array(coord).dot(np.linalg.inv(np.array(cell))).tolist()

def translate_strus(input_strus, input_stru_type, output_path = "."):
    """
    Translate the structure from one format to ABACUS stru.
    
    input_strus: str/list, the input structure.
    input_stru_type: str, the input structure type.
    output_stru_type: str, the output structure type.
    
    Return:
    output_strus: list, the output structure.
    """
    import dpdata

    def gen_path_name(base_path):
        if not os.path.exists(base_path):
            return base_path
        idx = 1
        while os.path.exists(f"{base_path}.{idx}"):
            idx += 1
        return f"{base_path}.{idx}"
    
    
    
    dpdata_formats = dpdata.format.Format.get_formats()
    if input_stru_type not in dpdata_formats and input_stru_type.lower() != "cif":
        print("ERROR: input_stru_type should be in cif, %s, but not %s" % (str(dpdata_formats),input_stru_type))
        return None
    
    if isinstance(input_strus,str):
        input_strus = [input_strus]
        
    output_folders = []
    idx = 0
    struinfo = {}
    try:
        for istru in input_strus:
            for iistru in glob.glob(istru):
                if input_stru_type in ["abacus/stru", "stru"]:
                    tpath = gen_path_name(os.path.join(output_path,"%06d" % idx))
                    os.makedirs(tpath,exist_ok=True)
                    os.system("cp %s %s" % (iistru, os.path.join(tpath,"STRU")))
                    output_folders.append(tpath)
                    idx += 1
                    print("Copy %s to %s" % (iistru, os.path.join(tpath,"STRU")))
                    with open(os.path.join(tpath,"struinfo.txt"),"w") as f:
                        f.write(istru)
                    struinfo[istru] = [os.path.basename(tpath)]
                else:    
                    if input_stru_type in dpdata_formats:
                        stru = dpdata.System(iistru,fmt=input_stru_type)
                    elif input_stru_type.lower() == "cif":
                        try:
                            from ase.io import read as ase_read
                            stru = ase_read(iistru)
                            stru = dpdata.System(stru, fmt="ase/structure")
                        except:
                            try:
                                from pymatgen.core import Structure
                                stru = Structure.from_file(iistru)
                                stru = dpdata.System(stru, fmt="pymatgen/structure")
                            except:
                                traceback.print_exc()
                                raise Exception("Cannot read cif file %s with ase or pymatgen" % iistru)
                    
                    print("Translating %s to ABACUS stru:" % iistru)
                    struinfo[istru] = []
                    for i in range(stru.get_nframes()):
                        tpath = os.path.join(output_path,"%06d" % idx)
                        os.makedirs(tpath,exist_ok=True)
                        stru.to("abacus/stru", os.path.join(tpath,"STRU"),i, pp_file=["" for _ in stru.data["atom_names"]])
                        output_folders.append(tpath)
                        idx += 1
                        print("    Save to %s" % os.path.join(tpath,"STRU"))
                        with open(os.path.join(tpath,"struinfo.txt"),"w") as f:
                            f.write(istru)
                        struinfo[istru].append(os.path.basename(tpath))
    except:
        traceback.print_exc()
        print("ERROR: %s to ABACUS STRU failed" % (input_stru_type))
        return None
    if len(struinfo) > 0:
        json.dump(struinfo, open(os.path.join(output_path,"struinfo.json"),"w"), indent=4)
    return output_folders

def read_pp_valence(pp_file):
    if not os.path.isfile(pp_file):
        print("ERROR: %s is not a file" % pp_file)
        return None
    
    with open(pp_file) as f: lines = f.readlines()
    for iline in lines:
        if "z_valence" in iline:
            #split the line by z_valence and get the second part, then split by ", and get the second part
            zv = re.split("z_valence",iline)[1].split("\"")[1]
            return float(zv)
        elif "Z valence" in iline:
            return int(iline.split()[0])
    return None

def get_element_name_from_file(filename):
        #the filename should be started with the element name and followed by character non-alpha
        def check_element(element):
            if element in MASS_DICT:
                return element
            else:
                return None
        filename = os.path.basename(filename)
        if filename == "":
            return None
        if len(filename) == 1:
            return check_element(filename)
        
        element_name = filename[:2]
        if element_name[-1].isalpha():
            return check_element(element_name)
        else:
            return check_element(filename[0])

def collect_pp(pp_path):
    # Read the pp_path and collect the pp files
    # pp_path: the path of the pp files
    
    if pp_path is None:
        return {}
    
    pp = {}
    if os.path.isdir(pp_path):
        if os.path.isfile(os.path.join(pp_path,"element.json")):
            for key,value in json.load(open(os.path.join(pp_path,"element.json"))).items():
                if os.path.isfile(os.path.join(pp_path,value)):
                    pp[key] = os.path.join(pp_path,value)
        else:
            allfiles = os.listdir(pp_path)
            for ifile in allfiles:
                if not os.path.isfile(os.path.join(pp_path,ifile)): continue
                element_name = get_element_name_from_file(ifile)
                if element_name is not None:
                    pp[element_name] = os.path.join(pp_path,ifile)
    else:
        print(f"Not find pp dir: \'{pp_path}\'\n\tcurrent path: {os.getcwd()}")
    return pp

def kspacing2kpt(kspacing, cell):
    """
    Convert kspacing to kpt.
    """
    if isinstance(kspacing, float):
        kspacing = [kspacing, kspacing, kspacing]
    elif isinstance(kspacing, str):
        a = kspacing.split()
        if len(a) == 1:
            kspacing = [float(a[0]), float(a[0]), float(a[0])]
        elif len(a) == 3:
            kspacing = [float(a[0]), float(a[1]), float(a[2])]
        else:
            raise ValueError("kspacing must be one or three floats")
    elif not isinstance(kspacing, list):
        raise TypeError("kspacing must be float or list")
    
    assert len(kspacing) == 3, "kspacing must be 3-dim"
    kpt = []
    V = abs(np.linalg.det(np.array(cell)))
    for i in range(3):
        kpt.append(int(np.linalg.norm(np.cross(cell[(i+1)%3],cell[(i+2)%3])) * 2 * np.pi / V / kspacing[i] + 1))
    return kpt

def kpt2kspacing(kpt,cell):
    """
    Convert kpt to kspacing.
    """
    if isinstance(kpt, int):
        kpt = [kpt, kpt, kpt]
    elif isinstance(kpt, str):
        a = kpt.split()
        if len(a) == 1:
            kpt = [int(a[0]), int(a[0]), int(a[0])]
        elif len(a) == 3:
            kpt = [int(a[0]), int(a[1]), int(a[2])]
        else:
            raise ValueError("kpt must be one or three ints")
    elif not isinstance(kpt, list):
        raise TypeError("kpt must be int or list")
    
    assert len(kpt) == 3, "kpt must be 3-dim"
    kspacing = []
    V = abs(np.linalg.det(np.array(cell)))
    for i in range(3):
        kspacing.append(np.linalg.norm(np.cross(cell[(i+1)%3],cell[(i+2)%3])) * 2 * np.pi / V / kpt[i])
    
    # remain 6 digits,
    kspacing = [round(i,6)+1e-5 for i in kspacing]
    return kspacing

def IsTrue(param):
    '''
    judge if a parameter is True
    
    If param is a string, then judge if it is "True" or "true" or "T" or "t" or "1", if yes return True, elif is "False" or "false" or "F" or "f" or "0", return False, else return None.
    If param is int or bool, return True if param is True, else return False.
    
    '''
    if isinstance(param,str):
        param = param.rstrip(".").lstrip(".")  # for the case of "True."
        if param.lower() in ["true","t","1"]:
            return True
        elif param.lower() in ["false","f","0"]:
            return False
        else:
            return None
    elif isinstance(param,(int,bool)):
        return True if param else False
    else:
        return None

def get_period(element):
    if isinstance(element,str):
        from abacustest import constant
        element = constant.PERIOD_DICT_NUMBER.get(element,0)
    if not isinstance(element,int):
        print("ERROR: element should be a str or int, but not %s" % str(element))
        return None
    
    if element in [1,2]:
        return 1
    elif element <= 10:
        return 2
    elif element <= 18:
        return 3
    elif element <= 36:
        return 4
    elif element <= 54:
        return 5
    elif element <= 86:
        return 6
    elif element <= 118:
        return 7
    else:
        print("ERROR: element should be in [1,118], but not %d" % element)
        return None

def perturb_cell(cell, perturb_ratio,coord=None):
    '''
    Perturb the cell and coord.
    
    cell: 3x3 list, the cell vectors.
    coord: Nx3 list, the coordinates.
    perturb_ratio: float, the perturb ratio; or list, the perturb ratio range.
    perturb_num: int, the number of perturbation.
    
    Return:
    new_cell: list, the perturbed cell.
    new_coord: list, the perturbed coord
    '''
    if isinstance(perturb_ratio, (float,int)):
        pratio = [0, abs(perturb_ratio)]
    elif isinstance(perturb_ratio, list):
        pratio = [abs(i) for i in perturb_ratio]
        pratio = [min(pratio),max(pratio)]
    else:
        print("ERROR: perturb_ratio should be a float or a list, but not %s" % str(perturb_ratio))
        return cell, coord

    cell = np.array(cell)

    r6 = np.random.uniform(pratio[0],pratio[1],6) * (np.random.randint(0,2,6)*2-1)
    
    perturb_matrix = np.array([[1+r6[0], 0.5*r6[3], 0.5*r6[4]],
                                  [0.5*r6[3], 1+r6[1], 0.5*r6[5]],
                                  [0.5*r6[4], 0.5*r6[5], 1+r6[2]]])
    new_cell = np.dot(cell, perturb_matrix).tolist()
    if coord is None:
        return new_cell, None
    else:
        coord = np.array(coord)
        new_coord = np.dot(coord, perturb_matrix).tolist()  
    return new_cell, new_coord   

def perturb_coord(coord:List[List], atom_pert_distance):
    '''
    Perturb the coordinates with the max distance in the range of atom_pert_distance.
    
    coord: Nx3 list, the coordinates.
    atom_pert_distance: float, the maximum perturb distance; or list, the perturb distance range.
    '''
    if isinstance(atom_pert_distance, (float,int)):
        atom_pert_distance = [0, abs(atom_pert_distance)]
    elif isinstance(atom_pert_distance, list):
        atom_pert_distance = [abs(i) for i in atom_pert_distance]
        atom_pert_distance = [min(atom_pert_distance),max(atom_pert_distance)]
    else:
        print("ERROR: atom_pert_distance should be a float or a list, but not %s" % str(atom_pert_distance))
        return coord

    new_coord = []
    for icoord in coord:
        icoord = np.array(icoord)
        random_vector = np.array([0,0,0])
        e = np.random.rand(3) 
        while np.linalg.norm(e) < 0.1:
            e = np.random.rand(3)
        e = e * 2 - 1
        random_unit_vector = e / np.linalg.norm(e)
        random_vector = random_unit_vector * np.random.uniform(atom_pert_distance[0],atom_pert_distance[1])
        new_coord.append((icoord + random_vector).tolist())
    return new_coord

def pert_vector(vectors:List[List[float]],max_angle):
    """
    perturb the vectors with a random angle.
    
    vectors: list, the vectors.
    max_angle: float, the max angle; or list, the angle range.
        
    Based on Rodrigues' Rotation Formula to rotate a vector.
    assume rotate angle is theta, rotate axis is k, the vector is v.
    v_rot = R(k,theta) * v
    R(k,theta) = I + sin(theta) * K + (1-cos(theta)) * K^2
    where K is the skew-symmetric matrix of k.
    K = [[0, -kz, ky],
         [kz, 0, -kx],
         [-ky, kx, 0]]

    R = [[kx^2*(1-cos(theta))+cos(theta), kx*ky*(1-cos(theta))-kz*sin(theta), kx*kz*(1-cos(theta))+ky*sin(theta)],
         [ky*kx*(1-cos(theta))+kz*sin(theta), ky^2*(1-cos(theta))+cos(theta), ky*kz*(1-cos(theta))-kx*sin(theta)],
         [kz*kx*(1-cos(theta))-ky*sin(theta), kz*ky*(1-cos(theta))+kx*sin(theta), kz^2*(1-cos(theta))+cos(theta)]]
    """
    if isinstance(max_angle, (float,int)):
        max_angle = [0, abs(max_angle)]
    elif isinstance(max_angle, list):
        max_angle = [abs(i) for i in max_angle]
        max_angle = [min(max_angle),max(max_angle)]
    else:
        print("ERROR: max_angle should be a float or a list, but not %s" % str(max_angle))
        return vectors
    
    #ref_vector = None
    #for ivector in vectors:
    #    if np.array(ivector).any() != 0:
    #        ref_vector = ivector
    #        break
    #if ref_vector is None:
    #    return vectors
    
    # generate a reandom ref_vector
    ref_vector = np.array([0,0,0])
    while np.linalg.norm(ref_vector) < 0.1:
        ref_vector = np.random.rand(3)
    ref_vector = ref_vector / np.linalg.norm(ref_vector)
    
    axis = np.random.rand(3)
    while np.linalg.norm(axis) < 0.1 or np.cross(axis,ref_vector).any() == 0:
        axis = np.random.rand(3)
    
    axis = np.cross(axis,ref_vector)
    axis = axis / np.linalg.norm(axis)
    angle = np.random.uniform(max_angle[0],max_angle[1]) * np.pi / 180
    skew_symmetric_matrix = np.array([[0, -axis[2], axis[1]],
                                      [axis[2], 0, -axis[0]],
                                      [-axis[1], axis[0], 0]])
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    R_matrix = np.eye(3) + sin_theta * skew_symmetric_matrix + (1-cos_theta) * np.dot(skew_symmetric_matrix,skew_symmetric_matrix)
    new_vectors = []
    for ivector in vectors:
        new_vectors.append(np.dot(R_matrix,ivector).tolist())
    return new_vectors
        