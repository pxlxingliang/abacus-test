import numpy as np
from typing import List

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
    perturb_ratio: float, the perturb ratio.
    perturb_num: int, the number of perturbation.
    
    Return:
    new_cell: list, the perturbed cell.
    new_coord: list, the perturbed coord
    '''
    perturb_ratio = abs(perturb_ratio)
    cell = np.array(cell)

    r6 = (np.random.rand(6) - 0.5 ) * 2 * perturb_ratio
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

def perturb_coord(coord:List[List], atom_pert_distance, atom_pert_style="normal"):
    '''
    atom_pert_style: str, the style of perturbation.
        - `'normal'`: the `distance` will be object to `chi-square distribution with 3 degrees of freedom` after normalization.
            The mean value of the distance is `atom_pert_fraction*side_length`
        - `'uniform'`: will generate uniformly random points in a 3D-balls with radius as `atom_pert_distance`.
            These points are treated as vector used by atoms to move.
            Obviously, the max length of the distance atoms move is `atom_pert_distance`.
        - `'const'`: The distance atoms move will be a constant `atom_pert_distance`.
    '''
    atom_pert_distance = abs(atom_pert_distance)

    new_coord = []
    for icoord in coord:
        icoord = np.array(icoord)
        random_vector = np.array([0,0,0])
        if atom_pert_style == "normal":
            e = np.random.rand(3) * 2 - 1 
            random_vector = (atom_pert_distance / np.sqrt(3)) * e
        elif atom_pert_style == "uniform":
            e = np.random.rand(3)
            while np.linalg.norm(e) < 0.1:
                e = np.random.rand(3) 
            e = e * 2 - 1
            random_unit_vector = e / np.linalg.norm(e)
            v0 = np.random.rand(1)
            v = np.power(v0, 1 / 3)
            random_vector = atom_pert_distance * v * random_unit_vector
        elif atom_pert_style == "const":
            e = np.random.rand(3) 
            while np.linalg.norm(e) < 0.1:
                e = np.random.rand(3)
            e = e * 2 - 1
            random_unit_vector = e / np.linalg.norm(e)
            random_vector = atom_pert_distance * random_unit_vector
        else:
            print(f"unsupported options atom_pert_style={atom_pert_style}")
        new_coord.append((icoord + random_vector).tolist())
    return new_coord

def pert_vector(vectors:List[List[float]],max_angle):
    """
    perturb the vectors with a random angle.
    
    vectors: list, the vectors.
    max_angle: float, the max angle.
        
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

    axis = np.random.rand(3)

    while np.linalg.norm(axis) < 0.1 or np.cross(axis,np.array(vectors[0])).any() == 0:
        axis = np.random.rand(3)
    axis = np.cross(axis,vectors[0])
    axis = axis / np.linalg.norm(axis)
    angle = np.random.rand() * max_angle * np.pi / 180
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
        