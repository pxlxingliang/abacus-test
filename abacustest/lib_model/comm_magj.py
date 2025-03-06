import numpy as np

def cal_norm_vector(v1, v2):
    """
    Calculate the normal vector of two vectors
    
    Parameters
    ----------
    v1 : list of three floats
        The first vector
    v2 : list of three floats
        The second vector
        
    Returns
    -------
    normal_v : list of three floats
        The normal vector of v1 and v2
    """
    v1 = np.array(v1) / np.linalg.norm(v1)
    v2 = np.array(v2) / np.linalg.norm(v2)
    dot = np.dot(v1, v2)
    if abs(dot - 1) < 1e-6:
        # the two vectors are parallel, then we need to find a vector that is not parallel to v1
        for i in [[1,0,0], [0,1,0], [0,0,1]]:
            if abs(np.dot(v1, i) - 1) > 1e-6:
                normal_v = np.cross(v1, i)
                return normal_v / np.linalg.norm(normal_v)
    else:
        normal_v = np.cross(v1, v2)
        return normal_v / np.linalg.norm(normal_v)

def cal_angle(v1, v2):
    """
    Calculate the angle between two vectors
    
    Parameters
    ----------
    v1 : list of three floats
        The first vector
    v2 : list of three floats
        The second vector
        
    Returns
    -------
    angle : float
        The angle between v1 and v2 in degree
    """
    v1 = np.array(v1) / np.linalg.norm(v1)
    v2 = np.array(v2) / np.linalg.norm(v2)
    dot = np.dot(v1, v2)
    return np.arccos(dot) * 180 / np.pi

def prepare_atom_mag(mag1, mag2, step, number, together_half=True):
    """
    Prepare the magnetic moment of the atoms
    
    Parameters
    ----------
    mag1 : list of three floats
        The magnetic moment of the first atom
    mag2 : list of three floats
        The magnetic moment of the second atom
    step : float
        The step size of the tilting angle
    number : int
        The number of tilting angles
    together_half : bool
        Whether the tilting angle of tilting atom1 and atom2 together is half of the tilting angle of tilting atom1 and tilting atom2 separately    
    
        
    Returns
    -------
    mag : list of lists
        The dimension is number * 3 * 2 * 3.
        Like: [[[mag1, mag2],[mag1, mag2], [mag1, mag2]],...]
        For each tilting angle, will return three pairs of magnetic moments, which are tilting atom1, tilting atom2 and tilting atom1 and atom2 together. 
    """
    
    mags = []
    norm1 = np.linalg.norm(mag1)
    norm2 = np.linalg.norm(mag2)
    v1 = mag1 / norm1
    v2 = mag2 / norm2
    angle0 = cal_angle(v1, v2) # the angle between mag1 and mag2
    normal_v = cal_norm_vector(v1, v2) # the normal vector of mag1 and mag2
    normal_v1 = cal_norm_vector(v1, normal_v) # the normal vector of mag1 and normal_v
    normal_v2 = cal_norm_vector(normal_v, v2) # the normal vector of mag2 and normal_v
    
    coef_together = 0.5 if together_half else 1
    
    for i in range(number):
        d_angle = (i+1) * step
        angle1 = angle0 + d_angle
        angle2 = angle0 + d_angle * coef_together
        mag = []
        # case1: tilting atom1
        mag1_tilt = v2 * np.cos(angle1 * np.pi / 180) * norm1 + normal_v2 * np.sin(angle1 * np.pi / 180) * norm1
        mag.append([mag1_tilt, mag2])
        
        # case2: tilting atom2
        mag2_tilt = v1 * np.cos(angle1 * np.pi / 180) * norm2 + normal_v1 * np.sin(angle1 * np.pi / 180) * norm2
        mag.append([mag1, mag2_tilt])
        
        # case3: tilting atom1 and atom2 together
        mag1_tilt = v2 * np.cos(angle2 * np.pi / 180) * norm1 + normal_v2 * np.sin(angle2 * np.pi / 180) * norm1
        mag2_tilt = v1 * np.cos(angle2 * np.pi / 180) * norm2 + normal_v1 * np.sin(angle2 * np.pi / 180) * norm2
        mag.append([mag1_tilt, mag2_tilt])
        
        mags.append(mag)
    
    if True:
        # check the result
        print("Original magnetic moments:")
        print("Atom1:", mag1)
        print("Atom2:", mag2)
        print("The angle between the two magnetic moments:", angle0)
        print("The normal vector of the two magnetic moments:", normal_v)
        print("The normal vector of mag1 and normal_v:", normal_v1)
        print("The normal vector of mag2 and normal_v:", normal_v2)
        print("The prepared magnetic moments:")
        for i, mag in enumerate(mags):
            for j, mag_pair in enumerate(mag):
                print(f"Angle {i}, Case {j+1}:", "%.3f" % cal_angle(mag_pair[0], mag_pair[1]), mag_pair)
    
    return np.array(mags).tolist()
        

