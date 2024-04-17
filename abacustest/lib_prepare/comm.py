import numpy as np

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
        