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