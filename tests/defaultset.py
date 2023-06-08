import os
from dflow.python import upload_packages
root_path = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]
upload_packages.append(os.path.join(root_path,"abacustest"))

STRU1 = '''ATOMIC_SPECIES
Si 1.000 Si.upf  #Element, Mass, Pseudopotential

NUMERICAL_ORBITAL
Si.orb

LATTICE_CONSTANT
10.2                    #Lattice constant

LATTICE_VECTORS
0.5 0.5 0.0             #Lattice vector 1
0.5 0.0 0.5             #Lattice vector 2
0.0 0.5 0.5             #Lattice vector 3

ATOMIC_POSITIONS
Cartesian               #Cartesian(Unit is LATTICE_CONSTANT)
Si                      #Name of element
0.0                     #Magnetic for this element.
2                       #Number of atoms
0.00 0.00 0.00 0 0 0    #x,y,z, move_x, move_y, move_z
0.25 0.25 0.25 1 1 1
'''

STRU2 = '''ATOMIC_SPECIES
Si 1.000 pplib/Si.upf1  #Element, Mass, Pseudopotential

NUMERICAL_ORBITAL
orblib/Si.orb1

LATTICE_CONSTANT
5.1                    #Lattice constant

LATTICE_VECTORS
0.5 0.5 0.0             #Lattice vector 1
0.5 0.0 0.5             #Lattice vector 2
0.0 0.5 0.5             #Lattice vector 3

ATOMIC_POSITIONS
Cartesian               #Cartesian(Unit is LATTICE_CONSTANT)
Si                      #Name of element
0.0                     #Magnetic for this element.
2                       #Number of atoms
0.00 0.00 0.00 0 0 0    #x,y,z, move_x, move_y, move_z
0.25 0.25 0.25 1 1 1
'''

INPUT1 = """INPUT_PARAMETERS
ecutwfc                 50
basis_type              lcao
"""

INPUT2 = """INPUT_PARAMETERS
basis_type              pw
"""

KPT1 = """K_POINTS
0
Gamma
4 4 4 0 0 0
"""

KPT2 = """K_POINTS
0
Gamma
2 2 2 0 0 0
"""