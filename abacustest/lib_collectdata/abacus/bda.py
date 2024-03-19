import os,sys,glob,re
from ..resultAbacus import ResultAbacus

class BdaAbacus(ResultAbacus):
    
    @ResultAbacus.register(bda_mag_moment="mag_moment of some metal element",
                           bda_bond_length="bond_length of some metal element")
    def GetBDAinfo(self):
        from pymatgen.core.structure import Structure
        from pymatgen.io.vasp.inputs import Poscar
        structure = Structure(lattice=self['cell'],
                              species=self['element_list'],
                              coords=self["coordinate"],
                              coords_are_cartesian=True)
        mag = tuple([{"tot":i} for i in self['atom_mag']])
        from ..comm_funcs.bda import BasicProperty
        output = BasicProperty(mag,None,None,Poscar(structure))
        mag = output.magnetic_moment('TM')
        bond = output.bond_length()
        if mag:
            self["bda_mag_moment"] = mag
        else:
            self["bda_mag_moment"] = None
        if bond:
            self["bda_bond_length"] = bond
        else: 
            self["bda_bond_length"] = None

        
        