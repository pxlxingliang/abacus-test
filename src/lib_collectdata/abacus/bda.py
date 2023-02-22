import os,sys,glob,re
from ..resultAbacus import ResultAbacus
import dpdata

class BdaAbacus(ResultAbacus):
    
    @ResultAbacus.register(bda_mag_moment="mag_moment of some metal element",
                           bda_bond_length="bond_length of some metal element")
    def GetBDAinfo(self):
        from pymatgen.core.structure import Structure
        from pymatgen.io.vasp.inputs import Poscar
        ls = dpdata.LabeledSystem(self.PATH,fmt="abacus/relax")
        atomname = ls.data["atom_names"]
        elementlist = [atomname[i] for i in ls.data["atom_types"]]
        structure = Structure(lattice=ls.data["cells"][-1],
                              species=elementlist,
                              coords=ls.data["coords"][-1],
                              coords_are_cartesian=True)
        mag = tuple([{"tot":i} for i in self['atom_mag']])
        from ..comm_funcs.bda import BasicProperty
        output = BasicProperty(mag,None,None,Poscar(structure))
        self["bda_mag_moment"] = output.magnetic_moment('TM')
        self["bda_bond_length"] = output.bond_length()

        
        