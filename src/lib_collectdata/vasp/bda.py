from ..resultVasp import ResultVasp
from ..comm_funcs.bda import BasicPropertyVasp
import os

class BDA(ResultVasp):
    @ResultVasp.register( bda_mag_moment="mag_moment of some metal element",
                          bda_bond_length="bond_length of some metal element")
    def GetBdaInfo(self):
        cwd = os.getcwd()
        os.chdir(self.PATH)
        output = BasicPropertyVasp('OUTCAR', 'CONTCAR')
        self["bda_mag_moment"] = output.magnetic_moment('TM')
        self["bda_bond_length"] = output.bond_length()
        os.chdir(cwd)
