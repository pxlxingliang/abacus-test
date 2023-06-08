import os,sys,json
from .result import Result
from . import comm

class ResultVasp(Result):
    _PARAM_DIC = {}

    def __init__(self,path = ".",resultREF="resultREF.json"):
        super().__init__()
        self.PATH = path   #the path of VASP job

        self.INCARf = os.path.join(self.PATH,'INCAR')
        self.KPOINTSf = os.path.join(self.PATH,'KPOINTS')
        self.POSCARf = os.path.join(self.PATH,'POSCAR')
        self.OUTCARf = os.path.join(self.PATH,'OUTCAR')
        self.OSZICARf = os.path.join(self.PATH,'OSZICAR')
        self.XMLf = os.path.join(self.PATH,'vasprun.xml')

        self.INCAR      = comm.ReadFile(self.INCARf,warn=False)
        self.KPOINTS    = comm.ReadFile(self.KPOINTSf,warn=False)
        self.POSCAR     = comm.ReadFile(self.POSCARf,warn=False)
        self.OUTCAR     = comm.ReadFile(self.OUTCARf,warn=True)
        self.OSZICAR    = comm.ReadFile(self.OSZICARf,warn=False)
        self.XMLROOT    = comm.ReadXmlFile(self.XMLf,warn=True)

        self.resultREF = {}
        if resultREF and os.path.isfile(resultREF):
            self.resultREF = json.load(open(resultREF))


