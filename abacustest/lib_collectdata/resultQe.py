import os,sys,json
from .result import Result
from . import comm

class ResultQe(Result):
    _PARAM_DIC = {}

    def __init__(self,path = ".",output = None,resultREF="resultREF.json"):
        super().__init__()
        self.PATH = path   #the path of QE job

        self.XMLf = os.path.join(self.PATH,"pwscf.xml")
        self.OUTPUTf = output if output != None else comm.FindOutput(self.PATH,keyinfo='Program PWSCF')

        self.OUTPUT = comm.ReadFile(self.OUTPUTf,warn=False)
        self.XMLROOT = comm.ReadXmlFile(self.XMLf,warn=True)

        self.resultREF = {}
        if resultREF and os.path.isfile(resultREF):
            self.resultREF = json.load(open(resultREF))
