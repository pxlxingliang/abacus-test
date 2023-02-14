import os,sys
from .result import Result
from . import comm

class ResultAbacus(Result):
    _PARAM_DIC = {}

    def __init__(self,path = ".",output = None):
        super().__init__()
        self.PATH = path   #the path of ABACUS job

        self.INPUTf = os.path.join(self.PATH,'INPUT')
        self.STRUf  = os.path.join(self.PATH,'STRU')
        self.KPTf   = os.path.join(self.PATH,'KPT')
        
        self.INPUT  = comm.ReadFile(self.INPUTf,warn=True)
        self.STRU   = comm.ReadFile(self.STRUf,warn=True)
        self.KPT    = comm.ReadFile(self.KPTf,warn=False)
        
        self.OUTPUTf = output if output != None else comm.FindOutput(self.PATH,keyinfo="Atomic-orbital Based Ab-initio")
        if self.OUTPUTf == None:
            print("WARNING: can not find the output of ABACUS in %s" % self.PATH)
            self.OUTPUT = []
        else:
            self.OUTPUT = comm.ReadFile(self.OUTPUTf,warn=False)
        
        self.SUFFIX,self.CALCULATION = self.SuffixCalculation(self.INPUT) 

        self.LOGf   = os.path.join(self.PATH,"OUT.%s/running_%s.log"%(self.SUFFIX,self.CALCULATION))
        self.LOG    = comm.ReadFile(self.LOGf,warn=True)

    def SuffixCalculation(self,INPUT):
        suffix = "ABACUS"
        calculation = "scf"
        for iline in INPUT:
            sline = iline.split("#")[0].split()
            if len(sline) >= 2 and sline[0].lower() == "suffix":
                suffix = sline[1].strip()
            elif len(sline) >= 2 and sline[0].lower() == "calculation":
                calculation = sline[1].strip()

        return suffix,calculation

