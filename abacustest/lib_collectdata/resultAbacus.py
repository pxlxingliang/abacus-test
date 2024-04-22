import os,sys,json
import traceback
from .result import Result
from . import comm

class ResultAbacus(Result):
    _PARAM_DIC = {}

    
    def __init__(self, path=".", output=None, resultREF=None, resultCriteria=None):
        """
        Initializes an instance of the ResultAbacus class.

        Args:
            path (str, optional): The path of the ABACUS job. Defaults to ".".
            output (str, optional): The screen output file. Defaults to None, and will be found automatically.
            resultREF (str, optional): The path to the resultRef.json file. Defaults to None. If None, resultREF = os.path.join(self.PATH, "resultRef.json").
        """
        super().__init__()
        self.PATH = path   # the path of ABACUS job

        self.INPUTf = os.path.join(self.PATH, 'INPUT')
        self.STRUf = os.path.join(self.PATH, 'STRU')
        self.KPTf = os.path.join(self.PATH, 'KPT')

        self.INPUT = comm.ReadFile(self.INPUTf, warn=True)
        self.STRU = comm.ReadFile(self.STRUf, warn=True)
        self.KPT = comm.ReadFile(self.KPTf, warn=False)
        
        # screen output
        self.OUTPUTf = output if output is not None else comm.FindOutput(self.PATH, keyinfo="Atomic-orbital Based Ab-initio")
        if self.OUTPUTf is None:
            print("WARNING: can not find the output of ABACUS in %s" % self.PATH)
            self.OUTPUT = []
        else:
            self.OUTPUT = comm.ReadFile(self.OUTPUTf, warn=False)

        self.SUFFIX, self.CALCULATION = self.SuffixCalculation(self.INPUT)

        # OUT.XXX/running_xxx.log
        self.LOGf = os.path.join(self.PATH, "OUT.%s/running_%s.log" % (self.SUFFIX, self.CALCULATION))
        self.LOG = comm.ReadFile(self.LOGf, warn=True)
        
        # time.json
        self.TIME = self.ReadTime()

        # read from abacus.json file
        self.JSON = self.ReadJson()

        # resultREF is a json file of a dict, where key is the param name, and value is the value of the param
        # such as : {"energy": -10000.1111111, "force": [[0.1,0.1,0.1],[0.1,0.1,0.1],[0.1,0.1,0.1]]}
        if resultREF is None:
            resultREF = os.path.join(self.PATH, "resultRef.json")
        self.resultREF = {}
        if os.path.isfile(resultREF):
            try:
                self.resultREF = json.load(open(resultREF))
            except:
                traceback.print_exc()
                print("WARNING: can not read resultRef.json")
                self.resultREF = {}
        

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
    
    def ReadTime(self): 
        if os.path.isfile(os.path.join(self.PATH,"time.json")):
            time_file = os.path.join(self.PATH,"time.json")
        elif os.path.isfile(os.path.join(self.PATH,f"OUT.{self.SUFFIX}","time.json")):
            time_file = os.path.join(self.PATH,f"OUT.{self.SUFFIX}","time.json")
        else:
            time_file = None
        
        if time_file:
            import json
            try:
                return json.load(open(time_file))
            except:
                traceback.print_exc()
                return None
        else:
            return None

    def ReadJson(self):
        return None
        if os.path.isfile(os.path.join(self.PATH,"abacus.json")):
            try:
                jsons = json.load(open(os.path.join(self.PATH,"abacus.json")))
                return jsons
            except:
                print("WARNING: can not read abacus.json")
                traceback.print_exc()
                return None
        return None
    
    def GetTime(self,class_name,func_name):
        # return the cpu_second of the function, and the calls of the function
        # if tunc_name is None, return the cpu_second of the class and the number of functions in the class
        '''
        time.json:
	    {
	        "total": 1.0,
	        "sub": [
	            {
	                "class_name": "wavefunc",
	                "sub": [
	                    {
	                        "name": "evc",
	                        "cpu_second": 0.000318,
	                        "calls": 2,
	                        "cpu_second_per_call": 0.000159,
	                        "cpu_second_per_total": 0.000318
	                    }
	                ]
	            }
	        ]
	    }
        '''
        if self.TIME == None:
            return (None,None)
        if class_name == "total":
            return (self.TIME.get("total",None),1)
        if "sub" not in self.TIME:
            return (None,None)
        for sub in self.TIME["sub"]:
            if sub.get("class_name",None) == class_name:
                if "sub" not in sub:
                    return (None,None)
                if func_name == None:
                    total_time = 0.0
                    for subsub in sub["sub"]:
                        total_time += subsub.get("cpu_second",0.0)
                    return (total_time,len(sub["sub"]))
                else:
                    for subsub in sub["sub"]:
                        if subsub.get("name",None) == func_name:
                            return (subsub.get("cpu_second",None),subsub.get("calls",None))
                break
        return (None,None)
