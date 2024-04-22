
from ..resultAbacus import ResultAbacus

class RefAbacus(ResultAbacus): 
    # This class is used to compare the results of the current calculation with the reference results.  
    
    @ResultAbacus.register(delta_energy="the difference between energy and the reference value. Unit in eV. Key in reference file is \"energy\"",
                           delta_energyPerAtom="delta_energy/natom, unit in eV")
    def GetDeltaEnergy(self): 
        if self.resultREF.get("energy",None) != None:
            if self["energy"] != None:
                self["delta_energy"] = self["energy"] - self.resultREF.get("energy")
                if self["natom"] != None:
                    self["delta_energyPerAtom"] = (self["energy"] - self.resultREF.get("energy")) / self["natom"]
                else:
                    self["delta_energyPerAtom"] = None
                return
        self["delta_energy"] = None
        self["delta_energyPerAtom"] = None
        
    @ResultAbacus.register(delta_force="the difference between force and the reference value. The reference value should be a list of 3*natom, each list is a force of each atom. Unit in eV/Angstrom",
                           delta_stress="the difference between stress and the reference value. The reference value should be a list 9 elements, unit in kbar"
                           )
    def GetDeltaForceStress(self): 
        df = None
        ds = None
        if self.resultREF.get("force",None) != None:
            if self["force"] != None:
                try:
                    df = [i - j for i,j in zip(self["force"],self.resultREF.get("force"))]
                except:
                    df = None
        if self.resultREF.get("stress",None) != None:
            if self["stress"] != None:
                try:
                    ds = [i - j for i,j in zip(self["stress"],self.resultREF.get("stress"))]
                except:
                    ds = None
        self["delta_force"] = df
        self["delta_stress"] = ds