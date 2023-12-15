import os,sys,glob,re
from ..resultAbacus import ResultAbacus
from .. import comm

class Abacus(ResultAbacus):
    
    @ResultAbacus.register(version="the version of ABACUS")
    def GetVersion(self):
        if len(self.LOG) > 0:
            for line_idx,line in enumerate(self.LOG):
                if "WELCOME TO ABACUS" in line:
                    version = line.split()[-1]
                    if version[0].lower() != 'v':
                        print("Unknow version of '%s'" % version)
                    self['version'] = version
                    return
                elif line[30:36] == "ABACUS":
                    version = line[36:].strip()
                    commit = "unknown"
                    for ii in range(30):
                        if line_idx + ii >= len(self.LOG):
                            break
                        if "Commit:" in self.LOG[line_idx + ii]:
                            commit = re.split(":",self.LOG[line_idx + ii].strip(),maxsplit=1)[-1].strip()
                            break
                    self['version'] = version + "(" + commit + ")"
                    return
                                              
    @ResultAbacus.register(ncore="the mpi cores")
    def GetNcore(self):
        for line in self.LOG:
            if "DSIZE =" in line:
                self['ncore'] = int(line.split()[-1])
                return
    
    @ResultAbacus.register(normal_end="if the job is normal ending")
    def GetNormalEnd(self):
        if len(self.LOG) == 0:
            self['normal_end'] = None
            return
        elif len(self.LOG) > 0 and " Total  Time  :" in self.LOG[-1]:
            self['normal_end'] = True
            return
        else:
            self['normal_end'] = False

            print("Job is not normal ending!!! The latest 10 lines is:")
            if len(self.LOG) < 10:
                print(''.join(self.LOG))
            else:
                print(''.join(self.LOG[-10:]))

    @ResultAbacus.register(INPUT="a dict to store the setting in OUT.xxx/INPUT")
    def GetInputParameter(self):
        def str2intfloat(ii):
            try:
                return int(ii)
            except:
                pass
            try:
                return float(ii)
            except:
                return ii

        INPUTf = os.path.join(self.PATH,"OUT.%s/INPUT" % self.SUFFIX)
        if os.path.isfile(INPUTf):
            with open(INPUTf) as f1:
                input_context = f1.readlines()
        else:
            input_context = self.INPUT
        
        readinput = False
        INPUT = {}
        for i,iline in enumerate(input_context):
            if iline.strip() == 'INPUT_PARAMETERS':
                readinput = True
            elif iline.strip() == '' or iline.strip()[0] in ['#']:
                continue
            elif readinput:
                sline = re.split('[ \t]',iline.split("#")[0].strip(),maxsplit=1)
                if len(sline) == 2:
                    INPUT[sline[0].lower()] = str2intfloat(sline[1].strip())
        self["INPUT"] = INPUT
    
    @ResultAbacus.register(kpt="list, the K POINTS setting in KPT file")
    def GetKptParam(self):
        if len(self.KPT) > 3:
            self["kpt"] = [int(i) for i in self.KPT[3].split()[:3]]
        else:
            self["kpt"] = None


    @ResultAbacus.register(nbands="number of bands",
                           nkstot = "total K point number",
                           ibzk = "irreducible K point number",
                           natom ="total atom number",
                           nelec = "total electron number",
                           fft_grid = "fft grid for charge/potential",
                           point_group="point group",
                           point_group_in_space_group="point group in space group")
    def GetLogParam(self):       
        natom = 0
        nelec = 0
        point_group = None
        point_group_in_space_group = None
        for i,line in enumerate(self.LOG):
            if "NBANDS =" in line:
                self['nbands'] = int(line.split()[2])
            elif 'nkstot =' in line:
                self['nkstot'] = int(line.split()[-1])
            elif 'nkstot_ibz =' in line:
                self['ibzk'] = int(line.split()[-1])
            elif 'number of atom for this type =' in line:
                natom += int(line.split()[-1])
            elif 'total electron number of element' in line:
                nelec += float(line.split()[-1])
            elif "[fft grid for charge/potential] =" in line:
                self['fft_grid'] = [float(i.strip()) for i in line.split('=')[1].split(',')]
            elif "POINT GROUP =" in line:
                point_group = line.split("=")[1].strip()
            elif "POINT GROUP IN SPACE GROUP =" in line:
                point_group_in_space_group = line.split("=")[1].strip()

        self["point_group"] = point_group
        self["point_group_in_space_group"] = point_group_in_space_group
        if natom > 0:
            self["natom"] = natom 
        if nelec > 0:
            self["nelec"] = nelec 

    @ResultAbacus.register(converge="if the SCF is converged",
                           total_mag="total magnetism (Bohr mag/cell)",
                           absolute_mag="absolute magnetism (Bohr mag/cell)",
                           energy = "the total energy (eV)",
                           volume = "the volume of cell, in A^3",
                           efermi = "the fermi energy (eV)",
                           energy_per_atom="the total energy divided by natom, (eV)")
    def GetLogResult(self):       
        total_mag = None
        absolute_mag = None
        efermi = None
        converge = None
        energy = None
        volume = None
        for i,line in enumerate(self.LOG):
            if 'charge density convergence is achieved' in line:
                converge = True
            elif 'convergence has NOT been achieved!' in line or\
                'convergence has not been achieved' in line:
                converge = False
            elif 'total magnetism (Bohr mag/cell)' in line:
                total_mag = float(line.split()[-1])
            elif 'absolute magnetism' in line:
                absolute_mag = float(line.split()[-1])
            elif "!FINAL_ETOT_IS" in line:
                energy = float(line.split()[1])
            elif "Volume (A^3) =" in line:
                volume = float(line.split()[-1])
            elif 'E_Fermi' in line:
                if 'E_Fermi_dw' in line:
                    if efermi == None:
                        efermi = float(line.split()[-1])
                        continue
                    if float(line.split()[-1]) > efermi:
                        efermi = float(line.split()[-1])
                else:
                    efermi = float(line.split()[-1])

        self["energy"] = energy
        self["converge"] = converge
        self["volume"] = volume
        
        if total_mag != None:
            self['total_mag'] = total_mag
        if absolute_mag != None:
            self['absolute_mag'] = absolute_mag
        if efermi != None:
            self["efermi"] = efermi

        if self["natom"] != None and self['energy'] != None:
            self["energy_per_atom"] = self['energy']/self["natom"]
    
    @ResultAbacus.register(force="list[3*natoms], force of the system, if is MD or RELAX calculation, this is the last one")
    def GetForceFromLog(self):
        force = None
        for i in range(len(self.LOG)):
            i = -1*i - 1
            line = self.LOG[i]
            if 'TOTAL-FORCE (eV/Angstrom)' in line:
                #head_pattern = re.compile(r'^\s*atom\s+x\s+y\s+z\s*$')
                value_pattern = re.compile(r'^\s*[A-Z][a-z]?[1-9][0-9]*\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*$')
                j = i
                noforce = False
                while not value_pattern.match(self.LOG[j]):
                    j += 1
                    if j >= i + 10:
                        print("Warning: can not find the first line of force")
                        noforce = True
                        break
                if noforce:
                    break
                
                while value_pattern.match(self.LOG[j]):
                    if force == None:
                        force = []
                    force += [float(ii) for ii in self.LOG[j].split()[1:4]]
                    j += 1
                break
        self['force'] = force
    
    @ResultAbacus.register(stress="list[9], stress of the system, if is MD or RELAX calculation, this is the last one",
                           virial="list[9], virial of the system,  = stress * volume, which is the last one.")
    def GetStessFromLog(self):
        stress = None
        for i in range(len(self.LOG)):
            i = -1*i - 1
            line = self.LOG[i]
            if 'TOTAL-STRESS (KBAR)' in line:
                value_pattern = re.compile(r'^\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*$')
                j = i 
                nostress = False
                # find the first line of stress
                while not value_pattern.match(self.LOG[j]):
                    j += 1
                    if j >= i + 10: # if can not find the first line of force in 10 lines, then stop
                        print("Warning: can not find the first line of force")
                        nostress = True
                        break
                if nostress:
                    break
                
                while value_pattern.match(self.LOG[j]):
                    if stress == None:
                        stress = []
                    stress += [float(ii) for ii in self.LOG[j].split()[:3]]
                    j += 1
                break
        self['stress'] = stress
        if stress != None and self["volume"] != None:
            self['virial'] = [i * self["volume"] * comm.KBAR2EVPERANGSTROM3 for i in stress]
        else:
            self['virial'] = None
        
    @ResultAbacus.register(largest_gradient="list, the largest gradient of each ION step. Unit in eV/Angstrom")
    def GetLargestGradientFromLog(self):
        lg = None
        for line in self.LOG:
            if "Largest gradient is" in line:
                if lg == None:
                    lg = []
                lg.append(float(line.split()[-1]))
        self['largest_gradient'] = lg
        
    '''
    @ResultAbacus.register(band_gap = "band gap of the system")
    def GetBandGapFromLog(self):
        def ErrorReturn(strinfo):
            print("WARNING: %s, skip the catch of band gap info")
            self['band_gap'] = None
            return
        
        nband = self['nbands']  
        if self["ibzk"] != None:    
            nk = self['ibzk']
        elif self["nkstot"] != None:
            nk = self["nkstot"]
        else:
            nk = None
            
        if nband == None or nk == None:
            ErrorReturn("no nbands or ibzk")
            return
                
        band_gap = None
        for i,line in enumerate(self.LOG):
            if 'STATE ENERGY(eV) AND OCCUPATIONS' in line:
                nspin = int(line.split()[-1])
                if nspin not in [1,2]:
                    ErrorReturn("NOT SUPPORT FOR NSPIN=%d now" % nspin)
                    return

                totalcb = None
                totalvb = None
                for ispin in range(nspin):
                    cb = None
                    vb = None
                    fermi = self['efermi']
                    if fermi == None:
                        ErrorReturn("can not get efermi")
                        return

                    for k in range(nk):
                        for m in range(nband):
                            ni = ((nband+2)*nk + 1) * ispin + (nband+2)*k + nspin + m + 1 + i
                            eband = float(self.LOG[ni].split()[1])
                            if eband > fermi:
                                if vb == None or eband < vb:
                                    vb = eband
                                if cb == None or float(self.LOG[ni-1].split()[1]) > cb:
                                    cb = float(self.LOG[ni-1].split()[1])
                                break
                    if totalcb == None or (cb != None and totalcb < cb): totalcb = cb
                    if totalvb == None or (vb != None and totalvb > vb): totalvb = vb
                band_gap = None if totalvb == None or totalcb == None else totalvb-totalcb
                break
                
        self['band_gap'] = band_gap
    '''
        
    @ResultAbacus.register(band_gap = "band gap of the system")
    def GetBandGapFromLog(self):
        def ErrorReturn(strinfo):
            print("WARNING: %s, skip the catch of band gap info" % strinfo)
            self['band_gap'] = None
            return
        
        nband = self['nbands']  
        if self["ibzk"] != None:    
            nk = self['ibzk']
        elif self["nkstot"] != None:
            nk = self["nkstot"]
        else:
            nk = None
            
        if nband == None or nk == None:
            ErrorReturn("no nbands or ibzk")
            return
        nelec = self['nelec']
        if nelec == None:
            ErrorReturn("no nelec")
            return
        occu_band = int(nelec/2)
        if occu_band == 0:
            ErrorReturn("no occu_band")
            return
        if occu_band >= nband:
            ErrorReturn("occu_band >= nband")
            return
                
        band_gap = None
        #print("nelec:",nelec,"occu_band:",occu_band,"nband:",nband,"nk:",nk)
        for i,line in enumerate(self.LOG):
            if 'STATE ENERGY(eV) AND OCCUPATIONS' in line:
                nspin = int(line.split()[-1])
                if nspin not in [1,2]:
                    ErrorReturn("NOT SUPPORT FOR NSPIN=%d now" % nspin)
                    return

                totalcb = None
                totalvb = None
                for ispin in range(nspin):
                    cb = None
                    vb = None
                    for k in range(nk):
                        eband1 = float(self.LOG[((nband+2)*nk + 1) * ispin + (nband+2)*k + nspin + occu_band + i].split()[1])  # the highest occupied band
                        eband2 = float(self.LOG[((nband+2)*nk + 1) * ispin + (nband+2)*k + nspin + occu_band + 1 + i].split()[1])  # the lowest unoccupied band
                        if cb == None or eband1 > cb:
                            cb = eband1
                        if vb == None or eband2 < vb:
                            vb = eband2
                        #print("k:",k,"eband1:",eband1,"eband2:",eband2,"cb:",cb,"vb:",vb)
                    if totalcb == None or (cb != None and totalcb < cb): totalcb = cb
                    if totalvb == None or (vb != None and totalvb > vb): totalvb = vb
                    #print("ispin:",ispin,"totalcb:",totalcb,"totalvb:",totalvb)
                if totalcb == None or totalvb == None:
                    band_gap = None
                else:
                    band_gap = totalvb-totalcb
                    if band_gap < 0: band_gap = 0
                break
                
        self['band_gap'] = band_gap

    @ResultAbacus.register(total_time="the total time of the job",
                           stress_time="the time to do the calculation of stress",
                           force_time = "the time to do the calculation of force",
                           scf_time = "the time to do SCF",
                           scf_time_each_step = "list, the time of each step of SCF",
                           step1_time = "the time of 1st SCF step",
                           scf_steps = "the steps of SCF")
    def GetTimeFromOutput(self):
        # first, check if self.time is already exist
        if self.TIME:
            total_time = self.GetTime("total",None)[0]
            # if PW basis, stress time = Stress_PW/cal_stress, force_time = cal_force_nl
            # if lcao basis, stress time = getForceStress, force_time = None
            stress_time_pw = self.GetTime("Stress_PW","cal_stress")[0]
            stress_time_lcao = self.GetTime("Force_Stress_LCAO","getForceStress")[0]
            force_time_pw = self.GetTime("Forces","cal_force_nl")[0]
            force_time_lcao = None
            if stress_time_pw != None:
                stress_time = stress_time_pw
            else:
                stress_time = stress_time_lcao
                
            if force_time_pw != None:
                force_time = force_time_pw
            else:
                force_time = force_time_lcao
        else:
            stress_time = None
            force_time = None
            total_time = None
            for i,line in enumerate(self.OUTPUT):
                if line[23:28] == 'total' or (len(line.split()) == 5 and line.split()[0] == 'total'):  # old version or new version
                    total_time = float(line.split()[1])
                elif line[23:33] == 'cal_stress' or (len(line.split()) == 6 and line.split()[0] == 'cal_stress'): # old version or new version
                    stress_time = float(line.split()[-5])
                elif line[23:35] == 'cal_force_nl' or (len(line.split()) == 6 and line.split()[0] == 'cal_force_nl'):
                    force_time = float(line.split()[-5])
                elif line[23:37] == 'getForceStress' or (len(line.split()) == 6 and line.split()[0] == 'getForceStress'):
                    stress_time = float(line.split()[-5])
                    force_time = None
                    
        self["total_time"] = total_time
        self['stress_time'] = stress_time
        self['force_time'] = force_time

        scftime = []
        for i,line in enumerate(self.OUTPUT):
            if line[1:5] == 'ITER':
                for j in range(i+1,len(self.OUTPUT)):
                    if self.OUTPUT[j][1:3] in ['CG','DA','GE','GV']:
                        scftime.append(float(self.OUTPUT[j].split()[-1]))
                break
        if len(scftime) > 0:
            import numpy as np
            self['scf_time'] = np.array(scftime).sum()
            self['step1_time'] = scftime[0]
            self['scf_steps'] = len(scftime)
            self['scf_time_each_step'] = scftime

    @ResultAbacus.register(atom_mag="list, the magnization of each atom")
    def GetAtomMag(self):
        mullikenf = os.path.join(os.path.split(self.LOGf)[0],"mulliken.txt")
        if not os.path.isfile(mullikenf):
            self['atom_mag'] = None
            return
        
        atom_mag = []
        with open(mullikenf) as f1: lines = f1.readlines()
        for line in lines:
            if line[:5] == "STEP:":
                atom_mag.append([])
            elif "Total Magnetism on atom" in line:
                if "(" in line and ")" in line:
                    # for nspin = 4
                    sline = line.split("(")[1].split(")")[0].split(",")
                    if len(sline) == 3:
                        atom_mag[-1].append([float(sline[0]),float(sline[1]),float(sline[2])])
                else:
                    atom_mag[-1].append(float(line.split()[-1]))
                    
        self['atom_mag'] = None if len(atom_mag) == 0 else atom_mag
    
    @ResultAbacus.register(drho="[], drho of each scf step",
                           drho_last="drho of the last scf step")
    def GetDrho(self):
        drho = []
        for line in self.LOG:
            if "Density error is" in line:
                drho.append(float(line.split()[-1]))
        
        if len(drho) == 0:
            self['drho'] = None
        else:
            self['drho'] = drho
            self["drho_last"] = drho[-1]
    
    
    @ResultAbacus.register(lattice_constant="unit in angstrom",
                           cell = "[[],[],[]], two-dimension list, unit in Angstrom. If is relax or md, will output the last one",
                           coordinate = "[[],..], two dimension list, is a cartesian type, unit in Angstrom. If is relax or md, will output the last one",
                           element = "list[], a list of the element name of all atoms",
                           label = "list[], a list of atom label of all atoms",
                           element_list = "same as element",
                           atomlabel_list = "same as label")
    def GetCell(self):    
        for line in self.LOG:
            if "lattice constant (Angstrom)" in line:
                self["lattice_constant"] = float(line.split()[-1])
                break
        
        cell = None
        for i in range(len(self.LOG)): 
            iline = -i - 1 
            line = self.LOG[iline]
            if "Lattice vectors: (Cartesian coordinate: in unit of a_0)" in line:
                cell = []
                for k in range(1, 4):
                    cell.append([float(x) * self["lattice_constant"] for x in self.LOG[iline + k].split()[0:3]])
                break
        self['cell'] = cell
        
        coordinate = None
        for i in range(len(self.LOG)): 
            iline = -i - 1 
            line = self.LOG[iline]
            if len(line.split()) >= 2 and line.split()[1] == "COORDINATES":  
                coordinate = []
                if line.split()[0] == "DIRECT":
                    for k in range(2, 2 + self["natom"]):
                        coordinate.append([float(x) for x in self.LOG[iline + k].split()[1:4]])
                    import numpy as np
                    coordinate = np.array(coordinate).dot(np.array(self["cell"])).tolist()
                elif line.split()[0] == "CARTESIAN":
                    for k in range(2, 2 + self["natom"]):
                        coordinate.append([float(x) for x in self.LOG[iline + k].split()[1:4]])
                else:
                    print("Unrecongnized coordinate type: %s" % (line))   
                break
        self['coordinate'] = coordinate   
              
        element_list = None
        atomlabel_list = None
        for line in self.LOG:
            if "atom label =" in line:
                label = line.split()[-1]
                element = label
                while element[-1].isdigit():
                    element = element[:-1]
            elif "number of atom for this type =" in line:
                if element_list == None: element_list = []
                if atomlabel_list == None: atomlabel_list = []
                natom = int(line.split()[-1])
                atomlabel_list += natom * [label]
                element_list += natom * [element]
                if len(atomlabel_list) >= self['natom']:
                    break
        self["element"] = element_list
        self["label"] = atomlabel_list
        self['element_list'] = element_list
        self['atomlabel_list'] = atomlabel_list
            
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


class AbacusRelax(ResultAbacus):
    
    @ResultAbacus.register(relax_converge="if the relax is converged")
    def GetRelaxConverge(self):
        #need read self.LOG
        if self.LOG:
            for i in range(len(self.LOG)):
                line = self.LOG[-i-1]
                if "Relaxation is converged!" in line:
                    self["relax_converge"] = True
                    return
                elif "Relaxation is not converged yet!" in line:
                    self["relax_converge"] = False
                    return
                elif "Ion relaxation is not converged yet" in line or \
                    "Lattice relaxation is not converged yet" in line:
                    self["relax_converge"] = False
                    return
                elif "Lattice relaxation is converged!" in line or \
                    "Ion relaxation is converged!" in line:
                    self["relax_converge"] = True
                    return
        self["relax_converge"] = None
    
    @ResultAbacus.register(relax_steps= "the total ION steps")
    def GetRelaxSteps(self):
        #need read self.LOG
        if self.LOG:
            for i in range(len(self.LOG)):
                line = self.LOG[-i-1]
                if "ALGORITHM --------------- ION=" in line:
                    index_ben = line.index("ION=") + 4
                    index_end = line.index("ELEC")
                    self["relax_steps"] = int(line[index_ben:index_end])
                    return
        self["relax_steps"] = None
