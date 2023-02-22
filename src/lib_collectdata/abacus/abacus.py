import os,sys,glob,re
from ..resultAbacus import ResultAbacus

class Abacus(ResultAbacus):
    
    @ResultAbacus.register(version="the version of ABACUS")
    def GetVersion(self):
        if len(self.LOG) > 0:
            for line in self.LOG:
                if "WELCOME TO ABACUS" in line:
                    version = line.split()[-1]
                    if version[0].lower() != 'v':
                        print("Unknow version of '%s'" % version)
                    self['version'] = version
                    
    @ResultAbacus.register(ncore="the mpi cores")
    def GetNcore(self):
        for line in self.LOG:
            if "DSIZE =" in line:
                self['ncore'] = int(line.split()[-1])
                return
    
    @ResultAbacus.register(normal_end="if the job is nromal ending")
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
        
        readinput = False
        INPUT = {}
        for i,iline in enumerate(self.INPUT):
            if iline.strip() == 'INPUT_PARAMETERS':
                readinput = True
            elif iline.strip() == '' or iline.strip()[0] in ['#']:
                continue
            elif readinput:
                sline = re.split('[ \t]',iline.split("#")[0],maxsplit=1)
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
                           converge="if the SCF is converged",
                           total_mag="total magnetism (Bohr mag/cell)",
                           absolute_mag="absolute magnetism (Bohr mag/cell)",
                           nkstot = "total K point number",
                           ibzk = "irreducible K point number",
                           natom ="total atom number",
                           nelec = "total electron number",
                           energy = "the total energy (eV)",
                           volume = "the volume of cell, in A^3",
                           fft_grid = "fft grid for charge/potential",
                           efermi = "the fermi energy (eV)",
                           energy_per_atom="the total energy divided by natom, (eV)")
    def GetLogParam(self):       
        natom = 0
        nelec = 0
        total_mag = None
        absolute_mag = None
        efermi = None
        for i,line in enumerate(self.LOG):
            if "NBANDS =" in line:
                self['nbands'] = int(line.split()[2])
            elif 'charge density convergence is achieved' in line:
                self['converge'] = True
            elif 'convergence has NOT been achieved!' in line or\
                'convergence has not been achieved' in line:
                self['converge'] = False
            elif 'total magnetism (Bohr mag/cell)' in line:
                total_mag = float(line.split()[-1])
            elif 'absolute magnetism' in line:
                absolute_mag = float(line.split()[-1])
            elif 'nkstot =' in line:
                self['nkstot'] = int(line.split()[-1])
            elif 'nkstot_ibz =' in line:
                self['ibzk'] = int(line.split()[-1])
            elif 'number of atom for this type =' in line:
                natom += int(line.split()[-1])
            elif 'total electron number of element' in line:
                nelec += float(line.split()[-1])
            elif "!FINAL_ETOT_IS" in line:
                self['energy'] = float(line.split()[1])
            elif "Volume (A^3) =" in line:
                self['volume'] = float(line.split()[-1])
            elif "[fft grid for charge/potential] =" in line:
                self['fft_grid'] = [float(i.strip()) for i in line.split('=')[1].split(',')]
            elif 'E_Fermi' in line:
                if 'E_Fermi_dw' in line:
                    if efermi == None:
                        efermi = float(line.split()[-1])
                        continue
                    if float(line.split()[-1]) > efermi:
                        efermi = float(line.split()[-1])
                else:
                    efermi = float(line.split()[-1])

        if natom > 0:
            self["natom"] = natom 
        if nelec > 0:
            self["nelec"] = nelec 

        if total_mag != None:
            self['total_mag'] = total_mag
        if absolute_mag != None:
            self['absolute_mag'] = absolute_mag
        if efermi != None:
            self["efermi"] = efermi

        if self["natom"] != None and self['energy'] != None:
            self["energy_per_atom"] = self['energy']/self["natom"]

    @ResultAbacus.register(stress="list[9], stress of the system, if is MD or RELAX calculation, this is the last one",
                           force="list[3*natoms], force of the system, if is MD or RELAX calculation, this is the last one")
    def GetForceStessFromLog(self):
        getforce = getstress = False
        stress = force = None
        for i in range(len(self.LOG)):
            if getforce and getstress:
                break
            i = -1*i - 1
            line = self.LOG[i]
            if not getstress and 'TOTAL-STRESS (KBAR)' in line:
                j = i + 4
                stress = []
                for k in range(3):
                    for m in self.LOG[j+k].split():
                        stress.append(float(m))
                getstress = True
            elif not getforce and 'TOTAL-FORCE (eV/Angstrom)' in line:
                j = i + 5
                force = []
                while len(self.LOG[j].split()) == 4:
                    for k in self.LOG[j].split()[1:4]:
                        force.append(float(k))
                    j += 1
                getforce = True

        self['stress'] = stress
        self['force'] = force
    
    @ResultAbacus.register(band_gap = "band gap of the system")
    def GetBandGapFromLog(self):
        def ErrorReturn(strinfo):
            print("WARNING: %s, skip the catch of band gap info")
            self['band_gap'] = None
            return
        
        band_gap = None
        for i,line in enumerate(self.LOG):
            if 'STATE ENERGY(eV) AND OCCUPATIONS' in line:
                nband = self['nbands']
                nk = self['ibzk']
                if nband == None or nk == None:
                    ErrorReturn("no nbands or ibzk")

                nspin = int(line.split()[-1])
                if nspin not in [1,2]:
                    ErrorReturn("NOT SUPPORT FOR NSPIN=%d now" % nspin)

                totalcb = -999999
                totalvb = 999999
                for ispin in range(nspin):
                    cb = -999999
                    vb = 999999
                    fermi = self['efermi']
                    if fermi == None:
                        ErrorReturn("can not get efermi")

                    for k in range(nk):
                        for m in range(nband):
                            ni = ((nband+2)*nk + 1) * ispin + (nband+2)*k + nspin + m + 1 + i
                            eband = float(self.LOG[ni].split()[1])
                            if eband > fermi:
                                if eband < vb:
                                    vb = eband
                                if float(self.LOG[ni-1].split()[1]) > cb:
                                    cb = float(self.LOG[ni-1].split()[1])
                                break
                    if totalcb < cb: totalcb = cb
                    if totalvb > vb: totalvb = vb
                band_gap = totalvb-totalcb
                
        self['band_gap'] = band_gap

    @ResultAbacus.register(total_time="the total time of the job",
                           stress_time="the time to do the calculation of stress",
                           force_time = "the time to do the calculation of force",
                           scf_time = "the time to do SCF",
                           scf_time_each_step = "list, the time of each step of SCF",
                           step1_time = "the time of 1st SCF step",
                           scf_steps = "the steps of SCF")
    def GetTimeFromOutput(self):
        scftime = []
        stress_time = None
        force_time = None
        for i,line in enumerate(self.OUTPUT):
            if line[1:5] == 'ITER':
                for j in range(i+1,len(self.OUTPUT)):
                    if self.OUTPUT[j][1:3] in ['CG','DA','GE','GV']:
                        scftime.append(float(self.OUTPUT[j].split()[-1]))
            elif line[23:28] == 'total':
                self['total_time'] = float(line.split()[1])
            elif line[23:33] == 'cal_stress':
                stress_time = float(line.split()[-5])
            elif line[23:35] == 'cal_force_nl':
                force_time = float(line.split()[-5])
            elif line[23:37] == 'getForceStress':
                stress_time = float(line.split()[-5])
                force_time = None
        
        self['stress_time'] = stress_time
        self['force_time'] = force_time
        
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
            if "Total Magnetism on atom" in line:
                atom_mag.append(float(line.split()[-1]))
        self['atom_mag'] = None if len(atom_mag) == 0 else atom_mag
    
    @ResultAbacus.register(drho="[], drho of each scf step")
    def GetDrho(self):
        drho = []
        for line in self.LOG:
            if "Density error is" in line:
                drho.append(float(line.split()[-1]))
        
        if len(drho) == 0:
            self['drho'] = None
        else:
            self['drho'] = drho
    
    
    @ResultAbacus.register(lattice_constant="unit in angstrom",
                           cell = "[[],[],[]], two-dimension list, unit in Angstrom. If is relax or md, will out the last one",
                           coordinate = "[[],..], two dimension list, is a cartessian type, unit in angstrom. If is relax or md, will out the last one",
                           element_list = "list[], a list of the element name of all atoms",
                           atomlabel_list = "list[], a list of atom label of all atoms")
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
        self['element_list'] = element_list
        self['atomlabel_list'] = atomlabel_list
            
