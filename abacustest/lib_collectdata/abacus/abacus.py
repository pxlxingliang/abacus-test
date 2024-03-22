import os,sys,glob,re,traceback
from ..resultAbacus import ResultAbacus
from .. import comm

class Abacus(ResultAbacus):
    
    @ResultAbacus.register(version="the version of ABACUS")
    def GetVersion(self):
        version = "unknown"
        commit = "unknown"
        if self.JSON:
            version = self.JSON.get("general_info",{}).get("version","unknown")
            commit = self.JSON.get("general_info",{}).get("commit","unknown")
        else:
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
                           nelec_dict = "dict of electron number of each species",
                           fft_grid = "fft grid for charge/potential",
                           point_group="point group",
                           point_group_in_space_group="point group in space group")
    def GetLogParam(self):       
        natom = 0
        nelec = 0
        elec_dict = {}
        point_group = None
        point_group_in_space_group = None
        for i,line in enumerate(self.LOG):
            if "NBANDS =" in line:
                self['nbands'] = int(line.split()[2])
            elif 'nkstot =' in line:
                self['nkstot'] = int(line.split()[-1])
            elif 'nkstot_ibz =' in line:
                self['ibzk'] = int(line.split()[-1])
            elif "            electron number of element" in line:
                elec_dict[line.split()[4]] = float(line.split()[-1])
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
        else:
            self["nelec"] = None
        
        if elec_dict:
            self["nelec_dict"] = elec_dict
        else:
            self["nelec_dict"] = None

    @ResultAbacus.register(converge="if the SCF is converged",
                           total_mag="total magnetism (Bohr mag/cell)",
                           absolute_mag="absolute magnetism (Bohr mag/cell)",
                           energy = "the total energy (eV)",
                           volume = "the volume of cell, in A^3",
                           efermi = "the fermi energy (eV). If has set nupdown, this will be a list of two values. The first is up, the second is down.",
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
                if 'E_Fermi_up' in line:
                    if efermi == None:
                        efermi = [None,None]
                    efermi[0] = float(line.split()[-1])
                elif 'E_Fermi_dw' in line:
                    if efermi == None:
                        efermi = [None,None]
                    efermi[1] = float(line.split()[-1])
                else:
                    efermi = float(line.split()[-1])

        self["energy"] = energy
        self["converge"] = converge
        self["volume"] = volume
        
        self['total_mag'] = total_mag
        self['absolute_mag'] = absolute_mag
        if efermi != None and isinstance(efermi,list):
            if abs(efermi[0] - efermi[1]) < 1e-6:
                efermi = efermi[0]
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
                           virial="list[9], virial of the system,  = stress * volume, which is the last one.",
                           pressure="the pressure of the system, unit in kbar.")
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
            self["pressure"] = (stress[0] + stress[4] + stress[8]) / 3.0
        else:
            self['virial'] = None
            self["pressure"] = None
        
    @ResultAbacus.register(largest_gradient="list, the largest gradient of each ION step. Unit in eV/Angstrom")
    def GetLargestGradientFromLog(self):
        lg = None
        for line in self.LOG:
            if "Largest gradient is" in line:
                if lg == None:
                    lg = []
                lg.append(float(line.split()[-1]))
        self['largest_gradient'] = lg

    
    @ResultAbacus.register(band = "Band of system. Dimension is [nspin,nk,nband].",
                           band_weight = "Band weight of system. Dimension is [nspin,nk,nband].")
    def GetBandFromLog(self): 
        nband = self['nbands']  
        if self["ibzk"] != None:    
            nk = self['ibzk']
        elif self["nkstot"] != None:
            nk = self["nkstot"]
        else:
            nk = None
            
        if nband == None or nk == None:
            print("no nbands or ibzk, and skip the catch of band info")
            return
        
        band = None
        band_weight = None
        for i,line in enumerate(self.LOG):
            if 'STATE ENERGY(eV) AND OCCUPATIONS' in line:
                nspin = int(line.split()[-1])
                band = []
                band_weight = []
                for ispin in range(nspin):
                    band.append([])
                    band_weight.append([])
                    for k in range(nk):
                        band[-1].append([])
                        band_weight[-1].append([])
                        for m in range(nband):
                            ni = ((nband+2)*nk + 1) * ispin + (nband+2)*k + nspin + m + 1 + i
                            eband = float(self.LOG[ni].split()[1])
                            wband = float(self.LOG[ni].split()[2])
                            band[-1][-1].append(eband)
                            band_weight[-1][-1].append(wband)
                break
        self['band'] = band
        self['band_weight'] = band_weight
        
    
    @ResultAbacus.register(band_gap = "band gap of the system")
    def GetBandGapFromLog(self):
        def ErrorReturn(strinfo):
            print("WARNING: %s, skip the catch of band gap info")
            self['band_gap'] = None
            return
        
        band = self['band']
        band_weight = self['band_weight']
        if band == None or band_weight == None:
            ErrorReturn("no band")
            return
        
        efermi = self['efermi']
        if efermi == None:
            ErrorReturn("no efermi")
            return
        
        self['band_gap'] = comm.cal_band_gap(band,efermi)
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
    '''
    
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
                    if self.OUTPUT[j][1:3] in ['CG','DA','GE','GV','BP']:
                        scftime.append(float(self.OUTPUT[j].split()[-1]))
                break
        if len(scftime) > 0:
            import numpy as np
            self['scf_time'] = np.array(scftime).sum()
            self['step1_time'] = scftime[0]
            self['scf_steps'] = len(scftime)
            self['scf_time_each_step'] = scftime

    @ResultAbacus.register(atom_mags="list of list, the magnization of each atom of each ion step.",
                           atom_mag = "list, the magnization of each atom. Only the last ION step.",
                           atom_elec="list of list of each atom. Each atom list is a list of each orbital, and each orbital is a list of each spin",
                           atom_orb_elec="list of list of each atom. Each atom list is a list of each orbital, and each orbital is a list of each spin",
                           )
    def GetAtomMag(self):
        mullikenf = os.path.join(os.path.split(self.LOGf)[0],"mulliken.txt")
        if not os.path.isfile(mullikenf):
            self['atom_mag'] = None
            self["atom_mags"] = None
            self["atom_electron"] = None
            return
        
        atom_mag = []
        atom_electron = []
        with open(mullikenf) as f1: lines = f1.readlines()
        for idx,line in enumerate(lines):
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
            elif "Zeta of" in line:
                two_spin = True if "Spin 2" in line else False
                atom_electron.append([])
                j = idx + 1
                while j < len(lines) and lines[j].strip() != "":
                    if "Zeta of" in lines[j]:
                        break
                    if "sum over m+zeta" in lines[j]:
                        atom_electron[-1].append([float(lines[j].split()[3])])
                        if two_spin:
                            atom_electron[-1][-1].append(float(lines[j].split()[4]))
                    j += 1
        
        if len(atom_mag) == 0:
            self['atom_mags'] = None
            self["atom_mag"] = None
        else:
            self['atom_mags'] = atom_mag
            self["atom_mag"] = atom_mag[-1]
        
        if not atom_electron:
            self["atom_elec"] = None
            self["atom_orb_elec"] = None
        else:
            self["atom_orb_elec"] = atom_electron
            atom_elec = []
            for i in atom_electron:
                t = 0
                for j in i:
                    t += sum(j)
                atom_elec.append(t)
            self["atom_elec"] = atom_elec
            
    
    @ResultAbacus.register(atom_mag_u="list of a dict, the magnization of each atom calculated by occupation number. Only the last SCF step.",
                           atom_elec_u = "list of a dict with keys are atom index, atom label, and electron of U orbital.")
    def GetAtomMag(self):
        atom_mag_u = None
        atom_elec_u = None
        if self.LOG:
            u_block = []
            start_line = end_line = None
            for i in  range(len(self.LOG)):
                i = -1*i - 1
                if "//=========================L(S)DA+U===========================//" in self.LOG[i]:
                    start_line = i + 1
                    break
                elif "//=======================================================//" in self.LOG[i]:
                    end_line = i
            print(start_line,end_line)
            if None not in [start_line, end_line]:
                u_block = self.LOG[start_line:end_line]
            else:
                print("Can not find the L(S)DA+U block")
                self['atom_mag_u'] = None
                self["atom_elec_u"] = None
                return
            
            labels = self["label"]
            atom_mag_u = []
            atom_elec_u = []
            for i,line in enumerate(u_block):
                if "atoms" in line:
                    atom_index = int(line.split()[-1])
                    atom_label = None if not labels else labels[atom_index]
                    atom_mag_u.append({"index":atom_index,"label":atom_label,"mag":None})
                    atom_elec_u.append({"index":atom_index,"label":atom_label,"elec":None})
                elif "eigenvalues" in line:
                    if atom_elec_u[-1]["elec"] == None:
                        atom_elec_u[-1]["elec"] = []
                    atom_elec_u[-1]["elec"].append(float(u_block[i+1].split()[-1]))
                elif "atomic mag:" in line:
                    atom_mag_u[-1]["mag"] = float(line.split()[-1])
        self['atom_mag_u'] = atom_mag_u
        self['atom_elec_u'] = atom_elec_u
    
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
    
    @ResultAbacus.register(denergy="[], denergy of each scf step",
                           denergy_last="denergy of the last scf step")
    def GetDenergy(self):
        denergy = None
        
        if self.OUTPUT:
            for i,line in enumerate(self.OUTPUT):
                if "ITER" in line and "ETOT(eV)" in line and "EDIFF(eV)" in line and "DRHO" in line and "TIME(s)" in line:
                    denergy = []
                    ncol = len(line.split())
                    ediff_idx = line.split().index("EDIFF(eV)")
                    j = i + 1
                    while j < len(self.OUTPUT):
                        if "----------------------------" in self.OUTPUT[j]:
                            break
                        if self.OUTPUT[j][1:3] in ['CG','DA','GE','GV','BP'] and len(self.OUTPUT[j].split()) == ncol:
                            denergy.append(float(self.OUTPUT[j].split()[ediff_idx]))
                        j += 1
                    break
        self["denergy"] = denergy
        if denergy != None:
            self["denergy_last"] = denergy[-1]
        else:
            self["denergy_last"] = None

    @ResultAbacus.register(lattice_constant="unit in angstrom",
                           cell = "[[],[],[]], two-dimension list, unit in Angstrom. If is relax or md, will output the last one.",
                           cell_init = "[[],[],[]], two-dimension list, unit in Angstrom. The initial cell",
                           coordinate = "[[],..], two dimension list, is a cartesian type, unit in Angstrom. If is relax or md, will output the last one",
                           coordinate_init = "[[],..], two dimension list, is a cartesian type, unit in Angstrom. The initial coordinate",
                           element = "list[], a list of the element name of all atoms",
                           label = "list[], a list of atom label of all atoms",
                           element_list = "same as element",
                           atomlabel_list = "same as label")
    def GetCell(self): 
        lc = 1   
        for line in self.LOG:
            if "lattice constant (Angstrom)" in line:
                self["lattice_constant"] = lc = float(line.split()[-1])
                break
        
        cell = None
        for i in range(len(self.LOG)): 
            iline = -i - 1 
            line = self.LOG[iline]
            if "Lattice vectors: (Cartesian coordinate: in unit of a_0)" in line:
                cell = []
                for k in range(1, 4):
                    cell.append([float(x) * lc for x in self.LOG[iline + k].split()[0:3]])
                break
        self['cell'] = cell
        
        cell_init = None
        for i in range(len(self.LOG)): 
            iline = i 
            line = self.LOG[iline]
            if "Lattice vectors: (Cartesian coordinate: in unit of a_0)" in line:
                cell_init = []
                for k in range(1, 4):
                    cell_init.append([float(x) * lc for x in self.LOG[iline + k].split()[0:3]])
                break
        self['cell_init'] = cell_init
        
        coordinate = None
        coordinate_init = None
        natom = self['natom']
        if cell != None and natom != None:
            try:
                for i in range(len(self.LOG)): 
                    iline = -i - 1 
                    line = self.LOG[iline]
                    if len(line.split()) >= 2 and line.split()[1] == "COORDINATES":  
                        coordinate = []
                        if line.split()[0] == "DIRECT":
                            for k in range(2, 2 + natom):
                                coordinate.append([float(x) for x in self.LOG[iline + k].split()[1:4]])
                            import numpy as np
                            coordinate = np.array(coordinate).dot(np.array(cell)).tolist()
                        elif line.split()[0] == "CARTESIAN":
                            for k in range(2, 2 + natom):
                                coordinate.append([float(x) for x in self.LOG[iline + k].split()[1:4]])
                        else:
                            print("Unrecongnized coordinate type: %s" % (line))   
                        break
                self['coordinate'] = coordinate  
            except:
                traceback.print_exc()
                self['coordinate'] = None 
            
            try:
                for i in range(len(self.LOG)): 
                    iline = i
                    line = self.LOG[iline]
                    if len(line.split()) >= 2 and line.split()[1] == "COORDINATES":  
                        coordinate_init = []
                        if line.split()[0] == "DIRECT":
                            for k in range(2, 2 + natom):
                                coordinate_init.append([float(x) for x in self.LOG[iline + k].split()[1:4]])
                            import numpy as np
                            coordinate_init = np.array(coordinate_init).dot(np.array(cell)).tolist()
                        elif line.split()[0] == "CARTESIAN":
                            for k in range(2, 2 + natom):
                                coordinate_init.append([float(x) for x in self.LOG[iline + k].split()[1:4]])
                        else:
                            print("Unrecongnized coordinate type: %s" % (line))   
                        break
                self['coordinate_init'] = coordinate_init 
            except:
                traceback.print_exc()
                self['coordinate_init'] = None   
        else:
            print("No cell or natom, skip the catch of coordinate info")
            self['coordinate'] = None
            self['coordinate_init'] = None
              
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
    
    @ResultAbacus.register(pdos="a dict, keys are 'energy' and 'orbitals', and 'orbitals' is a list of dict which is (index,species,l,m,z,data), dimension of data is nspin*ne",
                           )
    def GetPDOS(self): 
        pdos_file = os.path.join(self.PATH,f"OUT.{self.SUFFIX}","PDOS")
        print("PDOS file path:",pdos_file)
        pdos = None
        if os.path.isfile(pdos_file):
            import xml.etree.ElementTree as ET
            tree = ET.parse(pdos_file)
            root = tree.getroot()
            nspin = int(root.find('nspin').text)
            energy = [float(i) for i in root.find('energy_values').text.split()]

            ne = len(energy)

            all_orbitals = []
            for iorb in root.findall('orbital'):
                data = [[] for i in range(nspin)]
                for i in iorb.find('data').text.split("\n"):
                    for j,jj in enumerate(i.split()):
                        if j == 1:
                            data[j].append(-1*float(jj))
                        else:
                            data[j].append(float(jj))
                if len(data[0]) != ne:
                    print("WARNING: PDOS len(data[0]) != ne")      
                all_orbitals.append({
                    "index": int(iorb.get('index')),
                    "atom_index": int(iorb.get('atom_index')),
                    "species": iorb.get('species'),
                    "l": int(iorb.get('l')),
                    "m": int(iorb.get('m')),
                    "z": int(iorb.get('z')),
                    "data": data
                })
            pdos = {"energy":energy,"orbitals":all_orbitals}
        self["pdos"] = pdos


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
                elif " STEP OF ION RELAXATION : " in line:
                    self["relax_steps"] = int(line.split()[-1])
                    return
        self["relax_steps"] = None


class AbacusDeltaSpin(ResultAbacus):
    
    @ResultAbacus.register(ds_lambda_step="a list of DeltaSpin converge step in each SCF step",
                           ds_lambda_rms="a list of DeltaSpin RMS in each SCF step")
    def GetDSLambdaStep(self):
        '''
Step (Outer -- Inner) =  16 -- 1           RMS = 1.057e-07
Step (Outer -- Inner) =  16 -- 2           RMS = 1.355e-07
Step (Outer -- Inner) =  16 -- 3           RMS = 1.176e-07
Step (Outer -- Inner) =  16 -- 4           RMS = 9.760e-08
        '''
        lambda_step = None 
        lambda_rms = None
        if self.OUTPUT:
            scf_step = []
            lambda_step = []
            lambda_rms = []
            for i in self.OUTPUT:
                if "Step (Outer -- Inner) =" in i:
                    ss = int(i.split()[5])
                    lambdas = int(i.split()[7])
                    rms = float(i.split()[10])
                    if ss not in scf_step:
                        scf_step.append(ss)
                        lambda_step.append(lambdas)
                        lambda_rms.append(rms)
                    else:
                        lambda_step[-1] = lambdas
                        lambda_rms[-1] = rms
        self["ds_lambda_step"] = lambda_step
        self["ds_lambda_rms"] = lambda_rms
            