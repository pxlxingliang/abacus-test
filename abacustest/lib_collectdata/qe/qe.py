import os,sys,glob,re
import traceback
from ..resultQe import ResultQe
from .. import comm
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element
import numpy as np

xfmlt = comm.XmlFindMultiLayerText
xfml = comm.XmlFindMultiLayer

class Qe(ResultQe):
    @ResultQe.register(normal_end="if the job is normal ending")
    def GetNormalEnd(self):
        if self.XMLROOT:
            exit_status = self.XMLROOT.find('exit_status')
            if exit_status != None:
                if exit_status.text == '0':
                    self["normal_end"] = True
                else:
                    self["normal_end"] = False
                return   
            self["normal_end"] = None
        elif self.OUTPUT:
            if len(self.OUTPUT) > 1 and "JOB DONE." in self.OUTPUT[-2]:
                self["normal_end"] = True
            else:
                self["normal_end"] = False
        else:
            self["normal_end"] = None
            
    @ResultQe.register(version="the version of QE",
                       ncore = "the mpi cores")
    def GetGeneralInfo(self):
        if self.XMLROOT != None:
            self['version'] = self.XMLROOT.find("general_info/creator").attrib['VERSION']
            self['ncore'] = self.XMLROOT.find("parallel_info/nprocs").text
        elif self.OUTPUT:
            for iline in self.OUTPUT:
                if "Program PWSCF" in iline:
                    self['version'] = iline.split()[2].strip().lstrip("v.")
                elif "Parallel version (MPI)" in iline:
                    try:
                        self['ncore'] = int(iline.split()[5].strip())
                    except:
                        traceback.print_exc()
                        self['ncore'] = None
                    break
        else:
            self['version'] = None
            self['ncore'] = None
    
    @ResultQe.register(natom="total atom number",
                       ecutwfc="energy cutoff in Ry",
                       ks_solver="the diagonalization method",
                       mixing_type="the charge mixing method in SCF",
                       mixing_beta="the coefficient beta of charge mixing in SCF",
                       scf_thr="the convergence threshold of SCF, unit in Ry")
    def GetInputParam(self):
        if self.XMLROOT != None:

            input_param = self.XMLROOT.find('input')        

            self['natom'] = comm.iint(input_param.find('atomic_structure').attrib['nat'])


            ecutwfc = comm.ifloat(xfmlt(input_param,['basis','ecutwfc']))
            if ecutwfc != None: ecutwfc *= 2
            self["ecutwfc"] = ecutwfc

            self['ks_solver'] = xfmlt(input_param,['electron_control','diagonalization'])
            self['mixing_type'] = xfmlt(input_param,['electron_control','mixing_mode'])
            self['mixing_beta'] = comm.ifloat(xfmlt(input_param,['electron_control','mixing_beta']))

            scf_thr = comm.ifloat(xfmlt(input_param,['electron_control','conv_thr']))
            if scf_thr != None: scf_thr *= 2
            self['scf_thr'] = scf_thr
        elif self.OUTPUT:
            ks_solver = None
            for iline in self.OUTPUT:
                if "number of atoms/cell" in iline:
                    self["natom"] = int(iline.split()[4])
                elif "kinetic-energy cutoff" in iline:
                    self["ecutwfc"] = float(iline.split()[4])
                elif "diagonalization" in iline:
                    ks_solver = iline.split()[0].lower()
                elif "mixing beta" in iline:
                    self["mixing_beta"] = float(iline.split()[3])
                elif "number of iterations used =" in iline:
                    self["mixing_type"] = iline.split()[6]
                elif "convergence threshold" in iline:
                    self["scf_thr"] = float(iline.split()[4])
            self["ks_solver"] = ks_solver
        else:
            self["natom"] = None
            self["ecutwfc"] = None
            self["ks_solver"] = None
            self["mixing_type"] = None
            self["mixing_beta"] = None
            self["scf_thr"] = None

    @ResultQe.register(smearing_method="smearing method",
                       smearing_sigma = "smearing sigma in Ry")
    def GetSmearing(self):
        if self.XMLROOT != None:
            input_param = self.XMLROOT.find('input')

            smearing_method = xfmlt(input_param,['bands','occupations'])
            if smearing_method == 'fixed':
                smearing_sigma = 0.0
            elif smearing_method == 'smearing':
                smearing_method = xfmlt(input_param,['bands','smearing'])
                smearing_sigma = comm.ifloat(xfml(input_param,['bands','smearing']).attrib.get('degauss',0.0))
                if smearing_sigma != None:
                    smearing_sigma *= 2.0

            self["smearing_method"] = smearing_method
            self["smearing_sigma"] = smearing_sigma
        else:
            self["smearing_method"] = None
            self["smearing_sigma"] = None
            
    
    @ResultQe.register(converge="if SCF is converged",
                       scf_steps="the steps of SCF",
                       energy="total energy, unit in eV",
                       energies = "list, the total energy of each ION step",
                       nbands="band number",
                       ibzk="irreducible K point number",
                       nkstot="the total k points",
                       kpt="list, the K POINTS setting",
                       force="list, the force of all atoms, [atom1x,atom1y,atom1z,atom2x,atom2y,atom2z...]. Unit in eV/Angstrom",
                       stress="list, the stress, [xx,xy,xz,yx,yy,yz,zx,zy,zz]. Unit in kbar.",
                       virial="list, the virial, [xx,xy,xz,yx,yy,yz,zx,zy,zz]. Unit in eV.",
                       pressure="the pressure, unit in kbar",
                       forces = "list, the force of each ION step",
                       stresses = "list, the stress of each ION step",
                       virials = "list, the virial of each ION step",
                       pressures = "list, the pressure of each ION step",
                       cell = "list, the cell, [a1,a2,a3,b1,b2,b3,c1,c2,c3]. Unit in Angstrom",
                       volume = "the volume of cell, unit in Angstrom^3",
                       coord = "list, the coordinate of all atoms, [atom1x,atom1y,atom1z,atom2x,atom2y,atom2z...]. Unit in Angstrom",
                       label = "the label of each atom",
                       atomlabel_list = "same as label",
                       element = "the element of each atom. Because QE do not output the element, we actually output the label.",
                       element_list = "same as element.",
                       total_mag="total magnization",
                       absolute_mag="total absolute magnization",
                       atom_mag = "list, the magnization of each atom",
                       nelec="total electron number",
                       energy_per_atom="total energy divided by natom, unit in eV",)
    def GetOutputParam(self):
        output = None
        if self.XMLROOT != None:
            output = self.XMLROOT.findall('output')
            if len(output) == 0:
                print("ERROR: can not find the output in pwscf.xml")
                return
            output = output[-1]

            element = None
            self['converge'] = comm.ibool(xfmlt(output,['convergence_info','scf_conv','convergence_achieved']))
            self['scf_steps'] = comm.iint(xfmlt(output,['convergence_info','scf_conv','n_scf_steps']))
            self['ibzk'] = comm.iint(xfmlt(output,['band_structure','nks']))
            self['nelec'] = comm.ifloat(xfmlt(output,['band_structure','nelec']))
            self['nbands'] = comm.iint(xfmlt(output,['band_structure','nbnd']))

            total_energy = comm.ifloat(xfmlt(output,['total_energy','etot']))
            if total_energy != None: total_energy *= comm.HARTREE2EV
            self['energy'] = total_energy

            mp = xfml(output,['band_structure','starting_k_points','monkhorst_pack'])
            if isinstance(mp, Element):
                nk1 = comm.iint(mp.get('nk1'))
                nk2 = comm.iint(mp.get('nk2'))
                nk3 = comm.iint(mp.get('nk3'))
                if nk1 and nk2 and nk3: 
                    self['nkstot'] = nk1 * nk2 * nk3
                    self['kpt'] = [nk1, nk2, nk3]

            # structure
            structure = output.findall('atomic_structure')
            if len(structure) == 0:
                structure = None
            else:
                structure = structure[-1]
            volume = None
            if structure != None:
                cell_a = xfmlt(structure,['cell','a1'])
                cell_b = xfmlt(structure,['cell','a2'])
                cell_c = xfmlt(structure,['cell','a3'])
                if cell_a and cell_b and cell_c:
                    self['cell'] = cell = [[float(i)* comm.BOHR2A for i in j.split()]  for j in [cell_a,cell_b,cell_c]]
                    # calculate the volume by cell
                    
                    self["volume"] = volume = np.linalg.det(cell)
                coord = []
                label = []
                for icoord in structure.findall('atomic_positions/atom'):
                    coord.append([float(i) * comm.BOHR2A for i in icoord.text.split()])
                    label.append(icoord.attrib['name'])
                if len(coord) == self['natom']:
                    self['coord'] = coord
                    self["label"] = label
                    element = label
                else:
                    self['coord'] = None
                    self["label"] = None
                    element = None
                    print("ERROR: the length of coord is not equal to natom")
                self["atomlabel_list"] = self["label"]
            else:
                self['cell'] = None
                self['coord'] = None

            if output.find('forces') != None:
                self['force'] = [float(i) * comm.HARTREE2EV / comm.BOHR2A for i in output.find('forces').text.split()]
            if output.find('stress') != None:
                # the read in stress is in Hartree/Bohr^3, convert to kbar
                # 1 kbar = 1e8 Pa = 1e8 N/m^2 = 1e8 J/m^3 = 1e8 * 2.2937126583579E17 Hartree/m^3 = 2.2937126583579E25 * 5.29177E-11**3 Hartree/Bohr^3 = 3.398927420868445E-6 Hartree/Bohr^3
                self['stress'] = [ float(i)/comm.KBAR2HARTREEPERBOHR3 for i in output.find('stress').text.split()]
                self["pressure"] = (self['stress'][0] + self['stress'][4] + self['stress'][8]) / 3.0
                # calculate the virial, unit in eV
                if volume != None:
                    self['virial'] = [ float(i) * volume * comm.HARTREE2EV / comm.BOHR2A ** 3   for i in output.find('stress').text.split()]

            tot_mag = comm.ifloat(xfmlt(output,['magnetization','total']))
            if tot_mag == None and output.find('magnetization/total_vec') != None:
                tot_mag = [float(i) for i in output.find('magnetization/total_vec').text.split()]
            self['total_mag'] = tot_mag if tot_mag != None else None
            self['absolute_mag'] = comm.ifloat(xfmlt(output,['magnetization','absolute']))
            atom_mags = output.findall('magnetization/Scalar_Site_Magnetic_Moments/SiteMagnetization')
            if atom_mags:
                atom_mag = []
                element = []
                for i in atom_mags:
                    atom_mag.append(float(i.text))
                    element.append(i.attrib['species'].strip())
                self['atom_mag'] = atom_mag
            else:
                self['atom_mag'] = None
                
            if self['energy'] == None or self['natom'] == None:
                self['energy_per_atom'] = None
            else:
                self['energy_per_atom'] = self['energy'] / self['natom']
            self["element"] = element
            self["element_list"] = element
            
            allsteps = self.XMLROOT.findall('step')
            energies = []
            forces = []
            stresses = []
            pressures = []
            virials = []
            for istep in allsteps:
                energies.append(comm.ifloat(xfmlt(istep,['total_energy','etot'])) * comm.HARTREE2EV)
                force = istep.find('forces')
                if force != None:
                    forces.append([float(i) * comm.HARTREE2EV / comm.BOHR2A for i in force.text.split()])
                stress = istep.find('stress')
                if stress != None:
                    stresses.append([float(i)/comm.KBAR2HARTREEPERBOHR3 for i in stress.text.split()])
                    pressures.append((stresses[-1][0] + stresses[-1][4] + stresses[-1][8]) / 3.0)
                    if volume != None:
                        virials.append([float(i) * volume * comm.HARTREE2EV / comm.BOHR2A ** 3   for i in stress.text.split()])
            self["energies"] = energies if energies else None
            self["forces"] = forces if forces else None
            self["stresses"] = stresses if stresses else None
            self["pressures"] = pressures if pressures else None
            self["virials"] = virials if virials else None
            
        elif self.OUTPUT:
            converge = None
            ibzk = None
            nbands = None
            nelec = None
            energies = []
            stresses = []
            pressures = []
            forces = []
            natom = None
            cell = None
            coord = None
            label = None
            element = None
            tot_mag = None
            abs_mag = None
            atom_mag = None
            scf_steps = None
            for idx,iline in enumerate(self.OUTPUT):
                if "convergence has been achieved" in iline:
                    converge = True
                elif "convergence NOT achieved" in iline:
                    converge = False
                elif "atomic species   valence    mass     pseudopotential" in iline:
                    j = idx + 1
                    element = {}
                    while self.OUTPUT[j].strip() !="":
                        species = self.OUTPUT[j].split()[0]
                        element[species] = self.OUTPUT[j].split("(")[0].split()[-1].strip()
                        j += 1
                elif "number of k points" in iline:
                    ibzk = int(iline.split()[4])
                elif "total energy" in iline:
                    energies.append(float(iline.split()[-2]) * comm.RY2EV)
                elif "number of Kohn-Sham states" in iline:
                    nbands = int(iline.split()[-1])
                elif "number of electrons" in iline:
                    nelec = float(iline.split()[-1])
                elif "number of atoms/cell" in iline:
                    natom = int(iline.split()[-1])
                elif "total magnetization" in iline:
                    sline = iline.split()
                    if len(sline) == 8:
                        tot_mag = [float(imag) for imag in sline[3:6]]
                    else:
                        tot_mag = float(sline[3])
                elif "absolute magnetization" in iline:
                    abs_mag = float(iline.split()[3])
                elif "iteration #" in iline:
                    if scf_steps == None:
                        scf_steps = 1
                    else:
                        scf_steps += 1
                elif "magn=" in iline:
                    if atom_mag == None:
                        atom_mag = []
                    atom_mag.append(float(iline.split()[-1]))
                elif "total stress" in iline:
                    stress = []
                    for i in range(3):
                        stress.extend([float(j) for j in self.OUTPUT[idx+i+1].split()[3:6]])
                    pressure = float(iline.split()[-1])
                    stresses.append(stress)
                    pressures.append(pressure)
                elif "Forces acting on atoms" in iline:
                    force = []
                    j = idx + 2
                    while len(self.OUTPUT[j].split()) == 8:
                        force.extend([float(i) * comm.RY2EV / comm.BOHR2A for i in self.OUTPUT[j].split()[-3:]])
                        j += 1
                    forces.append(force)
                elif "CELL_PARAMETERS" in iline:
                    cell = []
                    lat0 = float(iline.split("=")[1].strip()[:-1])
                    for i in range(3):
                        cell.append([float(j) * lat0 for j in self.OUTPUT[idx+i+1].split()])
                elif "ATOMIC_POSITIONS" in iline:
                    coord = []
                    label = []
                    if natom == None:
                        print("ERROR: can not find natom when catch coord")
                        continue
                    for i in range(natom):
                        icoord = self.OUTPUT[self.OUTPUT.index(iline)+i+1].split()
                        label.append(icoord[0])
                        coord.append([float(j) for j in icoord[1:]])
                elif "site n.     atom                  positions (alat units)" in iline:
                    label = []
                    j = idx + 1
                    while self.OUTPUT[j].strip() != "":
                        label.append(self.OUTPUT[j].split()[1])
                        j += 1
                        
            self["converge"] = converge
            self["energy"] = None if len(energies) == 0 else energies[-1]
            self["energies"] = None if len(energies) == 0 else energies
            self["nbands"] = nbands
            self["ibzk"] = ibzk
            self["nelec"] = nelec
            self["stress"] = None if len(stresses) == 0 else stresses[-1]
            self["stresses"] = None if len(stresses) == 0 else stresses
            self["pressure"] = None if len(pressures) == 0 else pressures[-1]
            self["pressures"] = None if len(pressures) == 0 else pressures
            self["force"] = None if len(forces) == 0 else forces[-1]
            self["forces"] = None if len(forces) == 0 else forces
            self["natom"] = natom
            self["cell"] = cell
            self["coord"] = coord
            self["label"] = label
            self["atomlabel_list"] = label
            self["element"] = None if label ==None or element == None else [element[i] for i in label]
            self["element_list"] = label
            self["total_mag"] = tot_mag
            self["absolute_mag"] = abs_mag
            self["atom_mag"] = atom_mag
            self["scf_steps"] = scf_steps
            self["nkstot"] = None
            self["kpt"] = None
            if energies and natom != None:
                self["energy_per_atom"] = energies[-1] / natom
            else:
                self["energy_per_atom"] = None
            
            if cell != None:
                self["volume"] = abs(np.linalg.det(cell))
            else:
                self["volume"] = None
            
            # calc the virial
            if stresses and self["volume"] != None:
                virials = []
                for stress in stresses:
                    virials.append([i * self["volume"] * comm.KBAR2EVPERANGSTROM3 for i in stress])
                self["virial"] = virials[-1]
                self["virials"] = virials
            else:
                self["virial"] = None
                self["virials"] = None
        else:
            self["converge"] = None
            self["scf_steps"] = None
            self["energy"] = None
            self["energies"] = None
            self["nbands"] = None
            self["ibzk"] = None
            self["nkstot"] = None
            self["kpt"] = None
            self["force"] = None
            self["forces"] = None
            self["stress"] = None
            self["stresses"] = None
            self["virial"] = None
            self["virials"] = None
            self["pressure"] = None
            self["pressures"] = None
            self["cell"] = None
            self["volume"] = None
            self["coord"] = None
            self["label"] = None
            self["atomlabel_list"] = None
            self["element"] = None
            self["element_list"] = None
            self["total_mag"] = None
            self["absolute_mag"] = None
            self["atom_mag"] = None
            self["nelec"] = None
            self["energy_per_atom"] = None
            
    @ResultQe.register(band="list, the band, [[spin1],[spin2],...], unit in eV",
                       efermi="the Fermi energy, unit in eV",
                       band_gap="the band gap, unit in eV")
    def GetBand(self):  
        band = None
        efermi = None
        band_gap = None
        if self.XMLROOT != None:
            band_all = []
            for ik in self.XMLROOT.findall('output/band_structure/ks_energies/eigenvalues'):
                band_all.append([float(i) * comm.HARTREE2EV for i in ik.text.split()])
            nbnd_up = comm.iint(xfmlt(self.XMLROOT,['output','band_structure','nbnd_up']))
            if nbnd_up != None:
                band = [[],[]]
                for ik in band_all:
                    band[0].append(ik[:nbnd_up])
                    band[1].append(ik[nbnd_up:])
            else:
                band = [band_all]     
                
            efermi = comm.ifloat(xfmlt(self.XMLROOT,['output','band_structure','fermi_energy']))
            if efermi != None:
                efermi *= comm.HARTREE2EV
        elif self.OUTPUT:
            band_all = []
            for idx,iline in enumerate(self.OUTPUT):
                if "the Fermi energy is" in iline:
                    efermi = float(iline.split()[-2])
                    break
                elif " k =" in iline:
                    j = idx + 2
                    band_all.append([])
                    while self.OUTPUT[j].strip() != "":
                        band_all[-1].extend([float(i) for i in self.OUTPUT[j].split()])
                        j += 1
            ibzk = self["ibzk"]
            if ibzk != None and len(band_all) == 2*ibzk:
                band = [band_all[:ibzk],band_all[ibzk:]]
            elif band:
                band = [band_all]
            
        if band != None and efermi != None:
            band_gap = comm.cal_band_gap(band,efermi)
        self["band"] = band
        self["efermi"] = efermi
        self["band_gap"] = band_gap      
    
    @ResultQe.register(drho="[], drho of each scf step. Unit in Ry",
                       drho_last="drho of the last scf step. Unit in Ry",
                       denergy="[], denergy of each scf step. Unit in eV",
                       denergy_last="denergy of the last scf step. Unit in eV")
    def Getdrhode(self):
        drho = None
        drho_last = None
        denergy = None
        denergy_last = None
        energy = None
        if self.OUTPUT:
            for iline in self.OUTPUT:
                if "total energy              =" in iline:
                    if energy == None:
                        energy = []
                    energy.append(iline.split()[-2])
                elif "estimated scf accuracy    <" in iline:
                    if drho == None:
                        drho = []
                    drho.append(float(iline.split()[-2]))
        
        if drho != None:
            drho_last = drho[-1]
        if energy != None:
            e_acc = 10**(-1 * len(energy[0].split(".")[1]))
            de = e_acc* comm.RY2EV
            denergy = [0.0]
            for i in range(1,len(energy)):
                e1 = float(energy[i])
                e2 = float(energy[i-1])
                if not (abs(e1 - e2) < e_acc):
                    de = (e1 - e2)*comm.RY2EV
                denergy.append(de)
            denergy_last = denergy[-1]
            
            
        self["drho"] = drho
        self["drho_last"] = drho_last
        self["denergy"] = denergy
        self["denergy_last"] = denergy_last
                    
            
    
    @ResultQe.register(total_time="the total running time (WALL)",
                       force_time="the time to do the calculations of force (WALL)",
                       stress_time="the time to do the calculation of stress (WALL)")
    def GetTimeParam(self):
        if self.XMLROOT != None:

            timing = self.XMLROOT.find('timing_info')

            self['total_time'] = comm.ifloat(xfmlt(timing,['total','wall']))

            all_time = timing.findall('partial')
            for itime in all_time:
                if itime.get('label') == 'forces':
                    self['force_time'] = comm.ifloat(itime.find('wall').text)
                elif itime.get('label') == 'stress':
                    self['stress_time'] = comm.ifloat(itime.find('wall').text)
        elif self.OUTPUT:
            # read the time from the output file
            # begin to read from the end of the file
            force_time = None
            stress_time = None
            total_time = None
            
            for iline in self.OUTPUT[::-1]:
                if "Writing all to output data dir" in iline:
                    break
                if "PWSCF        :" in iline:
                    # may have case "PWSCF        :      7h43m CPU      8h 0m WALL"
                    total_time = comm.strtime2sec(iline.split("CPU")[-1].split("WALL")[0])
                elif "force" in iline:
                    force_time = comm.strtime2sec(iline.split()[-2])
                elif "stress" in iline:
                    stress_time = comm.strtime2sec(iline.split()[-2])
            self["total_time"] = total_time
            self["force_time"] = force_time
            self["stress_time"] = stress_time
        else:
            self["total_time"] = None
            self["force_time"] = None
            self["stress_time"] = None

    @ResultQe.register(scf_time="the SCF time (WALL)",
                       scf_time_each_step="list, the SCF time of each step (WALL)",
                       scf_time_step1 = "the time of the first step (WALL)")
    def GetSCFTime(self):
        if self.OUTPUT:
            start_time = None # the time before SCF
            scf_times = [] # the time of each SCF step
            scf_loop = False
            for i,line in enumerate(self.OUTPUT):
                if "Self-consistent Calculation" in line:
                    for j in range(i-1,-1,-1):
                        if "total cpu time spent up to now is" in self.OUTPUT[j]:
                            start_time = float(self.OUTPUT[j].split()[-2])
                            break
                    scf_loop = True
                elif "End of self-consistent calculation" in line:
                    break
                elif "total cpu time spent up to now is" in line and scf_loop and start_time != None:
                    scf_times.append(float(line.split()[-2]) - start_time)
                    start_time = float(line.split()[-2])
            self["scf_time"] = None if not scf_times else sum(scf_times)
            self["scf_time_each_step"] = None if not scf_times else scf_times
            self["scf_time_step1"] = None if not scf_times else scf_times[0]
        else:
            self["scf_time"] = None
            self["scf_time_each_step"] = None
            self["scf_time_step1"] = None
            
                    
                    
    
    @ResultQe.register(relax_converge="if the relax is converged")
    def GetRelaxConverge(self):
        if self.OUTPUT:
            for line in self.OUTPUT[::-1]:
                if "The maximum number of steps has been reached." in line:
                    self["relax_converge"] = False
                    return
                elif "bfgs converged in" in line:
                    self["relax_converge"] = True
                    return
            self["relax_converge"] = None
            
        elif self.XMLROOT != None:
            output = self.XMLROOT.find('output')
            rlx_conv = comm.ibool(xfmlt(output,['convergence_info','opt_conv','convergence_achieved']))
            # in QE, after the rlx is converged, it will do another SCF, but in some case the scf is not done, and this value will be false
            # we need to check if the ION steps is less than input/control_variables/nstep, if yes, we should set this value to True
            if not rlx_conv:
                nstep = comm.iint(xfmlt(self.XMLROOT,["input",'control_variables','nstep']))
                if nstep == None:
                    nstep = 50
                if self["relax_steps"] != None and self["relax_steps"] < nstep:
                    rlx_conv = True
            
            self["relax_converge"] = rlx_conv
        else:
            self["relax_converge"] = None

    @ResultQe.register(relax_steps="the total ION steps")
    def GetRelaxStepsFromXml(self):
        if self.XMLROOT != None:
            self["relax_steps"] = len(self.XMLROOT.findall('step'))
        elif self.OUTPUT:
            for iline in self.OUTPUT[::-1]:
                if "number of bfgs steps" in iline:
                    self["relax_steps"] = int(iline.split()[-1])
                    return
                
        self["relax_steps"] = None
        
        

