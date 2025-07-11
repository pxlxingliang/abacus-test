import os,sys,glob,re
from ..resultVasp import ResultVasp
from .. import comm
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element

xfmlt = comm.XmlFindMultiLayerText
xfml = comm.XmlFindMultiLayer

class Vasp(ResultVasp):
    
    @ResultVasp.register(version="the vasp version",
                         ncore = "mpi cores")
    def GetGeneralInfo(self):
        if self.XMLROOT != None:
            self['version'] = self.XMLROOT.find("./generator/i[@name='version']").text
        
        for line in self.OUTCAR:
            if "running on" in line and "total cores" in line:
                self['ncore'] = int(line.split()[2])
    
    @ResultVasp.register(normal_end="if the job is normal ending")
    def GetNormalEnd(self):
        if len(self.OUTCAR) == 0:
            self['normal_end'] = None
            return
        elif len(self.OUTCAR) > 0 and "Voluntary context switches:" in self.OUTCAR[-1]:
            self['normal_end'] = True
            return
        else:
            self['normal_end'] = False

            print("Job is not normal ending!!! The latest 10 lines is:")
            if len(self.OUTCAR) < 10:
                print(''.join(self.OUTCAR))
            else:
                print(''.join(self.OUTCAR[-10:]))

    @ResultVasp.register(kpt="list, the K POINTS setting",
                         nkstot = "total K point number",
                         ibzk = "irreducible K point number")
    def GetKPT(self):
        if self.XMLROOT != None:
            tree = "./kpoints/generation/v[@name='divisions']"
            kpts = self.XMLROOT.find(tree)
            if kpts != None:
                nkstot = 1
                kpt = []
                for i in kpts.text.split():
                    nkstot *= int(i)
                    kpt.append(int(i))
                self['kpt'] = kpt
                self['nkstot'] = nkstot
            else:
                self['kpt'] = None
                self['nkstot'] = None

            tree = "./kpoints/varray[@name='kpointlist']/v"
            ibzk = self.XMLROOT.findall(tree)
            if ibzk != None:
                self['ibzk'] = len(ibzk)
            else:
                self['ibzk'] = None
        else:
            ibzk = None
            kpt = None
            nkstot = None
            for i,line in enumerate(self.OUTCAR):
                if "k-points in BZ" in line:
                    ibzk = int(line.split()[3])
                elif " generate k-points for:" in line:
                    kpt = [int(j) for j in self.OUTCAR[i].split()[-3:]]
                elif "k-points in reciprocal lattice and weights" in line:
                    nkstot = 0
                    for j in range(i+1,len(self.OUTCAR)):
                        if self.OUTCAR[j].strip() == "":
                            break
                        nkstot += 1
                    break
                    
            self["ibzk"] = ibzk        
            self['kpt'] = kpt
            self['nkstot'] = nkstot        
                
    @ResultVasp.register(nbands="number of bands",
                         nelec = "total electron number",
                         spin = "the spin number",
                         encut = "eV, the energy cutoff",
                         ismear = "the smearing method",
                         sigma = "the SIGMA setting of smearing, in eV",
                         nelm = "value of NELM, the setted maximum SCF steps",
                         natom = "total atom number",
                         volume = "volume(A^3). if is relax or md, will return the volume of last ION step")
    def GetInputSetting(self):
        if self.XMLROOT != None:
            tree = ".//separator[@name='electronic']/i[@name='NBANDS']"
            self['nbands'] = comm.iint(comm.XmlGetText(self.XMLROOT.find(tree)))

            tree = ".//separator[@name='electronic']/i[@name='NELECT']"
            self['nelec'] = comm.ifloat(comm.XmlGetText(self.XMLROOT.find(tree)))

            tree = ".//separator[@name='electronic spin']/i[@name='ISPIN']"
            self['spin'] = comm.iint(comm.XmlGetText(self.XMLROOT.find(tree)))
            
            tree = ".//separator[@name='electronic']/i[@name='ENMAX']"
            self['encut'] = comm.ifloat(comm.XmlGetText(self.XMLROOT.find(tree)))

            tree = ".//separator[@name='electronic smearing']/i[@name='ISMEAR']"
            self['ismear'] = comm.iint(comm.XmlGetText(self.XMLROOT.find(tree)))

            tree = ".//separator[@name='electronic smearing']/i[@name='SIGMA']"
            self['sigma'] = comm.ifloat(comm.XmlGetText(self.XMLROOT.find(tree)))
            
            tree = ".//separator[@name='electronic convergence']/i[@name='NELM']"
            self['nelm'] = comm.iint(comm.XmlGetText(self.XMLROOT.find(tree)))

            tree = "./atominfo/atoms"
            self['natom'] = comm.iint(comm.XmlGetText(self.XMLROOT.find(tree)))
            
            tree = "./structure/crystal/i[@name='volume']"
            self['volume'] = comm.ifloat(comm.XmlGetText(self.XMLROOT.findall(tree),idx=-1))
            
        else:
            volume = None
            for line in self.OUTCAR:
                sline = line.split()
                if "number of bands    NBANDS" in line:
                    self['nbands'] = int(sline[-1])
                elif "number of ions     NIONS =" in line:
                    self['natom'] = int(sline[-1])
                elif "ISPIN  =" in line:
                    self['spin'] = int(sline[2])
                elif "ENCUT  =  " in line:
                    self['encut'] = float(sline[4])
                elif "NELECT =" in line:
                    self["nelec"] = float(sline[2])
                elif "ISMEAR =" in line:
                    self['ismear'] = sline[2].rstrip(";")
                    self["sigma"] = float(sline[5])
                elif "NELM   =" in line:
                    self['nelm'] = int(sline[2][:-1])
                elif "volume of cell" in line:
                    volume = float(sline[-1])
            self["volume"] = volume

    @ResultVasp.register(ldautype = "value of LDAUTYPE, the type of plus U",
                         ldaul = "list, value of LDAUL, the l-quantum number of each element",
                         ldauu = "list, value of LDAUU, the U setting",
                         ldauj = "list, value of LDAUJ, the J setting")
    def GetLdaUSetting(self):
        for i,line in enumerate(self.OUTCAR):
            sline = line.split()
            if "LDA+U is selected," in line:
                self['ldautype'] = int(sline[-1])
                self['ldaul'] = [int(j) for j in self.OUTCAR[i+1].split("=")[-1].split()]
                self['ldauu'] = [float(j) for j in self.OUTCAR[i+2].split("=")[-1].split()]
                self['ldauj'] = [float(j) for j in self.OUTCAR[i+3].split("=")[-1].split()]
                break

    @ResultVasp.register(denergy_last = 'The energy difference in SCF, eV',
                         denergy = "The last energy difference in SCF, eV")
    def GetDE(self):
        des = []
        for i in range(len(self.OSZICAR)):
            line = self.OSZICAR[i]
            if line.strip() == "": continue
            if line[3] == ":" and len(line.split()) == 7:
                de = float(line.split()[3])
                des.append(de)
        if len(des) > 0:
            self['denergy_last'] = des[-1]
            self['denergy'] = des
        else:
            self['denergy_last'] = None
            self['denergy'] = None
   
    @ResultVasp.register(scf_steps = 'the steps of SCF, if is relax or md job, only last ION step is read',
                         converge = "if the SCF is converged. If scf_steps is smaller than NELM, will be converged, else is not converged")
    def GetSCFInfo(self):
        for i in range(len(self.OUTCAR)):
            j = -1*i - 1
            line = self.OUTCAR[j]
            if 'Iteration' in line:
                self['scf_steps'] = int(line.split('(')[1].split(')')[0])
                if self['scf_steps'] < self['nelm'] and self['efermi'] is not None:
                    # in some case the SCF is not converged but terminated, so we need to check the fermi energy
                    self['converge'] = True
                else:
                    self['converge'] = False
                break
    
    @ResultVasp.register(energy = 'eV,the total energy, if is relax or md job, will return the energy of last ION step',
                         energy_per_atom = 'eV, the energy divided by natom, if is relax or md job, will return the energy of last ION step')
    def GetEnergy(self):
        for i in range(len(self.OUTCAR)):
            line = self.OUTCAR[-i-1]
            if "energy  without entropy=" in line:
                self['energy'] = float(line.split()[-4])
                if self['natom'] != None:
                    self['energy_per_atom'] = self['energy'] / self['natom']
                else:
                    self['energy_per_atom'] = None
                return
    
    @ResultVasp.register(force = 'list, eV/angstrom, the force of all atoms, [atom1x,atom1y,atom1z,atom2x,atom2y,atom2z...]',
                         stress = 'list, kBar, the stress, [xx,xy,xz,yx,yy,yz,zx,zy,zz]',
                         virial='list, eV, the virial, [xx,xy,xz,yx,yy,yz,zx,zy,zz]',
                         pressure="kBar, the pressure, 1/3*trace(stress)",
                         cell = 'list[list], Angstrom, the vector of cell. If is RELAX or MD job, will output the last cell.',
                         )
    def GetForceStress(self):
        force = None
        stress = None
        virial = None
        cell = None
        for i,line in enumerate(self.OUTCAR):
            if 'TOTAL-FORCE (eV/Angst)' in line:
                j = i+2
                force = []
                while self.OUTCAR[j][:3] != " --":
                    force += [float(k) for k in self.OUTCAR[j].split()[3:6]]
                    j += 1
            elif '  in kB' in line:
                s = [float(i) for i in line.split()[2:8]]
                stress = [s[0],s[3],s[5],s[3],s[1],s[4],s[5],s[4],s[2]]
                v = [float(i) for i in self.OUTCAR[i-1].split()[1:7]]
                virial = [v[0],v[3],v[5],v[3],v[1],v[4],v[5],v[4],v[2]]
            elif "VOLUME and BASIS-vectors are now :" in line:
                cell = []
                for j in range(3):
                    cell.append([float(k) for k in self.OUTCAR[i+j+5].split()[0:3]])
        self['force'] = force
        self['stress']  = stress
        self['virial'] = virial
        self["cell"] = cell
        if stress:
            self['pressure'] = 1.0/3.0*(stress[0]+stress[4]+stress[8])
        else:
            self['pressure'] = None
        
    @ResultVasp.register(total_time = 'Total CPU time (s)',
                         scf_time = 'the total SCF times, s',
                         stress_time = 'the time of calculating stress')                         
    def GetTimeInfo(self):
        stresst = None
        scft = 0
        for line in self.OUTCAR:
            if 'STRESS:  cpu time' in line:
                stresst = float(line.split()[-1])
            elif 'LOOP:  cpu time' in line:
                scft += float(line.split()[-1])
            elif 'Total CPU time used (sec):' in line:
                self['total_time'] = float(line.split()[-1])

        if stresst != None:
            self['stress_time'] = stresst
        if scft > 0:
            self['scf_time'] = scft

    @ResultVasp.register(total_mag = 'total magnization',
                         absolute_mag="absolute magnetism, the summation of the magnetic moment of all atoms",
                         atom_mag = 'list, the magnization of each atom')
    def GetMagInfo(self):
        getatommag = False
        atommagsx = []
        atommagsy = []
        atommagsz = []
        
        for i in range(len(self.OUTCAR)):
            i = -i-1
            line = self.OUTCAR[i]
            if line[:19] == " number of electron":
                self['total_mag'] = float(line.split()[-1])
                break
            elif line[:18] == ' magnetization (x)':
                j = i + 4
                atommag = []
                while self.OUTCAR[j][:3] != "---":
                    atommag.append(float(self.OUTCAR[j].split()[-1]))
                    j += 1
                atommagsx.append(atommag)
            elif line[:18] == ' magnetization (y)':
                j = i + 4
                atommag = []
                while self.OUTCAR[j][:3] != "---":
                    atommag.append(float(self.OUTCAR[j].split()[-1]))
                    j += 1
                atommagsy.append(atommag)
            elif line[:18] == ' magnetization (z)':
                j = i + 4
                atommag = []
                while self.OUTCAR[j][:3] != "---":
                    atommag.append(float(self.OUTCAR[j].split()[-1]))
                    j += 1
                atommagsz.append(atommag)
        
        if len(atommagsx) > 0 and len(atommagsy) > 0 and len(atommagsz) > 0:
            assert len(atommagsx) == len(atommagsy) == len(atommagsz), "The atom magnetism in x, y, z direction should have the same length."
            atomx = atommagsx[0]
            atomy = atommagsy[0]
            atomz = atommagsz[0]
            atommag = []
            abs_mag = 0
            for i in range(len(atomx)):
                atommag.append([atomx[i],atomy[i],atomz[i]])
                abs_mag += (atomx[i]**2 + atomy[i]**2 + atomz[i]**2)**0.5
            self["atom_mag"] = atommag
            self["absolute_mag"] = abs_mag
        elif len(atommagsx) > 0:
            self["atom_mag"] = atommagsx[0]
            self["absolute_mag"] = sum([abs(i) for i in atommagsx[0]])
        else:
            self["atom_mag"] = None
            self["absolute_mag"]
                
    
    @ResultVasp.register(atom_name = 'list, the element name of each atom' ,
                         atom_type = 'list, the element name of each atomtype',
                         element_list = 'list, the element name of all atoms',
                         element = "list[], a list of the element name of all atoms",
                           label = "list[], a list of atom label of all atoms",
                           atomlabel_list = "same as label",
                         efermi     = 'the fermi energy, eV')
    def GetXMLInfo(self):
        if self.XMLROOT != None:
            atom_name = comm.XmlGetText(self.XMLROOT.findall("./atominfo/array[@name='atoms']/set/rc/c[1]"))
            if isinstance(atom_name,list):
                atom_name = [i.strip() for i in atom_name]
            self["atom_name"] = atom_name
            
            atom_type = comm.XmlGetText(self.XMLROOT.findall("./atominfo/array[@name='atomtypes']/set/rc/c[2]"))
            if isinstance(atom_type,list):
                atom_type = [i.strip() for i in atom_type]
            self["atom_type"] = atom_type
            
            self['efermi'] = comm.XmlGetText(self.XMLROOT.findall("./calculation/dos/i[@name='efermi'][last()]"),func=float,idx = -1)
            element = []
            for i in self.XMLROOT.findall("./atominfo/array[@name='atoms']/set/rc/c[1]"):
                element.append(i.text)
            self["element"] = element
            self["label"] = element
            self["element_list"] = element
            self["atomlabel_list"] = element
        else:
            # read from OUTCAR
            atom_name = None
            atom_type = None
            efermi = None
            for i in range(len(self.OUTCAR)):
                if "VRHFIN" in self.OUTCAR[i]:
                    if atom_name == None:
                        atom_name = []
                    atom_name.append(self.OUTCAR[i].split(":")[0].split("=")[1].strip())
                elif "E-fermi" in self.OUTCAR[i]:
                    efermi = float(self.OUTCAR[i].split()[2])
            if atom_name != None:
                atom_type = []
                for i in atom_name:
                    if i not in atom_type:
                        atom_type.append(i)
            self['atom_type'] = atom_type
            self['atom_name'] = atom_name
            self["element"] = atom_name
            self["label"] = atom_name
            self["element_list"] = atom_name
            self["atomlabel_list"] = atom_name
            self['efermi'] = efermi

    @ResultVasp.register(band = '[[[]]], list with three dimension spin*kpoint*band')
    def GetBandInfo(self):
        if self.XMLROOT != None:
            band = []
            eigen = self.XMLROOT.findall('./calculation/eigenvalues/array')
            if eigen == None:
                self['band'] = None
            else:
                array = eigen[-1]
                for spin in array.find('set').findall('set'):
                    band.append([])
                    for kpoint in spin.findall('set'):
                        band[-1].append([])
                        for iband in kpoint.findall('r'):
                            band[-1][-1].append(float(iband.text.split()[0]))
                self['band'] = band
        else:
            # try to read band from OUTCAR
            band = []
            nspin = self["spin"]
            nk = self["nkstot"]
            nband = self["nbands"]
            if nspin == None or nk == None:
                self['band'] = None
                return
            print("nspin=%d, nk=%d, nband=%d" % (nspin,nk,nband))
            for i in range(len(self.OUTCAR)):
                if "spin component" in self.OUTCAR[i]:
                    band.append([])
                    for j in range(nk):
                        nk_start = i + 1 + (nband+3)*j + 3
                        nk_end = nk_start + nband
                        band[-1].append([float(k.split()[1]) for k in self.OUTCAR[nk_start:nk_end]])
                        
                    nspin -= 1
                    if nspin == 0:
                        break
            if len(band) == 0:
                self['band'] = None
            else:
                self['band'] = band
    
    @ResultVasp.register(band_gap = 'eV, the band gap')
    def GetBandGap(self):
        if self['band'] == None or self['efermi'] == None:
            self['band_gap'] = None
        else:
            vb = None
            cb = None
            fermi = self['efermi']
            self["band_gap"] = comm.cal_band_gap(self['band'],fermi)
    
    @ResultVasp.register(band_plot="Plot the band structure. Return the file name of the plot.")
    def PlotBandFromLog(self):  
        band = self['band']  
        efermi = self['efermi']
        if band == None:
            print("no band, and skip the plot of band")
            self['band_plot'] = None
            return
        band_plot = os.path.join(self.PATH,"band.png")
        comm.plot_band(band, band_plot, efermi)
        self['band_plot'] = band_plot     
    
    @ResultVasp.register(point_group = 'point group',
                         point_group_in_space_group = "point group in space group")
    def GetPointGroup(self):
        pg = None
        pgsg = None
        for i in self.OUTCAR:
            if pg != None and pgsg != None:
                break
            
            if "The dynamic configuration has the point symmetry" in i:
                pg = i.strip()[48:-1].strip()
            elif " The point group associated with its full space group is" in i:
                pgsg = i.strip()[56:-1].strip()
        self['point_group'] = pg
        self['point_group_in_space_group'] = pgsg
        
    @ResultVasp.register(relax_steps = 'the steps of ION',
                         relax_converge = 'if the relax is converged',)
    def GetRelaxInfo(self):
        if self.XMLROOT != None:
            self["relax_steps"] = len(self.XMLROOT.findall("calculation"))
            
            #max_ion = comm.iint(self.XMLROOT.findtext("./parameters/ionic[@name='NSW']",default=64))
            max_ion = comm.iint(self.XMLROOT.findtext("./parameters/separator[@name='ionic']/i[@name='NSW']",default=64))
            if max_ion != None and self["relax_steps"] != None:
                if max_ion > self["relax_steps"]:
                    self["relax_converge"] = True
                else:
                    self["relax_converge"] = False
            else:
                self["relax_converge"] = None
        elif self.OUTCAR:
            nsw = 1
            for i in self.OUTCAR:
                if "   NSW    =" in i:
                    nsw = int(i.split()[2])
            relax_steps = 0
            for i in self.OUTCAR[::-1]:
                if "--------------------------------------- Iteration" in i:
                    relax_steps = int(i[49:56])
                    break
            if relax_steps == 0 or relax_steps == nsw:
                self["relax_converge"] = False
            else:
                self["relax_converge"] = True
            self["relax_steps"] = relax_steps