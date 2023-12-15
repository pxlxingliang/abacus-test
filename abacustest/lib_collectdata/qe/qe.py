import os,sys,glob,re
from ..resultQe import ResultQe
from .. import comm
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element

xfmlt = comm.XmlFindMultiLayerText
xfml = comm.XmlFindMultiLayer

class Qe(ResultQe):
    @ResultQe.register(normal_end="if the job is normal ending")
    def GetNormalEnd(self):
        if self.XMLROOT == None:
            self["normal_end"] = None
            return
        
        exit_status = self.XMLROOT.find('exit_status')
        if exit_status != None:
            if exit_status.text == '0':
                self["normal_end"] = True
            else:
                self["normal_end"] = False
            return   
        self["normal_end"] = None
            
    @ResultQe.register(version="the version of QE",
                       ncore = "the mpi cores")
    def GetGeneralInfo(self):
        if self.XMLROOT != None:
            self['version'] = self.XMLROOT.find("general_info/creator").attrib['VERSION']
            self['ncore'] = self.XMLROOT.find("parallel_info/nprocs").text
    
    @ResultQe.register(natom="total atom number",
                       ecutwfc="energy cutoff in Ry",
                       ks_solver="the diagonalization method",
                       mixing_type="the charge mixing method in SCF",
                       mixing_beta="the coefficient beta of charge mixing in SCF",
                       scf_thr="the convergence threshold of SCF, unit in Ry")
    def GetInputParamFromXml(self):
        if self.XMLROOT == None:
            return

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

    @ResultQe.register(smearing_method="smearing method",
                       smearing_sigma = "smearing sigma in Ry")
    def GetSmearingFromXml(self):
        if self.XMLROOT == None:
            return

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
    
    @ResultQe.register(converge="if SCF is converged",
                       scf_steps="the steps of SCF",
                       energy="total energy, unit in eV",
                       nbands="band number",
                       ibzk="irreducible K point number",
                       nkstot="the total k points",
                       kpt="list, the K POINTS setting",
                       force="list, the force of all atoms, [atom1x,atom1y,atom1z,atom2x,atom2y,atom2z...]. Unit in eV/Angstrom",
                       stress="list, the stress, [xx,xy,xz,yx,yy,yz,zx,zy,zz]. Unit in kbar.",
                       virial="list, the virial, [xx,xy,xz,yx,yy,yz,zx,zy,zz]. Unit in eV.",
                       cell = "list, the cell, [a1,a2,a3,b1,b2,b3,c1,c2,c3]. Unit in Angstrom",
                       volume = "the volume of cell, unit in Angstrom^3",
                       coord = "list, the coordinate of all atoms, [atom1x,atom1y,atom1z,atom2x,atom2y,atom2z...]. Unit in Angstrom",
                       label = "the label of each atom",
                       total_mag="total magnization",
                       absolute_mag="total absolute magnization",
                       nelec="total electron number",
                       energy_per_atom="total energy divided by natom, unit in eV",)
    def GetOutputParamFromXml(self):
        if self.XMLROOT == None:
            return

        output = self.XMLROOT.findall('output')
        if len(output) == 0:
            return
        output = output[-1]

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
                import numpy as np
                self["volume"] = volume = np.linalg.det(cell)
            coord = []
            label = []
            for icoord in structure.findall('atomic_positions/atom'):
                coord.append([float(i) * comm.BOHR2A for i in icoord.text.split()])
                label.append(icoord.attrib['name'])
            if len(coord) == self['natom']:
                self['coord'] = coord
                self["label"] = label
            else:
                self['coord'] = None
                self["label"] = None
                print("ERROR: the length of coord is not equal to natom")
        else:
            self['cell'] = None
            self['coord'] = None
        
        if output.find('forces') != None:
            self['force'] = [float(i) * comm.HARTREE2EV / comm.BOHR2A for i in output.find('forces').text.split()]
        if output.find('stress') != None:
            # the read in stress is in Hartree/Bohr^3, convert to kbar
            # 1 kbar = 1e8 Pa = 1e8 N/m^2 = 1e8 J/m^3 = 1e8 * 2.2937126583579E17 Hartree/m^3 = 2.2937126583579E25 * 5.29177E-11**3 Hartree/Bohr^3 = 3.398927420868445E-6 Hartree/Bohr^3
            self['stress'] = [ float(i)/comm.KBAR2HARTREEPERBOHR3 for i in output.find('stress').text.split()]
            # calculate the virial, unit in eV
            if volume != None:
                self['virial'] = [ float(i) * volume * comm.HARTREE2EV / comm.BOHR2A ** 3   for i in output.find('stress').text.split()]
    
        self['total_mag'] = comm.ifloat(xfmlt(output,['magnetization','total']))
        self['absolute_mag'] = comm.ifloat(xfmlt(output,['magnetization','absolute']))
        if self['energy'] == None or self['natom'] == None:
            self['energy_per_atom'] = None
        else:
            self['energy_per_atom'] = self['energy'] / self['natom']
    
    @ResultQe.register(total_time="the total running time",
                       force_time="the time to do the calculations of force",
                       stress_time="the time to do the calculation of stress")
    def GetTimeParamFromXml(self):
        if self.XMLROOT == None:
            return

        timing = self.XMLROOT.find('timing_info')

        self['total_time'] = comm.ifloat(xfmlt(timing,['total','wall']))
        
        all_time = timing.findall('partial')
        for itime in all_time:
            if itime.get('label') == 'forces':
                self['force_time'] = comm.ifloat(itime.find('wall').text)
            elif itime.get('label') == 'stress':
                self['stress_time'] = comm.ifloat(itime.find('wall').text)

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
        if self.XMLROOT == None:
            return

        self["relax_steps"] = len(self.XMLROOT.findall('step'))

