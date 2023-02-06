import os,sys,glob,re
from ..resultQe import ResultQe
from .. import comm
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element

xfmlt = comm.XmlFindMultiLayerText
xfml = comm.XmlFindMultiLayer

class Qe(ResultQe):

    @ResultQe.register(version="the version of QE",
                       ncore = "the mpi cores")
    def GetGeneralInfo(self):
        if self.XMLROOT != None:
            self['version'] = self.XMLROOT.find("general_info/creator[@name='PWSCF']").attrib['VERSION']
            self['ncore'] = self.XMLROOT.find("parallel_info/nprocs").text
    
    @ResultQe.register(natom="total atom number",
                       nbands="band number",
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
        self['nbands'] = comm.iint(xfmlt(input_param,['bands','nbnd']))

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
                       ibzk="irreducible K point number",
                       nkstot="the total k points",
                       kpt="list, the K POINTS setting",
                       force="list, the force of all atoms, [atom1x,atom1y,atom1z,atom2x,atom2y,atom2z...]",
                       stress="list, the stress, [xx,xy,xz,yx,yy,yz,zx,zy,zz]",
                       total_mag="total magnization",
                       absolute_mag="total absolute magnization",
                       nelec="total electron number",
                       energy_per_atom="total energy divided by natom, unit in eV",)
    def GetOutputParamFromXml(self):
        if self.XMLROOT == None:
            return

        output = self.XMLROOT.find('output')

        self['converge'] = comm.ibool(xfmlt(output,['convergence_info','scf_conv','convergence_achieved']))
        self['scf_steps'] = comm.iint(xfmlt(output,['convergence_info','scf_conv','n_scf_steps']))
        self['ibzk'] = comm.iint(xfmlt(output,['band_structure','nks']))
        self['nelec'] = comm.ifloat(xfmlt(output,['band_structure','nelec']))
        
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
                
        if output.find('forces'):
            self['force'] = [float(i) for i in output.find('forces').text.split()]
        if output.find('stress'):
            self['stress'] = [float(i) for i in output.find('stress').text.split()]
    
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



