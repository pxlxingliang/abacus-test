import os
import subprocess
import numpy as np
from ase.calculators.calculator import FileIOCalculator, ReadError, CalculatorSetupError, CalculationFailed


# ABACUS stress output is in kbar.
# ASE convention expects stress in eV/Angstrom^3.
KBAR2EVPERA3 = 3.398927420868445E-6 * 27.211396132 / 0.52917721092**3


class AbacusCalculator(FileIOCalculator):
    implemented_properties = [
        'energy', 'forces', 'stress',
        'dipole', 'charges', 'magmom', 'magmoms',
    ]

    def __init__(
        self,
        *,
        atoms=None,
        command=None,
        directory='.',
        pp_dict=None,
        orb_dict=None,
        paw_dict=None,
        kpts=None,
        input_param=None,
    ):
        input_param = input_param or {}

        super().__init__(
            restart=None,
            label='abacus',
            atoms=atoms,
            command=command,
            directory=directory,
            pp_dict=pp_dict,
            orb_dict=orb_dict,
            paw_dict=paw_dict,
            kpts=kpts,
            input_param=input_param,
        )

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        pp_dict = self.parameters.get('pp_dict') or {}
        orb_dict = self.parameters.get('orb_dict') or {}
        paw_dict = self.parameters.get('paw_dict') or {}
        kpts_param = self.parameters.get('kpts')

        # ---------- build ABACUS INPUT dict ----------
        input_params = dict(self.parameters.get('input_param') or {})
        input_params['calculation'] = 'scf'
        input_params['suffix'] = "ABACUS"
        input_params["stru_file"] = None
        input_params["kpoint_file"] = None
        input_params["pseudo_dir"] = None
        input_params["orbital_dir"] = None

        if properties:
            if 'stress' in properties:
                input_params['cal_stress'] = 1
            if 'forces' in properties:
                input_params['cal_force'] = 1

        # ---------- convert ASE Atoms -> AbacusSTRU ----------
        from abacustest.lib_prepare.stru import AbacusSTRU

        stru = AbacusSTRU.from_ase(atoms)

        if not pp_dict and not paw_dict:
            raise CalculatorSetupError('Either pp_dict or paw_dict is required')

        if pp_dict:
            pp_abspath = {elem: os.path.abspath(path) for elem, path in pp_dict.items()}
            stru.set_pp(pp_abspath)

        if orb_dict:
            orb_abspath = {elem: os.path.abspath(path) for elem, path in orb_dict.items()}
            stru.set_orb(orb_abspath)

        if paw_dict:
            paw_abspath = {elem: os.path.abspath(path) for elem, path in paw_dict.items()}
            stru.set_paw(paw_abspath)

        # ---------- write files ----------
        stru.write(os.path.join(self.directory, 'STRU'))

        from abacustest.lib_prepare.abacus import WriteInput, WriteKpt
        WriteInput(input_params, os.path.join(self.directory, 'INPUT'))

        kpt_path = os.path.join(self.directory, 'KPT')
        if kpts_param is not None:
            kpts_list = list(map(float, kpts_param))
            if len(kpts_list) == 6:
                kpt_list = kpts_list
            elif len(kpts_list) == 3:
                kpt_list = kpts_list + [0, 0, 0]
            else:
                raise CalculatorSetupError(f"kpts must be length 3 or 6, got {kpts_param}")
            WriteKpt(kpoint_list=kpt_list, file_name=kpt_path)

    def read_results(self):
        from abacustest.lib_collectdata.abacus.abacus import Abacus

        result = Abacus(path=self.directory)

        energy = result['energy']
        if energy is None:
            raise ReadError('ABACUS calculation did not produce a valid energy')

        self.results = {'energy': energy}

        try:
            forces = result['force']
            self.results['forces'] = np.array(forces, dtype=float).reshape(-1, 3)
        except (KeyError, Exception):
            pass

        try:
            stress = result['stress']
            if len(stress) >= 9:
                s = np.array(stress[:9], dtype=float)
                self.results['stress'] = np.array([
                    s[0], s[4], s[8], s[5], s[2], s[1]
                ]) * KBAR2EVPERA3 * -1
        except (KeyError, Exception):
            pass

        #try:
        #    dipole = result['dipole']
        #    self.results['dipole'] = np.array(dipole, dtype=float)
        #except (KeyError, Exception):
        #    pass

        #try:
        #    charges = result['atom_elec_mul']
        #    self.results['charges'] = np.array(charges, dtype=float)
        #except (KeyError, Exception):
        #    pass

        try:
            magmom = result['total_mag']
            self.results['magmom'] = magmom
        except (KeyError, Exception):
            pass

        try:
            atom_mag = result['atom_mag']
            self.results['magmoms'] = np.array(atom_mag, dtype=float)
        except (KeyError, Exception):
            pass

    def todict(self, skip_default=True):
        d = super().todict(skip_default=skip_default)
        for key in ('pp_dict', 'orb_dict', 'paw_dict', 'kpts','input_param'):
            if key in self.parameters:
                d[key] = self.parameters[key]
        return d
