from abacustest.lib_tools.tool import Tool
import json, os
from abacustest.lib_prepare.abacus import WriteKpt, WriteInput, ReadInput, AbacusStru

from typing import List, Dict, Any
from typing import Union
from pathlib import Path
from abacustest.constant import EV2RY
from abacustest.lib_prepare.comm import collect_pp
import warnings


class Vasp2AbacusTool(Tool):
    @staticmethod
    def add_args(parser):
        parser.description = "Transform VASP input files to ABACUS input files. ENCUT and EDIFF will not be transformed."
        parser.add_argument('-j', '--job', default=[], action="extend", nargs="*", help='the paths of VASP jobs, should contain INCAR, POSCAR, or KPOINTS')
        parser.add_argument("--pp", default=None, type=str, help="the path of pseudopotential library, or read from enviroment variable ABACUS_PP_PATH")
        parser.add_argument("--orb", default=None, type=str, help="the path of orbital library, or read from enviroment variable ABACUS_ORB_PATH")
        parser.add_argument("--input", default=None, type=str, help="A template input file, used to set the default parameters in ABACUS INPUT. The setting in this file will overwrite the parameters transferred from VASP INCAR.")
        return parser

    def run(self, params):
        jobs = params.job
        pp_path = params.pp
        orb_path = params.orb

        if pp_path is None:
            pp_path = os.environ.get("ABACUS_PP_PATH", None)
        if orb_path is None:
            orb_path = os.environ.get("ABACUS_ORB_PATH", None)

        vasp2abacus = Vasp2Abacus(jobs, pp_path, orb_path, params.input)
        vasp2abacus.run()


class Vasp2Abacus():
    """Transform VASP input files to ABACUS input files.
    Parameters
    ----------
    jobs : List[Union[str,Path]]
        The paths of VASP jobs, will transfer INCAR to ABACUS INPUT, POSCAR to ABACUS STRU, and KPOINTS to ABACUS KPT.
        If file does not exist, then related transfer will be skipped.
    pp_path : Union[str,Path], optional
        The path of pseudopotential library, used to set the pseudopotentials in ABACUS STRU.
        If not set, no pseudopotentials will be set.
    orb_path : Union[str,Path], optional
        The path of orbital library, used to set the orbitals in ABACUS STRU.
        If not set, no orbitals will be set.
    input_file: Union[str,Path], optional
        The path of the template input file, used to set the default parameters in ABACUS INPUT.
        The setting in this file will overwrite the parameters transferred from VASP INCAR.
        If not set, the below default parameters will be used:
        - basis_type: "pw" if lcao is False, "lcao" if lcao is True
        - ecutwfc: 80 Ry for pw, 100 Ry for lcao.
    """
    def __init__(self,
                 jobs: List[Union[str,Path]],
                 pp_path: Union[str,Path]=None,
                 orb_path: Union[str,Path]=None,
                 input_file: Union[str,Path]=None,
                 ):
        self.jobs = jobs
        self.pp_path = pp_path
        self.orb_path = orb_path

        self.pp_lib = None
        self.orb_lib = None
        self.recommand_ecutwfc = None
        if pp_path is not None:
            self.pp_lib = collect_pp(pp_path)
            if os.path.isfile(os.path.join(pp_path, "ecutwfc.json")):
                self.recommand_ecutwfc = json.load(open(os.path.join(pp_path, "ecutwfc.json"), "r"))
        if orb_path is not None:
            self.orb_lib = collect_pp(orb_path)

        if input_file is not None:
            self.input_template = ReadInput(input_file)
        else:
            self.input_template = {}

    def run(self):
        transfer_jobs = []
        for job in self.jobs:
            if not os.path.exists(job):
                raise ValueError(f"Job path {job} does not exist.")
            print(f"\nProcessing job {job}...")

            if os.path.isfile(os.path.join(job, "INCAR")):
                input_param, magmom = self.transfer_incar(os.path.join(job, "INCAR"))
            else:
                input_param = None
                magmom = None

            if os.path.isfile(os.path.join(job, "POSCAR")):
                stru = self.transfer_poscar(os.path.join(job, "POSCAR"))
                elements = stru.get_element(total=False, number=False)
                natoms = stru.get_natoms()

                if magmom is not None:
                    if len(magmom) not in [natoms, 3*natoms]:
                        raise ValueError(f"Magmom length {len(magmom)} does not match the number of atoms {natoms}.")
                    if len(magmom) == 3 * natoms:
                        if input_param.get("nspin", 1) != 4:
                            warnings.warn(f"Warning: the magmom length {len(magmom)} is 3 times the number of atoms, but nspin is {input_param.get('nspin', 1)}, set nspin to 4.")
                            input_param["nspin"] = 4
                        magmom = [magmom[i:i+3] for i in range(0, len(magmom), 3)]
                    else:
                        if input_param.get("nspin", 1) != 2:
                            warnings.warn(f"Warning: the magmom length {len(magmom)} is equal to the number of atoms, but nspin is {input_param.get('nspin', 1)}, set nspin to 2.")
                            input_param["nspin"] = 2
                    stru.set_atommag(magmom)

                pps, orbs = self.set_pporb(elements)
                if pps is not None:
                    pp_names = self.link_pporb(job, pps)
                    stru.set_pp(pp_names)
                if orbs is not None:
                    orb_names = self.link_pporb(job, orbs)
                    stru.set_orb(orb_names)

                stru.write(os.path.join(job, "STRU"))
            else:
                stru, elements = None, None

            if input_param is not None:
                input_param.update(self.input_template)
                if input_param.get("basis_type", "pw") == "pw":
                    recommand_ecutwfc = self.find_recommand_ecutwfc(elements)
                else:
                    recommand_ecutwfc = None
                input_param = self.set_default_input(input_param, recommand_ecutwfc)

                WriteInput(input_param, os.path.join(job, "INPUT"))

            if os.path.isfile(os.path.join(job, "KPOINTS")):
                kpt, model = self.transfer_kpoints(os.path.join(job, "KPOINTS"))
                WriteKpt(kpt, os.path.join(job, "KPT"), model)

            transfer_jobs.append(job)

        return transfer_jobs

    def link_pporb(self, job_path, pporbs):
        for pporb in pporbs:
            tpporb = os.path.join(job_path, os.path.basename(pporb))
            if os.path.isfile(tpporb):
                os.remove(tpporb)
            if os.path.isfile(pporb):
                os.symlink(pporb, tpporb)
        return [os.path.basename(pporb) for pporb in pporbs]

    def set_pporb(self, elements):
        pps = None
        orbs = None

        if self.pp_lib is not None:
            element_pp_lack = [e for e in elements if e not in self.pp_lib]
            if len(element_pp_lack) > 0:
                warnings.warn(f"Warning: the following elements do not have pseudopotentials in the library {self.pp_path}: {', '.join(element_pp_lack)}")

            pps = [self.pp_lib.get(e, "") for e in elements]
        if self.orb_lib is not None:
            element_orb_lack = [e for e in elements if e not in self.orb_lib]
            if len(element_orb_lack) > 0:
                warnings.warn(f"Warning: the following elements do not have orbitals in the library {self.pp_path}: {', '.join(element_orb_lack)}")
            orbs = [self.orb_lib.get(e, "") for e in elements]

        return pps, orbs

    def set_default_input(self, input_param, recommand_ecutwfc=None):
        if "basis_type" not in input_param:
            input_param["basis_type"] = "pw"

        if "ecutwfc" not in input_param:
            if input_param["basis_type"] == "lcao":
                input_param["ecutwfc"] = 100.0
            elif recommand_ecutwfc is not None:
                input_param["ecutwfc"] = recommand_ecutwfc
            else:
                input_param["ecutwfc"] = 80.0

        if "ks_solver" not in input_param:
            if input_param["basis_type"] == "lcao":
                input_param["ks_solver"] = "genelpa"
            else:
                input_param["ks_solver"] = "dav_subspace"

        if "scf_thr" not in input_param:
            if input_param["basis_type"] == "lcao":
                input_param["scf_thr"] = 1e-7
            else:
                input_param["scf_thr"] = 1e-8
        return input_param

    def find_recommand_ecutwfc(self, elements):
        if elements is None:
            return None
        if self.recommand_ecutwfc is None:
            recommand_e_nonone = []
        else:
            recommand_e = [self.recommand_ecutwfc.get(e, None) for e in elements]
            recommand_e_nonone = [e for e in recommand_e if e is not None]

        if len(recommand_e_nonone) > 0:
            print(f"Recommanded ecutwfc for elements {elements} is {recommand_e}, set to {max(recommand_e_nonone)} Ry.")
            return max(recommand_e_nonone)
        else:
            warnings.warn(f"Warning: no recommand ecutwfc found for elements {elements}")
            return None

    def transfer_poscar(self, poscar: Union[str,Path]):
        if not os.path.exists(poscar):
            raise ValueError(f"POSCAR file {poscar} does not exist.")
        stru = AbacusStru.FromDpdata(poscar, "poscar")

        return stru

    def transfer_incar(self, incar: Union[str,Path]):
        from pymatgen.io.vasp import Incar
        vasp_param = Incar.from_file(incar)
        vasp_param = {k.upper(): vasp_param[k] for k in vasp_param}

        input_param = {}

        useless_keys = ["SYSTEM"]
        for key in useless_keys:
            vasp_param.pop(key, None)

        oneone_keys = {
            "SYMPREC": "symmetry_prec",
            "KPAR": "kpar",
            "NELM": "scf_nmax",
            "ISPIN": "nspin"
        }
        for key in oneone_keys:
            if key in vasp_param:
                input_param[oneone_keys[key]] = vasp_param.pop(key)

        if "LSORBIT" in vasp_param or "LNONCOLLINEAR" in vasp_param:
            input_param["lspinorb"] = vasp_param.pop("LSORBIT", False)
            if vasp_param.pop("LNONCOLLINEAR", False) or input_param["lspinorb"]:
                input_param["nspin"] = 4
                input_param["noncolin"] = True

        nsw = vasp_param.pop("NSW", 0)
        has_isif = "ISIF" in vasp_param
        if nsw in [-1, 0]:
            ibrion = vasp_param.pop("IBRION", -1)
        else:
            ibrion = vasp_param.pop("IBRION", 0)
        if ibrion == 0:
            isif = vasp_param.pop("ISIF", 0)
        else:
            isif = vasp_param.pop("ISIF", 2)

        if has_isif:
            input_param["cal_force"] = 1
            input_param["cal_stress"] = 1
        input_param.update(self.transfer_incar_calculation(ibrion, isif))

        if "EDIFFG" in vasp_param:
            if vasp_param["EDIFFG"] < 0:
                input_param["force_thr_ev"] = -1 * vasp_param.pop("EDIFFG")
            else:
                warnings.warn(f"Warning: the EDIFFG value {vasp_param['EDIFFG']} is not supported in ABACUS, use 0.01 eV/Ang.")
                input_param["force_thr_ev"] = 0.01
                vasp_param.pop("EDIFFG", None)

        input_param.update(self.transfer_incar_smear(vasp_param.pop("ISMEAR", None),
                                                     vasp_param.pop("SIGMA", None)))

        magmom = vasp_param.pop("MAGMOM", None)

        input_param.update(self.transfer_incar_u(vasp_param.pop("LDAU", None),
                                                 vasp_param.pop("LDAUTYPE", None),
                                                 vasp_param.pop("LDAUL", None),
                                                 vasp_param.pop("LDAUU", None),
                                                 vasp_param.pop("LDAUJ", None)))

        if len(vasp_param) > 0:
            print(f"Warning: the following keys in INCAR are not supported in ABACUS, ignored:\n\t {', '.join(vasp_param.keys())}")

        return input_param, magmom

    def transfer_incar_calculation(self, ibrion, isif):
        input_param = {}
        if ibrion == -1:
            input_param["calculation"] = "scf"
        elif ibrion == 0:
            input_param["calculation"] = "md"
        elif ibrion in [1, 2]:
            if ibrion == 1:
                warnings.warn("Warning: the IBRION value 1 is not supported in ABACUS, use cg algorithm.")
            input_param["relax_method"] = "cg"

            if isif < 3:
                input_param["calculation"] = "relax"
            elif isif < 9:
                input_param["calculation"] = "cell-relax"
                if isif == 3:
                    pass
                elif isif == 4:
                    input_param["fixed_axes"] = "volume"
                elif isif == 8:
                    input_param["fixed_axes"] = "shape"
                else:
                    warnings.warn(f"Warning: the ISIF value {isif} is not supported in ABACUS, use ISIF=3.")
            else:
                warnings.warn(f"Warning: the ISIF value {isif} is not supported in ABACUS, use ISIF=3.")
                input_param["calculation"] = "cell-relax"
        else:
            warnings.warn(f"Warning: the IBRION value {ibrion} is not supported in ABACUS, do SCF calculation.")
            input_param["calculation"] = "scf"

        return input_param

    def transfer_incar_smear(self, ismear, sigma):
        if ismear is None:
            return {}
        input_param = {}

        if ismear == 0:
            input_param["smearing_method"] = "gauss"
        elif ismear == -14:
            input_param["smearing_method"] = "fd"
        elif ismear == 1:
            input_param["smearing_method"] = "mp"
        elif ismear == 2:
            input_param["smearing_method"] = "mp2"
        else:
            warnings.warn(f"Warning: the ISMEAR value {ismear} is not supported in ABACUS, use gauss.")
            input_param["smearing_method"] = "gauss"

        if sigma is not None:
            input_param["smearing_sigma"] = sigma * EV2RY

        return input_param

    def transfer_incar_u(self, ldau, utype, ul, uu, uj):
        if ldau is None or not ldau:
            return {}

        if utype == 3:
            warnings.warn("Warning: the U type 3 is not supported in ABACUS, skip the settings of U")
            return {}

        if ul is None or uu is None:
            warnings.warn("Warning: the U parameters are not fully set, skip the settings of U")
            return {}

        if uj is None:
            uj = [0 for _ in uu]

        uu = [uu[i]-uj[i] for i in range(len(uu))]

        input_param = {
            "dft_plus_u": 1,
            "orbital_corr": " ".join([str(i) for i in ul]),
            "hubbard_u": " ".join([str(i) for i in uu]),
            "onsite_radius": 3.0
        }
        return input_param

    def transfer_kpoints(self, kpoint_file: Union[str,Path]):
        if not os.path.exists(kpoint_file):
            raise ValueError(f"KPOINTS file {kpoint_file} does not exist.")

        with open(kpoint_file, "r") as f:
            lines = f.readlines()

        assert len(lines) >= 5, "KPOINTS file should have at least 5 lines."

        nks = int(lines[1].strip())
        model = lines[2].strip().lower()
        if nks == 0:
            if model.startswith("g") or model.startswith("m"):
                kpt = [int(i) for i in lines[3].strip().split()] + [float(i) for i in lines[4].strip().split()[:3]]
                model = "gamma" if model.startswith("g") else "mp"
            else:
                raise ValueError(f"The automatic mesh generation method {model} is not supported, should be 'Gamma' or 'MP'.")
        elif nks > 0:
            if model.startswith("l"):
                if lines[3].strip().lower().startswith("c") or lines[3].strip().lower().startswith("k"):
                    model = "line_cartesian"
                else:
                    model = "line"

                k_pairs = []
                pair_lines = [line.strip() for line in lines[4:] if line.strip()]
                assert len(pair_lines) % 2 == 0, "Line k-point pairs should be in pairs."
                for i in range(0, len(pair_lines), 2):
                    k1_lines = pair_lines[i].split()
                    k2_lines = pair_lines[i+1].split()
                    k1 = [float(k) for k in k1_lines[:3]]
                    k2 = [float(k) for k in k2_lines[:3]]
                    if len(k1_lines) > 3:
                        k1.append("# " + " ".join(k1_lines[3:]))
                    if len(k2_lines) > 3:
                        k2.append("# " + " ".join(k2_lines[3:]))
                    k_pairs.append([k1, k2])

                kpt = []
                pre_k2_coord = None
                pre_k2 = None
                for k1, k2 in k_pairs:
                    k1_coord = "_".join([f"{k:.6f}" for k in k1[:3]])
                    k2_coord = "_".join([f"{k:.6f}" for k in k2[:3]])

                    if pre_k2_coord is not None and pre_k2_coord != k1_coord:
                        pre_k2.insert(3, 1)
                        kpt.append(pre_k2)
                    k1.insert(3, nks)
                    kpt.append(k1)

                    pre_k2 = k2
                    pre_k2_coord = k2_coord
                pre_k2.insert(3, 1)
                kpt.append(pre_k2)
            else:
                if model.startswith("c") or model.startswith("k"):
                    model = "cartesian"
                else:
                    model = "direct"
                kpt = []
                for line in lines[3:]:
                    if line.strip() == "":
                        continue
                    kpt.append([float(k) for k in line.strip().split()[:4]])
        else:
            raise ValueError(f"Invalid k-point setting in KPOINTS file: {kpoint_file}. The number of k-points should be non-negative.")

        return kpt, model
