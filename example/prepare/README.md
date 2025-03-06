Enter each folder, and execute `abacustest prepare -p prepare.json`, and then you can get the new inputs in folder abacustest.

If there has several new inputs are generated for each example, then a new subfolder will be created named by 00000, 00001, 00002, ...
The screen output will list for each new input the corresponding STRU file and KPT file/setting, and INPUT setting, like:
```
...
INPUT: 3/INPUT
Invariant INPUT setting: {'calculation': 'scf', 'ntype': 1, 'nbands': 129, 'symmetry': 1, 'scf_thr': 1e-08, 'scf_nmax': 100, 'cal_force': 1, 'cal_stress': 1, 'basis_type': 'pw', 'smearing_method': 'gauss', 'smearing_sigma': 0.002, 'mixing_type': 'broyden', 'mixing_beta': 0.2, 'mixing_gg0': 1.5, 'ks_solver': 'cg', 'pw_seed': 1}
KPT: [5, 6, 7, 8, 9]
STRU: 3/STRU
{'calculation': 'scf', 'ntype': 1, 'nbands': 129, 'symmetry': 1, 'ecutwfc': 50, 'scf_thr': 1e-08, 'scf_nmax': 100, 'cal_force': 1, 'cal_stress': 1, 'basis_type': 'pw', 'smearing_method': 'gauss', 'smearing_sigma': 0.002, 'mixing_type': 'broyden', 'mixing_beta': 0.2, 'mixing_gg0': 1.5, 'ks_solver': 'cg', 'pw_seed': 1}
label 'Pt': link the pseudopotential file './Pt_ONCV_PBE-1.0.upf' defined in 3/STRU
label 'Pt': the orbital file './Pt_gga_9au_100Ry_4s2p2d1f.orb' defined in 3/STRU is not found.
abacustest/0/00000:['0/STRU', [5, 5, 5, 0, 0, 0], {'ecutwfc': 50}]
abacustest/0/00001:['0/STRU', [5, 5, 5, 0, 0, 0], {'ecutwfc': 60}]
abacustest/0/00002:['0/STRU', [5, 5, 5, 0, 0, 0], {'ecutwfc': 70}]
abacustest/0/00003:['0/STRU', [5, 5, 5, 0, 0, 0], {'ecutwfc': 80}]
abacustest/0/00004:['0/STRU', [5, 5, 5, 0, 0, 0], {'ecutwfc': 90}]
abacustest/0/00005:['0/STRU', [6, 6, 6, 0, 0, 0], {'ecutwfc': 50}]
abacustest/0/00006:['0/STRU', [6, 6, 6, 0, 0, 0], {'ecutwfc': 60}]
...
```