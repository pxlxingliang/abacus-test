---
name: abacus-input-parameter
description: "ABACUS INPUT parameter reference guide. Use when: user needs to configure ABACUS calculations, understand INPUT parameters, or troubleshoot convergence issues."
metadata: { "openclaw": { "emoji": "📖", "requires": {} } }
---

# ABACUS INPUT Parameter Reference

Comprehensive reference for ABACUS INPUT file parameters, organized by category.

## Quick Access by Category

| Category | Parameters | Reference |
|----------|------------|-----------|
| **System Variables** | `suffix`, `calculation`, `esolver_type`, `symmetry`, `symmetry_prec`, `symmetry_autoclose`, `kpar`, `bndpar`, `latname`, `init_wfc`, `init_chg`, `init_vel`, `mem_saver`, `diago_proc`, `nbspline`, `kspacing`, `min_dist_coef`, `device`, `precision` | [`system-variables.md`](references/system-variables.md) |
| **Input Files** | `stru_file`, `kpoint_file`, `pseudo_dir`, `orbital_dir`, `read_file_dir`, `restart_load`, `restart_save`, `wannier_card` | [`input-files.md`](references/input-files.md) |
| **Plane Wave** | `ecutwfc`, `ecutrho`, `nx/ny/nz`, `ndx/ndy/ndz`, `pw_seed`, `pw_diag_thr`, `pw_diag_nmax`, `pw_diag_ndim`, `fft_mode`, `erf_ecut`, `erf_height`, `erf_sigma`, `diago_smooth_ethr`, `use_k_continuity` | [`plane-wave.md`](references/plane-wave.md) |
| **LCAO** | `lmaxmax`, `nb2d`, `lcao_ecut`, `lcao_dk`, `lcao_dr`, `lcao_rmax`, `search_radius`, `search_pbc`, `bx/by/bz`, `elpa_num_thread`, `num_stream` | [`lcao.md`](references/lcao.md) |
| **Electronic Structure** | `basis_type`, `ks_solver`, `nbands`, `nbands_mul`, `nelec`, `nelec_delta`, `nupdown`, `dft_functional`, `xc_temperature`, `nspin`, `lspinorb`, `noncolin`, `soc_lambda`, `smearing_method`, `smearing_sigma`, `mixing_type`, `mixing_beta`, `mixing_ndim`, `mixing_gg0`, `scf_nmax`, `scf_thr`, `gamma_only`, `chg_extrap` | [`electronic-structure.md`](references/electronic-structure.md) |
| **Geometry Relaxation** | `relax_method`, `relax_new`, `relax_nmax`, `relax_cg_thr`, `relax_scale_force`, `relax_bfgs_*`, `cal_force`, `cal_stress`, `force_thr`, `force_thr_ev`, `stress_thr`, `press1/2/3`, `fixed_axes`, `fixed_ibrav`, `fixed_atoms`, `cell_factor` | [`geometry-relaxation.md`](references/geometry-relaxation.md) |
| **Output** | `out_level`, `out_alllog`, `out_chg`, `out_pot`, `out_wfc_*`, `out_dos`, `out_band`, `out_mul`, `out_stru`, `out_mat_*`, `restart_save`, `restart_load`, `out_interval`, `out_ndigits`, `towannier90` | [`output.md`](references/output.md) |
| **Molecular Dynamics** | `md_type`, `md_nstep`, `md_dt`, `md_tfirst`, `md_tlast`, `md_restart`, `md_thermostat`, `md_tfreq`, `md_tchain`, `md_pmode`, `md_pcouple`, `md_pfirst`, `md_plast`, `msst_*`, `lj_*`, `dp_*`, `dump_*` | [`molecular-dynamics.md`](references/molecular-dynamics.md) |
| **DFT+U** | `dft_plus_u`, `orbital_corr`, `hubbard_u`, `yukawa_potential`, `yukawa_lambda`, `uramping`, `omc`, `onsite_radius` | [`dft-u.md`](references/dft-u.md) |
| **vdW Correction** | `vdw_method`, `vdw_s6`, `vdw_s8`, `vdw_a1`, `vdw_a2`, `vdw_d`, `vdw_abc`, `vdw_C6_file`, `vdw_R0_file`, `vdw_cutoff_*`, `vdw_cn_thr` | [`vdw-correction.md`](references/vdw-correction.md) |
| **Exact Exchange** | `dft_functional`, `exx_hybrid_alpha`, `exx_hse_omega`, `exx_separate_loop`, `exx_hybrid_step`, `exx_*_threshold`, `exx_distribute_type`, `exx_real_number`, `exx_opt_orb_*` | [`exact-exchange.md`](references/exact-exchange.md) |
| **TDDFT** | `td_vext`, `td_propagator`, `td_ttype`, `td_stype`, `td_gauss_*`, `td_trape_*`, `td_trigo_*`, `td_heavi_*`, `out_dipole`, `out_current`, `ocp`, `ocp_set` | [`tddft.md`](references/tddft.md) |
| **Berry Phase/Wannier90** | `berry_phase`, `gdir`, `towannier90`, `wannier_card`, `wannier_method`, `wannier_spin`, `out_wannier_*` | [`berry-phase.md`](references/berry-phase.md) |
| **DeePKS** | `deepks_out_labels`, `deepks_scf`, `deepks_model`, `deepks_equiv`, `deepks_bandgap`, `deepks_v_delta`, `bessel_descriptor_*` | [`deepks.md`](references/deepks.md) |
| **OFDFT** | `of_kinetic`, `of_method`, `of_conv`, `of_tole`, `of_tolp`, `of_*_weight`, `of_wt_*`, `of_lkt_a`, `of_read_kernel`, `of_full_pw` | [`ofdft.md`](references/ofdft.md) |
| **Electric Field** | `efield_flag`, `efield_dir`, `efield_amp`, `efield_pos_*`, `dip_cor_flag`, `gate_flag`, `zgate`, `block_*` | [`electric-field.md`](references/electric-field.md) |
| **SDFT** | `method_sto`, `nbands_sto`, `nche_sto`, `emin_sto`, `emax_sto`, `seed_sto`, `initsto_*`, `npart_sto` | [`sdft.md`](references/sdft.md) |
| **DOS** | `out_dos`, `dos_edelta_ev`, `dos_sigma`, `dos_scale`, `dos_emin_ev`, `dos_emax_ev`, `dos_nche` | [`dos.md`](references/dos.md) |
| **NAOs** | `bessel_nao_ecut`, `bessel_nao_rcut`, `bessel_nao_smooth`, `bessel_nao_sigma`, `bessel_nao_tolerence` | [`naos.md`](references/naos.md) |
| **Debugging** | `t_in_h`, `vl_in_h`, `vnl_in_h`, `vh_in_h`, `vion_in_h`, `test_force`, `test_stress`, `test_skip_ewald` | [`debugging.md`](references/debugging.md) |
| **Conductivities** | `cal_cond`, `cond_che_thr`, `cond_dw`, `cond_wcut`, `cond_dt`, `cond_smear`, `cond_fwhm`, `cond_nonlocal` | [`conductivities.md`](references/conductivities.md) |
| **Solvation** | `imp_sol`, `eb_k`, `tau`, `sigma_k`, `nc_k` | [`solvation.md`](references/solvation.md) |
| **QO Analysis** | `qo_switch`, `qo_basis`, `qo_strategy`, `qo_screening_coeff`, `qo_thr` | [`qo-analysis.md`](references/qo-analysis.md) |
| **PEXSI** | `pexsi_npole`, `pexsi_inertia`, `pexsi_nmax`, `pexsi_ordering`, `pexsi_method`, `pexsi_mu*`, `pexsi_thr`, `pexsi_temp` | [`pexsi.md`](references/pexsi.md) |
| **LR-TDDFT** | `xc_kernel`, `lr_solver`, `lr_thr`, `nocc`, `nvirt`, `lr_nstates`, `abs_wavelen_range`, `abs_broadening`, `ri_hartree_benchmark` | [`lr-tddft.md`](references/lr-tddft.md) |
| **RDMFT** | `rdmft`, `rdmft_power_alpha` | [`rdmft.md`](references/rdmft.md) |

---

## Common Parameters for Quick Setup

### Minimal SCF Calculation

```
INPUT_PARAMETERS
suffix ABACUS
calculation scf
ecutwfc 50
kspacing 0.14
scf_thr 1e-8
scf_nmax 100
```

### Structure Relaxation

```
INPUT_PARAMETERS
calculation relax
ecutwfc 50
kspacing 0.14
relax_method cg
force_thr_ev 0.01
relax_nmax 50
```

### Molecular Dynamics

```
INPUT_PARAMETERS
calculation md
md_type nvt
md_nstep 1000
md_dt 1.0
md_tfirst 300
md_tlast 300
```

### Hybrid Functional (HSE06)

```
INPUT_PARAMETERS
dft_functional hse
exx_hybrid_alpha 0.25
exx_hse_omega 0.11
```

### Spin-Polarized Calculation

```
INPUT_PARAMETERS
nspin 2
mixing_beta 0.4
```

### SOC Calculation

```
INPUT_PARAMETERS
lspinorb true
nspin 4
symmetry -1
```

---

## Parameter Types

| Type | Description | Example |
|------|-------------|---------|
| `String` | Text value | `calculation scf` |
| `Integer` | Whole number | `nspin 2` |
| `Real` | Decimal number | `ecutwfc 50.0` |
| `Boolean` | True/False or 0/1 | `symmetry 1` |

---

## How to Use

1. **Browse by category**: Click on a category link above to see all parameters in that group
2. **Search for specific parameters**: Each reference file contains detailed descriptions
3. **Check defaults**: Default values are provided for each parameter
4. **Note conditions**: Some parameters have conditional defaults or requirements

---

## Related Skills

- **Input preparation**: [`abacustest-prepare`](../abacustest-prepare/SKILL.md)
- **Structure inputs**: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
- **Job submission**: [`abacustest-submit`](../abacustest-submit/SKILL.md)
- **Result extraction**: [`abacustest-extract-dft-results`](../abacustest-extract-dft-results/SKILL.md)

---

## External Resources

- **ABACUS Documentation**: https://abacus.oss.cn-north-1.aliyuncs.com/
- **INPUT Parameters**: https://github.com/deepmodeling/abacus-develop/blob/develop/docs/advanced/input.md
- **STRU File**: https://github.com/deepmodeling/abacus-develop/blob/develop/docs/advanced/input_files/stru.md
- **KPT File**: https://github.com/deepmodeling/abacus-develop/blob/develop/docs/advanced/input_files/kpt.md
