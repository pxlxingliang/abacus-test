# Electronic Structure

Parameters controlling electronic structure calculations, including basis set, diagonalization, spin, smearing, and SCF convergence.

## Basis Set & Solver

### basis_type
- **Type**: String
- **Default**: `pw`
- **Options**:
  - `pw` - Plane-wave basis set
  - `lcao` - Localized atomic orbitals
  - `lcao_in_pw` - LCAO expanded in PW (non-SCF)

### ks_solver
- **Type**: String
- **Default**: Depends on basis
- **Options**:
  - `cg` - Conjugate gradient (PW default)
  - `bpcg` - Block-parallel CG (GPU-friendly)
  - `dav` - Davidson algorithm
  - `dav_subspace` - Davidson without orthogonalization (recommended for efficiency)
  - `lapack` - LAPACK (serial only)
  - `genelpa` - ELPA for LCAO
  - `scalapack_gvx` - ScaLAPACK
  - `cusolver` - CUDA solver (single GPU)
  - `cusolvermp` - CUDA multi-GPU
  - `elpa` - ELPA (CPU/GPU)

### nbands
- **Type**: Integer
- **Default**: Conditional (see below)
- **Description**: Number of Kohn-Sham orbitals. Recommended to set explicitly for smearing calculations.
- **Default Conditions**:
  - `nspin=1`: max(nbands_mul × occupied_bands, occupied_bands + 10)
  - `nspin=2`: max(nbands_mul × nelec_spin, nelec_spin + 10)
  - `nspin=4`: max(nbands_mul × nelec, nelec + 20)

### nbands_mul
- **Type**: Real
- **Default**: `1.2`
- **Description**: Multiplier for default nbands formula. Use larger values (e.g., 2.0) for complex systems near Fermi surface.

### nelec
- **Type**: Real
- **Default**: `0.0`
- **Description**: Total number of electrons. `0.0` = calculated from valence electrons (neutral system).

### nelec_delta
- **Type**: Real
- **Default**: `0.0`
- **Description**: Additional electrons: total = nelec + nelec_delta.

### nupdown
- **Type**: Real
- **Default**: `0.0`
- **Description**: Spin constraint (n_up - n_down). Range: [-nelec, nelec]. Used in constrained DFT.

## DFT Functional

### dft_functional
- **Type**: String
- **Default**: From pseudopotential file
- **Description**: Exchange-correlation functional. Supports LDA, GGA, meta-GGA, and hybrid functionals.
- **Common Values**:
  - LDA: `lda`, `pwl`
  - GGA: `pbe`, `pbesol`, `revpbe`, `blyp`, `bp86`, `pw91`
  - meta-GGA: `scan` (requires LIBXC)
  - Hybrid: `pbe0`, `hse`, `hf`, `b3lyp` (requires LIBXC)
- **LIBXC Format**: Combine components with `+`, e.g., `GGA_X_RPBE+GGA_C_PBE`

### xc_temperature
- **Type**: Real
- **Default**: `0.0`
- **Description**: Temperature for temperature-dependent XC functionals (KSDT, etc.).

## Spin & SOC

### nspin
- **Type**: Integer
- **Default**: `1`
- **Options**:
  - `1` - Spin degeneracy (non-spin-polarized)
  - `2` - Collinear spin-polarized
  - `4` - Non-collinear (auto-set for SOC)

### lspinorb
- **Type**: Boolean
- **Default**: `False`
- **Description**: Enable spin-orbit coupling.
- **Requirements**:
  - nspin auto-set to 4
  - Symmetry disabled
  - Full-relativistic pseudopotentials (has_so=true)

### noncolin
- **Type**: Boolean
- **Default**: `False`
- **Description**: Allow non-collinear magnetic moments.
- **Note**: Cannot use with gamma_only=true

### soc_lambda
- **Type**: Real
- **Default**: `1.0`
- **Range**: [0.0, 1.0]
- **Description**: Modulate SOC strength. 0.0 = scalar-relativistic, 1.0 = full SOC.

## Smearing

### smearing_method
- **Type**: String
- **Default**: `gauss`
- **Options**:
  - `fixed` - Fixed occupations (insulators only)
  - `gauss`/`gaussian` - Gaussian smearing
  - `mp` - Methfessel-Paxton (metals)
  - `mp2` - 2nd-order MP (metals)
  - `mv`/`cold` - Marzari-Vanderbilt
  - `fd` - Fermi-Dirac

### smearing_sigma
- **Type**: Real
- **Default**: `0.015` Ry
- **Description**: Smearing energy width.

### smearing_sigma_temp
- **Type**: Real
- **Default**: `2 * smearing_sigma / kB`
- **Description**: Smearing temperature (alternative to smearing_sigma).

## Charge Mixing

### mixing_type
- **Type**: String
- **Default**: `broyden`
- **Options**:
  - `plain` - Simple mixing
  - `pulay` - Pulay method
  - `broyden` - Modified Broyden (faster convergence)

### mixing_beta
- **Type**: Real
- **Default**: `0.8` (nspin=1), `0.4` (nspin=2,4)
- **Description**: Mixing parameter. Lower = more stable but slower.
- **Recommendations**:
  - `0.8` for nspin=1
  - `0.4` for nspin=2,4
  - `0.1` or less for difficult convergence

### mixing_beta_mag
- **Type**: Real
- **Default**: `4*mixing_beta` (max 1.6)
- **Description**: Mixing parameter for magnetic density.

### mixing_ndim
- **Type**: Integer
- **Default**: `8`
- **Description**: Mixing dimensions (history steps) for Pulay/Broyden.

### mixing_restart
- **Type**: Real
- **Default**: `0`
- **Description**: SCF restart threshold when drho < mixing_restart.

### mixing_gg0
- **Type**: Real
- **Default**: `1.0`
- **Description**: Kerker preconditioning parameter. `0` = disabled. Recommended for metals.

### mixing_gg0_mag
- **Type**: Real
- **Default**: `0.0`
- **Description**: Kerker preconditioning for magnetic density (not recommended unless necessary).

### mixing_angle
- **Type**: Real
- **Default**: `-10.0`
- **Description**: Angle mixing for non-collinear calculations. Use `1.0` for angle mixing.

## SCF Convergence

### scf_nmax
- **Type**: Integer
- **Default**: `100`
- **Description**: Maximum electronic iterations.

### scf_thr
- **Type**: Real
- **Default**: `1e-9` (PW), `1e-7` (LCAO)
- **Description**: Charge density convergence threshold.

### scf_ene_thr
- **Type**: Real
- **Default**: `-1.0` (disabled)
- **Description**: Energy convergence threshold.

### scf_thr_type
- **Type**: Integer
- **Default**: `1` (PW), `2` (LCAO)
- **Options**:
  - `1` - Coulomb norm: Δρ_G = ½∬Δρ(r)Δρ(r')/|r-r'| dr dr'
  - `2` - Real-space norm: Δρ_R = (1/Ne)∫|Δρ(r)| dr

### scf_os_stop
- **Type**: Boolean
- **Default**: `False`
- **Description**: Stop SCF early on oscillation detection.

### scf_os_thr
- **Type**: Real
- **Default**: `-0.01`
- **Description**: Slope threshold for oscillation detection.

### scf_os_ndim
- **Type**: Integer
- **Default**: `mixing_ndim`
- **Description**: Number of iterations for oscillation detection.

## Initialization & Extrapolation

### gamma_only
- **Type**: Integer (0/1)
- **Default**: `0`
- **Description**: Gamma-point only algorithm (faster, no KPT file needed).

### chg_extrap
- **Type**: String
- **Default**: `first-order` (relax), `second-order` (MD), `atomic` (else)
- **Options**: `atomic`, `first-order`, `second-order`
- **Description**: Charge density extrapolation for geometry relaxation/MD.

### pw_diag_thr
- **Type**: Real
- **Default**: `0.01`
- **Description**: Diagonalization threshold for first SCF iteration.

### pw_diag_nmax
- **Type**: Integer
- **Default**: `40`
- **Description**: Maximum diagonalization iterations.

### pw_diag_ndim
- **Type**: Integer
- **Default**: `4`
- **Description**: Davidson workspace dimension.

### diago_smooth_ethr
- **Type**: Boolean
- **Default**: `False`
- **Description**: Use smooth threshold strategy for empty states.

### use_k_continuity
- **Type**: Boolean
- **Default**: `False`
- **Description**: Use k-point continuity for wavefunction initialization (requires diago_smooth_ethr=1).

### initsto_ecut
- **Type**: Real
- **Default**: `0.0`
- **Description**: Initial stochastic wavefunction cutoff (SDFT).

---

## Related References

- [System Variables](system-variables.md) - Core configuration
- [Plane Wave](plane-wave.md) - PW-specific parameters
- [LCAO](lcao.md) - LCAO-specific parameters
- [SDFT](sdft.md) - Stochastic DFT parameters
