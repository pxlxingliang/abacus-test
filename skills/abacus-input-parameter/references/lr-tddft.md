# LR-TDDFT (Linear Response TDDFT)

Parameters for linear response time-dependent DFT (Casida equation) calculations.

## XC Kernel

### xc_kernel
- **Type**: String
- **Default**: `LDA`
- **Options**: `RPA`, `LDA`, `PBE`, `HSE`, `HF`
- **Description**: Exchange-correlation kernel for LR-TDDFT.

### lr_init_xc_kernel
- **Type**: String
- **Default**: `"default"`
- **Options**:
  - `"default"` - Calculate f_xc from ground-state charge density
  - `"file"` - Read f_xc from .cube files (LDA-type only)
  - `"from_charge_file"` - Calculate f_xc from charge density in .cube files

## Solver

### lr_solver
- **Type**: String
- **Default**: `dav`
- **Options**:
  - `dav`/`dav_subspace`/`cg` - Iterative diagonalization (Davidson/CG)
  - `lapack` - Full matrix diagonalization (LAPACK)
  - `spectrum` - Calculate absorption spectrum only (no Casida solution)

### lr_thr
- **Type**: Real
- **Default**: `1e-2`
- **Description**: Convergence threshold for iterative diagonalization.

## Orbital Selection

### nocc
- **Type**: Integer
- **Default**: `nbands`
- **Description**: Number of occupied orbitals (up to HOMO).
- **Note**: Auto-set to nelec/2 if illegal value.

### nvirt
- **Type**: Integer
- **Default**: `1`
- **Description**: Number of virtual orbitals (from LUMO).

### lr_nstates
- **Type**: Integer
- **Default**: `0`
- **Description**: Number of 2-particle states to solve.

### lr_unrestricted
- **Type**: Boolean
- **Default**: `False`
- **Description**: Use unrestricted LR-TDDFT (doubled matrix size).
- **Options**:
  - `True` - Always unrestricted
  - `False` - Unrestricted only for open-shell

## Absorption Spectrum

### abs_wavelen_range
- **Type**: Real Real
- **Default**: `0.0 0.0`
- **Description**: Wavelength range for absorption spectrum (min max).

### abs_broadening
- **Type**: Real
- **Default**: `0.01`
- **Description**: Broadening factor η for absorption spectrum.

## Output

### out_wfc_lr
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output LR-TDDFT eigenstates and eigenvectors.
- **Files**: `Excitation_Energy.dat`, `Excitation_Amplitude_${rank}.dat`

## RI-LVL Benchmark

### ri_hartree_benchmark
- **Type**: String
- **Default**: `none`
- **Options**:
  - `aims` - Use FHI-aims output files (RI-LVL tensors + KS eigenpairs)
  - `abacus` - Use ABACUS RI-LVL tensors with ABACUS KS eigenpairs
  - `none` - Standard Poisson equation + grid integration
- **Note**: Molecular systems, single processor, large supercell required.

### aims_nbasis
- **Type**: Integer array (ntype)
- **Default**: `{}` (use ABACUS basis)
- **Description**: FHI-aims basis set size per atom type.

---

## Quick Examples

### Standard LR-TDDFT (LDA)
```
xc_kernel LDA
lr_solver dav
lr_thr 1e-2
nocc 10
nvirt 20
lr_nstates 50
```

### Absorption Spectrum Only
```
lr_solver spectrum
abs_wavelen_range 200 800
abs_broadening 0.05
out_wfc_lr 1
```

### Hybrid Functional LR-TDDFT
```
xc_kernel HSE
lr_solver dav_subspace
lr_unrestricted 1
```

### RI-LVL Benchmark (FHI-aims)
```
ri_hartree_benchmark aims
aims_nbasis 5 14 14
```

---

## Related References

- [TDDFT](tddft.md) - Real-time TDDFT
- [Electronic Structure](electronic-structure.md) - SCF settings
- [Output](output.md) - Output control
