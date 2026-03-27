# Plane Wave Parameters

Parameters specific to plane-wave basis set calculations.

## Energy Cutoffs

### ecutwfc
- **Type**: Real
- **Default**: `50 Ry` (PW), `100 Ry` (LCAO)
- **Unit**: Rydberg
- **Description**: Plane-wave energy cutoff for wavefunctions.
- **Note**: Even for LCAO, ecutwfc is needed for pseudopotential and force calculations.

### ecutrho
- **Type**: Real
- **Default**: `4*ecutwfc`
- **Unit**: Rydberg
- **Description**: Energy cutoff for charge density and potential.
- **Notes**:
  - Norm-conserving PP: stick to default
  - Ultrasoft PP: 8-12 × ecutwfc recommended
  - GGA functionals in vacuum: may need higher values

## FFT Grid

### nx, ny, nz
- **Type**: Integer
- **Default**: `0` (calculated from ecutrho)
- **Description**: FFT grid points in x, y, z directions.
- **Note**: Must specify all three or none.

### ndx, ndy, ndz
- **Type**: Integer
- **Default**: `0` (calculated from ecutwfc)
- **Description**: FFT grid for dense charge density (ultrasoft PP).
- **Note**: Must be used with nx, ny, nz.

### fft_mode
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - FFTW_ESTIMATE
  - `1` - FFTW_MEASURE
  - `2` - FFTW_PATIENT
  - `3` - FFTW_EXHAUSTIVE
- **Description**: FFTW planning mode.

## Diagonalization (PW)

### pw_diag_thr
- **Type**: Real
- **Default**: `0.01`
- **Description**: Diagonalization threshold for first SCF iteration.
- **Note**: For nscf with PW, should be ≤ 1e-3.

### pw_diag_nmax
- **Type**: Integer
- **Default**: `40`
- **Description**: Maximum diagonalization iterations (cg/dav/bpcg).

### pw_diag_ndim
- **Type**: Integer
- **Default**: `4`
- **Description**: Davidson workspace dimension (number of wavefunction packets, min 2).

### pw_seed
- **Type**: Integer
- **Default**: `0`
- **Description**: Random seed for wavefunction initialization (PW only).

## Smooth Cutoff (Variable-Cell MD)

### erf_ecut
- **Type**: Real
- **Default**: `0.0`
- **Description**: Constant energy cutoff for kinetic functional modification.

### erf_height
- **Type**: Real
- **Default**: `0.0`
- **Description**: Height of energy step for G² > erf_ecut.

### erf_sigma
- **Type**: Real
- **Default**: `0.1`
- **Description**: Width of energy step.
- **Formula**: G² → G² + erf_height × (1 + erf((G² - erf_ecut)/erf_sigma))
- **Reference**: M. Bernasconi et al., J. Phys. Chem. Solids 56, 501 (1995)

## Diagonalization Improvements

### diago_smooth_ethr
- **Type**: Boolean
- **Default**: `False`
- **Description**: Use smooth threshold strategy (larger threshold for empty states).
- **Effect**: Improves efficiency without affecting ground-state properties.

### use_k_continuity
- **Type**: Boolean
- **Default**: `False`
- **Description**: Use k-point continuity for wavefunction initialization.
- **Requirements**: Must use with diago_smooth_ethr=1
- **Benefit**: Reduces initialization cost for dense k-point sampling.

---

## Quick Examples

### Standard SCF
```
ecutwfc 50
ecutrho 200
fft_mode 0
```

### High-Precision Calculation
```
ecutwfc 80
ecutrho 320
pw_diag_thr 1e-4
```

### Variable-Cell MD
```
erf_ecut 60
erf_height 0.5
erf_sigma 0.1
ref_cell_factor 1.05
```

### Efficient Large System
```
diago_smooth_ethr 1
use_k_continuity 1
fft_mode 1
```

---

## Related References

- [Electronic Structure](electronic-structure.md) - General SCF parameters
- [LCAO](lcao.md) - LCAO-specific parameters
- [System Variables](system-variables.md) - basis_type parameter
