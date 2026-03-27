# PEXSI (Pole Expansion Selected Inversion)

Parameters for PEXSI solver - linear scaling electronic structure method.

## PEXSI Control

### pexsi_npole
- **Type**: Integer
- **Default**: `40`
- **Description**: Number of poles in pole expansion (must be even).

### pexsi_inertia
- **Type**: Boolean
- **Default**: `True`
- **Description**: Use inertia counting at the beginning.

### pexsi_nmax
- **Type**: Integer
- **Default**: `80`
- **Description**: Maximum PEXSI iterations after each inertia counting.

## Communication & Storage

### pexsi_comm
- **Type**: Boolean
- **Default**: `True`
- **Description**: Construct PSelInv communication pattern.

### pexsi_storage
- **Type**: Boolean
- **Default**: `True`
- **Description**: Use symmetric storage for selected inversion.

## Ordering

### pexsi_ordering
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - Parallel ordering (ParMETIS)
  - `1` - Sequential ordering (METIS)
  - `2` - Multiple minimum degree ordering

### pexsi_row_ordering
- **Type**: Integer
- **Default**: `1`
- **Options**:
  - `0` - No row permutation
  - `1` - Make diagonal entries larger than off-diagonal

### pexsi_nproc
- **Type**: Integer
- **Default**: `1`
- **Description**: Number of processors for ParMETIS (only if pexsi_ordering=0).

## Matrix Properties

### pexsi_symm
- **Type**: Boolean
- **Default**: `True`
- **Description**: Matrix is symmetric.

### pexsi_trans
- **Type**: Boolean
- **Default**: `False`
- **Description**: Factorize transpose of matrix.

## Method & Parallelization

### pexsi_method
- **Type**: Integer
- **Default**: `1`
- **Options**:
  - `1` - Cauchy Contour Integral
  - `2` - Moussa optimized method

### pexsi_nproc_pole
- **Type**: Integer
- **Default**: `1`
- **Description**: Point parallelization. Recommended: 2.

## Physical Parameters

### pexsi_temp
- **Type**: Real
- **Default**: `0.015` Ry
- **Description**: Temperature in Fermi-Dirac distribution (same as smearing_sigma for FD).

### pexsi_gap
- **Type**: Real
- **Default**: `0`
- **Description**: Spectral gap (can be 0 in most cases).

### pexsi_delta_e
- **Type**: Real
- **Default**: `20`
- **Description**: Upper bound for spectral radius of S⁻¹H.

## Chemical Potential

### pexsi_mu_lower
- **Type**: Real
- **Default**: `-10`
- **Description**: Initial lower bound for chemical potential μ.

### pexsi_mu_upper
- **Type**: Real
- **Default**: `10`
- **Description**: Initial upper bound for μ.

### pexsi_mu
- **Type**: Real
- **Default**: `0`
- **Description**: Initial guess for μ.

### pexsi_mu_thr
- **Type**: Real
- **Default**: `0.05`
- **Description**: Stopping criterion for μ in inertia counting.

### pexsi_mu_expand
- **Type**: Real
- **Default**: `0.3`
- **Description**: Interval expansion if μ not in initial range.

### pexsi_mu_guard
- **Type**: Real
- **Default**: `0.2`
- **Description**: Safeguard criterion to re-invoke inertia counting.

## Convergence

### pexsi_elec_thr
- **Type**: Real
- **Default**: `0.001`
- **Description**: Stopping criterion (electron number vs exact).

### pexsi_zero_thr
- **Type**: Real
- **Default**: `1e-10`
- **Description**: Threshold for considering CCS matrix element as zero.

---

## Quick Examples

### Standard PEXSI
```
ks_solver elpa
pexsi_npole 40
pexsi_inertia 1
pexsi_nmax 80
```

### Large System (Memory Efficient)
```
pexsi_ordering 0
pexsi_nproc 4
pexsi_nproc_pole 2
pexsi_method 2
```

### High Precision
```
pexsi_npole 60
pexsi_elec_thr 1e-4
pexsi_zero_thr 1e-12
```

---

## Related References

- [Electronic Structure](electronic-structure.md) - ks_solver parameter
- [LCAO](lcao.md) - LCAO-specific settings
- [System Variables](system-variables.md) - General settings
