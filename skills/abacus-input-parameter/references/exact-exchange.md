# Exact Exchange (Hybrid Functionals)

Parameters for hybrid functional calculations with exact exchange (EXX).

## Hybrid Functional Control

### dft_functional
- **Type**: String
- **Default**: From pseudopotential
- **Hybrid Options**:
  - `hf` - Pure Hartree-Fock
  - `pbe0` - PBE0 hybrid
  - `hse` - HSE06 (requires LIBXC)
  - `b3lyp` - B3LYP (requires LIBXC)
  - `scan0` - SCAN0 (requires LIBXC)

### exx_hybrid_alpha
- **Type**: Real
- **Default**: `1.0` (hf), `0.25` (others)
- **Description**: Fock exchange fraction: E_X = αE_X^Fock + (1-α)E_X^DFT

### exx_hse_omega
- **Type**: Real
- **Default**: `0.11`
- **Description**: Range-separation parameter for HSE: 1/r = erfc(ωr)/r + erf(ωr)/r

## EXX Iteration Control

### exx_separate_loop
- **Type**: Boolean
- **Default**: `True`
- **Options**:
  - `False` - GGA loop first, then hybrid loop with H_exx update
  - `True` - Two-step method (inner: density matrix, outer: H_exx)

### exx_hybrid_step
- **Type**: Integer
- **Default**: `100`
- **Description**: Maximum outer-loop iterations for Fock exchange.

### exx_mixing_beta
- **Type**: Real
- **Default**: `1.0`
- **Description**: Mixing parameter for density matrix in outer loop.

### exx_lambda
- **Type**: Real
- **Default**: `0.3`
- **Description**: Compensates G=0 divergence in lcao_in_pw method.

## Acceleration Thresholds

### exx_pca_threshold
- **Type**: Real
- **Default**: `1e-4`
- **Description**: PCA threshold for auxiliary basis reduction.
- **Trade-off**: Larger = faster but less accurate.

### exx_c_threshold
- **Type**: Real
- **Default**: `1e-4`
- **Description**: Threshold for neglecting small C^k_ij matrix components.

### exx_v_threshold
- **Type**: Real
- **Default**: `1e-1`
- **Description**: Threshold for truncating V_ab matrix (double-center integrals).

### exx_dm_threshold
- **Type**: Real
- **Default**: `1e-4`
- **Description**: Threshold for truncating density matrix.

### exx_c_grad_threshold
- **Type**: Real
- **Default**: `1e-4`
- **Description**: Threshold for ∇C^k_ij (force/stress).

### exx_v_grad_threshold
- **Type**: Real
- **Default**: `1e-1`
- **Description**: Threshold for ∇V_ab (force/stress).

### exx_schwarz_threshold
- **Type**: Real
- **Default**: `0` (unused)
- **Description**: Cauchy-Schwarz threshold for four-center integrals.

### exx_cauchy_threshold
- **Type**: Real
- **Default**: `1e-7`
- **Description**: Cauchy-Schwarz threshold for Fock exchange matrix.

### exx_cauchy_force_threshold
- **Type**: Real
- **Default**: `1e-7`
- **Description**: Cauchy-Schwarz threshold for force matrix.

### exx_cauchy_stress_threshold
- **Type**: Real
- **Default**: `1e-7`
- **Description**: Cauchy-Schwarz threshold for stress matrix.

### exx_ccp_threshold
- **Type**: Real
- **Default**: `1e-8` (unused)
- **Description**: On-site Coulomb potential cutoff.

### exx_ccp_rmesh_times
- **Type**: Real
- **Default**: `5` (hf/pbe0/scan0), `1.5` (hse), `1` (others)
- **Description**: Radial mesh multiplier for Coulomb potential.

## Parallelization

### exx_distribute_type
- **Type**: String
- **Default**: `htime`
- **Options**:
  - `order` - Distribute by atom pair order
  - `htime` - Balance time per processor (recommended if memory sufficient)
  - `kmeans1`, `kmeans2` - K-means clustering for memory reduction (large systems)

### exx_real_number
- **Type**: Boolean
- **Default**: `True` (gamma_only), `False` (otherwise)
- **Description**: Use double (True) vs complex (False) data type in LibRI.
- **Note**: True improves hybrid functional SCF speed.

## Optimized Atomic Basis Functions

### exx_opt_orb_lmax
- **Type**: Integer
- **Default**: `0`
- **Description**: Maximum l for spherical Bessel functions in opt-ABF generation.

### exx_opt_orb_ecut
- **Type**: Real
- **Default**: `0`
- **Description**: PW cutoff for optimizing radial ABFs (reasonable: 60).

### exx_opt_orb_tolerence
- **Type**: Real
- **Default**: `0`
- **Description**: Threshold for spherical Bessel zeros (reasonable: 1e-12).

## Symmetry & Output

### exx_symmetry_realspace
- **Type**: Boolean
- **Default**: `True`
- **Options**:
  - `False` - Rotate only D(k) from irreducible k-points
  - `True` - Rotate both D(k) and Hexx(R)

### out_ri_cv
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output coefficient tensor C(R) and Coulomb matrix V(R).

### rpa_ccp_rmesh_times
- **Type**: Real
- **Default**: `10`
- **Description**: Radial mesh multiplier for RPA/GW bare Coulomb matrix.

---

## Quick Examples

### HSE06 Calculation
```
dft_functional hse
exx_hybrid_alpha 0.25
exx_hse_omega 0.11
exx_separate_loop 1
```

### PBE0 with Acceleration
```
dft_functional pbe0
exx_hybrid_alpha 0.25
exx_pca_threshold 1e-4
exx_c_threshold 1e-4
exx_v_threshold 1e-1
```

### Hybrid Functional (Fast)
```
dft_functional hse
exx_real_number 1
exx_distribute_type htime
exx_ccp_rmesh_times 1.5
```

---

## Related References

- [Electronic Structure](electronic-structure.md) - dft_functional, scf parameters
- [Output](output.md) - Output parameters
- [System Variables](system-variables.md) - symmetry, gamma_only
