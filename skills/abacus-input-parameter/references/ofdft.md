# OFDFT (Orbital-Free Density Functional Theory)

Parameters for orbital-free DFT calculations.

## Kinetic Energy Density Functional

### of_kinetic
- **Type**: String
- **Default**: `wt`
- **Options**:
  - `wt` - Wang-Teter
  - `tf` - Thomas-Fermi
  - `vw` - von Weizsäcker
  - `tf+` - TFλvW (λ set by of_vw_weight)
  - `lkt` - Luo-Karasiev-Trickey

## Optimization Method

### of_method
- **Type**: String
- **Default**: `tn`
- **Options**:
  - `cg1` - Polak-Ribiere (standard CG)
  - `cg2` - Hager-Zhang (generally faster)
  - `tn` - Truncated Newton

## Convergence

### of_conv
- **Type**: String
- **Default**: `energy`
- **Options**:
  - `energy` - Total energy change < of_tole
  - `potential` - Norm of potential < of_tolp
  - `both` - Both criteria must be satisfied

### of_tole
- **Type**: Real
- **Default**: `2e-6`
- **Description**: Energy change tolerance.

### of_tolp
- **Type**: Real
- **Default**: `1e-5`
- **Description**: Potential norm tolerance.

## KEDF Parameters

### of_tf_weight
- **Type**: Real
- **Default**: `1.0`
- **Description**: Weight of Thomas-Fermi KEDF.

### of_vw_weight
- **Type**: Real
- **Default**: `1.0`
- **Description**: Weight of von Weizsäcker KEDF.

### of_wt_alpha
- **Type**: Real
- **Default**: `5/6`
- **Description**: α parameter for Wang-Teter KEDF.

### of_wt_beta
- **Type**: Real
- **Default**: `5/6`
- **Description**: β parameter for Wang-Teter KEDF.

### of_wt_rho0
- **Type**: Real
- **Default**: `0.0`
- **Description**: Average electron density.

### of_hold_rho0
- **Type**: Boolean
- **Default**: `False`
- **Description**: Fix ρ0 even if volume changes.
- **Note**: Auto-set to True if of_wt_rho0 ≠ 0.

### of_lkt_a
- **Type**: Real
- **Default**: `1.3`
- **Description**: Parameter a for LKT KEDF.

## Kernel Settings

### of_read_kernel
- **Type**: Boolean
- **Default**: `False`
- **Description**: Read WT kernel from file.

### of_kernel_file
- **Type**: String
- **Default**: `WTkernel.txt`
- **Description**: WT kernel file name.

## FFT Settings

### of_full_pw
- **Type**: Boolean
- **Default**: `True`
- **Description**: Use all plane waves (ignore ecut).

### of_full_pw_dim
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - Either odd or even
  - `1` - Odd only
  - `2` - Even only
- **Note**: Even dimensions may cause FFT errors. Use `1` if nbspline ≠ -1.

---

## Quick Examples

### Standard Wang-Teter OFDFT
```
esolver_type ofdft
of_kinetic wt
of_method tn
of_conv energy
of_tole 2e-6
```

### TFλvW with Custom Weight
```
esolver_type ofdft
of_kinetic tf+
of_vw_weight 0.2
```

### LKT Functional
```
esolver_type ofdft
of_kinetic lkt
of_lkt_a 1.3
```

### Pre-computed Kernel
```
of_read_kernel 1
of_kernel_file my_kernel.txt
```

---

## Related References

- [System Variables](system-variables.md) - esolver_type
- [Plane Wave](plane-wave.md) - FFT settings
