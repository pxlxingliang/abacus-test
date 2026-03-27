# LCAO Parameters

Parameters specific to Localized Atomic Orbital (LCAO) basis set calculations.

## LCAO Basis Control

### lmaxmax
- **Type**: Integer
- **Default**: `2`
- **Description**: Maximum angular momentum channels. If ≠ 2, overrides the value from LCAO datasets.

### nb2d
- **Type**: Integer
- **Default**: `0` (auto)
- **Description**: 2D processor array arrangement for wavefunction matrix distribution.
- **Auto Settings**:
  - size ≤ 500: nb2d = 1
  - 500 < size ≤ 1000: nb2d = 32
  - size > 1000: nb2d = 64

## Two-Center Integrals

### lcao_ecut
- **Type**: Real
- **Default**: `ecutwfc`
- **Unit**: Ry
- **Description**: Energy cutoff for two-center integrals.

### lcao_dk
- **Type**: Real
- **Default**: `0.01` Bohr⁻¹
- **Description**: k-space grid spacing for two-center integrals.

### lcao_dr
- **Type**: Real
- **Default**: `0.01` Bohr
- **Description**: r-space grid spacing for two-center integrals.

### lcao_rmax
- **Type**: Real
- **Default**: `30` Bohr
- **Description**: Maximum distance for two-center integration table.

## Neighbor Search

### search_radius
- **Type**: Real
- **Default**: `-1` (auto)
- **Description**: Radius for finding neighboring atoms. Auto-determined from orbital and beta projector cutoffs.

### search_pbc
- **Type**: Boolean
- **Default**: `True`
- **Description**: Include periodic images in neighbor search.

## Grid Integral

### bx, by, bz
- **Type**: Integer
- **Default**: `0` (auto)
- **Description**: Grid grouping for matrix operations (x, y, z directions).

## ELPA Solver

### elpa_num_thread
- **Type**: Integer
- **Default**: `-1` (all threads)
- **Description**: Number of threads per ELPA calculation.

## GPU Acceleration

### num_stream
- **Type**: Integer
- **Default**: `4`
- **Description**: Number of CUDA streams for LCAO computation.
- **Note**: Most devices work well with ≥ 2 streams.

---

## Quick Examples

### Standard LCAO Calculation
```
basis_type lcao
ecutwfc 100
lcao_ecut 100
```

### Large System (Efficient Parallelization)
```
basis_type lcao
nb2d 64
ks_solver elpa
```

### GPU-Accelerated LCAO
```
basis_type lcao
device gpu
ks_solver cusolver
num_stream 4
```

### High-Precision Two-Center Integrals
```
lcao_ecut 150
lcao_dk 0.005
lcao_dr 0.005
lcao_rmax 40
```

---

## Related References

- [Plane Wave](plane-wave.md) - PW-specific parameters
- [Electronic Structure](electronic-structure.md) - ks_solver, basis_type
- [System Variables](system-variables.md) - device parameter
