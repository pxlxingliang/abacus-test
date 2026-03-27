# Geometry Relaxation

Parameters for structure and cell relaxation calculations.

## Relaxation Control

### relax_method
- **Type**: String
- **Default**: `cg`
- **Options**:
  - `cg` - Conjugate gradient
  - `bfgs` - BFGS algorithm
  - `bfgs_trad` - Traditional BFGS
  - `cg_bfgs` - CG first, then switch to BFGS
  - `sd` - Steepest descent
  - `fire` - FIRE algorithm (use with calculation=md, md_type=fire)

### relax_new
- **Type**: Boolean
- **Default**: `True`
- **Description**: Use new CG implementation (2022+).

### relax_nmax
- **Type**: Integer
- **Default**: `1` (SCF), `50` (relax/cell-relax)
- **Description**: Maximum ionic iterations. `0` = dry run (check inputs only).

### relax_cg_thr
- **Type**: Real
- **Default**: `0.5`
- **Description**: Force threshold for switching from CG to BFGS (cg_bfgs method).

### relax_scale_force
- **Type**: Real
- **Default**: `0.5`
- **Description**: Scale factor for first CG step. Smaller = safer for large systems.

### relax_bfgs_init
- **Type**: Real
- **Default**: `0.5`
- **Description**: Initial movement sum for BFGS.

### relax_bfgs_w1
- **Type**: Real
- **Default**: `0.01`
- **Description**: Wolfe condition parameter 1 for BFGS.

### relax_bfgs_w2
- **Type**: Real
- **Default**: `0.5`
- **Description**: Wolfe condition parameter 2 for BFGS.

### relax_bfgs_rmax
- **Type**: Real
- **Default**: `0.8`
- **Description**: Maximum atomic movement in BFGS.

### relax_bfgs_rmin
- **Type**: Real
- **Default**: `1e-5`
- **Description**: Minimum atomic movement. Below this = convergence failure.

## Force & Stress

### cal_force
- **Type**: Boolean
- **Default**: `False` (scf), `True` (relax/cell-relax/md)
- **Description**: Calculate forces at end of electronic iteration.

### cal_stress
- **Type**: Boolean
- **Default**: `True` (cell-relax), `False` (otherwise)
- **Description**: Calculate stress at end of electronic iteration.

### force_thr
- **Type**: Real
- **Default**: `0.001` Ry/Bohr
- **Description**: Force convergence threshold (Ry/Bohr).

### force_thr_ev
- **Type**: Real
- **Default**: `0.0257112` eV/Å
- **Description**: Force convergence threshold (eV/Å). Recommended: 0.04 eV/Å for LCAO.

### force_thr_ev2
- **Type**: Real
- **Default**: `0.0`
- **Description**: Forces below this value are set to 0.

### stress_thr
- **Type**: Real
- **Default**: `0.5`
- **Description**: Stress convergence threshold (compared with largest tensor component).

### press1, press2, press3
- **Type**: Real
- **Default**: `0`
- **Description**: External pressures along three axes (positive = compressive).

## Cell Relaxation Constraints

### fixed_axes
- **Type**: String
- **Default**: `None`
- **Options**:
  - `None` - All axes can relax
  - `volume` - Fixed volume relaxation
  - `shape` - Fixed shape (only lattice constant changes)
  - `a`, `b`, `c` - Fix specific axis
  - `ab`, `ac`, `bc` - Fix two axes
- **Note**: `shape` and `volume` require relax_new=True

### fixed_ibrav
- **Type**: Boolean
- **Default**: `False`
- **Description**: Preserve lattice type during relaxation.

### fixed_atoms
- **Type**: Boolean
- **Default**: `False`
- **Description**: Preserve atomic direct coordinates during variable-cell relaxation.

### cell_factor
- **Type**: Real
- **Default**: `1.2`
- **Description**: Factor for pseudopotential table construction. Should exceed maximum cell contraction.

---

## Quick Examples

### Standard Structure Relaxation
```
calculation relax
relax_method cg
force_thr_ev 0.01
relax_nmax 50
```

### Cell Relaxation with Fixed Volume
```
calculation cell-relax
fixed_axes volume
stress_thr 0.5
relax_nmax 100
```

### CG-BFGS Mixed Method
```
calculation relax
relax_method cg_bfgs
relax_cg_thr 0.5
force_thr_ev 0.02
```

### FIRE Relaxation (via MD)
```
calculation md
md_type fire
md_nstep 1000
```

---

## Related References

- [Molecular Dynamics](molecular-dynamics.md) - MD parameters including FIRE
- [System Variables](system-variables.md) - calculation parameter
- [Output](output.md) - Output structure files during relaxation
