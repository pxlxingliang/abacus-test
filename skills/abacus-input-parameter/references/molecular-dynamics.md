# Molecular Dynamics

Parameters for ab initio molecular dynamics simulations.

## MD Control

### md_type
- **Type**: String
- **Default**: `nvt`
- **Options**:
  - `nve` - NVE ensemble (velocity Verlet)
  - `nvt` - NVT ensemble (thermostat required)
  - `npt` - NPT ensemble (Nose-Hoover)
  - `langevin` - Langevin thermostat
  - `fire` - FIRE relaxation algorithm
  - `msst` - MSST shock wave method

### md_nstep
- **Type**: Integer
- **Default**: `10`
- **Description**: Total MD steps.

### md_dt
- **Type**: Real
- **Default**: `1.0` fs
- **Description**: Time step.

### md_tfirst, md_tlast
- **Type**: Real
- **Default**: No default
- **Description**: Initial and final temperature (K).
- **Notes**:
  - If unset or < 0, init_vel auto-set to true
  - If different, temperature ramps from tfirst to tlast
  - md_tlast only used in NVT/NPT

### md_restart
- **Type**: Boolean
- **Default**: `False`
- **Description**: Restart MD from `${read_file_dir}/Restart_md.dat`.

### md_restartfreq
- **Type**: Integer
- **Default**: `5`
- **Description**: Frequency of restart file output.

### md_dumpfreq
- **Type**: Integer
- **Default**: `1`
- **Description**: Frequency of MD_dump output.

### md_seed
- **Type**: Integer
- **Default**: `-1`
- **Description**: Random seed for MD. `< 0` = no srand().

## Thermostats (NVT)

### md_thermostat
- **Type**: String
- **Default**: `nhc`
- **Options**:
  - `nhc` - Nose-Hoover chain
  - `anderson` - Anderson thermostat
  - `berendsen` - Berendsen thermostat
  - `rescaling` - Velocity rescaling (method 1)
  - `rescale_v` - Velocity rescaling (method 2)

### md_tfreq
- **Type**: Real
- **Default**: `1/40/md_dt`
- **Description**: Temperature oscillation frequency. Range: 1/(40×md_dt) to 1/(100×md_dt).

### md_tchain
- **Type**: Integer
- **Default**: `1`
- **Description**: Number of thermostats in Nose-Hoover chain.

### md_nraise
- **Type**: Integer
- **Default**: `1`
- **Description**: Thermostat parameter:
  - Anderson: collision frequency = 1/md_nraise
  - Berendsen: rise time τ = md_nraise × md_dt
  - Rescale_v: rescale every md_nraise steps

### md_tolerance
- **Type**: Real
- **Default**: `100.0` K
- **Description**: Temperature tolerance for velocity rescaling.

## Barostat (NPT)

### md_pmode
- **Type**: Integer/String
- **Default**: `0`
- **Options**:
  - `iso` - Isotropic cell fluctuation
  - `aniso` - Anisotropic diagonal fluctuation
  - `tri` - Full triangular cell (6 DOF)
  - `0` - FFT grids fixed, G/K vectors change
  - `2` - FFT grids change per step

### md_pcouple
- **Type**: String
- **Default**: `none`
- **Options**: `none`, `xyz`, `xy`, `xz`, `yz`
- **Description**: Coupled lattice vectors in NPT.

### md_pfirst, md_plast
- **Type**: Real
- **Default**: `-1.0`
- **Description**: Target pressure (same unit as stress). If different, pressure ramps.

### md_pfreq
- **Type**: Real
- **Default**: `1/400/md_dt`
- **Description**: Pressure oscillation frequency.

### md_pchain
- **Type**: Integer
- **Default**: `1`
- **Description**: Number of barostat thermostats.

## MSST Method

### msst_direction
- **Type**: Integer
- **Default**: `2`
- **Options**: `0` (x), `1` (y), `2` (z)
- **Description**: Shock wave direction.

### msst_vel
- **Type**: Real
- **Default**: `0.0`
- **Description**: Shock wave velocity.

### msst_qmass
- **Type**: Real
- **Default**: No default (must set > 0)
- **Description**: Inertia of extended system variable.

### msst_vis
- **Type**: Real
- **Default**: `0.0`
- **Description**: Artificial viscosity.

### msst_tscale
- **Type**: Real
- **Default**: `0.01`
- **Description**: Initial temperature reduction percentage.

### ref_cell_factor
- **Type**: Real
- **Default**: `1.0`
- **Description**: Reference cell size factor (1.02-1.10 typical).

## Langevin Thermostat

### md_damp
- **Type**: Real
- **Default**: `1.0`
- **Description**: Damping parameter for Langevin method.

## Lennard-Jones Potential

### lj_rule
- **Type**: Integer
- **Default**: `2`
- **Options**:
  - `1` - Geometric average
  - `2` - Arithmetic average

### lj_eshift
- **Type**: Boolean
- **Default**: `False`
- **Description**: Shift LJ potential to zero at cutoff.

### lj_rcut
- **Type**: Real
- **Default**: No default
- **Description**: LJ cutoff radius.

### lj_epsilon
- **Type**: Real (array)
- **Default**: No default
- **Description**: LJ ε matrix (N(N+1)/2 components).

### lj_sigma
- **Type**: Real (array)
- **Default**: No default
- **Description**: LJ σ matrix (N(N+1)/2 components).

## Deep Potential MD

### pot_file
- **Type**: String
- **Default**: `graph.pb`
- **Description**: DP potential file.

### dp_rescaling
- **Type**: Real
- **Default**: `1.0`
- **Description**: Temperature-dependent DP rescaling factor.

### dp_fparam
- **Type**: Real (array)
- **Default**: `{}`
- **Description**: Frame parameters for DP.

### dp_aparam
- **Type**: Real (array)
- **Default**: `{}`
- **Description**: Atomic parameters for DP.

## Output Control

### dump_force
- **Type**: Boolean
- **Default**: `True`
- **Description**: Output forces to MD_dump.

### dump_vel
- **Type**: Boolean
- **Default**: `True`
- **Description**: Output velocities to MD_dump.

### dump_virial
- **Type**: Boolean
- **Default**: `True`
- **Description**: Output virials to MD_dump.

### cal_syns
- **Type**: Boolean
- **Default**: `False`
- **Description**: Calculate asynchronous overlap matrix (Hefei-NAMD).

### dmax
- **Type**: Real
- **Default**: `0.01`
- **Description**: Maximum atomic displacement per step (for cal_syns=True).

### init_vel
- **Type**: Boolean
- **Default**: `False`
- **Description**: Read initial velocity from STRU file.

---

## Quick Examples

### NVT Simulation (Nose-Hoover)
```
calculation md
md_type nvt
md_nstep 1000
md_dt 1.0
md_tfirst 300
md_tlast 300
md_thermostat nhc
md_tfreq 0.0005
```

### NPT Simulation
```
calculation md
md_type npt
md_nstep 5000
md_tfirst 300
md_pfirst 0.0
md_pmode iso
```

### Langevin Dynamics
```
calculation md
md_type langevin
md_damp 0.1
md_tfirst 300
```

### FIRE Relaxation
```
calculation md
md_type fire
md_nstep 500
```

### MSST Shock Wave
```
calculation md
md_type msst
msst_direction 2
msst_vel 7.0
msst_qmass 100.0
```

---

## Related References

- [Geometry Relaxation](geometry-relaxation.md) - Relaxation parameters
- [Output](output.md) - MD output control
- [System Variables](system-variables.md) - calculation parameter
