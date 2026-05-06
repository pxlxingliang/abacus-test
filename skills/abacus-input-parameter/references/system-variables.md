# System Variables

Core ABACUS system configuration parameters.

## Parameters

### suffix
- **Type**: String
- **Default**: `ABACUS`
- **Description**: Subdirectory name for output. ABACUS generates `OUT.${suffix}` in the working directory.

### calculation
- **Type**: String
- **Default**: `scf`
- **Options**:
  - `scf` - Self-consistent electronic structure calculations
  - `nscf` - Non-self-consistent calculations (requires charge density file)
  - `relax` - Structure relaxation (ionic iterations)
  - `cell-relax` - Cell relaxation calculations
  - `md` - Molecular dynamics simulations
  - `get_pchg` - Obtain partial charge densities (LCAO only)
  - `get_wf` - Obtain wave functions (LCAO only)
  - `get_S` - Obtain overlap matrix (LCAO with multiple k-points)
  - `gen_bessel` - Generate Bessel projectors for DeePKS (LCAO only)
  - `test_memory` - Estimate memory consumption
  - `test_neighbour` - Obtain neighboring atom information (LCAO only)

### esolver_type
- **Type**: String
- **Default**: `ksdft`
- **Options**:
  - `ksdft` - Kohn-Sham DFT
  - `ofdft` - Orbital-free DFT
  - `sdft` - Stochastic DFT
  - `tddft` - Real-time TDDFT
  - `lj` - Lennard-Jones potential
  - `dp` - Deep Potential (DeePMD)
  - `ks-lr` - KS-DFT + LR-TDDFT
  - `lr` - LR-TDDFT with given KS orbitals

### symmetry
- **Type**: Integer
- **Default**: `0` or `1` (conditional)
- **Options**:
  - `-1` - No symmetry (recommended for non-collinear + SOC)
  - `0` - Time reversal symmetry only
  - `1` - Full symmetry analysis
- **Default Conditions**:
  - `0` if: calculation=md/nscf/get_pchg/get_wf/get_S, gamma_only=True, dft_functional=hse/hf/pbe0/scan0/opt_orb, rpa=True, or efield_flag=1
  - `1` otherwise

### symmetry_prec
- **Type**: Real
- **Default**: `1.0e-6`
- **Description**: Accuracy for symmetry judgment. Enlarge if lattice parameters or atom positions are not accurate.

### symmetry_autoclose
- **Type**: Boolean
- **Default**: `True`
- **Options**:
  - `True` - Automatically set symmetry=0 on error and continue
  - `False` - Quit with error message

### kpar
- **Type**: Integer
- **Default**: `1`
- **Description**: Divide processors into kpar groups for k-point distribution. Must be ≤ number of k-points and MPI processes.

### bndpar
- **Type**: Integer
- **Default**: `1`
- **Description**: Divide processors into bndpar groups for band distribution (stochastic orbitals). Must be > 0.

### latname
- **Type**: String
- **Default**: `none`
- **Options**: `none`, `sc`, `fcc`, `bcc`, `hexagonal`, `trigonal`, `st`, `bct`, `so`, `baco`, `fco`, `bco`, `sm`, `bacm`, `triclinic`
- **Description**: Bravais lattice type. When not `none`, lattice vectors are generated automatically.

### init_wfc
- **Type**: String
- **Default**: `atomic`
- **Options**:
  - `atomic` - From atomic pseudo wave functions
  - `atomic+random` - Add small random numbers on atomic pseudo-wavefunctions
  - `file` - From binary files `WAVEFUNC*.dat`
  - `random` - Random numbers
  - `nao` - From numerical atomic orbitals
  - `nao+random` - Add small random numbers on numerical atomic orbitals

### init_chg
- **Type**: String
- **Default**: `atomic`
- **Options**:
  - `atomic` - Summation of atomic densities
  - `file` - Read from `charge-density.dat` or cube files
  - `wfc` - Calculate from wavefunctions
  - `auto` - Try file first, fallback to atomic

### init_vel
- **Type**: Boolean
- **Default**: `False`
- **Options**:
  - `True` - Read velocity from STRU file
  - `False` - Assign Gaussian distributed random velocities

### mem_saver
- **Type**: Boolean (0/1)
- **Default**: `0`
- **Description**: Memory saving technique for nscf calculations with many k-points.

### diago_proc
- **Type**: Integer
- **Default**: `0`
- **Description**: Number of processes for diagonalization. `0` = all MPI processes.

### nbspline
- **Type**: Integer
- **Default**: `-1`
- **Description**: Order of Cardinal B-spline for Structure Factor calculation. `-1` = disabled.

### kspacing
- **Type**: Real or 3 Reals
- **Default**: `0.0`
- **Description**: Minimum k-point spacing (1/Bohr). When > 0, KPT file is unnecessary. Suggested < 0.25.
- **Format**:
  - Single value: `kspacing 0.14` (same spacing in all directions)
  - Three values: `kspacing 0.1 0.1 1.0` (different spacing for a/b/c directions)
- **Use Case**: For 2D materials with vacuum layer, set larger spacing in vacuum direction to reduce k-points:
  ```
  kspacing 0.1 0.1 1.0  # Dense sampling in-plane, 1 k-point in vacuum direction
  ```

### min_dist_coef
- **Type**: Real
- **Default**: `0.2`
- **Description**: Factor for minimum allowed interatomic distance check. Reduce for high-pressure calculations.

### device
- **Type**: String
- **Default**: `cpu`
- **Options**:
  - `cpu` - CPU computation
  - `gpu` - GPU computation (CUDA/ROCm)

### precision
- **Type**: String
- **Default**: `double`
- **Options**:
  - `single` - Single precision
  - `double` - Double precision

---

## Related References

- [Input Files](input-files.md) - STRU, KPT, INPUT file configuration
- [Electronic Structure](electronic-structure.md) - SCF, diagonalization, mixing parameters
- [Plane Wave](plane-wave.md) - PW-specific parameters
