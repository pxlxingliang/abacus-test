# Vibration - Vibration Frequency Analysis

Calculate vibration frequencies of selected atoms using finite difference method.

---

## Overview

The `vibration` model provides a complete workflow for vibrational frequency calculations:

| Mode | Command | Description |
|------|---------|-------------|
| **prepare** | `abacustest model vibration prepare` | Generate displaced structures for frequency calculation |
| **post** | `abacustest model vibration post` | Calculate frequencies and thermochemical properties |
| **Direct Run** | ❌ | Not supported - use prepare + post workflow |

---

## Use When

- Calculate vibrational frequencies 
- Compute zero-point energy (ZPE)
- Obtain thermochemical properties (H, S, G) at various temperatures
- Study adsorbate vibrations on surfaces
- Analyze molecular vibrations and normal modes

---

## Input Directory Structure

**Important:** The `-j` parameter accepts **ABACUS calculation input directories**, not structure files.

```
project/
├── job1/
│   ├── INPUT      # ABACUS parameters (calculation = scf, force = 1)
│   ├── STRU       # Molecular/cluster structure
│   └── pp/        # Pseudopotentials
├── job2/
│   ├── INPUT
│   ├── STRU
│   └── pp/
└── ...
```

Each job directory must contain:
- Complete ABACUS inputs with force calculation enabled (`force = 1` in INPUT)
- Optimized molecular/cluster structure (forces should be near zero)

---

## Parameters

### Prepare Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--job` | `-j` | **ABACUS input directories** (supports multiple for batch) | Current directory |
| `--stepsize` | `-s` | Step size for finite difference (Å) | `0.01` Å |
| `--index` | `-i` | Atom indices to calculate (1-based, space-separated) | All atoms |
| `--image` | - | Docker image for Bohrium | `registry.dp.tech/dptech/abacus:LTSv3.10.1` |
| `--machine` | - | Machine type for Bohrium | `c32_m64_cpu` |
| `--abacus_command` | - | ABACUS execution command | `OMP_NUM_THREADS=1 mpirun -np 16 abacus \| tee out.log` |

**Example:**
```bash
# All atoms, default step size
abacustest model vibration prepare -j job1

# Specific atoms only (faster)
abacustest model vibration prepare -j job1 -i 1 2 3 4 5

# Larger step size
abacustest model vibration prepare -j job1 -s 0.02

# Multiple jobs (batch)
abacustest model vibration prepare -j job1 job22 job33

# Custom machine
abacustest model vibration prepare -j job1 --machine c64_m128_cpu
```

### Post Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--job` | `-j` | **Job directories** with vibration results (supports multiple) | Current directory |
| `--temperature` | `-t` | Temperature(s) for thermochemistry (K) | `298.15` K |

**Temperature format:**
- Single temperature: `-t 298.15` (calculate at 298.15 K)
- Temperature range: `-t 100 500 40` (40 points from 100K to 500K)

**Example:**
```bash
# Default temperature (298.15 K)
abacustest model vibration post -j job1

# Single custom temperature
abacustest model vibration post -j job1 -t 500

# Temperature range
abacustest model vibration post -j job1 -t 100 1000 10

# Multiple jobs
abacustest model vibration post -j job1 job2 job3 -t 298.15
```

---

## Output Files

After running `abacustest model vibration post`, the following files are generated:

### Directory Structure

```
results/                    # Output directory (specified by --job or save_path)
├── 000000/                 # One subdirectory per input job
│   └── vib/                # Vibration analysis results
│       ├── SCF/            # All displacement calculation outputs
│       │   ├── eq/         # Equilibrium structure calculation
│       │   ├── disp_1_x+   # Atom 1, +x displacement
│       │   ├── disp_1_x-   # Atom 1, -x displacement
│       │   ├── disp_1_y+   # Atom 1, +y displacement
│       │   ├── disp_1_y-   # Atom 1, -y displacement
│       │   ├── disp_1_z+   # Atom 1, +z displacement
│       │   └── disp_1_z-   # Atom 1, -z displacement
│       ├── cache/          # Force calculation cache files
│       │   ├── cache.eq.json
│       │   ├── cache.0x+.json
│       │   ├── cache.0x-.json
│       │   └── ...
│       └── cache.*.traj    # Trajectory files
├── metrics_vibration.json  # Main results summary
└── setting.json            # Run configuration
```

### Key Output Files

#### `metrics_vibration.json`

Main results file containing vibrational frequencies and thermochemical properties:

```json
{
    "000000/": {
        "frequencies": [
            182.11,
            182.25,
            182.30
        ],
        "zero_point_energy": 0.0339,
        "thermo_corr": {
            "298.15K": {
                "entropy": 0.00030,
                "free_energy": -0.00744
            }
        }
    }
}
```

**Fields:**
| Field | Unit | Description |
|-------|------|-------------|
| `frequencies` | cm⁻¹ | Vibrational frequencies (3N-6 for non-linear molecules) |
| `zero_point_energy` | eV | Zero-point energy (ZPE) |
| `thermo_corr` | eV, eV/K | Thermochemical corrections at specified temperatures |
| `thermo_corr.entropy` | eV/K | Entropy (S) |
| `thermo_corr.free_energy` | eV | Gibbs free energy correction (G) |


#### `vib/SCF/` Subdirectories

Each displacement directory contains complete ABACUS output:
- `INPUT`, `STRU`, `KPT` - Input files for that displacement
- `OUT.ABACUS/` - ABACUS output directory
- `out.log`, `running_scf.log` - Calculation logs
- `abacus.json` - Calculation metadata

#### `vib/cache/` Files

JSON cache files storing forces from each displacement calculation:
- `cache.eq.json` - Forces at equilibrium geometry
- `cache.0x+.json`, `cache.0x-.json` - Forces for ±x displacements of atom 1
- `cache.0y+.json`, `cache.0y-.json` - Forces for ±y displacements
- `cache.0z+.json`, `cache.0z-.json` - Forces for ±z displacements

These are used internally to compute the Hessian matrix.

---

## Tips

| Recommendation | Reason |
|----------------|--------|
| Optimize geometry first | Frequencies require minimum energy structure |
| Check for imaginary modes | Negative frequencies indicate transition state or unconverged geometry |
| Use small step size | 0.01-0.02 Å typical for accurate Hessian |
| Selective atoms for large systems | Calculate only relevant atoms (e.g., adsorbate) to save cost |
| Verify force convergence | Tight force convergence needed for accurate frequencies |

---

## Common Issues

| Issue | Solution |
|-------|----------|
| Imaginary frequencies | Geometry not fully optimized; re-optimize with tighter convergence |
| No vibration.json | Post step failed; check if all displacements completed |
| Wrong number of modes | Verify atom count; 3N-6 for molecules (3N-5 linear) |
| Force calculation failed | Ensure `force = 1` in INPUT file |
| Thermochemistry missing | Check temperature format; ensure frequencies are real |

---

## Related

- Main router: [`abacustest-models`](../SKILL.md)
- Phonon (solids): [`phonon.md`](phonon.md)
- Supercell: [`supercell.md`](supercell.md)
- Input generation: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
