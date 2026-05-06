# WorkFunc - Work Function

Calculate surface work function by analyzing electrostatic potential profile and vacuum level.

---

## Overview

The `workfunc` model provides a complete workflow for work function calculations:

| Mode | Command | Description |
|------|---------|-------------|
| **prepare** | `abacustest model workfunc prepare` | Prepare work function calculation with vacuum direction detection |
| **post** | `abacustest model workfunc post` | Extract work function from electrostatic potential |
| **Direct Run** | ❌ | Not supported - use prepare + post workflow |

---

## Use When

- Calculate surface work function (Φ = E_vacuum - E_Fermi)
- Analyze planar-averaged electrostatic potential profile
- Study catalytic activity and surface reactivity
- Determine vacuum level position
- Compare work functions between surfaces
- Analyze adsorption effects on work function

---

## Workflow: Prepare + Submit + Post

### Step 1: Prepare Work Function Calculation

```bash
# Prepare with automatic vacuum direction detection
abacustest model workfunc prepare -j job1 job2 --vacuum auto

# Specify vacuum direction manually
abacustest model workfunc prepare -j job1 --vacuum c --dipole-corr
```

**What `prepare` does:**
- Detects vacuum direction from structure (or uses specified direction)
- Sets up electrostatic potential calculation
- Configures dipole correction for asymmetric slabs (if enabled)
- Prepares submission configuration for remote computing
- Supports both ABACUS and VASP calculations

### Step 2: Submit to Remote Cluster

```bash
# Review generated files, then submit (after user confirmation)
abacustest submit -p setting.json
```

### Step 3: Post-process Results

```bash
# Extract work function and plot potential profile
abacustest model workfunc post -j job1 job2 --dft abacus

# Specify vacuum direction manually
abacustest model workfunc post -j job1 --vacuum c --dft abacus
```

---

## Input Directory Structure

**Important:** The `-j` parameter accepts **ABACUS calculation input directories**, not structure files.

```
project/
├── job1/
│   ├── INPUT      # ABACUS parameters (calculation = scf, out_pot = 1)
│   ├── STRU       # Surface structure with vacuum layer
│   └── pp/        # Pseudopotentials
├── job2/
│   ├── INPUT
│   ├── STRU
│   └── pp/
└── ...
```

Each job directory must contain:
- Complete ABACUS inputs
- Surface/slab structure with sufficient vacuum (>15 Å recommended)
- Electrostatic potential output enabled (`out_pot = 1` in INPUT)

---

## Parameters

### Prepare Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--job` | `-j` | **ABACUS input directories** (supports multiple for batch) | Current directory |
| `--vacuum` | - | Vacuum direction: `a`, `b`, `c`, or `auto` | `auto` |
| `--dipole-corr` | - | Enable dipole correction for asymmetric slabs | Disabled |
| `--image` | - | Docker image for Bohrium | Auto-detect |
| `--machine` | - | Machine type for Bohrium | Auto-detect |
| `--dft-command` | - | DFT execution command | Auto-detect |
| `--dft` | - | DFT software: `abacus` or `vasp` | `abacus` |
| `--potcar_path` | - | VASP POTCAR path (for VASP mode) | - |

**Example:**
```bash
# Single job with auto vacuum detection
abacustest model workfunc prepare -j pt-111

# Multiple jobs (batch)
abacustest model workfunc prepare -j pt-111 au-111 ag-111

# Manual vacuum direction with dipole correction
abacustest model workfunc prepare -j asymmetric-slab --vacuum c --dipole-corr

# VASP mode
abacustest model workfunc prepare -j vasp-job --dft vasp --potcar_path /path/to/potcars
```

### Post Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--job` | `-j` | **Job directories** with work function results (supports multiple) | Current directory |
| `--vacuum` | - | Vacuum direction: `a`, `b`, `c`, or `auto` | `auto` |
| `--dft` | - | DFT software: `abacus` or `vasp` | `abacus` |

**Example:**
```bash
# Post-process single job
abacustest model workfunc post -j job1

# Post-process multiple jobs
abacustest model workfunc post -j job1 job2 job3

# Specify vacuum direction
abacustest model workfunc post -j job1 --vacuum c

# VASP results
abacustest model workfunc post -j vasp-job --dft vasp
```

---

## Output Files

After post-processing:

```
job1/
├── work_function.json    # Work function results
├── potential.png         # Electrostatic potential profile plot
├── potential.dat         # Raw potential data (distance, potential)
└── run.sh                # Execution script (from prepare)
```

**work_function.json format:**
```json
{
  "work_function": 4.85,
  "vacuum_level": 8.23,
  "fermi_energy": 3.38,
  "vacuum_direction": "c",
  "dipole_corrected": false
}
```

**Plots generated:**
- `potential.png` - Planar-averaged electrostatic potential vs. distance along vacuum direction
  - Shows vacuum level plateau
  - Marks Fermi level position
  - Indicates work function value

**potential.dat format:**
```
# Distance(A)  Potential(eV)
0.00   5.23
0.10   5.25
...
15.00  8.23   # Vacuum level
...
30.00  5.21
```

---

## ⚠️ Important: Before Calculation

### Vacuum Layer Requirements

- **Minimum vacuum**: 15 Å recommended
- **Purpose**: Avoid periodic image interactions
- **Check**: Ensure potential reaches plateau in vacuum region

### Dipole Correction

- **When to use**: Asymmetric slabs (adsorption on one side only)
- **Purpose**: Remove artificial electric field from periodic boundary
- **Enable**: `--dipole-corr` in prepare step

### Confirm with User

**Before submitting:**
```bash
# 1. Prepare
abacustest model workfunc prepare -j surface_struct

# 2. Confirm with user
# "Work function calculation ready.
#  Vacuum direction: c (auto-detected)
#  Dipole correction: enabled/disabled
#  Submitting to Bohrium will incur costs. Proceed?"

# 3. Submit after confirmation
abacustest submit -p setting.json
```

---

## Tips

| Recommendation | Reason |
|----------------|--------|
| Vacuum >15 Å | Avoid periodic image interaction |
| Use dipole correction | For asymmetric slabs (one-sided adsorption) |
| Check potential plateau | Ensure vacuum level is well-defined |
| Verify vacuum direction | Auto-detection usually correct, but verify |
| Converge k-mesh | Work function sensitive to Fermi level accuracy |

---

## Common Issues

| Issue | Solution |
|-------|----------|
| No vacuum plateau | Increase vacuum layer thickness (>15 Å) |
| Wrong vacuum direction | Specify manually: `--vacuum c` |
| Work function negative | Check Fermi level and vacuum level alignment |
| Potential oscillates | Increase vacuum; check structure |
| Dipole correction fails | Ensure asymmetric slab; check INPUT settings |

---

## Related

- Main router: [`abacustest-models`](../SKILL.md)
- DOS/PDOS: [`dos.md`](dos.md)
- Input generation: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
- Band structure: [`band.md`](band.md)
