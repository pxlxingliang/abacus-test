# DOS - Density of States / PDOS

Analyze and plot total and partial density of states from completed ABACUS calculations.

---

## Overview

The `dos-pdos` model provides **post-processing only** for DOS and PDOS analysis:

| Mode | Command | Description |
|------|---------|-------------|
| **prepare** | ❌ | Not supported - this is post-processing only |
| **post** | ❌ | Not applicable - runs directly without subcommands |
| **Direct Run** | ✅ | `abacustest model dos-pdos [args]` - Analyze and plot DOS/PDOS |

**Important:** This is a **post-processing tool** that works on **a completed ABACUS calculation** that has DOS output enabled. It does NOT require a band structure calculation - any SCF calculation with `out_dos` set will work.

---

## Prerequisites: Enable DOS Output in ABACUS

Before running calculations, you must enable DOS output in the INPUT file:

### `out_dos` Parameter

| Value | Description | Output Files |
|-------|-------------|--------------|
| **0** | No DOS output (default) | None |
| **1** | Output DOS only | `DOS1` + `DOS1_smearing.dat` |
| **2** | Output DOS + PDOS (LCAO basis only) | DOS files + `PDOS` |

**Default:** `out_dos 0` (no DOS output)


## Use When

- Plot total density of states (DOS) from completed ABACUS calculations
- Analyze orbital contributions (s, p, d, f orbitals)
- Study element-specific PDOS (projected DOS)
- Identify impurity/defect states near Fermi level
- Compare orbital contributions between elements
- Analyze chemical bonding and electronic structure

---

## Workflow: ABACUS Calculation → DOS Analysis

### Step 1: Run ABACUS Calculation with DOS Enabled

```bash
# Ensure INPUT file has out_dos set
# INPUT file should contain:
#   out_dos 1    # or 2 for PDOS

# Run your ABACUS calculation
# This can be done via abacustest or directly
abacustest submit -p setting.json

# Or run ABACUS directly
abacus | tee out.log
```

**Note:** Any ABACUS calculation type works (SCF, relax, cell-relax, MD) as long as `out_dos` is set in the INPUT file. You do NOT need to run a band structure calculation.

### Step 2: Analyze DOS/PDOS (Direct Run)

```bash
# Analyze DOS from completed calculation directory
abacustest model dos-pdos -j results/job1

# Custom energy range and plot type
abacustest model dos-pdos -j results/job1 --range -10 10 --plot-type species
```

**Important:** The `-j` parameter should point to the directory containing completed ABACUS calculation results with DOS output files (typically `results/job1/` after submit).

---

## Input Directory Structure

**Important:** The `-j` parameter accepts **ABACUS calculation directories** with completed calculations that have DOS output.

**After `abacustest submit`, results are in `results/` directory:**
```
project/
├── job1/             # Original input directory
│   ├── INPUT         # Must have out_dos set
│   ├── STRU
│   └── ...
├── setting.json       # Submission configuration
└── results/           # Directory with actual calculation outputs
    └── job1/
        ├── INPUT           # ABACUS input (with out_dos)
        ├── STRU            # Structure file
        ├── KPT             # K-points
        ├── running_scf.log # ABACUS output log
        ├── OUT.ABACUS/     # ABACUS output directory
        │   ├── DOS1    # DOS output (nspin=1)
        │   ├── PDOS      # PDOS 
        │   └── ...
        └── ...
```

---

## Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--job` | `-j` | **ABACUS calculation directories** (supports multiple for batch) | Current directory |
| `--range` | - | Energy range for plots (efermi ± eV) | `-10 10` (eV) |
| `--plot-type` | - | Plot type: `species`, `shell`, `orbital`, or `atom` | `species` |
| `--atom-index` | - | Atom index (1-based) for `atom` plot type | All atoms |
| `--suffix` | - | Suffix for output filenames | - |
| `--no-save-data` | - | Skip saving data files (.dat) | Save data |
| `--no-save-plot` | - | Skip saving plot files (.png) | Save plots |

**Plot types:**
- `species` - DOS by element (e.g., C, O, Fe)
- `shell` - DOS by orbital shell (e.g., s, p, d, f)
- `orbital` - DOS by specific orbital (e.g., p_x, d_xy)
- `atom` - DOS for specific atom (requires `--atom-index`)

**Examples:**
```bash
# Single job with default settings
abacustest model dos-pdos -j results/job1

# Multiple jobs (batch processing)
abacustest model dos-pdos -j results/job1 results/job2 results/job3

# Custom energy range (efermi ± 15 eV)
abacustest model dos-pdos -j results/job1 --range -15 15

# Plot by orbital shell (s, p, d, f)
abacustest model dos-pdos -j results/job1 --plot-type shell

# Plot specific atom (atom index 5)
abacustest model dos-pdos -j results/job1 --plot-type atom --atom-index 5

# Custom suffix for output files
abacustest model dos-pdos -j results/job1 --suffix spin_up

# Save only plots, not data files
abacustest model dos-pdos -j results/job1 --no-save-data
```

---

## Output Files

After running dos-pdos (in the directory specified by `-j`, typically `results/job1/`):

```
results/job1/
├── DOS.png            # Total DOS plot
├── DOS.dat            # Total DOS data (energy, DOS)
├── PDOS.png           # Partial DOS plot
├── PDOS.dat           # Partial DOS data
├── PDOS_{element}.dat # Per-element PDOS files
└── ...
```

**With custom suffix:**
```
results/job1/
├── DOS_spin_up.png
├── DOS_spin_up.dat
├── PDOS_spin_up.png
└── PDOS_spin_up.dat
```

---

## Tips

| Recommendation | Reason |
|----------------|--------|
| Set `out_dos 1` in INPUT | Required for DOS output |
| Use `out_dos 2` for PDOS | Requires LCAO basis (`basis_type lcao`) |
| Use dense k-mesh | Better DOS resolution and smoother curves |
| Use Gaussian smearing | Smooth DOS curves for better visualization |
| Plot PDOS by element | Identify which elements contribute to specific states |
| Check energy range | Adjust `--range` to show relevant energy window |

---

## Common Issues

| Issue | Solution |
|-------|----------|
| No DOS.dat file generated | Check if `out_dos` was set in INPUT file |
| DOS plot looks noisy | Use denser k-mesh; increase smearing width |
| Wrong Fermi level | Ensure SCF calculation converged properly |
| Missing PDOS for element | Check if `out_dos 2` and `basis_type nao` are set |
| Empty PDOS files | Verify calculation completed; check OUT.ABACUS/ for doss* files |
| "No DOS data found" error | Ensure calculation has DOS output files in OUT.ABACUS/ |

---

## Related

- Main router: [`abacustest-models`](../SKILL.md)
- Input generation: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
- Band structure: [`band.md`](band.md)
- Work function: [`workfunc.md`](workfunc.md)
- ABACUS documentation: `out_dos` parameter reference
