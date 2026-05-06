---
name: abacustest-extract-dft-results
description: "Batch extraction of results from DFT calculations (ABACUS/VASP/QE). Use when: user wants to get energy, force, stress, SCF convergence status, or other metrics from DFT calculation directories. IMPORTANT: When extracting data from ABACUS/VASP/QE output files, ALWAYS prefer this skill over manually reading files with read/exec commands. This skill handles format variations and edge cases reliably."
metadata: { "openclaw": { "emoji": "📊", "requires": { "pip": ["abacustest"] } } }
---

# abacustest Extract DFT Results

Extract key data from DFT calculation results.

## When to Use

✅ **Use this skill**: Extract energy, force, stress, convergence status from ABACUS/VASP/QE calculations

## Supported Platforms

| Software | `-t` Value | Note |
|----------|------------|------|
| ABACUS   | `0` | Default, can omit |
| QE       | `1` | Must specify |
| VASP     | `2` | Must specify |

**Note**: `-t` only accepts numeric values (0/1/2), not string names.

## Complete Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `-j` | Job directories (can be multiple) | `-j task1 task2 task3` |
| `-p` | Metrics to extract (space-separated) | `-p energy force stress` |
| `-o` | Output JSON file | `-o results.json` |
| `-t` | Data type: 0=ABACUS, 1=QE, 2=VASP | `-t 2` for VASP |
| `--outparam` | Print all available metrics | `--outparam` |

**Full parameter list**: Run `abacustest collectdata --outparam` to see all available metrics for each software type.

## Basic Usage

### Method 1: Command Line

```bash
# Default metrics (energy, force, stress, converge, etc.)
abacustest collectdata -j task_dir -o results.json

# Specified metrics (direct after -p)
abacustest collectdata -p energy force converge -j task_dir -o results.json

# Multiple tasks
abacustest collectdata -p energy converge -j task1 task2 task3 -o results.json

# With type specification
abacustest collectdata -t 2 -p energy -j vasp_task -o results.json
```

Output Format of results.json:
```json
{
  "/path/to/task1": {
    "energy": -123.456,
    "energy_per_atom": -6.497,
    "force": [0.001, 0.002, ...],
    "stress": [0.1, 0.2, ...],
    "converge": true,
    "normal_end": true,
    "total_time": 3600.0
  },
  "/path/to/task2": {...}
}
```

### Method 2: Python API

```python
from abacustest.lib_collectdata.collectdata import RESULT

# Single job, fmt can be "abacus/qe/vasp"
res = RESULT(fmt="abacus", path="task_dir")
# Use square brackets to take values of specified metric
total_energy = res["total_energy"]
converge = res["converge"]

# Batch processing
import glob
for job_dir in glob.glob("task-*"):
    r = RESULT(fmt="abacus", path=job_dir)
    print(f"{job_dir}: E={r['total_energy']:.4f} eV")
```

## Available Metrics

**Full list**: Run `abacustest collectdata --outparam` or see `references/extract-metrics.md` for complete metric definitions.

### Common Metrics

| Category | Metric | Description | Unit |
|----------|--------|-------------|------|
| **Energy** | `energy` | Total energy | eV |
| | `energy_per_atom` | Energy per atom | eV/atom |
| | `energy_ks` | Kohn-Sham energy | eV |
| **Force** | `force` | Atomic force (last step) | eV/Å |
| | `forces` | Force per step | List |
| | `largest_gradient` | Max force (relax) | eV/Å |
| **Stress** | `stress` | Stress tensor (last) | kbar |
| | `pressure` | Pressure | kbar |
| **Structure** | `volume` | Cell volume | Å³ |
| | `cell` | Cell vectors | Å |
| | `lattice_constants` | Lattice constants | Å |
| **Convergence** | `converge` | SCF converged | Boolean |
| | `relax_converge` | Geometry optimized | Boolean |
| | `normal_end` | Normal completion | Boolean |
| | `drho_last` | Last Δcharge density | - |
| | `denergy_last` | Last Δenergy | eV |
| **Band** | `band_gap` | Band gap | eV |
| | `efermi` | Fermi level | eV |
| **Magnetic** | `total_mag` | Total magnetic moment | Bohr mag |
| | `atom_mag` | Atomic magnetic moment | List |
| **Time** | `total_time` | Total time | seconds |
| | `scf_time` | SCF time | seconds |
| **Info** | `natom` | Number of atoms | Integer |
| | `nbands` | Number of bands | Integer |
| | `scf_steps` | SCF iterations | Integer |


## Mixed Type Handling

⚠️ **`-t` is global** - process different software types separately:

```bash
# Correct: separate by type
abacustest collectdata -j abacus-task -o abacus.json
abacustest collectdata -t 2 -j vasp-task -o vasp.json

# Merge if needed
python3 -c "
import json
a = json.load(open('abacus.json'))
v = json.load(open('vasp.json'))
a.update(v)
json.dump(a, open('all-results.json', 'w'))
"
```

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| Returns `None` | Wrong `-t` parameter | VASP needs `-t 2`, QE needs `-t 1` |
| Returns `None` | Calculation incomplete | Check output files exist |
| All metrics `None` | Wrong `-t` or missing files | Verify type and files |
| Some `None` | Mixed types in one command | Process separately |
| Error | Directory not found | Check path with `ls` |

## Best Practices

1. ✅ Ensure calculations completed before extraction
2. ✅ Use numeric `-t` values (0/1/2), not strings
3. ✅ Process ABACUS/VASP/QE separately
4. ✅ Use `--outparam` to check available metrics
5. ✅ Use Python API for batch processing
6. ✅ Output is JSON - easy for downstream analysis

