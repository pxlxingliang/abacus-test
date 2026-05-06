# EOS - Equation of State

Calculate equation of state by generating volume-scaled structures, fitting E-V curves, and extracting bulk modulus and equilibrium properties.

---

## Overview

The `eos` model provides a complete workflow for equation of state calculations:

| Mode | Command | Description |
|------|---------|-------------|
| **prepare** | `abacustest model eos prepare` | Generate volume-scaled structures for EOS calculation |
| **post** | `abacustest model eos post` | Fit E-V curves and extract bulk modulus, equilibrium volume |
| **Direct Run** | ❌ | Not supported - use prepare + post workflow |

---

## Use When

- Calculate bulk modulus (B₀) and its derivative (B₀')
- Determine equilibrium lattice constant and volume
- Fit E-V curves using Birch-Murnaghan equation
- Compare structural stability between phases
- Generate equation of state data for thermodynamic modeling

---

## Workflow: Prepare + Submit + Post

### Step 1: Prepare EOS Calculations

```bash
# Generate volume-scaled structures (default: 9 points from 0.9 to 1.1)
abacustest model eos prepare -j job1 job2 job3

# Custom volume range and step
abacustest model eos prepare -j job1 -s 0.92 -e 1.08 -d 0.02
```

**What `prepare` does:**
- Creates volume-scaled subdirectories (eos-0.900, eos-0.925, etc.)
- Generates INPUT and STRU files for each volume point
- Prepares submission configuration for batch computing

### Step 2: Submit to Remote Cluster

```bash
# Review generated files, then submit (after user confirmation)
abacustest submit -p setting.json
```

**What happens during submit:**
- Creates `results/` directory (or path specified by `save_path` in setting.json)
- Submits calculations to Bohrium cluster
- Upon completion, calculation outputs are stored in `results/<job_id>/`

### Step 3: Post-process Results

```bash
cd results # cd to results directory

# Post-process from results directory after calculations complete
abacustest model eos post -j job1 -s result.json

# Generate combined plot for multiple jobs
abacustest model eos post -j job1 job2 --sumfile all_eos.png
```

**Important:** After `abacustest submit`, the calculation results are in the `results/` directory, not the original job directory. Change directory to results or Use `-j results/<job_id>` for post-processing.

---

## Input Directory Structure

**Important:** The `-j` parameter accepts **ABACUS calculation input directories**, not structure files.

```
project/
├── job1/
│   ├── INPUT      # ABACUS parameters (calculation = scf)
│   ├── STRU       # Atomic structure
│   └── pp/        # Pseudopotentials
├── job2/
│   ├── INPUT
│   ├── STRU
│   └── pp/
└── job3/
    └── ...
```

Each job directory must contain complete ABACUS inputs. The prepare step will create `eos-*/` subdirectories inside each job.

---

## Parameters

### Prepare Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--job` | `-j` | **ABACUS input directories** (supports multiple for batch) | Current directory |
| `--start` | `-s` | Minimum volume scaling factor (< 1) | `0.9` |
| `--end` | `-e` | Maximum volume scaling factor (> 1) | `1.1` |
| `--step` | `-d` | Volume step size | `0.025` |
| `--kspacing` | - | Use kspacing for KPT generation (0=no, 1=yes) | `0` |
| `--relax` | - | Calculation type: -1=from INPUT, 0=scf, 1=relax, 2=cell-relax | `-1` |
| `--clean` | - | Clean eos folder before running (0=no, 1=yes) | `1` |

**Number of volume points:** `(end - start) / step + 1`

**Example:**
```bash
# 9 volume points: 0.90, 0.925, 0.95, ..., 1.10
abacustest model eos prepare -j job1 -s 0.9 -e 1.1 -d 0.025

# 17 volume points with finer sampling
abacustest model eos prepare -j job1 -s 0.9 -e 1.1 -d 0.0125
```

### Post Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--job` | `-j` | **Job directories** with eos subdirectories (supports multiple) | Current directory |
| `--save` | `-s` | Output JSON filename for EOS results | `result.json` |
| `--result` | `-r` | Input result file with volume and energy data | - |
| `--ref` | - | Reference result file for delta comparison | - |
| `--sumfile` | - | Combined EOS plot filename for all jobs | `alleos.png` |

**Example:**
```bash
# Post-process single job
abacustest model eos post -j job1 -s job1-result.json

# Post-process multiple jobs with combined plot
abacustest model eos post -j job1 job2 job3 --sumfile comparison.png

# Compare with reference data
abacustest model eos post -j job1 --ref reference-result.json
```

---

## Output Files

**After submit (calculation complete):**
```
project/
├── job1/             # Original input subdirectories
│   ├── eos-0/
│   ├── eos-1/
│   └── ...
├── setting.json       # Submission configuration
├── metrics.json       # Extracted results from all calculations
└── results/           # Directory with actual calculation outputs
    └── job1/
        ├── eos-0/     # Completed SCF calculations
        │   ├── INPUT
        │   ├── STRU
        │   ├── OUT.ABACUS/
        │   └── running_scf.log
        ├── eos-1/
        └── ...
```

**After post-processing:**
```
results/
├── result.json        # volume & energy arrays
├── metrics.json        # key results of each sub-jobs 
├── plot_data.csv        # key results of each jobs
└── job1.png        # E-V curve plot (if regex module installed)
```

**metrics.json format** :
```json
{
    "job1/eos-0": {
        "volume": 190.124,
        "energy_per_atom": -1144.2915,
        "natom": 8,
        "converge": true,
        "normal_end": true,
        "scf_steps": 8,
        "total_time": 23.69
    }
}
```

**result.json format** (generated by post):
```json
{
    "job1": {
        "volume": [185.37, 180.62, ...],
        "energy_per_atom": [-1144.29, -1144.28, ...]
    }
}
```

**result.json format:**
```json
{
  "job1": { "volume": [112.5,...], "energy_per_atom": [-112144.887,...] },
  "job2": { "volume": [115.6,...], "energy_per_atom": [-112145.012,...] },
}
```

**plot_data.csv format** (generated by post):
```
Example job1:
x,185.37,180.62,175.87,171.11,190.12,194.88,199.63,204.38,209.14
y,-1144.288974,-1144.280261,-1144.264167,-1144.239419,-1144.291504,-1144.288838,-1144.281593,-1144.270159,-1144.254994
fitx,171.11,171.50,...,208.75,209.14
fity,-1144.239499,-1144.241826,...,-1144.256433,-1144.255112
v0,189.94
e0,-1144.291670
b0,7.295308
bp,4.894958
res,1.124934e-04
```
x/y: the calculated volume and energy
fitx/fity: the fitted volume and energy
v0: the fitted volume with minimum energy
e0: the fitted minimum energy
b0/bp: the fitted bulk_modulus and its derivative
res: the residuals between calculated and fitted values


**Plots generated:**
- `eos_fit.png` - E-V curve with Birch-Murnaghan fit
- `alleos.png` - Combined E-V curves for multiple jobs (if `--sumfile` used)

---

## Examples

### Example 1: Bulk Modulus of BCC Iron

```bash
# Prepare EOS calculation (9 volume points)
abacustest model eos prepare -j fe-bcc -s 0.92 -e 1.08 -d 0.02

# Submit after user confirmation
abacustest submit -p setting.json

# Wait for calculations to complete (results stored in results/ directory)

# change directory to results
cd rersults

# Post-process from results directory
abacustest model eos post -j fe-bcc -s fe-result.json
```

### Example 2: Compare BCC vs FCC Stability

```bash
abacustest model eos prepare -j fe-bcc fe-fcc
abacustest submit -p setting.json
abacustest model eos post -j fe-bcc fe-fcc -s result.json --sumfile result-bcc-fcc.png
```

---

## Tips

| Recommendation | Reason |
|----------------|--------|
| Volume range: ±8-12% | Too narrow → poor fit; too wide → anharmonic effects |
| At least 7-9 points | Better fitting accuracy |
| Use consistent KPT | Same k-point sampling across all volumes |
| Check convergence | Especially for compressed volumes (higher pressure) |
| Use `--relax 1` | Allow atomic positions to relax at each volume |

---

## Common Issues

| Issue | Solution |
|-------|----------|
| Volume range too small (<5%) | Use ±10% range for better fit |
| Step too large (>5%) | Use 2-3% steps for smoother curve |
| Some volumes not converged | Re-run unconverged points individually |
| Forgot to run post | Run `abacustest model eos post -j job_dir` |
| Negative bulk modulus | Check if volume range is appropriate; may indicate structural instability |

---

## Related

- Main router: [`abacustest-models`](../SKILL.md)
- Elastic constants: [`elastic.md`](elastic.md)
- Convergence tests: [`conv.md`](conv.md)
- Input generation: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
