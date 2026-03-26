# Conv - Convergence Tests

Automatically test convergence of DFT parameters (cutoff energy, k-point spacing, etc.).

---

## Overview

The `conv` model provides a complete workflow for convergence testing:

| Mode | Command | Description |
|------|---------|-------------|
| **prepare** | `abacustest model conv prepare` | Generate abacustest configuration for batch parameter testing |
| **post** | `abacustest model conv post` | Analyze convergence and plot results |
| **Direct Run** | ❌ | Not supported - use prepare + post workflow |

---

## Use When

- Test ecutwfc (cutoff energy) convergence
- Test k-spacing (k-point density) convergence
- Determine optimal calculation parameters
- Verify calculation accuracy before production runs
- Balance accuracy vs. computational cost

---

## Workflow: Prepare + Submit + Post

### Step 1: Prepare Convergence Test

```bash
# Test ecutwfc convergence (40, 60, 80, 100 Ry)
abacustest model conv prepare -j job1 -k ecutwfc -v 40 60 80 100

# Test kspacing convergence
abacustest model conv prepare -j job1 -k kspacing -v 0.2 0.15 0.1 0.05
```

**What `prepare` does:**
- Creates `setting.json` configuration file in the job directory
- The `setting.json` contains `mix_input` section specifying parameter values to test
- Does NOT create subdirectories yet (that happens during submit/prepare step)
- Configuration tells abacustest to run the same job multiple times with different parameter values

**Example `setting.json` generated:**
```json
{
    "save_path": "results",
    "bohrium_group_name": "convergence-test",
    "prepare": {
        "example_template": ["."],
        "mix_input": {
            "ecutwfc": ["40", "60", "80", "100"]
        }
    },
    "run_dft": {
        "command": "OMP_NUM_THREADS=1 mpirun -n 16 abacus | tee out.log",
        "image": "registry.dp.tech/dptech/abacus:LTSv3.10.1",
        "bohrium": {
            "scass_type": "c32_m64_cpu"
        }
    },
    "post_dft": {
        "metrics": {
            "dft_type": "abacus",
            "metrics_name": [
                "energy_per_atom",
                "force",
                {"ecutwfc": "{INPUT}['ecutwfc']"}
            ]
        }
    }
}
```

### Step 2: Submit to Remote Cluster (or Prepare Locally)

```bash
# Option A: Submit directly to Bohrium (creates subdirectories and submits)
abacustest submit -p setting.json

# Option B: Prepare inputs locally first (creates subdirectories without submitting)
abacustest prepare -p setting.json -s .

# After prepare, subdirectories are created:
# 00000/ (ecutwfc=40), 00001/ (ecutwfc=60), 00002/ (ecutwfc=80), 00003/ (ecutwfc=100)
```

**What `submit` or `prepare` does:**
- Reads `setting.json` and `mix_input` configuration
- Creates numbered subdirectories (00000, 00001, 00002, ...)
- Each subdirectory contains:
  - Symlinks to original STRU, KPT, pseudopotentials, orbitals
  - Modified INPUT file with the specific parameter value
- Submits calculations to Bohrium (if using `submit`)

**Directory structure after prepare:**
```
job1/
├── 00000/           # ecutwfc = 40 Ry
│   ├── INPUT        # Modified: ecutwfc = 40
│   ├── STRU         # Symlink to parent
│   ├── KPT          # Symlink to parent
│   └── pp/          # Symlinks to pseudopotentials
├── 00001/           # ecutwfc = 60 Ry
│   ├── INPUT        # Modified: ecutwfc = 60
│   └── ...
├── 00002/           # ecutwfc = 80 Ry
├── 00003/           # ecutwfc = 100 Ry
└── setting.json     # Configuration file
```

### Step 3: Post-process Results

```bash
# Analyze convergence and generate plot
abacustest model conv post -j job1 -k ecutwfc

# Include extra metrics in plot
abacustest model conv post -j job1 -k ecutwfc -e "energy_per_atom:Energy (eV/atom)"
```

---

## Input Directory Structure

**Important:** The `-j` parameter accepts **ABACUS calculation input directories**, not structure files.

```
project/
├── job1/
│   ├── INPUT          # Base ABACUS parameters
│   ├── STRU           # Atomic structure
│   ├── KPT            # K-points (optional)
│   └── pp/            # Pseudopotentials (symlinks ok)
└── ...
```

---

## Parameters

### Prepare Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--jobs` | `-j` | **ABACUS input directories** (supports multiple for batch) | Current directory |
| `--key` | `-k` | INPUT parameter to test (e.g., `ecutwfc`, `kspacing`) | Required |
| `--value` | `-v` | Parameter values to test (space-separated) | Required |

**Example:**
```bash
# Test ecutwfc: 40, 60, 80, 100 Ry (4 calculations)
abacustest model conv prepare -j job1 -k ecutwfc -v 40 60 80 100

# Test kspacing: 0.2, 0.15, 0.1, 0.05 (4 calculations)
abacustest model conv prepare -j job1 -k kspacing -v 0.2 0.15 0.1 0.05

# Multiple jobs (batch)
abacustest model conv prepare -j job1 job2 -k ecutwfc -v 40 60 80 100
```

**Common parameters to test:**
- `ecutwfc` - Plane wave cutoff energy (Ry)
- `kspacing` - K-point spacing (1/Å)
- `scf_thr` - SCF convergence threshold
- `mixing_beta` - Charge mixing parameter

### Post Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--key` | `-k` | Parameter that was tested | Required |
| `--jobs` | `-j` | **Job directories** with convergence results (supports multiple) | Current directory |
| `--metric` | `-m` | metrics.json file path | Auto-detect |
| `--extra` | `-e` | Extra metrics to plot (format: "name:title") | - |

**Example:**
```bash
# Post-process ecutwfc convergence
abacustest model conv post -j job1 -k ecutwfc

# Post-process kspacing convergence
abacustest model conv post -j job1 -k kspacing

# Include extra metric in plot
abacustest model conv post -j job1 -k ecutwfc -e "force_max:Max Force (eV/A)"

# Multiple jobs
abacustest model conv post -j job1 job2 -k ecutwfc
```

---

## Output Files

After post-processing:

```
job1/
├── 00000/           # SCF calculation with ecutwfc=40
├── 00001/           # SCF calculation with ecutwfc=60
├── ...
├── setting.json     # Configuration file (from prepare)
├── metrics.json     # Convergence data (from submit/post_dft)
├── conv_result.json # Processed convergence results
└── conv_plot.png    # Convergence curve plot
```

**metrics.json format** (generated by submit):
```json
{
    "00000": { "ecutwfc": 40, "energy_per_atom": -112144.887, "kspacing": 0.1 },
    "00001": { "ecutwfc": 60, "energy_per_atom": -112145.012, "kspacing": 0.1 },
    "00002": { "ecutwfc": 80, "energy_per_atom": -112145.089, "kspacing": 0.1 },
    "00003": { "ecutwfc": 100, "energy_per_atom": -112145.095, "kspacing": 0.1 }
}
```

**conv_result.json format** (generated by post):
```json
{
    "job1/00000": { "ecutwfc": 40, "energy_per_atom": -112144.887 },
    "job1/00001": { "ecutwfc": 60, "energy_per_atom": -112145.012 },
    "convergence": {
        "parameter": "ecutwfc",
        "converged_value": 80,
        "convergence_threshold": 0.001,
        "recommended": 100
    }
}
```

**Plots generated:**
- `conv_plot.png` - Convergence curve showing parameter value vs. metric (energy/atom)

---

## Examples

### Example 1: Ecutwfc Convergence Test

```bash
# Step 1: Prepare convergence test configuration
abacustest model conv prepare -j fe-bcc -k ecutwfc -v 40 60 80 100

# Step 2: Submit calculations (creates 00000-00003 subdirectories)
abacustest submit -p setting.json

# Step 3: Post-process and analyze
abacustest model conv post -j fe-bcc -k ecutwfc

# View convergence result
python -c "
import json
d = json.load(open('fe-bcc/metrics.json'))
print('Results:', d)
"
```

### Example 2: K-spacing Convergence Test

```bash
# Prepare: test 4 kspacing values
abacustest model conv prepare -j si-structure -k kspacing -v 0.2 0.15 0.1 0.05

# Submit
abacustest submit -p setting.json

# Post-process
abacustest model conv post -j si-structure -k kspacing
```

### Example 3: Batch Convergence Tests

```bash
# Test ecutwfc convergence for 3 materials
abacustest model conv prepare -j mat1 mat2 mat3 -k ecutwfc -v 40 60 80 100

# Submit all (each gets its own setting.json)
cd mat1 && abacustest submit -p setting.json && cd ..
cd mat2 && abacustest submit -p setting.json && cd ..
cd mat3 && abacustest submit -p setting.json && cd ..

# Post-process all
abacustest model conv post -j mat1 mat2 mat3 -k ecutwfc
```

### Example 4: Multiple Metrics in Plot

```bash
# Include both energy and force in convergence plot
abacustest model conv post -j job1 -k ecutwfc \
    -e "energy_per_atom:Energy (eV/atom)" \
    -e "force_max:Max Force (eV/A)"
```

---

## Tips

| Recommendation | Reason |
|----------------|--------|
| Test ecutwfc first | Affects all subsequent calculations |
| Then test k-spacing | K-point density significantly impacts accuracy |
| Use reasonable ranges | Don't waste resources on extreme values |
| Check energy convergence | Typical threshold: 1 meV/atom |
| Test on representative structure | Use primitive cell for faster testing |

**Typical convergence criteria:**
- Energy: < 1 meV/atom change
- Force: < 0.01 eV/Å change
- Stress: < 0.1 GPa change

---

## Common Issues

| Issue | Solution |
|-------|----------|
| No subdirectories created | Run `abacustest submit -p setting.json` or `abacustest prepare -p setting.json` after model conv prepare |
| No metrics.json | Calculations not completed; check submission status |
| Plot looks strange | Verify parameter values are in correct range |
| Calculations not converged | Increase SCF iterations or adjust mixing |
| Wrong parameter tested | Ensure `-k` matches INPUT parameter name exactly |

---

## Related

- Main router: [`abacustest-models`](../SKILL.md)
- Input generation: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
- EOS: [`eos.md`](eos.md)
- Band structure: [`band.md`](band.md)
