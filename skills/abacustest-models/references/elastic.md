# Elastic - Elastic Constants

Calculate full elastic tensor and derive mechanical properties (bulk, shear, Young's modulus, Poisson's ratio).

---

## Overview

The `elastic` model provides a complete workflow for elastic constant calculations:

| Mode | Command | Description |
|------|---------|-------------|
| **prepare** | `abacustest model elastic prepare` | Generate strained structures for elastic tensor calculation |
| **post** | `abacustest model elastic post` | Extract elastic tensor and calculate mechanical properties |
| **Direct Run** | ❌ | Not supported - use prepare + post workflow |

---

## Use When

- Calculate full elastic tensor (C_ij, 6×6 Voigt notation)
- Derive bulk modulus (B), shear modulus (G), Young's modulus (E)
- Calculate Poisson's ratio (ν)
- Analyze mechanical stability by crystal system
- Compare mechanical properties between materials
- Study elastic anisotropy

---

## Workflow: Prepare + Submit + Post

### Step 1: Prepare Elastic Calculations

```bash
# Generate strained structures (13-17 deformations)
abacustest model elastic prepare -j job1 job2 --norm 0.01 --shear 0.01
```

**What `prepare` does:**
- Applies normal strains (±0.5%, ±1%) along lattice vectors
- Applies shear strains (±0.5%, ±1%) for off-diagonal tensor elements
- Generates 24 strained structures depending on crystal symmetry
- Creates INPUT and STRU files for each deformed structure
- Prepares submission configuration for batch computing

### Step 2: Submit to Remote Cluster

```bash
# Review generated files, then submit (after user confirmation - expensive!)
abacustest submit -p setting.json
```

**What happens during submit:**
- Creates `results/` directory (or path specified by `save_path` in setting.json)
- Submits original and deformed structure calculations to Bohrium cluster
- Upon completion, calculation outputs are stored in `results/<job_id>/deformed_*/`

### Step 3: Post-process Results

```bash
# Option A: Change to results directory (recommended)
cd results
abacustest model elastic post -j job1

# Option B: Use full path from project directory
abacustest model elastic post -j results/job1

# Post-process multiple jobs
abacustest model elastic post -j job1 job2 job3
```

**Important:** After `abacustest submit`, the calculation results are in the `results/` directory. Either change to the `results/` directory first, or use the full path `results/job1` for post-processing.

---

## Input Directory Structure

**Important:** The `-j` parameter accepts **ABACUS calculation input directories**, not structure files.

```
project/
├── job1/
│   ├── INPUT      # ABACUS parameters (calculation = scf, stress = 1)
│   ├── STRU       # Atomic structure
│   ├── KPT        # K-points (optional, can be in INPUT)
│   └── pp/        # Pseudopotentials (or symlinks)
├── job2/
│   ├── INPUT
│   ├── STRU
│   └── pp/
└── ...
```

Each job directory must contain complete ABACUS inputs with `stress = 1` in the INPUT file.

---

## Parameters

### Prepare Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--job` | `-j` | **ABACUS input directories** (supports multiple for batch) | Current directory |
| `--norm` | - | Maximum strain for normal mode | `0.01` (1%) |
| `--shear` | - | Maximum strain for shear mode | `0.01` (1%) |
| `--norelax` | - | Skip atomic relaxation for deformed structures | Relax enabled |
| `--image` | - | Docker image for Bohrium | `registry.dp.tech/dptech/abacus:LTSv3.10.1` |
| `--machine` | - | Machine type for Bohrium | `c32_m64_cpu` |
| `--abacus_command` | - | ABACUS execution command | `OMP_NUM_THREADS=1 mpirun -np 16 abacus \| tee out.log` |

**Example:**
```bash
# Single job with default 1% strain
abacustest model elastic prepare -j job1

# Multiple jobs with custom strain
abacustest model elastic prepare -j job1 job2 --norm 0.015 --shear 0.015

# Skip relaxation (faster, but less accurate)
abacustest model elastic prepare -j job1 --norelax

# Custom machine for expensive calculations
abacustest model elastic prepare -j job1 --machine c64_m128_cpu
```

### Post Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--job` | `-j` | **Job directories** with elastic calculation results (supports multiple) | Current directory |

**Example:**
```bash
# Post-process single job
abacustest model elastic post -j job1

# Post-process multiple jobs
abacustest model elastic post -j job1 job2 job3
```

---

## Output Files

**After submit (calculation complete):**
```
project/
├── job1/             # Original input directory
│   ├── deformed_00/   # Deformed structure inputs
│   ├── deformed_01/
│   ├── ...
│   ├── deformed_23/
│   └── org/
├── setting.json       # Submission configuration
└── results/           # Directory with actual calculation outputs
    ├── job1/          # Completed calculations for job1
    │   ├── deformed_00/
    │   ├── deformed_01/
    │   ├── ...
    │   ├── deformed_23/
    │   └── org/
    ├── job2/          # Completed calculations for job2
    └── ...
```

**After post-processing:**
```
results/
├── metrics.json              # Key results for all calculation tasks
├── metrics_elastic.json      # Full elastic tensor and mechanical properties
├── metrics_elastic.csv       # CSV format of elastic results
└── metrics_elastic_*__*.csv  # Additional CSV files (elastic tensor, etc.) for each calculation task
```

---

## Checking Results: metrics.json

Before analyzing elastic properties, **always check `metrics.json`** to verify all calculations completed successfully:

```bash
# Check if all calculations converged
cat results/metrics.json | python -c "
import json, sys
data = json.load(sys.stdin)
for job, subjobs in data.items():
    for subjob, metrics in subjobs.items():
        converge = metrics.get('converge', False)
        normal_end = metrics.get('normal_end', False)
        status = '✓' if (converge and normal_end) else '✗'
        print(f'{status} {job}/{subjob}: converge={converge}, normal_end={normal_end}')
"
```

**metrics.json format:**
```json
{
    "job1/org": {
        "normal_end": true,
        "converge": true,
        "scf_steps": 9,
        "relax_steps": 1,
        "relax_converge": true,
        "denergy_last": -2.5e-06,
        "drho_last": 8.85e-08,
        "energy": -9154.332,
        "stress": [-0.917, 0, 0, 0, -0.918, 0, 0, 0, -0.917]
    },
    "job1/deformed_00": {
        "normal_end": true,
        "converge": true,
        "scf_steps": 9,
        "energy": -9154.327,
        "stress": [8.59, 0, 0, 0, 3.0, 0, 0, 0, 3.0]
    },
    "job1/deformed_01": { ... }
}
```

**Key fields to check:**
- `converge`: SCF calculation converged (should be `true` for all)
- `normal_end`: Calculation completed normally (should be `true` for all)
- `relax_converge`: Atomic relaxation converged (should be `true` if relaxation enabled)
- `scf_steps`: Number of SCF iterations (watch for unusually high values)

**If any calculation shows `converge: false` or `normal_end: false`:**
- That data point may be unreliable
- Elastic tensor fitting may be inaccurate
- Consider re-running failed calculations before proceeding

---

## Elastic Results: metrics_elastic.json

**metrics_elastic.json format:**
```json
{
    "job1/": {
        "elastic_tensor": [
            [92.69, 37.04, 37.05, 0, 0, 0],
            [37.27, 92.90, 37.26, 0, 0, 0],
            [37.14, 37.14, 92.78, 0, 0, 0],
            [0, 0, 0, 49.81, 0, 0],
            [0, 0, 0, 0, 49.83, 0],
            [0, 0, 0, 0, 0, 49.82]
        ],
        "bulk_modulus": 55.70,
        "shear_modulus": 41.03,
        "young_modulus": 98.82,
        "poisson_ratio": 0.204
    }
}
```

**Key mechanical properties:**
- `elastic_tensor`: 6×6 elastic stiffness tensor C_ij (GPa, Voigt notation)
- `bulk_modulus`: Bulk modulus B (GPa)
- `shear_modulus`: Shear modulus G (GPa)
- `young_modulus`: Young's modulus E (GPa)
- `poisson_ratio`: Poisson's ratio ν

**Note:** Post-processing requires the `pymatgen` Python module. Output files are generated in the `results/` directory.


## Examples

### Example 1: Elastic Constants Calculation

```bash
# Prepare elastic calculation 
abacustest model elastic prepare -j job1 --norm 0.01 --shear 0.01

# Submit after explicit user confirmation (expensive!)
abacustest submit -p setting.json

# Wait for calculations to complete (results in results/ directory)

# Change to results directory for post-processing
cd results

# Post-process from results directory
abacustest model elastic post -j job1

# Step 1: Check metrics.json for convergence
python -c "
import json
data = json.load(open('metrics.json'))
for job, metrics in data.items():
    converge = metrics.get('converge', False)
    normal_end = metrics.get('normal_end', False)
    status = '✓' if (converge and normal_end) else '✗'
    print(f'{status} {job}: converge={converge}, normal_end={normal_end}')
"

# Step 2: If all converged, view elastic results
python -c "
import json
d = json.load(open('metrics_elastic.json'))
for job, results in d.items():
    print(f'{job}:')
    print(f'  Bulk modulus: {results[\"bulk_modulus\"]:.2f} GPa')
    print(f'  Shear modulus: {results[\"shear_modulus\"]:.2f} GPa')
    print(f'  Young\\'s modulus: {results[\"young_modulus\"]:.2f} GPa')
    print(f'  Poisson\\'s ratio: {results[\"poisson_ratio\"]:.3f}')
"
```

### Example 2: Batch Processing Multiple jobs

```bash
# Prepare elastic calculations for 3 materials
abacustest model elastic prepare -j job1 job2 job3

# Submit all
abacustest submit -p setting.json

# Change to results directory
cd results

# Post-process all
abacustest model elastic post -j job1 job2 job3

# Compare results
python -c "
import json
d = json.load(open('metrics_elastic.json'))
print('Material | Bulk (GPa) | Shear (GPa) | Young (GPa) | Poisson')
print('-' * 60)
for job, r in d.items():
    print(f'{job:8} | {r[\"bulk_modulus\"]:10.2f} | {r[\"shear_modulus\"]:9.2f} | {r[\"young_modulus\"]:9.2f} | {r[\"poisson_ratio\"]:.3f}')
"
```

### Example 3: Larger Strain for Soft Materials

```bash
# Use 2% strain for better signal in soft materials
abacustest model elastic prepare -j job1 --norm 0.02 --shear 0.02

# Submit and post-process as usual
abacustest submit -p setting.json
cd results && abacustest model elastic post -j job1
```

### Example 4: Skip Relaxation for Speed

```bash
# Faster calculation (no atomic relaxation), less accurate
abacustest model elastic prepare -j job1 --norelax

# Submit and post-process
abacustest submit -p setting.json
cd results && abacustest model elastic post -j job1
```

### Example 5: Checking Convergence Before Analysis

```bash
# After post-processing, always check convergence first
cd results

# Check all calculations converged
python -c "
import json
data = json.load(open('metrics.json'))
failed = []
for job, metrics in data.items():
    if not (metrics.get('converge', False) and metrics.get('normal_end', False)):
        failed.append(job)
        
if failed:
    print('Warning: The following calculations did not converge:')
    for job in failed:
        print(f'  - {job}')
    print('Elastic results may be unreliable!')
else:
    print('All calculations converged successfully.')
"

# Only proceed with elastic analysis if all converged
```

---

## ⚠️ Important: Computational Cost

**Elastic calculations are computationally expensive:**

- **24 strained structures** must be calculated
- Each structure requires stress tensor calculation
- Total cost: 24× single SCF calculation
- **Always confirm with user before submitting!**

**Correct workflow:**
```bash
# 1. Prepare
abacustest model elastic prepare -j struct_dir

# 2. Confirm with user
# "Elastic calculation ready: 24 strained structures.
#  This is computationally expensive and will incur significant costs on Bohrium.
#  Estimated cost: [provide estimate based on machine type].
#  Proceed with submission?"

# 3. Submit after explicit confirmation
abacustest submit -p setting.json
```

---

## Tips

| Recommendation | Reason |
|----------------|--------|
| Use small strains (1-2%) | Stay in linear elastic regime |
| Check mechanical stability | Some materials are intrinsically unstable |
| Enable stress in INPUT | Required: `stress = 1` in INPUT file |
| Use higher precision | Elastic constants sensitive to convergence |
| Do cell-relaxation before elastic | Stability is better with relaxed cell |

---

## Common Issues

| Issue | Solution |
|-------|----------|
| Stress not calculated | Ensure `stress = 1` in INPUT file |
| Negative bulk modulus | Material may be mechanically unstable |
| Poor linear fit | Use smaller strain magnitude (0.5-1%) |
| Calculation too expensive | Use `--norelax` for faster but less accurate results |

---

## Related

- Main router: [`abacustest-models`](../SKILL.md)
- EOS (alternative B₀): [`eos.md`](eos.md)
- Input generation: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)

