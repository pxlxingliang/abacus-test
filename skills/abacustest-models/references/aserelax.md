# ASERelax - ASE Relaxation

Perform structure relaxation using ASE optimizers with ABACUS as the DFT calculator.

---

## Overview

The `aserelax` model provides three different usage modes:

| Mode | Command | Description |
|------|---------|-------------|
| **prepare** | `abacustest model aserelax prepare` | Prepare ASE+ABACUS relaxation scripts for a batch of calculations. Typically followed by `abacustest submit` for remote high-throughput computing on Bohrium. |
| **post** | `abacustest model aserelax post` | Post-process and plot results from optimized structures. Also supports plotting for native ABACUS relax tasks. |
| **Direct Run** | `abacustest model aserelax [args]` | Directly call ASE for relaxation. Requires `ase` and `ase-abacus` installed in the current environment. |

---

## Use When

- Flexible relaxation workflows with ASE optimizers
- High-throughput batch relaxation on remote clusters
- Use specific ASE optimizers not available in ABACUS
- Custom relaxation constraints
- Post-processing ABACUS native relax results
- Local relaxation (direct mode, requires ASE installed)

---

## Workflow 1: Prepare + Submit (Remote High-Throughput)

This is the **recommended workflow** for batch calculations on Bohrium or other remote clusters.

```bash
# Step 1: Prepare ASE relaxation scripts for all structures
abacustest model aserelax prepare -j job1 job2 job3

# Step 2: Submit to remote cluster (confirm with user first!)
abacustest submit -p setting.json

# Step 3: Post-process after calculations complete
abacustest model aserelax post -j job1 job2 job3
```

**What `prepare` does:**
- Generates Python scripts for ASE+ABACUS relaxation
- Sets up calculation parameters for each structure
- Prepares submission configuration for high-throughput computing

**Important:** The `-j` parameter accepts **ABACUS calculation input directories**, not structure files. Each directory must contain complete ABACUS inputs (INPUT, STRU, pseudopotentials, etc.).

**Input directory structure:**
```
project/
├── job1/
│   ├── INPUT      # ABACUS calculation parameters
│   ├── STRU       # Atomic structure
│   └── pp/        # Pseudopotentials (or symlinks)
├── job2/
│   ├── INPUT
│   ├── STRU
│   └── pp/
└── job3/
    ├── INPUT
    ├── STRU
    └── pp/
```

**Batch processing:** You can specify multiple job directories:
```bash
abacustest model aserelax prepare -j job1 job2 job3 job4 job5
# Or use wildcard (shell expansion)
abacustest model aserelax prepare -j job_*
```

### Prepare Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--jobs` | `-j` | **ABACUS input directories** (supports multiple for batch processing) | Current directory |
| `--rundftcommand` | `-c` | Command to execute aserelax | `abacustest model aserelax -o BFGS` |
| `--image` | `-i` | Docker image with ABACUS/ASE-ABACUS/abacustest | - |
| `--machine` | - | Machine type for remote execution | `c32_m128_cpu` |
| `--run` | `-r` | If run the test immediately (0=no, 1=yes) | `0` |

**Note:** Each job directory must contain complete ABACUS inputs: `INPUT`, `STRU`, and pseudopotentials (can be symlinks).

**Examples:**
```bash
# Single job directory
abacustest model aserelax prepare -j job1

# Multiple job directories (batch processing)
abacustest model aserelax prepare -j job1 job2 job3

# Using wildcard (shell expansion)
abacustest model aserelax prepare -j relax_job_*

# Custom command and machine
abacustest model aserelax prepare -j job1 job2 -c "abacustest model aserelax -o FIRE" --machine c64_m256_cpu
```

---

## Workflow 2: Post-process Only (Native ABACUS Relax)

Use this when you already have ABACUS relaxation results and want to analyze/plot them.

```bash
# Post-process existing ABACUS relax results from multiple job directories
abacustest model aserelax post -j job1 job2 job3
```

**What `post` does:**
- Reads relaxation trajectories from ABACUS output
- Generates energy vs. step plots
- Extracts final relaxed structures
- Creates summary reports

**Important:** The `-j` parameter accepts **ABACUS calculation job directories** (not structure files). Each directory should contain the completed calculation output.

**Supported input:**
- ASE+ABACUS relaxation results (from Workflow 1)
- Native ABACUS relax results (`calculation = relax` in INPUT)
- VASP relaxation results

### Post Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--jobs` | `-j` | **ABACUS job directories** with completed calculations (supports multiple) | Current directory |
| `--output` | `-o` | Output plot filename | `aserelax.png` |
| `--result` | `-r` | Save plot data to JSON file | `result.json` |
| `--type` | `-t` | Job type for data reading: `ase`, `abacus`, or `vasp` | `ase` |
| `--metric` | - | Post-process a specific metrics.json file | - |
| `--noplot` | - | Skip generating the plot | Plot enabled |

**Note:** Each job directory should contain the completed calculation output (e.g., `running_scf.log`, `OUTPUT/` directory, or `metrics.json`).

**Examples:**
```bash
# Post-process single job
abacustest model aserelax post -j job1

# Post-process multiple jobs (batch)
abacustest model aserelax post -j job1 job2 job3

# Post-process ABACUS native relax results
abacustest model aserelax post -j job1 job2 -t abacus -o energy_plot.png

# Post-process specific metrics file
abacustest model aserelax post --metric custom_metrics.json --noplot
```

---

## Workflow 3: Direct Run (Local ASE Relaxation)

Directly perform relaxation using ASE in the current environment.

```bash
# Direct relaxation (requires ase and ase-abacus installed)
abacustest model aserelax -j job1 -o BFGS --fmax 0.03

# Batch process multiple job directories
abacustest model aserelax -j job1 job2 job3 -o FIRE --fmax 0.01
```

**Requirements:**
- `ase` package installed
- `ase-abacus` calculator installed
- ABACUS executable available in PATH

**Install ase-abacus:**
```bash
# One-line installation from GitLab
pip install git+https://gitlab.com/1041176461/ase-abacus.git
```

This installs both `ase` (as a dependency) and `ase-abacus` in a single command.

**Important:** The `-j` parameter accepts **ABACUS calculation input directories**. Each directory must contain complete ABACUS inputs (INPUT, STRU, pseudopotentials).

**When to use:**
- Quick local relaxations
- Testing and development
- Small-scale calculations

**Quick Setup:**
```bash
pip install git+https://gitlab.com/1041176461/ase-abacus.git
```

### Direct Run Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--optimize` | `-o` | ASE optimizer method | `BFGS` |
| `--abacus` | `-a` | Path to ABACUS executable | `abacus` (from PATH) |
| `--fmax` | - | Force convergence threshold (eV/Å) | From INPUT or `0.0257112` |
| `--omp` | - | Number of OpenMP threads | `1` |
| `--mpi` | - | Number of MPI processes | All cores / OMP |
| `--cellrelax` | - | Relax cell (0=no, 1=yes) | Read from INPUT |
| `--job` | `-j` | **ABACUS input directories** (supports multiple for batch) | Current directory |

**Note:** Each job directory must contain complete ABACUS inputs: `INPUT`, `STRU`, and pseudopotentials.

**Supported Optimizers:**
`CG`, `BFGS`, `BFGSLineSearch`, `Berny`, `FIRE`, `GPMin`, `GoodOldQuasiNewton`, `LBFGS`, `LBFGSLineSearch`, `MDMin`, `ODE12r`, `QuasiNewton`, `RestartError`

**Examples:**
```bash
# Single job directory
abacustest model aserelax -j job1 -o FIRE --fmax 0.01 --omp 4

# Multiple job directories (batch processing)
abacustest model aserelax -j job1 job2 job3 -o LBFGS --fmax 0.02

# Cell relaxation
abacustest model aserelax -j bulk_job -o LBFGS --cellrelax 1 --mpi 8
```

---

## Output Files

After post-processing:

```
struct_dir/
├── ase_relax/
│   ├── relaxed_STRU       # Relaxed structure file
│   ├── relaxation_log.json # Relaxation trajectory
│   └── energy_vs_step.png  # Energy convergence plot
└── summary.json            # Summary of all jobs
```

**Plots generated:**
- `energy_vs_step.png` - Energy vs. relaxation step
- `force_vs_step.png` - Max force vs. step (if available)
- `structure_comparison.png` - Before/after structure overlay

---

## Examples

### Example 1: Batch Relaxation on Bohrium

```bash
# Prepare 10 job directories for remote calculation
# Each directory contains complete ABACUS inputs (INPUT, STRU, pp/)
abacustest model aserelax prepare -j job1 job2 job3

# Review generated scripts and setting.json
# Then submit (after user confirmation)
abacustest submit -p setting.json

# After calculations complete, post-process all jobs
abacustest model aserelax post -j job1 job2 job3 -o batch_relax.png
```

### Example 2: Post-process Native ABACUS Relax

```bash
# You already ran: abacustest model inputs --jtype relax
# And submitted with: abacustest submit -p setting.json

# Now analyze the results from multiple job directories
abacustest model aserelax post -j relax_job1 relax_job2 relax_job3 -t abacus
```

### Example 3: Local Direct Relaxation (Batch)

```bash
# Quick local relaxation for multiple job directories (requires ASE installed)
abacustest model aserelax -j job1 job2 job3 -o FIRE --fmax 0.01 --omp 4
```

### Example 4: Custom Machine and Command

```bash
# Prepare with custom machine and optimizer for batch jobs
abacustest model aserelax prepare -j job1 job2 job3 \
    -c "abacustest model aserelax -o LBFGS --fmax 0.02" \
    --machine c64_m256_cpu
```

### Example 5: Using Wildcard for Batch Processing

```bash
# If you have many job directories with common prefix
abacustest model aserelax prepare -j relax_job_*

# Post-process all completed jobs
abacustest model aserelax post -j relax_job_* -t abacus
```

---

## Common Issues

| Issue | Solution |
|-------|----------|
| `ase` not found | Install: `pip install ase` |
| `ase-abacus` not found | Install: `pip install git+https://gitlab.com/1041176461/ase-abacus.git` |
| ABACUS executable not found | Add to PATH or use `-a /path/to/abacus` |
| Relaxation not converging | Try `-o FIRE` or reduce `--fmax` |
| Cell relaxation fails | Ensure stress calculation is enabled in INPUT |
| metrics.json not found | Run prepare/submit first, or use `-t abacus` for native relax |

---

## Related

- Main router: [`abacustest-models`](../SKILL.md)
- Input generation: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
- Native ABACUS relax: Use `inputs --jtype relax`
- Submission: [`abacustest-submit`](../abacustest-submit/SKILL.md)
