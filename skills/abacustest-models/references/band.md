# Band - Band Structure

Calculate and plot electronic band structure with automatic high-symmetry k-path generation.

---

## Overview

The `band` model provides a complete workflow for band structure calculations:

| Mode | Command | Description |
|------|---------|-------------|
| **prepare** | `abacustest model band prepare` | Prepare SCF + NSCF + Band calculation inputs |
| **post** | `abacustest model band post` | Plot band structure and calculate band gap |
| **Direct Run** | ❌ | Not supported - use prepare + post workflow |

---

## Use When

- Calculate band gap (direct or indirect)
- Plot band structure diagrams (E-k dispersion)
- Analyze electronic band dispersion
- Determine semiconductor type (direct/indirect gap)
- Identify VBM (valence band maximum) and CBM (conduction band minimum)
- Study band alignment and electronic properties

---

## Workflow: Prepare + Submit + Post

### Step 1: Prepare Band Calculation

```bash
# Prepare SCF + NSCF + Band inputs for multiple structures
abacustest model band prepare -j job1 job2 job3
```

**What `prepare` does:**
- Generates three-stage calculation: SCF → NSCF → Band
- Automatically determines high-symmetry k-path from crystal structure
- Creates INPUT, STRU, and KPT files for each stage
- Prepares submission configuration for remote computing

### Step 2: Submit to Remote Cluster

```bash
# Review generated files, then submit (after user confirmation)
abacustest submit -p setting.json
```

**What happens during submit:**
- Creates `results/` directory (or path specified by `save_path` in setting.json)
- Submits SCF → NSCF → Band calculations to Bohrium cluster
- Upon completion, calculation outputs are stored in `results/<job_id>/`
- Generates `metrics.json` and `abacus.json` with extracted results

### Step 3: Post-process and Plot

```bash
# Post-process from results directory after calculations complete
abacustest model band post -j results/000000

# Custom energy range
abacustest model band post -j results/000000 --range -5 5
```

**Important:** After `abacustest submit`, the calculation results are in the `results/` directory. Use `-j results/<job_id>` for post-processing.

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

Each job directory must contain complete ABACUS inputs. The prepare step will create `01_scf/`, `02_nscf/`, `03_band/` subdirectories.

---

## Parameters

### Prepare Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--job` | `-j` | **ABACUS input directories** (supports multiple for batch) | Current directory |
| `--rundftcommand` | `-c` | Command to execute ABACUS | `OMP_NUM_THREADS=1 mpirun -np 16 abacus \| tee out.log` |
| `--image` | `-i` | Docker image for Bohrium | `registry.dp.tech/dptech/abacus-stable:LTSv3.10` |
| `--machine` | - | Machine type for Bohrium | `c32_m64_cpu` |
| `--run` | `-r` | Run immediately (0=no, 1=yes) | `0` |

**Example:**
```bash
# Single job
abacustest model band prepare -j job_si

# Multiple jobs (batch)
abacustest model band prepare -j job_si job_Ge job_GaN

# Custom machine
abacustest model band prepare -j job1 --machine c64_m128_cpu
```

### Post Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--job` | `-j` | **Job directories** with band calculation results (supports multiple) | Current directory |
| `--input` | - | NSCF input filename | `INPUT.nscf` |
| `--kpt` | - | K-points filename | `KPT.nscf` |
| `--range` | - | Energy range for plot (efermi ± eV) | `-5 5` (eV) |

**Example:**
```bash
# Post-process single job
abacustest model band post -j job1

# Post-process multiple jobs
abacustest model band post -j job1 job2 job3

# Custom energy range (efermi ± 5 eV)
abacustest model band post -j job1 --range -5 5

# Custom input files
abacustest model band post -j job1 --input INPUT.band --kpt KPT.band
```

---

## Output Files

**After submit (calculation complete):**
```
project/
├── job1/             # Original input directory
│   ├── INPUT.scf
│   ├── INPUT.nscf
│   └── ...
├── setting.json       # Submission configuration
└── results/           # Directory with actual calculation outputs
    └── job1/
        ├── INPUT.scf      # SCF input
        ├── INPUT.nscf     # NSCF input
        ├── STRU           # Structure file
        ├── KPT.scf        # SCF k-points
        ├── KPT.nscf       # NSCF k-points
        ├── scf.log        # SCF output log
        ├── nscf.log       # NSCF output log
        ├── OUT.ABACUS/    # ABACUS output directory
        └── ...
```

**After post-processing:**
```
results/job1/
└── band.png         # Band structure plot (generated in results directory)
```

**Console output from post-processing:**
```
E_fermi=6.76eV, BandGap=0.00eV
The band structure is plotted in the band.png file
```

**Note:** Currently `band.png` is generated in the `results/job1/` directory. The band gap and Fermi level are printed to console. Additional output files (band.dat, band_result.json) may require further configuration.

---

## Examples

### Example 1: Silicon Band Structure

```bash
# Prepare band calculation
abacustest model band prepare -j job_si

# Submit after user confirmation
abacustest submit -p setting.json

# Wait for calculations to complete (results in results/ directory)

# Post-process from results directory
abacustest model band post -j results/job_si

```

### Example 2: Batch Processing Multiple Materials

```bash
# Prepare band calculations for 3 materials
abacustest model band prepare -j job_si job_Ge job_GaN

# Submit all
abacustest submit -p setting.json

# Post-process all
cd results
abacustest model band post -j job*
```

### Example 3: Custom Energy Range

```bash
# Plot band structure with wider energy range (efermi ± 8 eV)
abacustest model band post -j job1 --range -8 8
```

### Example 4: Spin-Polarized Band Structure

```bash
# Prepare INPUT with nspin = 2 before running prepare
# Then run standard workflow
abacustest model band prepare -j job1
abacustest submit -p setting.json
abacustest model band post -j results/job1
```

---

## Calculation Stages

| Stage | Purpose | K-points |
|-------|---------|----------|
| **SCF** | Self-consistent charge density | Standard mesh |
| **NSCF** | Non-SCF with dense k-mesh | Dense uniform mesh |

**High-symmetry k-path** is automatically determined from crystal structure using SeeK-path conventions.

---

## Tips

| Recommendation | Reason |
|----------------|--------|
| Optimize structure first | Relaxed geometry gives accurate band structure |
| Use dense k-mesh for NSCF | Better band resolution and accurate gap |
| Check gap type | Direct vs indirect affects optical properties |
| Verify k-path | High-symmetry points should match literature |
| Spin-polarized systems | Set `nspin = 2` in INPUT for magnetic materials |

---

## Common Issues

| Issue | Solution |
|-------|----------|
| No band gap in result | Check if system is metallic; verify NSCF converged |
| Wrong k-path | Ensure STRU has correct lattice and atomic positions |
| Band plot looks strange | Check energy range; may need `--range -10 10` |
| SCF not converged | Increase `scf_nmax1` or adjust mixing parameters in INPUT |
| Missing 03_band directory | Prepare step may have failed; check error logs |

---

## Related

- Main router: [`abacustest-models`](../SKILL.md)
- DOS/PDOS: [`dos.md`](dos.md)
- Input generation: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
- Convergence tests: [`conv.md`](conv.md)
