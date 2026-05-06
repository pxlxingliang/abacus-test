# VASP2ABACUS - VASP to ABACUS Conversion

Transform VASP input files to ABACUS input files.

---

## Overview

The `vasp2abacus` model is a **direct-run utility**:

| Mode | Command | Description |
|------|---------|-------------|
| **prepare** | ❌ | Not supported |
| **post** | ❌ | Not supported |
| **Direct Run** | ✅ | `abacustest model vasp2abacus [args]` - Convert VASP → ABACUS |

**Note:** This is a **direct-run model** that converts VASP inputs to ABACUS format immediately.

**What gets converted:**
- ✅ POSCAR → STRU (structure)
- ✅ INCAR → INPUT (most parameters except ENCUT and EDIFF)
- ✅ KPOINTS → KPT (k-point mesh)

---

## Use When

- Convert existing VASP inputs to ABACUS format
- Migrate workflows from VASP to ABACUS
- Use VASP structures with ABACUS calculations
- Compare VASP and ABACUS results on same structure
- Leverage existing VASP inputs for ABACUS calculations

---

## Usage: Direct Run

```bash
# Convert single VASP job
abacustest model vasp2abacus -j vasp-job

# Convert multiple VASP jobs (batch)
abacustest model vasp2abacus -j vasp-job1 vasp-job2 vasp-job3

# With pseudopotential path
abacustest model vasp2abacus -j vasp-job --pp /path/to/pp_library

# With template INPUT file
abacustest model vasp2abacus -j vasp-job --input template.INPUT
```

---

## Input Directory Structure

**Input:** VASP job directory with standard VASP files

```
vasp-job/
├── INCAR      # VASP calculation parameters
├── POSCAR     # VASP structure file
├── KPOINTS    # VASP k-points 
└── POTCAR     # VASP pseudopotentials (not converted)
```

**Output:** ABACUS job directory

```
vasp-job/
├── INPUT      # Generated ABACUS input file
├── STRU       # Generated ABACUS structure file
├── KPT        # Generated ABACUS k-points (if KPOINTS existed)
└── pp/        # Pseudopotentials (must be provided separately)
```

---

## Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--job` | `-j` | **VASP job directories** (supports multiple for batch) | Current directory |
| `--pp` | - | Path to pseudopotential library, or `ABACUS_PP_PATH` env var | - |
| `--orb` | - | Path to orbital library, or `ABACUS_ORB_PATH` env var | - |
| `--input` | - | Template INPUT file for default ABACUS parameters | - |

**Notes:**
- `--pp` and `--orb` can use environment variables `ABACUS_PP_PATH` and `ABACUS_ORB_PATH`
- Template INPUT file parameters override converted INCAR parameters
- ENCUT and EDIFF from INCAR are NOT converted; must be set in template or manually


## Tips

| Recommendation | Reason |
|----------------|--------|
| Review converted INPUT | Verify all parameters converted correctly |
| Set ecutwfc manually | ENCUT not automatically converted |
| Set scf_thr manually | EDIFF not automatically converted |
| Check STRU format | Verify species blocks and atomic positions |
| Test with single job first | Verify conversion before batch processing |

**Parameters that need manual attention:**
- `ecutwfc` (from INCAR ENCUT) - set in INPUT or template
- `scf_thr` (from INCAR EDIFF) - set in INPUT or template
- Pseudopotentials - must provide `--pp` path or use `ABACUS_PP_PATH`
- Orbitals - must provide `--orb` path or use `ABACUS_ORB_PATH`

---

## Common Issues

| Issue | Solution |
|-------|----------|
| INPUT file missing parameters | Use `--input` template to set defaults |
| Pseudopotential not found | Provide `--pp` path or set `ABACUS_PP_PATH` |
| STRU format error | Check POSCAR format; ensure valid VASP structure |
| KPOINTS not converted | Ensure KPOINTS file exists in VASP directory |
| Wrong atomic species | Verify POTCAR matches POSCAR species order |

---

## Related

- Main router: [`abacustest-models`](../SKILL.md)
- Input generation: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
- Supercell: [`supercell.md`](supercell.md)
- Pseudopotentials: See `--pp` and `ABACUS_PP_PATH`
