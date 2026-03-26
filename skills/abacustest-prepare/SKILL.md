---
name: abacustest-prepare
description: "Prepare ABACUS input files from templates with parameter mixing and structure perturbations. Use when: user wants to generate multiple INPUT/KPT combinations, create perturbed structures for MD, or perform high-throughput input preparation from existing ABACUS templates."
metadata: { "openclaw": { "emoji": "🔧", "requires": { "pip": ["abacustest"] } } }
---

# abacustest Prepare

Prepare ABACUS input files (STRU, INPUT, KPT) from templates with parameter mixing and structure perturbations.

## ⚠️ Important: Skill Overlaps

This skill is part of a family of input preparation tools. Choose the right one for your task:

| Task | Recommended Skill |
|------|-------------------|
| **Convert CIF/POSCAR → ABACUS STRU** | [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md) ⭐ **Preferred** |
| **Mix parameters (ecutwfc, kspacing, etc.)** | `abacustest-prepare` (this skill) |
| **Generate perturbed structures** | `abacustest-prepare` (this skill) |
| **Convert ABACUS → VASP/QE/CP2K** | [`abacustest-abacus2VaspQeCp2k`](../abacustest-abacus2VaspQeCp2k/SKILL.md) ⭐ **See that skill** |

### Why Use `abacustest-prepare-inputs` for Structure Conversion?

While this skill (`abacustest-prepare`) can convert CIF/POSCAR to STRU, the [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md) skill is **recommended** because:

- ✅ Simpler command-line interface (no param.json needed for basic use)
- ✅ Built-in recommended parameters for different calculation types (SCF, relax, MD, band)
- ✅ Automatic DFT+U, magnetic, and SOC configuration
- ✅ Direct pseudopotential/orbital library download support
- ✅ Better for one-off structure conversions

**Use `abacustest-prepare` for structure conversion when:**
- You need to combine structure conversion with parameter mixing
- You need to combine structure conversion with structure perturbations
- You're already using param.json for other prepare features

### Why Use `abacustest-abacus2VaspQeCp2k` for Format Conversion?

For converting ABACUS inputs to VASP/QE/CP2K formats, see the dedicated [`abacustest-abacus2VaspQeCp2k`](../abacustest-abacus2VaspQeCp2k/SKILL.md) skill, which provides:
- ✅ Detailed conversion parameter documentation
- ✅ Pseudopotential library setup guidance
- ✅ Software-specific conversion limitations
- ✅ Complete examples for each target format

---

## When to Use

✅ **Use this skill**:
- "Generate inputs with different cutoff energies and k-points"
- "Create 20 perturbed structures for MD initialization"
- "Batch prepare high-throughput calculation inputs from templates"
- "Combine parameter mixing with structure perturbations"

❌ **Do not use**:
- Simple CIF/POSCAR → STRU conversion → Use [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
- ABACUS → VASP/QE/CP2K conversion → Use [`abacustest-abacus2VaspQeCp2k`](../abacustest-abacus2VaspQeCp2k/SKILL.md)
- Extract calculation results → Use `abacustest-extract-dft-results`

---

## Basic Usage

### Method 1: Using param.json (Recommended)

Create a `param.json` configuration file, then run:
```bash
abacustest prepare -p param.json -s output_dir
```

### Method 2: Command-line Arguments

```bash
abacustest prepare --stru structure.vasp --stru_format poscar -s output_dir
```

---

## Feature 1: Structure Conversion (POSCAR/CIF → STRU)

⚠️ **Note**: For simple structure file conversion, [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md) is recommended. Use this skill when combining conversion with other features.

### Example: Convert VASP POSCAR

**param.json**:
```json
{
  "prepare": {
    "strus": ["POSCAR"],
    "stru_format": "poscar",
    "pp_dict": {"Fe": "Fe.upf", "C": "C.upf"},
    "orb_dict": {"Fe": "Fe.orb", "C": "C.orb"},
    "pp_path": "/path/to/pseudopotentials",
    "orb_path": "/path/to/orbitals"
  }
}
```

**Execute**:
```bash
abacustest prepare -p param.json -s fe-c-structs
```

**Output**:
```
fe-c-structs/
└── 000000/
    ├── STRU
    ├── INPUT (template)
    └── KPT (template)
```

### Example: Convert CIF Files

**param.json**:
```json
{
  "prepare": {
    "strus": ["structure.cif"],
    "stru_format": "cif"
  }
}
```

### Parameter Reference

| Parameter | Type | Description |
|-----------|------|-------------|
| `strus` | list | Structure file paths (supports dpdata-readable formats) |
| `stru_format` | str | Format: `poscar`, `cif`, `deepmd/npy`, etc. |
| `input_template` | str | INPUT file template path |
| `kpt_template` | str | KPT file template path |
| `pp_dict` | dict | Pseudopotential mapping: `{"Fe": "Fe.upf"}` |
| `orb_dict` | dict | Orbital mapping: `{"Fe": "Fe.orb"}` |
| `pp_path` | str | Pseudopotential library path |
| `orb_path` | str | Orbital library path |

**Note on `pp_path`/`orb_path`**: These directories should contain an `element.json` file defining element names and filenames (e.g., `{"Fe": "Fe.upf"}`). If not present, abacustest uses the first two letters of filenames as element names.

---

## Feature 2: Parameter Mixing (High-Throughput Input Generation)

Generate multiple INPUT/KPT combinations from templates for convergence testing or high-throughput screening.

### Example: Mix Cutoff Energy and K-points

**param.json**:
```json
{
  "prepare": {
    "example_template": ["base_struct"],
    "mix_input": {
      "ecutwfc": [50, 60, 70],
      "kspacing": [0.1, 0.12]
    },
    "mix_kpt": [2, [3, 3, 3], [4, 4, 4]]
  }
}
```

**Execute**:
```bash
abacustest prepare -p param.json -s htf-inputs
```

**Output** (generates 3×2×3 = 18 input sets):
```
htf-inputs/
└── base_struct/
    ├── 000000/  (ecutwfc=50, kspacing=0.1, kpt=2)
    ├── 000001/  (ecutwfc=50, kspacing=0.1, kpt=3x3x3)
    ├── 000002/  ...
    └── ...
```

### mix_input Format

**Full combinatorial** (all combinations):
```json
{
  "mix_input": {
    "ecutwfc": [50, 60, 70],
    "kspacing": [0.1, 0.12]
  }
}
```
Generates 3×2 = 6 combinations.

**Paired values** (one-to-one correspondence using `|`):
```json
{
  "mix_input": {
    "ecutwfc|kspacing": ["50|0.1", "60|0.12", "70|0.15"]
  }
}
```
Generates 3 combinations (50/0.1, 60/0.12, 70/0.15).

### mix_kpt Format

Three ways to define K-points (Gamma-centered Monkhorst-Pack):

| Format | Example | Meaning |
|--------|---------|---------|
| Single integer | `2` | 2 2 2 0 0 0 |
| Three integers | `[3, 3, 3]` | 3 3 3 0 0 0 |
| Six values | `[4, 4, 4, 1, 1, 1]` | 4 4 4 1 1 1 (with shift) |

**Example**:
```json
{
  "mix_kpt": [2, [3, 3, 3], [4, 4, 4, 1, 1, 1]]
}
```
Generates 3 K-point configurations.

---

## Feature 3: Structure Perturbations

Generate perturbed structures for MD initialization, phonon calculations, or uncertainty quantification.

### Example: Generate 10 Perturbed Structures

**param.json**:
```json
{
  "prepare": {
    "example_template": ["initial_struct"],
    "pert_stru": {
      "pert_number": 10,
      "cell_pert_frac": 0.01,
      "atom_pert_dist": 0.1,
      "mag_rotate_angle": 5,
      "mag_tilt_angle": 10,
      "mag_norm_dist": 0.05
    }
  }
}
```

**Execute**:
```bash
abacustest prepare -p param.json -s perturbed-structs
```

**Output**:
```
perturbed-structs/
└── initial_struct/
    ├── 000000/  (perturbed structure 1)
    ├── 000001/  (perturbed structure 2)
    └── ...
```

### pert_stru Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `pert_number` | int | Number of perturbed structures | - |
| `cell_pert_frac` | float | Max cell perturbation fraction (0.01 = 1%) | `null` (no cell perturbation) |
| `atom_pert_dist` | float | Max atom displacement (Å) | `null` (no atom perturbation) |
| `mag_rotate_angle` | float | Magnetic moment rotation angle (degrees) | `null` (no rotation) |
| `mag_tilt_angle` | float | Magnetic moment tilt angle (degrees) | `null` (no tilt) |
| `mag_norm_dist` | float | Magnetic moment norm displacement (μ_B) | `null` (no displacement) |

**Notes**:
- Values can be lists `[min, max]` for random ranges: `"atom_pert_dist": [0.1, 0.15]`
- Magnetic perturbations require spin-constrained atoms (`sc = 1` in STRU)
- `pert_number` is the final count; each structure is independently perturbed

---

## Feature 4: Format Conversion (ABACUS → VASP/QE/CP2K)

⚠️ **Moved**: For ABACUS to VASP/QE/CP2K conversion, see the dedicated **[`abacustest-abacus2VaspQeCp2k`](../abacustest-abacus2VaspQeCp2k/SKILL.md)** skill.

That skill provides comprehensive documentation on:
- Pseudopotential library setup for target software
- Parameter conversion limitations
- Software-specific settings (vasp_setting, qe_setting, cp2k_setting)
- Complete examples for each conversion type

**Quick reference**:
```json
{
  "prepare": {
    "example_template": ["abacus_input"],
    "abacus2vasp": true,
    "potcar": {"Fe": "Fe.pot"},
    "vasp_setting": {"ENCUT": 500}
  }
}
```

For full documentation, see [`abacustest-abacus2VaspQeCp2k`](../abacustest-abacus2VaspQeCp2k/SKILL.md).

---

## Complete Parameter Reference

### Structure Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `strus` | list | Structure file paths |
| `stru_format` | str | Structure format (poscar/cif/deepmd/npy) |
| `example_template` | list | ABACUS input template directories |
| `input_template` | str | INPUT template file |
| `kpt_template` | str | KPT template file |

### Pseudopotential/Orbital Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `pp_dict` | dict | Pseudopotential dict `{"Fe": "Fe.upf"}` |
| `orb_dict` | dict | Orbital dict `{"Fe": "Fe.orb"}` |
| `pp_path` | str | Pseudopotential library path |
| `orb_path` | str | Orbital library path |

### Parameter Mixing

| Parameter | Type | Description |
|-----------|------|-------------|
| `mix_input` | dict | INPUT parameter combinations |
| `mix_kpt` | list | K-point combinations |
| `mix_stru` | list | STRU file combinations |

### Structure Perturbation

| Parameter | Type | Description |
|-----------|------|-------------|
| `pert_stru` | dict | Perturbation parameters |
| `pert_number` | int | Number of perturbed structures |
| `cell_pert_frac` | float | Cell perturbation fraction |
| `atom_pert_dist` | float | Atom perturbation distance (Å) |

### Format Conversion

| Parameter | Type | Description |
|-----------|------|-------------|
| `abacus2vasp` | bool | Convert to VASP |
| `abacus2qe` | bool | Convert to QE |
| `abacus2cp2k` | bool | Convert to CP2K |

See [`abacustest-abacus2VaspQeCp2k`](../abacustest-abacus2VaspQeCp2k/SKILL.md) for detailed conversion settings.

### Other Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `link_example_template_extra_files` | bool | Link extra files from templates | `true` |
| `extra_files` | list | Additional files to link | `[]` |
| `dpks_descriptor` | str | DeepKS descriptor file | `null` |

---

## Execution Order

When multiple features are combined, abacustest prepare runs in this order:

1. **Structure conversion** (`strus` → ABACUS inputs)
2. **Structure perturbation** (generate perturbed sub-examples)
3. **Parameter mixing** (INPUT/KPT combinations for each sub-example)
4. **Format conversion** (ABACUS → VASP/QE/CP2K)

**Rules**:
- When both `example_template` and `strus` are specified, `example_template` is ignored
- When `example_template` contains INPUT/KPT and `input_template`/`kpt_template` are also specified, the template values override

---

## Example Workflows

### Workflow 1: High-Throughput Cutoff Testing

**Goal**: Test ecutwfc = 40, 50, 60, 70, 80 Ry for one structure

**param.json**:
```json
{
  "prepare": {
    "example_template": ["base_struct"],
    "mix_input": {
      "ecutwfc": [40, 50, 60, 70, 80]
    }
  }
}
```

### Workflow 2: MD Initial Structures

**Goal**: Generate 20 perturbed structures for AIMD

**param.json**:
```json
{
  "prepare": {
    "example_template": ["equilibrium_struct"],
    "pert_stru": {
      "pert_number": 20,
      "atom_pert_dist": 0.05
    }
  }
}
```

### Workflow 3: Combined Conversion + Mixing

**Goal**: Convert CIF files and test multiple parameters

**param.json**:
```json
{
  "prepare": {
    "strus": ["structure1.cif", "structure2.cif"],
    "stru_format": "cif",
    "pp_path": "/path/to/pp",
    "mix_input": {
      "ecutwfc": [50, 60],
      "kspacing": [0.1, 0.12]
    }
  }
}
```

Generates 2 structures × 2 ecutwfc × 2 kspacing = 8 input sets.

### Workflow 4: Multi-Software Comparison

**Goal**: Prepare ABACUS, VASP, and QE inputs from one structure

**param.json**:
```json
{
  "prepare": {
    "example_template": ["reference_struct"],
    "abacus2vasp": true,
    "abacus2qe": true,
    "potcar": {"Fe": "Fe.pot"},
    "vasp_setting": {"ENCUT": 500},
    "qe_setting": {"system": {"ecutwfc": 50}}
  }
}
```

See [`abacustest-abacus2VaspQeCp2k`](../abacustest-abacus2VaspQeCp2k/SKILL.md) for full conversion documentation.

---

## Common Errors

### Error 1: Structure File Path Incorrect

```json
{
  "prepare": {
    "strus": ["wrong_path/POSCAR"]
  }
}
```

**Fix**:
```bash
# Verify file exists first
ls -la correct_path/POSCAR
```

```json
{
  "prepare": {
    "strus": ["correct_path/POSCAR"]
  }
}
```

### Error 2: Pseudopotential Files Not Found

```json
{
  "prepare": {
    "pp_dict": {"Fe": "Fe.upf"}
  }
}
# Error: Fe.upf does not exist
```

**Fix**:
```bash
# Use pp_path to specify pseudopotential library
ls /path/to/pseudopotentials/Fe_*.upf
```

```json
{
  "prepare": {
    "pp_path": "/path/to/pseudopotentials"
  }
}
```

### Error 3: Mix Parameter Syntax

```json
{
  "prepare": {
    "mix_input": {
      "ecutwfc": [50, 60, 70],
      "kspacing": [0.1, 0.12, 0.15]
    }
  }
}
# Generates 3×3 = 9 combinations (correct for full combinatorial)
```

```json
{
  "prepare": {
    "mix_input": {
      "ecutwfc|kspacing": ["50|0.1", "60|0.12", "70|0.15"]
    }
  }
}
# Generates 3 combinations (paired one-to-one)
```

---

## Tips

| Recommendation | Reason |
|----------------|--------|
| Use `abacustest-prepare-inputs` for simple CIF→STRU | Simpler CLI, recommended parameters |
| Test with one example first | Parameter mixing can generate many files |
| Use relative paths | Avoids issues when moving directories |
| Verify pseudopotential library | Missing PP files cause ABACUS failures |
| Check perturbation magnitudes | Too large perturbations break structures |

---

## Related Skills

- **Structure → ABACUS**: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md) ⭐ Preferred for simple conversions
- **ABACUS → VASP/QE/CP2K**: [`abacustest-abacus2VaspQeCp2k`](../abacustest-abacus2VaspQeCp2k/SKILL.md)
- **Submit jobs**: [`abacustest-submit`](../abacustest-submit/SKILL.md)
- **Extract results**: [`abacustest-extract-dft-results`](../abacustest-extract-dft-results/SKILL.md)
- **Specialized models**: [`abacustest-models`](../abacustest-models/SKILL.md)

---

## External Resources

- **abacustest GitHub**: https://github.com/pxlxingliang/abacus-test
- **ABACUS Docs**: https://abacus.oss.cn-north-1.aliyuncs.com/
- **dpdata** (structure formats): https://github.com/deepmodeling/dpdata
