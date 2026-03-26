---
name: abacustest-abacus2VaspQeCp2k
description: "Batch convert ABACUS input files to VASP/QE/CP2K formats. Use when: user wants to convert ABACUS calculations to other software, including both structure (STRU→POSCAR) and input parameters (INPUT→INCAR/pw.inp)."
metadata: { "openclaw": { "emoji": "🔄", "requires": { "pip": ["abacustest"] } } }
---

# abacustest Convert

Batch convert ABACUS input files (STRU + INPUT + KPT) to VASP, QE, or CP2K formats.

## When to Use

✅ **Use this skill**: Batch convert ABACUS inputs to VASP/QE/CP2K for cross-software comparison

❌ **Do not use**: Extract results (`abacustest-extract-dft-results`)

## Command

```bash
abacustest prepare -p param.json -s output_dir
```

**Note**: Conversion uses `abacustest prepare` command with specific parameters in `param.json`.

---

## ⚠️ Limitations: Parameter Conversion

**Different software have different input parameters - not all parameters can be automatically converted.**

- **What can be converted**: Common parameters (calculation type, cutoff, k-points, spin, etc.)
- **What cannot be converted**: Software-specific parameters (mixing schemes, advanced algorithms, etc.)
- **Unconverted parameters**: Will be printed to screen during execution for your reference

**Recommendation**: After conversion, manually review the generated files and add/adjust parameters as needed.


## Conversion: ABACUS → VASP

Convert ABACUS input files to VASP format (INCAR, POSCAR, KPOINTS, POTCAR).

### ⚠️ Before You Start: Pseudopotential Library

**Question**: Do you have a VASP pseudopotential library (e.g., PAW_PBE, POTCAR files)?

- **Yes**: Provide the path (e.g., `/path/to/vasp/PAW_PBE`)
- **No**: POTCAR will NOT be automatically configured; you need to prepare it manually

**Note**: ABACUS uses norm-conserving pseudopotentials, while VASP uses PAW pseudopotentials. They are not directly compatible.

### Configuration

**Option 1: Using pseudopotential library (recommended)**

```json
{
  "prepare": {
    "example_template": ["abacus_calc"],
    "abacus2vasp": true,
    "potcar": "/path/to/vasp/PAW_PBE"
  }
}
```

**Option 2: Manual POTCAR mapping**

```json
{
  "prepare": {
    "example_template": ["abacus_calc"],
    "abacus2vasp": true,
    "potcar": {"Fe": "/path/Fe.pot", "C": "/path/C.pot"}
  }
}
```

**Option 3: No POTCAR (manual preparation)**

```json
{
  "prepare": {
    "example_template": ["abacus_calc"],
    "abacus2vasp": true
  }
}
```

### Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `abacus2vasp` | Enable conversion | `true` |
| `example_template` | ABACUS input directories | `["calc1", "calc2"]` |
| `potcar` | POTCAR library path OR file mapping | `"/path/PAW_PBE"` or `{"Fe": "Fe.pot"}` |
| `vasp_setting` | Additional VASP INCAR parameters (dict, optional) | `{"ENCUT": 500}` |

**`vasp_setting` Details**:

- **Type**: Dictionary (dict)
- **Purpose**: Additional/custom VASP INCAR parameter settings
- **Location**: Under `prepare` in param.json
- **Default behavior**: If not specified, abacustest automatically converts ABACUS INPUT parameters to VASP equivalents
- **Override behavior**: If a parameter is set in `vasp_setting`, it will override the automatically converted value

**Example** (full param.json structure):
```json
{
  "prepare": {
    "example_template": ["abacus_calc"],
    "abacus2vasp": true,
    "potcar": "/path/to/PAW_PBE",
    "vasp_setting": {
      "ENCUT": 500,      // Override automatic ENCUT conversion
      "ISMEAR": 0,       // Override automatic smearing setting
      "EDIFF": 1e-5      // Override automatic convergence threshold
    }
  }
}
```

**Special Parameter: `emax_coef`**

`vasp_setting` supports a special parameter `emax_coef` for automatic ENCUT calculation:

```python
# How emax_coef works:
encut_recommended = max([POTCAR[i].ENCUT for i in all_elements])
ENCUT = emax_coef * encut_recommended

# Example: Fe (ENCUT=240.6 eV) + C (ENCUT=400.0 eV)
# encut_recommended = max(240.6, 400.0) = 400.0 eV
# If emax_coef = 1.3:
# ENCUT = 1.3 * 400.0 = 520.0 eV
```

**ENCUT Priority** (highest to lowest):
1. `vasp_setting.ENCUT` → Use explicit value directly
2. `vasp_setting.emax_coef` → Calculate: `emax_coef × max(POTCAR ENCUT)`
3. Not specified → Auto-convert: `ecutwfc (ABACUS) × Ry2eV`

**Important Notes**:

1. **`potcar` as library path**: abacustest will automatically find and configure POTCAR files from the VASP pseudopotential library based on elements in the structure.

2. **`potcar` not provided**: POTCAR will NOT be automatically configured. You need to prepare POTCAR files manually after conversion.

3. **For ABACUS/VASP comparison**: 
   - Recommend NOT setting `vasp_setting` (use automatic conversion)
   - Or only set `ENCUT` to account for pseudopotential differences
   - See "Why Only Set ENCUT for Comparison?" section below

### Why Only Set ENCUT for Comparison?

**ABACUS uses norm-conserving pseudopotentials, VASP uses PAW pseudopotentials.**

- PAW pseudopotentials are softer → VASP's kinetic energy cutoff is naturally **lower** than ABACUS
- Typical conversion: `ENCUT_VASP ≈ ENCUT_ABACUS × 0.6~0.8`
- Example: ABACUS ecutwfc=80 Ry (≈1088 eV) → VASP ENCUT=500-700 eV may be sufficient

### Output

```
output_dir/
└── abacus_calc/
    ├── INCAR      # VASP parameters (auto-converted or from vasp_setting)
    ├── POSCAR     # Structure
    ├── KPOINTS    # K-point mesh (if KPT file exists)
    └── POTCAR/    # Pseudopotentials (if potcar provided)
```

**Note on KPOINTS**:

- **If ABACUS uses `kspacing`**: abacustest converts it to `KSPACING` in INCAR (no KPOINTS file generated)
  - **Unit conversion**: ABACUS `kspacing` is in **1/Bohr**, VASP `KSPACING` is in **1/Å**
  - **Conversion**: `KSPACING (1/Å) = kspacing (1/Bohr) × 1.8897`
- **If ABACUS uses KPT file**: abacustest generates a KPOINTS file with the same k-point mesh

### Parameter Mapping (Automatic Conversion)

| ABACUS | VASP | Note |
|--------|------|------|
| `calculation` | `IBRION` | scf→-1, relax→2, cell-relax→3 |
| `ecutwfc` | `ENCUT` | Ry→eV (auto-adjusted for PAW) |
| `scf_thr` | `EDIFF` | Auto: scf_thr/1e-2 (PW) or /1e-1 (LCAO) |
| `kspacing` | `KPOINTS` | K-point density |
| `nspin` | `ISPIN` | Spin polarization |
| `mag` | `MAGMOM` | Magnetic moments |
| `dft_plus_u` | `LDAU` | DFT+U settings |

**Note**: If `vasp_setting` is specified, those parameters will use your values instead of automatic conversion.

---

## Conversion: ABACUS → QE

Convert ABACUS input files to Quantum ESPRESSO format (pw.in).

### Pseudopotentials

**abacustest uses pseudopotentials defined in the ABACUS STRU.**

- Pseudopotential files (.upf) from the ABACUS directory are copied to `output_dir/abacus_calc`
- The pw.in file references these pseudopotential files automatically
- **Note**: Ensure your ABACUS pseudopotentials are in UPF format (QE standard format)

---

### Configuration

```json
{
  "prepare": {
    "example_template": ["abacus_calc"],
    "abacus2qe": true,
    "qe_setting": {...}  
  }
}
```

### Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `abacus2qe` | Enable conversion | `true` |
| `example_template` | ABACUS input directories | `["calc1"]` |
| `qe_setting` | Additional QE input parameters (dict, optional) | See below |

**`qe_setting` Details**:

- **Type**: Dictionary (dict) with nested structure matching QE input format
- **Purpose**: Additional/custom QE parameters
- **Location**: Under `prepare` in param.json
- **Default behavior**: If not specified, abacustest automatically converts ABACUS INPUT parameters to QE equivalents
- **Override behavior**: If a parameter is set in `qe_setting`, it will override the automatically converted value

**Structure** (matching QE namelists):
```json
{
  "qe_setting": {
    "version": "7.0",           // Special parameter (see below)
    "control": {                // CONTROL namelist
      "calculation": "scf",
      "restart_mode": "from_scratch"
    },
    "system": {                 // SYSTEM namelist
      "ibrav": 0,
      "ecutwfc": 50,
      "ecutrho": 400
    },
    "electrons": {              // ELECTRONS namelist
      "conv_thr": 1e-8
    }
  }
}
```

**Special Parameter: `version`**

- **Type**: String
- **Purpose**: QE version number (e.g., `"6.5"`, `"7.0"`, `"7.1"`)
- **Why it matters**: Different QE versions use different formats for certain parameters, especially DFT+U (Hubbard) settings

**Example: DFT+U format changes by version**:
```
# QE < 7.0 (old format):
Hubbard_U(1) = 5.3
Hubbard_alpha(1) = 0.0

# QE >= 7.0 (new format):
Hubbard_U(Fe-3d) = 5.3
Hubbard_alpha(Fe-3d) = 0.0
```

abacustest uses the `version` parameter to generate the correct format for your QE version.

**Supported namelists in `qe_setting`**:
- `control` - CONTROL namelist
- `system` - SYSTEM namelist
- `electrons` - ELECTRONS namelist
- `ions` - IONS namelist
- `cell` - CELL namelist
- Block definitions: `"HUBBARD (ortho-atom)": ["U Fe1-3d 5.3"]`

QE input parameter reference: https://www.quantum-espresso.org/Doc/INPUT_PW.html

### Output

```
output_dir/
└── abacus_calc/
    ├── input      # QE input (pwscf format)
    └── *.upf      # Pseudopotentials (copied from ABACUS directory)
```

## Conversion: ABACUS → CP2K

Convert ABACUS input files to CP2K format (input.inp).

### Configuration

```json
{
  "prepare": {
    "example_template": ["abacus_calc"],
    "abacus2cp2k": true,
    "cp2k_setting": {
      "GLOBAL": {
        "RUN_TYPE": "ENERGY_FORCE"
      },
      "FORCE_EVAL": {
        "DFT": {
          "SCF": {
            "EPS_SCF": 1e-6
          },
          "QS": {
            "EPS_DEFAULT": 1e-10
          }
        }
      }
    }
  }
}
```

### Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `abacus2cp2k` | Enable conversion | `true` |
| `example_template` | ABACUS input directories | `["calc1"]` |
| `cp2k_setting` | CP2K input parameters | Nested dictionary |

**Note**: `cp2k_setting` structure matches CP2K input format. Keys are section names, values are subsections or parameters.

### Output

```
output_dir/
└── abacus_calc/
    └── input.inp  # CP2K input
```

### Parameter Mapping

| ABACUS | CP2K | Note |
|--------|------|------|
| `calculation` | `RUN_TYPE` | scf→ENERGY_FORCE |
| `ecutwfc` (PW) | `CUTOFF` | ×1000 (Ry→Ry) |
| `ecutwfc` (LCAO) | `REL_CUTOFF` | ×100 |
| `scf_thr` | `EPS_SCF` | SCF convergence |
| `force_thr` | `RMS_FORCE` | Force threshold |

---

## Common Errors

| Error | Fix |
|-------|-----|
| POTCAR not found | Check `potcar` paths or use full path |
| Wrong ENCUT after conversion | Set `emax_coef: 1.3` or define explicitly |
| QE version mismatch | Specify `qe_setting.version` |
| CP2K input format error | Verify `cp2k_setting` structure |

## Example: Cross-Software Benchmark

**Goal**: Compare ABACUS and VASP results for same structure

```bash
# 1. Run ABACUS calculation
# Directory: abacus-calc/ contains INPUT, STRU, KPT

# 2. Convert to VASP
cat > param.json << 'EOF'
{
  "prepare": {
    "example_template": ["abacus-calc"],
    "abacus2vasp": true,
    "potcar": {"Fe": "/path/Fe.pot"},
    "vasp_setting": {"ENCUT": 500, "EDIFF": 1e-5}
  }
}
EOF
abacustest prepare -p param.json -s vasp-calc

# 3. Run VASP calculation in vasp-calc/abacus-calc/

# 4. Compare results
abacustest collectdata -j abacus-calc -o abacus-result.json
abacustest collectdata -j vasp-calc/abacus-calc -o vasp-result.json
```

## Best Practices

1. ✅ Verify POTCAR/pseudopotential compatibility
2. ✅ Check converted parameters (especially ENCUT)
3. ✅ Use same k-point density across software
4. ✅ Ensure DFT+U settings are equivalent
5. ✅ Validate structure integrity after conversion

