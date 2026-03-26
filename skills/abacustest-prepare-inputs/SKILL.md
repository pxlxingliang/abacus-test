---
name: abacustest-prepare-inputs
description: "Batch configure structure files (CIF/POSCAR/STRU, etc.) as ABACUS input files for various calculation types (SCF, relax, cell-relax, MD, band), and download ABACUS recommended pseudopotential and orbital libraries (APNS-v1). Use when: user wants to batch prepare ABACUS inputs from structure files, download recommended pseudopotentials and orbital libraries, set up DFT+U/magnetic calculations with recommended parameters."
metadata: { "openclaw": { "emoji": "📝", "requires": { "pip": ["abacustest"] } } }
---

# abacustest Model: Inputs

This skill has two main functions:
1. **Download recommended pseudopotentials and orbitals** for ABACUS
2. **Generate ABACUS input files** from structure files with recommended parameters

---

## Part 1: Download Recommended Pseudopotentials & Orbitals

ABACUS requires pseudopotential files (`.upf`) for all calculations, and orbital files (`.orb`) for LCAO basis calculations. This skill can download the official APNS (Atomic Pseudopotential & Numerical Basis Set) library.

### ⚠️ IMPORTANT: Always Confirm Before Downloading

**DO NOT download automatically.** The pseudopotential/orbital library is large (~several MB). Always confirm with the user first:

```
⚠️ Pseudopotential/Orbital Library Download

The APNS v1 library contains:
- Pseudopotentials (.upf files) for most elements (H to Rn)
- Orbital files (.orb files) for LCAO calculations
- Total size: ~50 MB
- Download time: Varies by network connection

Before proceeding, please confirm:

□ Do you already have pseudopotentials/orbitals installed?
  → If yes, skip download and set ABACUS_PP_PATH/ABACUS_ORB_PATH environment variables

□ Do you want to download the recommended APNS v1 library?
  → If yes, I will run: abacustest model inputs --download-pporb apns-v1
  → After download, you'll need to move files to a shared location and set environment variables

□ Or do you have a preferred pseudopotential library?
  → If yes, please provide the path

Please confirm before I proceed with the download.
```

### When to Offer Download

Only suggest downloading when:
1. User explicitly asks to download pseudopotentials/orbitals
2. User has no existing PP/orbital setup (no `--pp/--orb` args, no env vars)
3. User confirms they want to proceed

### Download Command (After Confirmation)

```bash
abacustest model inputs --download-pporb apns-v1
```

### What Gets Downloaded

- **Pseudopotentials** (`pp/`): `.upf` files for most elements 
- **Orbitals** (`orb/`): `.orb` files for LCAO calculations
- **Version**: APNS v1 (ABACUS recommended)

### Setup After Download

After downloading, move the files to a shared location and set environment variables:

```bash
# Move to a shared location (example)
mv ABACUS-APNS-PPORBs-v1 /shared/ABACUS-APNS-v1/

# Add to your shell config (~/.bashrc or ~/.zshrc)
export ABACUS_PP_PATH=/shared/ABACUS-APNS-v1/pp
export ABACUS_ORB_PATH=/shared/ABACUS-APNS-v1/orb

# Apply changes
source ~/.bashrc  # or source ~/.zshrc
```

### Verify Setup

```bash
# Check environment variables
echo $ABACUS_PP_PATH
echo $ABACUS_ORB_PATH

# List available pseudopotentials
ls $ABACUS_PP_PATH/*.upf | head -5
```

Once set up, you won't need to specify `--pp` or `--orb` paths for future calculations.

---

## Part 2: Generate ABACUS Input Files

Convert structure files (CIF, POSCAR, STRU) to complete ABACUS input directories with recommended parameters.

### ⚠️ Before You Start: Check Pseudopotential & Orbital Paths

ABACUS requires pseudopotential files for all calculations. Before generating inputs, follow this checklist:

```
┌─────────────────────────────────────────────────────────────┐
│  Step 1: Did you provide --pp/--orb command-line arguments? │
│  └─ Yes → Use provided paths, skip to Step 4               │
│  └─ No  → Continue to Step 2                               │
├─────────────────────────────────────────────────────────────┤
│  Step 2: Check if environment variables are set:            │
│  $ echo $ABACUS_PP_PATH                                     │
│  $ echo $ABACUS_ORB_PATH                                    │
│  └─ Set → Proceed using environment variable paths          │
│  └─ Not set → Continue to Step 3                            │
├─────────────────────────────────────────────────────────────┤
│  Step 3: No paths found! Options:                           │
│  a) Download recommended PP/orbital:                        │
│     abacustest model inputs --download-pporb apns-v1        │
│  b) Ask user: "Where are your pseudopotentials located?"    │
│  c) Skip PP/orbital linking (manual editing required)       │
├─────────────────────────────────────────────────────────────┤
│  Step 4: Generate inputs with confirmed paths               │
└─────────────────────────────────────────────────────────────┘
```

### Basic Usage

```bash
# From CIF file
abacustest model inputs -f structure.cif

# From VASP POSCAR
abacustest model inputs -f POSCAR --ftype poscar

# From ABACUS STRU
abacustest model inputs -f stru.cif --ftype abacus

# Batch process multiple files
abacustest model inputs -f *.cif
```

### Generated Output

For each structure, a directory is created containing:
```
000000/
├── INPUT           # Calculation parameters
├── STRU            # Structure file
├── KPT             # K-point settings
├── *.upf           # Pseudopotential links
├── *.orb           # Orbital links (if LCAO)
└── run.sh          # Run script
```

---

## Calculation Types & Recommended Parameters

Use `--jtype` to specify the calculation type. The skill automatically sets appropriate parameters for each type.

| `--jtype` | Description | Key Parameters |
|-----------|-------------|----------------|
| `scf` | Self-consistent field (default) | `ecutwfc 80`, `kspacing 0.14`, `scf_thr 1e-8` |
| `relax` | Ionic relaxation | `force_thr_ev 0.01`, `relax_method cg` |
| `cell-relax` | Cell relaxation | `stress_thr 0.5`, `cal_stress 1` |
| `md` | Molecular dynamics | `md_type nvt`, `md_dt 1.0`, `md_tfirst 300` |
| `band` | Band structure | Generates SCF + NSCF + Band sequence |

**Example**:
```bash
# SCF calculation (default)
abacustest model inputs -f Fe.cif

# Relaxation
abacustest model inputs -f Fe.cif --jtype relax

# Molecular dynamics
abacustest model inputs -f Fe.cif --jtype md
```

---

## Parameter Reference

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `-f` | Structure file(s) | `-f Fe.cif` |
| `--ftype` | File type | `cif`, `poscar`, `abacus` |

### Calculation Settings

| Parameter | Description | Options | Default |
|-----------|-------------|---------|---------|
| `--jtype` | Calculation type | `scf`, `relax`, `cell-relax`, `md`, `band` | `scf` |
| `--lcao` | Use LCAO basis | - | No (PW) |
| `--kpt` | K-point setting | `6 6 6` or single value | Auto |

### Spin & Magnetic Settings

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--nspin` | Spin type: `1`=no spin, `2`=spin-polarized, `4`=non-collinear+SOC | `--nspin 2` |
| `--soc` | Enable spin-orbit coupling | `--soc` (requires `--nspin 4`) |
| `--init_mag` | Initial magnetic moments (element value pairs) | `--init_mag Fe 2.0 Mn -1.5` |
| `--afm` | Antiferromagnetic ordering (alternates signs) | `--afm` |

**Examples**:
```bash
# Spin-polarized with magnetic moment
abacustest model inputs -f Fe.cif --nspin 2 --init_mag Fe 2.5

# Antiferromagnetic (alternating +2.0, -2.0)
abacustest model inputs -f Fe2.cif --nspin 2 --init_mag Fe 2.0 --afm

# SOC calculation (non-collinear)
abacustest model inputs -f Pt.cif --nspin 4 --soc
```

### DFT+U Settings

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--dftu` | Enable DFT+U correction | `--dftu` |
| `--dftu_param` | Hubbard U values (element value pairs) | `--dftu_param Fe 4.0 Ti 1.0` |

**Example**:
```bash
# DFT+U with default U values (d=4eV, f=6eV)
abacustest model inputs -f FeO.cif --dftu

# Specify custom U values
abacustest model inputs -f FeO.cif --dftu --dftu_param Fe 4.0 O 1.0
```

---

## Pseudopotential & Orbital Paths

### Why Paths Matter

- **Pseudopotentials (`.upf`)**: Required for ALL ABACUS calculations. Define the electron-ion interaction.
- **Orbitals (`.orb`)**: Required only for LCAO basis calculations (`--lcao`). Define the atomic orbital basis functions.

### How Paths Are Resolved

The skill looks for pseudopotential/orbital paths in this order:

1. **Command-line arguments** (highest priority):
   ```bash
   abacustest model inputs -f Fe.cif --pp /my/pp/path --orb /my/orb/path
   ```

2. **Environment variables**:
   ```bash
   export ABACUS_PP_PATH=/shared/ABACUS-APNS-v1/pp
   export ABACUS_ORB_PATH=/shared/ABACUS-APNS-v1/orb
   abacustest model inputs -f Fe.cif
   ```

3. **No path provided** (lowest priority):
   - If neither `--pp/--orb` nor environment variables are set, the skill will generate STRU files **without** pseudopotential/orbital paths.
   - You must manually add them to the STRU file before running ABACUS.

### Recommendation

**Set up environment variables once** (see Part 1) so you don't need to specify paths every time:

```bash
# Add to ~/.bashrc or ~/.zshrc
export ABACUS_PP_PATH=/shared/ABACUS-APNS-v1/pp
export ABACUS_ORB_PATH=/shared/ABACUS-APNS-v1/orb
```

---

## Other Options

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--copy_pp_orb` | Copy files instead of symlinking | `--copy_pp_orb` |
| `--folder-syntax` | Custom output directory naming | `'{x[:-4]}'` |
| `--overwrite` | Overwrite existing directories | `--overwrite` |
| `--input` | Use custom INPUT template | `--input my_input.txt` |

**Examples**:
```bash
# Custom directory names (Fe.cif → Fe/)
abacustest model inputs -f *.cif --folder-syntax '{x[:-4]}'

# Copy PP files (for portability)
abacustest model inputs -f Fe.cif --copy_pp_orb

# Overwrite existing output
abacustest model inputs -f Fe.cif --overwrite
```

---

## Common Errors

### Missing Pseudopotentials

```bash
# Warning: some elements have no pseudopotentials
abacustest model inputs -f Fe.cif
```

**Fix**: Set `ABACUS_PP_PATH` or use `--pp`:
```bash
abacustest model inputs -f Fe.cif --pp /path/to/pp
```

### DFT+U Parameter Format

```bash
# Error: parameters should be in pairs
abacustest model inputs -f FeO.cif --dftu --dftu_param Fe
```

**Fix**: Provide element-value pairs:
```bash
abacustest model inputs -f FeO.cif --dftu --dftu_param Fe 4.0
```

### SOC Without nspin=4

```bash
# Warning: SOC requires nspin=4
abacustest model inputs -f Pt.cif --soc --nspin 2
```

**Fix**:
```bash
abacustest model inputs -f Pt.cif --soc --nspin 4
```


