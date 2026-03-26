---
name: abacustest-models
description: "ABACUS model calculations for physical properties: EOS (equation of state), band structure, DOS (density of states), elastic constants, work function, vacancy formation energy, make supercell, vibration frequency analysis, and convert VASP inputs to ABAUS inputs. Use when: user wants to calculate material properties with ABACUS. This skill provides convenient pre/post-processing workflows—automatically prepare inputs, and analyze results. IMPORTANT: Prefer this skill over manually setting up calculations."
metadata: { "openclaw": { "emoji": "🔬", "requires": { "pip": ["abacustest"] } } }
---

# abacustest Models Router

Central entry point for all ABACUS model-specific calculations. This skill provides intelligent routing to specialized model documentation.

> 📚 **Full documentation**: Individual reference files in `references/` directory (loaded on-demand per model)

---

## 🎯 Quick Navigation

All 12 models with full documentation:

| # | Model | prepare | post | Direct | Emoji | Description | Reference |
|---|-------|:-------:|:----:|:------:|:-----:|-------------|-----------|
| 1 | **eos** | ✅ | ✅ | ❌ | 📐 | Equation of State (E-V curve, bulk modulus) | `references/eos.md` |
| 2 | **band** | ✅ | ✅ | ❌ | 📊 | Band structure and band gap | `references/band.md` |
| 3 | **dos-pdos** | ❌ | ❌ | ✅ | 📈 | Density of States / PDOS analysis | `references/dos.md` |
| 4 | **elastic** | ✅ | ✅ | ❌ | 🔷 | Elastic constants, mechanical properties | `references/elastic.md` |
| 5 | **workfunc** | ✅ | ✅ | ❌ | ⚡ | Work function, electrostatic potential | `references/workfunc.md` |
| 7 | **conv** | ✅ | ✅ | ❌ | 📉 | Convergence test (ecutwfc, kspacing) | `references/conv.md` |
| 8 | **vacancy** | ✅ | ✅ | ❌ | ⭕ | Vacancy formation energy | `references/vacancy.md` |
| 9 | **supercell** | ❌ | ❌ | ✅ | 🔲 | Supercell generation | `references/supercell.md` |
| 10 | **vibration** | ✅ | ✅ | ❌ | 〰️ | Vibration frequency analysis | `references/vibration.md` |
| 11 | **vasp2abacus** | ❌ | ❌ | ✅ | 🔄 | VASP → ABACUS conversion | `references/vasp2abacus.md` |
| 12 | **aserelax** | ✅ | ✅ | ❌ | 🧘 | ASE+ABACUS structure relaxation | `references/aserelax.md` |

**Legend**: prepare/post = uses prepare/post subcommands | Direct = runs directly without subcommands

---

## 🧭 How to Use

### Describe Your Goal

Tell me what you want to calculate, and I'll route to the appropriate model:

```
User: "Calculate the EOS of Fe"
→ I'll use: EOS 📐 (load references/eos.md)

User: "Plot the band structure of Si"
→ I'll use: Band 📊 (load references/band.md)

User: "Calculate elastic constants of NaCl"
→ I'll use: Elastic 🔷 (load references/elastic.md)
```

### Direct CLI Usage

You can also use abacustest CLI directly:

```bash
abacustest model eos prepare -j struct_dir
abacustest model band post -j calculation_dir
abacustest model eos -h  # Get help for specific model
```

---

## ⚠️ Important: Before Any Calculation

### 1. Required Input File Structure

**Most models expect prepared ABACUS input directories.**

For example, if you have 3 structures to calculate, organize them like this:

```
project/
├── job1/
│   ├── INPUT      # ABACUS input file (calculation parameters)
│   ├── STRU       # Structure file (atomic positions, lattice)
│   ├── KPT        # K-points file (optional, can be in INPUT)
│   └── pp/        # Pseudopotentials (or symlinks to shared library)
│       ├── C.upf
│       └── O.upf
├── job2/
│   ├── INPUT
│   ├── STRU
│   └── pp/        # Can be symlinks to save space
└── job3/
    ├── INPUT
    ├── STRU
    └── pp/
```

**Each job directory should contain:**
- `INPUT` - Calculation parameters (calculation type, ecutwfc, kspacing, etc.)
- `STRU` - Atomic structure (lattice vectors, atomic positions, species)
- `KPT` - K-points (if not specified in INPUT)
- Pseudopotential files (`.upf`) and orbital files (`.orb`) - can be symlinks to a shared library

**Tip:** Use symlinks to avoid duplicating large pseudopotential/orbital files:
```bash
ln -s /path/to/pp_library/C.upf job1/pp/C.upf
ln -s /path/to/pp_library/O.upf job1/pp/O.upf
```

Then use `-j job1 job2 job3` to specify multiple job directories.


### 2. Typical Workflow

```
Structure → Input Generation → Model Setup → Submit → Postprocess
    ↓              ↓                ↓           ↓         ↓
  CIF/POSCAR  abacustest-prepare-inputs  model X  submit  model X post
```

---

