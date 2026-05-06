# Supercell - Supercell Generation

Generate supercells from unit cell structures.

---

## Overview

The `supercell` model is a **direct-run utility**:

| Mode | Command | Description |
|------|---------|-------------|
| **prepare** | ❌ | Not supported |
| **post** | ❌ | Not supported |
| **Direct Run** | ✅ | `abacustest model supercell [args]` - Generate supercell directly |

**Note:** This is a **direct-run model** that generates supercells immediately without prepare/post workflow.

---

## Use When

- Generate a supercell of a ABACUS structure STRU

---

## Usage: Direct Run

```bash
# Generate 2×2×2 supercell
abacustest model supercell 2 2 2 -i input.stru -o output.stru

# Generate 3×3×1 supercell (surface slab)
abacustest model supercell 3 3 1 -i slab.stru -o slab_3x3x1.stru

# Multiple structures (batch)
abacustest model supercell 2 2 2 -i struct1.stru -o struct1_2x2x2.stru
```

---

## Input/Output

**Input:** Single structure file (STRU format)
**Output:** Single supercell structure file (STRU format)

```
# Before
unit_cell.stru  →  abacustest model supercell 2 2 2 -i unit_cell.stru -o supercell.stru
                       ↓
# After
supercell.stru (2×2×2 = 8× atoms)
```

---

## Parameters

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `sc` | (positional) | Supercell size in a, b, c directions (3 integers) | Required |
| `--input` | `-i` | Input structure file (STRU format) | Required |
| `--output` | `-o` | Output supercell file name | `{input}_{a}_{b}_{c}` |

**Positional arguments:**
- Three integers specifying supercell multiplication factors along a, b, c lattice vectors

**Examples:**
```bash
# 2×2×2 supercell
abacustest model supercell 2 2 2 -i input.stru

# 3×3×1 supercell (expand in-plane, keep c unchanged)
abacustest model supercell 3 3 1 -i slab.stru

# 2×1×1 supercell (expand only along a)
abacustest model supercell 2 1 1 -i chain.stru

# Custom output filename
abacustest model supercell 2 2 2 -i input.stru -o large_cell.stru

# Default output: input_2_2_2.stru
abacustest model supercell 2 2 2 -i input.stru
```

---

## Output File

**Default naming:** `{input_filename}_{a}_{b}_{c}`

```
# Input: si.stru
# Command: abacustest model supercell 2 2 2 -i si.stru
# Output: si_2_2_2.stru

# Input: slab.stru  
# Command: abacustest model supercell 3 3 1 -i slab.stru -o slab_large.stru
# Output: slab_large.stru
```

**Output STRU format:**
- Same format as input STRU
- Lattice vectors multiplied by supercell factors
- Atomic positions replicated accordingly
- Atom count = original × (a × b × c)

---

## Tips

| Recommendation | Reason |
|----------------|--------|
| Check atom count | Verify: N_super = N_unit × (a×b×c) |
| Surface calculations | Keep vacuum direction unchanged (usually 1) |

---

## Common Issues

| Issue | Solution |
|-------|----------|
| Output file not created | Check input file path; ensure STRU format is valid |
| Wrong atom count | Verify supercell factors; check input structure |
| Lattice vectors wrong | Ensure input STRU has correct lattice definition |
| File format error | Input must be valid ABACUS STRU format |

---

## Related

- Main router: [`abacustest-models`](../SKILL.md)
- Phonon calculations: [`phonon.md`](phonon.md)
- Vacancy (needs supercell): [`vacancy.md`](vacancy.md)
- Vibration (molecules/clusters): [`vibration.md`](vibration.md)
- Input generation: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
