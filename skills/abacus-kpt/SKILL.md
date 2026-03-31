---
name: abacus-kpt
description: "Write, modify, and validate ABACUS KPT files. Use when: user needs to create KPT files for SCF/NSCF calculations, configure k-point grids, set up band structure paths, or troubleshoot k-point sampling issues."
metadata: { "openclaw": { "emoji": "📐", "requires": { "knowledge": ["Brillouin zone", "Monkhorst-Pack sampling", "band structure calculations"] } } }
---

# ABACUS KPT File Format

The `KPT` file defines the k-point sampling for ABACUS calculations. It can either specify k-points explicitly or generate k-mesh automatically using the Monkhorst-Pack method.

## ⚠️ When to Use This Skill

✅ **Use this skill**:
- "Create a KPT file for SCF calculation with 8×8×8 k-mesh"
- "Set up k-points for band structure calculation along Γ-X-L-W-K"
- "Convert explicit k-points to MP grid"
- "Add k-point shift for metallic systems"
- "Fix KPT file format errors"

❌ **Do not use**:
- Generate k-points automatically → Use `kspacing` in INPUT file
- Convert CIF/POSCAR → Use `abacustest-prepare-inputs` skill
- High-throughput k-point testing → Use `abacustest-prepare` skill

---

## File Structure Overview

A KPT file has three modes:

| Mode | Use Case | Format |
|------|----------|--------|
| **Auto MP Grid** | SCF/NSCF calculations | `0` + `Gamma/MP` + mesh parameters |
| **Explicit K-points** | Special sampling | `N` + `Direct/Cartesian` + coordinates + weights |
| **Line Mode** | Band structure | `N` + `Line` + high-symmetry path |

### Basic Structure

```
K_POINTS          # Keyword (can be K_POINTS, KPOINTS, or K)
<N>               # Number of k-points (0 = auto-generate)
<Type>            # Gamma, MP, Direct, Cartesian, Line, Line_Cartesian
<parameters...>   # Mesh size, coordinates, or path
```

---

## Mode 1: Auto-Generate K-Mesh (SCF/NSCF)

Automatically generates Monkhorst-Pack k-point grids for SCF and NSCF calculations.

### Format

```
K_POINTS
0
Gamma
nx ny ny sx sy sz
```

| Field | Type | Description |
|-------|------|-------------|
| Line 1 | keyword | `K_POINTS`, `KPOINTS`, or `K` |
| Line 2 | int | `0` = auto-generate k-mesh |
| Line 3 | str | `Gamma` (Γ-centered) or `MP` (Monkhorst-Pack) |
| Line 4 | 6 values | `nx ny nz` (grid divisions) + `sx sy sz` (shifts) |

### Monkhorst-Pack Methods

| Type | Description | When to Use |
|------|-------------|-------------|
| `Gamma` | Γ-centered MP grid | Insulators, semiconductors, most systems |
| `MP` | Standard MP grid | Metals, special sampling needs |

### Grid Divisions (nx, ny, nz)

Number of k-points along each reciprocal lattice vector:
- **Larger values** = denser sampling = better accuracy = more computational cost
- **Typical values**: 4–12 for SCF, 15–30 for DOS

### Grid Shifts (sx, sy, sz)

Offset of the k-grid (real numbers, typically 0 or 0.5):
- `0 0 0`: No shift (standard)
- `0.5 0.5 0.5`: Shifted grid (sometimes better for metals)

### Examples

**Standard 6×6×6 Γ-centered grid:**
```
K_POINTS
0
Gamma
6 6 6 0 0 0
```

**8×8×8 MP grid with no shift:**
```
K_POINTS
0
MP
8 8 8 0 0 0
```

**Shifted 4×4×4 grid (for metals):**
```
K_POINTS
0
MP
4 4 4 0.5 0.5 0.5
```

**Asymmetric grid (layered materials):**
```
K_POINTS
0
Gamma
12 12 4 0 0 0
```

**Slab/Surface (vacuum in c-direction):**
```
K_POINTS
0
Gamma
12 12 1 0 0 0
```

---

## Mode 2: Explicit K-Points

Manually specify k-point coordinates and weights for special sampling.

### Format

```
K_POINTS
N
Direct
k1x k1y k1z w1
k2x k2y k2z w2
...
kNx kNy kNz wN
```

| Field | Type | Description |
|-------|------|-------------|
| Line 1 | keyword | `K_POINTS`, `KPOINTS`, or `K` |
| Line 2 | int | Number of explicit k-points |
| Line 3 | str | `Direct` (fractional) or `Cartesian` (Bohr⁻¹) |
| Lines 4+ | 4 values | `kx ky kz weight` |

### Coordinate Types

| Type | Unit | Description |
|------|------|-------------|
| `Direct` | fractional | Fractional coordinates in reciprocal space (0–1) |
| `Cartesian` | Bohr⁻¹ | Cartesian coordinates in reciprocal space |

### Weights

- **Normalized**: Sum of all weights should equal 1.0
- **Symmetry**: ABACUS can use symmetry to reduce k-points; weights account for this

### Example: 2×2×2 MP Grid (Explicit)

```
K_POINTS
8
Direct
0.0 0.0 0.0 0.125
0.5 0.0 0.0 0.125
0.0 0.5 0.0 0.125
0.5 0.5 0.0 0.125
0.0 0.0 0.5 0.125
0.5 0.0 0.5 0.125
0.0 0.5 0.5 0.125
0.5 0.5 0.5 0.125
```

---

## Mode 3: Line Mode (Band Structure)

Specifies high-symmetry k-point paths for band structure calculations.

### Format

```
K_POINTS
N
Line
k1x k1y k1z n1  # Label1
k2x k2y k2z n2  # Label2
...
kNx kNy kNz nN  # LabelN
```

| Field | Type | Description |
|-------|------|-------------|
| Line 1 | keyword | `K_POINTS`, `KPOINTS`, or `K` |
| Line 2 | int | Number of high-symmetry points |
| Line 3 | str | `Line` (Direct) or `Line_Cartesian` |
| Lines 4+ | 4+ values | `kx ky kz n #Label` |

### Parameters

| Field | Type | Description |
|-------|------|-------------|
| `kx ky kz` | real | High-symmetry point coordinates (fractional) |
| `n` | int | Number of interpolation points to next point |
| `#Label` | str (optional) | High-symmetry point name (for plotting) |

### Important Rules

1. **Last point**: Set `n=1` for the final high-symmetry point (path ends here)
2. **Interpolation**: ABACUS interpolates `n` k-points between consecutive points
3. **Path continuity**: Each point connects to the next; plan your path carefully

### High-Symmetry Point Labels

Common labels for cubic systems (FCC):

| Label | Coordinates | Description |
|-------|-------------|-------------|
| Γ (G) | (0, 0, 0) | Brillouin zone center |
| X | (0.5, 0, 0.5) | [100] boundary |
| L | (0.5, 0.5, 0.5) | [111] boundary |
| W | (0.5, 0.25, 0.75) | Edge intersection |
| K | (0.375, 0.375, 0.75) | [110] boundary |
| U | (0.625, 0.25, 0.625) | Edge midpoint |

### Example: FCC Band Path (Γ-X-U-K-Γ-L-W-X)

For standard paths, use [SeeK-path](https://www.materialscloud.org/work/tools/seekpath).

```
K_POINTS
8
Line
0.000 0.000 0.000 20  # G
0.500 0.000 0.500 20  # X
0.625 0.250 0.625  1  # U
0.375 0.375 0.750 20  # K
... (continue to G, L, W, X)
0.500 0.000 0.500  1  # X (last point, n=1)
```

**Note**: Last point must have `n=1`. Use 20–40 points per segment for smooth bands.

### Example: Simple Path (Γ-X-M-Γ-R)

```
K_POINTS
5
Line
0.0 0.0 0.0 30  # G
0.5 0.0 0.0 30  # X
0.5 0.5 0.0 30  # M
0.0 0.0 0.0 30  # G
0.5 0.5 0.5  1  # R
```

### ⚠️ Tips for Band Paths

- **Interpolation points**: Use 20–40 points per segment for smooth bands
- **Path planning**: Include all relevant high-symmetry points for your lattice type
- **Reference**: Use [SeeK-path](https://www.materialscloud.org/work/tools/seekpath) to find standard paths for your structure

---

## ⚠️ INPUT Parameters That Override KPT

### gamma_only

For LCAO calculations, you can use `gamma_only 1` in INPUT file:

```
# INPUT file
gamma_only 1
```

**Effect**: ABACUS ignores the KPT file and uses only the Γ point.

**When to use**:
- Large supercells (Γ-point sampling sufficient)
- LCAO basis calculations
- Quick tests

**⚠️ Warning**: If `gamma_only 1`, the KPT file is **overwritten**. Make sure to set `gamma_only 0` for multi-k calculations.

### kspacing

When `kspacing` is set in INPUT file:

```
# INPUT file
kspacing 0.14  # unit: 1/Bohr
```

**Effect**: ABACUS **automatically generates** an appropriate k-mesh based on the reciprocal lattice and ignores the KPT file.

**How it works**:
- ABACUS calculates the required k-point density from `kspacing` and the reciprocal lattice vectors
- A Γ-centered MP grid is generated automatically
- The generated k-mesh **overwrites** the KPT file

**Single value (isotropic)**:
```
kspacing 0.14  # Same spacing in all directions
```

**Three values (anisotropic)**:
```
kspacing 0.1 0.1 1.0  # Different spacing for a/b/c directions
```

**Example: 2D Material (Graphene)**

For 2D materials with vacuum layer in c-direction, use larger spacing in vacuum direction to get 1 k-point:

```
# INPUT
kspacing 0.1 0.1 1.0  # Dense in-plane, 1 k-point in vacuum

# KPT (will be overwritten!)
K_POINTS
0
Gamma
12 12 1 0 0 0  # ❌ This is ignored; ABACUS generates ~12×12×1 grid
```

**Benefits for 2D systems**:
- Automatic determination of in-plane k-point density
- Single k-point in vacuum direction (no waste)
- Consistent sampling across different supercell sizes

**When to use `kspacing`**:
- ✅ High-throughput calculations (consistent k-point density across structures)
- ✅ Structure relaxations (automatic adjustment with cell changes)
- ✅ Quick setup (no need to manually determine grid size)
- ✅ 2D materials/vacuum systems (e.g., `0.1 0.1 1.0` for single k-point in vacuum)

**When to use explicit KPT**:
- ✅ Band structure calculations (Line mode required)
- ✅ Special k-point sampling (explicit coordinates needed)
- ✅ Reproducing published results (exact k-points required)
- ✅ NSCF DOS/band calculations (specific dense grid needed)

**⚠️ Warning**: Do not rely on KPT file contents when `kspacing` is set in INPUT. The KPT file will be regenerated by ABACUS.

---

## K-Point Convergence Guidelines

### SCF Calculations

| System Type | Recommended Grid | Notes |
|-------------|------------------|-------|
| Metals | 12×12×12 or denser | Need dense sampling near E_F |
| Semiconductors | 8×8×8 to 12×12×12 | Moderate density sufficient |
| Insulators | 6×6×6 to 8×8×8 | Can use coarser grids |
| Layered materials | 12×12×4 | Dense in-plane, sparse out-of-plane |
| Nanowires | 4×12×12 | Sparse along wire axis |
| Molecules (supercell) | 2×2×2 or Γ-only | Large cell = small BZ |

### DOS Calculations

| Property | Recommended Grid |
|----------|------------------|
| TDOS | 1.5–2× SCF grid density |
| PDOS | 2× SCF grid density |

### Band Structure

| Segment | Recommended Points |
|---------|-------------------|
| Standard path | 20–30 points |
| Complex dispersion | 40–50 points |
| Quick test | 10–15 points |

---

## Common Errors and Fixes

### Error 1: Wrong Number of K-Points

```
K_POINTS
0
Gamma
6 6 6  # ❌ Missing shift values
```

**Fix**:
```
K_POINTS
0
Gamma
6 6 6 0 0 0  # ✅ All 6 values required
```

### Error 2: Line Mode Last Point

```
K_POINTS
5
Line
0.0 0.0 0.0 30
0.5 0.0 0.0 30
0.5 0.5 0.0 30
0.0 0.0 0.0 30
0.5 0.5 0.5 30  # ❌ Last point should have n=1
```

**Fix**:
```
K_POINTS
5
Line
0.0 0.0 0.0 30
0.5 0.0 0.0 30
0.5 0.5 0.0 30
0.0 0.0 0.0 30
0.5 0.5 0.5 1  # ✅ n=1 for last point
```

### Error 3: Weights Don't Sum to 1

```
K_POINTS
4
Direct
0.0 0.0 0.0 0.5
0.5 0.0 0.0 0.5
0.0 0.5 0.0 0.5
0.0 0.0 0.5 0.5  # ❌ Sum = 2.0
```

**Fix**:
```
K_POINTS
4
Direct
0.0 0.0 0.0 0.25
0.5 0.0 0.0 0.25
0.0 0.5 0.0 0.25
0.0 0.0 0.5 0.25  # ✅ Sum = 1.0
```

### Error 4: Gamma-Only Conflict

```
# INPUT
gamma_only 1

# KPT
K_POINTS
0
Gamma
8 8 8 0 0 0  # ❌ Ignored when gamma_only=1
```

**Fix**: Either remove KPT or set `gamma_only 0`.

---

## Tips

| Recommendation | Reason |
|----------------|--------|
| Use `Gamma` for most systems | Γ-centered grids are more efficient |
| Test k-point convergence | Energy/forces change with grid density |
| Use asymmetric grids for low-D systems | Save computation on vacuum direction |
| Plan band paths carefully | Include all relevant high-symmetry points |
| Use 20–30 points per band segment | Smooth curves without excessive cost |
| Check SeeK-path for standard paths | Avoid missing important k-points |

---

## Related Skills

- **STRU file**: [`abacus-stru`](../abacus-stru/SKILL.md)
- **INPUT parameters**: [`abacus-input-parameter`](../abacus-input-parameter/SKILL.md)
- **Input preparation**: [`abacustest-prepare-inputs`](../abacustest-prepare-inputs/SKILL.md)
- **Parameter mixing**: [`abacustest-prepare`](../abacustest-prepare/SKILL.md)

---

## External Resources

- **ABACUS KPT Docs**: https://abacus.deepmodeling.com/en/latest/advanced/input_files/kpt.html
- **SeeK-path** (high-symmetry paths): https://www.materialscloud.org/work/tools/seekpath
- **Brillouin Zone Database**: https://www.cryst.ehu.es/cryst/help/definitions.html
- **Monkhorst-Pack Method**: Monkhorst & Pack, Phys. Rev. B 13, 5188 (1976)
