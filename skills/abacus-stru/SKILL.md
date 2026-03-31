---
name: abacus-stru
description: "Write, modify, and validate ABACUS STRU structure files. Use when: user needs to create a new STRU file from scratch, modify atomic positions/lattice parameters, add magnetic moments or velocity constraints, or troubleshoot STRU format issues."
metadata: { "openclaw": { "emoji": "🔬", "requires": { "knowledge": ["ABACUS STRU format", "crystallography basics"] } } }
---

# ABACUS STRU File Format

The `STRU` file contains the crystal structure information for ABACUS calculations, including lattice geometry, atomic positions, pseudopotential files, and numerical orbital files (for LCAO calculations).

## ⚠️ When to Use This Skill

✅ **Use this skill**:
- "Create a STRU file for Si diamond structure"
- "Add magnetic moments to Fe atoms in this STRU"
- "Convert lattice vectors from Angstrom to Bohr"
- "Fix the STRU format error in this file"
- "Add velocity constraints for MD simulation"

❌ **Do not use**:
- Convert CIF/POSCAR → STRU → Use `abacustest-prepare-inputs` skill
- Modify INPUT parameters → Use `abacus-input-parameter` skill
- Generate high-throughput inputs → Use `abacustest-prepare` skill

---

## File Structure Overview

A STRU file consists of several sections, each starting with a keyword:

```
ATOMIC_SPECIES
[element definitions]

NUMERICAL_ORBITAL      # LCAO only
[orbital files]

LATTICE_CONSTANT
[lattice scaling factor in Bohr]

LATTICE_VECTORS        # Or use latname in INPUT
[3x3 lattice vector matrix]

LATTICE_PARAMETERS     # Only with latname
[lattice parameters]

ATOMIC_POSITIONS
[coordinate type]
[element blocks with atomic positions]
```

---

## Section 1: ATOMIC_SPECIES

Defines chemical elements, atomic masses, and pseudopotential files.

### Format

```
ATOMIC_SPECIES
Label  Mass  PseudoFile  [PseudoType]
```

| Field | Type | Description | Required |
|-------|------|-------------|----------|
| `Label` | str | Element label (e.g., `Si`, `Fe`) | ✅ |
| `Mass` | float | Atomic mass (amu), used in MD only | ✅ |
| `PseudoFile` | str | Pseudopotential filename | ✅ |
| `PseudoType` | str | Pseudopotential format | ❌ (default: `auto`) |

### PseudoType Options

| Value | Format | Source |
|-------|--------|--------|
| `upf` | .UPF | Quantum ESPRESSO |
| `upf201` | .UPF (new format) | Quantum ESPRESSO |
| `vwr` | .vwr | Vanderbilt |
| `blps` | bulk-derived local | Princeton BLPS |
| `auto` | auto-detect | Default |

### Examples

**Single element:**
```
ATOMIC_SPECIES
Si 28.085 Si_ONCV_PBE-1.0.upf upf201
```

**Multiple elements:**
```
ATOMIC_SPECIES
Fe 55.845 Fe.upf
O 15.999 O.upf
```

**With explicit paths:**
```
ATOMIC_SPECIES
Mg 24.305 /path/to/pseudo/Mg.upf
O 15.999 /path/to/pseudo/O.upf
```

### ⚠️ Important Notes

- **Mass value**: Only used in molecular dynamics; any reasonable value works for SCF
- **Pseudopotential path**: If not specified, file is assumed in working directory
- **Path resolution**: Pseudopotential paths in STRU are **relative to `pseudo_dir` in INPUT**. For example:
  - INPUT: `pseudo_dir a/b/`
  - STRU: `Si 14 c/d/Si.upf`
  - ABACUS looks for: `a/b/c/d/Si.upf`
- **XC functional**: All pseudopotentials must use the same XC functional, or set `dft_functional` explicitly in INPUT

### Common Pseudopotential Sources

1. [Quantum ESPRESSO](http://www.quantum-espresso.org/pseudopotentials/)
2. [SG15-ONCV](http://quantum-simulation.org/potentials/sg15_oncv/upf/)
3. [PseudoDojo](http://www.pseudo-dojo.org/)
4. [BLPS Library](https://github.com/PrincetonUniversity/BLPSLibrary)

---

## Section 2: NUMERICAL_ORBITAL

Defines numerical atomic orbital files. Ignored in plane-wave (`basis_type pw`) calculations for most cases.

### Format

```
NUMERICAL_ORBITAL
OrbitalFile1
OrbitalFile2
...
```

One orbital file per element type, in the same order as `ATOMIC_SPECIES`.

### Example

```
NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb
```

### Orbital File Naming Convention

Typical format: `{Element}_{xc}_{cutoff}Ry_{orbitals}.orb`

Example: `Si_gga_8au_60Ry_2s2p1d.orb`
- `gga`: XC functional
- `8au`: Cutoff radius (8 a.u.)
- `60Ry`: Energy cutoff (60 Ry)
- `2s2p1d`: Orbital basis (2 s, 2 p, 1 d)

### ⚠️ Important Notes

- **LCAO only**: This section is ignored for `basis_type pw` **except** in special cases:
  - **DFT+U calculations** (`dft_plus_u 1`): Orbital files needed for projection even with `basis_type pw`
  - **Projected magnetic moments** (`onsite_radius` set): Orbital files needed for atomic-projected magnetism
  - See [`abacus-input-parameter`](../abacus-input-parameter/SKILL.md) for `dft_plus_u` and `onsite_radius` details
- **Path resolution**: Orbital file paths are **relative to `orbital_dir` in INPUT** (same as pseudopotentials)
- **Order matters**: Orbital files must match element order in `ATOMIC_SPECIES`
- **Download**: Orbitals available from [ABACUS official](http://abacus.ustc.edu.cn/pseudo/list.htm)

---

## Section 3: LATTICE_CONSTANT

Defines the lattice scaling factor in **Bohr** units.

### Format

```
LATTICE_CONSTANT
<Value>
```

### Unit Conversion

| From | To Bohr | Example |
|------|---------|---------|
| Angstrom | × 1.889726 | 5.43 Å → 10.26 Bohr |
| Bohr | × 1.0 | 10.26 Bohr → 10.26 Bohr |

### Example

```
LATTICE_CONSTANT
10.26  # 5.43 Angstrom
```

---

## Section 4: LATTICE_VECTORS

Defines the 3×3 lattice vector matrix. **Vectors are scaled by LATTICE_CONSTANT**.

### Format

```
LATTICE_VECTORS
v1_x  v1_y  v1_z  # latvec1
v2_x  v2_y  v2_z  # latvec2
v3_x  v3_y  v3_z  # latvec3
```

### Example: Si Diamond (FCC primitive cell)

```
LATTICE_CONSTANT
10.2

LATTICE_VECTORS
0.5 0.5 0.0  # latvec1
0.5 0.0 0.5  # latvec2
0.0 0.5 0.5  # latvec3
```

Actual lattice vectors = LATTICE_CONSTANT × matrix values

### ⚠️ Alternative: Use `latname` in INPUT

Instead of manually specifying `LATTICE_VECTORS`, you can use the `latname` parameter in INPUT file:

```
# In INPUT file:
latname fcc
```

Then **remove** the `LATTICE_VECTORS` section from STRU. See [LATTICE_PARAMETERS](#section-5-lattice_parameters) for complex lattices.

---

## Section 5: LATTICE_PARAMETERS

Used **only** with `latname` in INPUT file. Provides additional lattice parameters for non-cubic systems.

### Format

```
LATTICE_PARAMETERS
param1  param2  ...
```

### Supported Lattice Types

| latname | Parameters | Description | Generated Vectors |
|---------|------------|-------------|-------------------|
| `sc` | none | Simple cubic | v1=(1,0,0), v2=(0,1,0), v3=(0,0,1) |
| `fcc` | none | Face-centered cubic | v1=(-0.5,0,0.5), v2=(0,0.5,0.5), v3=(-0.5,0.5,0) |
| `bcc` | none | Body-centered cubic | v1=(0.5,0.5,0.5), v2=(-0.5,0.5,0.5), v3=(-0.5,-0.5,0.5) |
| `hexagonal` | c/a | c/a ratio | v1=(1,0,0), v2=(-0.5,√3/2,0), v3=(0,0,c/a) |
| `trigonal` | cos(γ) | Angle cosine | See formula below |
| `st` | c/a | Simple tetragonal | v1=(1,0,0), v2=(0,1,0), v3=(0,0,c/a) |
| `bct` | c/a | Body-centered tetragonal | v1=(0.5,-0.5,c/a), v2=(0.5,0.5,c/a), v3=(-0.5,-0.5,c/a) |
| `so` | b/a, c/a | Simple orthorhombic | v1=(1,0,0), v2=(0,b/a,0), v3=(0,0,c/a) |
| `baco` | b/a, c/a | Base-centered orthorhombic | v1=(0.5,b/a/2,0), v2=(-0.5,b/a/2,0), v3=(0,0,c/a) |
| `fco` | b/a, c/a | Face-centered orthorhombic | v1=(0.5,0,c/a/2), v2=(0.5,b/a/2,0), v3=(0,b/a/2,c/a/2) |
| `bco` | b/a, c/a | Body-centered orthorhombic | v1=(0.5,b/a/2,c/a/2), v2=(-0.5,b/a/2,c/a/2), v3=(-0.5,-b/a/2,c/a/2) |
| `sm` | b/a, c/a, cos(αβ) | Simple monoclinic | See formula below |
| `bacm` | b/a, c/a, cos(αβ) | Base-centered monoclinic | See formula below |
| `triclinic` | b/a, c/a, cos(ab), cos(ac), cos(bc) | Triclinic | See formula below |

### Trigonal Lattice Formula

For `latname = "trigonal"`, parameter x = cos(γ):

```
tx = sqrt((1-x)/2)
ty = sqrt((1-x)/6)
tz = sqrt((1+2*x)/3)

v1 = (tx, -ty, tz)
v2 = (0, 2*ty, tz)
v3 = (-tx, -ty, tz)
```

### Example: Hexagonal Lattice

**INPUT file:**
```
latname hexagonal
```

**STRU file:**
```
LATTICE_CONSTANT
5.0

LATTICE_PARAMETERS
1.63  # c/a ratio
```

---

## Section 6: ATOMIC_POSITIONS

Defines atomic positions and associated parameters (magnetism, constraints, velocities).

### Format

```
ATOMIC_POSITIONS
<CoordinateType>
<ElementLabel>
<Magnetism>
<NumAtoms>
<x> <y> <z> [m <move_x> <move_y> <move_z>] [mag <magmom>] [v <vx> <vy> <vz>] ...
```

### CoordinateType Options

| Value | Unit | Description |
|-------|------|-------------|
| `Direct` | fractional | Fractional coordinates (0-1) |
| `Cartesian` | LATTICE_CONSTANT | Cartesian in lattice constant units |
| `Cartesian_au` | Bohr | Cartesian in Bohr (LATTICE_CONSTANT=1.0) |
| `Cartesian_angstrom` | Angstrom | Cartesian in Angstrom (LATTICE_CONSTANT=1.889726) |
| `Cartesian_angstrom_center_xy` | Angstrom | Centered at (0.5, 0.5, 0.0) |
| `Cartesian_angstrom_center_xz` | Angstrom | Centered at (0.5, 0.0, 0.5) |
| `Cartesian_angstrom_center_yz` | Angstrom | Centered at (0.0, 0.5, 0.5) |
| `Cartesian_angstrom_center_xyz` | Angstrom | Centered at (0.5, 0.5, 0.5) |

### Element Block Structure

For each element type:

1. **Line 1**: Element label (must match `ATOMIC_SPECIES`)
2. **Line 2**: Default magnetic moment (μ_B) for this element
3. **Line 3**: Number of atoms of this element
4. **Lines 4+**: Atomic positions with optional parameters

### Atomic Position Parameters

After coordinates, add parameters using keywords:

| Keyword | Values | Description |
|---------|--------|-------------|
| `m` (or none) | 0/1 0/1 0/1 | Movement constraints for relaxation |
| `v` / `vel` / `velocity` | vx vy vz | Initial velocity |
| `mag` / `magmom` | m (collinear) or mx my mz (non-collinear) | Atomic magnetic moment |
| `angle1` | degrees | Polar angle (non-collinear, z-axis) |
| `angle2` | degrees | Azimuthal angle (non-collinear, xy-plane) |

### ⚠️ Non-Collinear Magnetism: mag + angle1/angle2

When **both** `mag` (3 components) **and** `angle1`/`angle2` are specified in non-collinear calculations (`nspin=4`):

1. ABACUS first calculates the **magnitude** from the `mag` vector: `|mag| = sqrt(mx² + my² + mz²)`
2. Then uses `angle1` and `angle2` to define the **direction**
3. Final magnetic moment = `|mag|` in the direction specified by angles

**Example:**
```
Fe
1.0
1
0.0 0.0 0.0 mag 1.0 1.0 1.0 angle1 90 angle2 0
```
- `mag 1.0 1.0 1.0` → magnitude = √3 ≈ 1.73 μ_B
- `angle1 90 angle2 0` → direction in xy-plane along x-axis
- Final moment: 1.73 μ_B along +x direction

### Examples

**Simple Si diamond (Direct coordinates):**
```
ATOMIC_POSITIONS
Direct
Si
0.0
2
0.00 0.00 0.00 0 0 0
0.25 0.25 0.25 1 1 1
```

**With magnetic moments (collinear):**
```
ATOMIC_POSITIONS
Direct
Fe
1.0
2
0.0 0.0 0.0 m 0 0 0 mag 2.5
0.5 0.5 0.5 m 1 1 1 mag -2.5
```

**With velocities (MD):**
```
ATOMIC_POSITIONS
Cartesian_angstrom
Mg
0.0
1
0.0 0.0 0.0 v 1.0 0.5 0.0
```

**Non-collinear magnetism:**
```
ATOMIC_POSITIONS
Direct
Fe
1.0
2
0.0 0.0 0.0 m 0 0 0 mag 1.0 angle1 90 angle2 0
0.5 0.5 0.5 m 1 1 1 mag 1.0 angle1 90 angle2 180
```

### Movement Constraints (`m` keyword)

Three values (0 or 1) control atom movement in relaxation:

- `0 0 0`: Atom fixed in all directions
- `1 1 1`: Atom free to move in all directions
- `1 0 0`: Atom moves only in x direction

**Example:**
```
0.00 0.00 0.00 0 0 0  # Fixed atom
0.25 0.25 0.25 1 1 1  # Free atom
```

### Magnetic Moment Settings

**Default behavior:**

| nspin | Default magmom if not specified |
|-------|--------------------------------|
| `nspin=2` (collinear) | Auto-set to 1.0 μ_B |
| `nspin=4` (non-collinear) | Auto-set to (1, 1, 1) |

**Override:** If any atom has explicit `mag` keyword, auto-setting is disabled.

---

## Complete Examples

### Example 1: Si Diamond (FCC, no latname)

```
ATOMIC_SPECIES
Si 28.085 Si_ONCV_PBE-1.0.upf upf201

NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb

LATTICE_CONSTANT
10.2

LATTICE_VECTORS
0.5 0.5 0.0
0.5 0.0 0.5
0.0 0.5 0.5

ATOMIC_POSITIONS
Direct
Si
0.0
2
0.00 0.00 0.00 0 0 0
0.25 0.25 0.25 1 1 1
```

### Example 2: Si Diamond (using latname)

**INPUT file:**
```
latname fcc
```

**STRU file:**
```
ATOMIC_SPECIES
Si 28.085 Si_ONCV_PBE-1.0.upf

NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb

LATTICE_CONSTANT
10.2

ATOMIC_POSITIONS
Direct
Si
0.0
2
0.00 0.00 0.00 0 0 0
0.25 0.25 0.25 1 1 1
```

### Example 3: FeO (Antiferromagnetic)

```
ATOMIC_SPECIES
Fe 55.845 Fe.upf
O 15.999 O.upf

LATTICE_CONSTANT
8.0

LATTICE_VECTORS
0.5 0.5 0.0
0.5 0.0 0.5
0.0 0.5 0.5

ATOMIC_POSITIONS
Direct
Fe
1.0
2
0.0 0.0 0.0 m 1 1 1 mag 4.0
0.5 0.5 0.5 m 1 1 1 mag -4.0

O
0.0
2
0.5 0.0 0.0 m 1 1 1
0.0 0.5 0.5 m 1 1 1
```

### Example 4: MgO (Cartesian Angstrom)

```
ATOMIC_SPECIES
Mg 24.305 Mg.upf
O 15.999 O.upf

LATTICE_CONSTANT
1.889726

LATTICE_VECTORS
4.21 0.0 0.0
0.0 4.21 0.0
0.0 0.0 4.21

ATOMIC_POSITIONS
Cartesian_angstrom
Mg
0.0
1
0.0 0.0 0.0 m 1 1 1

O
0.0
1
2.105 2.105 2.105 m 1 1 1
```

---

## Common Errors and Fixes

### Error 1: Element Label Mismatch

```
ATOMIC_SPECIES
Si 28.085 Si.upf

ATOMIC_POSITIONS
Direct
si  # ❌ Wrong: case mismatch
```

**Fix:**
```
ATOMIC_POSITIONS
Direct
Si  # ✅ Match ATOMIC_SPECIES exactly
```

### Error 2: Missing NUMERICAL_ORBITAL for LCAO

```
# INPUT: basis_type lcao
# STRU: (missing NUMERICAL_ORBITAL section)
```

**Fix:**
```
NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb
```

### Error 3: Wrong Coordinate Type

```
LATTICE_CONSTANT
10.2

ATOMIC_POSITIONS
Cartesian  # ❌ Values are in fractional
0.0 0.0 0.0
0.25 0.25 0.25
```

**Fix:**
```
ATOMIC_POSITIONS
Direct  # ✅ Use Direct for fractional
0.0 0.0 0.0
0.25 0.25 0.25
```

### Error 4: LATTICE_VECTORS with latname

```
# INPUT: latname fcc
# STRU:
LATTICE_VECTORS  # ❌ Conflict with latname
0.5 0.5 0.0
...
```

**Fix:** Remove `LATTICE_VECTORS` section when using `latname`.

---

## Tips

| Recommendation | Reason |
|----------------|--------|
| Use `Direct` coordinates | Easier to read and modify |
| Match element labels exactly | Case-sensitive matching |
| Set `m 0 0 0` for substrate atoms | Fix bottom layers in slab calculations |
| Use `Cartesian_angstrom` for MD | Easier to set velocities in Å/fs |
| Verify pseudopotential paths | Missing PP files cause ABACUS crash |
| Use `latname` for common lattices | Avoids manual vector errors |

---

## Related Skills

- **INPUT parameters**: `abacus-input-parameter`
- **Structure conversion**: `abacustest-prepare-inputs`
- **High-throughput inputs**: `abacustest-prepare`
- **Format conversion**: `abacustest-abacus2VaspQeCp2k`

---

## External Resources

- **ABACUS STRU Docs**: https://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html
- **Pseudopotential Library**: http://abacus.ustc.edu.cn/pseudo/list.htm
- **Crystallography Reference**: https://www.cryst.ehu.es/cryst/help/definitions.html
