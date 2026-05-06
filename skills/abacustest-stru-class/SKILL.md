---
name: abacustest-stru-class
description: "Programmatically manipulate ABACUS structure files using the AbacusSTRU Python class from abacustest package. Use when: reading/writing STRU files, modifying crystal structures, format conversion (VASP/CIF/ASE/etc.), batch processing structures, fixing atoms, setting magnetic moments, and other programmatic tasks."
metadata: { "openclaw": { "emoji": "🐍", "requires": { "pip": ["abacustest"] } } }
---

# AbacusSTRU Class Usage Guide

`AbacusSTRU` is the core Python class in the `abacustest` package for **programmatically** manipulating ABACUS crystal structures. It provides complete functionality for structure reading, writing, modification, and conversion.

## ⚠️ When to Use This Skill

✅ **Use this skill**:
- "Read and modify STRU files using Python code"
- "Batch convert CIF to STRU format"
- "Fix bottom layers of a surface slab in code"
- "Set magnetic moments for antiferromagnetic systems"
- "Convert structures from ASE/Pymatgen to ABACUS"
- "Apply random perturbations for MD initial structures"

❌ **Do not use**:
- Manually write/edit STRU file format → Use `abacus-stru` skill
- Batch prepare ABACUS inputs via command line → Use `abacustest-prepare-inputs` skill
- Modify INPUT parameters → Use `abacus-input-parameter` skill

---

## Quick Start

### Import Classes

```python
from abacustest.lib_prepare.stru import AbacusSTRU, AbacusATOM
```

### Read Structures

```python
# From STRU file
stru = AbacusSTRU.read("STRU", fmt="stru")

# From POSCAR
stru = AbacusSTRU.read("POSCAR", fmt="poscar")

# From CIF
stru = AbacusSTRU.read("structure.cif", fmt="cif")
```

### Write Structures

```python
# Write to STRU
stru.write("STRU", fmt="stru")

# Write to POSCAR
stru.write("POSCAR", fmt="poscar")

# Write to CIF
stru.write("structure.cif", fmt="cif")
```

### Basic Property Access

```python
print(stru.natoms)           # Total number of atoms
print(stru.cell)             # Cell vectors (Å)
print(stru.labels)           # List of atom labels
print(stru.coords)           # Cartesian coordinates (Å)
print(stru.elements)         # List of element symbols
```

---

## Core Classes

### 1. AbacusSTRU - Structure Class

Represents a complete ABACUS crystal structure.

#### Creation Methods

**Method 1: Read from file**
```python
stru = AbacusSTRU.read("STRU", fmt="stru")
```

**Method 2: Convert from other formats**
```python
# From ASE
stru = AbacusSTRU.from_ase(ase_atoms)

# From Pymatgen
stru = AbacusSTRU.from_pymatgen(pmg_structure)

# From DPData
stru = AbacusSTRU.from_dpdata(dpdata_system, index=0)

# From Phonopy
stru = AbacusSTRU.from_phonopy(phonopy_atoms)
```

**Method 3: Create from scratch**
```python
cell = [[4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
        [0.0, 0.0, 4.0]]

atoms = [
    AbacusATOM(label="Si", coord=(0.0, 0.0, 0.0)),
    AbacusATOM(label="Si", coord=(2.0, 2.0, 2.0)),
]

stru = AbacusSTRU(cell=cell, atoms=atoms)
```

#### Writing to Files

```python
# Basic write
stru.write("STRU", fmt="stru")

# Specify coordinate type
stru.write("STRU", fmt="stru", direct=True)   # Direct coordinates
stru.write("STRU", fmt="stru", direct=False)  # Cartesian coordinates

# Empty atom handling
stru.write("STRU", fmt="stru", empty2x=True)  # empty → X
```

#### Export to Other Formats

```python
# ASE Atoms
ase_atoms = stru.to("ase")

# Pymatgen Structure
pmg_struct = stru.to("pymatgen")

# Phonopy Atoms
phonopy_atoms = stru.to("phonopy")
```

---

### 2. AbacusATOM - Atom Class

Represents all properties of a single atom.

#### Creating an Atom

```python
atom = AbacusATOM(
    label="Fe",                    # Atom label
    coord=(0.0, 0.0, 0.0),         # Coordinate (Å)
    element="Fe",                  # Element symbol (auto-inferred)
    mass=55.845,                   # Atomic mass (auto-obtained)
    pp="Fe.upf",                   # Pseudopotential file
    orb="Fe.orb",                  # Orbital file
    type_mag=2.0,                  # Atom type magnetic moment
    move=(True, True, True),       # Movement constraints
    mag=3.0,                       # Atomic magnetic moment (overrides type_mag)
    velocity=(0.0, 0.0, 0.0),      # Initial velocity
)
```

#### Accessing Properties

```python
print(atom.label)              # Atom label
print(atom.coord)              # Coordinate
print(atom.element)            # Element symbol
print(atom.mass)               # Atomic mass
print(atom.atommag)            # Actual magnetic moment
print(atom.noncolinear)        # Whether non-collinear magnetic
```

#### Setting Magnetic Moments

```python
# Collinear magnetism
atom.set_atommag(mag=3.0)

# Non-collinear magnetism (with angles)
atom.set_atommag(mag=2.0, angle1=90, angle2=45)

# Non-collinear magnetism (3D vector)
atom.set_atommag(mag=(1.0, 0.0, 0.0))
```

#### Setting Magnetic Constraints (DeltaSpin/CDFT)

```python
# For DeltaSpin or CDFT calculations, constrain atomic magnetic moments

# Collinear case (scalar constraint)
atom.constrain = True  # Constrain this atom's magnetic moment

# Non-collinear case (vector constraint)
atom.constrain = (True, True, True)   # Constrain all 3 components
atom.constrain = (True, False, True)  # Constrain x and z, free y

# Set constraint with lambda parameter (for DeltaSpin)
atom.constrain = True
atom.lambda_ = 0.1  # Lambda parameter for constraint strength

# Batch set constraints on structure
stru = AbacusSTRU.read("STRU")
for atom in stru.atoms:
    if atom.element == "Fe":
        atom.constrain = True  # Constrain all Fe atoms
```

**Important Notes for Constrain:**
- `constrain=True` is used in DeltaSpin/CDFT calculations to fix atomic magnetic moments
- For collinear calculations (`nspin=2`), `constrain` is a boolean
- For non-collinear calculations (`nspin=4`), `constrain` can be a tuple of 3 booleans for x/y/z components
- The `lambda_` parameter controls the constraint strength in DeltaSpin calculations
- Constrained atoms will have their magnetic moments fixed during SCF iterations

---

## Structure Operations

### Modifying the Cell

```python
# Directly set cell
stru.cell = [[5.0, 0.0, 0.0],
             [0.0, 5.0, 0.0],
             [0.0, 0.0, 5.0]]

# Rotate structure
import numpy as np
rot_mat = np.array([[0, 1, 0], 
                    [-1, 0, 0], 
                    [0, 0, 1]])  # 90° rotation around z-axis
stru.rotate(rot_mat)

# Permute lattice vectors
stru.permute_lat_vec(mode="bca")  # (a,b,c) → (b,c,a)
stru.permute_lat_vec(mode="cab")  # (a,b,c) → (c,a,b)
```

### Modifying Atomic Coordinates

```python
# Batch set coordinates (Cartesian, in Å)
new_coords = [(0.0, 0.0, 0.0), (2.0, 2.0, 2.0)]
stru.coords = new_coords

# Batch set coordinates (direct coordinates)
stru.coords_direct = [(0, 0, 0), (0.5, 0.5, 0.5)]

# Modify single atom
stru[0].coord = (1.0, 1.0, 1.0)

# Add atoms
new_atom = AbacusATOM(label="C", coord=(0.5, 0.5, 0.5))
stru.append(new_atom)
stru.insert(0, new_atom)  # Insert at specific position

# Delete atoms
del stru[0]

# Create subset
subset = stru.create_subset([0, 1, 2])  # Keep only atoms with indices 0,1,2
```

### Setting Magnetic Moments

```python
# Batch set magnetic moments (collinear)
stru.atom_mags = [1.0, -1.0, 0.5, 2.0]

# Batch set magnetic moments (non-collinear, 3D vectors)
stru.atom_mags = [
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
]

# Set single atom
stru.atoms[0].set_atommag(2.0, angle1=90, angle2=0)
```

### Fixing Atoms (for Structure Relaxation)

```python
# Fix by indices
stru.fix_atom_by_index([0, 1, 2], move=(False, False, False))

# Fix by coordinate range (fix atoms with z < 5Å)
stru.fix_atom_by_coord(
    min=0.0,
    max=5.0,
    cartesian=True,
    direction=2,  # z direction
    move=(False, False, False)
)

# Fix selected atoms, release others
stru.fix_atom_by_index([0, 1], move=(False, False, False), only=True)
```

### Setting Pseudopotential/Orbital Files

```python
# Batch set pseudopotentials
pp_dict = {"Si": "Si.upf", "C": "C.upf", "O": "O.upf"}
stru.set_pp(pp_dict, key_type="element")  # or key_type="label"

# Batch set orbitals
orb_dict = {"Si": "Si.orb", "C": "C.orb"}
stru.set_orb(orb_dict)

# Batch set PAW files
paw_dict = {"Si": "Si.paw"}
stru.set_paw(paw_dict)
```

### Atom Sorting

```python
# Sort atoms by type (same types grouped together)
indices = stru.sort(keep_first_order=True)
# indices returns original atom indices before rearrangement
```

### Build Supercell

```python
# Create supercell with replication factors
stru_super = stru.supercell([2, 2, 2])  # 2x2x2 supercell
stru_super = stru.supercell([3, 1, 2])  # Asymmetric 3x1x2 supercell

# The method automatically:
# - Scales the cell vectors
# - Replicates all atoms with proper coordinate shifts
# - Preserves all atom properties (mag, move, velocity, etc.)
# - Preserves metadata and dpks settings
```

### Getting High-Symmetry k-Point Path

```python
# Use SeekPath to get high-symmetry k-points and path
point_coords, path = stru.get_kline(
    with_time_reversal=True,
    recipe="hpkot",
    symprec=1e-5
)

# point_coords: {"G": [0,0,0], "X": [0.5,0,0], ...}
# path: [("G", "X"), ("X", "W"), ...]
```

---

## Common Use Cases

### Use Case 1: Batch Format Conversion

```python
from pathlib import Path

# Convert all POSCAR to STRU
for poscar in Path(".").glob("*/POSCAR"):
    stru = AbacusSTRU.read(str(poscar), fmt="poscar")
    stru.write(poscar.parent / "STRU", fmt="stru")
```

### Use Case 2: Fix Surface Bottom Layers

```python
# Read structure
stru = AbacusSTRU.read("STRU")

# Fix atoms with z < 5Å
stru.fix_atom_by_coord(min=0.0, max=5.0, cartesian=True, direction=2)

# Write new structure
stru.write("STRU_fixed", fmt="stru")
```

### Use Case 3: Set Antiferromagnetic System

```python
stru = AbacusSTRU.read("STRU")

# Antiferromagnetic setup (alternating positive/negative moments)
for i, atom in enumerate(stru.atoms):
    if i % 2 == 0:
        atom.set_atommag(3.0)   # Up
    else:
        atom.set_atommag(-3.0)  # Down

stru.write("STRU_AFM", fmt="stru")
```

### Use Case 4: Create ABACUS Input from CIF

```python
# Read CIF
stru = AbacusSTRU.read("structure.cif", fmt="cif")

# Set pseudopotentials
stru.set_pp({"Si": "Si.upf", "O": "O.upf"})

# Set orbitals
stru.set_orb({"Si": "Si.orb", "O": "O.orb"})

# Write STRU
stru.write("STRU", fmt="stru")
```

### Use Case 5: Structure Perturbation (for MD Initial Structure)

```python
import numpy as np
import copy

def perturb_structure(stru, amplitude=0.01):
    """Apply random perturbations to atomic positions"""
    new_stru = copy.deepcopy(stru)
    for atom in new_stru.atoms:
        displacement = np.random.uniform(-amplitude, amplitude, 3)
        atom.coord = tuple(c + d for c, d in zip(atom.coord, displacement))
    return new_stru

# Usage
stru = AbacusSTRU.read("STRU")
stru_perturbed = perturb_structure(stru, amplitude=0.02)
stru_perturbed.write("STRU_md", fmt="stru")
```

### Use Case 6: Build Supercell

```python
# Using the built-in supercell() method
stru = AbacusSTRU.read("STRU")

# Create 2x2x2 supercell
stru_super = stru.supercell([2, 2, 2])
stru_super.write("STRU_2x2x2", fmt="stru")

# Create asymmetric supercell (e.g., 3x1x2)
stru_super = stru.supercell([3, 1, 2])
stru_super.write("STRU_3x1x2", fmt="stru")
```

### Use Case 7: DeltaSpin/CDFT Calculation Setup

```python
# For DeltaSpin or Constrained DFT (CDFT) calculations
stru = AbacusSTRU.read("STRU")

# Collinear case: constrain magnetic moments of specific atoms
for atom in stru.atoms:
    if atom.element == "Fe":
        atom.set_atommag(mag=3.0)  # Set initial magnetic moment
        atom.constrain = True      # Constrain the moment

# Non-collinear case: constrain specific components
for atom in stru.atoms:
    if atom.element == "Co":
        atom.set_atommag(mag=2.0, angle1=90, angle2=0)
        atom.constrain = (True, True, False)  # Constrain x,y; free z

# Set lambda parameter for constraint strength (DeltaSpin)
for atom in stru.atoms:
    if atom.constrain:
        atom.lambda_ = 0.1  # Adjust lambda as needed

stru.write("STRU_deltaconstrain", fmt="stru")
```

**Note:** For DeltaSpin/CDFT calculations, also set appropriate parameters in INPUT file:
- `dft_plus_u 1` - Enable DFT+U (if needed)
- `sc_mag_switch 1` - Enable DeltaSpin constraint
- See `abacus-input-parameter` skill for detailed INPUT settings

---

## Unit Conventions

| Operation | Unit | Description |
|-----------|------|-------------|
| Coordinates/cell in Python API | Angstrom (Å) | Code uniformly uses Å |
| STRU file read/write | Bohr | Automatically converted |
| Direct coordinates | Dimensionless | Fractional coordinates (0-1) |

---

## Magnetic Moment Settings

| Parameter | Description | Priority |
|-----------|-------------|----------|
| `type_mag` | Atom type-level magnetic moment (defined once per type in STRU) | Low |
| `mag` | Single atom-level magnetic moment | High (overrides type_mag) |
| `angle1/angle2` | Angle settings for non-collinear magnetism | Used with mag |
| `constrain` | Constraint flag for DeltaSpin/CDFT calculations | N/A |
| `lambda_` | Lambda parameter for constraint strength (DeltaSpin) | N/A |

**Non-collinear magnetism**: When both `mag` (3 components) and `angle1/angle2` are set:
1. ABACUS first calculates the magnitude of `mag`: `|mag| = √(mx² + my² + mz²)`
2. Then uses `angle1/angle2` to define the direction
3. Final magnetic moment = `|mag|` along the direction specified by angles

**Magnetic constraints (DeltaSpin/CDFT)**:
- `constrain=True` (collinear) or `constrain=(True, True, True)` (non-collinear) fixes the magnetic moment during SCF
- Used in DeltaSpin and Constrained DFT (CDFT) calculations
- `lambda_` parameter controls the constraint strength (typical values: 0.01-1.0)
- For non-collinear calculations, `constrain` can be a tuple to constrain specific components (e.g., `(True, False, True)` constrains x and z, leaves y free)

---

## Common Errors and Fixes

### Error 1: Element Label Mismatch

```python
# ❌ Wrong: Atom labels in STRU don't match ATOMIC_SPECIES
# Causes pseudopotential setup failure

# ✅ Fix: Ensure label consistency
stru.set_pp({"Si": "Si.upf"}, key_type="label")
```

### Error 2: Atom Types Not Sorted

```python
# ❌ Wrong: Same type atoms are scattered
# ABACUS requires same-type atoms to be consecutive

# ✅ Fix: Sort before writing
stru.sort(keep_first_order=True)
stru.write("STRU", fmt="stru")
```

### Error 3: Magnetic Moment Count Mismatch

```python
# ❌ Wrong: Number of moments doesn't match atom count
stru.atom_mags = [1.0, -1.0]  # Only 2 moments for 4 atoms

# ✅ Fix: Ensure counts match
stru.atom_mags = [1.0, -1.0, 1.0, -1.0]  # 4 moments for 4 atoms
```

### Error 4: Pseudopotential Path Issues

```python
# ⚠️ Note: PP paths in STRU are relative to pseudo_dir in INPUT
# Example:
# INPUT: pseudo_dir = a/b/
# STRU:  Si 14 c/d/Si.upf
# ABACUS actually looks for: a/b/c/d/Si.upf

# ✅ Fix: Use absolute paths or correctly set pseudo_dir
stru.set_pp({"Si": "/absolute/path/Si.upf"}, key_type="label")
```

---

## Tips

| Tip | Description |
|-----|-------------|
| Use `Direct` coordinates for output | More readable and easier to modify |
| Match atom labels exactly | Case-sensitive matching |
| Set `move=(False,False,False)` for substrate atoms | Fix bottom layers for surface calculations |
| Use `Cartesian_angstrom` for MD | Easier to set velocities |
| Verify pseudopotential paths | Missing files cause ABACUS crash |
| Use `latname` to avoid manual cell setup | Reduces errors |
| Use `constrain=True` for DeltaSpin/CDFT | Fix magnetic moments of specific atoms |
| Set `lambda_` for constraint strength | Typical values: 0.01-1.0 |

---

## Related Skills

- **STRU file format**: `abacus-stru`
- **Batch input preparation**: `abacustest-prepare-inputs`
- **Model calculations**: `abacustest-models`
- **INPUT parameters**: `abacus-input-parameter`
- **Format conversion**: `abacustest-abacus2VaspQeCp2k`

---

## External Resources

- **ABACUS Documentation**: https://abacus.deepmodeling.com/
- **Pseudopotential Library**: http://abacus.ustc.edu.cn/pseudo/list.htm
- **abacustest GitHub**: https://github.com/DeepModeling/abacustest
