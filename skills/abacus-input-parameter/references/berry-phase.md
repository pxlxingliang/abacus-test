# Berry Phase and Wannier90 Interface

Parameters for Berry phase (polarization) calculations and Wannier90 interface.

## Berry Phase

### berry_phase
- **Type**: Boolean
- **Default**: `False`
- **Description**: Enable Berry phase calculation for polarization.

### gdir
- **Type**: Integer
- **Default**: `3`
- **Options**:
  - `1` - Polarization along lattice vector a₁
  - `2` - Polarization along lattice vector a₂
  - `3` - Polarization along lattice vector a₃
- **Description**: Direction of polarization calculation.

## Wannier90 Interface

### towannier90
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - Do not generate Wannier90 files
  - `1` - Generate Wannier90 interface files

### wannier_card
- **Type**: String
- **Default**: `"none"`
- **Description**: Wannier90 input file name (e.g., `seedname.win`).

### nnkpfile
- **Type**: String
- **Default**: `seedname.nnkp`
- **Description**: File from `wannier90 -pp` command.

### wannier_method
- **Type**: Integer
- **Default**: `1`
- **Options**:
  - `1` - lcao_in_pw method (increase ecutwfc for accuracy)
  - `2` - Grid integration (Gauss-Legendre radial, Lebedev-Laikov spherical)
- **Note**: LCAO basis only.

### wannier_spin
- **Type**: String
- **Default**: `up`
- **Options**: `up`, `down`
- **Description**: Spin direction for nspin=2 calculations.

## Wannier90 Output Files

### out_wannier_mmn
- **Type**: Boolean (0/1)
- **Default**: `1`
- **Description**: Write `*.mmn` file (overlap matrices).

### out_wannier_amn
- **Type**: Boolean (0/1)
- **Default**: `1`
- **Description**: Write `*.amn` file (projection matrices).

### out_wannier_eig
- **Type**: Boolean (0/1)
- **Default**: `1`
- **Description**: Write `*.eig` file (eigenvalues).

### out_wannier_unk
- **Type**: Boolean (0/1)
- **Default**: `0`
- **Description**: Write `UNK.*` files (wavefunctions).

### out_wannier_wvfn_formatted
- **Type**: Boolean (0/1)
- **Default**: `1`
- **Description**: UNK file format:
  - `0` - Binary format
  - `1` - ASCII (text) format

---

## Quick Examples

### Berry Phase (Polarization)
```
berry_phase 1
gdir 3
```

### Wannier90 Interface (LCAO)
```
basis_type lcao
towannier90 1
wannier_card mysystem.win
nnkpfile mysystem.nnkp
out_wfc_lcao 2
```

### Wannier90 with Spin
```
nspin 2
towannier90 1
wannier_spin up
out_wannier_mmn 1
out_wannier_amn 1
```

### Full Wannier90 Workflow
```
basis_type lcao
towannier90 1
wannier_card seedname.win
wannier_method 1
out_wannier_mmn 1
out_wannier_amn 1
out_wannier_eig 1
out_wannier_unk 1
out_wannier_wvfn_formatted 1
```

---

## Related References

- [Electronic Structure](electronic-structure.md) - nspin, basis_type
- [LCAO](lcao.md) - LCAO-specific parameters
- [Output](output.md) - General output parameters
