# van der Waals Correction

Parameters for dispersion correction (DFT-D2, DFT-D3).

## vdW Method Selection

### vdw_method
- **Type**: String
- **Default**: `none`
- **Options**:
  - `none` - No vdW correction
  - `d2` - Grimme DFT-D2
  - `d3_0` - Grimme DFT-D3(0) (zero-damping)
  - `d3_bj` - Grimme DFT-D3(BJ) (BJ-damping)
- **Auto-Setting**: ABACUS ≥ 3.8.3 auto-sets D3 parameters for common functionals.
- **Note**: Must specify dft_functional explicitly for auto-setting to work.

## DFT-D2 Parameters

### vdw_s6
- **Type**: Real
- **Default**: `0.75` (PBE for d2)
- **Description**: Scale factor for dispersion energy.
- **Recommended Values** (d2):
  - PBE: 0.75
  - BLYP: 1.2
  - B-P86: 1.05
  - TPSS: 1.0
  - B3LYP: 1.05

### vdw_d
- **Type**: Real
- **Default**: `20`
- **Description**: Damping rate for DFT-D2.

### vdw_C6_file
- **Type**: String
- **Default**: `default` (built-in)
- **Description**: Custom C6 parameter file.
- **Format**: Element C6 (one per line)
```
H  0.1
Si 9.0
```

### vdw_C6_unit
- **Type**: String
- **Default**: `Jnm6/mol`
- **Options**: `Jnm6/mol`, `eVA` (eV·Å⁶)

### vdw_R0_file
- **Type**: String
- **Default**: `default` (built-in)
- **Description**: Custom R0 parameter file.
- **Format**: Element R0 (one per line)

### vdw_R0_unit
- **Type**: String
- **Default**: `A` (Ångström)
- **Options**: `A`, `Bohr`

## DFT-D3 Parameters

### vdw_s8
- **Type**: Real
- **Default**: Auto from dft_functional
- **Description**: Scale factor for D3(0) and D3(BJ).
- **Reference**: https://github.com/dftd3/simple-dftd3/blob/main/assets/parameters.toml

### vdw_a1
- **Type**: Real
- **Default**: Auto from dft_functional
- **Description**: Damping function parameter for D3.

### vdw_a2
- **Type**: Real
- **Default**: Auto from dft_functional
- **Description**: Damping function parameter for D3(BJ).

### vdw_abc
- **Type**: Boolean
- **Default**: `False`
- **Description**: Include three-body terms for DFT-D3.

## Cutoff Settings

### vdw_cutoff_type
- **Type**: String
- **Default**: `radius`
- **Options**:
  - `radius` - Spherical cutoff (vdw_cutoff_radius)
  - `period` - Supercell cutoff (vdw_cutoff_period)

### vdw_cutoff_radius
- **Type**: Real
- **Default**: `56.6918` (d2), `95` (d3_0/d3_bj) Bohr
- **Description**: Cutoff sphere radius.

### vdw_radius_unit
- **Type**: String
- **Default**: `Bohr`
- **Options**: `A`, `Bohr`

### vdw_cutoff_period
- **Type**: Integer Integer Integer
- **Default**: `3 3 3`
- **Description**: Supercell extent in lattice vector directions.

### vdw_cn_thr
- **Type**: Real
- **Default**: `40` Bohr
- **Description**: Coordination number cutoff radius.

### vdw_cn_thr_unit
- **Type**: String
- **Default**: `Bohr`
- **Options**: `A`, `Bohr`

---

## Quick Examples

### DFT-D3(BJ) with PBE
```
dft_functional pbe
vdw_method d3_bj
```

### DFT-D2 with Custom Parameters
```
vdw_method d2
vdw_s6 1.05
vdw_d 20
```

### DFT-D3 with Three-Body Terms
```
dft_functional b3lyp
vdw_method d3_bj
vdw_abc 1
```

### Custom C6 Parameters
```
vdw_method d2
vdw_C6_file my_c6.txt
vdw_C6_unit eVA
```

### Large Cutoff for Accuracy
```
vdw_method d3_bj
vdw_cutoff_type radius
vdw_cutoff_radius 120
vdw_radius_unit bohr
```

---

## Special Cases

### wB97X-D3BJ
```
dft_functional HYB_GGA_WB97X_V
vdw_method d3_bj
```

### wB97X-D3(0)
```
dft_functional HYB_GGA_WB97X_D3
vdw_method d3_0
```

---

## Related References

- [Electronic Structure](electronic-structure.md) - dft_functional
- [System Variables](system-variables.md) - General settings
