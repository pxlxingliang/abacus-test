# Density of States (DOS)

Parameters for density of states calculations.

## DOS Output

### out_dos
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - No DOS output
  - `1` - DOS only
  - `2` - (reserved)
  - `lcao-only` - DOS + PDOS (projected DOS)
- **Note**: See [Output](output.md) for details.

## DOS Parameters

### dos_edelta_ev
- **Type**: Real
- **Default**: `0.01` eV
- **Description**: Energy step size for DOS output.

### dos_sigma
- **Type**: Real
- **Default**: `0.07` eV
- **Description**: Gaussian smearing width for DOS.

### dos_scale
- **Type**: Real
- **Default**: `0.01`
- **Description**: Energy range expansion factor.
- **Formula**: Range = (emax - emin) × (1 + dos_scale), centered at (emax + emin)/2
- **Note**: Ignored if dos_emin_ev or dos_emax_ev is set.

### dos_emin_ev
- **Type**: Real
- **Default**: Minimal eigenenergy of Ĥ
- **Description**: Minimum energy for DOS output.
- **Note**: If set, dos_scale is ignored.

### dos_emax_ev
- **Type**: Real
- **Default**: Maximal eigenenergy of Ĥ
- **Description**: Maximum energy for DOS output.
- **Note**: If set, dos_scale is ignored.

### dos_nche
- **Type**: Integer
- **Default**: `100`
- **Description**: Chebyshev expansion order for DOS.

---

## Quick Examples

### Standard DOS
```
out_dos 1
dos_edelta_ev 0.01
dos_sigma 0.05
```

### PDOS (LCAO)
```
basis_type lcao
out_dos lcao-only
dos_sigma 0.07
```

### Custom Energy Range
```
out_dos 1
dos_emin_ev -10
dos_emax_ev 10
dos_edelta_ev 0.005
```

### High-Resolution DOS
```
out_dos 1
dos_edelta_ev 0.001
dos_sigma 0.02
dos_nche 200
```

---

## Related References

- [Output](output.md) - out_dos parameter
- [LCAO](lcao.md) - PDOS requires LCAO basis
- [SDFT](sdft.md) - Stochastic DOS calculations
