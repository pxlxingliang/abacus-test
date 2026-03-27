# DFT+U Correction

Parameters for DFT+U calculations (Hubbard correction for correlated electrons).

## DFT+U Control

### dft_plus_u
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - No DFT+U correction
  - `1` - Radius-adjustable localized projections (with onsite_radius)
  - `2` - First zeta of NAOs as projections (old method, testing)

### orbital_corr
- **Type**: Integer
- **Default**: `-1` (per atom type)
- **Options** (per atom type):
  - `-1` - No +U correction
  - `1` - p-electron orbits
  - `2` - d-electron orbits
  - `3` - f-electron orbits
- **Format**: Specify for each atom type in order (e.g., `2 -1` for Fe-O system with U on Fe only)

### hubbard_u
- **Type**: Real
- **Default**: `0.0` eV
- **Description**: Hubbard U parameter (actually U_eff = U - J in Dudarev scheme).
- **Format**: Specify for each atom type (unless using yukawa_potential).

### yukawa_potential
- **Type**: Boolean
- **Default**: `False`
- **Description**: Use local screened Coulomb potential to calculate U and J.
- **Note**: When True, hubbard_u does not need to be specified.

### yukawa_lambda
- **Type**: Real
- **Default**: Calculated on-the-fly (average of system)
- **Description**: Screening length of Yukawa potential.
- **Recommendation**: Stick to default unless有特殊 reason.

## U-Ramping

### uramping
- **Type**: Real
- **Default**: `-1.0` (disabled)
- **Description**: U-ramping for difficult convergence.
- **Mechanism**:
  - Start SCF with U = 0 eV (normal LDA/PBE)
  - When SCF restarts (drho < mixing_restart), increase U by uramping eV
  - Repeat until U reaches target hubbard_u
- **Recommendation**: For uramping = 1.0 eV, use mixing_restart ≈ 5e-4

## Occupation Matrix Control

### omc
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - No control; calculate from wavefunctions each SCF step
  - `1` - Read initial density matrix from `initial_onsite.dm`, then update
  - `2` - Use same density matrix from `initial_onsite.dm` throughout
- **Note**: Create `initial_onsite.dm` from a previous DFT+U calculation's `onsite.dm` file.

### onsite_radius
- **Type**: Real
- **Default**: `3.0` Bohr
- **Description**: Modulate single-zeta portion of NAOs for DFT+U projections.
- **Algorithm**: Smooth truncation with normalization:
  - g(r;σ) = 1 - exp(-(r-r_c)²/(2σ²)) for r < r_c
  - Optimized to minimize error with derivative term

---

## Quick Examples

### Simple DFT+U on Fe d-orbitals
```
dft_plus_u 1
orbital_corr 2 -1
hubbard_u 4.0 0.0
```

### DFT+U with U-Ramping (difficult convergence)
```
dft_plus_u 1
orbital_corr 2
hubbard_u 5.0
uramping 0.5
mixing_restart 5e-4
```

### Yukawa Potential (auto U)
```
dft_plus_u 1
orbital_corr 2
yukawa_potential 1
```

### Fixed Occupation Matrix
```
dft_plus_u 1
orbital_corr 2
hubbard_u 4.0
omc 2
```

---

## Related References

- [Electronic Structure](electronic-structure.md) - mixing parameters, scf convergence
- [LCAO](lcao.md) - NAO-related parameters
- [Output](output.md) - Output density matrix (out_dm)
