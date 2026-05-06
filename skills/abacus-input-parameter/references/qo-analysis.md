# Quasiatomic Orbital (QO) Analysis

Parameters for quasiatomic orbital analysis.

## QO Control

### qo_switch
- **Type**: Boolean
- **Default**: `0`
- **Description**: Enable QO analysis output.

### qo_basis
- **Type**: String
- **Default**: `szv`
- **Options**:
  - `pswfc` - Pseudowavefunctions from pseudopotential files
  - `hydrogen` - Hydrogen-like atomic basis (with Slater screening)
  - `szv` - First zeta of NAOs (recommended for ABACUS-LCAO)
- **Notes**:
  - `pswfc`: Requires norm-conserving PP with pseudowavefunctions (SG15 not supported)
  - `szv`: Recommended for ABACUS-LCAO calculations

### qo_strategy
- **Type**: String [String...] (optional)
- **Default**: `energy-valence` (hydrogen), `all` (pswfc/szv)
- **Options** (for hydrogen basis):
  - `minimal-nodeless` - Only nodeless orbitals up to highest occupied
  - `minimal-valence` - Only highest n orbitals
  - `full` - All possible orbitals up to max n
  - `energy-full` - Aufbau principle orbitals
  - `energy-valence` - Highest n and n-1 layer orbitals
- **Options** (for pswfc/szv):
  - `all` - All orbitals from file
  - `s`, `p`, `d`, `f` - Specific angular momentum
  - `spd` - s, p, d orbitals (any combination)
- **Note**: One parameter applies to all atom types; multiple parameters can be given for different types.

### qo_screening_coeff
- **Type**: Real [Real...] (optional)
- **Default**: `0.1`
- **Description**: Rescale radial orbital shape.
- **For pswfc**: Screening factor e^(-η|r|) multiplied to mimic electron behavior.
- **For hydrogen**: Apply Slater screening (Z_eff = Z - σ).
- **Note**: One value applies to all atom types; multiple values for different types.

### qo_thr
- **Type**: Real
- **Default**: `1.0e-6`
- **Description**: Convergence threshold for orbital generation cutoff.
- **Effect**: Lower threshold = larger cutoff radius.

---

## Quick Examples

### Standard QO Analysis (LCAO)
```
qo_switch 1
qo_basis szv
qo_strategy all
```

### QO with Pseudowavefunctions
```
qo_switch 1
qo_basis pswfc
qo_strategy spd
qo_thr 1e-7
```

### Hydrogen-like Basis with Slater Screening
```
qo_switch 1
qo_basis hydrogen
qo_strategy energy-valence
qo_screening_coeff 0.5
```

### Custom Strategy per Element
```
qo_switch 1
qo_basis hydrogen
qo_strategy energy-valence minimal-valence
qo_screening_coeff 0.3 0.5
```

---

## Related References

- [LCAO](lcao.md) - Numerical atomic orbitals
- [Output](output.md) - Output control
