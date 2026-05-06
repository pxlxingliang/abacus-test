# Implicit Solvation Model

Parameters for implicit solvation calculations (continuum solvation).

## Solvation Control

### imp_sol
- **Type**: Boolean
- **Default**: `False`
- **Description**: Enable implicit solvation correction.

## Solvent Properties

### eb_k
- **Type**: Real
- **Default**: `80`
- **Description**: Relative permittivity (dielectric constant) of bulk solvent.
- **Example**: 80 for water.

### tau
- **Type**: Real
- **Default**: `1.0798e-05`
- **Description**: Effective surface tension parameter.
- **Includes**: Cavitation, dispersion, and repulsion interactions not captured by electrostatics.

### sigma_k
- **Type**: Real
- **Default**: `0.6`
- **Description**: Width of diffuse cavity determined by electronic structure.

### nc_k
- **Type**: Real
- **Default**: `0.00037`
- **Description**: Electron density threshold for dielectric cavity formation.

---

## Quick Examples

### Water Solvation
```
imp_sol 1
eb_k 80
```

### Custom Solvent
```
imp_sol 1
eb_k 37
tau 1.5e-5
sigma_k 0.5
```

---

## Related References

- [Electronic Structure](electronic-structure.md) - SCF settings
- [Output](output.md) - Output control
