# Debugging Parameters

Parameters useful for debugging and testing ABACUS calculations.

## Hamiltonian Components

### t_in_h
- **Type**: Boolean (0/1)
- **Default**: `1`
- **Description**: Include kinetic energy term in Hamiltonian.

### vl_in_h
- **Type**: Boolean (0/1)
- **Default**: `1`
- **Description**: Include local pseudopotential term in Hamiltonian.

### vnl_in_h
- **Type**: Boolean (0/1)
- **Default**: `1`
- **Description**: Include non-local pseudopotential term in Hamiltonian.

### vh_in_h
- **Type**: Boolean (0/1)
- **Default**: `1`
- **Description**: Include Hartree potential term in Hamiltonian.

### vion_in_h
- **Type**: Boolean (0/1)
- **Default**: `1`
- **Description**: Include local ionic potential term in Hamiltonian.

## Force & Stress Debugging

### test_force
- **Type**: Boolean (0/1)
- **Default**: `0`
- **Description**: Output detailed force components.

### test_stress
- **Type**: Boolean (0/1)
- **Default**: `0`
- **Description**: Enable colorful terminal output for stress.

## Energy Debugging

### test_skip_ewald
- **Type**: Boolean (0/1)
- **Default**: `0`
- **Description**: Skip Ewald energy calculation.

---

## Quick Examples

### Debug Hamiltonian Construction
```
vl_in_h 0
vnl_in_h 1
vh_in_h 1
```

### Detailed Force Analysis
```
test_force 1
cal_force 1
```

### Skip Ewald (Testing)
```
test_skip_ewald 1
```

---

## Related References

- [Output](output.md) - Output control parameters
- [Electronic Structure](electronic-structure.md) - SCF parameters
