# RDMFT (Reduced Density Matrix Functional Theory)

Parameters for reduced density matrix functional theory calculations.

## RDMFT Control

### rdmft
- **Type**: Boolean
- **Default**: `False`
- **Description**: Enable RDMFT calculation.

## Power Functional

### rdmft_power_alpha
- **Type**: Real
- **Default**: `0.656`
- **Description**: Alpha parameter for power functional.
- **Formula**: g(occupation) = occupation^α
- **Note**: Also applies to other EXX-type/hybrid functionals in RDMFT.

---

## Quick Examples

### Standard RDMFT (Power Functional)
```
rdmft 1
rdmft_power_alpha 0.656
```

### Custom Alpha
```
rdmft 1
rdmft_power_alpha 0.7
```

---

## Related References

- [Electronic Structure](electronic-structure.md) - General SCF settings
- [Exact Exchange](exact-exchange.md) - Hybrid functional parameters
