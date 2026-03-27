# Electronic Conductivities

Parameters for calculating electronic conductivities.

## Conductivity Control

### cal_cond
- **Type**: Boolean
- **Default**: `False`
- **Description**: Enable electronic conductivity calculation.

## Chebyshev Expansion

### cond_che_thr
- **Type**: Real
- **Default**: `1e-8`
- **Description**: Error threshold for Chebyshev expansion.

### cond_dtbatch
- **Type**: Integer
- **Default**: `0` (auto)
- **Description**: Expansion batch size for exp(iH×dt×cond_dtbatch).
- **Note**: `0` = auto-set for > 100 expansion orders. Faster but more memory.

## Frequency Settings

### cond_dw
- **Type**: Real
- **Default**: `0.1`
- **Unit**: Energy
- **Description**: Frequency interval (dω) for frequency-dependent conductivities.

### cond_wcut
- **Type**: Real
- **Default**: `10.0`
- **Unit**: Energy
- **Description**: Cutoff frequency for conductivities.

## Time Integration

### cond_dt
- **Type**: Real
- **Default**: `0.02`
- **Unit**: Time
- **Description**: Time interval (dt) for Onsager coefficient integration.

## Smearing

### cond_smear
- **Type**: Integer
- **Default**: `1`
- **Options**:
  - `1` - Gaussian smearing
  - `2` - Lorentzian smearing

### cond_fwhm
- **Type**: Real
- **Default**: `0.4`
- **Description**: Full width at half maximum.
- **Formula**:
  - Gaussian: FWHM = 2√(2 ln 2) × σ
  - Lorentzian: FWHM = 2γ

## Velocity Matrix

### cond_nonlocal
- **Type**: Boolean
- **Default**: `True`
- **Description**: Include nonlocal potential correction in velocity matrix.
- **Formula**:
  - `True`: m v̂ = p̂ + (im/ℏ)[V̂_NL, r̂]
  - `False`: m v̂ ≈ p̂

---

## Quick Examples

### Standard Conductivity
```
cal_cond 1
cond_dw 0.1
cond_wcut 10.0
cond_smear 1
cond_fwhm 0.4
```

### High-Precision Conductivity
```
cal_cond 1
cond_che_thr 1e-10
cond_dt 0.01
cond_dw 0.05
```

### Lorentzian Smearing
```
cal_cond 1
cond_smear 2
cond_fwhm 0.5
```

---

## Related References

- [SDFT](sdft.md) - Stochastic DFT for large systems
- [Electronic Structure](electronic-structure.md) - General SCF settings
- [Output](output.md) - Output control
