# Electric Field and Dipole Correction

Parameters for external electric field and dipole correction in slab calculations.

## Electric Field

### efield_flag
- **Type**: Boolean
- **Default**: `False`
- **Description**: Add saw-like potential simulating electric field.
- **Note**: When enabled, symmetry is auto-set to 0.

### efield_dir
- **Type**: Integer
- **Default**: `2`
- **Options**:
  - `0` - Parallel to reciprocal lattice vector b₁
  - `1` - Parallel to reciprocal lattice vector b₂
  - `2` - Parallel to reciprocal lattice vector b₃
- **Description**: Direction of electric field.

### efield_pos_max
- **Type**: Real
- **Default**: Auto (center of vacuum - width/20)
- **Range**: [0, 1)
- **Description**: Position of maximum saw-like potential (crystal axis).

### efield_pos_dec
- **Type**: Real
- **Default**: Auto (width of vacuum / 10)
- **Range**: (0, 1)
- **Description**: Zone where potential decreases.

### efield_amp
- **Type**: Real
- **Default**: `0.0`
- **Description**: Electric field amplitude.
- **Note**: Potential increases with slope efield_amp, then decreases. Discontinuity must be in empty region.

## Dipole Correction

### dip_cor_flag
- **Type**: Boolean
- **Default**: `False`
- **Description**: Add dipole correction to ionic potential.
- **Usage**: Slab geometry for surface calculations only.
- **Note**: Discontinuity must fall in empty space.

---

## Quick Examples

### Electric Field Along z
```
efield_flag 1
efield_dir 2
efield_amp 0.01
```

### Dipole Correction for Slab
```
dip_cor_flag 1
efield_amp 0
efield_dir 2
```

### Combined Field + Dipole
```
efield_flag 1
dip_cor_flag 1
efield_dir 2
efield_amp 0.005
symmetry 0
```

---

## Gate Field (Compensating Charge)

### gate_flag
- **Type**: Boolean
- **Default**: `False`
- **Description**: Add compensating charge via charged plate for charged cells.

### zgate
- **Type**: Real
- **Default**: `0.5`
- **Range**: [0, 1)
- **Description**: Position of charged plate (crystal axis).

### block
- **Type**: Boolean
- **Default**: `False`
- **Description**: Add potential barrier to prevent electron spillover.

### block_down
- **Type**: Real
- **Default**: `0.45`
- **Range**: [0, block_up)
- **Description**: Lower bound of potential barrier.

### block_up
- **Type**: Real
- **Default**: `0.55`
- **Range**: (block_down, 1)
- **Description**: Upper bound of potential barrier.

### block_height
- **Type**: Real
- **Default**: `0.1`
- **Description**: Height of potential barrier.

---

## Related References

- [System Variables](system-variables.md) - symmetry parameter
- [Output](output.md) - Output potential (out_pot)
