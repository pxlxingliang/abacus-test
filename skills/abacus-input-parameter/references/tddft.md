# TDDFT (Time-Dependent Density Functional Theory)

Parameters for real-time TDDFT simulations.

## TDDFT Control

### td_edm
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - New method (original formula)
  - `1` - Old method (ground state formula)
- **Description**: Energy density matrix calculation method.

### td_propagator
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - Crank-Nicolson
  - `1` - 4th-order Taylor expansion
  - `2` - Enforced Time-Reversal Symmetry (ETRS)

## External Electric Field

### td_vext
- **Type**: Boolean
- **Default**: `False`
- **Description**: Add external laser field.

### td_vext_dire
- **Type**: String
- **Default**: `1`
- **Description**: Electric field direction(s).
- **Format**: Space-separated directions (1=x, 2=y, 3=z)
- **Example**: `1 2` = x and y directions simultaneously

### td_stype
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - Length gauge
  - `1` - Velocity gauge

### td_ttype
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - Gaussian
  - `1` - Trapezoid
  - `2` - Trigonometric
  - `3` - Heaviside
  - `4` - HHG

### td_tstart, td_tend
- **Type**: Integer
- **Default**: `1`, `100`
- **Description**: Electric field start/end steps.

## Length Gauge

### td_lcut1, td_lcut2
- **Type**: Real
- **Default**: `0.05`, `0.05`
- **Description**: Interval for length gauge field.
- **Field Profile**:
  - E = E0 for cut1 < x < cut2
  - E = -E0/(cut1+1-cut2) elsewhere

## Gaussian Field

### td_gauss_freq
- **Type**: Real
- **Default**: `22.13` fs⁻¹
- **Description**: Frequency of Gaussian field.

### td_gauss_phase
- **Type**: Real
- **Default**: `0.0`
- **Description**: Phase delay.

### td_gauss_sigma
- **Type**: Real
- **Default**: `30.0` fs
- **Description**: Pulse width.

### td_gauss_t0
- **Type**: Real
- **Default**: `100`
- **Description**: Time center (steps).

### td_gauss_amp
- **Type**: Real
- **Default**: `0.25` V/Å
- **Description**: Amplitude.

**Formula**: amp × cos(2π×freq×(t-t0)+phase) × exp(-(t-t0)²/(2σ²))

## Trapezoid Field

### td_trape_freq
- **Type**: Real
- **Default**: `1.60` fs⁻¹

### td_trape_phase
- **Type**: Real
- **Default**: `0.0`

### td_trape_t1, td_trape_t2, td_trape_t3
- **Type**: Real
- **Default**: `1875`, `5625`, `7500`
- **Description**: Time intervals.

### td_trape_amp
- **Type**: Real
- **Default**: `2.74` V/Å

**Formula**:
- E = amp×cos(2π×freq×t+phase)×t/t1 for t < t1
- E = amp×cos(2π×freq×t+phase) for t1 < t < t2
- E = amp×cos(2π×freq×t+phase)×(1-(t-t2)/(t3-t2)) for t2 < t < t3
- E = 0 for t > t3

## Trigonometric Field

### td_trigo_freq1, td_trigo_freq2
- **Type**: Real
- **Default**: `1.164656`, `0.029116` fs⁻¹

### td_trigo_phase1, td_trigo_phase2
- **Type**: Real
- **Default**: `0.0`, `0.0`

### td_trigo_amp
- **Type**: Real
- **Default**: `2.74` V/Å

**Formula**: amp × cos(2π×freq1×t+phase1) × sin(2π×freq2×t+phase2)²

## Heaviside Field

### td_heavi_t0
- **Type**: Real
- **Default**: `100`
- **Description**: Switch time.

### td_heavi_amp
- **Type**: Real
- **Default**: `2.74` V/Å

**Formula**: E = amp for t < t0, E = 0 for t > t0

## Output

### out_dipole
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output dipole moment.

### out_current
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output current (velocity gauge).

### out_current_k
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output current for all k-points.

### out_efield
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output E-field (V/Å).

### out_vecpot
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output vector potential (a.u.) to `At.dat`.

### init_vecpot_file
- **Type**: Boolean
- **Default**: `False`
- **Description**: Initialize vector potential from `At.dat` file.

## Fixed Occupation

### ocp
- **Type**: Boolean
- **Default**: `False`
- **Description**: Fix band occupations from ocp_set.

### ocp_set
- **Type**: String
- **Default**: `None`
- **Description**: Occupation string (space-separated, N*x shorthand supported).
- **Example**: `1 10*1 0 1` = 13 bands, 12th unoccupied

---

## Quick Examples

### Gaussian Pulse
```
td_vext 1
td_stype 0
td_ttype 0
td_gauss_freq 0.1
td_gauss_amp 0.1
td_gauss_sigma 50
td_gauss_t0 100
out_dipole 1
```

### Trapezoid Pulse
```
td_vext 1
td_ttype 1
td_trape_freq 0.05
td_trape_amp 0.05
td_trape_t1 500
td_trape_t2 1500
td_trape_t3 2000
```

### Multi-Direction Field
```
td_vext 1
td_vext_dire 1 2
td_gauss_amp 0.1
td_gauss_phase 0 1.57
```

---

## Related References

- [Electronic Structure](electronic-structure.md) - SCF parameters
- [LR-TDDFT](lr-tddft.md) - Linear response TDDFT
- [Output](output.md) - General output parameters
