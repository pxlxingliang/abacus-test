# NAOs (Numerical Atomic Orbitals)

Parameters for numerical atomic orbital generation and manipulation.

## Bessel NAOs

### bessel_nao_ecut
- **Type**: Real
- **Default**: `ecutwfc`
- **Unit**: Ry
- **Description**: Energy cutoff for spherical Bessel functions.
- **Formula**: Number of Bessel functions ≈ √(bessel_nao_ecut) × bessel_nao_rcut / π

### bessel_nao_tolerence
- **Type**: Real
- **Default**: `1.0e-12`
- **Description**: Tolerance for finding zeros of spherical Bessel functions.

### bessel_nao_rcut
- **Type**: Real
- **Default**: `6.0` Bohr
- **Description**: Cutoff radius and common node of spherical Bessel functions.

### bessel_nao_smooth
- **Type**: Boolean
- **Default**: `True`
- **Description**: Smooth NAOs near cutoff radius.
- **Formula**: 1 - exp(-(r - r_cut)² / (2σ²))

### bessel_nao_sigma
- **Type**: Real
- **Default**: `0.1` Bohr
- **Description**: Smoothing range for NAOs.

---

## Quick Examples

### Standard NAO Generation
```
bessel_nao_ecut 60
bessel_nao_rcut 6.0
bessel_nao_smooth 1
bessel_nao_sigma 0.1
```

### High-Precision NAOs
```
bessel_nao_ecut 100
bessel_nao_rcut 8.0
bessel_nao_tolerence 1e-14
```

### Generate Projectors (DeePKS)
```
calculation gen_bessel
bessel_nao_ecut 60
bessel_nao_lmax 2
```

---

## Related References

- [DeePKS](deepks.md) - Bessel projectors for DeePKS
- [LCAO](lcao.md) - LCAO basis using NAOs
- [Plane Wave](plane-wave.md) - ecutwfc reference
