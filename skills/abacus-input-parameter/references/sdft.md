# SDFT (Stochastic Density Functional Theory)

Parameters for stochastic DFT calculations - efficient for large systems.

## SDFT Method

### method_sto
- **Type**: Integer
- **Default**: `2`
- **Options**:
  - `1` - Calculate T_n(ĥ)|χ⟩ twice (less memory, slower)
  - `2` - Calculate T_n(ĥ)|χ⟩ once (more memory, faster)
  - Other - Use method 2
- **Description**: Stochastic DFT algorithm.

### nbands_sto
- **Type**: Integer or String
- **Default**: `256`
- **Options**:
  - `> 0` - Number of stochastic orbitals (perform SDFT)
  - `0` - Perform Kohn-Sham DFT
  - `all` - All complete basis (Chebyshev method, same as KSDFT without stochastic error)
- **Note**: For mixed stochastic-deterministic DFT, also set `nbands` for KS orbitals.
- **Accuracy**: Stochastic error scales as 1/√N_χ

### nche_sto
- **Type**: Integer
- **Default**: `100`
- **Description**: Chebyshev expansion order.

## Energy Bounds

### emin_sto
- **Type**: Real
- **Default**: `0.0`
- **Description**: Trial lower bound of Hamiltonian eigenvalues.

### emax_sto
- **Type**: Real
- **Default**: `0.0`
- **Description**: Trial upper bound of Hamiltonian eigenvalues.

## Random Seed

### seed_sto
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `≥ 0` - Stochastic orbitals: exp(i2πθ(G)), θ ~ Uniform(0,1)
  - `≤ -1` - Stochastic orbitals: ±1 with equal probability
  - `0` - Seed from time(NULL)

## Initialization

### initsto_ecut
- **Type**: Real
- **Default**: `0.0`
- **Description**: Cutoff for stochastic wavefunction initialization box.
- **Requirement**: Must be > ecutwfc (otherwise disabled).
- **Benefit**: Results independent of core count; consistent G coefficients.

### initsto_freq
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `> 0` - Update stochastic orbitals every initsto_freq steps (MD)
  - `0` - Never change stochastic orbitals

## Post-Processing

### npart_sto
- **Type**: Integer
- **Default**: `1`
- **Description**: Memory reduction factor for SDFT post-processing (DOS, conductivities).
- **Effect**: Memory cost = 1/npart_sto × original.

---

## Quick Examples

### Standard SDFT
```
esolver_type sdft
nbands_sto 256
nche_sto 100
method_sto 2
```

### Mixed Stochastic-Deterministic
```
esolver_type sdft
nbands 50
nbands_sto 200
```

### High Accuracy (KSDFT-equivalent)
```
esolver_type sdft
nbands_sto all
nche_sto 200
```

### SDFT-MD with Orbital Updates
```
esolver_type sdft
calculation md
initsto_freq 10
```

### Large System (Memory Efficient)
```
nbands_sto 512
npart_sto 4
method_sto 1
```

---

## Related References

- [Electronic Structure](electronic-structure.md) - nbands, esolver_type
- [Conductivities](conductivities.md) - SDFT-based conductivity calculations
- [DOS](dos.md) - Density of states from SDFT
