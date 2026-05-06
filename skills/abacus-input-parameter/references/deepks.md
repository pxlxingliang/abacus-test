# DeePKS (Deep Learning Kohn-Sham)

Parameters for DeePKS method - machine learning enhanced DFT calculations.

## DeePKS Output

### deepks_out_labels
- **Type**: Boolean
- **Default**: `False`
- **Description**: Print labels and descriptors for DeePKS training.
- **Output Files**: Start with "deepks" in `OUT.${suffix}/`
- **Note**: Requires numerical descriptor file (`jle.orb`) in STRU under `NUMERICAL_DESCRIPTOR` tag.

### deepks_scf
- **Type**: Boolean
- **Default**: `False`
- **Description**: Perform self-consistent field iteration in DeePKS.
- **Requirement**: Traced model file needed.

### deepks_equiv
- **Type**: Boolean
- **Default**: `False`
- **Description**: Use equivariant version of DeePKS.
- **Note**: Under development, internal use only.

### deepks_model
- **Type**: String
- **Default**: `None`
- **Description**: Path to trained, traced neural network model from deepks-kit.
- **Source**: https://github.com/deepmodeling/deepks-kit

## Bessel Projectors

### bessel_descriptor_lmax
- **Type**: Integer
- **Default**: `2`
- **Description**: Maximum angular momentum of Bessel projectors.
- **Generation**: Use `calculation = gen_bessel` to generate projectors.

### bessel_descriptor_ecut
- **Type**: Real
- **Default**: Same as ecutwfc
- **Description**: Energy cutoff of Bessel functions.

### bessel_descriptor_tolerence
- **Type**: Real
- **Default**: `1.0e-12`
- **Description**: Tolerance for finding Bessel function zeros.

### bessel_descriptor_rcut
- **Type**: Real
- **Default**: `6.0` Bohr
- **Description**: Cutoff radius of Bessel functions.

### bessel_descriptor_smooth
- **Type**: Boolean
- **Default**: `False`
- **Description**: Smooth Bessel functions at cutoff radius.

### bessel_descriptor_sigma
- **Type**: Real
- **Default**: `0.1` Bohr
- **Description**: Smoothing parameter at cutoff.

## DeePKS Labels

### deepks_bandgap
- **Type**: Boolean
- **Default**: `False`
- **Description**: Include bandgap label for training.

### deepks_v_delta
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - No V_delta output
  - `1` - Output h_base.npy, v_delta.npy, h_tot.npy, v_delta_precalc.npy
  - `2` - Output h_base.npy, v_delta.npy, h_tot.npy, phialpha.npy, grad_evdm.npy (no v_delta_precalc.npy)
- **Note**: Option 2 recommended for large systems (v_delta_precalc.npy becomes very large).

### deepks_out_unittest
- **Type**: Boolean
- **Default**: `False`
- **Description**: Generate files for DeePKS unit test.
- **Note**: Run with 1 process only.

---

## Quick Examples

### DeePKS Label Generation
```
deepks_out_labels 1
deepks_model /path/to/model.pt
```

### DeePKS-SCF Calculation
```
deepks_scf 1
deepks_model /path/to/traced_model.pt
```

### Generate Bessel Projectors
```
calculation gen_bessel
bessel_descriptor_lmax 2
bessel_descriptor_rcut 6.0
bessel_descriptor_ecut 60
```

### DeePKS with V_delta (Large System)
```
deepks_out_labels 1
deepks_v_delta 2
```

---

## STRU File Setup

For DeePKS, the STRU file needs:
```
NUMERICAL_ORBITAL
H_gga_8au_60Ry_2s1p.orb
O_gga_7au_60Ry_2s2p1d.orb
NUMERICAL_DESCRIPTOR
jle.orb
```

---

## Related References

- [System Variables](system-variables.md) - calculation parameter
- [NAOs](naos.md) - Numerical atomic orbital parameters
- [Output](output.md) - General output parameters
