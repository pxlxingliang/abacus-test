# Output Parameters

Control what information and files ABACUS outputs during calculations.

## Output Level & Logging

### out_level
- **Type**: String
- **Default**: `ie`
- **Options**:
  - `ie` - Electronic iteration level (useful SCF info)
  - `i` - Geometry relaxation level (adds relaxation info)
  - `m` - Molecular dynamics level (simplified output)

### out_alllog
- **Type**: Boolean
- **Default**: `False`
- **Description**: Write individual logs from all MPI ranks.
  - `True`: `running_${calculation}_${rank+1}.log`
  - `False`: Only `running_${calculation}.log` from rank 0

### printe
- **Type**: Integer
- **Default**: `scf_nmax`
- **Description**: Print band energies every printe steps.

## Charge Density Output

### out_chg
- **Type**: Integer [Integer] (optional)
- **Default**: `0 3`
- **Options**:
  - First integer:
    - `0` - No output
    - `1` - Output charge density (cube files)
    - `2` - Also output initial charge density
    - `-1` - Disable auto-backup restart file
  - Second integer: Precision (default 3, recommend 10 for high precision)
- **Output Files**:
  - nspin=1: `SPIN1_CHG.cube`
  - nspin=2: `SPIN1_CHG.cube`, `SPIN2_CHG.cube`
  - nspin=4: `SPIN1_CHG.cube` to `SPIN4_CHG.cube`

### out_pot
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - No output
  - `1` - Total local potential
  - `2` - Electrostatic potential only
  - `3` - Also output initial potential

### out_dm
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output density matrix (LCAO).

### out_dm1
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output density matrix in sparse format.

## Wavefunction Output

### out_wfc_pw
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - No output
  - `1` - Text format (`WAVEFUNC${K}.txt`)
  - `2` - Binary format (`WAVEFUNC${K}.dat`)

### out_wfc_r
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output real-space wavefunctions.

### out_wfc_lcao
- **Type**: Integer
- **Default**: `False`
- **Options**:
  - `0` - No output
  - `1` - Text format
  - `2` - Binary format
- **Files**: `WFC_NAO_GAMMA{K}_ION{step}.dat` or `WFC_NAO_K{K}_ION{step}.dat`

## Electronic Properties

### out_dos
- **Type**: Integer
- **Default**: `0`
- **Options**:
  - `0` - No output
  - `1` - DOS only
  - `2` - (reserved)
  - `lcao-only` - DOS + PDOS

### out_band
- **Type**: Boolean [Integer] (optional)
- **Default**: `False`
- **Description**: Output band structure (eV). Second parameter = precision (default 8).

### out_proj_band
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output projected band structure.

### out_bandgap
- **Type**: Boolean
- **Default**: `False`
- **Description**: Print bandgap per SCF iteration to log.

### out_elf
- **Type**: Integer [Integer] (optional)
- **Default**: `0 3`
- **Description**: Output electron localization function.
  - First: 0=none, 1=output
  - Second: Precision (default 3)

## Structure Output

### out_stru
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output STRU files per ionic step (`STRU_ION${istep}_D`).

### out_mul
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output Mulliken population analysis.

## Matrix Output

### out_mat_hs
- **Type**: Boolean [Integer] (optional)
- **Default**: `False 8`
- **Description**: Output Hamiltonian and overlap matrices (upper triangular).

### out_mat_hs2
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output H(R) and S(R) matrices in sparse format.

### out_mat_tk
- **Type**: Boolean [Integer] (optional)
- **Default**: `False [8]`
- **Description**: Output kinetic matrices T(k).

### out_mat_r
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output position matrix representation.

### out_mat_t
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output kinetic energy matrix T(R) (LCAO).

### out_mat_dh
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output Hamiltonian derivatives dH(R).

### out_mat_xc
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output XC matrices in KS orbital representation (for GW).

### out_eband_terms
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output band energy terms separately.

### out_hr_npz / out_dm_npz
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output H(R)/DM(R) in npz format (internal use).

### dm_to_rho
- **Type**: Boolean
- **Default**: `False`
- **Description**: Read DM(R) npz and create electron density (serial only).

## Restart & Backup

### restart_save
- **Type**: Boolean
- **Default**: `False`
- **Description**: Save charge density per ionic step for restart.
  - `auto`: Save to `OUT.${suffix}/restart/`
  - Other: Save to `${read_file_dir}/restart/`

### restart_load
- **Type**: Boolean
- **Default**: `False`
- **Description**: Load charge density from previous calculation for restart.

### rpa
- **Type**: Boolean
- **Default**: `False`
- **Description**: Generate output files for RPA calculations.

## Frequency & Format Control

### out_freq_elec
- **Type**: Integer
- **Default**: `scf_nmax`
- **Description**: Output frequency for charge density and wavefunctions (electronic iterations).

### out_interval
- **Type**: Integer
- **Default**: `1`
- **Description**: Output interval for MD calculations.

### out_app_flag
- **Type**: Boolean
- **Default**: `True`
- **Description**: Append mode for matrix outputs during MD.

### out_ndigits
- **Type**: Integer
- **Default**: `8`
- **Description**: Decimal precision for output data.

### out_element_info
- **Type**: Boolean
- **Default**: `False`
- **Description**: Output element information (pseudopotential/orbital info).

## Special Outputs

### nbands_istate
- **Type**: Integer
- **Default**: `5`
- **Description**: Number of bands around Fermi level for get_wf/get_pchg.

### bands_to_print
- **Type**: String
- **Default**: `none`
- **Description**: Flexible band selection string (e.g., `1 4*0 5*1 0`).

### if_separate_k
- **Type**: Boolean
- **Default**: `False`
- **Description**: Write partial charge densities for each k-point separately.

### wannier_card
- **Type**: String
- **Default**: `"none"`
- **Description**: Wannier90 input file name.

### towannier90
- **Type**: Integer
- **Default**: `0`
- **Description**: Generate Wannier90 interface files.

---

## Quick Examples

### Standard SCF Output
```
out_level ie
out_chg 0
out_dos 0
```

### Full Analysis Output
```
out_chg 1
out_dos 2
out_band 1
out_mul 1
out_stru 1
```

### Restart-Capable Calculation
```
restart_save 1
restart_load 1
read_file_dir ./previous_calc/OUT.ABACUS
```

### MD Trajectory Output
```
out_interval 10
dump_force 1
dump_vel 1
dump_virial 1
```

### Wannier90 Interface
```
towannier90 1
wannier_card seedname.win
out_wfc_lcao 2
```

---

## Related References

- [System Variables](system-variables.md) - suffix parameter
- [Electronic Structure](electronic-structure.md) - SCF parameters
- [DOS](dos.md) - DOS-specific output settings
