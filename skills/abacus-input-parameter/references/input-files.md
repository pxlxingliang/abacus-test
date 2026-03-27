# Input Files

Parameters related to input file paths and restart settings.

## File Paths

### stru_file
- **Type**: String
- **Default**: `STRU`
- **Description**: Structure file name.
- **Contents**: Atomic species, pseudopotentials, orbitals, cell, positions.
- **Note**: Ignored when calculation=md and md_restart=true.

### kpoint_file
- **Type**: String
- **Default**: `KPT`
- **Description**: K-points file name.
- **Notes**:
  - Not needed for gamma_only=1 (auto-generated)
  - Required for multiple k-points

### pseudo_dir
- **Type**: String
- **Default**: `""`
- **Description**: Pseudopotential directory.
- **Usage**: Combined with STRU filenames to form full paths.
- **Example**: `pseudo_dir "../"` + `Si.upf` → `../Si.upf`

### orbital_dir
- **Type**: String
- **Default**: `""`
- **Description**: Orbital file directory.
- **Usage**: Combined with STRU filenames to form full paths.

### read_file_dir
- **Type**: String
- **Default**: `OUT.${suffix}`
- **Description**: Directory for reading initial files (charge density, etc.).
- **Example**: `'./'` = working directory.

## Restart

### restart_load
- **Type**: Boolean
- **Default**: `False`
- **Description**: Load charge density from previous calculation.
- **Requirements**:
  - restart_save=true in previous run
  - Correct read_file_dir
  - Charge density file exists
- **Note**: For EXX calculations, Hexx(R) files are also read.

### restart_save
- **Type**: Boolean
- **Default**: `False`
- **Description**: Save charge density per ionic step for restart.
- **Output**:
  - `auto`: `OUT.${suffix}/restart/`
  - Other: `${read_file_dir}/restart/`
- **Note**: For EXX, Hexx(R) files also saved.

## Wannier90

### wannier_card
- **Type**: String
- **Default**: `"none"`
- **Description**: Wannier90 input file name.

---

## Quick Examples

### Standard Setup
```
stru_file STRU
kpoint_file KPT
pseudo_dir ./pp/
orbital_dir ./orb/
```

### Restart from Previous Calculation
```
read_file_dir ./previous_calc/OUT.ABACUS
restart_load 1
```

### Save for Future Restart
```
restart_save 1
read_file_dir ./
```

### Custom File Names
```
stru_file my_structure.stru
kpoint_file my_kpoints.kpt
```

---

## Related References

- [System Variables](system-variables.md) - suffix, calculation
- [Output](output.md) - Restart output settings
- [Electronic Structure](electronic-structure.md) - init_chg parameter
