# Band Structure Calculation Post-processing Guide

## Usage Command

```bash
abacustest model band post -j ./
```

## Output Results

### 1. band.png

Band structure plot showing electronic band dispersion relations along a specific k-path.

### 2. metrics_band.json

Contains key metrics from band structure analysis, with the following main content:

```json
{
  ".": {
    "band_gap": Band gap value (eV),
    "cbm": {
      "band_index": Conduction band minimum band indices,
      "kpoint_index": Conduction band minimum k-point indices,
      "kpoint_labels": Conduction band minimum k-point labels,
      "kpoint_coord": Conduction band minimum k-point coordinates,
      "energy": Conduction band minimum energy (eV)
    },
    "vbm": {
      "band_index": Valence band maximum band indices,
      "kpoint_index": Valence band maximum k-point indices,
      "kpoint_labels": Valence band maximum k-point labels,
      "kpoint_coord": Valence band maximum k-point coordinates,
      "energy": Valence band maximum energy (eV)
    }
  }
}
```

## band_index Explanation

`band_index` represents the band numbers (starting from 1) that lie at the Conduction Band Minimum (CBM) or Valence Band Maximum (VBM) for different spin channels and at different k-points.

Data structure is a 3-dimensional array: `[nspin][nkpoint][nband]`

- **First dimension**: Spin channel index (nspin=1 for non-spin-polarized, nspin=2 for spin-polarized)
- **Second dimension**: Index of k-points where CBM/VBM are found
- **Third dimension**: List of band numbers that lie at CBM/VBM at that k-point

For example, `[[[1, 2, 3], [1, 2, 3]]]` in `vbm.band_index` indicates:
- Spin channel 1 (non-spin-polarized)
- At two different k-points
- Bands 1, 2, and 3 all lie at the valence band maximum
