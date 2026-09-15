# MOI1 (Algiers–Cairo) vs MOI2 (Gibraltar–Lod)

Built from `mo_era5_mslp.nc` (lat 28–42 °N, lon 10 °W–35 °E, 1981–2025) via Cell 7, which now
carries a `VARIANT` switch (MOI1/MOI2) and a `CONSTRUCTION` switch
(`z_first` / `diff_first`). Four index files exist; nothing was overwritten.

**Data check first.** Rebuilding MOI2 in `diff_first` form from the *new* file reproduces the
original `MO_index.nc` (built from `era5-slp.nc`) to **max |difference| = 1.8e-07**, so the new
download is consistent with the old one and the two variants are comparable. Each variant's
pole decomposition is verified against its own index file at **max |difference| = 0.0** before
anything is read off it.

## Pole geography

| Pole | Distance from Morocco |
|---|---|
| MOI2 west — Gibraltar (36.0 °N, 5.0 °W) | **34 km** |
| MOI1 west — Algiers (36.5 °N, 3.0 °E) | **471 km** |
| MOI1 east — Cairo (30.0 °N, 31.5 °E) | 3,096 km |
| MOI2 east — Lod (32.0 °N, 34.5 °E) | 3,330 km |

The two indices correlate at **r = 0.935** monthly over 1981–2024 (DJF 0.945, SON 0.952,
MAM 0.918) — they are near-substitutes despite the 437 km displacement of the western pole.

## Results, canonical construction (CRU `diff_first` MO + Martín-Vide WeMO)

| Season | Model | MOI2 canonical | **MOI1 canonical** |
|---|---|---|---|
| SON | M1 indices | 0.653 | 0.633 |
| SON | M3 all five | 0.665 | 0.638 |
| **DJF** | **M1 indices** | **0.808** | **0.784** |
| **DJF** | **M3 all five** | **0.843** | **0.820** |
| MAM | M1 indices | 0.370 | **0.424** |
| MAM | M3 all five | 0.623 | 0.615 |

DJF: r(MO, precip) −0.912 → **−0.854**; ΔR²_indices +0.185 → +0.162; ΔR²_circulation +0.035 →
+0.036; MO's VIF 4.40 → **3.31**. MOI1 is slightly weaker in DJF and SON, slightly stronger in
MAM, and marginally less collinear with Z850 throughout.

## The framing does not change — and that is the finding

Moving the western pole 437 km further from Morocco does **not** break the local-barometer
behaviour:

(1991–2020 index baseline, matching the canonical run; files
`output_canonical_moi{1,2}/mo_pole_diagnostic_MOI{1,2}.csv`)

| DJF | MOI2 | MOI1 |
|---|---|---|
| r(index, precip) | −0.917 | −0.854 |
| r(west pole alone, precip) | −0.908 | −0.839 |
| r(index, its own west pole) | +0.925 | **+0.935** |
| r(west pole, Z850) | +0.879 | +0.831 |
| east pole added to west-pole regression | +0.046 (p=0.002) | +0.040 (p=0.024) |

Both variants are dominated by their western pole to the same degree, and in SON the eastern
pole is again indistinguishable from zero (MOI1 r = −0.019, p = 0.92; MOI2 r = +0.017, p = 0.93),
while in MAM it again subtracts skill (−0.028 adjusted R², p = 0.97). So the earlier conclusion
was not an artifact of Gibraltar's unusual proximity to Morocco. Western-Mediterranean sea-level
pressure anomalies are spatially coherent over ~1000 km, so Algiers carries nearly the same
information as Gibraltar.

**The right description is therefore "western Mediterranean sea-level pressure", not
"Gibraltar barometer".** Either MO variant is, for Moroccan winter rainfall, a western-basin SLP
index whose eastern pole contributes little — MOI1 just makes that harder for a reviewer to
dismiss as a geographic accident. The choice between variants is close to immaterial for the
results: use MOI2 to match the series you validated against (`moi2.output.dat`), or MOI1 if you
prefer a pole that is not sitting on Morocco's doorstep, and report the other as a sensitivity.

**Not yet validated.** There is no published MOI1 series in `Article3_data/` to check the
Algiers–Cairo reconstruction against; `moi1.output.dat` from the same CRU/UEA source would be
needed to repeat the Cell 9 validation for this variant.

## Output folders

| Build | MO file | Folder |
|---|---|---|
| `published` | `MO_index.nc` (MOI2, diff-first) | `output/` |
| `canonical` | `MO_index.nc` | `output_canonical/` |
| `zfirst` | `MO_index_MOI2_zfirst.nc` | `output_zfirst/` |
| `moi1_canonical` | `MO_index_MOI1_difffirst.nc` | `output_moi1_canonical/` |
| `moi1_zfirst` | `MO_index_MOI1_zfirst.nc` | `output_moi1_zfirst/` |

All five use the same precipitation, Z850, NAO and Niño3.4 inputs; `published` uses the
original WeMO, the other four use `WeMO_index_zfirst.nc`.
