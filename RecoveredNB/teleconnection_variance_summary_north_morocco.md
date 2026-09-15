# Teleconnection vs circulation control — northern Morocco (Western Sahara excluded)

> **Status: retired sensitivity test.** The analysis proceeds on the full Morocco domain
> (`output/`); this run is kept only as the record that excluding Western Sahara changes nothing.

Second run of [`teleconnection_variance.py`](teleconnection_variance.py), invoked as
`python teleconnection_variance.py --domain north_morocco`; outputs in
[`output_north_morocco/`](output_north_morocco/). Period 1991–2024, seasons SON / DJF / MAM, DJF primary.
Companion to [`teleconnection_variance_summary.md`](teleconnection_variance_summary.md), which covers
the full-Morocco domain (`output/`).

Domain: GADM `MAR_adm0`, Morocco excluding Western Sahara — **639 grid cells** at 0.25° spanning
27.75–35.75 °N, 13.00–1.25 °W, against 1017 cells spanning 20.75–35.75 °N for the original domain.
Everything else — inputs, baseline, seasons, models, tests — is held fixed, so the two runs differ
only in which cells enter the two domain means.

## Western Sahara was not biasing the result

Excluding it makes every headline number **slightly stronger, never weaker**, and changes no conclusion. The arid south was diluting the signal, not creating it.

| DJF, raw | Full Morocco | Northern Morocco |
|---|---|---|
| M1 indices R²adj | 0.808 | **0.828** |
| M2 Z850 R²adj | 0.658 | **0.696** |
| M3 all five R²adj | 0.843 | **0.852** |
| ΔR²_indices | +0.184 | +0.156 |
| ΔR²_circulation | +0.035 | +0.024 |
| MO vs precipitation | −0.912 | **−0.921** |
| Z850 vs precipitation | −0.818 | **−0.840** |
| NAO vs precipitation | −0.676 | **−0.698** |

Every predictor keeps its sign, its rank and its Benjamini-Hochberg verdict in all three seasons. The one marginal case sits in MAM, where WeMO moves from p=0.048 to p=0.058 raw — it failed BH correction in both runs anyway. Niño3.4 remains non-significant everywhere. The 4×4 index correlation matrix is identical between runs by construction, since the indices do not depend on the domain.

## DJF results (primary), northern Morocco

| Model | R² | R²adj | F | p |
|---|---|---|---|---|
| M1 indices | 0.849 | 0.828 | F(4,28)=39.49 | 4.0e-11 |
| M2 Z850 | 0.706 | 0.696 | F(1,31)=74.37 | 9.7e-10 |
| M3 all five | 0.875 | 0.852 | F(5,27)=37.90 | 2.2e-11 |

ΔR²_indices = **+0.156**, ΔR²_circulation = **+0.024**. Detrended: +0.180 / +0.015. Correlations with precipitation: MO −0.92, Z850 −0.84, NAO −0.70, WeMO −0.39, Niño3.4 −0.01; all but Niño3.4 survive BH at q=0.05 across all 15 tests. Lag-1 autocorrelation of DJF precipitation is +0.17 (p=0.36), so the season-years remain effectively independent.

## The result a reviewer will attack — unchanged, and slightly sharper

In M3, **MO is still the only index with a significant coefficient** (−11.19 mm/month per σ, SE 2.55, p<0.001); NAO (−0.02 ± 1.76), WeMO (−0.18 ± 1.33) and Niño3.4 (1.75 ± 1.31) all sit at or near zero. MO now correlates with Z850 at r=0.83 (was 0.78), and the two highest VIFs rise accordingly to 4.79 (MO) and 4.02 (Z850) — still under the conventional threshold of 5, but closer to it. Tightening the domain to where the rain actually falls makes the Gibraltar–Lod gradient and the local 850 hPa height *more* alike, which sharpens rather than softens the caveat: ΔR²_indices is very largely MO adding skill, and MO is itself a circulation diagnostic drawn from the same reanalysis as Z850. **Superseded:** MO is the CRU/UEA MOI2 (Gibraltar-Lod) variant and its western pole sits 34 km off the Moroccan coast, so this is a local-barometer effect rather than a teleconnection — see the pole decomposition in [`teleconnection_variance_summary.md`](teleconnection_variance_summary.md). The domain comparison below stands; the causal reading of MO does not. NAO, WeMO and ENSO add essentially nothing once MO is in the model, in either domain.

Season structure is also unchanged: SON behaves like DJF (MO r=−0.77, Z850 r=−0.76), and MAM is still the one season where circulation dominates the indices (ΔR²_circulation=+0.358, up from +0.307).

## Notes specific to this run

- **n is 33 for DJF here too.** The PERSIANN-CDR February 1992 gap is global to the product, not domain-specific, so DJF 1992 drops out identically. SON and MAM are n=34. The same `DOCUMENTED_GAPS` guard applies. The gap was verified against a second, independently downloaded export (bit-identical over the common grid, same five flagged months); a census of 1990-12 to 2024-12 returns **zero partial-coverage months** in this domain, so no renormalised partial fill enters any season here either. Details in [`teleconnection_variance_summary.md`](teleconnection_variance_summary.md).
- **The script was parameterized rather than copied.** `--domain {morocco,north_morocco}` selects the shapefile and the output folder; `morocco` remains the default and was re-run to confirm it reproduces `output/` byte-for-byte, so the first set of results is untouched.
- **Both domains produce the same six CSVs** (`results_table`, `correlations`, `index_correlation_matrix`, `vif`, `coefficients`, `diagnostics`), in their own folders, with the domain recorded in the printed header.
- Sanity checks: the northern domain climatology is **268 mm/yr** against 189 mm/yr for the full domain — the expected wetting once the Saharan south is removed — with the wettest months November, December and March. Domain-mean precipitation variability rises correspondingly (monthly anomaly σ 16.95 vs 11.97 mm/month), as does Z850 variability (17.61 vs 14.38 gpm). Both grids place 639 cells inside the polygon, and the shapefile is EPSG:4326 on −180/180 longitudes.
