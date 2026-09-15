# Teleconnection vs circulation control on Moroccan seasonal precipitation

Notes on the analysis produced by [`teleconnection_variance.py`](teleconnection_variance.py);
outputs in [`output/`](output/). Period 1991–2024, seasons SON / DJF / MAM, DJF primary.

## One thing broke the n=34 requirement

**PERSIANN-CDR February 1992 is missing over the entire grid** (all cells flagged −99, not NaN — so an unguarded mean would have silently returned a large negative "anomaly"). DJF 1992 therefore cannot be formed from three months and **DJF runs at n=33**. SON and MAM are both n=34.

I didn't let this hard-crash the run, because that would deliver nothing. Instead the completeness check aborts on *any* undocumented shortfall and passes only gaps listed in a `DOCUMENTED_GAPS` dict at the top of the file, where DJF 1992 is recorded with its reason. The shortfall prints at three stages and `n` is written into every output row. If you'd rather force `n==34`, empty that dict and the run stops. Feb 1990 is also fully missing but falls outside the window (only Dec 1990 is needed for DJF 1991, and it's present).

### Verification of the gap and of every partial fill

The gap is in the product as distributed, not in the local copy. Two exports held in
`Article3_data/` — `CDR.nc` (Dec 2023, 68x73 cells, 20.25-37 N / 18-0 W) and
`CDR-1983-2025.nc` (Sep 2025, 78x86 cells, 19.25-38.5 N / 20 W-1.25 E) — were pulled
independently about two years apart over different spatial windows. Over their common grid
they are **bit-identical**: max |difference| = 0.0 across 2,417,468 pixel-months, and they
disagree on the missing-value flag in **zero** pixel-months. Both flag exactly the same five
months: 1984-02, 1984-05 and 1984-06 (partial) and 1990-02 and 1992-02 (whole grid).

The partial months also carry a structure no failed transfer would produce: the western edge
of the missing block steps eastward with latitude (-6.25 deg at 38.5 N, -3.75 deg through the
mid-latitudes, -1.25 deg at 20.5 N), i.e. a curved satellite viewing-geometry boundary in the
ISCCP B1 / GridSat-B1 infrared source, not a truncated file. Note that neither the NCEI
product landing page nor the NCAR Climate Data Guide entry publishes a gap inventory for
PERSIANN-CDR, so this rests on reproducibility across independent downloads plus the physical
structure of the gap, not on a published list of missing months.

**No partial fill was ever averaged.** All -99 cells are converted to NaN before any spatial
mean is taken, and the cos-lat weights are renormalised over the valid cells of each time step.
More to the point, a census of the window that actually feeds the analysis (1990-12 to 2024-12)
returns **zero partial-coverage months** in either domain: the only affected month is 1992-02 at
100% missing, which is dropped as an incomplete season rather than averaged. The three partial
months predate the window by seven years and fall outside the 1991-2020 baseline as well. The
single residual effect of the gap is that the February climatology rests on 29 years rather
than 30 (all other calendar months: 30).

## DJF results (primary)

| Model | R² | R²adj | F | p |
|---|---|---|---|---|
| M1 indices | 0.832 | 0.808 | F(4,28)=34.65 | 1.8e-10 |
| M2 Z850 | 0.669 | 0.658 | F(1,31)=62.65 | 6.2e-09 |
| M3 all five | 0.867 | 0.843 | F(5,27)=35.27 | 5.1e-11 |

ΔR²_indices = **+0.184**, ΔR²_circulation = **+0.035**. Detrended: +0.209 / +0.025 — the increments are not trend artifacts. Correlations with precipitation: MO −0.91, Z850 −0.82, NAO −0.68, WeMO −0.39, Niño3.4 −0.04. All but Niño3.4 survive BH at q=0.05 across all 15 tests. Lag-1 autocorrelation of DJF precipitation is +0.15 (p=0.41), so the season-years are effectively independent.

## MO is a local barometer, not a teleconnection

**This supersedes the earlier framing of this section.** MO here is the CRU/UEA **MOI2**
variant: the normalised pressure difference between Gibraltar's Northern Frontier
(36.1 N, 5.3 W) and Lod Airport, Israel (32.0 N, 34.5 E). (MOI1 is the Algiers-Cairo pair.)
The poles were recovered by exhaustive search against `MO_index.nc`; the western one sits on
the ERA5 grid point at **36.0 N, 5.0 W**, in the Strait of Gibraltar, **34 km from the
Moroccan coastline** - the nominal CRU station coordinate is closer still, 22 km - 53 km from the nearest precipitation grid
cell entering the domain mean, with 51 of the 1017 domain cells inside 200 km. The eastern
pole, the ERA5 point nearest Lod Airport (32.0 N, 34.5 E), is 3,330 km away.

Decomposing the index into its two poles ([`mo_pole_diagnostic.py`](mo_pole_diagnostic.py),
results in [`output/mo_pole_diagnostic_MOI2.csv`](output/mo_pole_diagnostic_MOI2.csv)) shows the
western pole carries essentially the whole relationship. **These numbers are on the 1981-2010
index baseline, matching the `published` run in `output/`;** the 1991-2020 canonical values are
in `output_canonical_moi2/` and differ slightly in MAM (see
[`mo_variant_comparison.md`](mo_variant_comparison.md)):

| Season | r(MO, precip) | r(z SLP_west, precip) | r(z SLP_east, precip) | r(SLP_west, Z850) |
|---|---|---|---|---|
| SON | −0.727 | **−0.735** | +0.017 | +0.795 |
| DJF | −0.912 | **−0.904** | −0.430 | +0.885 |
| MAM | −0.380 | **−0.483** | −0.232 | +0.756 |

In SON the pattern is exactly what a local barometer predicts: the western pole reproduces MO
(−0.735 vs −0.727) and the eastern pole is indistinguishable from zero (+0.017, p=0.93). In
MAM the western pole is *stronger* than MO (−0.483 vs −0.380), so differencing against the
east actively degrades the signal. In DJF the eastern pole is not zero (−0.430) and does add
to a regression on the western pole (+0.051 adjusted R2, p=0.001) — but note its **sign**: it
correlates negatively with precipitation, the same direction as the western pole. A genuine
east-west seesaw would require the opposite sign. What the eastern pole is doing is removing a
pressure component common to both ends of the basin (the two poles correlate at +0.67 in DJF)
that is irrelevant to Moroccan rainfall — it acts as a suppressor variable, not as an
independent eastern-Mediterranean signal.

The western pole also correlates with the Morocco-mean 850 hPa height at **+0.885** in DJF.
MO and Z850 are therefore close to the same measurement, made at slightly different places and
levels; MO wins because a point barometer at the Strait is better placed for Moroccan frontal
rainfall than a height field averaged over the whole country including the Sahara.

**Consequence for the manuscript.** MO cannot be described as a teleconnection index in this
analysis. Sea-level pressure sitting directly over a region in winter is anti-correlated with
that region's rainfall almost tautologically — high pressure means no fronts — so
r(MO, precip) = −0.91 is close to a definition, not a finding. The earlier result that
"the four indices add ΔR2adj = +0.184 beyond the circulation diagnostic" is really
"a better-placed local barometer beats a domain-averaged height field by +0.184", which is a
statement about diagnostic placement, not about remote forcing. The genuinely remote indices —
NAO, WeMO, Niño3.4 — all collapse to ≈0 in M3 with standard errors larger than their
estimates (NAO +0.37 ± 1.24, WeMO −0.02 ± 0.93, Niño3.4 +1.13 ± 0.96), and that, not the MO
result, is the defensible teleconnection finding here.

**The Cell 9 validation is sound.** `moi2.output.dat` is the published MOI2 (Gibraltar-Lod)
series, the same variant as `MO_index.nc`, so the r = 0.91 monthly / 0.98 DJF-mean agreement
reported there is a like-for-like validation of the reconstruction. The residual difference is
what you would expect from substituting ERA5 grid points for the Gibraltar and Lod station
barometers. Nothing in this section questions that the index is faithfully reconstructed - the
question is what the index measures.

## Deviations from the spec, all deliberate

- **Two extra CSVs.** `coefficients.csv` holds the coefficients/SEs/t/p/CIs that (d) asks you to report — the `results_table.csv` schema you specified has no room for a variable-length coefficient list. `diagnostics.csv` holds the lag-1 statistics from (g).
- **Detrended models live in `results_table.csv`** under a `detrended` boolean column (18 rows = 3 seasons × 3 models × 2), rather than in a separate file.
- **ΔR² values appear only on the M3 rows**, since both increments are properties of that comparison; M1/M2 rows carry NaN there.
- **BH is reported twice** — `p_bh_all` across all 15 tests, and `p_bh_within_season` across each season's 5 — because "across all index-precipitation tests" is ambiguous and the choice changes WeMO's MAM verdict.
- **Detrending regresses on season-year, not array position**, so the DJF gap doesn't distort the trend estimate.

Sanity checks that passed: Morocco-mean climatology is 189 mm/yr with a Nov–Mar maximum, Z850 heights are 1502–1540 gpm (confirming the ÷g conversion), 1017 grid cells inside the polygon on both grids, and both files are on −180/180 longitudes (the 0–360 branch is implemented but not exercised — stated in the printed output). Note the shapefile includes Western Sahara, which is why the domain mean is drier than Morocco proper.
