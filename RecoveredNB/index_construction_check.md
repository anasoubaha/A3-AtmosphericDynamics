# Construction order of MO and WeMO: does it inflate the rainfall link?

Evidence from [`index_construction_check.py`](index_construction_check.py); numbers in
[`output/index_construction_check.csv`](output/index_construction_check.csv).
Both indices are currently built as *difference the poles in Pa, then standardise*
(Cells 5–6). The question is whether *standardise each pole, then difference* is the
canonical form, and what the order does to the results.

## The mechanism is real and large

Differencing raw pressure weights both poles equally in hPa, so the higher-variance pole
dominates. In DJF over 1991–2020:

| Index | West pole | East pole | Variance ratio W/E |
|---|---|---|---|
| MO | Gibraltar **3.84 hPa** | Lod 1.74 hPa | **4.88** |
| WeMO | Cádiz 3.84 hPa | Padua **4.91 hPa** | 0.61 |

So the concern is well founded — but it cuts in opposite directions for the two indices.
For MO the storm-track pole dominates; for WeMO it is the *eastern* (Padua) pole that does.

## But the two published indices are defined differently

The published MOI2 daily series carries an exact normalisation fingerprint: **mean −0.0001 and
sd 1.0000 over 1961–1990**, against 1.0128 (1951–80) and 0.9635 (1981–2010). Its sd by calendar
month ranges 0.43 (Aug) to 1.34 (Dec), so it was normalised **once over the whole series, not
per calendar month**. A difference of two separately standardised series would not land on
sd = 1.0000 by construction. CRU/Palutikof's MOI2 is therefore **difference first, then
normalise** — the order already in use here.

Martín-Vide & López-Bustins (2006) define WeMOi explicitly as the difference of the
**standardised** values: standardise first. The two published indices genuinely differ in
construction, and testing each form against its own published series confirms it:

| | vs published, current (`diff_then_std`) | vs published, `z_first` |
|---|---|---|
| **MO** vs `moi2.output.dat` | **r = 0.894** | r = 0.832 |
| **WeMO** vs `Web_WeMOi-2020.txt` | r = 0.856 | **r = 0.891** |

(both series monthly-standardised before correlating; raw-form correlations show the same
ordering). Per season the WeMO improvement is substantial: DJF 0.962 → **0.983**,
SON 0.859 → **0.926**, MAM 0.851 → **0.908**.

**Verdict: MO is already canonical; WeMO is not.**

## What the order does to the rainfall correlations

| Index | Season | current (`diff_then_std`) | `z_first` | west pole alone |
|---|---|---|---|---|
| MO | SON | −0.727 | −0.547 | −0.735 |
| MO | **DJF** | **−0.912** | **−0.597** | −0.904 |
| MO | MAM | −0.380 | −0.246 | −0.483 |
| WeMO | SON | −0.500 | −0.652 | −0.720 |
| WeMO | **DJF** | −0.392 | **−0.567** | −0.905 |
| WeMO | MAM | −0.342 | −0.534 | −0.563 |

Two things follow.

**For MO, construction order accounts for most of the headline correlation.** Under `z_first`
the DJF link falls from −0.912 to −0.597 (R² 0.83 → 0.36), and the index's attachment to its
western pole falls from r = 0.925 to r = 0.452. So the local-barometer behaviour is largely a
*consequence of the differencing order*. Crucially, that order is the published definition's,
not an implementation error: **MOI2 as defined is dominated by its Gibraltar pole in winter**,
because differencing raw pressures weights a 3.84 hPa pole against a 1.74 hPa one. The critique
therefore lands on the published index, not on this reconstruction — which also means a
reviewer cannot object that the index is not the one cited. It is.

**For WeMO, adopting the canonical form strengthens the index** against precipitation in every
season (DJF −0.392 → −0.567). Its current non-canonical form is *understating* WeMO, so this is
a correction that would change the regression results.

## The western poles are effectively one barometer

MO's western pole (Gibraltar, 36.0 °N 5.0 °W) and WeMO's (Cádiz, 36.25 °N 6.0 °W) are **94 km
apart**, 34 km and 50 km from Morocco respectively, and their standardised DJF series correlate
at **r = 0.9929**. The western pole alone gives −0.904 (MO file) and −0.905 (WeMO file) against
Moroccan DJF rainfall — the same number twice, because it is the same measurement. Any framing
that treats MO and WeMO as two independent teleconnection predictors of Moroccan rainfall has
to reckon with the fact that both are anchored on a barometer sitting just off the Moroccan
coast.

## Two further deviations from CRU, for the record

The MO index here standardises **per calendar month over 1981–2010**; CRU normalises **once
over 1961–1990**. The per-calendar-month choice is defensible and probably preferable (it
removes the annual variance cycle that would otherwise distort pooled multi-season
correlations — it is why the Cell 9 validation rises from r = 0.67 to r = 0.91 when both series
are put on the same footing), but it should be stated in the Methods rather than left implicit,
since it means the series is not numerically identical to the distributed MOI2.

## What it does to the regression — three parallel builds

Cells 5 and 7 now build each index as `z(west) − z(east)` with a hard pole-coverage guard,
writing to **new** files (`WeMO_index_zfirst.nc`, `MO_index_MOI2_zfirst.nc`); the originals are
untouched. `teleconnection_variance.py --indices {published,zfirst,canonical}` runs any of
three builds into its own folder:

| Build | MO | WeMO | Output |
|---|---|---|---|
| `published` | diff-then-std | diff-then-std | `output/` |
| `zfirst` | z-first | z-first | `output_zfirst/` |
| `canonical` | diff-then-std (= CRU MOI2) | z-first (= Martín-Vide WeMOi) | `output_canonical/` |

Adjusted R², raw series:

| Season | Model | published | z-first both | canonical |
|---|---|---|---|---|
| SON | M1 indices | 0.636 | 0.486 | 0.653 |
| SON | M3 all five | 0.655 | 0.615 | 0.665 |
| **DJF** | **M1 indices** | **0.808** | **0.591** | **0.808** |
| **DJF** | **M3 all five** | **0.843** | **0.817** | **0.843** |
| MAM | M1 indices | 0.272 | 0.318 | 0.370 |
| MAM | M3 all five | 0.579 | 0.634 | 0.623 |

(M2, Z850 alone, is identical in all three by construction: 0.509 / 0.658 / 0.377.)

Increments in DJF: ΔR²_indices +0.184 → +0.158 → +0.185; ΔR²_circulation +0.035 → **+0.226** →
+0.035.

**Two things follow.**

*The canonical build barely moves DJF.* Fixing WeMO changes the DJF headline not at all
(M1 0.808, M3 0.843, ΔR²_indices +0.185, all within rounding of the published build), because
MO already carries that variance. The WeMO correction mainly improves **spring**: MAM M1 rises
0.272 → 0.370 and M3 0.579 → 0.623. So the published DJF result is not an artifact of the WeMO
construction error.

*The `zfirst` build is the sensitivity test, not a candidate for the paper.* Standardising MO's
poles first is a departure from CRU's MOI2 definition, so it should not be used for a result
that cites Palutikof. But as a diagnostic it is decisive: with MO's poles equalised, DJF M1
falls from 0.808 to 0.591, ΔR²_circulation rises from +0.035 to **+0.226**, DJF max VIF drops
from 4.16 to **2.23**, and in M3 the largest coefficient becomes **Z850** (−7.91 ± 1.33,
p < 0.001) rather than MO (−4.14 ± 1.14, p = 0.001). In other words, roughly a quarter of the
DJF variance that the published build attributes to "teleconnection indices" is attributable to
MOI2's differencing order concentrating the index on its Gibraltar pole. That is a property of
the published index definition, and it belongs in the discussion.

## The Algiers–Cairo (MOI1) variant is blocked on data

`era5-slp.nc` spans **lat 31.75–36.25 °N, lon 5.25 °W–34.75 °E** — it contains the MOI2 poles
but neither MOI1 pole. Run as written with `method='nearest'`, the Cell 7 code would have
snapped silently:

- Algiers (36.5 °N, 3.0 °E) → 36.25 °N, **28 km** south — tolerable;
- Cairo (30.0 °N, 31.5 °E) → **31.75 °N**, 31.5 °E — **194 km north, in the Mediterranean off
  the Nile Delta**. The resulting "MOI1" would have been meaningless.

Cell 7 now raises instead of snapping. To compute MOI1, request ERA5 monthly-mean
`mean_sea_level_pressure` on a box such as **lat 28–38 °N, lon 1–34 °E**, 1950–2024 (a superset
covering both variants), from the same CDS dataset already cited in Cell 2.

The one other file that covers both poles, `surfacepressure_80N0S50W80E.nc`, is not a
substitute: it is *surface* pressure, it starts in 1991 (excluding the 1981–2010 baseline), and
its grid points at the MOI1 poles sit at roughly **435 m** (Algiers) and **222 m** (Cairo)
elevation, where surface pressure is not interchangeable with MSLP. At genuinely sea-level
points the substitution is harmless — checked at Gibraltar and Lod, where standardised `sp` and
`msl` correlate at r = 0.9999 and the resulting MO differs by r = 0.9999 — but that check does
not license using it 435 m up. Note also that the nearest sea-level grid cells to the CRU
stations are around (36.75 °N, 2.75 °E) for Algiers and (30.5 °N, 31.25 °E) for Cairo, which is
where the poles should be placed once a proper MSLP file is available.
