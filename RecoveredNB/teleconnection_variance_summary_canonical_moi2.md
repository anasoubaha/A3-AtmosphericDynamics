# Restructured predictor blocks — MOI2 (Gibraltar–Lod)

Run: `python teleconnection_variance.py --indices canonical_moi2` → [`output_canonical_moi2/`](output_canonical_moi2/).
Canonical construction (MO normalised after differencing per Palutikof; WeMO standardised
per station then differenced per Martín-Vide), both indices renormalised on **1991–2020**.
Full Morocco domain (1017 cells, verified), 1991–2024, DJF n=33 / SON n=34 / MAM n=34.

MO is placed in the **circulation** block, not the teleconnection block, following the pole
diagnostic: the index reproduces its own western pole (r = 0.925 MOI2 / 0.935 MOI1 in DJF),
that pole correlates 0.83–0.88 with Morocco-mean Z850, and the eastern pole is ≈0 in SON and
subtracts skill in MAM. Those figures are in `mo_pole_diagnostic_MOI2.csv` in this
same folder, computed on the same 1991–2020 index baseline as the tables below.

## DJF adjusted R², all five models

| Model | Predictors | R²adj |
|---|---|---|
| T | NAO + WeMO + Niño3.4 | 0.536 |
| C | Z850 + MO | 0.861 |
| F | all five | 0.857 |
| T_conv | NAO + MO + WeMO + Niño3.4 | 0.820 |
| C_conv | Z850 | 0.658 |

## Both pairs of increments (adjusted R², raw series)

| Season | dR2_teleconnections (F−C) | dR2_circulation (F−T) | dR2_indices_conv (F−C_conv) | dR2_circulation_conv (F−T_conv) |
|---|---|---|---|---|
| SON | +0.099 | +0.270 | +0.153 | +0.013 |
| **DJF** | **-0.004** | **+0.322** | +0.199 | +0.037 |
| MAM | +0.231 | +0.312 | +0.211 | +0.292 |

Detrended DJF: dR2_teleconnections = -0.006,
dR2_circulation = +0.352.

The two groupings are not interchangeable and must not be quoted together: with MO counted as
circulation the teleconnection block adds **-0.004** in DJF; with MO counted as an index
the same full model appears to show indices adding +0.199. The difference is entirely
where MO sits.

## Significant coefficients in F (DJF, mm/month per σ)

| Predictor | Coefficient | SE | t | p |
|---|---|---|---|---|
| **MO** | -7.94 | 1.65 | -4.82 | 0.000 |
| **Z850** | -4.10 | 1.42 | -2.89 | 0.008 |
| NAO | +0.67 | 1.19 | +0.56 | 0.578 |
| WeMO | -0.37 | 1.04 | -0.36 | 0.723 |
| Nino34 | +1.27 | 0.90 | +1.42 | 0.168 |

Retaining significance at p < 0.05: **MO, Z850**. Not significant: NAO, WeMO, Nino34.

## Collinearity

Maximum VIF in F (DJF) = **4.40**, carried by **MO**. Full DJF set:
MO 4.40, Z850 3.27, NAO 2.30, WeMO 1.77, Nino34 1.31. All below the conventional threshold of 5.

## Does any conclusion differ between MOI2 and MOI1?

**No conclusion differs; two details do.** Under both variants the circulation block dominates
(dR2_circulation +0.322 MOI2 / +0.299 MOI1 in DJF), the teleconnection block adds essentially
nothing in winter (−0.004 / +0.034), MO is the largest coefficient in F, Niño3.4 is the only
predictor failing Benjamini-Hochberg in every season, and all VIFs stay under 5.

Two differences, neither changing an interpretation:

1. **WeMO's DJF coefficient.** Under MOI2 it is −0.37 ± 1.04 (p=0.72, not significant); under
   MOI1 it is −2.85 ± 1.04 (p=0.011, significant). MOI2's Gibraltar pole is 94 km from WeMO's
   Cádiz pole and the two correlate at r=0.993, so MOI2 absorbs the variance WeMO would
   otherwise carry. With the Algiers pole 471 km further east, WeMO keeps its own share. If
   WeMO's winter role matters to the argument, say which variant produced it.
2. **Where the collinearity sits.** MOI2: max VIF 4.40 on MO. MOI1: 3.42 on Z850, with MO at
   3.16 — lower overall, as expected from the more distant western pole.

DJF adjusted R²(F) is 0.857 (MOI2) against 0.835 (MOI1); C is 0.861 against 0.801.

## Note on R²adj(F) < R²adj(C) in DJF

The full model's adjusted R² (0.857) is marginally **below** the circulation block alone
(0.861). This is legitimate, not an error: adding NAO, WeMO and Niño3.4 to Z850+MO raises raw
R² (0.877 vs 0.867) but not enough to pay for three degrees of freedom at n=33. The script
prints this rather than failing, and enforces the raw-R² ordering separately. It is also the
cleanest statement of the result — in winter the three teleconnection indices contribute
nothing that the two circulation diagnostics have not already captured.
