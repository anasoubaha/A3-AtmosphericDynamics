# Restructured predictor blocks — MOI1 (Algiers–Cairo)

Run: `python teleconnection_variance.py --indices canonical_moi1` → [`output_canonical_moi1/`](output_canonical_moi1/).
Canonical construction (MO normalised after differencing per Palutikof; WeMO standardised
per station then differenced per Martín-Vide), both indices renormalised on **1991–2020**.
Full Morocco domain (1017 cells, verified), 1991–2024, DJF n=33 / SON n=34 / MAM n=34.

MO is placed in the **circulation** block, not the teleconnection block, following the pole
diagnostic: the index reproduces its own western pole (r = 0.925 MOI2 / 0.935 MOI1 in DJF),
that pole correlates 0.83–0.88 with Morocco-mean Z850, and the eastern pole is ≈0 in SON and
subtracts skill in MAM. Those figures are in `mo_pole_diagnostic_MOI1.csv` in this
same folder, computed on the same 1991–2020 index baseline as the tables below.

## DJF adjusted R², all five models

| Model | Predictors | R²adj |
|---|---|---|
| T | NAO + WeMO + Niño3.4 | 0.536 |
| C | Z850 + MO | 0.801 |
| F | all five | 0.835 |
| T_conv | NAO + MO + WeMO + Niño3.4 | 0.798 |
| C_conv | Z850 | 0.658 |

## Both pairs of increments (adjusted R², raw series)

| Season | dR2_teleconnections (F−C) | dR2_circulation (F−T) | dR2_indices_conv (F−C_conv) | dR2_circulation_conv (F−T_conv) |
|---|---|---|---|---|
| SON | +0.137 | +0.241 | +0.123 | +0.006 |
| **DJF** | **+0.034** | **+0.299** | +0.176 | +0.037 |
| MAM | +0.213 | +0.306 | +0.205 | +0.245 |

Detrended DJF: dR2_teleconnections = +0.025,
dR2_circulation = +0.325.

The two groupings are not interchangeable and must not be quoted together: with MO counted as
circulation the teleconnection block adds **+0.034** in DJF; with MO counted as an index
the same full model appears to show indices adding +0.176. The difference is entirely
where MO sits.

## Significant coefficients in F (DJF, mm/month per σ)

| Predictor | Coefficient | SE | t | p |
|---|---|---|---|---|
| **MO** | -6.07 | 1.50 | -4.04 | 0.000 |
| **WeMO** | -2.85 | 1.04 | -2.75 | 0.011 |
| **Z850** | -4.22 | 1.56 | -2.70 | 0.012 |
| NAO | +0.30 | 1.27 | +0.23 | 0.817 |
| Nino34 | +1.21 | 0.98 | +1.23 | 0.230 |

Retaining significance at p < 0.05: **MO, WeMO, Z850**. Not significant: NAO, Nino34.

## Collinearity

Maximum VIF in F (DJF) = **3.42**, carried by **Z850**. Full DJF set:
Z850 3.42, MO 3.16, NAO 2.26, WeMO 1.51, Nino34 1.35. All below the conventional threshold of 5.

## Does any conclusion differ between MOI1 and MOI2?

**No conclusion differs; two details do.** Under both variants the circulation block dominates
(dR2_circulation +0.299 MOI1 / +0.322 MOI2 in DJF), the teleconnection block adds little in
winter (+0.034 / −0.004), MO is the largest coefficient in F, Niño3.4 is the only predictor
failing Benjamini-Hochberg in every season, and all VIFs stay under 5.

Two differences, neither changing an interpretation:

1. **WeMO's DJF coefficient.** Under MOI1 it is −2.85 ± 1.04 (p=0.011, significant); under
   MOI2 it is −0.37 ± 1.04 (p=0.72, not significant). MOI2's Gibraltar pole sits 94 km from
   WeMO's Cádiz pole (r=0.993 between them), so MOI2 absorbs variance that WeMO would
   otherwise carry; MOI1's Algiers pole is 471 km further east and leaves WeMO its own share.
   This run is the one in which WeMO retains an independent winter signal.
2. **Where the collinearity sits.** MOI1: max VIF 3.42 on Z850, MO at 3.16. MOI2: 4.40 on MO.
   MOI1 is the less collinear configuration.

DJF adjusted R²(F) is 0.835 (MOI1) against 0.857 (MOI2); C is 0.801 against 0.861. MOI1 does
not reproduce the MOI2 result that R²adj(F) falls marginally below R²adj(C) — here the full
model is the better of the two on adjusted R², by +0.034.
