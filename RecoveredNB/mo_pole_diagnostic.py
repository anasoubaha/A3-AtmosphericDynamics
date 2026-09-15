"""
Is the MO index acting as a teleconnection, or as a local barometer?

MO here is the CRU/UEA MOI2 (Gibraltar-Lod) variant: a normalised mean-sea-level-pressure
difference between a western and an eastern pole. If the western pole alone
reproduces the MO-precipitation correlation and the eastern pole contributes
nothing, then MO is not measuring an east-west Mediterranean seesaw over Morocco;
it is measuring the pressure sitting on top of Morocco, which is anti-correlated
with Moroccan rainfall almost by definition.

This script decomposes MO into its two poles and correlates each against the same
seasonal precipitation series used by teleconnection_variance.py.

Run:  python mo_pole_diagnostic.py [--domain {morocco,north_morocco}]
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import shapely
import scipy.stats as stats
import statsmodels.api as sm

HERE = Path(__file__).resolve().parent
DATA = (HERE / ".." / ".." / "Article3_data").resolve()
# mo_era5_mslp.nc (lat 28-42N, lon 10W-35E, 1981-2025) contains all four poles, so
# one file serves both variants and both baselines.
SLP_FILE = DATA / "mo_era5_mslp.nc"
DOMAINS = {
    "morocco": (Path("/Users/workad/Workspace/shapfiles/Morocco/Morocco.shp"), "output"),
    "north_morocco": (Path("/Users/workad/Workspace/shapfiles/maradm/MAR_adm/MAR_adm0.shp"),
                      "output_north_morocco"),
}
# Pole grid points per variant, recovered by exhaustive search against the index
# files. Each (variant, baseline) pair is checked against the index file it should
# reproduce (max |difference| ~ 0, asserted below) before anything is inferred.
VARIANTS = {
    "MOI2": {"west": (36.00, -5.00), "east": (32.00, 34.50)},
    "MOI1": {"west": (36.50, 3.00), "east": (30.00, 31.50)},
}
# The standardisation baseline must be the one the analysis it accompanies used,
# otherwise the seasonal aggregates differ: standardisation is per calendar month,
# so a change of baseline re-weights the three months of a season relative to each
# other and shifts seasonal correlations (noticeably in MAM). 1991-2020 matches the
# canonical headline run; 1981-2010 reproduces the earlier runs in ./output*.
BASELINES = {
    "1991-2020": {
        "clim": ("1991", "2020"),
        "index_file": {"MOI2": "MO_index_MOI2_difffirst_9120.nc",
                       "MOI1": "MO_index_MOI1_difffirst_9120.nc"},
        "outdir": {"MOI2": "output_canonical_moi2", "MOI1": "output_canonical_moi1"},
    },
    "1981-2010": {
        "clim": ("1981", "2010"),
        "index_file": {"MOI2": "MO_index.nc", "MOI1": "MO_index_MOI1_difffirst.nc"},
        "outdir": {"MOI2": "output", "MOI1": "output"},
    },
}
YEARS = range(1991, 2025)
SEASONS = {"SON": (9, 10, 11), "DJF": (12, 1, 2), "MAM": (3, 4, 5)}

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--domain", choices=sorted(DOMAINS), default="morocco")
parser.add_argument("--variant", choices=sorted(VARIANTS), default="MOI2")
parser.add_argument("--baseline", choices=sorted(BASELINES), default="1991-2020",
                    help="index standardisation baseline; must match the run this "
                         "diagnostic is reported alongside")
args = parser.parse_args()
POLE_WEST = VARIANTS[args.variant]["west"]
POLE_EAST = VARIANTS[args.variant]["east"]
_BASE = BASELINES[args.baseline]
CLIM = _BASE["clim"]
INDEX_FILE = _BASE["index_file"][args.variant]
SHAPEFILE, _ = DOMAINS[args.domain]
OUTDIR = HERE / (_BASE["outdir"][args.variant] if args.domain == "morocco"
                 else DOMAINS[args.domain][1])
OUTDIR.mkdir(exist_ok=True)


def zscore_monthly(s, clim=CLIM):
    """Standardise per calendar month over `clim` - the index's own convention."""
    c = s.loc[clim[0]:clim[1]]
    return (s - s.index.month.map(c.groupby(c.index.month).mean())) \
        / s.index.month.map(c.groupby(c.index.month).std())


def seasonal_mean(s, months):
    """DJF(Y) = Dec(Y-1)+Jan(Y)+Feb(Y); all three months required."""
    sub = s[s.index.month.isin(months)]
    key = pd.Series(sub.index.year + (sub.index.month == 12), index=sub.index)
    g = sub.groupby(key)
    out = g.mean()
    out[g.count() < 3] = np.nan
    return out


# ---------------------------------------------------------------------------
# 1. The two poles, and proof that they are the ones behind MO_index.nc
# ---------------------------------------------------------------------------
print(f"variant: {args.variant}  west={POLE_WEST} east={POLE_EAST}")
print(f"index standardisation baseline: {args.baseline}  "
      f"(verifying against {INDEX_FILE})")
print(f"output: {OUTDIR}")
slp = xr.open_dataset(SLP_FILE).drop_vars(["expver", "number"], errors="ignore")
time = pd.to_datetime(slp.valid_time.values)
west = pd.Series(slp["msl"].sel(latitude=POLE_WEST[0], longitude=POLE_WEST[1],
                                method="nearest").values, index=time)
east = pd.Series(slp["msl"].sel(latitude=POLE_EAST[0], longitude=POLE_EAST[1],
                                method="nearest").values, index=time)

mo_reconstructed = zscore_monthly(west - east)
mo_published = xr.open_dataset(DATA / INDEX_FILE)["MO_std"].to_series()
mo_published.index = mo_published.index.to_period("M").to_timestamp()
check = pd.concat({"file": mo_published, "recon": mo_reconstructed}, axis=1,
                  join="inner").dropna()
max_diff = float((check["file"] - check["recon"]).abs().max())
print(f"MO reconstruction check: n = {len(check)}, max |difference| = {max_diff:.2e}")
assert max_diff < 1e-6, "the pole pair does not reproduce MO_index.nc"

# Each pole standardised the same way, so the three series are on one scale.
z_west = zscore_monthly(west)
z_east = zscore_monthly(east)

# ---------------------------------------------------------------------------
# 2. The same precipitation predictand the main analysis uses
# ---------------------------------------------------------------------------
geom = gpd.read_file(SHAPEFILE).union_all()
precip = xr.open_dataset(DATA / "CDR-1983-2025.nc").rename({"datetime": "time"})["precip"]
precip = precip.where(precip != -99.0)                    # PERSIANN-CDR missing flag
lon2d, lat2d = np.meshgrid(precip.lon.values.astype(float), precip.lat.values.astype(float))
mask = xr.DataArray(shapely.contains_xy(geom, lon2d, lat2d), dims=("lat", "lon"),
                    coords={"lat": precip.lat, "lon": precip.lon})
weights = np.cos(np.deg2rad(precip.lat)) * mask
p_series = precip.weighted(weights).mean(dim=("lat", "lon"), skipna=True).to_series()
p_series.index = p_series.index.to_period("M").to_timestamp()
base = p_series.loc["1991":"2020"]
p_anom = p_series - p_series.index.month.map(base.groupby(base.index.month).mean())

z850 = xr.open_dataset(DATA / "Z850.nc").drop_vars(["number", "expver"], errors="ignore")
z850 = z850.rename({"valid_time": "time"})["z"].squeeze("pressure_level", drop=True) / 9.80665
lon2d, lat2d = np.meshgrid(z850.longitude.values.astype(float),
                           z850.latitude.values.astype(float))
mask_z = xr.DataArray(shapely.contains_xy(geom, lon2d, lat2d), dims=("latitude", "longitude"),
                      coords={"latitude": z850.latitude, "longitude": z850.longitude})
wz = np.cos(np.deg2rad(z850.latitude)) * mask_z
z_series = z850.weighted(wz).mean(dim=("latitude", "longitude"), skipna=True).to_series()
z_series.index = z_series.index.to_period("M").to_timestamp()
basez = z_series.loc["1991":"2020"]
z850_anom = z_series - z_series.index.month.map(basez.groupby(basez.index.month).mean())

# ---------------------------------------------------------------------------
# 3. Season by season: MO vs its two poles
# ---------------------------------------------------------------------------
rows = []
for season, months in SEASONS.items():
    frame = pd.DataFrame({
        "precip": seasonal_mean(p_anom, months),
        "MO": seasonal_mean(mo_reconstructed, months),
        "SLP_west": seasonal_mean(z_west, months),
        "SLP_east": seasonal_mean(z_east, months),
        "Z850": seasonal_mean(z850_anom, months),
    })
    frame = frame.loc[frame.index.isin(YEARS)].dropna()

    print(f"\n=== {season}  (n = {len(frame)}) ===")
    for name in ["MO", "SLP_west", "SLP_east", "Z850"]:
        r, p = stats.pearsonr(frame[name], frame["precip"])
        rows.append({"variant": args.variant, "baseline": args.baseline,
                     "domain": args.domain, "season": season, "n": len(frame),
                     "series": name, "r_vs_precip": r, "p": p, "r2": r ** 2})
        print(f"   r({name:<8}, precip) = {r:+.3f}   p = {p:.2e}   R2 = {r**2:.3f}")

    # How much of MO's own variance each pole carries, and how alike the western
    # pole and the local circulation diagnostic are.
    extra = {
        "r(MO, SLP_west)": stats.pearsonr(frame["MO"], frame["SLP_west"])[0],
        "r(MO, SLP_east)": stats.pearsonr(frame["MO"], frame["SLP_east"])[0],
        "r(SLP_west, SLP_east)": stats.pearsonr(frame["SLP_west"], frame["SLP_east"])[0],
        "r(SLP_west, Z850)": stats.pearsonr(frame["SLP_west"], frame["Z850"])[0],
    }
    for k, v in extra.items():
        print(f"   {k:<24} = {v:+.3f}")
        rows.append({"variant": args.variant, "baseline": args.baseline,
                     "domain": args.domain, "season": season, "n": len(frame),
                     "series": k, "r_vs_precip": np.nan, "p": np.nan, "r2": np.nan,
                     "value": v})

    # Does the eastern pole add anything once the western pole is in the model?
    X1 = sm.add_constant(frame[["SLP_west"]])
    X2 = sm.add_constant(frame[["SLP_west", "SLP_east"]])
    f1, f2 = sm.OLS(frame["precip"], X1).fit(), sm.OLS(frame["precip"], X2).fit()
    print(f"   R2adj: west only = {f1.rsquared_adj:.3f}, west+east = {f2.rsquared_adj:.3f}"
          f"  (east pole adds {f2.rsquared_adj - f1.rsquared_adj:+.3f});"
          f" east-pole coefficient p = {f2.pvalues['SLP_east']:.3f}")
    rows.append({"variant": args.variant, "baseline": args.baseline,
                 "domain": args.domain, "season": season, "n": len(frame),
                 "series": "R2adj_west_only", "r_vs_precip": np.nan, "p": np.nan,
                 "r2": f1.rsquared_adj})
    rows.append({"variant": args.variant, "baseline": args.baseline,
                 "domain": args.domain, "season": season, "n": len(frame),
                 "series": "R2adj_west_plus_east", "r_vs_precip": np.nan,
                 "p": f2.pvalues["SLP_east"], "r2": f2.rsquared_adj})

out = OUTDIR / f"mo_pole_diagnostic_{args.variant}.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(f"\nwritten: {out}")
