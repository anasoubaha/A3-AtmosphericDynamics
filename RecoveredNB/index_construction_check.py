"""
Does the construction order of MO and WeMO inflate their link to Moroccan rainfall?

Both indices are built in Cells 5-6 as "difference the two poles in Pa, then
standardise the difference". The published definitions are not obviously the same:
Martin-Vide & Lopez-Bustins (2006) define WeMOi as the difference of the two
standardised pressures, while CRU/Palutikof describe the MOI only as a "normalized
pressure difference".

Order matters mechanically. Differencing raw pressure weights both poles equally in
hPa, so the higher-variance pole dominates; standardising first equalises them.

This script (1) measures the variance imbalance at each pole, (2) fingerprints how
each PUBLISHED series was normalised, (3) tests which construction order reproduces
each published series, and (4) reports what each order does to the correlation with
Moroccan seasonal precipitation.

Run:  python index_construction_check.py      -> output/index_construction_check.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import shapely
import scipy.stats as stats

HERE = Path(__file__).resolve().parent
DATA = (HERE / ".." / ".." / "Article3_data").resolve()
SLP_MO = Path("/Users/workad/Workspace/A4-DroughtForecasting-Preprocessing/data/era5-slp.nc")
SHAPEFILE = Path("/Users/workad/Workspace/shapfiles/Morocco/Morocco.shp")
OUTDIR = HERE / "output"
YEARS = range(1991, 2025)
SEASONS = {"SON": (9, 10, 11), "DJF": (12, 1, 2), "MAM": (3, 4, 5)}
# pole grid points as selected by Cells 5-6 (MO verified against MO_index.nc)
POLES = {
    "MO":   {"file": SLP_MO,                        "west": (36.00, -5.00), "east": (32.00, 34.50),
             "west_name": "Gibraltar", "east_name": "Lod"},
    "WeMO": {"file": DATA / "wemo_era5_mslp.nc",    "west": (36.25, -6.00), "east": (45.50, 11.50),
             "west_name": "Cadiz", "east_name": "Padua"},
}


def zmon(s, a="1981", b="2010"):
    """Standardise per calendar month over a base period."""
    c = s.loc[a:b]
    return (s - s.index.month.map(c.groupby(c.index.month).mean())) \
        / s.index.month.map(c.groupby(c.index.month).std())


def seasonal(s, months):
    sub = s[s.index.month.isin(months)]
    key = pd.Series(sub.index.year + (sub.index.month == 12), index=sub.index)
    g = sub.groupby(key)
    out = g.mean()
    out[g.count() < 3] = np.nan
    return out


def pole_series(path, lat, lon):
    d = xr.open_dataset(path).drop_vars(["expver", "number"], errors="ignore")
    return pd.Series(d["msl"].sel(latitude=lat, longitude=lon, method="nearest").values,
                     index=pd.to_datetime(d.valid_time.values))


def morocco_precip_anomaly():
    geom = gpd.read_file(SHAPEFILE).union_all()
    p = xr.open_dataset(DATA / "CDR-1983-2025.nc").rename({"datetime": "time"})["precip"]
    p = p.where(p != -99.0)
    lon2d, lat2d = np.meshgrid(p.lon.values.astype(float), p.lat.values.astype(float))
    mask = xr.DataArray(shapely.contains_xy(geom, lon2d, lat2d), dims=("lat", "lon"),
                        coords={"lat": p.lat, "lon": p.lon})
    s = p.weighted(np.cos(np.deg2rad(p.lat)) * mask).mean(dim=("lat", "lon"), skipna=True).to_series()
    s.index = s.index.to_period("M").to_timestamp()
    base = s.loc["1991":"2020"]
    return s - s.index.month.map(base.groupby(base.index.month).mean())


def published_mo():
    raw = pd.read_fwf(DATA / "moi2.output.dat", colspecs=[(0, 4), (5, 7), (8, 10), (11, 18)],
                      names=["year", "month", "day", "MO"])
    daily = pd.Series(raw["MO"].values, index=pd.to_datetime(raw[["year", "month", "day"]]))
    return daily, daily.resample("MS").mean()


def published_wemo():
    tab = pd.read_csv(DATA / "Web_WeMOi-2020.txt", sep="\t", header=None, skiprows=9,
                      usecols=range(13), names=["year"] + list(range(1, 13)), dtype=str)
    tab = tab.apply(lambda c: pd.to_numeric(c.str.replace(",", ".", regex=False), errors="coerce"))
    tab = tab.dropna(subset=["year"]).set_index("year")
    lg = tab.stack().replace(-999, np.nan).dropna()
    return pd.Series(lg.values, index=pd.to_datetime(dict(
        year=lg.index.get_level_values(0).astype(int),
        month=lg.index.get_level_values(1).astype(int), day=1))).sort_index()


rows = []
precip = morocco_precip_anomaly()

# ---------------------------------------------------------------------------
# 1. How was each PUBLISHED series normalised?
# ---------------------------------------------------------------------------
print("=" * 78)
print("1. NORMALISATION FINGERPRINT OF THE PUBLISHED SERIES")
print("=" * 78)
mo_daily, mo_monthly = published_mo()
print("published MOI2 (daily), mean and sd over candidate base periods:")
for a, b in [("1951", "1980"), ("1961", "1990"), ("1981", "2010")]:
    c = mo_daily.loc[a:b]
    flag = "  <== mean 0, sd 1: this is the base period" if abs(c.mean()) < 1e-3 and abs(c.std() - 1) < 1e-3 else ""
    print(f"   {a}-{b}: mean = {c.mean():+.4f}, sd = {c.std():.4f}{flag}")
    rows.append({"index": "MO", "check": "published_base_period", "detail": f"{a}-{b}",
                 "mean": c.mean(), "sd": c.std()})
by_month = mo_daily.groupby(mo_daily.index.month).std()
print(f"   sd by calendar month ranges {by_month.min():.2f} (Aug) to {by_month.max():.2f} (Dec)"
      f" -> normalised ONCE over the whole series, not per calendar month")
print("\n=> published MOI2 = (P_west - P_east), then normalised over 1961-1990: DIFFERENCE FIRST.")
print("   Martin-Vide & Lopez-Bustins (2006) define WeMOi as the difference of the")
print("   standardised pressures: STANDARDISE FIRST. The two indices differ by definition.")

# ---------------------------------------------------------------------------
# 2-4. Variance imbalance, agreement with the published series, precipitation
# ---------------------------------------------------------------------------
for name, cfg in POLES.items():
    west = pole_series(cfg["file"], *cfg["west"]) / 100.0      # hPa
    east = pole_series(cfg["file"], *cfg["east"]) / 100.0
    forms = {"diff_then_std": zmon(west - east),               # as built in Cells 5-6
             "z_first": zmon(west) - zmon(east),               # difference of standardised values
             "west_pole_only": zmon(west)}
    published = mo_monthly if name == "MO" else published_wemo()

    print("\n" + "=" * 78)
    print(f"{name}: west = {cfg['west_name']} {cfg['west']}, east = {cfg['east_name']} {cfg['east']}")
    print("=" * 78)

    djf = lambda s: s[s.index.month.isin([12, 1, 2])].loc["1991":"2020"]
    sw, se = djf(west).std(), djf(east).std()
    print(f"  DJF pole variability: west {sw:.2f} hPa, east {se:.2f} hPa "
          f"-> variance ratio west/east = {(sw / se) ** 2:.2f}")
    rows.append({"index": name, "check": "djf_pole_sd_hPa", "detail": "west", "sd": sw})
    rows.append({"index": name, "check": "djf_pole_sd_hPa", "detail": "east", "sd": se})

    print("  agreement with the published series (monthly overlap):")
    for form, s in forms.items():
        j = pd.concat({"p": published, "r": s}, axis=1, join="inner").dropna()
        r = stats.pearsonr(j["p"], j["r"])[0]
        r_std = stats.pearsonr(zmon(j["p"]), zmon(j["r"]))[0]
        print(f"    {form:<15} n = {len(j)}  r = {r:.4f}   r(both standardised) = {r_std:.4f}")
        rows.append({"index": name, "check": "vs_published", "detail": form, "n": len(j),
                     "r": r, "r_both_standardised": r_std})

    print("  correlation with Moroccan precipitation:")
    for season, months in SEASONS.items():
        line = f"    {season}: "
        for form, s in forms.items():
            j = pd.concat({"p": seasonal(precip, months), "i": seasonal(s, months)}, axis=1)
            j = j.loc[j.index.isin(YEARS)].dropna()
            r, p = stats.pearsonr(j["p"], j["i"])
            line += f"{form} {r:+.3f}   "
            rows.append({"index": name, "check": "vs_precip", "detail": form, "season": season,
                         "n": len(j), "r": r, "p": p})
        print(line)

    j = pd.concat({k: seasonal(v, (12, 1, 2)) for k, v in forms.items()}, axis=1)
    j = j.loc[j.index.isin(YEARS)].dropna()
    r_cur = stats.pearsonr(j["diff_then_std"], j["west_pole_only"])[0]
    r_can = stats.pearsonr(j["z_first"], j["west_pole_only"])[0]
    print(f"  DJF r(index, west pole): diff_then_std = {r_cur:+.3f}, z_first = {r_can:+.3f}")
    rows.append({"index": name, "check": "djf_r_index_westpole", "detail": "diff_then_std", "r": r_cur})
    rows.append({"index": name, "check": "djf_r_index_westpole", "detail": "z_first", "r": r_can})

pd.DataFrame(rows).to_csv(OUTDIR / "index_construction_check.csv", index=False)
print(f"\nwritten: {OUTDIR / 'index_construction_check.csv'}")
