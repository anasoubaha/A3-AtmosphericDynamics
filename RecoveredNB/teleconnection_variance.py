"""
Variance in Moroccan seasonal precipitation explained by teleconnection indices,
and the increment they add beyond a regional circulation diagnostic.

Period      : 1991-2024 (34 season-years), seasons SON / DJF / MAM.
              DJF of year Y = Dec(Y-1) + Jan(Y) + Feb(Y).
Primary     : DJF (pre-specified).
Predictand  : domain-mean PERSIANN-CDR precipitation anomaly (mm/month).
Predictors  : NAO, MO, WeMO, Nino3.4 (teleconnection indices) and Z850
              (domain-mean 850 hPa geopotential height anomaly, gpm).
Domain      : selected with --domain (see DOMAINS below); each domain writes to
              its own output folder so the two runs never overwrite each other.

Run:  python teleconnection_variance.py                        -> ./output/
      python teleconnection_variance.py --domain north_morocco -> ./output_north_morocco/
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
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.outliers_influence import variance_inflation_factor

# ===========================================================================
# 0. Configuration
# ===========================================================================
HERE = Path(__file__).resolve().parent
DATA = (HERE / ".." / ".." / "Article3_data").resolve()
# Spatial domains. "morocco" is the original run; "north_morocco" repeats it
# over Morocco excluding Western Sahara, to test whether the arid, almost
# rainless south is driving the result. Everything else is held fixed, so the
# two runs differ ONLY in which grid cells enter the two domain means.
DOMAINS = {
    "morocco": {
        "name": "Morocco",
        "long_name": "Morocco including Western Sahara",
        "shapefile": Path("/Users/workad/Workspace/shapfiles/Morocco/Morocco.shp"),
        "outdir": "output",
    },
    "north_morocco": {
        "name": "northern Morocco",
        "long_name": "Morocco excluding Western Sahara (GADM MAR_adm0)",
        "shapefile": Path("/Users/workad/Workspace/shapfiles/maradm/MAR_adm/MAR_adm0.shp"),
        "outdir": "output_north_morocco",
    },
}
# Index construction variants. "published" reproduces the original run (both
# indices built as difference-then-standardise). "zfirst" swaps in the canonical
# difference-of-standardised-poles form, which is how Martin-Vide & Lopez-Bustins
# define WeMOi; see index_construction_check.md. Each writes to its own folder.
INDEX_SETS = {
    "published": {"MO": "MO_index.nc", "WeMO": "WeMO_index.nc", "suffix": ""},
    "zfirst": {"MO": "MO_index_MOI2_zfirst.nc", "WeMO": "WeMO_index_zfirst.nc",
               "suffix": "_zfirst"},
    # Each index in the form its OWN publisher defines: CRU's MOI2 is normalised
    # AFTER differencing (so MO_index.nc is already canonical), while WeMOi is the
    # difference of standardised poles. This is the build to use for the paper.
    "canonical": {"MO": "MO_index.nc", "WeMO": "WeMO_index_zfirst.nc",
                  "suffix": "_canonical"},
    # The MOI1 (Algiers-Cairo) pole pair, whose western pole sits 471 km from
    # Morocco instead of 34 km, paired with the canonical WeMO in both cases.
    "moi1_canonical": {"MO": "MO_index_MOI1_difffirst.nc", "WeMO": "WeMO_index_zfirst.nc",
                       "suffix": "_moi1_canonical"},
    "moi1_zfirst": {"MO": "MO_index_MOI1_zfirst.nc", "WeMO": "WeMO_index_zfirst.nc",
                    "suffix": "_moi1_zfirst"},
    # Headline builds: canonical construction (MO normalised AFTER differencing per
    # Palutikof; WeMO standardised per station BEFORE differencing per Martin-Vide),
    # both rebuilt on the 1991-2020 baseline that every other anomaly here uses.
    "canonical_moi2": {"MO": "MO_index_MOI2_difffirst_9120.nc",
                       "WeMO": "WeMO_index_zfirst_9120.nc", "suffix": "_canonical_moi2",
                       "mo_order": "diff_first (normalised after differencing)",
                       "wemo_order": "z_first (standardised per pole, then differenced)",
                       "baseline": "1991-2020",
                       "previous": {"MO": "MO_index.nc", "WeMO": "WeMO_index_zfirst.nc"},
                       "file_tag": "_MOI2"},
    "canonical_moi1": {"MO": "MO_index_MOI1_difffirst_9120.nc",
                       "WeMO": "WeMO_index_zfirst_9120.nc", "suffix": "_canonical_moi1",
                       "mo_order": "diff_first (normalised after differencing)",
                       "wemo_order": "z_first (standardised per pole, then differenced)",
                       "baseline": "1991-2020",
                       "previous": {"MO": "MO_index_MOI1_difffirst.nc",
                                    "WeMO": "WeMO_index_zfirst.nc"},
                       "file_tag": "_MOI1"},
}
_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument("--domain", choices=sorted(DOMAINS), default="morocco",
                     help="spatial domain for the precipitation and Z850 means")
_parser.add_argument("--indices", choices=sorted(INDEX_SETS), default="published",
                     help="which MO/WeMO construction to use")
_args = _parser.parse_args()
DOMAIN = DOMAINS[_args.domain]
INDEX_SET = INDEX_SETS[_args.indices]
SHAPEFILE = DOMAIN["shapefile"]
OUTDIR = HERE / (DOMAIN["outdir"] + INDEX_SET["suffix"])
OUTDIR.mkdir(exist_ok=True)

YEARS = range(1991, 2025)              # 34 season-years, inclusive of 2024
EXPECTED_N = len(YEARS)                # 34
CLIM = ("1991-01-01", "2020-12-31")    # anomaly baseline (WMO 1991-2020)
SEASONS = {"SON": (9, 10, 11), "DJF": (12, 1, 2), "MAM": (3, 4, 5)}
INDICES = ["NAO", "MO", "WeMO", "Nino34"]     # teleconnection indices (model M1)
PREDICTORS = INDICES + ["Z850"]               # full predictor set (model M3)
PRIMARY_SEASON = "DJF"
ALPHA = 0.05                                  # BH false-discovery rate q
G0 = 9.80665                                  # m s-2, for geopotential -> gpm
CDR_MISSING = -99.0                           # PERSIANN-CDR missing-value flag

# Season-years that CANNOT be formed because a source month is missing outright.
# Any deviation from EXPECTED_N that is not listed here aborts the run (see
# check_season_count). DJF 1992 needs Feb 1992, which is absent from
# PERSIANN-CDR over the entire grid (all cells flagged -99), as is Feb 1990.
DOCUMENTED_GAPS = {
    ("DJF", 1992): "PERSIANN-CDR 1992-02 is missing over the whole grid (-99)",
}

print("=" * 78)
print("Moroccan seasonal precipitation: teleconnection vs circulation control")
print("=" * 78)
print(f"data directory : {DATA}")
print(f"domain         : {DOMAIN['long_name']}")
print(f"index build    : {_args.indices}  (MO: {INDEX_SET['MO']}, WeMO: {INDEX_SET['WeMO']})")
if "mo_order" in INDEX_SET:
    print(f"  MO construction   : {INDEX_SET['mo_order']}")
    print(f"  WeMO construction : {INDEX_SET['wemo_order']}")
    print(f"  index baseline    : {INDEX_SET['baseline']}")
print(f"shapefile      : {SHAPEFILE}")
print(f"output         : {OUTDIR}")
print(f"period         : {min(YEARS)}-{max(YEARS)}  (expected n = {EXPECTED_N} per season)")
print(f"anomaly base   : {CLIM[0][:7]} to {CLIM[1][:7]}")


# ===========================================================================
# 1. Helpers
# ===========================================================================
def to_180(ds, lon_name):
    """Return `ds` on a -180..180 longitude axis, stating the convention found.

    Both input grids here are already -180..180, but the conversion is done
    explicitly so the code is safe against a re-download on a 0..360 axis
    (ERA5 is served in both conventions depending on the request).
    """
    lon = ds[lon_name]
    if float(lon.max()) > 180.0:
        print(f"    longitude convention: 0-360  -> converting to -180..180")
        ds = ds.assign_coords({lon_name: (((lon + 180) % 360) - 180)})
        ds = ds.sortby(lon_name)
    else:
        print(f"    longitude convention: -180..180 (no conversion needed)"
              f"  [{float(lon.min()):.2f} to {float(lon.max()):.2f}]")
    return ds


def polygon_mask(ds, lat_name, lon_name, geom):
    """Boolean (lat, lon) mask: True where the grid-cell CENTRE is inside `geom`.

    Centre-in-polygon is the standard cheap masking rule; at 0.25 deg it
    resolves the domain with several hundred cells, so fractional-area weighting of the
    coastal/border cells would shift the domain mean only marginally.
    """
    lon2d, lat2d = np.meshgrid(ds[lon_name].values.astype(float),
                               ds[lat_name].values.astype(float))
    mask = shapely.contains_xy(geom, lon2d, lat2d)
    return xr.DataArray(mask, dims=(lat_name, lon_name),
                        coords={lat_name: ds[lat_name], lon_name: ds[lon_name]})


def area_mean(da, mask, lat_name, lon_name):
    """Cosine-latitude-weighted spatial mean inside `mask`.

    cos(lat) makes each cell's weight proportional to its true area on the
    sphere. skipna=True renormalises the weights over the valid cells of each
    time step, so a partially missing field is still averaged correctly; a
    fully missing field returns NaN (weights sum to zero) and is caught later.
    """
    weights = np.cos(np.deg2rad(da[lat_name])) * mask
    return da.weighted(weights).mean(dim=(lat_name, lon_name), skipna=True)


def monthly_anomaly(series, clim=CLIM):
    """Subtract the calendar-month climatology computed over `clim`."""
    base = series.loc[clim[0]:clim[1]]
    clim_mean = base.groupby(base.index.month).mean()
    return series - series.index.month.map(clim_mean)


def seasonal_mean(series, months):
    """Mean over the three months of a season, indexed by season-year.

    December is assigned to the FOLLOWING year so that DJF(Y) = Dec(Y-1),
    Jan(Y), Feb(Y). A season-year is returned as NaN unless all three of its
    months are present and non-missing - no 2-month seasons enter the analysis.
    """
    sub = series[series.index.month.isin(months)]
    season_year = pd.Series(sub.index.year + (sub.index.month == 12), index=sub.index)
    grouped = sub.groupby(season_year)
    out = grouped.mean()
    out[grouped.count() < 3] = np.nan       # completeness requirement
    out.index.name = "season_year"
    return out


def check_season_count(season, frame):
    """Fail loudly unless n == EXPECTED_N or the shortfall is documented."""
    missing = sorted(set(YEARS) - set(frame.index))
    undocumented = [y for y in missing if (season, y) not in DOCUMENTED_GAPS]
    if undocumented:
        raise RuntimeError(
            f"{season}: n = {len(frame)} (expected {EXPECTED_N}). "
            f"Undocumented missing season-years: {undocumented}. "
            f"Inspect the inputs before trusting any result below."
        )
    for y in missing:
        print(f"    ! {season} {y} dropped - {DOCUMENTED_GAPS[(season, y)]}")
    if len(frame) != EXPECTED_N:
        print(f"    ! {season}: n = {len(frame)} (expected {EXPECTED_N}); "
              f"shortfall is documented above and carried into every result.")


def standardize(frame, columns):
    """z-score the given columns over the analysis sample (ddof=1).

    Applied to predictors only: it puts every regression coefficient in
    'mm per one standard deviation of predictor', so the coefficients are
    directly comparable. R2, F and the p-values are unaffected by this.
    """
    out = frame.copy()
    out[columns] = (frame[columns] - frame[columns].mean()) / frame[columns].std(ddof=1)
    return out


def detrend_on_year(frame, columns):
    """Remove a least-squares linear trend in TIME from each column.

    Regressing on the season-year (not on position) keeps the trend estimate
    correct when a season-year is missing from the record.
    """
    out = frame.copy()
    years = frame.index.values.astype(float)
    for c in columns:
        slope, intercept = np.polyfit(years, frame[c].values, 1)
        out[c] = frame[c].values - (slope * years + intercept)
    return out


def lag1_autocorr(series):
    """Lag-1 autocorrelation over consecutive season-years only.

    Pairs whose season-years are not one year apart (i.e. either side of a
    dropped season) are excluded rather than silently treated as adjacent.
    """
    y = series.index.values
    x = series.values
    pairs = [(x[i], x[i + 1]) for i in range(len(x) - 1) if y[i + 1] - y[i] == 1]
    a = np.array(pairs)
    r, p = stats.pearsonr(a[:, 0], a[:, 1])
    return r, p, len(pairs)


# ===========================================================================
# 2. Load inputs
# ===========================================================================
print("\n" + "-" * 78)
print("2. INPUTS")
print("-" * 78)

domain_gdf = gpd.read_file(SHAPEFILE)
assert domain_gdf.crs is not None and domain_gdf.crs.to_epsg() == 4326, \
    "shapefile must be EPSG:4326 to match the NetCDF grids"
geom = domain_gdf.union_all()
print(f"  shapefile      : {len(domain_gdf)} feature(s), EPSG:4326, "
      f"bounds lon [{domain_gdf.total_bounds[0]:.2f}, {domain_gdf.total_bounds[2]:.2f}] "
      f"lat [{domain_gdf.total_bounds[1]:.2f}, {domain_gdf.total_bounds[3]:.2f}]")

# --- 2a. ERA5 850 hPa geopotential -------------------------------------------
print("  Z850 (ERA5 850 hPa geopotential)")
z850_ds = xr.open_dataset(DATA / "Z850.nc", engine="netcdf4")
z850_ds = z850_ds.drop_vars(["number", "expver"], errors="ignore")
z850_ds = z850_ds.rename({"valid_time": "time"})
z850_ds = to_180(z850_ds, "longitude")
# ERA5 stores geopotential (m2 s-2); divide by g to get geopotential height (gpm)
z850 = z850_ds["z"].squeeze("pressure_level", drop=True) / G0
print(f"    grid {z850.sizes['latitude']} x {z850.sizes['longitude']} @0.25deg, "
      f"n_months = {z850.sizes['time']} "
      f"({str(z850.time.values[0])[:7]} to {str(z850.time.values[-1])[:7]})")

# --- 2b. PERSIANN-CDR precipitation ------------------------------------------
print("  Precipitation (PERSIANN-CDR)")
precip_ds = xr.open_dataset(DATA / "CDR-1983-2025.nc", engine="netcdf4")
precip_ds = precip_ds.rename({"datetime": "time"})
precip_ds = to_180(precip_ds, "lon")
precip = precip_ds["precip"]
# PERSIANN-CDR flags missing data as -99 rather than NaN; convert before any
# averaging, otherwise -99 is treated as a (very dry) valid observation.
n_flagged = int((precip.values == CDR_MISSING).sum())
precip = precip.where(precip != CDR_MISSING)
print(f"    grid {precip.sizes['lat']} x {precip.sizes['lon']} @0.25deg, "
      f"n_months = {precip.sizes['time']} "
      f"({str(precip.time.values[0])[:7]} to {str(precip.time.values[-1])[:7]})")
print(f"    {n_flagged:,} cell-months flagged {CDR_MISSING:g} -> set to NaN")

# --- 2c. Teleconnection indices ----------------------------------------------
print("  Teleconnection indices")
raw_indices = {
    "NAO":    ("NAO_index.nc", "NAO"),                      # NOAA CPC, standardised
    "MO":     (INDEX_SET["MO"], "MO_std"),                  # ERA5-derived, standardised
    "WeMO":   (INDEX_SET["WeMO"], "WeMO_std"),              # ERA5-derived, standardised
    "Nino34": ("NOAAPSL-HadISST-Nino34.nc", "value"),       # NOAA PSL, SST anomaly (degC)
}
monthly = {}
for name, (fname, var) in raw_indices.items():
    da = xr.open_dataset(DATA / fname, engine="netcdf4")[var]
    s = da.to_series()
    s.index = s.index.to_period("M").to_timestamp()   # snap to month start
    monthly[name] = s
    print(f"    {name:<7}: n_months = {len(s)} "
          f"({s.index[0]:%Y-%m} to {s.index[-1]:%Y-%m}), "
          f"mean = {s.mean():+.3f}, sd = {s.std():.3f}")
# NAO/MO/WeMO arrive already standardised and Nino3.4 is already an SST anomaly,
# so no further de-seasonalisation is applied to the indices themselves.


# ===========================================================================
# 3. Build the two domain-mean monthly series
# ===========================================================================
print("\n" + "-" * 78)
print(f"3. DOMAIN-MEAN SERIES - {DOMAIN['long_name']} (cos-lat weighted)")
print("-" * 78)

mask_z = polygon_mask(z850, "latitude", "longitude", geom)
mask_p = polygon_mask(precip, "lat", "lon", geom)
print(f"  cells inside the domain: Z850 {int(mask_z.sum())}, "
      f"precip {int(mask_p.sum())}")
EXPECTED_CELLS = {"morocco": 1017, "north_morocco": 639}[_args.domain]
if int(mask_p.sum()) != EXPECTED_CELLS or int(mask_z.sum()) != EXPECTED_CELLS:
    raise RuntimeError(f"domain cell count changed: expected {EXPECTED_CELLS}, got "
                       f"Z850 {int(mask_z.sum())}, precip {int(mask_p.sum())}")
print(f"  cell-count check: matches the {_args.domain} run ({EXPECTED_CELLS}) OK")
assert int(mask_z.sum()) > 0 and int(mask_p.sum()) > 0, "empty domain mask"

z_series = area_mean(z850, mask_z, "latitude", "longitude").to_series()
p_series = area_mean(precip, mask_p, "lat", "lon").to_series()
z_series.index = z_series.index.to_period("M").to_timestamp()
p_series.index = p_series.index.to_period("M").to_timestamp()

# Anomalies against the 1991-2020 calendar-month climatology.
monthly["Z850"] = monthly_anomaly(z_series)
precip_anom = monthly_anomaly(p_series)
print(f"  Z850   monthly anomaly: n = {monthly['Z850'].notna().sum()}, "
      f"sd = {monthly['Z850'].std():.2f} gpm")
print(f"  precip monthly anomaly: n = {precip_anom.notna().sum()} valid of "
      f"{len(precip_anom)}, sd = {precip_anom.std():.2f} mm/month")


# ===========================================================================
# 4. Seasonal aggregation
# ===========================================================================
print("\n" + "-" * 78)
print("4. SEASONAL MEANS")
print("-" * 78)

panels = {}
for season, months in SEASONS.items():
    cols = {name: seasonal_mean(s, months) for name, s in monthly.items()}
    cols["precip"] = seasonal_mean(precip_anom, months)
    frame = pd.DataFrame(cols)
    frame = frame.loc[frame.index.isin(YEARS)].dropna()      # complete cases only
    print(f"  {season}: n = {len(frame)} season-years "
          f"({frame.index.min()}-{frame.index.max()})")
    check_season_count(season, frame)
    assert not frame.isna().any().any(), f"{season}: NaN survived into the panel"
    panels[season] = frame


# ===========================================================================
# 5. Correlations, FDR control, and the index correlation matrix
# ===========================================================================
print("\n" + "-" * 78)
print("5. CORRELATIONS")
print("-" * 78)

corr_rows = []
matrix_rows = []
for season, frame in panels.items():
    for predictor in PREDICTORS:
        r, p = stats.pearsonr(frame[predictor], frame["precip"])
        corr_rows.append({"season": season, "predictor": predictor, "n": len(frame),
                          "r": r, "r2": r ** 2, "p_raw": p})
    # 4x4 matrix among the teleconnection indices (collinearity of M1)
    cmat = frame[INDICES].corr()
    for a in INDICES:
        for b in INDICES:
            matrix_rows.append({"season": season, "index_1": a, "index_2": b,
                                "r": cmat.loc[a, b]})

corr = pd.DataFrame(corr_rows)
# Benjamini-Hochberg across the whole family of predictor-precipitation tests
# (5 predictors x 3 seasons = 15). The per-season version (5 tests) is reported
# alongside it because the DJF analysis is pre-specified and stands on its own.
corr["p_bh_all"] = multipletests(corr["p_raw"], alpha=ALPHA, method="fdr_bh")[1]
corr["significant_bh_all"] = corr["p_bh_all"] < ALPHA
corr["p_bh_within_season"] = corr.groupby("season")["p_raw"].transform(
    lambda p: multipletests(p, alpha=ALPHA, method="fdr_bh")[1])
corr["significant_bh_within_season"] = corr["p_bh_within_season"] < ALPHA

for season in SEASONS:
    print(f"  {season} (vs precipitation anomaly)")
    for _, row in corr[corr.season == season].iterrows():
        flag = "*" if row.significant_bh_all else " "
        print(f"    {row.predictor:<7} r = {row.r:+.3f}  p = {row.p_raw:.4f}  "
              f"p_BH(15) = {row.p_bh_all:.4f} {flag}")
print(f"  (* = passes Benjamini-Hochberg at q = {ALPHA} across all 15 tests)")

# ---------------------------------------------------------------------------
# 5b. Baseline check: same construction order, previous reference period
# ---------------------------------------------------------------------------
# Changing the reference period rescales each index. Because the standardisation
# is per CALENDAR MONTH, it is not a single global rescaling: the three months of
# a season can be re-weighted slightly relative to one another, so a seasonal-mean
# correlation may shift in the third decimal. A shift larger than that, or any
# shift in a sign, would mean something other than rescaling changed.
if "previous" in INDEX_SET:
    print("\n  baseline check vs the previous reference period "
          f"(same construction order):")
    worst = 0.0
    for idx_name, prev_file in INDEX_SET["previous"].items():
        prev = xr.open_dataset(DATA / prev_file, engine="netcdf4")[
            "MO_std" if idx_name == "MO" else "WeMO_std"].to_series()
        prev.index = prev.index.to_period("M").to_timestamp()
        for season, months in SEASONS.items():
            new_r = float(corr[(corr.season == season)
                               & (corr.predictor == idx_name)]["r"].iloc[0])
            both = pd.concat({"i": seasonal_mean(prev, months),
                              "p": panels[season]["precip"]}, axis=1).dropna()
            old_r = stats.pearsonr(both["i"], both["p"])[0]
            delta = new_r - old_r
            worst = max(worst, abs(delta))
            print(f"    {idx_name:<5} {season}: r(new {INDEX_SET['baseline']}) = {new_r:+.4f}  "
                  f"r(previous, {prev_file}) = {old_r:+.4f}  delta = {delta:+.4f}")
    print(f"    largest |delta r| = {worst:.4f}")

    # A seasonal-mean correlation can move even when nothing is wrong, because the
    # standardisation is per CALENDAR MONTH. The test that separates arithmetic from
    # a bug is done one calendar month at a time: if the new index is an affine map
    # of the old within a month, the within-month correlation with precipitation is
    # mathematically invariant, and any seasonal-mean shift is only the three months
    # being re-weighted relative to each other.
    print("\n  within-calendar-month test (the diagnostic that separates rescaling "
          "from a bug):")
    for idx_name, prev_file in INDEX_SET["previous"].items():
        var = "MO_std" if idx_name == "MO" else "WeMO_std"
        prev = xr.open_dataset(DATA / prev_file, engine="netcdf4")[var].to_series()
        prev.index = prev.index.to_period("M").to_timestamp()
        new_idx = monthly[idx_name]
        joined = pd.concat({"old": prev, "new": new_idx, "precip": precip_anom},
                           axis=1, join="inner").dropna()
        joined = joined.loc[f"{min(YEARS)}":f"{max(YEARS)}"]
        for month in SEASONS[PRIMARY_SEASON]:
            k = joined[joined.index.month == month]
            slope, intercept = np.polyfit(k["old"], k["new"], 1)
            resid = float(np.abs(k["new"] - (slope * k["old"] + intercept)).max())
            r_old = stats.pearsonr(k["old"], k["precip"])[0]
            r_new = stats.pearsonr(k["new"], k["precip"])[0]
            verdict = ("affine -> r invariant" if resid < 1e-5
                       else "NOT affine (poles rescaled separately)")
            print(f"    {idx_name:<5} month {month:>2}: scale = {slope:.4f}, "
                  f"max residual = {resid:.1e}, r_old = {r_old:+.6f}, "
                  f"r_new = {r_new:+.6f}, delta = {r_new - r_old:+.1e}  [{verdict}]")
    print("    A diff_first index is affine within a month, so its within-month r is")
    print("    unchanged and only the seasonal re-weighting moves; a z_first index is")
    print("    rebuilt pole by pole, which is more than a rescaling. Neither is a bug.")
    print("    Coefficients in mm/sigma move with the index sd by design.")


# ===========================================================================
# 6. Regression models, variance increments, VIF, persistence
# ===========================================================================
print("\n" + "-" * 78)
print("6. REGRESSION MODELS")
print("-" * 78)

# Predictor blocks. MO sits in the CIRCULATION block, not the teleconnection
# block: the pole diagnostic showed it reproduces its own western pole (r ~ 0.93 in
# both variants) and that pole correlates 0.84-0.89 with Morocco-mean Z850, so it is
# a western-Mediterranean SLP diagnostic rather than a remote index. The previous
# grouping is kept verbatim as T_conv / C_conv for one comparison sentence.
TELECONNECTIONS = ["NAO", "WeMO", "Nino34"]
CIRCULATION = ["Z850", "MO"]
MODELS = {
    "T": TELECONNECTIONS,          # teleconnection indices
    "C": CIRCULATION,              # regional circulation diagnostics
    "F": PREDICTORS,               # full
    "T_conv": INDICES,             # previous grouping: NAO+MO+WeMO+Nino34
    "C_conv": ["Z850"],            # previous grouping: Z850 alone
}
assert set(TELECONNECTIONS) | set(CIRCULATION) == set(PREDICTORS), "blocks must partition F"
assert not (set(TELECONNECTIONS) & set(CIRCULATION)), "blocks must not overlap"


def fit(frame, predictors):
    """OLS of the precipitation anomaly on standardised predictors."""
    X = sm.add_constant(frame[predictors], has_constant="add")
    return sm.OLS(frame["precip"], X).fit()


results_rows, coef_rows, vif_rows, diag_rows = [], [], [], []

for season, raw_frame in panels.items():
    print(f"\n  === {season} (n = {len(raw_frame)}) ===")

    # lag-1 persistence of the predictand (affects the effective sample size)
    r1, p1, npairs = lag1_autocorr(raw_frame["precip"])
    print(f"    lag-1 autocorrelation of precipitation: r = {r1:+.3f} "
          f"(p = {p1:.3f}, {npairs} consecutive pairs)")
    diag_rows.append({"season": season, "n": len(raw_frame), "lag1_r": r1,
                      "lag1_p": p1, "lag1_n_pairs": npairs})

    for detrended in (False, True):
        frame = detrend_on_year(raw_frame, raw_frame.columns) if detrended else raw_frame
        frame = standardize(frame, PREDICTORS)
        assert not frame.isna().any().any(), "NaN entering the regression"

        fits = {}
        for model_name, predictors in MODELS.items():
            res = fit(frame, predictors)
            fits[model_name] = res
            results_rows.append({
                "season": season, "model": model_name,
                "predictors": "+".join(predictors),
                "detrended": detrended, "n": int(res.nobs),
                "df_model": int(res.df_model), "df_resid": int(res.df_resid),
                "r2": res.rsquared, "r2_adj": res.rsquared_adj,
                "f_stat": res.fvalue, "f_pvalue": res.f_pvalue,
            })
            for term in res.params.index:
                coef_rows.append({
                    "season": season, "model": model_name, "detrended": detrended,
                    "term": term, "coef_mm_per_sd": res.params[term],
                    "std_err": res.bse[term], "t": res.tvalues[term],
                    "p": res.pvalues[term],
                    "ci_low": res.conf_int().loc[term, 0],
                    "ci_high": res.conf_int().loc[term, 1],
                })

        # Variance increments on ADJUSTED R2 (the blocks differ in predictor count).
        # Two groupings, deliberately named so they cannot be confused:
        #   NEW  (MO counted as circulation):
        #     dR2_teleconnections = F - C       what NAO+WeMO+Nino34 add over Z850+MO
        #     dR2_circulation     = F - T       what Z850+MO add over NAO+WeMO+Nino34
        #   CONVENTIONAL (MO counted as an index), for comparison only:
        #     dR2_indices_conv      = F - C_conv    four indices over Z850 alone
        #     dR2_circulation_conv  = F - T_conv    Z850 over the four indices
        r2adj = {k: v.rsquared_adj for k, v in fits.items()}
        increments = {
            "dR2_teleconnections": r2adj["F"] - r2adj["C"],
            "dR2_circulation": r2adj["F"] - r2adj["T"],
            "dR2_indices_conv": r2adj["F"] - r2adj["C_conv"],
            "dR2_circulation_conv": r2adj["F"] - r2adj["T_conv"],
        }
        for row in results_rows:
            if (row["season"] == season and row["detrended"] == detrended
                    and row["model"] == "F"):
                row.update(increments)

        tag = "detrended" if detrended else "raw      "
        print(f"    [{tag}] R2adj: T = {r2adj['T']:+.3f} | C = {r2adj['C']:+.3f} | "
              f"F = {r2adj['F']:+.3f} | T_conv = {r2adj['T_conv']:+.3f} | "
              f"C_conv = {r2adj['C_conv']:+.3f}")
        print(f"    {'':<11} new grouping : dR2_teleconnections = "
              f"{increments['dR2_teleconnections']:+.3f}, "
              f"dR2_circulation = {increments['dR2_circulation']:+.3f}")
        print(f"    {'':<11} conventional : dR2_indices_conv = "
              f"{increments['dR2_indices_conv']:+.3f}, "
              f"dR2_circulation_conv = {increments['dR2_circulation_conv']:+.3f}")

        # --- sanity checks, printed rather than silently asserted ---
        n_used = {k: int(v.nobs) for k, v in fits.items()}
        print(f"    {'':<11} n per model: {n_used}")
        assert len(set(n_used.values())) == 1, f"models disagree on n: {n_used}"
        # F is one fit reused by both groupings; this records that fact in the log
        print(f"    {'':<11} R2adj(F) used by both groupings = {r2adj['F']:.6f} "
              f"(single fit, identical by construction)")
        raw_r2 = {k: v.rsquared for k, v in fits.items()}
        for block in ("T", "C"):
            if raw_r2["F"] + 1e-12 < raw_r2[block]:
                raise RuntimeError(f"{season}: raw R2(F) {raw_r2['F']:.4f} < "
                                   f"R2({block}) {raw_r2[block]:.4f} - blocks misconfigured")
        for block in ("T", "C"):
            if r2adj["F"] < r2adj[block]:
                print(f"    {'':<11} NOTE: adjusted R2(F) {r2adj['F']:.3f} < "
                      f"R2adj({block}) {r2adj[block]:.3f} - legitimate, the extra "
                      f"predictors do not pay for their degrees of freedom")

        # VIF for every predictor of M3 (computed on the fitted design matrix,
        # constant included; VIF > 5 is the usual flag for collinearity)
        if not detrended:
            X = sm.add_constant(frame[PREDICTORS], has_constant="add")
            for i, term in enumerate(X.columns):
                if term == "const":
                    continue
                vif_rows.append({"season": season, "predictor": term,
                                 "vif": variance_inflation_factor(X.values, i)})
            season_vifs = {r["predictor"]: r["vif"] for r in vif_rows
                           if r["season"] == season}
            carrier = max(season_vifs, key=season_vifs.get)
            print(f"    max VIF in F = {season_vifs[carrier]:.2f} ({carrier})")

results = pd.DataFrame(results_rows)
coefficients = pd.DataFrame(coef_rows)
vif = pd.DataFrame(vif_rows)
diagnostics = pd.DataFrame(diag_rows)
index_matrix = pd.DataFrame(matrix_rows).pivot_table(
    index=["season", "index_1"], columns="index_2", values="r").reindex(columns=INDICES)


# ===========================================================================
# 7. Write outputs
# ===========================================================================
print("\n" + "-" * 78)
print("7. OUTPUT FILES")
print("-" * 78)
# The MOI1 and MOI2 runs are supplied side by side, so their filenames carry the
# variant; builds that are not variant-specific keep their bare names.
FILE_TAG = INDEX_SET.get("file_tag", "")
for stem, obj, index in [
    ("results_table", results, False),
    ("correlations", corr, False),
    ("index_correlation_matrix", index_matrix, True),
    ("vif", vif, False),
    ("coefficients", coefficients, False),
    ("diagnostics", diagnostics, False),
]:
    name = f"{stem}{FILE_TAG}.csv"
    obj.to_csv(OUTDIR / name, index=index)
    print(f"  {name:<36} {len(obj)} rows")


# ===========================================================================
# 8. Plain-language summary of the primary (DJF) analysis
# ===========================================================================
season = PRIMARY_SEASON
frame = panels[season]
res = results[(results.season == season) & (~results.detrended)].set_index("model")
res_dt = results[(results.season == season) & (results.detrended)].set_index("model")
cdf = corr[corr.season == season].set_index("predictor")
cf = coefficients[(coefficients.season == season) & (~coefficients.detrended)
                  & (coefficients.model == "F")].set_index("term").drop("const")
sig = cf[cf["p"] < 0.05]
season_vif = vif[vif.season == season].set_index("predictor")["vif"]
n = len(frame)
r1, _, _ = lag1_autocorr(frame["precip"])
inc = res.loc["F"]

print("\n" + "=" * 78)
print(f"PLAIN-LANGUAGE SUMMARY - {DOMAIN['name']}, {season} "
      f"{min(YEARS)}-{max(YEARS)} (n = {n})")
print("=" * 78)
print(f"""
Over {min(YEARS)}-{max(YEARS)} the {season} {DOMAIN['name']}-mean precipitation anomaly comprises
{n} season-years. MO is treated here as a circulation diagnostic alongside Z850,
not as a teleconnection index.

Blocks (adjusted R2):
  T      NAO+WeMO+Nino34          {res.loc['T', 'r2']:.3f} raw, {res.loc['T', 'r2_adj']:.3f} adjusted,
         F({int(res.loc['T', 'df_model'])},{int(res.loc['T', 'df_resid'])}) = {res.loc['T', 'f_stat']:.2f}, p = {res.loc['T', 'f_pvalue']:.4f}
  C      Z850+MO                  {res.loc['C', 'r2']:.3f} raw, {res.loc['C', 'r2_adj']:.3f} adjusted,
         F({int(res.loc['C', 'df_model'])},{int(res.loc['C', 'df_resid'])}) = {res.loc['C', 'f_stat']:.2f}, p = {res.loc['C', 'f_pvalue']:.4f}
  F      all five                 {res.loc['F', 'r2']:.3f} raw, {res.loc['F', 'r2_adj']:.3f} adjusted,
         F({int(res.loc['F', 'df_model'])},{int(res.loc['F', 'df_resid'])}) = {res.loc['F', 'f_stat']:.2f}, p = {res.loc['F', 'f_pvalue']:.4f}
  T_conv NAO+MO+WeMO+Nino34       {res.loc['T_conv', 'r2_adj']:.3f} adjusted
  C_conv Z850                     {res.loc['C_conv', 'r2_adj']:.3f} adjusted

With MO counted as circulation, the three teleconnection indices add
{inc['dR2_teleconnections']:+.3f} adjusted R2 on top of Z850+MO, while the circulation pair adds
{inc['dR2_circulation']:+.3f} on top of the indices. Under the previous grouping the same full
model gives dR2_indices_conv = {inc['dR2_indices_conv']:+.3f} and dR2_circulation_conv = {inc['dR2_circulation_conv']:+.3f}.

Predictors with p < 0.05 in F: {', '.join(f"{t} ({cf.loc[t, 'coef_mm_per_sd']:+.2f} +/- {cf.loc[t, 'std_err']:.2f} mm/month per sd)" for t in sig.index) or 'none'}.
Largest standardised coefficient: {cf['coef_mm_per_sd'].abs().idxmax()}. Maximum VIF in F is
{season_vif.max():.2f}, carried by {season_vif.idxmax()}.

After Benjamini-Hochberg control at q = {ALPHA} across all 15 predictor-precipitation
tests, the predictors that remain significant are:
{', '.join(cdf[cdf.significant_bh_all].index.tolist()) or 'none'}.

Robustness: after linear detrending the same blocks give adjusted R2 of
{res_dt.loc['T', 'r2_adj']:.3f} (T), {res_dt.loc['C', 'r2_adj']:.3f} (C) and {res_dt.loc['F', 'r2_adj']:.3f} (F), with dR2_teleconnections =
{res_dt.loc['F', 'dR2_teleconnections']:+.3f} and dR2_circulation = {res_dt.loc['F', 'dR2_circulation']:+.3f}. The predictand shows lag-1
autocorrelation of {r1:+.2f} from one {season} to the next, so the {n} season-years
can be treated as {'effectively independent' if abs(r1) < 0.3 else 'only partly independent'} for inference.
""")
print("=" * 78)
