#!/usr/bin/env python3
"""
Celestial Theodolite — Pikes Peak Occlusion EL Statistics
=========================================================

This script analyzes the elevation angle (EL) of stars at the moment they are
occluded by Pikes Peak. Two geometric predictions are calculated from the
surveyed peak parameters, then compared against the observed EL data.

  - FE (red/flat earth):  theta_FE = arctan(delta_peak / d)
      The geometric peak elevation as a simple rise-over-run angle.
      No curvature correction applied.

  - GE (green/globe earth):  theta_GE = theta_FE - inscribed_angle
      The globe model's predicted geometric alignment at time of
      occlusion. The inscribed angle is gamma/2 (half the central
      angle), representing the angular "drop" due to earth curvature
      over the observer-to-peak distance.

Each dataset is analyzed for its own internal consistency (how tightly
the observed EL values cluster), and then compared against its model's
predicted angle to show the residuals.

Usage:
    python3 pikes_peak_stats.py
"""

import math
import statistics


# =============================================================================
#  PEAK PARAMETERS — surveyed / known values
# =============================================================================

PEAK_NAME = "Pikes Peak"         # name of the peak (used in output labels)
OBSERVER  = "Dr. Heffron"        # observer who recorded the data

h_peak    = 4300                 # peak elevation above sea level (meters)
h_obs     = 2073                 # observer elevation above sea level (meters)
d         = 50551                # observer-to-peak distance (meters)
r         = 6371000              # earth radius (meters)


# =============================================================================
#  PREDICTED ANGLES — derived from peak parameters
# =============================================================================

# delta_peak: the height difference between peak and observer (meters)
delta_peak = h_peak - h_obs

# theta_FE: flat earth predicted elevation angle (degrees)
# Simple rise-over-run — arctan(height / distance)
theta_fe = math.degrees(math.atan(delta_peak / d))

# gamma: the central angle subtended by the observer-to-peak arc (radians)
# On a sphere, the arc length d relates to the central angle by: gamma = d / r
gamma = d / r

# inscribed_angle: half the central angle (degrees)
# Per the inscribed angle theorem, the inscribed angle is always
# half the central angle it subtends. This is the angular "drop"
# the globe model predicts over the observer-to-peak distance.
inscribed_angle = math.degrees(gamma) / 2

# theta_GE: globe earth predicted elevation angle (degrees)
# The globe prediction is the flat earth angle minus the inscribed angle.
# ge = fe - inscribed_angle
theta_ge = theta_fe - inscribed_angle


# =============================================================================
#  OBSERVATION DATA
# =============================================================================

# Each observation is a tuple:
#   (row, date, time, star_name, LSTa_green, LSTa_red)
#
# - row:        spreadsheet row number (for cross-referencing)
# - date:       date of observation (MM/DD/YY)
# - time:       local time of the observed occlusion (UTC-7)
# - star_name:  catalog name of the star that was occluded
# - LSTa_green: Local Sidereal Time alpha (degrees) — from GE (green) sheet
# - LSTa_red:   Local Sidereal Time alpha (degrees) — from FE (red) sheet

observations = [
    (3,  "01/27/25", "18:43:07", "39 Aquarii",              48.864752, 48.317538),
    (4,  "01/27/25", "19:31:40", "HD 217533",               61.035482, 60.488268),
    (5,  "01/30/25", "18:31:22", "39 Aquarii (under moon)",  48.876033, 48.328819),
    (6,  "01/30/25", "18:36:02", "HD 211410",               50.045893, 49.498679),
    (7,  "01/30/25", "19:19:19", "HD 217429",               60.896352, 60.349138),
    (8,  "01/30/25", "19:21:19", "HD 217721",               61.397721, 60.850507),
    (9,  "01/30/25", "19:52:44", "HD 221665",               69.273391, 68.726177),
    (10, "01/30/25", "19:58:27", "Omega1 Aquarii",          70.706470, 70.159256),
    (11, "01/30/25", "20:08:28", "HD 223559",               73.217493, 72.670279),
    (12, "01/30/25", "20:10:51", "HD 223774",               73.814957, 73.267743),
    (13, "01/30/25", "20:45:31", "HD 2345",                 82.505351, 81.958137),
    (14, "01/31/25", "18:27:25", "39 Aquarii",              48.871437, 48.324223),
    (15, "01/31/25", "19:13:47", "TYC 5826-746-1",         60.494839, 59.947625),
    (16, "01/31/25", "19:15:57", "HD 217533",               61.037989, 60.490775),
    (17, "01/31/25", "19:54:23", "Omega1 Aquarii",          70.672628, 70.125414),
    (18, "01/31/25", "20:04:33", "HD 223559",               73.221253, 72.674039),
    (19, "01/31/25", "20:06:57", "HD 223774",               73.822896, 73.275682),
]


# =============================================================================
#  OBSERVED EL VALUES (decimal degrees)
# =============================================================================

# EL = elevation angle of the star at the moment of occlusion.
# This is the angle from the observer's horizontal to the star
# at the instant it disappears behind the peak.

# FE (red sheet) — observed EL using flat earth geometry
# Each value corresponds to one star occlusion event.
el_fe = [
    2.571259, 2.576620, 2.563023, 2.572243, 2.576189,
    2.565196, 2.557369, 2.537649, 2.551488, 2.568583,
    2.555965, 2.566379, 2.542652, 2.574791, 2.562355,
    2.548747, 2.562789,
]

# GE (green sheet) — observed EL using globe earth geometry
# ge = fe - inscribed_angle applied per observation
el_ge = [
    2.116652, 2.122010, 2.108416, 2.117636, 2.121579,
    2.110586, 2.102761, 2.083042, 2.096874, 2.113973,
    2.101355, 2.111771, 2.088042, 2.120181, 2.107747,
    2.094134, 2.108180,
]


# =============================================================================
#  HELPER FUNCTIONS
# =============================================================================

def deg_to_dms(deg):
    """
    Convert decimal degrees to a human-readable D°M'S" string.
    Handles negative values by prefixing a minus sign.

    Example:  2.522507 -> "2°31'21.0"
             -0.187850 -> "-0°11'16.3"

    Parameters:
        deg (float): Angle in decimal degrees

    Returns:
        str: Formatted string in degrees, arcminutes, arcseconds
    """
    sign = "-" if deg < 0 else ""
    deg = abs(deg)
    d = int(deg)                # whole degrees
    rem = (deg - d) * 60        # remaining fractional degrees -> arcminutes
    m = int(rem)                # whole arcminutes
    s = (rem - m) * 60          # remaining fractional arcminutes -> arcseconds
    return f"{sign}{d}°{m}'{s:.1f}\""


def deg_to_dms_signed(deg):
    """
    Same as deg_to_dms but always shows a +/- sign.
    Used for residuals where the direction matters.

    Example:  +0.039452 -> "+0°2'22.0"
              -0.187850 -> "-0°11'16.3"
    """
    sign = "-" if deg < 0 else "+"
    deg = abs(deg)
    d = int(deg)
    rem = (deg - d) * 60
    m = int(rem)
    s = (rem - m) * 60
    return f"{sign}{d}°{m}'{s:.1f}\""


def run_stats(label, values):
    """
    Calculate sample standard deviation step-by-step, printing every
    intermediate value so it can be traced and verified by hand.

    Standard deviation formula (sample, n-1):

             ┌──────────────────────────┐
             │  Σ (xi - x̄)²            │
        s =  │ ─────────────            │
             │    n - 1                 │
             └──────────────────────────┘

    Steps:
        1. Compute the mean (x̄ = sum / n)
        2. For each value, compute the squared deviation: (xi - x̄)²
        3. Sum all squared deviations (SS = sum of squares)
        4. Divide SS by (n - 1) to get the sample variance (s²)
        5. Take the square root of the variance to get std dev (s)

    Why n-1 instead of n?
        We use n-1 (Bessel's correction) because our observations
        are a *sample* of all possible occlusions, not the entire
        population. Dividing by n-1 corrects for the bias that comes
        from estimating the population mean using the sample mean.

    Parameters:
        label (str):    descriptive name for the dataset
        values (list):  list of EL values in decimal degrees

    Returns:
        float: the calculated standard deviation in decimal degrees
    """
    n = len(values)

    # ── Step 1: Mean ──
    # x̄ = (x1 + x2 + ... + xn) / n
    total = sum(values)
    mean = total / n

    # ── Step 2: Squared deviations ──
    # For each observation, compute how far it is from the mean,
    # then square it. Squaring ensures negative deviations don't
    # cancel out positive ones.
    sq_devs = [(x - mean) ** 2 for x in values]

    # ── Step 3: Sum of squared deviations (SS) ──
    # Total squared "distance" of all points from the mean.
    ss = sum(sq_devs)

    # ── Step 4: Sample variance ──
    # Divide by (n - 1) — Bessel's correction for sample data.
    variance = ss / (n - 1)

    # ── Step 5: Standard deviation ──
    # Square root brings us back to the original units (degrees).
    stddev = math.sqrt(variance)

    # ── Print everything ──
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")

    print(f"  Step 1 — Mean:")
    print(f"    sum     = {total:.6f}")
    print(f"    n       = {n}")
    print(f"    mean    = {total:.6f} / {n} = {mean:.6f}°")

    print(f"\n  Step 2 — Squared deviations (xi - mean)²:")
    for i, (x, sq) in enumerate(zip(values, sq_devs)):
        print(f"    [{i+1:2d}]  ({x:.6f} - {mean:.6f})² = {sq:.10f}")

    print(f"\n  Step 3 — Sum of squared deviations:")
    print(f"    SS = {ss:.10f}")

    print(f"\n  Step 4 — Sample variance (SS / (n-1)):")
    print(f"    s² = {ss:.10f} / {n - 1} = {variance:.10f}")

    print(f"\n  Step 5 — Standard deviation (√variance):")
    print(f"    s  = √{variance:.10f} = {stddev:.6f}°")

    # ── Convert to arcminutes and arcseconds ──
    # 1 degree = 60 arcminutes = 3600 arcseconds
    arcmin = stddev * 60
    arcsec = stddev * 3600
    print(f"    s  = {arcmin:.2f} arcmin = {arcsec:.1f} arcsec")
    print(f"{'='*60}")

    return stddev


# =============================================================================
#  MAIN
# =============================================================================

if __name__ == "__main__":

    # ── Print peak parameters and predicted angles ──
    dates = sorted(set(obs[1] for obs in observations))
    stars = [obs[3] for obs in observations]

    print(f"\n{'='*60}")
    print(f"  {PEAK_NAME} — Peak Parameters")
    print(f"{'='*60}")
    print(f"  Observer      = {OBSERVER}")
    print(f"  h_peak        = {h_peak} m")
    print(f"  h_obs         = {h_obs} m")
    print(f"  Δh_peak       = {delta_peak} m")
    print(f"  distance      = {d} m")
    print(f"  r (earth)     = {r} m")
    print(f"{'='*60}")

    print(f"\n{'='*60}")
    print(f"  {PEAK_NAME} — Predicted Angles")
    print(f"{'='*60}")
    print(f"  θ_FE            = arctan(Δpeak / d)")
    print(f"                  = arctan({delta_peak} / {d})")
    print(f"                  = {theta_fe:.6f}°  ({deg_to_dms(theta_fe)})")
    print(f"")
    print(f"  γ (central)     = d / r")
    print(f"                  = {d} / {r}")
    print(f"                  = {math.degrees(gamma):.6f}°  ({deg_to_dms(math.degrees(gamma))})")
    print(f"")
    print(f"  inscribed angle = γ / 2")
    print(f"                  = {math.degrees(gamma):.6f} / 2")
    print(f"                  = {inscribed_angle:.6f}°  ({deg_to_dms(inscribed_angle)})")
    print(f"")
    print(f"  θ_GE            = θ_FE − inscribed_angle")
    print(f"                  = {theta_fe:.6f} − {inscribed_angle:.6f}")
    print(f"                  = {theta_ge:.6f}°  ({deg_to_dms(theta_ge)})")
    print(f"{'='*60}")

    print(f"\n  {len(observations)} stars over {len(dates)} days")
    print(f"  Dates: {', '.join(dates)}")

    # ── Internal consistency (spread around own mean) ──

    sd_fe = run_stats("FE (Red Sheet) — Flat Earth EL", el_fe)
    sd_ge = run_stats("GE (Green Sheet) — Globe Earth EL", el_ge)

    # ── Cross-check against Python's statistics.stdev() ──
    sd_fe_check = statistics.stdev(el_fe)
    sd_ge_check = statistics.stdev(el_ge)

    match_fe = "MATCH" if abs(sd_fe - sd_fe_check) < 1e-10 else "MISMATCH"
    match_ge = "MATCH" if abs(sd_ge - sd_ge_check) < 1e-10 else "MISMATCH"

    print(f"\n{'='*60}")
    print(f"  Cross-check: manual vs statistics.stdev()")
    print(f"{'='*60}")
    print(f"  FE: manual = {sd_fe:.6f}  |  stdev() = {sd_fe_check:.6f}  |  {match_fe}")
    print(f"  GE: manual = {sd_ge:.6f}  |  stdev() = {sd_ge_check:.6f}  |  {match_ge}")
    print(f"{'='*60}")

    # =================================================================
    #  RESIDUALS — observed vs predicted, per star
    # =================================================================
    #
    # For each star occlusion, compute how far the observed EL is from
    # the model's predicted angle:
    #
    #   FE residual = observed_FE_EL - θ_FE  (predicted flat earth angle)
    #   GE residual = observed_GE_EL - θ_GE  (predicted globe earth angle)
    #
    # A residual near zero means the observation matches the prediction.
    # The std dev of the residuals tells us how tightly the observations
    # cluster around each model's predicted angle — the model with the
    # smaller residual std dev is the better fit.
    #
    # NOTE: because theta_FE and theta_GE are single predicted values
    # (not the sample mean), these residuals measure accuracy against
    # the model, not just internal consistency.

    # ── FE residuals: each observed FE EL minus the predicted θ_FE ──
    residuals_fe = [obs_el - theta_fe for obs_el in el_fe]

    # ── GE residuals: each observed GE EL minus the predicted θ_GE ──
    residuals_ge = [obs_el - theta_ge for obs_el in el_ge]

    # ── Print FE residuals per star ──
    print(f"\n{'='*70}")
    print(f"  FE Residuals — observed EL vs predicted θ_FE ({deg_to_dms(theta_fe)})")
    print(f"{'='*70}")
    for star, obs_el, res in zip(stars, el_fe, residuals_fe):
        print(f"    {star:<28s}  {obs_el:.6f}° − {theta_fe:.6f}° = {res:+.6f}°  ({deg_to_dms_signed(res)})")

    # ── FE residual statistics ──
    fe_res_mean = statistics.mean(residuals_fe)
    fe_res_sd   = statistics.stdev(residuals_fe)
    print(f"  {'─'*65}")
    print(f"    Mean residual  = {fe_res_mean:+.6f}°  ({deg_to_dms_signed(fe_res_mean)})")
    print(f"    Std dev        = {fe_res_sd:.6f}°  ({deg_to_dms(fe_res_sd)})")
    print(f"{'='*70}")

    # ── Print GE residuals per star ──
    print(f"\n{'='*70}")
    print(f"  GE Residuals — observed EL vs predicted θ_GE ({deg_to_dms(theta_ge)})")
    print(f"{'='*70}")
    for star, obs_el, res in zip(stars, el_ge, residuals_ge):
        print(f"    {star:<28s}  {obs_el:.6f}° − {theta_ge:.6f}° = {res:+.6f}°  ({deg_to_dms_signed(res)})")

    # ── GE residual statistics ──
    ge_res_mean = statistics.mean(residuals_ge)
    ge_res_sd   = statistics.stdev(residuals_ge)
    print(f"  {'─'*65}")
    print(f"    Mean residual  = {ge_res_mean:+.6f}°  ({deg_to_dms_signed(ge_res_mean)})")
    print(f"    Std dev        = {ge_res_sd:.6f}°  ({deg_to_dms(ge_res_sd)})")
    print(f"{'='*70}")

    # ── Summary comparison ──
    print(f"\n{'='*70}")
    print(f"  Summary — Residuals vs Predicted Angles")
    print(f"{'='*70}")
    print(f"  FE predicted  = {theta_fe:.6f}°  ({deg_to_dms(theta_fe)})")
    print(f"  FE obs mean   = {statistics.mean(el_fe):.6f}°  ({deg_to_dms(statistics.mean(el_fe))})")
    print(f"  FE mean resid = {fe_res_mean:+.6f}°  ({deg_to_dms_signed(fe_res_mean)})")
    print(f"  FE resid s    = {fe_res_sd:.6f}°  ({deg_to_dms(fe_res_sd)})")
    print(f"")
    print(f"  GE predicted  = {theta_ge:.6f}°  ({deg_to_dms(theta_ge)})")
    print(f"  GE obs mean   = {statistics.mean(el_ge):.6f}°  ({deg_to_dms(statistics.mean(el_ge))})")
    print(f"  GE mean resid = {ge_res_mean:+.6f}°  ({deg_to_dms_signed(ge_res_mean)})")
    print(f"  GE resid s    = {ge_res_sd:.6f}°  ({deg_to_dms(ge_res_sd)})")
    print(f"{'='*70}")
