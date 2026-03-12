#!/usr/bin/env python3
"""
Celestial Theodolite — Blodgett Peak Occlusion EL Statistics
=============================================================

This script analyzes the elevation angle (EL) of stars at the moment they are
occluded by Blodgett Peak. Two geometric predictions are calculated from the
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
    python3 blodgett_peak_stats.py
"""

import math
import statistics


# =============================================================================
#  PEAK PARAMETERS — surveyed / known values
# =============================================================================

PEAK_NAME = "Blodgett Peak"      # name of the peak (used in output labels)
OBSERVER  = "Dr. Heffron"        # observer who recorded the data

h_peak    = 2870                 # peak elevation above sea level (meters)
h_obs     = 2073                 # observer elevation above sea level (meters)
d         = 35729                # observer-to-peak distance (meters)
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
    (20, "02/07/25", "19:06:56", "LP Aquarii",  65.266912, 65.676904),
    (21, "02/07/25", "19:20:24", "HD 216718",   68.642796, 69.052788),
    (22, "02/07/25", "20:52:00", "HD 2324",     91.605492, 92.015484),
    (23, "02/07/25", "21:15:52", "HD 4915",     97.588494, 97.998486),
    (24, "02/07/25", "21:28:20", "25 Ceti",    100.713694, 101.123686),
    (25, "02/07/25", "21:39:00", "HD 7449",    103.387661, 103.797653),
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
    1.205888, 1.197605, 1.196577, 1.202692, 1.204086, 1.212806,
]

# GE (green sheet) — observed EL using globe earth geometry
# ge = fe - inscribed_angle applied per observation
el_ge = [
    0.884569, 0.876289, 0.875260, 0.881374, 0.882775, 0.891488,
]


# =============================================================================
#  HELPER FUNCTIONS
# =============================================================================

def deg_to_dms(deg):
    """
    Convert decimal degrees to a human-readable D°M'S" string.
    Handles negative values by prefixing a minus sign.

    Example:  1.277571 -> "1°16'39.3"
             -0.321319 -> "-0°19'16.7"

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
