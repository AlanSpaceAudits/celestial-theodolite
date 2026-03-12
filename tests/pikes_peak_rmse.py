#!/usr/bin/env python3
"""
Celestial Theodolite — Pikes Peak RMSE Analysis
================================================

RMSE (Root Mean Square Error) measures how far observations are from
a model's predicted angle, capturing BOTH:
  - systematic offset (bias) — is the model consistently too high/low?
  - scatter — how spread out are the observations?

This is the key metric for comparing models. Std dev only measures
scatter (it subtracts the mean, zeroing out the offset). RMSE subtracts
the PREDICTED value, so a model that consistently misses by 11 arcminutes
gets penalized for that miss.

                 ┌─────────────────────────────┐
                 │  Σ (observed_i - predicted)² │
    RMSE =      │  ─────────────────────────── │
                 │            n                 │
                 └─────────────────────────────┘

    Where:
        observed_i = each star's measured EL at time of occlusion
        predicted  = the model's geometric prediction (θ_FE or θ_GE)
        n          = number of observations

    Steps:
        1. For each observation, compute the residual: observed - predicted
        2. Square each residual
        3. Sum all squared residuals
        4. Divide by n to get the Mean Squared Error (MSE)
        5. Take the square root to get RMSE

A model with small RMSE means its geometric prediction closely
matches what was actually observed at the time of occlusion.

Usage:
    python3 pikes_peak_rmse.py
"""

import math
import statistics


# =============================================================================
#  GRAPH TOGGLE — set to 1 to show the summary graph, 0 to skip
# =============================================================================

SHOW_GRAPH = 1


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

# FE (red sheet) — observed EL using flat earth geometry
el_fe = [
    2.571259, 2.576620, 2.563023, 2.572243, 2.576189,
    2.565196, 2.557369, 2.537649, 2.551488, 2.568583,
    2.555965, 2.566379, 2.542652, 2.574791, 2.562355,
    2.548747, 2.562789,
]

# GE (green sheet) — observed EL using globe earth geometry
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
    """
    sign = "-" if deg < 0 else ""
    deg = abs(deg)
    d = int(deg)
    rem = (deg - d) * 60
    m = int(rem)
    s = (rem - m) * 60
    return f"{sign}{d}°{m}'{s:.1f}\""


def deg_to_dms_signed(deg):
    """
    Same as deg_to_dms but always shows a +/- sign.
    Used for residuals where the direction matters.
    """
    sign = "-" if deg < 0 else "+"
    deg = abs(deg)
    d = int(deg)
    rem = (deg - d) * 60
    m = int(rem)
    s = (rem - m) * 60
    return f"{sign}{d}°{m}'{s:.1f}\""


def run_rmse(label, observed, predicted, stars):
    """
    Calculate RMSE step-by-step, printing every intermediate value.

    Steps:
        1. Compute residual for each star: observed_i - predicted
        2. Square each residual: (observed_i - predicted)²
        3. Sum all squared residuals (SS)
        4. Divide SS by n to get Mean Squared Error (MSE)
        5. Take √MSE to get RMSE

    Parameters:
        label (str):      descriptive name for the dataset
        observed (list):  list of observed EL values (decimal degrees)
        predicted (float): model's predicted angle (decimal degrees)
        stars (list):     star names for labeling

    Returns:
        float: RMSE in decimal degrees
    """
    n = len(observed)

    # ── Step 1 & 2: Residuals and squared residuals ──
    # residual = how far the observation is from the prediction
    # squaring ensures all values are positive and penalizes
    # large misses more than small ones
    residuals  = [obs - predicted for obs in observed]
    sq_resids  = [res ** 2 for res in residuals]

    # ── Step 3: Sum of squared residuals ──
    ss = sum(sq_resids)

    # ── Step 4: Mean Squared Error ──
    # Divide by n (not n-1). RMSE uses n because we're measuring
    # distance from a known predicted value, not estimating a
    # population parameter from the sample.
    mse = ss / n

    # ── Step 5: Root Mean Square Error ──
    rmse = math.sqrt(mse)

    # ── Print everything ──
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"  Predicted: {predicted:.6f}°  ({deg_to_dms(predicted)})")
    print(f"{'='*70}")

    print(f"\n  Step 1 & 2 — Residuals and squared residuals:")
    print(f"  {'Star':<28s}  {'Observed':>12s}  {'Residual':>12s}  {'Residual²':>14s}")
    print(f"  {'─'*28}  {'─'*12}  {'─'*12}  {'─'*14}")
    for star, obs, res, sq in zip(stars, observed, residuals, sq_resids):
        print(f"  {star:<28s}  {obs:>11.6f}°  {res:>+11.6f}°  {sq:>14.10f}")

    print(f"\n  Step 3 — Sum of squared residuals:")
    print(f"    SS  = {ss:.10f}")

    print(f"\n  Step 4 — Mean Squared Error (SS / n):")
    print(f"    MSE = {ss:.10f} / {n} = {mse:.10f}")

    print(f"\n  Step 5 — Root Mean Square Error (√MSE):")
    print(f"    RMSE = √{mse:.10f} = {rmse:.6f}°  ({deg_to_dms(rmse)})")
    print(f"{'='*70}")

    return rmse


# =============================================================================
#  MAIN
# =============================================================================

if __name__ == "__main__":

    stars = [obs[3] for obs in observations]
    dates = sorted(set(obs[1] for obs in observations))

    # ── Print peak parameters and predicted angles ──
    print(f"\n{'='*70}")
    print(f"  {PEAK_NAME} — RMSE Analysis")
    print(f"{'='*70}")
    print(f"  Observer      = {OBSERVER}")
    print(f"  h_peak        = {h_peak} m")
    print(f"  h_obs         = {h_obs} m")
    print(f"  Δh_peak       = {delta_peak} m")
    print(f"  distance      = {d} m")
    print(f"  r (earth)     = {r} m")
    print(f"  Observations  = {len(observations)} stars over {len(dates)} days")
    print(f"  Dates         = {', '.join(dates)}")
    print(f"{'='*70}")

    print(f"\n{'='*70}")
    print(f"  Predicted Angles")
    print(f"{'='*70}")
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
    print(f"{'='*70}")

    # ── Run RMSE for each model ──

    fe_rmse = run_rmse(
        "FE (Red Sheet) — RMSE vs predicted θ_FE",
        el_fe, theta_fe, stars,
    )

    ge_rmse = run_rmse(
        "GE (Green Sheet) — RMSE vs predicted θ_GE",
        el_ge, theta_ge, stars,
    )

    # ── Summary ──
    print(f"\n{'='*70}")
    print(f"  Summary — RMSE Comparison")
    print(f"{'='*70}")
    print(f"  FE predicted  = {theta_fe:.6f}°  ({deg_to_dms(theta_fe)})")
    print(f"  FE obs mean   = {statistics.mean(el_fe):.6f}°  ({deg_to_dms(statistics.mean(el_fe))})")
    print(f"  FE RMSE       = {fe_rmse:.6f}°  ({deg_to_dms(fe_rmse)})")
    print(f"")
    print(f"  GE predicted  = {theta_ge:.6f}°  ({deg_to_dms(theta_ge)})")
    print(f"  GE obs mean   = {statistics.mean(el_ge):.6f}°  ({deg_to_dms(statistics.mean(el_ge))})")
    print(f"  GE RMSE       = {ge_rmse:.6f}°  ({deg_to_dms(ge_rmse)})")
    print(f"")
    print(f"  GE/FE ratio   = {ge_rmse / fe_rmse:.1f}x")
    print(f"{'='*70}")

    # =================================================================
    #  GRAPH — RMSE Summary (optional)
    # =================================================================
    #
    # Side-by-side bar chart comparing FE and GE:
    #   - Bar height = observed mean EL
    #   - Dashed line = model's predicted angle
    #   - Error bars = ±RMSE (how far observations sit from prediction)
    #   - Annotations for the delta, inscribed angle, and sigma assessment

    if SHOW_GRAPH:
        import matplotlib.pyplot as plt

        # ── Computed values for the graph ──
        fe_obs_mean = statistics.mean(el_fe)
        ge_obs_mean = statistics.mean(el_ge)
        fe_resid    = fe_obs_mean - theta_fe         # signed offset from prediction
        ge_resid    = ge_obs_mean - theta_ge         # signed offset from prediction
        # Global σ from Gaussian fit across all peaks (star-intersections analysis)
        # These represent the overall measurement uncertainty of the system,
        # not just per-peak scatter. Derived from 29 observations via MLE.
        fe_sd       = 0.085992                       # global FE σ (degrees)
        ge_sd       = 0.090953                       # global GE σ (degrees)

        # ±2σ converted to meters at the observer-to-peak distance
        fe_2sig_m   = d * math.tan(math.radians(2 * fe_sd))
        ge_2sig_m   = d * math.tan(math.radians(2 * ge_sd))

        # ── Create figure with two subplots ──
        fig, (ax_fe, ax_ge) = plt.subplots(1, 2, figsize=(14, 7))
        fig.suptitle(
            f"{PEAK_NAME} — RMSE Summary\n"
            f"{len(observations)} stars, {len(dates)} days  |  "
            f"Observer: {OBSERVER}",
            fontsize=15, fontweight="bold",
        )

        # ── Shared y-axis range ──
        # Both panels use the same scale for fair visual comparison.
        y_max = max(theta_fe, fe_obs_mean, theta_ge, ge_obs_mean) * 1.35
        for ax in (ax_fe, ax_ge):
            ax.set_ylim(0, y_max)
            ax.set_ylabel("Angle (°)", fontsize=12)

        # ────────────────────────────────────────────────────
        #  LEFT PANEL — Flat Earth
        # ────────────────────────────────────────────────────

        ax_fe.set_title(f"Flat Earth  (error bars = ±2σ = ±{fe_2sig_m:.0f} m)", fontsize=12)

        # Bar: predicted θ_FE (light green)
        ax_fe.bar(0, theta_fe, width=0.5, color="#a8e6a0", edgecolor="none",
                  label="FE predicted θ", zorder=2)

        # Bar: observed mean (darker green), with ±2σ error bars
        ax_fe.bar(0.55, fe_obs_mean, width=0.5, color="#4caf50", edgecolor="none",
                  label="Observed mean", zorder=2)
        ax_fe.errorbar(0.55, fe_obs_mean, yerr=2 * fe_sd, fmt="none",
                       ecolor="black", capsize=6, capthick=1.5, linewidth=1.5, zorder=3)

        # Dashed line: predicted reference
        ax_fe.axhline(y=theta_fe, color="red", linestyle="--", linewidth=1.5,
                      label="predicted ref", zorder=4)

        # ── Annotations ──

        # Delta label — positioned just above the predicted ref dashed line
        ax_fe.annotate(
            f"Δ={deg_to_dms_signed(fe_resid)} ({fe_resid:+.4f}°)",
            xy=(0.275, theta_fe + y_max * 0.02), fontsize=10, color="#2e7d32",
            ha="center", va="bottom",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#2e7d32", alpha=0.9),
        )

        # RMSE label
        ax_fe.text(
            0.275, y_max * 0.92,
            f"RMSE = {deg_to_dms(fe_rmse)}",
            fontsize=12, fontweight="bold", color="white", ha="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#2e7d32", edgecolor="none"),
        )

        # Mean residual (bias)
        ax_fe.text(
            0.275, y_max * 0.84,
            f"Mean residual = {deg_to_dms_signed(fe_resid)}",
            fontsize=10, color="#2e7d32", ha="center",
        )

        ax_fe.set_xticks([0, 0.55])
        ax_fe.set_xticklabels(["Predicted\nθ_FE", "Observed\nMean"], fontsize=10)
        ax_fe.legend(loc="upper right", fontsize=9)

        # ────────────────────────────────────────────────────
        #  RIGHT PANEL — Globe Earth
        # ────────────────────────────────────────────────────

        ax_ge.set_title(f"Globe Earth  (error bars = ±2σ = ±{ge_2sig_m:.0f} m)", fontsize=12)

        # Bar: predicted θ_GE (light pink)
        ax_ge.bar(0, theta_ge, width=0.5, color="#f8bbd0", edgecolor="none",
                  label="GE predicted θ", zorder=2)

        # Bar: observed mean (magenta), with ±2σ error bars
        ax_ge.bar(0.55, ge_obs_mean, width=0.5, color="#e91e90", edgecolor="none",
                  label="Observed mean", zorder=2)
        ax_ge.errorbar(0.55, ge_obs_mean, yerr=2 * ge_sd, fmt="none",
                       ecolor="black", capsize=6, capthick=1.5, linewidth=1.5, zorder=3)

        # Dashed line: predicted reference
        ax_ge.axhline(y=theta_ge, color="red", linestyle="--", linewidth=1.5,
                      label="predicted ref", zorder=4)

        # ── Annotations ──

        # Delta label — positioned just above the predicted ref dashed line
        ax_ge.annotate(
            f"Δ={deg_to_dms_signed(ge_resid)} ({ge_resid:+.4f}°)",
            xy=(0.275, theta_ge + y_max * 0.02), fontsize=10, color="#ad1457",
            ha="center", va="bottom",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#ad1457", alpha=0.9),
        )

        # RMSE label
        ax_ge.text(
            0.275, y_max * 0.92,
            f"RMSE = {deg_to_dms(ge_rmse)}",
            fontsize=12, fontweight="bold", color="white", ha="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#c62828", edgecolor="none"),
        )

        # Mean residual (bias)
        ax_ge.text(
            0.275, y_max * 0.84,
            f"Mean residual = {deg_to_dms_signed(ge_resid)}",
            fontsize=10, color="#c62828", ha="center",
        )

        # Drop label — inscribed angle shown between the two bars
        drop_y = (theta_ge + ge_obs_mean) / 2
        drop_m = d * math.radians(inscribed_angle)
        ax_ge.text(
            0.275, drop_y,
            f"drop = {deg_to_dms(inscribed_angle)} ({drop_m:.1f} m)",
            fontsize=9, color="#555555", ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#999999", alpha=0.85),
        )

        ax_ge.set_xticks([0, 0.55])
        ax_ge.set_xticklabels(["Predicted\nθ_GE", "Observed\nMean"], fontsize=10)
        ax_ge.legend(loc="upper right", fontsize=9)

        # ── Ratio annotation at bottom ──
        fig.text(
            0.5, 0.02,
            f"GE RMSE / FE RMSE = {ge_rmse / fe_rmse:.1f}x",
            fontsize=12, fontweight="bold", ha="center", color="#333333",
        )

        plt.tight_layout(rect=[0, 0.05, 1, 0.95])
        plt.savefig("pikes_peak_rmse.png", dpi=150, bbox_inches="tight")
        print(f"\n  Graph saved: pikes_peak_rmse.png")
        plt.show()
