#!/usr/bin/env python3
"""
Celestial Theodolite — Master Analysis Script
===============================================

Single entry point for all Celestial Theodolite analysis. Reads observation
data from celtheo_data.csv and dynamically computes all derived values
including global σ (measurement uncertainty across all peaks).

The Celestial Theodolite records the precise time a star is occluded by a
mountain peak. Two geometric models predict where the star should be at
the moment of occlusion:

  FE (flat earth):   θ_FE = arctan(Δh / d)
      Simple rise-over-run angle from observer to peak.

  GE (globe earth):  θ_GE = θ_FE − inscribed_angle
      Same angle minus the "drop" due to earth curvature.
      The inscribed angle = γ/2, where γ = d/r is the central
      angle subtended by the observer-to-peak arc on the sphere.

Each observation records the star's elevation angle (EL) at occlusion.
The FE sheet records EL with no curvature correction; the GE sheet
applies the inscribed angle correction. By comparing observed EL values
against each model's predicted angle, we measure how well each model
fits reality.

Supports two analysis modes:
  1. RMSE  — Root Mean Square Error vs predicted angles
             Measures total miss (bias + scatter combined).
  2. Stats — Standard deviation (internal consistency) + residuals
             Measures scatter alone, then shows per-star residuals.

Data flow:
  celtheo_data.csv  →  load_data()  →  compute_predictions()  →  analysis
                                     →  compute_global_sigma()  →  error bars

Usage:
    python3 celtheo.py                              # interactive menu
    python3 celtheo.py --mode rmse --peak 1         # RMSE for peak #1
    python3 celtheo.py --mode stats --peak all       # stats for all peaks
    python3 celtheo.py --mode rmse --peak 3 --graph  # with graph

Interactive menu shorthand:
    1-3   = mode 1 (RMSE), peak 3 (Cheyenne Mountain)
    2-5   = mode 2 (Stats), peak 5 (All peaks)
    1-2g  = mode 1 (RMSE), peak 2 (Blodgett), with graph
"""

import argparse
import csv
import math
import os
import statistics
import sys


# =============================================================================
#  CONSTANTS
# =============================================================================

R_EARTH  = 6371000          # mean earth radius (meters) — used for central angle calc
CSV_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "celtheo_data.csv")


# =============================================================================
#  HELPERS — Angle Format Conversion
# =============================================================================
#
# All internal calculations use decimal degrees. These helpers convert
# to degrees-arcminutes-arcseconds (D°M'S") for human-readable output.
# 1° = 60 arcminutes = 3600 arcseconds.

def deg_to_dms(deg):
    """
    Convert decimal degrees to D°M'S" string.

    Examples:
        1.723358  → "1°43'24.1\""
        0.214164  → "0°12'51.0\""
       -0.074600  → "-0°4'28.6\""
    """
    sign = "-" if deg < 0 else ""
    deg = abs(deg)
    d = int(deg)                # whole degrees
    rem = (deg - d) * 60        # leftover → arcminutes
    m = int(rem)                # whole arcminutes
    s = (rem - m) * 60          # leftover → arcseconds
    return f"{sign}{d}°{m}'{s:.1f}\""


def deg_to_dms_signed(deg):
    """
    Same as deg_to_dms but always shows +/- sign.
    Used for residuals where the direction of the offset matters.

    Examples:
        +0.039452  → "+0°2'22.0\""
        -0.187850  → "-0°11'16.3\""
    """
    sign = "-" if deg < 0 else "+"
    deg = abs(deg)
    d = int(deg)
    rem = (deg - d) * 60
    m = int(rem)
    s = (rem - m) * 60
    return f"{sign}{d}°{m}'{s:.1f}\""


# =============================================================================
#  DATA LOADING
# =============================================================================
#
# All observation data lives in celtheo_data.csv. Each row is one star
# occlusion event. Peak parameters (h_peak, h_obs, distance) are
# denormalized into every row so the CSV is self-contained.
#
# CSV columns:
#   peak      — human-readable peak name (e.g., "Pikes Peak")
#   slug      — filename-safe identifier (e.g., "pikes_peak") for PNG output
#   observer  — who recorded the data (e.g., "Dr. Heffron")
#   h_peak    — peak elevation above sea level (meters)
#   h_obs     — observer elevation above sea level (meters)
#   distance  — observer-to-peak distance (meters)
#   row       — spreadsheet row number (for cross-referencing source data)
#   date      — date of observation (MM/DD/YY)
#   time      — local time of occlusion
#   star      — catalog name of the occluded star
#   el_fe     — observed elevation angle, FE (red) sheet (decimal degrees)
#   el_ge     — observed elevation angle, GE (green) sheet (decimal degrees)
#
# To add a new peak: just append rows to the CSV. Everything else
# (predictions, global σ, error bars) recomputes automatically.

def load_data(csv_path=CSV_FILE):
    """
    Load all observation data from CSV, grouped by peak.

    Returns:
        peaks (dict):      peak_name → {slug, observer, h_peak, h_obs, d, obs: [...]}
        peak_order (list):  peak names in the order they first appear in the CSV
    """
    peaks = {}          # peak_name → peak data dict
    peak_order = []     # preserves CSV ordering for consistent menu numbering
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["peak"]
            # First time seeing this peak — create its entry
            if name not in peaks:
                peaks[name] = {
                    "slug":     row["slug"],        # for output filenames
                    "observer": row["observer"],
                    "h_peak":   int(row["h_peak"]),
                    "h_obs":    int(row["h_obs"]),
                    "d":        int(row["distance"]),
                    "obs":      [],                 # list of observation dicts
                }
                peak_order.append(name)
            # Append this observation to the peak's list
            peaks[name]["obs"].append({
                "row":   int(row["row"]),
                "date":  row["date"],
                "time":  row["time"],
                "star":  row["star"],
                "el_fe": float(row["el_fe"]),       # FE observed EL (degrees)
                "el_ge": float(row["el_ge"]),       # GE observed EL (degrees)
                "h_peak": int(row["h_peak"]),       # per-obs peak elevation (m)
                "h_obs":  int(row["h_obs"]),        # per-obs observer elevation (m)
                "d":      int(row["distance"]),     # per-obs distance (m)
            })
    return peaks, peak_order


# =============================================================================
#  GEOMETRIC PREDICTIONS
# =============================================================================
#
# Both models start from the same surveyed inputs:
#   h_peak  = peak elevation (m above sea level)
#   h_obs   = observer elevation (m above sea level)
#   d       = observer-to-peak distance (m)
#
# FE model: the star should be at arctan(Δh/d) — simple triangle.
# GE model: the star should be lower by the inscribed angle (γ/2),
#           which is the angular "drop" the globe predicts over distance d.
#
# The two predictions are mutually exclusive. Whichever one the
# observations cluster around is the better geometric model.

def compute_predictions(h_peak, h_obs, d):
    """
    Derive the predicted occlusion angles for both models.

    Parameters:
        h_peak (int): peak elevation above sea level (meters)
        h_obs  (int): observer elevation above sea level (meters)
        d      (int): observer-to-peak distance (meters)

    Returns:
        dict with:
            delta_peak      — height difference: h_peak - h_obs (meters)
            theta_fe        — FE predicted angle (degrees)
            gamma           — central angle subtended by arc d on sphere (radians)
            inscribed_angle — γ/2, the curvature "drop" (degrees)
            theta_ge        — GE predicted angle (degrees)
    """
    # Height difference between peak and observer
    delta_peak = h_peak - h_obs

    # FE prediction: simple rise-over-run
    #   θ_FE = arctan(Δh / d)
    theta_fe = math.degrees(math.atan(delta_peak / d))

    # Central angle: the angle at earth's center subtended by the arc
    # from observer to peak.  γ = arc_length / radius = d / r
    gamma = d / R_EARTH

    # Inscribed angle: half the central angle (inscribed angle theorem).
    # This is the angular amount the globe model predicts the line of
    # sight "drops" over the distance d due to curvature.
    #   inscribed_angle = γ / 2   (converted to degrees)
    inscribed_angle = math.degrees(gamma) / 2

    # GE prediction: FE angle minus the curvature drop
    #   θ_GE = θ_FE − inscribed_angle
    theta_ge = theta_fe - inscribed_angle

    return {
        "delta_peak":      delta_peak,
        "theta_fe":        theta_fe,
        "gamma":           gamma,           # radians
        "inscribed_angle": inscribed_angle, # degrees
        "theta_ge":        theta_ge,
    }


# =============================================================================
#  GLOBAL SIGMA — Dynamic Measurement Uncertainty
# =============================================================================
#
# σ (sigma) quantifies the measurement uncertainty of the Celestial
# Theodolite instrument, pooled across all peaks and both models.
#
# Method: averaged-residual Gaussian fit.
#   1. For each observation, compute both residuals:
#        fe_resid = obs_fe − θ_FE
#        ge_resid = obs_ge − θ_GE
#   2. Average them: avg_resid = (fe_resid + ge_resid) / 2
#      This combines both models' information per observation,
#      so peaks where FE and GE residuals partially cancel
#      (e.g., FE high + GE low) don't inflate σ.
#   3. Fit a Gaussian to all averaged residuals:
#        μ = mean, σ = √(Σ(r−μ)² / N)
#
# Why averaged residuals?
#   Using separate σ_FE and σ_GE then averaging lets a model with
#   consistently large but tight residuals get artificially small σ.
#   Using raw pooled residuals lets distant peaks (with large biases
#   in both models) inflate σ for all peaks. Averaging the FE and GE
#   residuals per observation balances these effects — the two models'
#   biases partially cancel, giving a σ that reflects the instrument's
#   measurement uncertainty without being dominated by any single peak.
#
# This σ is used for the ±2σ error bars on the RMSE graphs.
# It updates automatically when new data is added to the CSV.

def compute_global_sigma(peaks):
    """
    Compute a single shared σ from averaged FE+GE residuals.

    For each observation, the FE and GE residuals are averaged together,
    then an MLE Gaussian is fit to all averaged residuals across all peaks.
    This gives a fair shared σ where both models contribute equally and
    peaks with opposing biases (FE high, GE low) don't inflate σ.

    Returns:
        sigma     (float): shared global σ (degrees)
        n         (int):   total number of observations
    """
    avg_resids = []   # will hold (fe_resid + ge_resid) / 2 for every obs

    for name, peak in peaks.items():
        for o in peak["obs"]:
            # Use each observation's own geometry for its prediction
            preds = compute_predictions(o["h_peak"], o["h_obs"], o["d"])
            fe_resid = o["el_fe"] - preds["theta_fe"]
            ge_resid = o["el_ge"] - preds["theta_ge"]
            avg_resids.append((fe_resid + ge_resid) / 2)

    n = len(avg_resids)

    # MLE Gaussian fit: μ = mean, σ = √(Σ(r−μ)² / N)
    mu = sum(avg_resids) / n
    sigma = math.sqrt(sum((r - mu)**2 for r in avg_resids) / n)
    return sigma, n


# =============================================================================
#  RMSE ANALYSIS
# =============================================================================
#
# RMSE (Root Mean Square Error) measures how far observations are from
# a model's predicted angle. Unlike standard deviation (which measures
# scatter around the mean), RMSE measures distance from the PREDICTION.
# This means it captures BOTH:
#   - systematic offset (bias): is the model consistently too high/low?
#   - scatter: how spread out are the observations?
#
#                 ┌─────────────────────────────┐
#                 │  Σ (observed_i - predicted)² │
#    RMSE =      │  ─────────────────────────── │
#                 │            n                 │
#                 └─────────────────────────────┘
#
# A model with small RMSE means its geometric prediction closely
# matches what was actually observed. This is the key comparison metric.

def run_rmse(peak_name, peak, preds, sigma, n_global, show_graph=False):
    """
    Run full RMSE analysis for a single peak: print parameters,
    predicted angles, step-by-step RMSE for both FE and GE, and summary.

    Parameters:
        peak_name  (str):   human-readable peak name
        peak       (dict):  peak data from load_data()
        preds      (dict):  output of compute_predictions()
        sigma      (float): shared global σ for error bars (degrees)
        n_global   (int):   total observations used to compute σ
        show_graph (bool):  if True, generate and display bar chart

    Returns:
        (fe_rmse, ge_rmse): RMSE values in decimal degrees
    """
    observer = peak["observer"]
    stars    = [o["star"] for o in peak["obs"]]
    dates    = sorted(set(o["date"] for o in peak["obs"]))
    el_fe    = [o["el_fe"] for o in peak["obs"]]
    el_ge    = [o["el_ge"] for o in peak["obs"]]

    # Peak-level params (from first observation) — used for derivation display
    h_peak          = peak["h_peak"]
    h_obs           = peak["h_obs"]
    d               = peak["d"]
    theta_fe        = preds["theta_fe"]
    theta_ge        = preds["theta_ge"]
    delta_peak      = preds["delta_peak"]
    gamma           = preds["gamma"]
    inscribed_angle = preds["inscribed_angle"]

    # Per-observation predictions (supports peaks with varying geometry)
    obs_preds = [compute_predictions(o["h_peak"], o["h_obs"], o["d"]) for o in peak["obs"]]
    pred_fe = [p["theta_fe"] for p in obs_preds]
    pred_ge = [p["theta_ge"] for p in obs_preds]
    inscribed_angles = [p["inscribed_angle"] for p in obs_preds]
    mean_pred_fe = statistics.mean(pred_fe)
    mean_pred_ge = statistics.mean(pred_ge)

    # ── Print peak parameters ──
    print(f"\n{'='*70}")
    print(f"  {peak_name} — RMSE Analysis")
    print(f"{'='*70}")
    print(f"  Observer      = {observer}")
    print(f"  h_peak        = {h_peak} m")
    print(f"  h_obs         = {h_obs} m")
    print(f"  Δh_peak       = {delta_peak} m")
    print(f"  distance      = {d} m")
    print(f"  r (earth)     = {R_EARTH} m")
    print(f"  Observations  = {len(peak['obs'])} stars over {len(dates)} days")
    print(f"  Dates         = {', '.join(dates)}")
    print(f"{'='*70}")

    # ── Print predicted angles (showing full derivation) ──
    print(f"\n{'='*70}")
    print(f"  Predicted Angles")
    print(f"{'='*70}")
    print(f"  θ_FE            = arctan(Δpeak / d)")
    print(f"                  = arctan({delta_peak} / {d})")
    print(f"                  = {theta_fe:.6f}°  ({deg_to_dms(theta_fe)})")
    print(f"")
    print(f"  γ (central)     = d / r")
    print(f"                  = {d} / {R_EARTH}")
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

    # ── RMSE calculations (step-by-step printed inside) ──
    fe_rmse = _calc_rmse("FE (Red Sheet) — RMSE vs predicted θ_FE", el_fe, pred_fe, stars)
    ge_rmse = _calc_rmse("GE (Green Sheet) — RMSE vs predicted θ_GE", el_ge, pred_ge, stars)

    # ── Summary comparison ──
    print(f"\n{'='*70}")
    print(f"  Summary — RMSE Comparison")
    print(f"{'='*70}")
    print(f"  FE predicted  = {mean_pred_fe:.6f}°  ({deg_to_dms(mean_pred_fe)})")
    print(f"  FE obs mean   = {statistics.mean(el_fe):.6f}°  ({deg_to_dms(statistics.mean(el_fe))})")
    print(f"  FE RMSE       = {fe_rmse:.6f}°  ({deg_to_dms(fe_rmse)})")
    print(f"")
    print(f"  GE predicted  = {mean_pred_ge:.6f}°  ({deg_to_dms(mean_pred_ge)})")
    print(f"  GE obs mean   = {statistics.mean(el_ge):.6f}°  ({deg_to_dms(statistics.mean(el_ge))})")
    print(f"  GE RMSE       = {ge_rmse:.6f}°  ({deg_to_dms(ge_rmse)})")
    print(f"")
    print(f"  GE/FE ratio   = {ge_rmse / fe_rmse:.1f}x")
    print(f"")
    print(f"  Global σ       = {sigma:.6f}°  (pooled from {n_global} obs × 2 models)")
    print(f"{'='*70}")

    if show_graph:
        _rmse_graph(peak_name, peak, pred_fe, pred_ge, inscribed_angles,
                    el_fe, el_ge, fe_rmse, ge_rmse, sigma, dates)

    return fe_rmse, ge_rmse


def _calc_rmse(label, observed, predicted, stars):
    """
    Calculate RMSE step-by-step, printing every intermediate value
    so the math can be traced and verified by hand.

    RMSE formula:
        1. residual_i = observed_i − predicted
        2. squared_i  = residual_i²
        3. SS  = Σ squared_i           (sum of squared residuals)
        4. MSE = SS / n                (mean squared error)
        5. RMSE = √MSE                (root mean square error)

    We divide by n (not n-1) because we're measuring distance from
    a KNOWN predicted value, not estimating a population parameter.

    Parameters:
        label     (str):   descriptive name for output header
        observed  (list):  observed EL values (decimal degrees)
        predicted (float): model's predicted angle (decimal degrees)
        stars     (list):  star names for per-row labeling

    Returns:
        float: RMSE in decimal degrees
    """
    n = len(observed)

    # Steps 1 & 2: compute residual and its square for each observation
    # predicted is a list (one per observation) to support varying geometry
    residuals = [obs - pred for obs, pred in zip(observed, predicted)]
    sq_resids = [res ** 2 for res in residuals]

    # Step 3: sum of squared residuals
    ss = sum(sq_resids)

    # Step 4: mean squared error
    mse = ss / n

    # Step 5: root mean square error
    rmse = math.sqrt(mse)

    # ── Print the full calculation ──
    pred_mean = statistics.mean(predicted)
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"  Predicted: {pred_mean:.6f}°  ({deg_to_dms(pred_mean)})")
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
#  RMSE GRAPH
# =============================================================================
#
# Side-by-side bar chart comparing FE and GE for a single peak:
#
#   LEFT panel (green) — Flat Earth:
#     Light bar  = predicted θ_FE
#     Dark bar   = observed mean EL
#     Dashed red = predicted reference line
#     Error bars = ±2σ (global measurement uncertainty)
#
#   RIGHT panel (pink/magenta) — Globe Earth:
#     Same layout but with θ_GE prediction.
#     Also shows the "drop" label (inscribed angle in degrees + meters).
#
# Error bars are ±2σ (95% confidence) using a single shared σ pooled
# from both models' residuals — representing instrument uncertainty.
# The angular σ is converted to meters in the subtitle using:
# Δh = d × tan(2σ), where d is the observer-to-peak distance.
# Both panels show identical error bars for a fair comparison.

def _rmse_graph(peak_name, peak, pred_fe, pred_ge, inscribed_angles,
                el_fe, el_ge, fe_rmse, ge_rmse, sigma, dates):
    """Generate RMSE summary bar chart and save as PNG."""
    import matplotlib.pyplot as plt

    observer        = peak["observer"]
    d               = statistics.mean([o["d"] for o in peak["obs"]])
    slug            = peak["slug"]          # used for output filename
    theta_fe        = statistics.mean(pred_fe)
    theta_ge        = statistics.mean(pred_ge)
    inscribed_angle = statistics.mean(inscribed_angles)
    h_peak          = peak["h_peak"]
    h_obs           = peak["h_obs"]
    delta_h         = h_peak - h_obs

    # Observed means and signed residuals (observed − predicted)
    fe_obs_mean = statistics.mean(el_fe)
    ge_obs_mean = statistics.mean(el_ge)
    fe_resid    = fe_obs_mean - theta_fe    # positive = obs above prediction
    ge_resid    = ge_obs_mean - theta_ge    # negative = obs below prediction

    # Convert ±2σ from degrees to meters at the peak distance.
    # At small angles: Δh ≈ d × tan(θ) ≈ d × θ_rad
    # This gives the physical height uncertainty the error bars represent.
    # Both models share the same σ (instrument uncertainty), so same meter value.
    sig_2_m = d * math.tan(math.radians(2 * sigma))

    # ── Create figure: two subplots side by side ──
    fig, (ax_fe, ax_ge) = plt.subplots(1, 2, figsize=(14, 7))
    fig.suptitle(
        f"{peak_name} — RMSE Summary\n"
        f"{len(peak['obs'])} stars, {len(dates)} days  |  "
        f"Observer: {observer}\n"
        f"dist = {d:,.0f} m  |  Δh = {delta_h:,} m",
        fontsize=15, fontweight="bold",
    )

    # Both panels share the same y-axis scale for fair visual comparison.
    # 1.35x headroom leaves space for annotations above the bars.
    y_max = max(theta_fe, fe_obs_mean, theta_ge, ge_obs_mean) * 1.35
    for ax in (ax_fe, ax_ge):
        ax.set_ylim(0, y_max)
        ax.set_ylabel("Angle (°)", fontsize=12)

    # ────────────────────────────────────────────────────
    #  LEFT PANEL — Flat Earth (green)
    # ────────────────────────────────────────────────────

    # Title includes ±2σ in meters for physical context
    ax_fe.set_title(f"Flat Earth  (error bars = ±2σ = ±{sig_2_m:.0f} m)", fontsize=12)

    # Light green bar: model's predicted angle θ_FE
    ax_fe.bar(0, theta_fe, width=0.5, color="#a8e6a0", edgecolor="none",
              label="FE predicted θ", zorder=2)

    # Dark green bar: observed mean EL, with ±2σ error bars
    ax_fe.bar(0.55, fe_obs_mean, width=0.5, color="#4caf50", edgecolor="none",
              label="Observed mean", zorder=2)
    ax_fe.errorbar(0.55, fe_obs_mean, yerr=2 * sigma, fmt="none",
                   ecolor="black", capsize=6, capthick=1.5, linewidth=1.5, zorder=3)

    # Dashed red line: predicted reference (makes the offset visually obvious)
    ax_fe.axhline(y=theta_fe, color="red", linestyle="--", linewidth=1.5,
                  label="predicted ref", zorder=4)

    # ── Annotations ──

    # Delta label: shows the signed offset between observed mean and prediction
    ax_fe.annotate(
        f"Δ = {deg_to_dms_signed(fe_resid)}  ({'+' if fe_resid >= 0 else '-'}{d * math.tan(math.radians(abs(fe_resid))):.1f} m)",
        xy=(0.275, theta_fe + y_max * 0.02), fontsize=10, color="#2e7d32",
        ha="center", va="bottom",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#2e7d32", alpha=0.9),
    )

    # RMSE badge: the key metric, prominently displayed
    ax_fe.text(
        0.275, y_max * 0.92,
        f"RMSE = {deg_to_dms(fe_rmse)}",
        fontsize=12, fontweight="bold", color="white", ha="center",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#2e7d32", edgecolor="none"),
    )

    ax_fe.set_xticks([0, 0.55])
    ax_fe.set_xticklabels(["Predicted\nθ_FE", "Observed\nMean"], fontsize=10)
    ax_fe.legend(loc="upper right", fontsize=9)

    # ────────────────────────────────────────────────────
    #  RIGHT PANEL — Globe Earth (pink/magenta)
    # ────────────────────────────────────────────────────

    ax_ge.set_title(f"Globe Earth  (error bars = ±2σ = ±{sig_2_m:.0f} m)", fontsize=12)

    # Solid pink bar: GE predicted angle (bottom portion, 0 to θ_GE)
    ax_ge.bar(0, theta_ge, width=0.5, color="#f8bbd0", edgecolor="none",
              label="GE predicted θ", zorder=2)

    # Hatched "drop" zone: the inscribed angle removed by the globe model
    # (from θ_GE up to θ_FE), with diagonal lines to show it's subtracted
    drop_m = d * math.radians(inscribed_angle)
    import matplotlib as mpl
    ax_ge.bar(0, inscribed_angle, width=0.5, bottom=theta_ge,
              color="#f8bbd0", edgecolor="#e91e90", linewidth=0.8,
              hatch="//", zorder=3)
    # Drop label inside the hatched zone
    drop_center_y = theta_ge + inscribed_angle / 2
    ax_ge.text(0, drop_center_y,
               f"drop = {deg_to_dms(inscribed_angle)}\n({drop_m:.1f} m)",
               fontsize=8, color="#ad1457", ha="center", va="center",
               bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#e91e90", alpha=0.85),
               zorder=5)
    mpl.rcParams['hatch.color'] = "#e91e90"

    # Magenta bar: observed mean EL, with ±2σ error bars
    ax_ge.bar(0.55, ge_obs_mean, width=0.5, color="#e91e90", edgecolor="none",
              label="Observed mean", zorder=2)
    ax_ge.errorbar(0.55, ge_obs_mean, yerr=2 * sigma, fmt="none",
                   ecolor="black", capsize=6, capthick=1.5, linewidth=1.5, zorder=3)

    # Dashed red line: predicted reference at θ_GE
    ax_ge.axhline(y=theta_ge, color="red", linestyle="--", linewidth=1.5,
                  label="predicted ref", zorder=4)

    # ── Annotations ──

    # Delta label — position above the hatched drop zone (θ_FE level)
    ax_ge.annotate(
        f"Δ = {deg_to_dms_signed(ge_resid)}  ({'+' if ge_resid >= 0 else '-'}{d * math.tan(math.radians(abs(ge_resid))):.1f} m)",
        xy=(0.275, theta_fe + y_max * 0.02), fontsize=10, color="#ad1457",
        ha="center", va="bottom",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#ad1457", alpha=0.9),
    )

    # RMSE badge
    ax_ge.text(
        0.275, y_max * 0.92,
        f"RMSE = {deg_to_dms(ge_rmse)}",
        fontsize=12, fontweight="bold", color="white", ha="center",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#c62828", edgecolor="none"),
    )

    ax_ge.set_xticks([0, 0.55])
    ax_ge.set_xticklabels(["Predicted\nθ_GE", "Observed\nMean"], fontsize=10)
    ax_ge.legend(loc="upper right", fontsize=9)

    # ── Ratio annotation at bottom center ──
    # Shows how many times worse the GE RMSE is compared to FE.
    fig.text(
        0.5, 0.02,
        f"GE RMSE / FE RMSE = {ge_rmse / fe_rmse:.1f}x",
        fontsize=12, fontweight="bold", ha="center", color="#333333",
    )

    # ── Save and display ──
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    png_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"{slug}_rmse.png")
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    print(f"\n  Graph saved: {png_path}")
    plt.show()


# =============================================================================
#  STATS ANALYSIS (Standard Deviation)
# =============================================================================
#
# This mode examines each peak's observations for:
#
# 1. Internal consistency (standard deviation):
#    How tightly do the observations cluster around their own mean?
#    This tells you about measurement precision independent of
#    whether the model prediction is right or wrong.
#
# 2. Residuals vs predicted angle:
#    For each star, how far was the observed EL from the model's
#    prediction? The mean residual shows systematic bias; the
#    residual std dev shows scatter around the prediction.
#
# Standard deviation formula (sample, n-1 / Bessel's correction):
#
#          ┌──────────────────────────┐
#          │  Σ (xi - x̄)²            │
#     s =  │ ─────────────            │
#          │    n - 1                 │
#          └──────────────────────────┘
#
# We use n-1 because our observations are a sample, not the full
# population. Bessel's correction prevents underestimating variance.

def run_stats(peak_name, peak, preds):
    """
    Run standard deviation + residual analysis for a single peak.
    Prints full step-by-step calculations for both FE and GE.

    Parameters:
        peak_name (str):   human-readable peak name
        peak      (dict):  peak data from load_data()
        preds     (dict):  output of compute_predictions()
    """
    observer = peak["observer"]
    stars    = [o["star"] for o in peak["obs"]]
    dates    = sorted(set(o["date"] for o in peak["obs"]))
    el_fe    = [o["el_fe"] for o in peak["obs"]]
    el_ge    = [o["el_ge"] for o in peak["obs"]]

    h_peak          = peak["h_peak"]
    h_obs           = peak["h_obs"]
    d               = peak["d"]
    theta_fe        = preds["theta_fe"]
    theta_ge        = preds["theta_ge"]
    delta_peak      = preds["delta_peak"]
    gamma           = preds["gamma"]
    inscribed_angle = preds["inscribed_angle"]

    # ── Peak parameters ──
    print(f"\n{'='*60}")
    print(f"  {peak_name} — Peak Parameters")
    print(f"{'='*60}")
    print(f"  Observer      = {observer}")
    print(f"  h_peak        = {h_peak} m")
    print(f"  h_obs         = {h_obs} m")
    print(f"  Δh_peak       = {delta_peak} m")
    print(f"  distance      = {d} m")
    print(f"  r (earth)     = {R_EARTH} m")
    print(f"{'='*60}")

    # ── Predicted angles (full derivation) ──
    print(f"\n{'='*60}")
    print(f"  {peak_name} — Predicted Angles")
    print(f"{'='*60}")
    print(f"  θ_FE            = arctan(Δpeak / d)")
    print(f"                  = arctan({delta_peak} / {d})")
    print(f"                  = {theta_fe:.6f}°  ({deg_to_dms(theta_fe)})")
    print(f"")
    print(f"  γ (central)     = d / r")
    print(f"                  = {d} / {R_EARTH}")
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

    print(f"\n  {len(peak['obs'])} stars over {len(dates)} days")
    print(f"  Dates: {', '.join(dates)}")

    # ── Internal consistency: how tightly do EL values cluster? ──
    sd_fe = _calc_stddev("FE (Red Sheet) — Flat Earth EL", el_fe)
    sd_ge = _calc_stddev("GE (Green Sheet) — Globe Earth EL", el_ge)

    # ── Cross-check: verify our manual calc matches Python's built-in ──
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

    # ── Residuals: each observation vs the model's prediction ──
    # residual = observed - predicted
    # Positive residual = observation is ABOVE the predicted angle
    # Negative residual = observation is BELOW the predicted angle
    # Per-observation predictions (supports peaks with varying geometry)
    obs_preds = [compute_predictions(o["h_peak"], o["h_obs"], o["d"]) for o in peak["obs"]]
    pred_fe = [p["theta_fe"] for p in obs_preds]
    pred_ge = [p["theta_ge"] for p in obs_preds]
    mean_pred_fe = statistics.mean(pred_fe)
    mean_pred_ge = statistics.mean(pred_ge)

    residuals_fe = [obs_el - pred for obs_el, pred in zip(el_fe, pred_fe)]
    residuals_ge = [obs_el - pred for obs_el, pred in zip(el_ge, pred_ge)]

    # ── FE residuals per star ──
    print(f"\n{'='*70}")
    print(f"  FE Residuals — observed EL vs predicted θ_FE ({deg_to_dms(mean_pred_fe)})")
    print(f"{'='*70}")
    for star, obs_el, pred, res in zip(stars, el_fe, pred_fe, residuals_fe):
        print(f"    {star:<28s}  {obs_el:.6f}° − {pred:.6f}° = {res:+.6f}°  ({deg_to_dms_signed(res)})")

    fe_res_mean = statistics.mean(residuals_fe)      # systematic bias
    fe_res_sd   = statistics.stdev(residuals_fe)     # scatter around prediction
    print(f"  {'─'*65}")
    print(f"    Mean residual  = {fe_res_mean:+.6f}°  ({deg_to_dms_signed(fe_res_mean)})")
    print(f"    Std dev        = {fe_res_sd:.6f}°  ({deg_to_dms(fe_res_sd)})")
    print(f"{'='*70}")

    # ── GE residuals per star ──
    print(f"\n{'='*70}")
    print(f"  GE Residuals — observed EL vs predicted θ_GE ({deg_to_dms(mean_pred_ge)})")
    print(f"{'='*70}")
    for star, obs_el, pred, res in zip(stars, el_ge, pred_ge, residuals_ge):
        print(f"    {star:<28s}  {obs_el:.6f}° − {pred:.6f}° = {res:+.6f}°  ({deg_to_dms_signed(res)})")

    ge_res_mean = statistics.mean(residuals_ge)
    ge_res_sd   = statistics.stdev(residuals_ge)
    print(f"  {'─'*65}")
    print(f"    Mean residual  = {ge_res_mean:+.6f}°  ({deg_to_dms_signed(ge_res_mean)})")
    print(f"    Std dev        = {ge_res_sd:.6f}°  ({deg_to_dms(ge_res_sd)})")
    print(f"{'='*70}")

    # ── Summary: side-by-side comparison ──
    print(f"\n{'='*70}")
    print(f"  Summary — Residuals vs Predicted Angles")
    print(f"{'='*70}")
    print(f"  FE predicted  = {mean_pred_fe:.6f}°  ({deg_to_dms(mean_pred_fe)})")
    print(f"  FE obs mean   = {statistics.mean(el_fe):.6f}°  ({deg_to_dms(statistics.mean(el_fe))})")
    print(f"  FE mean resid = {fe_res_mean:+.6f}°  ({deg_to_dms_signed(fe_res_mean)})")
    print(f"  FE resid s    = {fe_res_sd:.6f}°  ({deg_to_dms(fe_res_sd)})")
    print(f"")
    print(f"  GE predicted  = {mean_pred_ge:.6f}°  ({deg_to_dms(mean_pred_ge)})")
    print(f"  GE obs mean   = {statistics.mean(el_ge):.6f}°  ({deg_to_dms(statistics.mean(el_ge))})")
    print(f"  GE mean resid = {ge_res_mean:+.6f}°  ({deg_to_dms_signed(ge_res_mean)})")
    print(f"  GE resid s    = {ge_res_sd:.6f}°  ({deg_to_dms(ge_res_sd)})")
    print(f"{'='*70}")


def _calc_stddev(label, values):
    """
    Calculate sample standard deviation step-by-step, printing every
    intermediate value so the math can be traced and verified by hand.

    Standard deviation formula (sample, Bessel-corrected):

             ┌──────────────────────────┐
             │  Σ (xi - x̄)²            │
        s =  │ ─────────────            │
             │    n - 1                 │
             └──────────────────────────┘

    Steps:
        1. Compute the mean (x̄ = Σx / n)
        2. For each value, compute the squared deviation: (xi - x̄)²
        3. Sum all squared deviations (SS)
        4. Divide SS by (n-1) to get sample variance (s²)
        5. Take √variance to get std dev (s)

    Why n-1?
        Bessel's correction. Our observations are a sample, not the
        entire population. Using n-1 corrects for the bias introduced
        by estimating the population mean from the sample mean.

    Parameters:
        label  (str):  descriptive name for the dataset
        values (list): observed EL values in decimal degrees

    Returns:
        float: sample standard deviation in decimal degrees
    """
    n = len(values)

    # Step 1: Mean
    total = sum(values)
    mean = total / n

    # Step 2: Squared deviations from the mean
    # Squaring ensures negatives don't cancel positives
    sq_devs = [(x - mean) ** 2 for x in values]

    # Step 3: Sum of squared deviations
    ss = sum(sq_devs)

    # Step 4: Sample variance (Bessel-corrected)
    variance = ss / (n - 1)

    # Step 5: Standard deviation
    stddev = math.sqrt(variance)

    # ── Print the full calculation ──
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

    # Also show in arcminutes and arcseconds for intuition
    arcmin = stddev * 60
    arcsec = stddev * 3600
    print(f"    s  = {arcmin:.2f} arcmin = {arcsec:.1f} arcsec")
    print(f"{'='*60}")

    return stddev


# =============================================================================
#  INTERACTIVE MENU
# =============================================================================
#
# When run with no arguments, the script shows a numbered menu.
# Input format: [mode]-[peak] with optional 'g' suffix for graph.
#
# Examples:
#   1-1   → RMSE for peak #1
#   2-3   → Stats for peak #3
#   1-5g  → RMSE for all peaks, with graph
#   q     → quit

def interactive_menu(peaks, peak_order):
    """
    Display the interactive menu and parse the user's selection.

    Returns:
        mode       (str):  "rmse" or "stats"
        selected   (list): list of peak names to analyze
        show_graph (bool): whether to generate graphs
    """
    total_obs = sum(len(p["obs"]) for p in peaks.values())

    print(f"\n{'='*60}")
    print(f"  Celestial Theodolite — Master Analysis")
    print(f"{'='*60}")
    print(f"  Data: {total_obs} observations across {len(peaks)} peaks")
    print(f"{'='*60}")

    # ── Show available modes ──
    print(f"\n  Analysis mode:")
    print(f"    1. RMSE (Root Mean Square Error)")
    print(f"    2. Stats (Standard Deviation)")

    # ── Show available peaks with observation counts ──
    print(f"\n  Select peak:")
    for i, name in enumerate(peak_order, 1):
        n_obs  = len(peaks[name]["obs"])
        n_days = len(set(o["date"] for o in peaks[name]["obs"]))
        print(f"    {i}. {name:<20s}  ({n_obs} obs, {n_days} day{'s' if n_days > 1 else ''})")
    print(f"    {len(peak_order) + 1}. All peaks")

    # ── Input format instructions ──
    print(f"\n  Enter as [mode]-[peak]  (e.g., 1-3 = RMSE for peak #3)")
    print(f"  Add 'g' for graph       (e.g., 1-3g = RMSE + graph)")

    while True:
        choice = input("\n  > ").strip().lower()

        # Allow quitting
        if choice in ("q", "quit", "exit"):
            sys.exit(0)

        # Check for trailing 'g' (graph flag)
        show_graph = choice.endswith("g")
        if show_graph:
            choice = choice[:-1]    # strip the 'g' before parsing numbers

        # Parse "mode-peak" format
        parts = choice.split("-")
        if len(parts) == 2:
            try:
                mode_idx = int(parts[0])    # 1 = RMSE, 2 = Stats
                peak_idx = int(parts[1])    # 1-N = specific peak, N+1 = all

                if mode_idx not in (1, 2):
                    print("  Mode must be 1 or 2.")
                    continue

                mode = "rmse" if mode_idx == 1 else "stats"

                if 1 <= peak_idx <= len(peak_order):
                    selected = [peak_order[peak_idx - 1]]
                elif peak_idx == len(peak_order) + 1:
                    selected = list(peak_order)
                else:
                    print(f"  Peak must be 1-{len(peak_order) + 1}.")
                    continue

                return mode, selected, show_graph

            except ValueError:
                pass

        print(f"  Format: [mode]-[peak]  e.g., 1-2 or 2-5g")

    return mode, selected, show_graph


# =============================================================================
#  MAIN — Entry Point
# =============================================================================
#
# Execution flow:
#   1. Load all observation data from CSV
#   2. Compute global σ from ALL observations (for error bars)
#   3. Get user selection (CLI args or interactive menu)
#   4. For each selected peak:
#      a. Compute predicted angles from peak parameters
#      b. Run the selected analysis mode (RMSE or Stats)
#      c. Optionally generate and display the graph

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Celestial Theodolite — Master Analysis",
        epilog="Run with no arguments for interactive menu.",
    )
    parser.add_argument("--mode", choices=["rmse", "stats"],
                        help="Analysis mode: rmse or stats")
    parser.add_argument("--peak", type=str,
                        help="Peak number (1-N) or 'all'")
    parser.add_argument("--graph", action="store_true",
                        help="Show summary graph (RMSE mode)")
    parser.add_argument("--csv", type=str, default=CSV_FILE,
                        help="Path to data CSV")
    args = parser.parse_args()

    # ── Step 1: Load all data from CSV ──
    peaks, peak_order = load_data(args.csv)

    # ── Step 2: Compute global σ dynamically from ALL data ──
    # This updates automatically when new rows are added to the CSV.
    sigma, n_global = compute_global_sigma(peaks)

    # ── Step 3: Get user selection ──
    if args.mode and args.peak:
        # CLI mode: arguments provided directly
        mode = args.mode
        show_graph = args.graph
        if args.peak.lower() == "all":
            selected = list(peak_order)
        else:
            try:
                idx = int(args.peak)
                selected = [peak_order[idx - 1]]
            except (ValueError, IndexError):
                # Fall back to substring matching (e.g., --peak "pikes")
                match = [n for n in peak_order if args.peak.lower() in n.lower()]
                if match:
                    selected = match
                else:
                    print(f"  Unknown peak: {args.peak}")
                    print(f"  Available: {', '.join(f'{i+1}={n}' for i, n in enumerate(peak_order))}")
                    sys.exit(1)
    else:
        # Interactive mode: show menu and get selection
        mode, selected, show_graph = interactive_menu(peaks, peak_order)

    # ── Print global σ (computed from all peaks combined) ──
    print(f"\n  Global σ (pooled from {n_global} obs × 2 models across {len(peaks)} peaks):")
    print(f"    σ = {sigma:.6f}°  ({deg_to_dms(sigma)})")

    # ── Step 4: Run analysis for each selected peak ──
    for peak_name in selected:
        peak  = peaks[peak_name]
        preds = compute_predictions(peak["h_peak"], peak["h_obs"], peak["d"])

        if mode == "rmse":
            run_rmse(peak_name, peak, preds, sigma, n_global, show_graph)
        else:
            run_stats(peak_name, peak, preds)
