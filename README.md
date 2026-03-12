# Celestial Theodolite — RMSE Analysis

A tool for testing the shape of the Earth using star occultations behind mountain peaks.

## What It Does

When a star disappears behind a mountain peak, the exact time tells us the star's precise position. Two geometric models predict where the star *should* be at that moment:

- **Flat Earth (FE):** Simple rise-over-run angle from observer to peak — `arctan(height_difference / distance)`
- **Globe Earth (GE):** Same angle minus a curvature correction — the "inscribed angle" that accounts for Earth's curve over the observer-to-peak distance

These two predictions are different. By comparing what we actually observe against what each model predicts, we measure which one fits reality better using RMSE (Root Mean Square Error).

## How It Works

1. An observer records the precise time a star is occluded by a known mountain peak
2. The star's elevation angle at that moment is calculated from its geocentric position
3. Both models' predicted angles are computed from the peak's surveyed elevation, the observer's position, and the distance between them
4. RMSE measures how far off each model's prediction is from the observations

The mountain peak acts as a natural crosshair — its elevation is surveyed, permanent, and verifiable by anyone. This makes the experiment fully repeatable.

## The Data

Observations are stored in `celtheo_data.csv`. Each row records:
- Peak name, elevation, and distance from observer
- Star name, date, and time of occultation
- Observed elevation angles under each model (FE and GE sheets)

Currently includes **64 observations** across **14 peaks** from multiple observers.

## Error Bars — How We Got Them

The error bars on the RMSE graphs show **plus or minus 2 sigma** — that's a 95% confidence window representing how precise the measurement system is.

Here's how we compute sigma:

1. For every observation, compute the residual (how far the observation was from the prediction) for both FE and GE
2. Fit a bell curve (Gaussian) to each set of residuals — this gives us sigma_FE and sigma_GE
3. **Average the two sigmas** to get one shared sigma

We use a single shared sigma for both models on purpose. If we let each model have its own error bars, the globe model would get artificially small ones — not because it's more precise, but because its residuals happen to cluster tightly (it's consistently wrong by the same amount). That would make it *look* more certain than it actually is.

By averaging, both models are evaluated against the same measurement uncertainty window. Neither model gets an unfair advantage. The error bars reflect the **instrument's precision**, not how well either model fits.

To convert from degrees to meters (so the graphs show a physical height uncertainty at the peak):

```
error_meters = distance * tan(2 * sigma)
```

Since sigma is the same everywhere but each peak is a different distance away, the meter value changes per peak — closer peaks have smaller error bars in meters, farther peaks have larger ones.

## Usage

```bash
# Interactive menu
python3 celtheo.py

# RMSE for a specific peak with graph
python3 celtheo.py --mode rmse --peak "pikes" --graph

# RMSE for all peaks with graphs
python3 celtheo.py --mode rmse --peak all --graph

# Stats mode (standard deviation + residuals)
python3 celtheo.py --mode stats --peak all
```

## Requirements

- Python 3.8+
- matplotlib (for graph output)
