# Direction Finding by RSSI and Phase Difference: A Simulation Study

**English** | [Türkçe](README.tr.md)

> **Simulation only.** Every number in this repository comes from simulation. No real radio measurement was made. The measurement step (two LoRa boards, RSSI against distance and antenna angle) is left for later; see [Real measurement](#real-measurement).

In search and rescue, the direction and distance of a lost person's transmitter have to be estimated from a receiver, for example one on an aircraft, and shown on a map. This project simulates and compares the basic methods:

1. **Distance from signal strength (RSSI)** with a log-distance path loss model
2. **Direction with a rotating directional antenna**: the bearing of the strongest RSSI
3. **Direction from phase difference** between two antennas (interferometer)

On top of these, RSSI is filtered with a Kalman filter, the position is found from three receivers (trilateration), and the result is shown on a map.

File, column and variable names are in Turkish: `yol_kaybi` = path loss, `yonlu_anten` = directional antenna, `faz_farki` = phase difference, `trilaterasyon` = trilateration, `saglama` = sanity checks, `ucgen` = triangle, `dogrusal` = collinear, `hata` = error.

## Common settings

| | Value |
|---|---|
| Frequency | 868 MHz (λ = 0.3456 m) |
| Transmit power | 14 dBm |
| Reference distance d0 | 1 m, RSSI(d0) = P_tx − FSPL(1 m) = −17.21 dBm |
| Random seed | 0 for every script, results are reproducible |

Units are stated in every function's docstring: distance m, power dBm, loss and noise dB, angles degrees unless written as radians.

## 1. Distance from RSSI

**Model.** `RSSI(d) = RSSI(d0) − 10 n log10(d / d0) + X`, with shadowing noise `X ~ N(0, σ)` in dB. Distance is found by inverting it: `d̂ = d0 · 10^((RSSI(d0) − RSSI) / (10 n))`.

Substituting gives `d̂ / d = 10^(X / (10 n)) = exp(s Z)` with `s = ln(10) σ / (10 n)` and `Z ~ N(0, 1)`. So:

- the error is **multiplicative and does not depend on distance**,
- `s` is the standard deviation of `ln(d̂ / d)`, not the mean relative error,
- the mean relative error is `E|d̂/d − 1| = exp(s²/2) · erf(s / √2)`,
- the median of `d̂ / d` is 1, but its mean is `exp(s²/2) > 1`: on average distance is over-estimated.

`src/yol_kaybi.py` sweeps 30 distances from 10 to 1000 m with 2000 noise draws each (`results/yol_kaybi_ozet.csv`):

| n | σ [dB] | Mean rel. error | 90% rel. error | std(ln d̂/d) simulated | s formula | Mean rel. error formula |
|---|---|---|---|---|---|---|
| 2.0 | 6 | 64.9% | 141.3% | 0.691 | 0.691 | 64.8% |
| 2.7 | 3 | 20.8% | 42.0% | 0.255 | 0.256 | 20.9% |
| 2.7 | 6 | 44.7% | 93.1% | 0.513 | 0.512 | 44.6% |
| 2.7 | 8 | 64.0% | 140.5% | 0.683 | 0.682 | 63.7% |
| 3.5 | 6 | 33.3% | 67.1% | 0.396 | 0.395 | 33.2% |

![Path loss error](results/yol_kaybi_hata.png)

Simulation and closed-form values agree to about three digits. The curves are flat in distance, as the model predicts. **RSSI ranging is coarse:** with typical values (n = 2.7, σ = 6 dB) a single packet is off by 45% on average and by more than 93% in one case out of ten. A larger n makes the error *smaller*, because the same dB noise corresponds to a smaller change in distance when the signal falls off faster.

## 2. Direction with a rotating directional antenna

**Model.** Antenna gain in dB is Gaussian: `G(θ) = −12 (θ / B)²`, B is the 3 dB beamwidth (`G(B/2) = −3 dB`). The gain is floored at −30 dB, because real antennas have side and back lobes. The antenna turns through 360° in 10° steps and one noisy RSSI is recorded per step. The true bearing is drawn uniformly at random, 3000 times per case.

Three estimators (`src/yonlu_anten.py`):

- **max RSSI:** the angle of the strongest sample. Resolution is limited to the step.
- **parabola (3 points):** the peak sample and its two neighbours. Gain in dB is exactly quadratic in angle, so this is exact without noise. It is equivalent to fitting a Gaussian in linear power.
- **Gaussian fit:** the same parabola, but least squares over all main lobe samples within about one beamwidth of the peak.

Bearing RMSE in degrees (`results/yon_hata.csv`):

| B | σ = 0 dB (max / parabola / fit) | σ = 2 dB | σ = 6 dB |
|---|---|---|---|
| 20° | 2.86 / 0.00 / 0.00 | 4.94 / 3.05 / 1.32 | 8.03 / 6.72 / 4.38 |
| 30° | 2.89 / 0.00 / 0.00 | 7.12 / 5.64 / 1.79 | 11.65 / 10.79 / 6.19 |
| 45° | 2.90 / 0.00 / 0.00 | 10.46 / 9.39 / 2.84 | 16.45 / 15.91 / 9.66 |
| 60° | 2.89 / 0.00 / 0.00 | 13.10 / 12.36 / 2.92 | 21.60 / 21.19 / 10.11 |
| 90° | 2.87 / 0.00 / 0.00 | 18.59 / 18.13 / 3.53 | 31.49 / 31.20 / 12.94 |

![Bearing error](results/yon_hata.png)

- Without noise, max RSSI gives 2.89°, which is exactly the quantisation error of a 10° step (step / √12). The two fitting methods are exact.
- With noise the 3-point parabola is only slightly better than max RSSI, because it uses just three samples. The least squares fit is clearly the best method: at B = 30° and σ = 2 dB it gives 1.8° against 7.1°.
- Wider beams give larger errors. A wide beam has a flatter peak, so the same dB noise moves the apparent peak further.

![One rotation](results/yon_ornek_tarama.png)

## 3. Direction from phase difference (interferometer)

**Model.** With antenna spacing d, the phase difference is `Δφ = 2π d sin θ / λ`, and the estimate is `θ = arcsin(λ Δφ / (2π d))`. Phase is measured in [−π, π). Phase noise σφ is added before wrapping. Small-noise error: `σθ ≈ σφ λ / (2π d cos θ)`.

**Ambiguity.** If `d > λ/2`, several angles give the same measured phase (`src/faz_farki.py`, `results/faz_belirsizlik.csv`):

| d / λ | d [cm] | Mean number of matching angles | Unambiguous sector |
|---|---|---|---|
| 0.25 | 8.6 | 1.00 | 180° |
| 0.5 | 17.3 | 1.00 | 180° |
| 1 | 34.6 | 2.00 | 60° |
| 2 | 69.1 | 4.01 | 29° (\|θ\| < 14.5°) |
| 4 | 138.3 | 8.01 | 14° |

![Phase ambiguity](results/faz_belirsizlik.png)

Without noise, d = λ/2 recovers every angle (largest error 5e-13°). With d = 2λ the principal solution is wrong for 84% of angles. For example, a source at 40° gives the same phase as −59.0°, −20.9° and 8.2°.

**Noise.** Error for σφ = 5° (`results/faz_hata.csv`, 5000 draws per angle):

| θ | d = λ/2: RMSE | median | formula | sign flipped | d = 2λ, correct k assumed: RMSE | formula |
|---|---|---|---|---|---|---|
| 0° | 1.60° | 1.12° | 1.59° | – | 0.40° | 0.40° |
| 30° | 1.82° | 1.20° | 1.84° | 0% | 0.46° | 0.46° |
| 60° | 3.24° | 2.17° | 3.18° | 0% | 0.80° | 0.80° |
| 70° | 20.7° | 3.15° | 4.65° | 1.7% | 1.16° | 1.16° |
| 75° | 51.8° | 4.20° | 6.15° | 10.8% | 1.58° | 1.54° |
| 80° | 85.4° | 6.11° | 9.17° | 28.4% | 2.70° | 2.29° |
| 85° | 112.2° | 14.1° | 18.3° | 46.4% | 3.60° | 4.57° |

![Phase error](results/faz_hata.png)

- Up to 60° the simulation matches the formula. The error grows as `1 / cos θ` towards 90°.
- **Unexpected result:** with d = λ/2, close to 90° the true phase is close to ±π. Noise then wraps it to the other end of the range, and the estimate jumps to the opposite side. At 80° this happens in 28% of trials and the RMSE explodes. The formula does not capture this.
- d = 2λ is four times more precise, but only when the ambiguity is solved (here the correct k is given to the estimator, which is not possible in practice). Precision and ambiguity have to be traded off.
- Ordinary LoRa modules cannot measure phase. This needs a coherent two-channel receiver, for example a two-channel SDR. This step is simulation only in any case.

## 4. Filtering RSSI with a Kalman filter

**Filter** (`src/kalman.py`). State: true RSSI. Model `x_k = x_(k−1) + w`, `w ~ N(0, Q)`. Measurement `z_k = x_k + v`, `v ~ N(0, R = σ²)`. Each step: `P = P + Q`, `K = P / (P + R)`, `x = x + K (z − x)`, `P = (1 − K) P`.

With Q = 0 and P0 = R the gain is `K_k = 1 / k`, so the filter equals the running mean of all packets. The code confirms this: the largest difference is 7e-14 dB and `P = R / N = 0.72` after 50 packets.

Static target at 100 m, n = 2.7, σ = 6 dB, 50 packets, 4000 runs (`results/kalman_ozet.csv`):

| Method | RSSI RMSE [dB] | Mean distance error | Median | 90% | Drop vs raw |
|---|---|---|---|---|---|
| Raw (single packet) | 5.95 | 43.5% | 32.6% | 89.1% | – |
| Moving average (10) | 1.87 | 12.8% | 10.8% | 26.2% | 70.5% |
| Kalman Q = 0 | 0.85 | 5.8% | 4.9% | 11.9% | 86.7% |
| Kalman Q = 0.1 | 1.01 | 6.8% | 5.7% | 14.1% | 84.3% |

Theory: raw 6.00 dB, moving average 6 / √10 = 1.90 dB, Kalman Q = 0 6 / √50 = 0.85 dB.

![Kalman](results/kalman.png)

The moving average stops improving after its window is full. Kalman with Q = 0 keeps improving as σ / √k. A small Q costs a little accuracy but lets the filter follow a slowly changing signal.

**Important caveat:** the simulation draws the noise independently for every packet. For a real static target, shadowing (an obstacle between transmitter and receiver) stays the same from packet to packet. Filtering can only remove the fast part (fading, measurement noise), not this bias. The real gain will be smaller than shown here.

## 5. Trilateration and map

**Method** (`src/trilaterasyon.py`). The position minimises `Σ (‖p − p_i‖ − d_i)²` (`scipy.optimize.least_squares`), started from five points. The target is at (350, 250) m, n = 2.7, σ = 6 dB, 2000 runs. Distances come from a single packet, or from 20 packets per receiver filtered with Kalman (Q = 0).

Three receiver layouts:

- triangle (0, 0), (800, 0), (400, 693)
- nearly collinear: middle receiver 40 m off the line
- collinear (0, 0), (400, 0), (800, 0)

Position error in metres (`results/trilaterasyon.csv`):

| Geometry | Packets | Median | 90% | RMSE | Estimates on the y < 0 side |
|---|---|---|---|---|---|
| Triangle | 1 | 254.5 | 649.5 | 391.4 | 15.6% |
| Triangle | 20 | 51.7 | 100.7 | 65.9 | 0.0% |
| Nearly collinear | 1 | 323.9 | 544.3 | 376.3 | 38.7% |
| Nearly collinear | 20 | 64.0 | 438.4 | 244.1 | 31.2% |
| Collinear | 1 | 379.6 | 708.3 | 459.7 | 48.2% |
| Collinear | 20 | 439.0 | 531.2 | 362.0 | 51.2% |

![Trilateration](results/trilaterasyon.png)

- Without noise, the triangle and the nearly collinear layout find the target exactly (error 1e-9 m).
- With collinear receivers the point and its mirror image across the line, (350, −250), have **exactly the same cost**. The data cannot tell them apart. The estimator then picks one of the equal-cost solutions at random. In 40 noise-free solves, 16 landed on the target and 24 on the mirror. Filtering does not help: the error stays around the 500 m mirror distance.
- Moving the middle receiver 40 m off the line is still not enough. After filtering, 31% of estimates are on the wrong side.
- In the triangle, filtering reduces the 90% error from 650 m to 101 m.

An early version of the code always kept the first solution found. Because of the order of the starting points, this was always the side where the target was, so collinear receivers looked better than the triangle. Picking at random among equal-cost solutions fixed this. The table above is from the corrected code.

**Map.** `results/harita.html` (folium) shows, for one 20-packet run per geometry: the receivers, a circle of each estimated distance, the true position, the estimate, and a circle with the simulated 90% error radius. Each geometry is a separate layer. Local metres are converted to latitude and longitude around an **arbitrary** origin (39.0° N, 35.0° E). It does not represent a real search area. The background tiles need an internet connection, so the same run is also drawn as a static image (`results/harita.png`, no background):

![Map view](results/harita.png)

In this run the triangle misses by 56 m and the target is inside the 101 m circle. Both the nearly collinear and the collinear layouts pick the mirror solution below the receivers, 442 m and 487 m away from the target. A search team sent to the centre of those circles would search on the wrong side of the line.

## Real measurement

**Not done.** Everything above is simulation only. The planned measurement:

- one LoRa board as transmitter, one as receiver, in open ground,
- at least 30 packets at 10, 20, 50 and 100 m, distances measured with a tape, saved as CSV (distance, RSSI, SNR) in `data/`,
- fit n and σ for that environment and compare with the simulation,
- with a directional antenna: record RSSI every 15° and use the estimators from section 2.

## Sanity checks

`src/saglama.py` checks the conditions from the project plan with `assert` and stops at the first failure. All 28 checks pass:

- wavelengths for 433, 868 and 121.5 MHz,
- FSPL 71.21 dB (868 MHz, 100 m) and 65.17 dB (433 MHz, 100 m),
- noise-free distance, bearing, angle of arrival and position are exact,
- the noise formulas match the simulation,
- `G(B/2) = −3 dB`,
- d = 2λ is ambiguous outside ±14.5°,
- Kalman with Q = 0 equals the running mean,
- collinear receivers produce the mirror solution,
- the same seed gives identical results.

## Limitations

- **Simulation does not replace measurement.** Multipath (reflections), irregular antenna patterns, antenna coupling, receiver calibration and hardware errors are not modelled.
- **Shadowing is drawn independently for every packet.** In reality it is correlated in time and space for a static target, so filtering helps less than shown in section 4.
- **RSSI ranging is coarse by nature.** At σ = 6 dB a single packet is off by tens of percent, and the error grows proportionally with distance.
- **RSSI bearing depends on beamwidth and reflections.** The Gaussian main lobe and the −30 dB floor are idealised.
- **The phase method is unambiguous only for d ≤ λ/2.** Even then it fails near ±90° (end-fire), where noise wraps the phase. Ordinary LoRa modules cannot measure phase.
- **Fixed model parameters.** The path loss model uses the same n and σ for every receiver. A real search area is not homogeneous.
- **Single target, 2D only.** The trilateration uses a single fixed target position and a flat 2D geometry. Receiver height and aircraft motion are not modelled.

## Reproducing

On Windows with Python 3.12:

```powershell
py -3.12 -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

Then:

```powershell
python src\yol_kaybi.py       # section 1
python src\yonlu_anten.py     # section 2
python src\faz_farki.py       # section 3
python src\kalman.py          # section 4
python src\trilaterasyon.py   # section 5 and results\harita.html
python src\saglama.py         # all sanity checks
```

Each script writes its tables and plots to `results/`. On the development machine every script took about 2 seconds, except `trilaterasyon.py` (about 70 seconds: 12,000 position solves, five least squares starts each).

## Project layout

```
src/        source code
data/       real measurement CSV files (empty: no measurement yet)
results/    tables, plots and the map
```

## License

MIT, see [LICENSE](LICENSE). The map in `results/harita.html` loads its street and satellite background tiles from Esri (attribution is shown on the map). An internet connection is needed to see the background.
