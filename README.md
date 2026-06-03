# FP2DG: Fractional-Polynomial Powered Descent Guidance

MATLAB implementation of the **FP2DG** algorithm for lunar powered descent guidance. The algorithm parameterizes the commanded thrust acceleration as a two-term fractional polynomial in time-to-go, then uses constrained optimization (`fmincon`, SQP) to find the best parameter triplet for a given mission.

The accompanying papers describe the theory in full:

- **"Optimizing Fractional-Polynomial Powered Descent Guidance Laws"** — covers the optimization-based approach implemented here
- **"Theory of Fractional-Polynomial Powered Descent Guidance"** — covers the analytical foundations of the guidance law family

---

## Table of Contents

1. [Background: The Guidance Law](#1-background-the-guidance-law)
2. [Repository Layout](#2-repository-layout)
3. [Quick Start: Running `userTest.m`](#3-quick-start-running-usertestm)
4. [Input Structs Reference](#4-input-structs-reference)
5. [Entry Points](#5-entry-points)
6. [Core Library](#6-core-library)
7. [Coordinate Functions](#7-coordinate-functions)
8. [Parameter Sweeps](#8-parameter-sweeps)
9. [Dispersion and Monte Carlo Studies](#9-dispersion-and-monte-carlo-studies)
10. [DIDO Optimal Control (Optional)](#10-dido-optimal-control-optional)
11. [Non-Dimensionalization](#11-non-dimensionalization)
12. [Constraints](#12-constraints)
13. [Two-Phase ODE Integration](#13-two-phase-ode-integration)
14. [File Call Graph](#14-file-call-graph)
15. [Output Variables Guide](#15-output-variables-guide)

---

## 1. Background: The Guidance Law

FP2DG commands thrust acceleration as:

```
a_T(t_go) = a_f* + c1 * t_go^gamma1 + c2 * t_go^gamma2
```

where `t_go` is time-to-go (counts **down** from `tgo` to 0 at landing), `a_f*` is the desired final acceleration (typically a hover-thrust-minus-gravity value), and `gamma1`, `gamma2` are fractional exponents — the free parameters.

**Coefficients `c1`, `c2`** are solved analytically from the two terminal boundary conditions — final position `r_f` and final velocity `v_f` — by inverting a 2×2 linear system. This is handled by `calculateCoeffs.m`.

**`kr`** is a composite parameter derived from the exponents:
```
kr = (gamma2 + 2)(gamma1 + 2)
```
It appears directly in the **closed-loop (state-feedback) form** of the guidance law used in simulation:
```
a_T = gamma1*(kr/(2*gamma1+4) - 1)*a_f*
    + (gamma1*kr/(2*gamma1+4) - gamma1 - 1)*g
    + ((gamma1+1)/t_go)*(1 - kr/(gamma1+2))*(v_f* - v)
    + (kr/t_go^2)*(r_f* - r - v*t_go)
```
This tracking form feeds back the current position `r` and velocity `v`, making the guidance law robust to perturbations without needing re-optimization.

**The optimization** searches over `[gamma1, gamma2, tgo]` to minimize a weighted cost:
```
J = beta * integral(|a_T| dt)  +  (1 - beta) * integral(|a_T|^2 dt)
```
- `beta = 1.0` — pure fuel-optimal (minimizes delta-V)
- `beta = 0.0` — minimum control effort (smoothest throttle profile)
- `0 < beta < 1` — a trade-off between the two

Optional constraints enforced during optimization: thrust saturation bounds, a glideslope cone, and a thrust-pointing angle limit.

---

## 2. Repository Layout

```
Thesis Repo/
|
|-- userTest.m                  <- Entry point: single run or beta sweep
|
|-- getParams.m                 <- Core pipeline: setup -> optimize -> simulate -> plot
|-- getParamsDIVERT.m           <- Variant of getParams for divert scenarios
|-- optimizationLoop.m          <- fmincon wrapper (SQP), enforces all constraints
|-- closedLoopSim.m             <- Static closed-loop simulation (no re-opt)
|-- simReOpt.m                  <- Simulation with periodic re-optimization + divert
|-- calculateCoeffs.m           <- Solves for guidance coefficients c1, c2
|-- plotting.m                  <- All figures for a single run
|-- simpsonComp13Integral.m     <- Composite Simpson 1/3 rule numerical integrator
|-- enuBasis.m                  <- Computes ENU unit vectors at a lat/lon point
|
|-- comparisonFigures.m         <- Utility: merge PNGs side-by-side for thesis figures
|-- comparisonFiguresVert.m     <- Same, vertical stacking
|
|-- CoordinateFunctions/
|   |-- PDI2MCMF.m              <- PDI state (lat/lon/alt/vel/FPA/heading) -> MCMF
|   |-- ENU2MCMF.m              <- ENU vector/position -> MCMF
|   `-- MCMF2ENU.m              <- MCMF vector/position -> ENU
|
|-- ParameterSweeps/
|   |-- gammaSweep.m            <- Entry point: 2D cost map over [gamma1, gamma2] at fixed tgo
|   `-- tgoSweep.m              <- Entry point: 1D cost curve over tgo at fixed [gamma1, gamma2]
|
|-- Dispersion and MonteCarlo/
|   |-- dispersionStudy.m       <- Entry point: full state dispersion study (parfor)
|   |-- accelMonteCarlo.m       <- Entry point: accelerometer scale factor Monte Carlo
|   |-- statsPlotting.m         <- Statistics plots called by dispersionStudy
|   `-- Seeds/
|       |-- alt_seeds.dat       <- Pre-generated standard-normal random samples
|       |-- lon_seeds.dat
|       |-- lat_seeds.dat
|       |-- V_seeds.dat
|       |-- fpa_seeds.dat
|       |-- azmth_seeds.dat
|       |-- mass_seeds.dat
|       `-- accel_seeds.dat
|
`-- DIDO/                       <- Optional: requires DIDO Toolbox (commercial)
    |-- FP2DG_problem.m         <- Entry point: run DIDO for comparison
    |-- FP2DG_cost.m
    |-- FP2DG_dynamics.m
    |-- FP2DG_events.m
    `-- FP2DG_path.m
```

---

## 3. Quick Start: Running `userTest.m`

`userTest.m` is the primary interactive entry point. Open it in MATLAB and run from the repo root (it uses `addpath` to find all dependencies automatically).

**Step 1 — Choose a beta value (or a vector for comparison sweeps):**
```matlab
betaVec = [0.7];           % single run
% betaVec = 1.0:-0.1:0.0; % sweep beta from 1 to 0, all overlaid on same figures
```

**Step 2 — Choose which features to enable:**
```matlab
glideSlopeEnabled     = false;  % enforce glideslope cone constraint
pointingEnabled       = false;  % enforce thrust-pointing angle constraint
reOptimizationEnabled = false;  % periodically re-optimize mid-flight
divertEnabled         = false;  % enable target-switch divert capability
```

> **Note:** Setting `divertEnabled = true` automatically disables glideslope and pointing constraints and enables re-optimization. The script switches to calling `getParamsDIVERT` instead of `getParams`.

**Step 3 — Run.** The script prints optimization results to the console and opens figures 1–15. If running a beta sweep, all runs overlay on the same figures in different colors.

**Step 4 — Read the summary table.** After all beta values are run, a formatted summary table is printed to the console showing gamma1, gamma2, kr, tgo, fuel costs, simulation error, and exit flags.

---

## 4. Input Structs Reference

These four structs are the primary inputs to `getParams` and appear in identical form in all five entry points.

### `PDIState` — Initial spacecraft state at PDI (Powered Descent Initiation)

| Field | Units | Description |
|---|---|---|
| `altitude` | m | Altitude above the lunar surface |
| `lonInitDeg` | deg | Longitude of the PDI location |
| `latInitDeg` | deg | Latitude of the PDI location |
| `inertialVelocity` | m/s | Inertial speed magnitude at PDI |
| `flightPathAngleDeg` | deg | Flight path angle (positive = above horizontal) |
| `azimuth` | deg | Heading angle, clockwise from North |

### `planetaryParams` — Moon parameters

| Field | Units | Description |
|---|---|---|
| `rPlanet` | m | Mean lunar radius (used for gravity and coordinate transforms) |
| `gPlanet` | m/s^2 | Lunar surface gravity |
| `gEarth` | m/s^2 | Earth gravity — used only for Isp unit conversion (Isp * g_Earth = exhaust velocity) |

### `vehicleParams` — Spacecraft parameters

| Field | Units | Description |
|---|---|---|
| `massInit` | kg | Initial (wet) mass at PDI |
| `dryMass` | kg | Dry mass (structure + payload, no propellant) |
| `isp` | s | Specific impulse |
| `maxThrust` | N | Maximum thrust magnitude |
| `minThrust` | N | Minimum thrust magnitude (engine-on lower bound) |

### `targetState` — Landing target and divert configuration

| Field | Units | Description |
|---|---|---|
| `landingLonDeg` | deg | Landing site longitude |
| `landingLatDeg` | deg | Landing site latitude |
| `rfLanding` | m, 3x1 ENU | Final position relative to landing site (typically [0;0;0]) |
| `vfLanding` | m/s, 3x1 ENU | Desired final velocity (typically [0;0;-1] for near-hover) |
| `afLanding` | m/s^2, 3x1 ENU | Desired final acceleration (typically [0;0;2g] for hover) |
| `divertEnabled` | bool | Enable mid-flight target switch |
| `altDivert` | m | Altitude at which the divert event fires |
| `divertPoints` | m, Nx3 ENU | Candidate divert target locations (East, North, Up per row) |

### `optimizationParams` — Optimizer and simulation configuration

| Field | Default | Description |
|---|---|---|
| `paramsX0` | [0.3, 0.4, 700] | Initial guess [gamma1, gamma2, tgo_seconds] |
| `nodeCount` | 301 | Number of equally-spaced nodes for cost/constraint evaluation (must be odd) |
| `glideSlopeEnabled` | false | Enforce glideslope cone constraint |
| `glideSlopeFinalTheta` | 45 | Glideslope half-angle (deg) below `glideSlopeLow` altitude |
| `glideSlopeHigh` | 500 | Altitude (m) above which glideslope is unconstrained (180 deg) |
| `glideSlopeLow` | 250 | Altitude (m) below which glideslope is at `glideSlopeFinalTheta` |
| `freeGlideNodes` | 1 | Terminal nodes exempted from glideslope (avoids singularity at landing) |
| `pointingEnabled` | false | Enforce thrust-pointing angle constraint |
| `maxTiltAccel` | 2 | Maximum angular acceleration of pointing limit (deg/s^2) |
| `minPointing` | 10 | Minimum pointing angle from vertical at landing (deg) |
| `updateOpt` | false | Enable periodic re-optimization during simulation |
| `updateFreq` | 30 | Re-optimization interval (s) |
| `updateStop` | 120 | Freeze re-optimization with this many seconds of tgo remaining |
| `gamma1eps` | 1e-2 | Minimum allowed gamma1 |
| `gamma2eps` | 1e-2 | Minimum gamma2 - gamma1 gap (enforced as linear inequality in fmincon) |

---

## 5. Entry Points

These are the five script files intended to be run directly. Each sets up the four input structs, calls the appropriate pipeline function(s), and saves or displays results.

---

### `userTest.m` — Interactive single run or beta sweep

**What it does:** Runs one full optimization + simulation cycle per beta value. Calls `getParams` (no divert) or `getParamsDIVERT` (divert). Generates up to 15 figures and prints a summary table.

**Key parameters at the top of the file:**
```matlab
betaVec = [0.7];             % cost trade-off value(s) to evaluate
paramsIC = [0.3, 0.4, 700];  % initial guess [gamma1, gamma2, tgo_sec]
glideSlopeEnabled     = false;
pointingEnabled       = false;
reOptimizationEnabled = false;
divertEnabled         = false;
```

**Divert setup:** When `divertEnabled = true`, the script builds a grid of candidate divert targets at configurable distances and headings around the nominal landing site (8 directions × N distances). The `divertPoints` matrix is Nx3 ENU. `getParamsDIVERT` will run a full trajectory simulation for every row of that matrix.

**Figures produced (by `plotting.m`):**

| Figure | Name | Content |
|---|---|---|
| 1 | Optim Throttle | Optimizer throttle fraction vs time |
| 2 | Optimization Ground Range vs Alt | Optimal trajectory in range-altitude space |
| 3 | Optimization 3D | Optimal trajectory in ENU (final 2 km), with glideslope cone overlay |
| 4 | Opt Vel | Optimal velocity magnitude vs time |
| 5 | Opt Accel | Optimal acceleration profile |
| 6 | Simulation 3D | Simulated trajectory in ENU (final 2 km) |
| 7 | Simulation Ground Range vs Alt | Simulated trajectory in range-altitude space |
| 8 | Sim Accel | Simulated acceleration profile |
| 9 | Sim Velocity | Simulated velocity magnitude vs time |
| 10 | Sim Throttle | Simulated throttle fraction vs time |
| 11 | Sim Mass | Vehicle mass vs time |
| 12 | Sim Alt | Altitude vs time |
| 13 | Pointing Analysis | Thrust pointing angle vs limit (optimizer trajectory) |
| 14 | Pointing Analysis (Sim) | Thrust pointing angle vs limit (simulation trajectory) |
| 15 | ReOpt Parameter History | gamma1 and gamma2 evolution across re-opt events (only when `updateOpt = true`) |

---

### `gammaSweep.m` — 2D cost map over [gamma1, gamma2]

**What it does:** Sweeps gamma1 and gamma2 over a 2D grid at a fixed `tgo`, evaluating the cost and checking thrust constraints at every point. No `fmincon` is called — this is a direct forward evaluation of the guidance law over the parameter space.

**Key parameters:**
```matlab
fixedTgo = 761.65;         % seconds, fixed for all evaluations
gammaRange = [0.01, 1.0];  % range swept for both gamma1 and gamma2
gridResolution = 1000;     % grid is gridResolution x gridResolution points
betaVal = 0;
```

Points where `gamma2 <= gamma1` (with a small tolerance) are set to `NaN` and shown as blank.

**Output:** Two figures — a filled contour plot of `J(gamma1, gamma2)` and a constraint-activity map (colored green/yellow/red for 0/1/2+ saturated nodes). The minimum feasible cost and its location are printed to the console.

**Calls:** `calculateCoeffs`, `simpsonComp13Integral` (locally via `evaluateTraj`). Does **not** call `getParams`.

---

### `tgoSweep.m` — 1D cost curve over tgo

**What it does:** Sweeps `tgo` at fixed `[gamma1, gamma2]`, evaluating cost and counting active thrust constraints. Useful for seeing how sensitive the cost is to flight time for a given guidance shape.

**Key parameters:**
```matlab
fixedGamma1 = 0.01;
fixedGamma2 = 1.0033;
tgoRange = [566.7, 780];   % seconds
betaVal = 0.0;
numPoints = 1000;
```

**Output:** Prints min/max cost and an approximate quadratic fit to the console. Figure-saving code is present but commented out.

**Calls:** `calculateCoeffs`, `simpsonComp13Integral` (locally via `evaluateTraj`). Does **not** call `getParams`.

---

### `dispersionStudy.m` — State dispersion study

**What it does:** Runs N cases (one per seed file row) in parallel via `parfor`, each with a different initial state drawn from pre-generated normally-distributed seeds. Captures landing accuracy and fuel use across the full dispersion.

**Dispersion sources and default 3-sigma ranges:**

| Variable | 3-sigma Range |
|---|---|
| Altitude | ±200 m |
| Longitude | ±0.25 deg |
| Latitude | ±0.25 deg |
| Inertial velocity | ±3 m/s |
| Flight path angle | ±0.1 deg |
| Azimuth | ±0.2 deg |
| Mass | ±0.66% of nominal |

Seeds are loaded from `Seeds/` (inside the `Dispersion and MonteCarlo` folder). Each seed file contains values drawn from N(0,1); they are scaled by `(3-sigma range) / 3` to produce properly-distributed 3-sigma dispersions.

If a case throws any error (optimizer failure, etc.), the `catch` block fills all numeric result fields with `NaN` and sets `exit_ok = false`, so the overall struct array stays consistent for the table-building code.

**Output:** `.mat` and `.csv` result files in a timestamped subdirectory of `Dispersion301/`. Produces histogram, boxchart, and 3D landing scatter plots via `statsPlotting.m`.

**Requires:** MATLAB Parallel Computing Toolbox (for `parfor`).

---

### `accelMonteCarlo.m` — Accelerometer scale factor Monte Carlo

**What it does:** Like `dispersionStudy.m`, but varies only one quantity: a multiplicative scale factor applied to the commanded thrust acceleration inside `closedLoopSim`. This models uncertainty in the onboard accelerometer — the vehicle thinks it is applying the commanded thrust, but the actual thrust is `scale * a_T`. Landing accuracy and fuel consumption are recorded across all cases.

**Key parameter:**
```matlab
accel_monte_carlo = 0.10;  % 10% 3-sigma accelerometer error
```

Seeds are loaded from `Seeds/accel_seeds.dat`. The scale factor applied to case `idx` is:
```matlab
accel_scale = accel_seeds(idx) * (accel_monte_carlo / 3);  % N(0,1) * sigma
accel_scale = 1 + accel_scale;                              % centered at 1.0
```

**Output:** `.mat` and `.csv` results in `AccelMonteCarlo301/`. Four figures: position error vs scale factor, 2D landing dispersion scatter with 3-sigma circle, fuel vs scale factor, and component histograms.

---

## 6. Core Library

---

### `getParams.m` — Main pipeline function

**Signature:**
```matlab
[gammaOpt, gamma2Opt, krOpt, tgoOptSec, aTOptim, exitflag, ...
 optFuelCost, simFuelCost, aTSim, finalPosSim, optHistory, ...
 ICstates, exitFlags, problemParams, nonDimParams, refVals, ...
 optTable, simTable, nActiveConstraints] = ...
    getParams(PDIState, planetaryParams, targetState, vehicleParams, ...
              optimizationParams, betaParam, doPlots, verboseOutput, ...
              dispersion, runSimulation [, monteCarloSeed])
```

> **Note on arg 9 (`dispersion`):** This is a legacy argument retained for call-site compatibility. It has no effect on behavior and can be passed as `false`.

> **Note on arg 11 (`monteCarloSeed`):** Optional. When provided, it is passed through to `closedLoopSim` as a multiplicative scalar on the commanded thrust acceleration, modeling accelerometer scale error.

**What it does, step by step:**

1. **Coordinate transform:** Calls `PDI2MCMF` to convert the PDI state into Cartesian position and velocity in the MCMF frame. Converts landing target (ENU) to MCMF via `ENU2MCMF`.

2. **Non-dimensionalization:** Divides all quantities by reference values (`L_ref`, `T_ref`, `V_ref`, etc.) so the optimizer works in a well-conditioned, unit-free space.

3. **Struct packing:** Bundles all parameters into `problemParams`, `nonDimParams`, and `refVals` for clean passing to downstream functions.

4. **Optimization:** Calls `optimizationLoop` with the initial guess. If exit flag is not 1, re-runs from the result of the first round to improve convergence.

5. **Simulation:** Selects the appropriate simulation mode:
   - `runSimulation = false` — skip simulation entirely
   - `monteCarloSeed` provided — call `closedLoopSim` with the seed
   - `updateOpt = false` — call `closedLoopSim` (static)
   - `updateOpt = true` — call `simReOpt` (re-optimizing)

6. **Plotting:** Calls `plotting.m` if `doPlots = true`.

7. **Output tables:** Returns `optTable = [landingError; optFuelCost; optCost]` and `simTable = [simLandingError; simFuelCost; simCost]`.

---

### `getParamsDIVERT.m` — Divert pipeline function

Same structure as `getParams`, but specifically for divert scenarios. Key differences:

- Forces `glideSlopeEnabled = false`, `pointingEnabled = false`, `updateOpt = true`.
- After the initial optimization, calls `simReOpt` in a **loop** over every row of `targetState.divertPoints`, simulating what happens if a divert command fires mid-flight to each candidate target.
- Stores each trajectory in `divertData{idx}` and produces dedicated divert plots: 3D trajectories, final landing scatter, throttle profiles, and gamma1 vs gamma2 parameter evolution across re-optimizations.
- Does not return `optTable`/`simTable`.

---

### `optimizationLoop.m` — fmincon wrapper

**Signature:**
```matlab
[optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag, nActiveConstraints] = ...
    optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimizationParams, refVals, verboseOutput)
```

**What it does:**

Runs `fmincon` (SQP algorithm) over the decision variables `[gamma1, gamma2, tgo_ND]` subject to:

- **Linear inequalities:** `gamma1 >= gamma1eps`, `gamma2 >= gamma1 + gamma2eps`, `tgo >= 0.01`
- **Box bounds:** `gamma1 in [gamma1eps, 5]`, `gamma2 in [0, 10]`, `tgo in [0.01, 15]`
- **Nonlinear constraints** (from local `nonLinearLimits`):
  - Thrust bounds at every node: `minThrust <= m(t) * |a_T(t)| <= maxThrust` (2 * nodeCount constraints)
  - Glideslope cone at every node if enabled (nodeCount constraints)
  - Pointing angle at every node if enabled (nodeCount constraints)

**Objective function** (from local `objectiveFunction`): evaluates the cost integral using `simpsonComp13Integral` after calling `calculateCoeffs` for the candidate `[gamma1, gamma2, tgo]`. Mass is propagated via the Tsiolkovsky rocket equation for the thrust constraints.

After optimization, the function reconstructs the full trajectory at the optimal parameters — position `rdOptim`, velocity `vdOptim`, thrust acceleration `aTOptim`, and mass `mOptim` — as dense arrays over `nodeCount` equally-spaced tgo nodes.

---

### `closedLoopSim.m` — Static simulation

**Signature:**
```matlab
[tTraj, stateTraj, aTList, flag_thrustGotLimited] = ...
    closedLoopSim(gamma, gamma2, tgo0, nonDimParams, refVals [, monteCarloSeed])
```

**What it does:** Simulates the vehicle flying the guidance law from `t = 0` to `t = tgo0` using `ode45`. The guidance parameters `[gamma1, gamma2, tgo0]` are held constant throughout — there is no re-optimization mid-flight. The state vector is `X = [r; v; m]` (7 elements, all non-dimensional).

At each integration step, the ODE function calls `computeAccel` (a local function) which:
1. Evaluates the closed-loop guidance law to get the commanded `a_T`
2. Applies the optional `monteCarloSeed` scale factor
3. Clamps thrust to `[minThrust, maxThrust]` and sets `flag_thrustGotLimited` if clamped

After integration, `aTList` is reconstructed by re-evaluating `computeAccel` at every saved state point.

**Two-phase integration:** See [Section 13](#13-two-phase-ode-integration).

---

### `simReOpt.m` — Re-optimizing simulation

**Signature:**
```matlab
[tTraj, stateTraj, aTList, flag_thrustGotLimited, optHistory, ICstates, exitFlags] = ...
    simReOpt(gamma0, gamma20, tgo0, problemParams, nonDimParams, refVals, ...
             optimizationParams, betaParam, verboseOutput [, divertPoint])
```

**What it does:** Simulates the vehicle in sequential segments. After each segment of length `updateFreq` seconds, calls `optimizationLoop` from the current state to update `[gamma1, gamma2, tgo]`. Re-optimization continues until `tgo <= updateStop`, after which the remaining flight uses fixed parameters with the same two-phase frozen-thrust integration as `closedLoopSim`.

**Divert logic:** If `problemParams.divertEnabled = true` and a `divertPoint` is provided, the ODE integrator monitors altitude. When the vehicle descends below `altDivert`, the landing target is switched to `divertPoint` and the guidance is re-optimized from the current state. The divert can fire during any segment, including the final no-reopt segment.

**`optHistory`** (Nx5 table): records `[t_elapsed, gamma1, gamma2, kr, tgo_ND]` at each re-optimization event.
**`ICstates`** (7xN table): records `[r; v; m]` at each re-optimization event.
**`exitFlags`** (Nx1): the `fmincon` exit flag at each re-optimization.

---

### `calculateCoeffs.m` — Guidance coefficient solver

**Signature:**
```matlab
[c1, c2, c1_num, c2_num] = calculateCoeffs(r, v, tgo, gamma1, gamma2, afStar, rfStar, vfStar, g)
```

Solves the 2×2 linear system derived from the requirement that the trajectory satisfies both position and velocity boundary conditions at `t_go = 0`. The result is two 3×1 coefficient vectors `c1` and `c2` (one per space dimension, all non-dimensional). The optional third and fourth outputs return the numerator vectors before division by the determinant `delta`, which can be used for diagnostics.

---

### `plotting.m` — Visualization

Called automatically by `getParams` when `doPlots = true`. Produces figures 1–15 as described in the `userTest.m` section. Each call overlays on top of existing open figures (matching by figure name), so running a beta sweep produces a single set of figures with all runs shown in different colors with automatic legend entries.

---

### `simpsonComp13Integral.m` — Numerical integrator

Evaluates the integral of `y` over vector `t` using composite Simpson's 1/3 rule. Requires an **odd** number of points (even number of intervals). This is why `nodeCount = 301` is the default.

---

## 7. Coordinate Functions

All live in `CoordinateFunctions/`. Added to the MATLAB path by every entry-point script.

Two coordinate frames are used throughout the codebase:

- **MCMF** (Moon-Centered Moon-Fixed): a Cartesian frame rotating with the Moon. The Moon's center is the origin; the x-axis points toward 0° lat, 0° lon; z-axis points toward the North Pole.
- **ENU** (East-North-Up): a local topocentric frame anchored at a specific lat/lon (always the landing site). East, North, and Up are tangent or normal to the lunar surface at that point.

All physics and ODE integration happen in MCMF. Results are converted back to ENU for plotting and error reporting.

---

### `CoordinateFunctions/PDI2MCMF.m`

**Signature:**
```matlab
[r_mcmf, v_mcmf] = PDI2MCMF(altitude_km, lonInitDeg, latInitDeg, ...
    landingLonDeg, landingLatDeg, inertialVel_mps, flightPathAngleDeg, azimuth_rad, radMoon)
```

Converts the PDI state (scalar altitude, latitude, longitude, speed magnitude, flight path angle, heading) into a 3D MCMF Cartesian position and velocity vector.

> **Units note:** `azimuth` must be in **radians** here. The `PDIState.azimuth` struct field is in degrees; `getParams.m` converts it with `* pi/180` before calling this function.

---

### `CoordinateFunctions/ENU2MCMF.m`

**Signature:**
```matlab
mcmf = ENU2MCMF(enu, anchorLatDeg, anchorLonDeg, isPosition, radBody)
```

Rotates a vector from the local ENU frame (anchored at the given lat/lon) into MCMF. If `isPosition = true`, also adds the MCMF position of the anchor point (the lunar surface at `radBody` radius). Use `isPosition = false` for velocities and accelerations, which transform by rotation only.

---

### `CoordinateFunctions/MCMF2ENU.m`

**Signature:**
```matlab
[enu, alt] = MCMF2ENU(X, landingLatDeg, landingLonDeg, isPosition, radBody)
```

Inverse of `ENU2MCMF`. Converts MCMF Cartesian to ENU relative to the landing site. If `isPosition = true`, subtracts the landing site MCMF position before rotating. The optional second output `alt` is the altitude above the surface at radius `radBody`.

---

### `enuBasis.m`

**Signature:**
```matlab
[E, N, U] = enuBasis(lat_rad, lon_rad)
```

Returns the three orthonormal unit vectors of the local ENU frame at a given point (all in MCMF coordinates). Used internally by all three coordinate conversion functions.

---

## 8. Parameter Sweeps

Both sweep scripts (`ParameterSweeps/gammaSweep.m` and `ParameterSweeps/tgoSweep.m`) are self-contained. They set up the same four input structs as `userTest.m`, perform their own non-dimensionalization and coordinate transforms inline, and then call a local `evaluateTraj` function that directly invokes `calculateCoeffs` and `simpsonComp13Integral`.

**No `fmincon` or ODE integration is used** — this is a fast, direct evaluation of the cost function and constraint activity across a parameter grid.

### `ParameterSweeps/gammaSweep.m`

**Purpose:** Map out the cost landscape over the guidance shape parameter space. Useful for understanding what the optimizer is searching through, choosing good initial guesses, and identifying the feasible region.

**Output:**
- Figure 1: Filled contour plot of `J(gamma1, gamma2)` with minimum marked
- Figure 2: `pcolor` map of active constraint count (0 = green / 1 = yellow / 2+ = red)
- Table `T_Gamma` in workspace (can be saved to CSV by uncommenting the save block)

### `ParameterSweeps/tgoSweep.m`

**Purpose:** Understand how cost varies with time-of-flight at fixed guidance shape. Useful for seeing whether the cost landscape is shallow (optimizer has flexibility on tgo) or steep (tgo is tightly constrained by the minimum-cost condition).

**Output:** Console summary (min/max cost, approximate quadratic fit coefficients). Figure-saving code is present but commented out.

---

## 9. Dispersion and Monte Carlo Studies

Both studies use `parfor` for parallelism and require the **MATLAB Parallel Computing Toolbox**. Both produce timestamped output directories with `.mat` and `.csv` result files.

### Seed files

The `Seeds/` folder (inside `Dispersion and MonteCarlo/`) contains pre-generated files of standard-normal random samples. Each file has one value per row. The studies scale these by `sigma = (3-sigma range) / 3` to produce the correct distribution.

Using pre-generated seeds instead of `randn` at runtime ensures the exact same set of dispersion cases is reproduced in every run — essential for reproducibility.

### `Dispersion and MonteCarlo/dispersionStudy.m`

For each case `idx`, the script:
1. Loads one row from each seed file
2. Constructs a perturbed `PDI` struct and `vehicle` struct by adding `seed * sigma` to each field
3. Calls `getParams(PDI, ..., runSimulation)` — full optimize + simulate pipeline
4. Stores results in `Results(idx)` struct array

If a case throws any error, the `catch` block fills all numeric fields with `NaN` and sets `exit_ok = false`, keeping the struct array shape consistent.

After the loop, `statsPlotting(Results, betaParam, optimizationParams)` produces histograms, boxcharts, and a 3D landing scatter plot for all converged cases (exitflag 1 or 2 only).

### `Dispersion and MonteCarlo/accelMonteCarlo.m`

Same structure as `dispersionStudy`, but the only perturbation is a multiplicative `accel_scale` factor loaded from `accel_seeds.dat`. This is passed as the 11th argument to `getParams` (`monteCarloSeed`), which forwards it to `closedLoopSim`, which multiplies all commanded thrust accelerations by this factor.

Additional result field: `Results(idx).accel_scale` — the scale factor applied to that case.

### `Dispersion and MonteCarlo/statsPlotting.m`

Called by `dispersionStudy` after the parallel loop completes. Filters to cases with `exitflag == 1` or `exitflag == 2`, then produces:
- Histograms of gamma1, gamma2, tgo, fuel_opt, fuel_sim, and four guidance coefficients
- Boxcharts of the same variables
- 3D scatter of final landing positions (East, North, Up error)

---

## 10. DIDO Optimal Control (Optional)

The `DIDO/` folder contains a separate workflow requiring the **DIDO Toolbox** (a commercial MATLAB pseudospectral optimal control solver, not included). This is used for comparison — DIDO solves the true fuel-optimal (or mixed-cost) problem without the FP2DG structure, providing a benchmark to evaluate how close FP2DG comes to true optimality.

### `DIDO/FP2DG_problem.m` — Entry point

Sets up the DIDO problem struct, runs the solver, extracts and redimensionalizes results, verifies the solution (Hamiltonian constancy check + ODE45 propagation check), and produces comparison figures. Mission parameters match `userTest.m`.

### `DIDO/FP2DG_cost.m`, `FP2DG_dynamics.m`, `FP2DG_events.m`, `FP2DG_path.m`

Callback functions called by the DIDO solver:

| File | Purpose |
|---|---|
| `FP2DG_cost` | Defines the objective functional |
| `FP2DG_dynamics` | Equations of motion (point-mass, constant gravity) |
| `FP2DG_events` | Terminal boundary conditions (initial and final state) |
| `FP2DG_path` | Path (thrust saturation) constraints |

---

## 11. Non-Dimensionalization

All physics inside `getParams`, `optimizationLoop`, `closedLoopSim`, and `simReOpt` is performed in non-dimensional units. The reference values are chosen so that typical trajectory quantities are O(1):

| Quantity | Symbol | Value | Notes |
|---|---|---|---|
| Length | `L_ref` | 10 000 m | Fixed |
| Time | `T_ref` | `sqrt(L_ref / g_planet)` | ~78.6 s (lunar) |
| Velocity | `V_ref` | `L_ref / T_ref` | ~127 m/s |
| Acceleration | `A_ref` | `g_planet` | 1.622 m/s^2 |
| Mass | `M_ref` | `massInit` | 15 103 kg |

The non-dimensional gravity vector `gConst` is approximated as constant, evaluated at the landing site:
```
gConst = -(rPlanet_ND)^2 * rfStar_ND / |rfStar_ND|^3
```
This is valid for the short-duration, relatively-flat powered descent arc.

The reference values are packed into the `refVals` struct and returned by `getParams`. To convert outputs back to physical units:

```matlab
traj_pos_m   = stateTraj(:, 1:3) * refVals.L_ref;
traj_vel_mps = stateTraj(:, 4:6) * refVals.V_ref;
traj_mass_kg = stateTraj(:, 7)   * refVals.M_ref;
accel_mps2   = aTSim              * refVals.A_ref;
time_sec     = tTraj              * refVals.T_ref;
```

---

## 12. Constraints

### Thrust saturation (always active)

At every `nodeCount` node, the thrust magnitude must satisfy:
```
minThrust <= m(t) * |a_T(t)| <= maxThrust
```
These are `2 * nodeCount` nonlinear inequality constraints in `fmincon`. The mass `m(t)` at each node is computed from the Tsiolkovsky rocket equation using the accumulated thrust integral from landing backwards.

### Glideslope cone (optional — `glideSlopeEnabled`)

The position trajectory must lie within a cone centered on the vertical axis through the landing site. The cone half-angle varies with altitude:

| Altitude | Half-angle |
|---|---|
| Above `glideSlopeHigh` (500 m) | 180 deg (no constraint) |
| Between `glideSlopeLow` and `glideSlopeHigh` | Linear interpolation |
| Below `glideSlopeLow` (250 m) | `glideSlopeFinalTheta` (45 deg) |

The bottom `freeGlideNodes` nodes are exempted to avoid numerical issues as the trajectory converges exactly to the landing site.

### Thrust pointing angle (optional — `pointingEnabled`)

The angle between the thrust vector and the local vertical must not exceed a time-varying limit:
```
phi(t) <= Theta(t_go) = phi_0 + 0.5 * maxTiltAccel * t_go^2
```
This means: far from landing (large `t_go`), tilting is freely allowed; near landing, the vehicle must be nearly vertical. `phi_0 = minPointing` sets the minimum angle at touchdown.

---

## 13. Two-Phase ODE Integration

The closed-loop guidance law has a singularity as `t_go -> 0` (the `1/t_go` and `1/t_go^2` terms diverge). To prevent numerical blow-up in the ODE solver, both `closedLoopSim` and the final segment of `simReOpt` use a two-phase integration strategy:

**Phase 1: Guidance active** (`t = 0` to `t = tgo0 - 0.2 s`)

The ODE function evaluates the closed-loop guidance law at every step, using the current state `[r, v]` and the current `t_go = tgo0 - t`. This is the true closed-loop behavior: the guidance reacts to the actual trajectory.

**Phase 2: Frozen thrust** (`t = tgo0 - 0.2 s` to `t = tgo0`)

The thrust acceleration computed at `t_go = 0.2 s` is held fixed. The ODE becomes a simple constant-thrust propagation problem, with no singularity. The 0.2 s threshold is hard-coded and converted to non-dimensional units as `minTime = 0.2 / T_ref`.

---

## 14. File Call Graph

```
userTest.m
|-- PDI2MCMF         <- uses enuBasis
|-- ENU2MCMF         <- uses enuBasis
|-- getParams
|   |-- PDI2MCMF
|   |-- ENU2MCMF
|   |-- optimizationLoop
|   |   |-- calculateCoeffs
|   |   `-- simpsonComp13Integral
|   |-- closedLoopSim              (if updateOpt = false)
|   |   `-- calculateCoeffs        (via computeAccel)
|   |-- simReOpt                   (if updateOpt = true)
|   |   |-- optimizationLoop
|   |   `-- calculateCoeffs        (via computeAccel)
|   |-- MCMF2ENU     <- uses enuBasis
|   `-- plotting
|       |-- MCMF2ENU
|       |-- enuBasis
|       `-- calculateCoeffs
`-- getParamsDIVERT                (if divertEnabled)
    |-- optimizationLoop
    |-- simReOpt (once per divert target)
    `-- MCMF2ENU

accelMonteCarlo.m
`-- getParams  (with monteCarloSeed -> closedLoopSim)

dispersionStudy.m
|-- getParams
`-- statsPlotting

gammaSweep.m / tgoSweep.m
|-- calculateCoeffs         (via local evaluateTraj)
`-- simpsonComp13Integral

DIDO/FP2DG_problem.m        (requires DIDO Toolbox)
|-- FP2DG_cost
|-- FP2DG_dynamics
|-- FP2DG_events
`-- FP2DG_path
```

---

## 15. Output Variables Guide

Key variables returned by `getParams`:

| Variable | Size | Units | Description |
|---|---|---|---|
| `gammaOpt` | scalar | — | Optimal gamma1 |
| `gamma2Opt` | scalar | — | Optimal gamma2 |
| `krOpt` | scalar | — | (gamma2+2)(gamma1+2) — composite guidance parameter |
| `tgoOptSec` | scalar | s | Optimal time-of-flight (dimensional) |
| `aTOptim` | 3 x nodeCount | ND | Thrust acceleration along optimizer trajectory |
| `exitflag` | scalar | — | fmincon exit flag (1 = fully converged, 2 = near-converged, <=0 = failed) |
| `optFuelCost` | scalar | kg | Propellant consumed along optimizer trajectory |
| `simFuelCost` | scalar | kg | Propellant consumed along simulation trajectory |
| `aTSim` | 3 x N | ND | Thrust acceleration along simulation trajectory |
| `finalPosSim` | 3x1 | m, ENU | Final position relative to landing site (simulation) |
| `optHistory` | table (K x 5) | ND | [t_elapsed, gamma1, gamma2, kr, tgo] at each re-opt event |
| `ICstates` | table (K x 7) | ND | [r_x, r_y, r_z, v_x, v_y, v_z, m] at each re-opt event |
| `exitFlags` | K x 1 | — | fmincon exit flags at each re-opt event |
| `problemParams` | struct | dimensional | All dimensional problem parameters |
| `nonDimParams` | struct | ND | All non-dimensional parameters (r0, v0, isp, thrust bounds, etc.) |
| `refVals` | struct | — | L_ref, T_ref, V_ref, A_ref, M_ref |
| `optTable` | 3 x 1 | m, kg, — | [landing_error; opt_fuel_kg; opt_cost] |
| `simTable` | 3 x 1 | m, kg, — | [sim_landing_error; sim_fuel_kg; sim_cost] |
| `nActiveConstraints` | scalar | — | Number of nearly-active constraints at the optimum |

**Redimensionalizing the simulation trajectory:**
```matlab
traj_pos_m   = stateTraj(:, 1:3) * refVals.L_ref;   % meters, MCMF
traj_vel_mps = stateTraj(:, 4:6) * refVals.V_ref;   % m/s, MCMF
traj_mass_kg = stateTraj(:, 7)   * refVals.M_ref;   % kg
accel_mps2   = aTSim              * refVals.A_ref;   % m/s^2, MCMF
time_sec     = tTraj              * refVals.T_ref;   % seconds

% Convert position to ENU relative to landing site
pos_ENU = MCMF2ENU(traj_pos_m', landingLatDeg, landingLonDeg, true, rPlanet);
```
