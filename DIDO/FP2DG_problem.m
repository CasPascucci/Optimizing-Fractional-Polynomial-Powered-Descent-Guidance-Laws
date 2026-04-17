clear all; clc;
addpath([pwd, '/../CoordinateFunctions']);
addpath([pwd, '/../']);

%% 1. Mission Parameters (identical to userTest.m)
beta = 0.96; % 1.0 = Fuel Optimal, 0.0 = Smoothest Throttle

PDIState = struct('altitude', 15240, ...       % m
                  'lonInitDeg', 41.85, ...     % deg
                  'latInitDeg', -71.59, ...    % deg
                  'inertialVelocity', 1693.8, ... % m/s
                  'flightPathAngleDeg', 0, ... % deg
                  'azimuth', 180);             % deg

planetaryParams = struct('rPlanet', 1736.01428 * 1000, ... % m
                         'gPlanet', 1.622, ...             % m/s^2
                         'gEarth', 9.81);                  % m/s^2

vehicleParams = struct('massInit', 15103.0, ... % kg
                       'dryMass', 6855, ...     % kg
                       'isp', 311, ...          % s
                       'maxThrust', 45000, ...  % N
                       'minThrust', 4500);      % N

targetState = struct('landingLonDeg', 41.85, ...
                     'landingLatDeg', -90.00, ...
                     'rfLanding', [0;0;0], ...
                     'vfLanding', [0;0;-1]);    % m/s

%% 2. Reference Values & Non-Dimensionalization
rPlanet = planetaryParams.rPlanet;
gPlanet = planetaryParams.gPlanet;
gEarth  = planetaryParams.gEarth;

L_ref = 10000;
T_ref = sqrt(L_ref / gPlanet);
A_ref = gPlanet;
V_ref = L_ref / T_ref;
M_ref = vehicleParams.massInit;

fprintf('Reference values: L=%.1f m, T=%.3f s, V=%.3f m/s, A=%.3f m/s^2, M=%.1f kg\n', ...
    L_ref, T_ref, V_ref, A_ref, M_ref);

%% 3. Coordinate Transforms to MCMF
altitude_km        = PDIState.altitude / 1000;
lonInitDeg         = PDIState.lonInitDeg;
latInitDeg         = PDIState.latInitDeg;
inertialVelocity   = PDIState.inertialVelocity;
flightPathAngleDeg = PDIState.flightPathAngleDeg;
azimuth            = PDIState.azimuth * pi / 180;

landingLonDeg = targetState.landingLonDeg;
landingLatDeg = targetState.landingLatDeg;

[r0Dim, v0Dim] = PDI2MCMF(altitude_km, lonInitDeg, latInitDeg, ...
                           landingLonDeg, landingLatDeg, ...
                           inertialVelocity, flightPathAngleDeg, azimuth, rPlanet);

rfDim = 10000 * ENU2MCMF(targetState.rfLanding/10000, landingLatDeg, landingLonDeg, true);
vfDim = ENU2MCMF(targetState.vfLanding, landingLatDeg, landingLonDeg, false);

%% 4. Non-Dimensional State Values
rPlanetND   = rPlanet / L_ref;
r0ND        = r0Dim / L_ref;
v0ND        = v0Dim / V_ref;
rfStarND    = rfDim / L_ref;
vfStarND    = vfDim / V_ref;
m0ND        = 1.0;
dryMassND   = vehicleParams.dryMass / M_ref;
ispND       = vehicleParams.isp * gEarth / V_ref;
maxThrustND = vehicleParams.maxThrust / (M_ref * A_ref);
minThrustND = vehicleParams.minThrust / (M_ref * A_ref);

fprintf('\nNon-dimensional initial state:\n');
fprintf('  r0 = [%.6f; %.6f; %.6f]\n', r0ND);
fprintf('  v0 = [%.6f; %.6f; %.6f]\n', v0ND);
fprintf('  m0 = %.6f\n', m0ND);
fprintf('Non-dimensional target state:\n');
fprintf('  rf = [%.6f; %.6f; %.6f]\n', rfStarND);
fprintf('  vf = [%.6f; %.6f; %.6f]\n', vfStarND);
fprintf('Thrust bounds (ND): [%.4f, %.4f]\n', minThrustND, maxThrustND);
fprintf('Isp (ND): %.4f\n', ispND);

%% 5. DIDO Problem Definition
constants.nonDim.rMoonND = rPlanetND;
constants.nonDim.ispND   = ispND;
constants.beta           = beta;

FP2DG.cost     = 'FP2DG_cost';
FP2DG.dynamics = 'FP2DG_dynamics';
FP2DG.events   = 'FP2DG_events';
FP2DG.path     = 'FP2DG_path';
FP2DG.constants = constants;

%% 6. Search Space
search.states = [-200, 200;   % r_x
                 -200, 200;   % r_y
                 -200, 200;   % r_z
                  -15,  15;   % v_x
                  -15,  15;   % v_y
                  -15,  15;   % v_z
                  0.3,  1.1]; % m

search.controls = [-5, 5;  % aT_x
                   -5, 5;  % aT_y
                   -5, 5]; % aT_z

FP2DG.search = search;

%% 7. Bounds
bounds.events = [r0ND(1),     r0ND(1);      % r_x(0)
                 r0ND(2),     r0ND(2);      % r_y(0)
                 r0ND(3),     r0ND(3);      % r_z(0)
                 v0ND(1),     v0ND(1);      % v_x(0)
                 v0ND(2),     v0ND(2);      % v_y(0)
                 v0ND(3),     v0ND(3);      % v_z(0)
                 m0ND,        m0ND;         % m(0)
                 rfStarND(1), rfStarND(1);  % r_x(tf)
                 rfStarND(2), rfStarND(2);  % r_y(tf)
                 rfStarND(3), rfStarND(3);  % r_z(tf)
                 vfStarND(1), vfStarND(1);  % v_x(tf)
                 vfStarND(2), vfStarND(2);  % v_y(tf)
                 vfStarND(3), vfStarND(3)]; % v_z(tf)

bounds.path = [minThrustND, maxThrustND];

bounds.initial.time = [0, 0];    % t0 fixed
bounds.final.time   = [1, 15];   % tf free, ~7 ND expected

FP2DG.bounds = bounds;

%% 8. Run DIDO
check(FP2DG);

algorithm.nodes = 60; % Lower doesn't improve cost, but landing accuracy for ode check starts to suffer
fprintf('\n=== Running DIDO (%d nodes) ===\n', algorithm.nodes);
tic;
[cost, primal, dual] = dido(FP2DG, algorithm);
runTime = toc;
fprintf("*****************\n");
fprintf('DIDO run time: %.2f s\n', runTime);
fprintf('Cost: %.6f\n', cost);

%% 9. Extract & Redimensionalize Results
t_DIDO  = primal.time;
r_DIDO  = primal.states(1:3, :);
v_DIDO  = primal.states(4:6, :);
m_DIDO  = primal.states(7, :);
aT_DIDO = primal.controls(1:3, :);

t_dim  = t_DIDO * T_ref;      % s
r_dim  = r_DIDO * L_ref;      % m
v_dim  = v_DIDO * V_ref;      % m/s
m_dim  = m_DIDO * M_ref;      % kg
aT_dim = aT_DIDO * A_ref;     % m/s^2

aT_mag_dim = vecnorm(aT_dim);
thrust_dim = m_dim .* aT_mag_dim;  % N
throttle   = thrust_dim / vehicleParams.maxThrust;

r_ENU = zeros(3, size(r_dim, 2));
for i = 1:size(r_dim, 2)
    r_ENU(:,i) = MCMF2ENU(r_dim(:,i), landingLatDeg, landingLonDeg, true, true);
end
alt_dim = vecnorm(r_dim, 2, 1) - rPlanet; % m

fuelCost = (m_dim(1) - m_dim(end)); % kg
finalPosENU = MCMF2ENU(r_dim(:,end), landingLatDeg, landingLonDeg, true, true);
landingError = norm(finalPosENU);

%% 10. Results Summary
fprintf('\n--- DIDO Results ---\n');
fprintf('Cost:            %.4f\n', cost);
fprintf('Time of Flight:  %.2f sec (%.4f ND)\n', t_dim(end) - t_dim(1), t_DIDO(end) - t_DIDO(1));
fprintf('Fuel Used:       %.2f kg\n', fuelCost);
fprintf('Final Mass:      %.2f kg\n', m_dim(end));
fprintf('Landing Error:   %.4f m\n', landingError);
fprintf('Final Pos (ENU): [%.4f, %.4f, %.4f] m\n', finalPosENU);
fprintf('Max Throttle:    %.4f\n', max(throttle));
fprintf('Min Throttle:    %.4f\n', min(throttle));

%% 11. Verify and Validate
H = dual.Hamiltonian;
H_mean = mean(H);
H_std  = std(H);
fprintf('\nV&V - Hamiltonian:\n');
fprintf('  Mean: %.6f, Std: %.6e, Max deviation: %.6e\n', H_mean, H_std, max(abs(H - H_mean)));

% ODE Check
fprintf('\nV&V - ODE45 Trajectory Check...\n');
aT_interp = griddedInterpolant(t_DIDO, aT_DIDO', 'spline');
X0_prop = [r_DIDO(:,1); v_DIDO(:,1); m_DIDO(1)];

odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[~, X_prop] = ode45(@(t, X) propagateDynamics(t, X, aT_interp, rPlanetND, ispND), ...
    [t_DIDO(1), t_DIDO(end)], X0_prop, odeopts);

r_prop_final = X_prop(end, 1:3)' * L_ref;
v_prop_final = X_prop(end, 4:6)' * V_ref;
m_prop_final = X_prop(end, 7) * M_ref;

fprintf('  Position error (DIDO vs ODE45): %.6f m\n', norm(r_prop_final - r_dim(:,end)));
fprintf('  Velocity error (DIDO vs ODE45): %.6f m/s\n', norm(v_prop_final - v_dim(:,end)));
fprintf('  Mass error: %.6f kg\n', abs(m_prop_final - m_dim(end)));
fprintf('  ODE45 fuel used: %.2f kg\n', M_ref - m_prop_final);

%% 13. Plots
figure('Name', 'DIDO Throttle Profile'); hold on; grid on;
plot(t_dim, throttle, 'b-', 'LineWidth', 2, 'DisplayName', 'DIDO Optimal Throttle');
yline(1.0, 'k--', 'LineWidth', 1, 'DisplayName', 'Max Throttle');
yline(vehicleParams.minThrust/vehicleParams.maxThrust, 'k--', 'LineWidth', 1, 'DisplayName', 'Min Throttle');
xlabel('Time (s)'); ylabel('Throttle Fraction');
title('DIDO Optimal Throttle Profile');
ylim([0, 1.2]);
legend('Location', 'best');
set(gca, 'FontSize', 20);

groundRange = vecnorm(r_ENU(1:2,:), 2, 1);
figure('Name', 'DIDO Altitude vs Ground Range');
plot(groundRange/1000, alt_dim/1000, 'b-', 'LineWidth', 2, 'DisplayName', 'DIDO Trajectory');
xlabel('Ground Range (km)'); ylabel('Altitude (km)');
title('DIDO Optimal Trajectory');
legend('Location', 'best'); grid on;
set(gca, 'FontSize', 20);

figure('Name', 'DIDO 3D Trajectory');
plot3(r_ENU(1,:)/1000, r_ENU(2,:)/1000, r_ENU(3,:)/1000, 'b-', 'LineWidth', 2, 'DisplayName', 'DIDO Trajectory');
xlabel('East (km)'); ylabel('North (km)'); zlabel('Up (km)');
title('DIDO 3D Trajectory (ENU)');
legend('Location', 'best'); grid on; axis equal;
set(gca, 'FontSize', 20);

v_ENU = zeros(3, size(v_dim, 2));
for i = 1:size(v_dim, 2)
    v_ENU(:,i) = MCMF2ENU(v_dim(:,i), landingLatDeg, landingLonDeg, false, true);
end
figure('Name', 'DIDO Velocity');
plot(t_dim, v_ENU(1,:), 'r-', t_dim, v_ENU(2,:), 'g-', t_dim, v_ENU(3,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('East', 'North', 'Up'); title('DIDO Velocity (ENU)');
grid on;
set(gca, 'FontSize', 20);

aT_ENU = zeros(3, size(aT_dim, 2));
for i = 1:size(aT_dim, 2)
    aT_ENU(:,i) = MCMF2ENU(aT_dim(:,i), landingLatDeg, landingLonDeg, false, true);
end
figure('Name', 'DIDO Acceleration');
plot(t_dim, aT_ENU(1,:), 'r-', t_dim, aT_ENU(2,:), 'g-', t_dim, aT_ENU(3,:), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
legend('East', 'North', 'Up'); title('DIDO Commanded Acceleration (ENU)');
grid on;
set(gca, 'FontSize', 20);

figure('Name', 'DIDO Mass'); hold on; grid on;
plot(t_dim, m_dim, 'b-', 'LineWidth', 2, 'DisplayName', 'DIDO Optimal Fuel Consumption');
yline(vehicleParams.dryMass, 'k--', 'LineWidth', 1, 'DisplayName', 'Dry Mass');
xlabel('Time (s)'); ylabel('Mass (kg)');
title('DIDO Mass History');
legend('Location', 'best');
set(gca, 'FontSize', 20);
propUsedDIDO = max(0, m_dim(1) - m_dim(end));
subtitle(sprintf('Propellant used = %.1f kg', propUsedDIDO));

figure('Name', 'DIDO Hamiltonian');
plot(t_DIDO, H, 'b-', 'LineWidth', 2, 'DisplayName', 'DIDO Solved Hamiltonian');
xlabel('Time (ND)'); ylabel('Hamiltonian');
title('DIDO Hamiltonian');
legend('Location', 'best'); grid on;
set(gca, 'FontSize', 20);

fprintf('\n=== DIDO Complete ===\n');
fprintf('Compare fuel used %.2f kg against FP2DG result from userTest.m\n', fuelCost);

%% Local Functions
function dXdt = propagateDynamics(t, X, aT_interp, rMoonND, ispND)
    r    = X(1:3);
    v    = X(4:6);
    mass = X(7);

    aT = aT_interp(t)'; % interpolate DIDO control
    rMag = norm(r);
    g = -(rMoonND^2) * r / (rMag^3);

    aTmag = norm(aT);
    dmdt = -(aTmag * mass) / ispND;

    dXdt = [v; aT + g; dmdt];
end
