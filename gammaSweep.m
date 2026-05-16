clear all; clc; close all;
addpath([pwd, '/CoordinateFunctions']);

% Key Parameters
fixedTgo = 761.65;
gammaRange = [0.01, 1.0];
betaVal = 0;
gridResolution = 1000;

glideSlopeEnabled = false; % keep all of these flags false here, constraints aren't used in these tests
pointingEnabled = false;
divertEnabled = false;

nodeCount = 301;
% Setup
PDIState = struct('altitude', 15240, ... % meters
                  'lonInitDeg', 41.85, ... % deg
                  'latInitDeg', -71.59, ... % deg
                  'inertialVelocity', 1693.8, ... % m/s
                  'flightPathAngleDeg', 0, ... % deg, angle below horizontal for initial trajectory
                  'azimuth', 180); % deg, heading angle, clockwise from North

planetaryParams = struct('rPlanet', 1736.01428 * 1000, ... % m
                         'gPlanet', 1.622, ... % m/s^2
                         'gEarth', 9.81); % m/s^2

vehicleParams = struct('massInit', 15103.0, ... % kg
                       'dryMass', 6855, ... % kg
                       'isp', 311, ... % seconds
                       'maxThrust', 45000, ... % Newtons
                       'minThrust', 4500); % N

targetState = struct('landingLonDeg', 41.85, ... % deg
                     'landingLatDeg', -90.00, ... % deg
                     'rfLanding', [0;0;0], ...  % meters, centered at lunar radius at above landing coordinates
                     'vfLanding', [0;0;-1], ...  % m/s
                     'afLanding', [0;0;2*planetaryParams.gPlanet], ... % m/s^2
                     'delta_t', 5, ... % seconds, not currently implemented
                     'divertEnabled', divertEnabled, ... % keep false here
                     'altDivert', 1000, ... % not used here
                     'divertPoints', [0,0,0]); % not used here

L_ref = 10000;
T_ref = sqrt(L_ref/planetaryParams.gPlanet);
refVals = struct('L_ref', L_ref, ...
                 'T_ref', T_ref, ...
                 'A_ref', planetaryParams.gPlanet, ...
                 'V_ref', L_ref/T_ref, ...
                 'M_ref', vehicleParams.massInit);

[r0Dim, v0Dim] = PDI2MCMF(PDIState.altitude/1000, PDIState.lonInitDeg, PDIState.latInitDeg, ...
    targetState.landingLonDeg, targetState.landingLatDeg, PDIState.inertialVelocity, ...
    PDIState.flightPathAngleDeg, PDIState.azimuth*pi/180, planetaryParams.rPlanet);

rfDim = refVals.L_ref * ENU2MCMF(targetState.rfLanding/10000, targetState.landingLatDeg, targetState.landingLonDeg, true);
vfDim = ENU2MCMF(targetState.vfLanding, targetState.landingLatDeg, targetState.landingLonDeg, false);
afDim = ENU2MCMF(targetState.afLanding, targetState.landingLatDeg, targetState.landingLonDeg, false);

rfStarND = rfDim/L_ref;
vfStarND = vfDim/refVals.V_ref;
afStarND = afDim/refVals.A_ref;
nonDim = struct('r0ND',r0Dim/L_ref, ...
                'v0ND',v0Dim/refVals.V_ref, ...
                'rfStarND',rfStarND, ...
                'vfStarND',vfStarND, ...
                'afStarND',afStarND, ...
                'gConst',-(planetaryParams.rPlanet/L_ref)^2 * rfStarND / (norm(rfStarND)^3), ...
                'm0ND', 1.0,...
                'ispND', vehicleParams.isp * 9.81 / refVals.V_ref,...
                'maxThrustND',vehicleParams.maxThrust / (vehicleParams.massInit * refVals.A_ref), ...
                'minThrustND',vehicleParams.minThrust / (vehicleParams.massInit * refVals.A_ref));

fprintf('Gamma Sweep, Fixed Tgo=%.0fs\n', fixedTgo);


g1Vec = linspace(gammaRange(1), gammaRange(2), gridResolution);
g2Vec = linspace(gammaRange(1), gammaRange(2), gridResolution);
[G1, G2] = meshgrid(g1Vec, g2Vec);

costGammaSweep = zeros(gridResolution,gridResolution);
activeGammaSweep = zeros(gridResolution,gridResolution);

for idx = 1:numel(G1)
    if abs(G1(idx) - G2(idx)) < 1e-4
        costGammaSweep(idx) = NaN;
        activeGammaSweep(idx) = NaN;
        continue;
    end
    if G2(idx) < G1(idx)
        costGammaSweep(idx) = NaN;
        activeGammaSweep(idx) = NaN;
        continue;
    end
    [cost, nActive] = evaluateTraj(G1(idx),G2(idx), fixedTgo, betaVal, nonDim, refVals, nodeCount);
    costGammaSweep(idx) = cost;
    activeGammaSweep(idx) = nActive;
end

maxConstraints = max(activeGammaSweep, [], 'all', 'omitnan');
if maxConstraints >= 2
    fprintf("Thrust Saturated at more than one point on at least one trajectory\n")
end

maxCost = max(costGammaSweep, [], 'all', 'omitnan');
feasibleMask = ~isnan(costGammaSweep) & (activeGammaSweep <= 1);
feasibleCost = costGammaSweep;
feasibleCost(~feasibleMask) = NaN;
minCost = min(feasibleCost, [], 'all', 'omitnan');
if any(feasibleMask, 'all')
    [minRow, minCol] = find(feasibleCost == minCost, 1);
    [maxRow, maxCol] = find(costGammaSweep == maxCost, 1);
    bestG1min = G1(minRow, minCol);
    bestG2min = G2(minRow, minCol);
    bestG1max = G1(maxRow, maxCol);
    bestG2max = G2(maxRow, maxCol);
    fprintf("Min feasible cost of %.1f @ G1=%.2f, G2=%.2f\n", minCost, bestG1min, bestG2min);
    fprintf("Max cost of %.1f @ G1=%.2f, G2=%.2f\n", maxCost, bestG1max, bestG2max);
else
    fprintf("No valid solutions found.\n");
end

folderName = sprintf('Gamma Sweep/Tgo_%d', round(fixedTgo));
if ~exist(folderName,"dir")
    mkdir(folderName);
end

% Figure 1: Cost Contour (replace surf with contourf)
f1 = figure('Name','Cost Sensitivity to Gamma'); hold on;
contourf(G1, G2, costGammaSweep, 30, 'LineColor', 'none','DisplayName','');
cb = colorbar;
cb.Label.String = 'Cost $J$';
cb.Label.Interpreter = 'latex';
% Draw the gamma2 = gamma1 constraint boundary
diagLine = linspace(gammaRange(1), gammaRange(2), 200);
plot(diagLine, diagLine, 'w--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
% Mark the minimum
hMin = plot(bestG1min, bestG2min, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
legend(hMin, 'Minimum', 'Interpreter', 'latex', 'Location', 'northwest');
xlabel('$\gamma_1$', 'Interpreter', 'latex');
ylabel('$\gamma_2$', 'Interpreter', 'latex');
title('Cost Sensitivity to $[\gamma_1, \gamma_2]$', 'Interpreter', 'latex');
subtitle(sprintf('$\\beta$ = %.2f, $t_{go}$ = %.2f s', betaVal, fixedTgo), 'Interpreter', 'latex');
set(gca, 'FontSize', 20);
saveas(f1, fullfile(folderName, 'Cost_Contour.png'));
saveas(f1, fullfile(folderName, 'Cost_Contour.fig'));

% Figure 2: Active Constraints (discrete colormap)
f2 = figure('Name','Active Constraints Map'); hold on;
pcolor(G1, G2, activeGammaSweep); shading flat;
colormap(f2, [0.2 0.6 0.2; 1 0.8 0; 0.8 0.1 0.1]);  % green/yellow/red for 0/1/2+
caxis([0, 2]);  % pin limits so 0 always maps to green regardless of data range
cb2 = colorbar('Ticks', [0, 1, 2], 'TickLabels', {'0', '1', '2+'});
cb2.Label.String = 'Active Thrust Constraints';
plot(diagLine, diagLine, 'w--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hMin2 = plot(bestG1min, bestG2min, 'b*', 'MarkerSize', 14, 'LineWidth', 2);
legend(hMin2, 'Minimum', 'Interpreter', 'latex', 'Location', 'northwest');
xlabel('$\gamma_1$', 'Interpreter', 'latex');
ylabel('$\gamma_2$', 'Interpreter', 'latex');
title('Active Thrust Constraints', 'Interpreter', 'latex');
subtitle(sprintf('$t_{go}$ = %d s', round(fixedTgo)), 'Interpreter', 'latex');
set(gca, 'FontSize', 20);
saveas(f2, fullfile(folderName, 'Constraints_Map.png'));
saveas(f2, fullfile(folderName, 'Constraints_Map.fig'));



T_Gamma = table(G1(:), G2(:), costGammaSweep(:), activeGammaSweep(:), ...
    'VariableNames', {'Gamma1', 'Gamma2', 'FuelCost_kg', 'ActiveConstraints'});
writetable(T_Gamma, fullfile(folderName, 'Gamma_Sweep_Data.csv'));

fprintf('Results saved to: %s\n', folderName);
%% Functions
function [cost, activeCount] = evaluateTraj(gamma1, gamma2, tgoSec, betaVal, nonDim, refVals, nodeCount)
    tgoND = tgoSec / refVals.T_ref;
    try
        [c1, c2] = calculateCoeffs(nonDim.r0ND, nonDim.v0ND, tgoND, gamma1, gamma2, nonDim.afStarND, nonDim.rfStarND, nonDim.vfStarND, nonDim.gConst);
    catch
        cost = NaN; activeCount = NaN; return;
    end
    tgospan = linspace(0, tgoND, nodeCount);
    aT = nonDim.afStarND + c1.*(tgospan.^gamma1) + c2.*(tgospan.^gamma2);
    aTmag = vecnorm(aT, 2, 1);

    Q = cumtrapz(tgospan, aTmag ./ nonDim.ispND);
    Q = Q(end) - Q;
    simpson1 = simpsonComp13Integral(tgospan,aTmag);
    simpson2 = simpsonComp13Integral(tgospan,dot(aT,aT));
    cost = betaVal*simpson1 + (1-betaVal)*simpson2;

    mCurrent = nonDim.m0ND * exp(-Q);
    thrustMag = mCurrent .* aTmag;

    tol = 1e-4;
    activeNodes = sum((thrustMag - nonDim.maxThrustND) > tol) + sum((nonDim.minThrustND - thrustMag) > tol);
    activeCount = activeNodes;
end    