clear all; clc; close all;
addpath([pwd, '/CoordinateFunctions']);


% Key Parameters 
fixedGamma1 = 0.01;
fixedGamma2 = 1.0033;
tgoRange = [566.7, 780];
betaVal = 0.0;
% IF console prints message stating that too many constraints are active,
% adjust tgoRange to a safer region, constraints at beginning of
% activeTgoSweep mean the range should be moved to higher values, and vice
% versa for constraints at end of activeTgoSweep

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
refVals = struct('L_ref', L_ref, 'T_ref', T_ref, 'A_ref', planetaryParams.gPlanet, ...
    'V_ref', L_ref/T_ref, 'M_ref', vehicleParams.massInit);

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

fprintf('Tgo Sweep, Fixed Params: G1=%.2f, G2=%.2f\n', fixedGamma1, fixedGamma2);

numPoints = 1000;

tgoVec = linspace(tgoRange(1),tgoRange(2),numPoints);
costTgoSweep = zeros(numPoints,1);
activeTgoSweep = zeros(numPoints,1);

for idx = 1: numPoints
    [cost, nActive] = evaluateTraj(fixedGamma1,fixedGamma2, tgoVec(idx), betaVal, nonDim, refVals, nodeCount);
    costTgoSweep(idx) = cost;
    activeTgoSweep(idx) = nActive;
end

% Check if all tests are within constraints
maxConstraints = max(activeTgoSweep, [], 'all', 'omitnan');
if maxConstraints >= 2
    fprintf("Thrust Saturated at more than one point on at least one trajectory\n")
end

% Find max and minimum fuel, and slope approximation between them
[maxCost, maxCostIdx] = max(costTgoSweep);
[minCost, minCostIdx] = min(costTgoSweep);
maxCostTime = tgoVec(maxCostIdx);
minCostTime = tgoVec(minCostIdx);
fprintf("Max cost of %.1f @ tgo = %.1f\n", maxCost, maxCostTime);
fprintf("Min cost of %.1f @ tgo = %.1f\n", minCost, minCostTime);

% Approximate as a quadratic
coeffs = polyfit(tgoVec, costTgoSweep, 2);
fprintf("Approx quadratic: %.4fx^2 + %.4fx + %.2f\n", coeffs);
fprintf("Approx slope: %.4f\n", (maxCost-minCost)/(maxCostTime-minCostTime));

%% Saving
folderName = sprintf('Tgo Sweep/G1_%.2f_G2_%.2f', fixedGamma1, fixedGamma2);
if ~exist(folderName,"dir")
    mkdir(folderName);
end

f1 = figure('Name','Cost Function Sensitivty to Tgo Sweep'); hold on;
plot(tgoVec, costTgoSweep, 'b-', 'LineWidth',2);
hMin = plot(minCostTime, minCost, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
legend(hMin, sprintf('Minimum ($t_{go}$ = %.1f s)', minCostTime), 'Interpreter', 'latex', 'Location', 'northeast');
xlabel('Time of Flight of Trajectory'); ylabel('Cost Function Evaluation of Trajectory');
title('Cost Function Sensitivty to $t_{go}$ Sweep','Interpreter','latex');
subtitle(sprintf("$\\beta$ = %.2f, $\\gamma_1$ = %.4f, $\\gamma_2$ = %.4f", betaVal, fixedGamma1, fixedGamma2), 'Interpreter','latex');
set(gca, 'FontSize', 20);



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

end