clear all;  clc; format short
%close all;
addpath([pwd, '/CoordinateFunctions']);

%% Key Parameters
beta = 0.9;  % Weighting: 1.0 = Fuel Optimal, 0.0 = Smoothest Throttle

% First IC set is for everything but GS
paramsIC = [0.3, 0.4, 700]; % Initial guess in optimization for gamma1, gamma2, tgo (dimensional seconds)
%paramsIC = [0.5, 0.8, 700]; % Use for GS, recommended to stay to higher values

glideSlopeEnabled = false;
pointingEnabled = false;
reOptimizationEnabled = false;
divertEnabled = false; % Will internally force reOpt On, glideSlope and pointing Off

%% 1. Initial State Definitions
% PDI (Powered Descent Initiation) State, Dimensional
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
                     'divertEnabled', divertEnabled); % Requires reopt on, will disable pointing and glideslope if not already done, and will enable Re-Opt

%% 2. Optimization Configuration
optimizationParams = struct;
optimizationParams.paramsX0 = paramsIC;

optimizationParams.nodeCount = 301; % Node Count for Optimization, must be odd for Simpson's rule

% Glideslope Constraints
optimizationParams.glideSlopeEnabled    = glideSlopeEnabled;
optimizationParams.glideSlopeFinalTheta = 45;  % deg
optimizationParams.glideSlopeHigh       = 500; % m
optimizationParams.glideSlopeLow        = 250; % m
optimizationParams.freeGlideNodes     = 1;  % Number of optimization nodes starting at landing site, counting backwards, to leave unconstrained for glideslope. Default is 1, as landing node requires infinite precision

% Pointing Constraints
optimizationParams.pointingEnabled = pointingEnabled;
optimizationParams.maxTiltAccel    = 2;  % deg/s^2
optimizationParams.minPointing     = 10; % deg

% Re-Optimization Settings
optimizationParams.updateOpt  = reOptimizationEnabled; 
optimizationParams.updateFreq = 30;   % s
optimizationParams.updateStop = 120;   % s (Time before landing to stop updates)

% Tolerances
optimizationParams.gamma1eps = 1e-2;
optimizationParams.gamma2eps = 1e-2;

% Divert
divertDistances = [1000, 2000, 3000]; % Distances in meters away from original site, the 8-point rings of divert will occur at each of these distances

divert1E = zeros(length(divertDistances),1); % each line here sets up the North and East Coordinates of the divert points
divert1N = divertDistances';
divert2E = divertDistances'.*cosd(45);
divert2N = divert2E;
divert3E = divertDistances';
divert3N = zeros(length(divertDistances),1);
divert4E = divertDistances'.*cosd(45);
divert4N = -divert4E;
divert5E = zeros(length(divertDistances),1);
divert5N = -divertDistances';
divert6E = -divertDistances'.*cosd(45);
divert6N = divert6E;
divert7E = -divertDistances';
divert7N = zeros(length(divertDistances),1);
divert8E = -divertDistances'.*cosd(45);
divert8N = -divert8E;


% Group the divert coordinates into a matrix. Will be 8n+1 x 3, where n is
% number of divert distances chosen above
targetState.divertPoints = [0, 0, 0;
                            divert1E, divert1N, zeros(length(divertDistances),1);
                            divert2E, divert2N, zeros(length(divertDistances),1);
                            divert3E, divert3N, zeros(length(divertDistances),1);
                            divert4E, divert4N, zeros(length(divertDistances),1);
                            divert5E, divert5N, zeros(length(divertDistances),1);
                            divert6E, divert6N, zeros(length(divertDistances),1);
                            divert7E, divert7N, zeros(length(divertDistances),1);
                            divert8E, divert8N, zeros(length(divertDistances),1)];
%targetState.divertPoints = [-150, -300, 0]; % meters % Single case, above block is multiple divert plots
targetState.altDivert = 1000; % m, altitude above ground to trigger divert scenario

%% 3. Execution Flags & Run
runSimulation = true;
doPlotting    = true; 
verboseOutput = true;
timePreParam = toc;

if ~ targetState.divertEnabled
[gammaOpt, gamma2Opt, krOpt, tgoOptSec, ~, ~, optFuelCost, simFuelCost, ...
 aTSim, finalPosSim, optHistory, ICstates, exitFlags, problemParams, ...
 nonDimParams, refVals, optTable, simTable] = getParams(PDIState, planetaryParams, targetState, ...
    vehicleParams, optimizationParams, beta, doPlotting, verboseOutput, false, runSimulation);
else
 [gammaOpt, gamma2Opt, krOpt, tgoOptSec, ~, ~, optFuelCost, simFuelCost, ...
 aTSim, finalPosSim, optHistory, ICstates, exitFlags, problemParams, ...
 nonDimParams, refVals] = getParamsDIVERT(PDIState, planetaryParams, targetState, ...
    vehicleParams, optimizationParams, beta, doPlotting, verboseOutput, false, runSimulation);
end

%% 4. Results Display
if optimizationParams.updateOpt
    fprintf("\n--- Initial Solution ---\n")
else
    fprintf('\n--- Final Results ---\n');
end
fprintf('Gamma1:      %.4f\n', gammaOpt);
fprintf('Gamma2:      %.4f\n', gamma2Opt);
fprintf('Kr:          %.4f\n', krOpt);
fprintf('Tgo (sec):   %.2f\n', tgoOptSec);
fprintf('Cost Value:  %.4f\n', optTable(3));
fprintf('Opt Fuel:    %.2f kg\n', optFuelCost);
fprintf('Sim Fuel:    %.2f kg\n', simFuelCost);
if exist("optTable","var")
    table(optTable,simTable, 'VariableNames',["Optimization", "Simulation"],'RowNames',["Landing Error", "Fuel Cost", "Cost Function"]);
end

%% 5. Single Segment Re-Run
% Use this block to isolate and troubleshoot specific re-optimization segments
if optimizationParams.updateOpt
    % segmentIdx = 63;
    % outputSingle = reOptReRun(segmentIdx, ICstates, optHistory, beta, ...
    %     problemParams, nonDimParams, optimizationParams, refVals, ...
    %     targetState.delta_t / refVals.T_ref, verboseOutput);
end