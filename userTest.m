clear all; clc; format short
addpath([pwd, '/CoordinateFunctions']);

%% Key Parameters
betaVec = [0.7]; % Called betaVec because a vector of beta values can be given, with each one being ran one after another on the same figures for comparison
%betaVec = 1.0:-0.05:0.0 % Example of alternate way to give a beta vector

paramsIC = [0.3, 0.4, 700]; % Good for everything but GS
% paramsIC = [0.5, 0.8, 700]; % Use for GS

glideSlopeEnabled    = false;
pointingEnabled      = false;
reOptimizationEnabled = false;
divertEnabled        = false; % will internally disable GS and PT, and will enable reOpt

%% 1. Initial State Definitions
PDIState = struct('altitude', 15240, ...
                  'lonInitDeg', 41.85, ...
                  'latInitDeg', -71.59, ...
                  'inertialVelocity', 1693.8, ...
                  'flightPathAngleDeg', 0, ...
                  'azimuth', 180);

planetaryParams = struct('rPlanet', 1736.01428 * 1000, ...
                         'gPlanet', 1.622, ...
                         'gEarth', 9.81);

vehicleParams = struct('massInit', 15103.0, ...
                       'dryMass', 6855, ...
                       'isp', 311, ...
                       'maxThrust', 45000, ...
                       'minThrust', 4500);

targetState = struct('landingLonDeg', 41.85, ...
                     'landingLatDeg', -90.00, ...
                     'rfLanding', [0;0;0], ...
                     'vfLanding', [0;0;-1], ...
                     'afLanding', [0;0;2*planetaryParams.gPlanet], ...
                     'divertEnabled', divertEnabled);

%% 2. Optimization Configuration
optimizationParams = struct;
optimizationParams.paramsX0       = paramsIC;
optimizationParams.nodeCount       = 301;
optimizationParams.glideSlopeEnabled    = glideSlopeEnabled;
optimizationParams.glideSlopeFinalTheta = 45;
optimizationParams.glideSlopeHigh       = 500;
optimizationParams.glideSlopeLow        = 250;
optimizationParams.freeGlideNodes       = 1;
optimizationParams.pointingEnabled = pointingEnabled;
optimizationParams.maxTiltAccel    = 2;
optimizationParams.minPointing     = 10;
optimizationParams.updateOpt  = reOptimizationEnabled;
optimizationParams.updateFreq = 30;
optimizationParams.updateStop = 120;
optimizationParams.gamma1eps  = 1e-2;
optimizationParams.gamma2eps  = 1e-2;

% Divert setup
divertDistances = [1000, 2000, 3000];
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
targetState.altDivert = 1000;

%% 3. Beta Sweep Loop
runSimulation = true;
doPlotting    = true;
verboseOutput = false;  % Recommend false during sweeps

nBeta   = length(betaVec);
results = struct(); % Preallocate results struct

for i = 1:nBeta
    beta = betaVec(i);
    if ~targetState.divertEnabled
        [gammaOpt, gamma2Opt, krOpt, tgoOptSec, ~, optExitFlag, optFuelCost, simFuelCost, ...
         aTSim, finalPosSim, optHistory, ICstates, exitFlags, problemParams, ...
         nonDimParams, refVals, optTable, simTable, nActive] = getParams(PDIState, planetaryParams, targetState, ...
            vehicleParams, optimizationParams, beta, doPlotting, verboseOutput, false, runSimulation);
    else
        [gammaOpt, gamma2Opt, krOpt, tgoOptSec, ~, optExitFlag, optFuelCost, simFuelCost, ...
         aTSim, finalPosSim, optHistory, ICstates, exitFlags, problemParams, ...
         nonDimParams, refVals] = getParamsDIVERT(PDIState, planetaryParams, targetState, ...
            vehicleParams, optimizationParams, beta, doPlotting, verboseOutput, false, runSimulation);
        optTable  = [NaN; NaN; NaN];
        simTable  = [NaN; NaN; NaN];
    end

    if ~targetState.divertEnabled

        % Store results for this beta
        results(i).beta        = beta;
        results(i).gammaOpt    = gammaOpt;
        results(i).gamma2Opt   = gamma2Opt;
        results(i).krOpt       = krOpt;
        results(i).tgoOptSec   = tgoOptSec;
        results(i).optFuelCost = optFuelCost;
        results(i).simFuelCost = simFuelCost;
        results(i).costValue   = optTable(3);
        results(i).simCost     = simTable(3);
        results(i).landingError = optTable(1);
        results(i).exitFlag       = optExitFlag;
        results(i).nActive        = nActive;
        results(i).exitFlags      = exitFlags;
        if ~isempty(finalPosSim)
            results(i).simLandingError = norm(finalPosSim);  % [m]
        else
            results(i).simLandingError = NaN;
        end
    end
end

%% 4. Summary Table
if ~targetState.divertEnabled
    fprintf('\n\n========== BETA SWEEP SUMMARY ==========\n');
    fprintf('%-8s %-10s %-10s %-8s %-12s %-8s %-12s %-8s %-12s %-12s %-8s %-8s\n', ...
            'Beta', 'Gamma1', 'Gamma2', 'Kr', 'Tgo (s)', 'Cost', 'Opt Fuel', 'Sim Cost', 'Sim Fuel', 'SimErr (m)', 'ExitFlag', 'NActive');
    fprintf('%s\n', repmat('-', 1, 106));
    
    for i = 1:nBeta
        r = results(i);
        fprintf('%-8.4f %-10.4f %-10.4f %-8.4f %-12.2f %-8.4f %-12.2f %-8.4f %-12.2f %-12.3f %-8d %-8d\n', ...
                r.beta, r.gammaOpt, r.gamma2Opt, r.krOpt, ...
                r.tgoOptSec, r.costValue, r.optFuelCost, r.simCost, r.simFuelCost, r.simLandingError, r.exitFlag, r.nActive);
    
    end
end
%% 5. Optional: Save results
% save('betaSweepResults.mat', 'results', 'betaVec');