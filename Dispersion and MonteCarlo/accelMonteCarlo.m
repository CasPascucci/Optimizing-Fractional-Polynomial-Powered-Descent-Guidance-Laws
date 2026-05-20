clear all; clc; %close all;
addpath(fileparts(mfilename('fullpath')) + "\..");
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'CoordinateFunctions'));
PDINom = struct;
PDINom.altitude           = 15240;
PDINom.lonInitDeg         = 41.85;
PDINom.latInitDeg         = -71.59;
PDINom.inertialVelocity   = 1693.8;
PDINom.flightPathAngleDeg = 0;
PDINom.azimuth            = 180;

planetaryParams = struct;
planetaryParams.rPlanet = 1736.01428 * 1000; % m
planetaryParams.gPlanet = 1.622; % m/s^2
planetaryParams.gEarth = 9.81; % m/s^2

vehicleNom = struct;
vehicleNom.massInit = 15103.0; % kg
vehicleNom.dryMass = vehicleNom.massInit - 8248; % kg
vehicleNom.isp = 311; % s
vehicleNom.maxThrust = 45000; % N
vehicleNom.minThrust = 4500; % N

targetState = struct;
targetState.landingLonDeg      = 41.85;
targetState.landingLatDeg      = -90.0;
targetState.rfLanding = [0;0;0];
targetState.vfLanding = [0;0;-1];
targetState.afLanding = [0;0;2*planetaryParams.gPlanet];
targetState.divertEnabled = false; % No divert in dispersion test
targetState.altDivert = 0;
targetState.divertPoints = 0;

targetState.delta_t   = 5; % seconds dim, for btt, not implemented

optimizationParams = struct;
optimizationParams.paramsX0 = [0.3, 0.4, 700];
optimizationParams.nodeCount = 301; %Count must be odd for Simpson
optimizationParams.gamma1eps = 1e-2;
optimizationParams.gamma2eps = 1e-2;

optimizationParams.glideSlopeFinalTheta = 45; %deg
optimizationParams.glideSlopeHigh = 500; %m
optimizationParams.glideSlopeLow = 250; %m
optimizationParams.freeGlideNodes = 1; %m
optimizationParams.glideSlopeEnabled = true;

optimizationParams.pointingEnabled = true;
optimizationParams.maxTiltAccel = 2; % deg/s^2
optimizationParams.minPointing = 10; %deg, floor for pointing constraint

optimizationParams.updateFreq = 10;
optimizationParams.updateStop = 120;
optimizationParams.updateOpt = false;

betaParam = 0.6;
runSimulation = true; % needs to be true
doPlotting = false; % disable this to not plot results
verboseOutput = false;

%% Setup Stats 3 Sigma Ranges
seedDir = fileparts(mfilename("fullpath"));
seedFilenames = fullfile(seedDir,'Seeds/accel_seeds.dat');
accel_seeds =load(seedFilenames);

% 3-sigma ranges
accel_monte_carlo = 0.10; % percent accel


caseCount = numel(accel_seeds);
[queue, done] = waitBarQueue(caseCount, 'Acceleration Monte Carlo');
Results = struct('k',{},'gamma',{},'gamma2',{},'kr',{},'tgo',{},'fuel_opt',{},'fuel_sim',{},...
    'final_error',{},'coeff1',{},'coeff2',{},'coeff3',{},'coeff4',{},'exit_ok',{},'msg',{}, 'exitflag',{}, 'accel_scale',{});

L_ref = 10000; A_ref = planetaryParams.gPlanet; T_ref = sqrt(L_ref/A_ref); V_ref = L_ref/T_ref; M_ref = 15103.0;
dispTime = tic;
%% Stats Loop - Requires Parallel Computing Toolbox
parfor (idx = 1:caseCount)
    try
        accel_scale = accel_seeds(idx);

        PDI = PDINom;
        vehicle = vehicleNom;


        [gammaOpt, gamma2Opt, krOpt, tgoOpt, ~, exitflag, optFuel, simFuel, ~, finalPosSim, ~, ~, ~] = getParams(PDI, planetaryParams, targetState, vehicle, optimizationParams, betaParam, doPlotting, verboseOutput, true, runSimulation, accel_scale);
        
        Results(idx).k = idx;
        Results(idx).gamma = gammaOpt;
        Results(idx).gamma2 = gamma2Opt;
        Results(idx).kr = krOpt;
        Results(idx).tgo = tgoOpt;
        Results(idx).fuel_opt = optFuel;
        Results(idx).fuel_sim = simFuel;
        Results(idx).final_error = finalPosSim;
        Results(idx).coeff1 = gammaOpt*(krOpt/(2*gammaOpt +4) -1);
        Results(idx).coeff2 = (gammaOpt*krOpt/(2*gammaOpt+4)-gammaOpt-1);
        Results(idx).coeff3 = ((gammaOpt+1)/tgoOpt)*(1-krOpt/(gammaOpt+2));
        Results(idx).coeff4 = (krOpt/tgoOpt^2);
        Results(idx).exit_ok = true;
        Results(idx).msg = '';
        Results(idx).exitflag = exitflag;
        Results(idx).accel_scale = accel_scale;
    catch ME
        Results(idx).k = idx;
        Results(idx).exit_ok = false;
        Results(idx).msg = ME.message;
    end
    send(queue,1);
end
elapsed = toc(dispTime)
done();

%% Save Results
dispDir = 'AccelMonteCarlo301';
if ~exist(dispDir,'dir')
    mkdir(dispDir);
end

% Timestamp for the run
timeRun = datestr(now,'yyyymmdd_HHMMSS');

% Build constraint tags based on which conditions are active
condTags = {};
if isfield(optimizationParams,'glideSlopeEnabled') && optimizationParams.glideSlopeEnabled
    condTags{end+1} = 'GS';
end
if isfield(optimizationParams,'pointingEnabled') && optimizationParams.pointingEnabled
    condTags{end+1} = 'PT';
end
if isempty(condTags)
    condTags = {'NONE'};
end
condTag = strjoin(condTags,'_');
betaVal = betaParam * 100;
suffix  = sprintf('%dbeta_%s_%s', round(betaVal), condTag, timeRun);

runTitle = sprintf('Monte Carlo Study (Beta: %.2f | Constraints: %s)', betaParam, strjoin(condTags, ', '));
runName = suffix;
runDir  = fullfile(dispDir, runName);
mkdir(runDir);

% Filenames
matfile = fullfile(runDir, ['results_dispersion_' runName '.mat']);
csvfile = fullfile(runDir, ['results_dispersion_' runName '.csv']);

% Save results
save(matfile, 'Results');

% Export summary CSV - might need non matlab app to examine results in
% future
T = table( (1:numel(Results)).', [Results.exit_ok].',[Results.gamma].', [Results.gamma2].', [Results.kr].', ...
           [Results.tgo].', [Results.fuel_opt].', [Results.fuel_sim].', vecnorm([Results.final_error].',2,2), [Results.coeff1].', [Results.coeff2].', ...
           [Results.coeff3].', [Results.coeff4].', ...
           'VariableNames', {'k','exit_ok','gamma','gamma2','kr','tgo','fuel_opt','fuel_sim','final_error_norm','coeff1','coeff2','coeff3','coeff4'} );
writetable(T, csvfile);

%% Generate Stats Plots and Tables
ok_mask = [Results.exit_ok];
R_ok = Results(ok_mask);
accel_scales = [R_ok.accel_scale];
final_errors = [R_ok.final_error]; % 3 x N
error_norms = vecnorm(final_errors, 2, 1);
fuel_sim_vals = [R_ok.fuel_sim];

% Figure 1 — Position Error Norm vs Acceleration Scale Factor
figure('Name', 'Position Error vs Accel Scale');
scatter(accel_scales, error_norms, 20, 'filled');
xlabel('Acceleration Scale Factor');
ylabel('Position Error Norm (m)');
title(sprintf('Position Error vs Acceleration Scale Factor (\\beta = %.2f)', betaParam));
grid on;
set(gca, 'FontSize', 20);

% Figure 2 — 2D Landing Dispersion (East vs North)
figure('Name', '2D Landing Dispersion');
scatter(final_errors(1,:), final_errors(2,:), 20, 'filled');
hold on;
sigma_E = std(final_errors(1,:));
sigma_N = std(final_errors(2,:));
r3sigma = 3 * max(sigma_E, sigma_N);
theta = linspace(0, 2*pi, 200);
plot(r3sigma*cos(theta), r3sigma*sin(theta), 'r--', 'LineWidth', 1.5);
hold off;
xlabel('East Error (m)');
ylabel('North Error (m)');
title(sprintf('Landing Dispersion (\\beta = %.2f) — 3\\sigma circle r=%.1f m', betaParam, r3sigma));
axis equal; grid on;
set(gca, 'FontSize', 20);

% Figure 3 — Fuel Consumed vs Acceleration Scale Factor
figure('Name', 'Fuel vs Accel Scale');
scatter(accel_scales, fuel_sim_vals, 20, 'filled');
xlabel('Acceleration Scale Factor');
ylabel('Fuel Consumed (kg)');
title(sprintf('Fuel Consumed vs Acceleration Scale Factor (\\beta = %.2f)', betaParam));
grid on;
set(gca, 'FontSize', 20);

% Figure 4 — Position Error Component Histograms
figure('Name', 'Position Error Histograms');
comp_labels = {'East', 'North', 'Up'};
for ci = 1:3
    subplot(3,1,ci);
    histogram(final_errors(ci,:), 30);
    xlabel([comp_labels{ci} ' Error (m)']);
    ylabel('Count');
    title(sprintf('%s Position Error (\\mu=%.2f m, \\sigma=%.2f m)', comp_labels{ci}, mean(final_errors(ci,:)), std(final_errors(ci,:))));
    grid on;
    set(gca, 'FontSize', 20);
end


% Stats Summary Table
vals = struct();
fields = {'gamma','gamma2','kr','tgo','fuel_opt','fuel_sim','coeff1','coeff2','coeff3','coeff4'};

for f = 1:numel(fields)
    x = [Results.(fields{f})];
    ok = [Results.exit_ok];
    x_ok = x(ok);

    vals.(fields{f}).mean = mean(x_ok);
    vals.(fields{f}).std  = std(x_ok,0);
    vals.(fields{f}).min  = min(x_ok);
    vals.(fields{f}).max  = max(x_ok);
    vals.(fields{f}).range3sigma = [vals.(fields{f}).mean - 3*vals.(fields{f}).std, ...
                                    vals.(fields{f}).mean + 3*vals.(fields{f}).std];
end

% Create summary table
names = fieldnames(vals);
meanVals = []; stdVals = []; minVals = []; maxVals = []; lowVals = []; highVals = [];

for i = 1:numel(names)
    meanVals(i,1) = vals.(names{i}).mean;
    stdVals(i,1)  = vals.(names{i}).std;
    minVals(i,1)  = vals.(names{i}).min;
    maxVals(i,1)  = vals.(names{i}).max;
    lowVals(i,1)  = vals.(names{i}).range3sigma(1);
    highVals(i,1) = vals.(names{i}).range3sigma(2);

end

SummaryTable = table(names, meanVals, stdVals, minVals, maxVals, lowVals, highVals, ...
    'VariableNames', {'Parameter','Mean','StdDev','Min','Max','-3σ','+3σ'});

disp(SummaryTable);


% Exit Flag/Conditions Summary
allExit = [Results.exitflag]';
ok      = [Results.exit_ok]';
N       = numel(allExit);

% Count per unique exit flag
[uniqFlags, ~, idxu] = unique(allExit);
counts   = accumarray(idxu, 1);
percents = 100 * counts / N;

% Exit Flag Meanings for table
meaning = strings(numel(uniqFlags),1);
for i = 1:numel(uniqFlags)
    switch uniqFlags(i)
        case 1
            meaning(i) = "Optimality and Constraint tolerances met";
        case 0
            meaning(i) = "Stopped by limit";
        case -2
            meaning(i) = "Infeasible / no feasible point";
        case -3
            meaning(i) = "Constraint tolerance met, went under ObjectiveLimit";
        case 2
            meaning(i) = "Under StepTolerance, Constraint tolerance met";
        otherwise
            meaning(i) = "Shouldn't be here with interior point";
    end
end

ExitFlagTable = table(uniqFlags, counts, percents, meaning, ...
    'VariableNames', {'ExitFlag','Count','Percent','Meaning'});

% Success/failure roll-up (based on exit_ok)
succCount = sum(ok);
failCount = N - succCount;
succPct   = 100 * succCount / N;
failPct   = 100 - succPct;

SuccessTable = table( ...
    categorical({'success';'failure'}), ...
    [succCount; failCount], ...
    [succPct;   failPct], ...
    'VariableNames', {'Outcome','Count','Percent'});



% Table figures
disp('Exit flag summary:'); disp(ExitFlagTable);
disp('Success summary:');   disp(SuccessTable);

%% Rerun section
%rerunDispersionCase(1614, PDINom, planetaryParams, targetState, vehicleNom, optimizationParams, beta, seeds)
%% Save both tables
exitCsv   = fullfile(runDir, 'exitflag_summary.csv');
succCsv   = fullfile(runDir, 'success_summary.csv');
writetable(ExitFlagTable, exitCsv);
writetable(SuccessTable,  succCsv);
summaryFile = fullfile(runDir, 'results_dispersion_summary.csv');
writetable(SummaryTable, summaryFile);

%% Progress Bar Function
function [queue, done] = waitBarQueue(caseCount, titleStr)
    queue = parallel.pool.DataQueue;
    p = 0;
    h = waitbar(0, sprintf('0 / %d', caseCount), 'Name', titleStr);

    afterEach(queue, @updateWait);

    function updateWait(~)
        p = p + 1;
        if isvalid(h)
            waitbar(p/caseCount, h, sprintf('%d / %d', p, caseCount));
        end
    end

    done = @() (delete(h));
end