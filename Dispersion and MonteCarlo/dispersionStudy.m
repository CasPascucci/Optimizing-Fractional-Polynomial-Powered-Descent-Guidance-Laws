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
targetState.divertPoints = zeros(1,3);

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
seedFilenames = struct( ...
    'alt', fullfile(seedDir,'Seeds/alt_seeds.dat'), ...
    'lon', fullfile(seedDir,'Seeds/lon_seeds.dat'), ...
    'lat', fullfile(seedDir,'Seeds/lat_seeds.dat'), ...
    'v', fullfile(seedDir,'Seeds/V_seeds.dat'), ...
    'fpa', fullfile(seedDir,'Seeds/fpa_seeds.dat'), ...
    'azmth', fullfile(seedDir,'Seeds/azmth_seeds.dat'), ...
    'mass', fullfile(seedDir,'Seeds/mass_seeds.dat'));

seeds = struct(...
    'alt_seeds', load(seedFilenames.alt),...
    'lon_seeds', load(seedFilenames.lon),...
    'lat_seeds', load(seedFilenames.lat),...
    'v_seeds', load(seedFilenames.v),...
    'fpa_seeds', load(seedFilenames.fpa),...
    'azmth_seeds', load(seedFilenames.azmth),...
    'mass_seeds', load(seedFilenames.mass));

% 3-sigma ranges
r_disp = 200; % m
lon_disp = 0.25; % deg
lat_disp = 0.25; % deg
v_disp = 3; % m/s
fpa_disp = 0.1; % deg
azmth_disp = 0.2; % deg
mass_disp = 0.0066; % fraction of nominal mass

caseCount = numel(seeds.alt_seeds);
[queue, done] = waitBarQueue(caseCount, 'Dispersion');
Results = struct('k',{},'gamma',{},'gamma2',{},'kr',{},'tgo',{},'fuel_opt',{},'fuel_sim',{},...
    'final_error',{},'coeff1',{},'coeff2',{},'coeff3',{},'coeff4',{},'exit_ok',{},'msg',{}, 'exitflag',{});

L_ref = 10000; A_ref = planetaryParams.gPlanet; T_ref = sqrt(L_ref/A_ref); V_ref = L_ref/T_ref; M_ref = 15103.0;
dispTime = tic;
%% Stats Loop - Requires Parallel Computing Toolbox
parfor (idx = 1:caseCount)
    try
        dalt = seeds.alt_seeds(idx) * (r_disp / 3);
        dlon = seeds.lon_seeds(idx) * (lon_disp / 3);
        dlat = seeds.lat_seeds(idx) * (lat_disp / 3);
        dv = seeds.v_seeds(idx) * (v_disp / 3);
        dfpa = seeds.fpa_seeds(idx) * (fpa_disp / 3);
        dazmth = seeds.azmth_seeds(idx) * (azmth_disp / 3);
        mass_mult = seeds.mass_seeds(idx) * (mass_disp / 3);

        PDI = PDINom;
        PDI.altitude_km = PDINom.altitude / 1000 + dalt/1000;
        PDI.lonInitDeg = PDINom.lonInitDeg + dlon;
        PDI.latInitDeg = PDINom.latInitDeg + dlat;
        PDI.inertialVelocity = PDINom.inertialVelocity + dv;
        PDI.flightPathAngleDeg = PDINom.flightPathAngleDeg + dfpa;
        PDI.azimuth = PDINom.azimuth + dazmth;

        vehicle = vehicleNom;
        vehicle.massInit = vehicleNom.massInit * (1+ mass_mult);
        vehicle.dryMass = vehicle.massInit - 8248; % this way makes the dry mass variable but fuel amount constant


        [gammaOpt, gamma2Opt, krOpt, tgoOpt, ~, exitflag, optFuel, simFuel, ~, finalPosSim, ~, ~, ~] = getParams(PDI, planetaryParams, targetState, vehicle, optimizationParams, betaParam, doPlotting, verboseOutput, true, runSimulation);
        
        Results(idx).k = idx;
        Results(idx).gamma = gammaOpt;
        Results(idx).gamma2 = gamma2Opt;
        Results(idx).kr = krOpt;
        Results(idx).tgo = tgoOpt;
        Results(idx).fuel_opt = optFuel;
        Results(idx).fuel_sim = simFuel; % New
        Results(idx).final_error = finalPosSim; % New
        Results(idx).coeff1 = gammaOpt*(krOpt/(2*gammaOpt +4) -1);
        Results(idx).coeff2 = (gammaOpt*krOpt/(2*gammaOpt+4)-gammaOpt-1);
        Results(idx).coeff3 = ((gammaOpt+1)/tgoOpt)*(1-krOpt/(gammaOpt+2));
        Results(idx).coeff4 = (krOpt/tgoOpt^2);
        Results(idx).exit_ok = true;
        Results(idx).msg = '';
        Results(idx).exitflag = exitflag;
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
dispDir = 'Dispersion301';
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

runTitle = sprintf('Dispersion Study (Beta: %.2f | Constraints: %s)', betaParam, strjoin(condTags, ', '));
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
statsPlotting(Results, betaParam, optimizationParams); % Only plots results with flag 1 or 2


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