function [tTraj, stateTraj, aTList, flag_thrustGotLimited, optHistory, ICstates, exitFlags] = simReOpt(gamma0,gamma20,tgo0, problemParams, nonDimParams, refVals, delta_t, optimizationParams, betaParam, verboseOutput, divertPoint)

if nargin < 11
    divertPoint = [];
end


% Original IC's
r0 = nonDimParams.r0ND;
v0 = nonDimParams.v0ND;
m0 = nonDimParams.m0ND;
rfStar = nonDimParams.rfStarND;
vfStar = nonDimParams.vfStarND;
afStar = nonDimParams.afStarND;
isp = nonDimParams.ispND;
rMoonND = nonDimParams.rMoonND;
gGuidance = nonDimParams.gConst;

X0 = [r0; v0; m0];

updateFreqND = optimizationParams.updateFreq / refVals.T_ref;
updateStopND = optimizationParams.updateStop / refVals.T_ref;

% Check if doing Divert
divertEnabled = problemParams.divertEnabled;
divertOccured = false;
if divertEnabled
    altDivertND = problemParams.altDivert /refVals.L_ref;
    eventDivert = @(t,y) divertTrigger(t, y, nonDimParams.rMoonND, altDivertND);
    
else
    odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
end
% Setup Results
kr0 = (gamma20+2)*(gamma0+2);

t_elapsed = 0;
gamma = gamma0;
gamma2 = gamma20;
kr = kr0;
tgo = tgo0;

tTraj = [];
stateTraj = [];
ICstates = [];
ICstates(:,1) = X0;
ICcounter = 2;
aTList = [];
flag_thrustGotLimited = false;
optHistory = [t_elapsed, gamma, gamma2, kr, tgo];
exitFlags = [];
exitFlags(1) = 99;



if verboseOutput
    fprintf('\n=== Starting Re-Optimization Simulation ===\n');
    fprintf('Initial tgo: %.2f s\n', tgo0 * refVals.T_ref);
    fprintf('Update frequency: %.2f s\n', optimizationParams.updateFreq);
    fprintf('Update stop time: %.2f s (%.2f s remaining)\n', ...
        optimizationParams.updateStop, tgo0 * refVals.T_ref);
    fprintf('==========================================\n\n');
end
segmentCount = 0;
minTime = 0.2/refVals.T_ref; % Time to stop sim at the end

% Main loop
    while tgo > updateStopND % keep looping until update stop time
        segmentCount = segmentCount + 1; % Start a new segment
        tgo_remaining = tgo - updateStopND;
        t_segment = min(tgo,min(updateFreqND, tgo_remaining));
    
        if t_segment <= 0
            break;
        end
    
        t_start = t_elapsed; % Start segment at total elapsed time
        t_end = t_elapsed + t_segment;
    
        if verboseOutput
            fprintf('--- Segment %d ---\n', segmentCount);
            fprintf('Time: %.2f to %.2f s (elapsed)\n', t_start * refVals.T_ref, t_end * refVals.T_ref);
            fprintf('Tgo at start: %.2f s\n', tgo * refVals.T_ref);
            fprintf('Current params: gamma=%.4f, gamma2=%.4f, kr=%.4f\n', gamma, gamma2, kr);
        end
        if divertEnabled && ~divertOccured
            odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', eventDivert);
            [tSeg, stateSeg, tEvent, ~, ~] = ode45(@(t, X) trajectorySegmentODE(t, X, t_elapsed, gamma, kr, tgo, ...
                isp, rMoonND, rfStar, vfStar, afStar, gGuidance, nonDimParams), ...
                [t_start, t_end], X0, odeoptions);
        else
            odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
            [tSeg, stateSeg] = ode45(@(t, X) trajectorySegmentODE(t, X, t_elapsed, gamma, kr, tgo, ...
                isp, rMoonND, rfStar, vfStar, afStar, gGuidance, nonDimParams), ...
                [t_start, t_end], X0, odeoptions);
            tEvent = [];
        end

        if ~isempty(tEvent) % If divert during this segment
            if verboseOutput
                fprintf('= Divert Triggered at t=%.2f s =\n', tEvent(end)*refVals.T_ref);
            end
            rfStar = divertPoint;
            nonDimParams.rfStarND = divertPoint;
            
            divertOccured = true;
        end
    
        tTraj = [tTraj;tSeg];
        stateTraj = [stateTraj;stateSeg];
    
        numPoints = length(tSeg);
        segAT = zeros(3,numPoints);
        for i = 1:numPoints
            r = stateSeg(i,1:3)';
            v = stateSeg(i,4:6)';
            m = stateSeg(i,7);
            t_since_reopt = tSeg(i) - t_elapsed;
            tgo_now = tgo - t_since_reopt;
            [aTi, limited] = computeAccel(r, v, m, gamma, kr, tgo_now, rfStar, vfStar, afStar, gGuidance, nonDimParams);
            segAT(:, i) = aTi;
            flag_thrustGotLimited = flag_thrustGotLimited || limited;
        end
        aTList = [aTList, segAT];
    
        X0 = stateSeg(end,:)'; % Next segment IC is the end of this segment
        ICstates(:,ICcounter) = X0;
        ICcounter = ICcounter + 1;

        t_elapsed = tSeg(end);
        tgo = tgo - t_segment;
    
        if verboseOutput
            fprintf('Segment complete. New tgo: %.2f s\n', tgo * refVals.T_ref);
        end
        if tgo <= updateStopND
            if verboseOutput
                fprintf('Reached update stop time. Exiting re-optimization loop.\n\n');
            end
            break;
        end
    

        newNonDimParams = nonDimParams;
        newNonDimParams.r0ND = X0(1:3);
        newNonDimParams.v0ND = X0(4:6);
        newNonDimParams.m0ND = X0(7);
        
        paramsX0 = [gamma, gamma2, tgo];
        
        if verboseOutput
            fprintf('\nRe-optimizing at t=%.2f s (tgo=%.2f s)...\n', t_elapsed * refVals.T_ref, tgo * refVals.T_ref);
        end
        
        ReOptTimer = tic;
        [optParams, ~, ~, ~, ~, ~, exitflag] = optimizationLoop(paramsX0, betaParam, ...
            problemParams, newNonDimParams, optimizationParams, refVals, delta_t, false, false);
        ReOptTime = toc(ReOptTimer);
        fprintf("ReOpt time: %.3fs\n",ReOptTime);
        if exitflag > 0
            gamma = optParams(1);
            gamma2 = optParams(2);
            kr = (gamma2+2)*(gamma+2);
            tgo = optParams(3);
        else
            % OPTIMIZATION FAILED
            if verboseOutput
                fprintf('WARNING: Optimization failed (exitflag=%d). Keeping previous parameters and natural tgo.\n', exitflag);
            end
        end
        
        optHistory = [optHistory; [t_elapsed, gamma, gamma2, kr, tgo]];
        exitFlags = [exitFlags; exitflag];
        
        if verboseOutput
            fprintf('New parameters: gamma=%.4f, gamma2=%.4f, kr=%.4f, tgo=%.2f s (exitflag=%d)\n\n', ...
                gamma, gamma2, kr, tgo * refVals.T_ref, exitflag);
        end
        if divertOccured
            if verboseOutput
                fprintf('Divert occurred, stopping reoptimization.\n');
            end
            break;
        end
    end
    
    %% Main loop concluded, and updateStop reached
    if tgo > 0  % safety check if there is time to do final segment
        segmentCount = segmentCount + 1;
    
        if verboseOutput
            fprintf('--- Final Segment %d (No Re-opt) ---\n', segmentCount);
            fprintf('Integrating final %.2f s to landing\n', tgo * refVals.T_ref);
        end

        t_start = t_elapsed;
        t_end    = t_elapsed + tgo;
        t_freeze = t_end - minTime;
        if t_freeze < t_start
            t_freeze = t_start;
        end

        t_reopt_start1 = t_elapsed;
        tgo_at_reopt1  = tgo;
        gamma1 = gamma;
        kr1    = kr;
        rfStar1 = rfStar;
    
        % Integrate to freeze, maybe divert
        if divertEnabled && ~divertOccured
            odeoptions = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',eventDivert);
            [tSeg1, stateSeg1, tEvent1, stateEvent1, ~] = ode45(@(t, X) trajectorySegmentODE( ...
                t, X, t_elapsed, gamma, kr, tgo, isp, rMoonND, rfStar, vfStar, afStar, gGuidance, nonDimParams), ...
                [t_start, t_freeze], X0, odeoptions);
        else
            odeoptions = odeset('RelTol',1e-6,'AbsTol',1e-6);
            [tSeg1, stateSeg1] = ode45(@(t, X) trajectorySegmentODE( ...
                t, X, t_elapsed, gamma, kr, tgo, isp, rMoonND, rfStar, vfStar, afStar, gGuidance, nonDimParams), ...
                [t_start, t_freeze], X0, odeoptions);
            tEvent1 = [];
            stateEvent1 = [];
        end
        n1 = length(tSeg1);
        segAT1 = zeros(3,n1);
        for i = 1:n1
            r = stateSeg1(i,1:3)';  v = stateSeg1(i,4:6)';  m = stateSeg1(i,7);
            tgo_now = tgo_at_reopt1 - (tSeg1(i) - t_reopt_start1);
            [aTi, limited] = computeAccel(r, v, m, gamma1, kr1, tgo_now, rfStar1, vfStar, afStar, gGuidance, nonDimParams);
            segAT1(:,i) = aTi;
            flag_thrustGotLimited = flag_thrustGotLimited || limited;
        end
        % Default: no post-divert segment
        tSeg2 = [];
        stateSeg2 = [];
        segAT2 = [];
    
        % If divert, reopt and update end/freeze
        if ~isempty(tEvent1)
            t_divert = tEvent1(end);
            X_divert = stateEvent1(end, :)';
    
            if verboseOutput
                fprintf('= Divert Triggered at t=%.2f s (final segment) =\n', t_divert * refVals.T_ref);
            end
    
            % Switch target
            rfStar = divertPoint;
            nonDimParams.rfStarND = divertPoint;
            divertOccured = true;
    
            % Update clocks to the divert time
            dt = t_divert - t_start;
            t_elapsed = t_divert;
            tgo = tgo - dt;
    
            % Re-optimize at divert
            newNonDimParams = nonDimParams;
            newNonDimParams.r0ND = X_divert(1:3);
            newNonDimParams.v0ND = X_divert(4:6);
            newNonDimParams.m0ND = X_divert(7);
    
            paramsX0 = [gamma, gamma2, tgo];
    
            [optParams, ~, ~, ~, ~, ~, exitflag] = optimizationLoop(paramsX0, betaParam, ...
                problemParams, newNonDimParams, optimizationParams, refVals, delta_t, false, false);
    
            if exitflag > 0
                gamma  = optParams(1);
                gamma2 = optParams(2);
                kr     = (gamma2+2)*(gamma+2);
                tgo    = optParams(3);
            else
                if verboseOutput
                    fprintf('WARNING: Optimization failed at divert (exitflag=%d). Keeping previous parameters and natural tgo.\n', exitflag);
                end
            end
    
            optHistory = [optHistory; [t_elapsed, gamma, gamma2, kr, tgo]];
            exitFlags  = [exitFlags; exitflag];

            t_end    = t_elapsed + tgo;
            t_freeze = t_end - minTime;
            if t_freeze < t_elapsed
                t_freeze = t_elapsed;
            end

            % Integrate post-divert until the updated freeze time
            if t_freeze > t_elapsed
                odeoptions = odeset('RelTol',1e-6,'AbsTol',1e-6);
                t_reopt_start2 = t_elapsed;
                tgo_at_reopt2  = tgo;
                gamma2_use = gamma;
                kr2_use    = kr;
                rfStar2    = rfStar;
                [tSeg2, stateSeg2] = ode45(@(t, X) trajectorySegmentODE( ...
                    t, X, t_elapsed, gamma, kr, tgo, isp, rMoonND, rfStar, vfStar, afStar, gGuidance, nonDimParams), ...
                    [t_elapsed, t_freeze], X_divert, odeoptions);
                n2 = length(tSeg2);
                segAT2_full = zeros(3,n2);
                for i = 1:n2
                    r = stateSeg2(i,1:3)';  v = stateSeg2(i,4:6)';  m = stateSeg2(i,7);
                    tgo_now = tgo_at_reopt2 - (tSeg2(i) - t_reopt_start2);
                    [aTi, limited] = computeAccel(r, v, m, gamma2_use, kr2_use, tgo_now, rfStar2, vfStar, afStar, gGuidance, nonDimParams);
                    segAT2_full(:,i) = aTi;
                    flag_thrustGotLimited = flag_thrustGotLimited || limited;
                end
                if n2 >= 2
                    segAT2 = segAT2_full(:,2:end);
                else
                    segAT2 = zeros(3,0);
                end
            end
        end
    
        % Determine freeze state
        if ~isempty(stateSeg2)
            freezeState = stateSeg2(end, :)';
            t_freeze_actual = tSeg2(end);
        else
            freezeState = stateSeg1(end, :)';
            t_freeze_actual = tSeg1(end);
        end
    
        % Compute frozen thrust at minTime
        r_freeze = freezeState(1:3);
        v_freeze = freezeState(4:6);
        m_freeze = freezeState(7);
    
        [aT_frozen, ~] = computeAccel(r_freeze, v_freeze, m_freeze, gamma, kr, ...
            minTime, rfStar, vfStar, afStar, gGuidance, nonDimParams);
    
        % Frozen thrust from freeze to end time
        odeoptions = odeset('RelTol',1e-6,'AbsTol',1e-6);
        if t_end > t_freeze_actual
            [tSeg3, stateSeg3] = ode45(@(t, X) trajectoryFrozenThrust(t, X, aT_frozen, isp, rMoonND), ...
                [t_freeze_actual, t_end], freezeState, odeoptions);
        else
            tSeg3 = t_freeze_actual;
            stateSeg3 = freezeState.';
        end
        n3 = length(tSeg3);
        if n3 >= 2
            segAT3 = repmat(aT_frozen, 1, n3-1);
        else
            segAT3 = zeros(3,0);
        end
        % Combine segments (avoid duplicates)
        tSeg = tSeg1;
        stateSeg = stateSeg1;
        segAT = segAT1;
    
        if ~isempty(tSeg2)
            tSeg = [tSeg; tSeg2(2:end)];
            stateSeg = [stateSeg; stateSeg2(2:end,:)];
            segAT = [segAT, segAT2];
        end
    
        if ~isempty(tSeg3)
            tSeg = [tSeg; tSeg3(2:end)];
            stateSeg = [stateSeg; stateSeg3(2:end,:)];
            segAT = [segAT, segAT3];
        end
    
        tTraj = [tTraj; tSeg];
        stateTraj = [stateTraj; stateSeg];
        aTList = [aTList, segAT];

        if verboseOutput
            fprintf('Final segment complete.\n\n');
        end
    end

    if verboseOutput
        fprintf('=== Simulation Complete ===\n');
        fprintf('Total elapsed time: %.2f s\n', tTraj(end) * refVals.T_ref);
        fprintf('Total segments: %d\n', segmentCount);
        fprintf('Re-optimizations: %d\n', size(optHistory, 1) - 1);
        fprintf('Thrust limited: %s\n', string(flag_thrustGotLimited));
        fprintf('===========================\n');
    end
end


function [aTi, limited] = computeAccel(r, v, m, gamma, kr, tgo, rfStar, vfStar, afStar, gGuidance, nonDimParams)
    % Compute thrust acceleration with limits
    minAccel = nonDimParams.minThrustND / m;
    maxAccel = nonDimParams.maxThrustND / m;
    
    % Guidance law
    aT1 = gamma * (kr/(2*gamma + 4) - 1) * afStar;
    aT2 = (gamma*kr/(2*gamma+4) - gamma - 1) * gGuidance;
    aT3 = ((gamma+1)/tgo) * (1 - kr/(gamma+2)) * (vfStar - v);
    aT4 = (kr/tgo^2) * (rfStar - r - v*tgo);
    aTi = aT1 + aT2 + aT3 + aT4;
    
    % Apply thrust limits
    normAT = norm(aTi);
    if normAT > maxAccel
        aTi = aTi / normAT * maxAccel;
        limited = true;
    elseif normAT < minAccel
        aTi = aTi / normAT * minAccel;
        limited = true;
    else
        limited = false;
    end
end

function dXdt = trajectorySegmentODE(t, X, t_reopt_start, gamma, kr, tgo_at_reopt, isp, rMoonND, rfStar, vfStar, afStar, gGuidance, nonDimParams)
    % ODE function for trajectory segment
    % t - current time (elapsed from simulation start)
    % t_reopt_start - time when current parameters were set
    % tgo_at_reopt - time-to-go when current parameters were set
    
    r = X(1:3);
    v = X(4:6);
    m = X(7);
    
    % Calculate current time-to-go
    t_since_reopt = t - t_reopt_start;
    tgo = tgo_at_reopt - t_since_reopt;
    
    % Ensure tgo doesn't go negative
    if tgo < 0.001
        tgo = 0.001;
    end
    
    % Gravitational acceleration
    g = -(rMoonND^2) * r / (norm(r)^3);
    
    % Compute thrust acceleration
    [aT, ~] = computeAccel(r, v, m, gamma, kr, tgo, rfStar, vfStar, afStar, gGuidance, nonDimParams);
    
    % Thrust magnitude and mass flow rate
    F_mag = norm(aT) * m;
    dm_dt = -F_mag / isp;
    
    % State derivatives
    dXdt = [v; aT + g; dm_dt];
end

function dXdt = trajectoryFrozenThrust(t, X, aT_frozen, isp, rMoonND)
    % ODE function for trajectory with frozen thrust acceleration
    % Used to avoid singularity when tgo approaches zero
    
    r = X(1:3);
    v = X(4:6);
    m = X(7);
    
    g = -(rMoonND^2) * r / (norm(r)^3);
    aT = aT_frozen;

    F_mag = norm(aT) * m;
    dm_dt = -F_mag / isp;

    dXdt = [v; aT + g; dm_dt];
end

function [value, isterminal, direction] = divertTrigger(~, y, rMoonND, altDivertND)
    r = y(1:3);
    currentAlt = norm(r) - rMoonND;
    value = currentAlt - altDivertND;
    isterminal = 1;
    direction = -1;
end