function outputSingle = reOptReRun(idx, ICStates, optHistory, betaParam, problemParams, nonDimParams, optimizationParams, refVals, delta_tND, verboseOutput)
       ICStates = table2array(ICStates);
       optHistory = table2array(optHistory);

       % --- synchronization Logic ---
       % ICStates(idx) is the vehicle state at the START of the segment.
       % optHistory(idx) are the guidance parameters active for that segment.
       
       IC = ICStates(idx,:);
       
       % Setup Non-Dim Params for this specific instant
       nonDimParams.r0ND = IC(1:3)';
       nonDimParams.v0ND = IC(4:6)';
       nonDimParams.m0ND = IC(7);
       
       X0 = [nonDimParams.r0ND; nonDimParams.v0ND; nonDimParams.m0ND];
       
       % Extract the parameters that drove this segment
       tStart = optHistory(idx,1);
       gamma  = optHistory(idx,2); 
       gamma2 = optHistory(idx,3); 
       kr     = optHistory(idx,4);
       tgo    = optHistory(idx,5);
       
       updateFreqND = optimizationParams.updateFreq / refVals.T_ref;
       tEnd = tStart + updateFreqND;

       % Standard Constants
       rfStar = nonDimParams.rfStarND;
       vfStar = nonDimParams.vfStarND;
       afStar = nonDimParams.afStarND;
       isp = nonDimParams.ispND;
       rMoonND = nonDimParams.rMoonND;
       gGuidance = nonDimParams.gConst;
       flag_thrustGotLimited = false;

       % --- RECONSTRUCT PLAN (Do not re-optimize) ---
       % We use the EXACT parameters from history to see what the guidance *thought* it was doing.
       % This prevents the "New Parameters vs Old State" off-by-one error.
       
       [c1, c2] = calculateCoeffs(nonDimParams.r0ND, nonDimParams.v0ND, tgo, gamma, gamma2, afStar, rfStar, vfStar, gGuidance);
       
       % Generate the arrays for the "Planned" trajectory
       nodeCount = optimizationParams.nodeCount;
       tgospan = linspace(0, tgo, nodeCount);
       
       % Calculate Coeffs functions (inline replication of optimizationLoop logic)
       phi1hat = (tgospan.^(gamma+2))./((gamma+1)*(gamma+2));
       phi2hat = (tgospan.^(gamma2+2))./((gamma2+1)*(gamma2+2));
       phi1bar = (tgospan.^(gamma+1))./(gamma+1);
       phi2bar = (tgospan.^(gamma2+1))./(gamma2+1);
       
       rdOptim = rfStar + c1*phi1hat + c2*phi2hat - vfStar.*tgospan + 0.5*(gGuidance+afStar).*tgospan.^2;
       vdOptim = vfStar - c1*phi1bar - c2*phi2bar -(gGuidance+afStar).*tgospan;
       aTOptim = afStar + c1*tgospan.^gamma + c2*tgospan.^gamma2;
       
       % Reconstruct Mass Plan
       aTNorm = vecnorm(aTOptim, 2, 1);
       Q = cumtrapz(tgospan, aTNorm./isp);
       Q = Q(end) - Q;
       mOptim = nonDimParams.m0ND .* exp(-Q);
       
       % Pack parameters for plotting
       optParams = [gamma, gamma2, tgo];
       optCost = 0; % Irrelevant here

       % --- RUN SIMULATION (The "Flown" Segment) ---
       odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
       [tSeg, stateSeg] = ode45(@(t, X) trajectorySegmentODE(t, X, tStart, gamma, kr, tgo, ...
            isp, rMoonND, rfStar, vfStar, afStar, gGuidance, nonDimParams), ...
            [tStart, tEnd], X0, odeoptions);

       numPoints = length(tSeg);
       segAT = zeros(3,numPoints);
        for i = 1:numPoints
            r = stateSeg(i,1:3)';
            v = stateSeg(i,4:6)';
            m = stateSeg(i,7);
            t_since_reopt = tSeg(i) - tStart;
            tgo_now = tgo - t_since_reopt;
            [aTi, limited] = computeAccel(r, v, m, gamma, kr, tgo_now, rfStar, vfStar, afStar, gGuidance, nonDimParams);
            segAT(:, i) = aTi;
            flag_thrustGotLimited = flag_thrustGotLimited || limited;
        end

       % --- PLOT ---
       plotReRunSegment(tSeg, stateSeg, optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, segAT, refVals, problemParams, nonDimParams, optimizationParams, flag_thrustGotLimited)
       
       % --- RETURN DATA ---
       outputSingle.tSeg = tSeg;
       outputSingle.stateSeg = stateSeg;
       outputSingle.optParams = optParams;
       outputSingle.optCost = optCost;
       outputSingle.aTOptim = aTOptim;
       outputSingle.mOptim = mOptim;
       outputSingle.rdOptim = rdOptim;
       outputSingle.vdOptim = vdOptim;
       outputSingle.segAT = segAT;
       outputSingle.flag_thrustGotLimited = flag_thrustGotLimited;

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