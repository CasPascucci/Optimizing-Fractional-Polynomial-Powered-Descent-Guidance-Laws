function [tTraj, stateTraj, aTList, flag_thrustGotLimited] = closedLoopSim(gamma, gamma2, tgo0, nonDimParams, refVals, monteCarloSeed)

    if nargin < 6
        monteCarloSeed = [];
    end

    kr = (gamma2+2)*(gamma+2);

    r0      = nonDimParams.r0ND;
    v0      = nonDimParams.v0ND;
    m0      = nonDimParams.m0ND;
    rMoonND = nonDimParams.rMoonND;
    isp     = nonDimParams.ispND;
    rfStar  = nonDimParams.rfStarND;
    vfStar  = nonDimParams.vfStarND;
    afStar  = nonDimParams.afStarND;
    gConst  = nonDimParams.gConst;

    X0       = [r0; v0; m0];
    minTime  = 0.2 / refVals.T_ref;          % freeze threshold (ND)
    t_freeze = max(0, tgo0 - minTime);
    odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

    % Phase 1: guidance active (t = 0 to t_freeze)
    [tPhase1, statePhase1] = ode45( ...
        @(t, X) guidanceODE(t, X, gamma, kr, tgo0, isp, rMoonND, ...
                            rfStar, vfStar, afStar, gConst, nonDimParams, monteCarloSeed), ...
        [0, t_freeze], X0, odeoptions);

    % Compute frozen thrust at the freeze point (tgo = minTime)
    X_freeze = statePhase1(end, :)';
    [aT_frozen, flag_thrustGotLimited] = computeAccel( ...
        X_freeze(1:3), X_freeze(4:6), X_freeze(7), gamma, kr, minTime, ...
        rfStar, vfStar, afStar, gConst, nonDimParams, monteCarloSeed);

    % Phase 2: frozen thrust (t = t_freeze to tgo0)
    [tPhase2, statePhase2] = ode45( ...
        @(t, X) frozenThrustODE(t, X, aT_frozen, isp, rMoonND), ...
        [t_freeze, tgo0], X_freeze, odeoptions);

    % Combine phases
    tTraj     = [tPhase1;          tPhase2(2:end)];
    stateTraj = [statePhase1;      statePhase2(2:end, :)];
    numStates = size(stateTraj, 1);
    aTList    = zeros(3, numStates);

    for idx = 1:numStates
        tgo = tgo0 - tTraj(idx);
        if tgo > minTime
            [aTi, limited] = computeAccel( ...
                stateTraj(idx,1:3)', stateTraj(idx,4:6)', stateTraj(idx,7), ...
                gamma, kr, tgo, rfStar, vfStar, afStar, gConst, nonDimParams, monteCarloSeed);
            flag_thrustGotLimited = flag_thrustGotLimited || limited;
        else
            aTi = aT_frozen;
        end
        aTList(:, idx) = aTi;
    end
end


%% Local functions

function dXdt = guidanceODE(t, X, gamma, kr, tgo0, isp, rMoonND, rfStar, vfStar, afStar, gConst, nonDimParams, monteCarloSeed)
    r    = X(1:3);
    v    = X(4:6);
    mass = X(7);
    tgo  = tgo0 - t;
    g    = -(rMoonND^2) * r / (norm(r)^3);
    [aT, ~] = computeAccel(r, v, mass, gamma, kr, tgo, rfStar, vfStar, afStar, gConst, nonDimParams, monteCarloSeed);
    dXdt = [v; aT + g; -norm(aT)*mass/isp];
end

function dXdt = frozenThrustODE(~, X, aT_frozen, isp, rMoonND)
    r    = X(1:3);
    v    = X(4:6);
    mass = X(7);
    g    = -(rMoonND^2) * r / (norm(r)^3);
    dXdt = [v; aT_frozen + g; -norm(aT_frozen)*mass/isp];
end

function [aT, limited] = computeAccel(r, v, m, gamma, kr, tgo, rfStar, vfStar, afStar, gConst, nonDimParams, monteCarloSeed)
    minAccel = nonDimParams.minThrustND / m;
    maxAccel = nonDimParams.maxThrustND / m;
    aT = gamma*(kr/(2*gamma+4)-1)*afStar ...
       + (gamma*kr/(2*gamma+4)-gamma-1)*gConst ...
       + ((gamma+1)/tgo)*(1-kr/(gamma+2))*(vfStar-v) ...
       + (kr/tgo^2)*(rfStar-r-v*tgo);
    if ~isempty(monteCarloSeed)
        aT = aT * monteCarloSeed;
    end
    n = norm(aT);
    if n > maxAccel
        aT = aT / n * maxAccel;
        limited = true;
    elseif n < minAccel
        aT = aT / n * minAccel;
        limited = true;
    else
        limited = false;
    end
end
