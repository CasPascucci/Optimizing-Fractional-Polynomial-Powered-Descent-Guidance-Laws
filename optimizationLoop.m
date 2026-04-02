function [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag] = optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimizationParams, refVals, delta_t, verboseOutput, dispersion)

    % Some Tweakable Parameters / Defaults
    if optimizationParams.gamma1eps < 0
        optimizationParams.gamma1eps = 1e-2;
    end
    if optimizationParams.gamma2eps < 0
        optimizationParams.gamma2eps = 1e-2;
    end
    % Initiailizing
    r0 = nonDimParams.r0ND;
    v0 = nonDimParams.v0ND;

    rfStar = nonDimParams.rfStarND;
    vfStar = nonDimParams.vfStarND;
    afStar = nonDimParams.afStarND;

    gConst = nonDimParams.gConst;

    maxThrust = nonDimParams.maxThrustND;
    minThrust = nonDimParams.minThrustND;
    isp = nonDimParams.ispND;

    [rot_MCMF2ENU, r0_MCMF_ND] = computeMCMF2ENUTransform(...
        problemParams.landingLatDeg, ...
        problemParams.landingLonDeg, ...
        nonDimParams.rMoonND);

    % Fmincon Constraints
    Aineq = [-1  0  0;      % -gamma <= -gamma1eps ---- gamma >= gamma1eps
              1 -1  0;      %  gamma -gamma2 <= -tol ---- gamma2 >= gamma + tol
              0  0 -1];     % tgo >= 0.01
    bineq = [-optimizationParams.gamma1eps; -optimizationParams.gamma2eps; -0.01];

    lb = [optimizationParams.gamma1eps, 0, 0.01];
    ub = [5, 10, 15];
    % lb = [0.3,0.757,0.01];
    % ub = [0.3,0.757,15];

    fminconOptions = optimoptions('fmincon', 'Display', 'none', 'MaxFunctionEvaluations', 10000, ...
    'FiniteDifferenceType','forward','MaxIterations', 1000, ...
    'Algorithm','sqp', 'OptimalityTolerance',1e-4);

    nodeCount = optimizationParams.nodeCount;

    obj = @(params) objectiveFunction(params, betaParam, afStar, rfStar, r0, vfStar, v0, gConst, nonDimParams, optimizationParams);
    nonlincon = @(params) nonLinearLimits(params, r0, v0, rfStar, vfStar, afStar, gConst, isp, minThrust, maxThrust, optimizationParams, problemParams, refVals, rot_MCMF2ENU, r0_MCMF_ND, nonDimParams.m0ND);
    
    [optParams, optCost, exitflag, output] = fmincon(obj, paramsX0, Aineq, bineq, [], [], lb, ub, nonlincon, fminconOptions);
    if verboseOutput
        fprintf('=== Optimization Results ===\n');
        fprintf('Exit flag: %d\n', exitflag);
        fprintf('Exit message: %s\n', output.message);
        fprintf('Final cost: %.6e\n', optCost);
        fprintf('Iterations: %d\n', output.iterations);
        fprintf('First-order optimality: %.6e (tolerance: %.1e)\n', output.firstorderopt, fminconOptions.OptimalityTolerance);
        fprintf('\nOptimal parameters:\n');
        fprintf('  gamma = %.6f\n', optParams(1));
        fprintf('  gamma2    = %.6f\n', optParams(2));
        fprintf('  tgo   = %.6f (ND)\n', optParams(3));
    
        [c, ~] = nonlincon(optParams);
        activeIneq = find(c >= -1e-6); % Nearly active
        numActiveIneq = length(activeIneq);
        activeIneqString = strjoin(string(activeIneq),', ');
        if abs(optParams(1) - optimizationParams.gamma1eps) < 1e-3
            numActiveIneq = numActiveIneq + 1;
        end
        if abs(optParams(2) - optParams(1)) - optimizationParams.gamma1eps < 1e-3
            numActiveIneq = numActiveIneq + 1;
        end
        fprintf('\nTotal active inequality constraints: %d/%d\n', numActiveIneq, length(c));
        if ~isempty(activeIneq)
            fprintf("Constraints #: %s\n", activeIneqString);
        end
        if length(activeIneq) ~= numActiveIneq
            fprintf("Linear constraint also active\n");
        end
        [maxViolation, idx] = max(c);
        fprintf('Max violation: %.3f at constraint %d\n', maxViolation, idx);
    end

    if exitflag <= 0
            fprintf('\n=== Optimization Failed (exitflag=%d) ===\n', exitflag);
            fprintf('Params: gamma=%.4f, gamma2=%.4f, tgo=%.4f\n', optParams(1), optParams(2), optParams(3));
            [c, ~] = nonlincon(optParams);
            violated = find(c > 1e-6);
            fprintf('Violated constraints (%d total):\n', length(violated));
            for vi = 1:min(10, length(violated))
                fprintf('  c(%d) = %.4f\n', violated(vi), c(violated(vi)));
            end
            fprintf('Max violation: c(%d) = %.6f\n', find(c==max(c),1), max(c));
    end

    


    gamma1 = optParams(1);
    gamma2 = optParams(2);
    if gamma1 < 0
        gamma1 = 0;
    end
    if gamma2 < 0
        gamma2 = 0;
    end

    [c1, c2] = calculateCoeffs(r0, v0, optParams(3), gamma1, gamma2, afStar, rfStar, vfStar, gConst);

    if verboseOutput
        fprintf('\nGuidance Coefficients:\n');
        fprintf('  c1 = [%.6e, %.6e, %.6e]\n', c1(1), c1(2), c1(3));
        fprintf('  c2 = [%.6e, %.6e, %.6e]\n', c2(1), c2(2), c2(3));
    end

    tgospan = linspace(0,optParams(3),nodeCount);
    aTOptim = afStar + c1*tgospan.^gamma1 + c2*tgospan.^gamma2;

    phi1hat = (tgospan.^(gamma1+2))./((gamma1+1)*(gamma1+2));
    phi2hat = (tgospan.^(gamma2+2))./((gamma2+1)*(gamma2+2));
    phi1bar = (tgospan.^(gamma1+1))./(gamma1+1);
    phi2bar = (tgospan.^(gamma2+1))./(gamma2+1);

    rdOptim = rfStar + c1*phi1hat + c2*phi2hat - vfStar.*tgospan + 0.5*(gConst+afStar).*tgospan.^2;
    vdOptim = vfStar + c1*phi1bar + c2*phi2bar -(gConst+afStar).*tgospan;

    aTNorm = vecnorm(aTOptim,2,1);

    Q = cumtrapz(tgospan,aTNorm./isp);
    Q = Q(end) - Q;

    mOptim = 1 .* exp(-Q);

    
end
%% Functions
function cost = objectiveFunction(params, betaParam, afStar, rfStar, r, vfStar, v, gConst, nonDimParams, optimizationParams)
    gamma1  = params(1);
    gamma2     = params(2);
    tgo   = params(3);
    nodeCount = optimizationParams.nodeCount;

    if gamma1 < 0
        gamma1 = 0;
    end
    if gamma2 < 0
        gamma2 = 0;
    end


    [c1, c2] = calculateCoeffs(r, v, tgo, gamma1, gamma2, afStar, rfStar, vfStar, gConst);

    tspan = linspace(0,tgo,nodeCount);

    aT = afStar + c1*tspan.^gamma1 + c2*tspan.^gamma2;
    
    aTmag = vecnorm(aT,2,1);

    simpson1 = simpsonComp13Integral(tspan,aTmag);
    simpson2 = simpsonComp13Integral(tspan,dot(aT,aT));

    cost = betaParam*simpson1 + (1-betaParam)*simpson2;

end

function [c, ceq] = nonLinearLimits(params, r0, v0, rfStar, vfStar, afStar, gConst, isp, minThrust, maxThrust, optimizationParams, problemParams, refVals, rot_MCMF2ENU, r0_MCMF_ND, m0)
    nodeCount = optimizationParams.nodeCount;
    glideSlopeFlag = optimizationParams.glideSlopeEnabled;
    pointingFlag = optimizationParams.pointingEnabled;
    gamma1  = params(1);
    gamma2     = params(2);
    tgo0   = params(3);

    [c1, c2] = calculateCoeffs(r0, v0, tgo0, gamma1, gamma2, afStar, rfStar, vfStar, gConst); % Calculate guidance coeffs
    tgospan = linspace(0,tgo0,nodeCount); % reverse of tspan in tgo time

    aT = afStar + c1*tgospan.^gamma1 + c2*tgospan.^gamma2;
    aTmag = vecnorm(aT,2,1);

    Q = cumtrapz(tgospan,aTmag./isp); % cumulative integral of fuel expenditure
    Q = Q(end) - Q;

    m = m0 .* exp(-Q); % Convert Q actually into mass units

    thrust = m .*(aTmag);

    nThrustCon = 2 * nodeCount; % 2x nodes to include upper and lower thrust bounds
    nGlideCon = 0;
    nPointCon = 0;
    
    if glideSlopeFlag % num constraint needs to be consistent, preallocate number of constraints for glideslope and pointing
        nGlideCon = floor(nodeCount * 0.15); % bottom 15% of trajectory is analyzed
    end
    if pointingFlag
        nPointCon = nodeCount;
    end

    c = zeros(nThrustCon + nGlideCon + nPointCon, 1);
    ceq = [];

    idx = 1;
    upper = thrust - (maxThrust);
    lower = minThrust - thrust;
    c(idx:idx+nodeCount-1) = upper(:);
    idx = idx + nodeCount;
    c(idx:idx+nodeCount-1) = lower(:);
    idx = idx + nodeCount;

    % Glideslope Constraint
    % if glideSlopeFlag
    %     freeGlideNodes = floor(nodeCount*0.01) + 1;
    %     tgospanGlide = tgospan(1:nGlideCon);
    %     phi1hat = (tgospanGlide.^(gamma1+2))./((gamma1+1)*(gamma1+2));
    %     phi2hat = (tgospanGlide.^(gamma2+2))./((gamma2+1)*(gamma2+2));
    %     rdOptim = rfStar + c1*phi1hat + c2*phi2hat - vfStar.*tgospanGlide + 0.5*(gConst+afStar).*tgospanGlide.^2;
    %     rdOptimTOPO = rot_MCMF2ENU * (rdOptim - r0_MCMF_ND);
    %     %rdOptimTOPO = MCMF2ENU(rdOptim,problemParams.landingLatDeg,problemParams.landingLonDeg,true,false);
    %     rUnitVec = rdOptimTOPO./vecnorm(rdOptimTOPO,2,1);
    % 
    % 
    %     if norm(rdOptimTOPO(:,1)) < 1e-10
    %         rUnitVec(:,1) = [0;0;1];
    %     end
    %     vertUnitVec = [0;0;1];
    % 
    %     rMoon = problemParams.rMoon;
    %     rdOptimDim = rdOptim*refVals.L_ref;
    %     alt_opt = vecnorm(rdOptimDim,2,1) - rMoon;
    %     frac = (alt_opt - 250)/250;
    %     frac(frac < 0) = 0;
    %     frac(frac > 1) = 1;
    %     theta = (frac*45) + 45;
    %     theta(1:freeGlideNodes) = 90;
    % 
    %     cGlide = cosd(theta) - dot(rUnitVec, repmat(vertUnitVec, 1, nGlideCon));
    %     c(idx:idx+nGlideCon-1) = cGlide;
    %     idx = idx + nGlideCon;
    % end
    if glideSlopeFlag
        nGlideCon = nodeCount;
        tgospanGlide = tgospan(1:nGlideCon);
        phi1hat = (tgospanGlide.^(gamma1+2))./((gamma1+1)*(gamma1+2));
        phi2hat = (tgospanGlide.^(gamma2+2))./((gamma2+1)*(gamma2+2));
        rdOptim = rfStar + c1*phi1hat + c2*phi2hat - vfStar.*tgospanGlide + 0.5*(gConst+afStar).*tgospanGlide.^2;
        rdOptimTOPO = rot_MCMF2ENU * (rdOptim - r0_MCMF_ND);
        rUnitVec = rdOptimTOPO./vecnorm(rdOptimTOPO,2,1);

        if norm(rdOptimTOPO(:,1)) < 1e-10
            rUnitVec(:,1) = [0;0;1];
        end
        vertUnitVec = [0;0;1];
        
        rMoon = problemParams.rMoon;
        rdOptimDim = rdOptim*refVals.L_ref;
        alt_opt = vecnorm(rdOptimDim,2,1) - rMoon;

        theta = zeros(1,nGlideCon);

        %Below glideSlopeLow meters, glideSlopeFinalTheta
        idxLow = alt_opt <= optimizationParams.glideSlopeLow;
        theta(idxLow) = optimizationParams.glideSlopeFinalTheta;

        idxMed = (alt_opt > 250) & (alt_opt <= optimizationParams.glideSlopeHigh);
        theta(idxMed) = optimizationParams.glideSlopeFinalTheta + optimizationParams.glideSlopeFinalTheta*((alt_opt(idxMed) - optimizationParams.glideSlopeLow) / optimizationParams.glideSlopeLow);

        theta(alt_opt > optimizationParams.glideSlopeHigh) = 180;
        
        % free up bottom constraints (too narrow / infinitely small)
        theta(1:optimizationParams.freeGlideNodes) = 180;

        cGlide = cosd(theta) - dot(rUnitVec, repmat(vertUnitVec, 1, nGlideCon));
        c(idx:idx+nGlideCon-1) = cGlide;
        idx = idx + nGlideCon;
    end
    if pointingFlag 
        aTTOPO = rot_MCMF2ENU * aT;
        %aTTOPO = MCMF2ENU(aT,problemParams.landingLatDeg,problemParams.landingLonDeg,false,false);
        aTmagTOPO = vecnorm(aTTOPO,2,1);
        thrustUnitVec = aTTOPO./aTmagTOPO;
        vertUnitVec = [0;0;1];
        dotProducts = thrustUnitVec'*vertUnitVec;
        phi = acosd(dotProducts); % Values at idx 1 are landing, idx end are PDI, Deg
        phi0 = optimizationParams.minPointing;
        phiA = 0.5 * (optimizationParams.maxTiltAccel * refVals.T_ref^2) * (tgospan(:).^2); % DEG
        Theta = zeros(nodeCount,1);
        idxCon1 = (phiA + phi0) >= 180;
        Theta(idxCon1) = 180;
        idxCon2 = ~idxCon1;
        Theta(idxCon2) = phiA(idxCon2) + phi0;

        cPoint = phi - Theta;
        c(idx:idx+nodeCount-1) = cPoint;
    end
end
function [rot, r0_MCMF] = computeMCMF2ENUTransform(landingLatDeg, landingLonDeg, Rm)
    
    lat0 = deg2rad(landingLatDeg);
    lon0 = deg2rad(landingLonDeg);
    
    clat = cos(lat0); slat = sin(lat0);
    clon = cos(lon0); slon = sin(lon0);

    E0 = [-slon;       clon;      0];
    N0 = [-slat*clon; -slat*slon; clat];
    U0 = [clat*clon;   clat*slon; slat];
    
    rot = [E0, N0, U0]';

    % Landing site position in MCMF
    r0_MCMF = Rm * U0;
end