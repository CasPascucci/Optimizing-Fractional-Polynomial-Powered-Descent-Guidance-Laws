function [gammaOpt, gamma2Opt, krOpt, tgoOpt, aTOptim, exitflag, optFuelCost, simFuelCost, aTSim, finalPosSim, optHistory, ICstates, exitFlags, problemParams, nonDimParams, refVals, optTable, simTable, nActiveConstraints] = ...
    getParams(PDIState, planetaryParams, targetState, vehicleParams, optimizationParams, betaParam, doPlots, verboseOutput, dispersion, runSimulation, monteCarloSeed)
    
    if nargin > 10
        monteCarlo = true; % 11th arg in is only for Accel Monte Carlo Sim
    else
        monteCarlo = false;
        monteCarloSeed = [];
    end

    %% 1. Initialization & Defaults
    % Planetary Constants
    rPlanet = planetaryParams.rPlanet;
    gPlanet = planetaryParams.gPlanet;

    % Reference & Scaling Constants
    gEarth = 9.81;
    L_ref = 10000;
    T_ref = sqrt(L_ref/gPlanet);
    A_ref = gPlanet;
    V_ref = L_ref / T_ref;
    M_ref = vehicleParams.massInit;

    % Extract State & Vehicle Parameters
    altitude_km        = PDIState.altitude / 1000;
    lonInitDeg         = PDIState.lonInitDeg;
    latInitDeg         = PDIState.latInitDeg;
    inertialVelocity   = PDIState.inertialVelocity;
    flightPathAngleDeg = PDIState.flightPathAngleDeg;
    azimuth            = PDIState.azimuth * pi / 180; 

    landingLonDeg      = targetState.landingLonDeg;
    landingLatDeg      = targetState.landingLatDeg;
    rfLanding          = targetState.rfLanding;
    vfLanding          = targetState.vfLanding;
    afLanding          = targetState.afLanding;
    delta_t            = targetState.delta_t;
    divertEnabled      = targetState.divertEnabled & optimizationParams.updateOpt;
    altDivert          = targetState.altDivert;
    divertPoints       = targetState.divertPoints;

    massInit           = vehicleParams.massInit;
    dryMass            = vehicleParams.dryMass;
    isp                = vehicleParams.isp;
    maxThrust          = vehicleParams.maxThrust;
    minThrust          = vehicleParams.minThrust;

    %% 2. Coordinate Transformation & Non-Dimensionalization
    [r0Dim, v0Dim] = PDI2MCMF(altitude_km, lonInitDeg, latInitDeg, ...
                                       landingLonDeg, landingLatDeg, ...
                                       inertialVelocity, flightPathAngleDeg, azimuth, rPlanet);

    rfDim = 10000 * ENU2MCMF(rfLanding/10000, landingLatDeg, landingLonDeg, true);
    vfDim = ENU2MCMF(vfLanding, landingLatDeg, landingLonDeg, false);
    afDim = ENU2MCMF(afLanding, landingLatDeg, landingLonDeg, false);

    rPlanetND   = rPlanet / L_ref;
    r0ND        = r0Dim / L_ref;
    v0ND        = v0Dim / V_ref;
    rfStarND    = rfDim / L_ref;
    vfStarND    = vfDim / V_ref;
    afStarND    = afDim / A_ref;
    gConst      = -(rPlanetND^2) * rfStarND / (norm(rfStarND)^3); 
    massInitND  = massInit / M_ref;
    dryMassND   = dryMass / M_ref;
    ispND       = isp * gEarth / (V_ref);
    maxThrustND = maxThrust / (M_ref * A_ref);
    minThrustND = minThrust / (M_ref * A_ref);
    delta_tND   = delta_t / T_ref;

    %% 3. Structs Packing
    problemParams = struct;
    problemParams.r0Dim = r0Dim; 
    problemParams.v0Dim = v0Dim; 
    problemParams.rMoon = rPlanet; 
    problemParams.gMoon = gPlanet; 
    problemParams.g0 = gEarth; 
    problemParams.rfDim = rfDim; 
    problemParams.vfDim = vfDim; 
    problemParams.afDim = afDim; 
    problemParams.massInitDim = massInit; 
    problemParams.dryMassDim = dryMass; 
    problemParams.ispDim = isp; 
    problemParams.maxThrustDim = maxThrust; 
    problemParams.minThrustDim = minThrust; 
    problemParams.landingLatDeg = landingLatDeg;
    problemParams.landingLonDeg = landingLonDeg;
    problemParams.divertEnabled = divertEnabled;
    problemParams.altDivert = altDivert;
    problemParams.divertPoints = divertPoints;

    refVals = struct;
    refVals.L_ref = L_ref;
    refVals.T_ref = T_ref;
    refVals.A_ref = A_ref;
    refVals.V_ref = V_ref;
    refVals.M_ref = M_ref;

    nonDimParams = struct;
    nonDimParams.rMoonND = rPlanetND;
    nonDimParams.r0ND = r0ND;
    nonDimParams.v0ND = v0ND;
    nonDimParams.rfStarND = rfStarND;
    nonDimParams.vfStarND = vfStarND;
    nonDimParams.afStarND = afStarND;
    nonDimParams.gConst = gConst;
    nonDimParams.m0ND = massInitND;
    nonDimParams.mMinND = dryMassND;
    nonDimParams.ispND = ispND;
    nonDimParams.maxThrustND = maxThrustND;
    nonDimParams.minThrustND = minThrustND;

    %% 4. Optimization
    
    paramsX0 = optimizationParams.paramsX0;
    paramsX0(3) = paramsX0(3)/refVals.T_ref;
    reopt = optimizationParams.updateOpt;
    fprintf("=== Starting Optimization ===\n");
    OptTimer = tic;
    [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag, nActiveConstraints] = ...
        optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimizationParams, refVals, delta_tND, verboseOutput, dispersion);
    optTime = toc(OptTimer);
    fprintf("Optimization Time: %.3fs\n", optTime);
    if exitflag ~= 1
        fprintf("\n First Optimization Converged to flag ~=1, rerunning optimization starting from first rounds parameters:\n");
        [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag, nActiveConstraints] = ...
        optimizationLoop(optParams, betaParam, problemParams, nonDimParams, optimizationParams, refVals, delta_tND, verboseOutput, dispersion);
    end

    gammaOpt = optParams(1);
    gamma2Opt = optParams(2);
    krOpt = (gamma2Opt+2)*(gammaOpt+2);
    tgoOpt = optParams(3) * T_ref;
    optFuelCost = (mOptim(end)-mOptim(1))*M_ref;

    %% 6. Simulation
    tTraj = [];
    stateTraj = [];
    aTSim = [];
    flag_thrustGotLimited = false;
    simFuelCost = [];
    finalPosSim = [];
    optHistory = [];
    ICstates = [];
    exitFlags = [];
    if monteCarlo
        [tTraj, stateTraj, aTSim, flag_thrustGotLimited] = ...
            closedLoopSim(gammaOpt, gamma2Opt, tgoOpt/T_ref, nonDimParams, refVals, monteCarloSeed);
        if ~isempty(optHistory)
            optHistory = array2table(optHistory);
            optHistory.Properties.VariableNames(1:5) = {'t_elapsedND','gamma1','gamma2','kr','tgoND'};
            ICstates = array2table(ICstates');
            ICstates.Properties.VariableNames(1:7) = {'r0_X', 'r0_Y', 'r0_Z', 'v0_X', 'v0_Y', 'v0_Z', 'm0'};
        end
        simFuelCost = M_ref * (stateTraj(1,7) - stateTraj(end,7));
        finalPosSim = MCMF2ENU(stateTraj(end,1:3)' * L_ref, landingLatDeg, landingLonDeg, true, true);
    elseif runSimulation % Only enter if simulation is required
        if ~reopt % If not reoptimizing call this variant
            % Static Simulation
            [tTraj, stateTraj, aTSim, flag_thrustGotLimited] = ...
                closedLoopSim(gammaOpt, gamma2Opt, tgoOpt/T_ref, nonDimParams, refVals);
            [c1_static, c2_static, c1_num, c2_num] = calculateCoeffs(nonDimParams.r0ND, nonDimParams.v0ND, optParams(3), ...
                                       gammaOpt, gamma2Opt, nonDimParams.afStarND, ...
                                       nonDimParams.rfStarND, nonDimParams.vfStarND, nonDimParams.gConst);
            optHistory = [c1_static, c2_static, c1_num, c2_num];
        else % If reoptimizing, call this one
            divertTrajectories = cell(1,size(problemParams.divertPoints, 1));
            if divertEnabled % Re-Optimization with Divert, as reopt is required for divert
                for idx = 1:size(problemParams.divertPoints, 1)
                    divertPoint = problemParams.divertPoints(idx,:)' ./ refVals.L_ref;
                    divertPoint = ENU2MCMF(divertPoint, landingLatDeg, landingLonDeg, true);
                    [tTraj, stateTraj, aTSim, flag_thrustGotLimited, optHistory, ICstates, exitFlags] = ...
                    simReOpt(gammaOpt, gamma2Opt, tgoOpt/T_ref, problemParams, nonDimParams, refVals, delta_tND, optimizationParams, betaParam, verboseOutput, divertPoint);
                    divertTrajectories{idx} = stateTraj(:,1:3);
                end
            else
                % regular Re-Optimization Simulation
                [tTraj, stateTraj, aTSim, flag_thrustGotLimited, optHistory, ICstates, exitFlags] = ...
                    simReOpt(gammaOpt, gamma2Opt, tgoOpt/T_ref, problemParams, nonDimParams, refVals, delta_tND, optimizationParams, betaParam, verboseOutput);
            end



            
            % Format History Output
            if ~isempty(optHistory)
                optHistory = array2table(optHistory);
                optHistory.Properties.VariableNames(1:5) = {'t_elapsedND','gamma1','gamma2','kr','tgoND'};
                ICstates = array2table(ICstates');
                ICstates.Properties.VariableNames(1:7) = {'r0_X', 'r0_Y', 'r0_Z', 'v0_X', 'v0_Y', 'v0_Z', 'm0'};
            end
        end
        
        simFuelCost = M_ref * (stateTraj(1,7) - stateTraj(end,7));
        finalPosSim = MCMF2ENU(stateTraj(end,1:3)' * L_ref, landingLatDeg, landingLonDeg, true, true);
    end

    if flag_thrustGotLimited
        fprintf("Simulation Trajectory Thrust Limited!\n")
    end

    % Get Pseudo cost from Simulation
    aTmagSim = vecnorm(aTSim,2,1);
    simCost = betaParam* trapz(tTraj,aTmagSim') + (1-betaParam)*trapz(tTraj, dot(aTSim,aTSim)');
    %% 7. Plotting
    if doPlots
        plotting(tTraj, stateTraj, optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, aTSim, ...
            refVals, problemParams, nonDimParams, optimizationParams, flag_thrustGotLimited, ...
            optHistory, ICstates, betaParam);
    end

    %% 8. Divert Plot
    if divertEnabled
        figure(); hold on;
        for idx = 1:size(problemParams.divertPoints, 1)
            traj = divertTrajectories{idx};
            plot3(traj(:,1),traj(:,2),traj(:,3))
        end
        xlabel("X");ylabel("Y");zlabel("Z");
        xlim([-0.2,0.2]);ylim([-0.2,0.2]);zlim([-173.65,-173.6]);

        figure(); hold on;
        for idx = 1:size(problemParams.divertPoints, 1)
            traj = divertTrajectories{idx} * refVals.L_ref;
            plot(traj(end,1),traj(end,2),'.','MarkerSize',12);
        end
        xlabel("X");ylabel("Y");
    end
optError = MCMF2ENU(rdOptim(:,1), landingLatDeg,landingLonDeg,true,false);
optErrorNorm = norm(rfLanding*refVals.L_ref - optError*refVals.L_ref);
optTable = [optErrorNorm; optFuelCost; optCost];
if runSimulation
    simError = MCMF2ENU(stateTraj(end,1:3)', landingLatDeg,landingLonDeg,true,false);
    simErrorNorm = norm(rfLanding*refVals.L_ref - simError*refVals.L_ref);
    simTable = [simErrorNorm; simFuelCost; simCost];
end
end