function [gammaOpt, gamma2Opt, krOpt, tgoOpt, aTOptim, exitflag, optFuelCost, simFuelCost, aTSim, finalPosSim, optHistory, ICstates, exitFlags, problemParams, nonDimParams, refVals] = ...
    getParamsDIVERT(PDIState, planetaryParams, targetState, vehicleParams, optimizationParams, betaParam, doPlots, verboseOutput, runSimulation)

    if targetState.divertEnabled % Enforce that these constraints don't apply during divert scenarios, and that reopt is enabled
        optimizationParams.glideSlopeEnabled = false;
        optimizationParams.pointingEnabled = false;
        optimizationParams.updateOpt = true;
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

    rfDim = 10000 * ENU2MCMF(rfLanding/10000, landingLatDeg, landingLonDeg, true, rPlanet/L_ref);
    vfDim = ENU2MCMF(vfLanding, landingLatDeg, landingLonDeg, false, rPlanet/L_ref);
    afDim = ENU2MCMF(afLanding, landingLatDeg, landingLonDeg, false, rPlanet/L_ref);

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
    [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag] = ...
        optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimizationParams, refVals, verboseOutput);
    if exitflag ~= 1
        fprintf("\n First Optimization Converged to flag ~=1, rerunning optimization starting from first rounds parameters:\n");
        [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag] = ...
        optimizationLoop(optParams, betaParam, problemParams, nonDimParams, optimizationParams, refVals, verboseOutput);
    end

    gammaOpt = optParams(1);
    gamma2Opt = optParams(2);
    krOpt = (gamma2Opt+2)*(gammaOpt+2);
    tgoOpt = optParams(3) * T_ref;
    optFuelCost = (mOptim(end)-mOptim(1))*M_ref;

    %% 5. Simulation
    tTraj = [];
    stateTraj = [];
    aTSim = [];
    flag_thrustGotLimited = false;
    simFuelCost = [];
    finalPosSim = [];
    optHistory = [];
    ICstates = [];
    exitFlags = [];

    divertData = cell(1,size(problemParams.divertPoints, 1));


    for idx = 1:size(problemParams.divertPoints, 1)
        divertPoint = problemParams.divertPoints(idx,:)' ./ refVals.L_ref;
        divertPoint = ENU2MCMF(divertPoint, landingLatDeg, landingLonDeg, true, rPlanetND);
        [tTraj, stateTraj, aTSim, flag_thrustGotLimited, optHistory, ICstates, exitFlags] = ...
        simReOpt(gammaOpt, gamma2Opt, tgoOpt/T_ref, problemParams, nonDimParams, refVals, optimizationParams, betaParam, verboseOutput, divertPoint);

        divertData{idx}.tTraj = tTraj;
        divertData{idx}.stateTraj = stateTraj;
        divertData{idx}.aTSim = aTSim;
        divertData{idx}.optHistory = optHistory;
    end

    
    % Format History Output
    if ~isempty(optHistory)
        optHistory = array2table(optHistory);
        optHistory.Properties.VariableNames(1:5) = {'t_elapsedND','gamma1','gamma2','kr','tgoND'};
        ICstates = array2table(ICstates');
        ICstates.Properties.VariableNames(1:7) = {'r0_X', 'r0_Y', 'r0_Z', 'v0_X', 'v0_Y', 'v0_Z', 'm0'};
    end
    
    simFuelCost = M_ref * (stateTraj(1,7) - stateTraj(end,7));
    finalPosSim = MCMF2ENU(stateTraj(end,1:3)' * L_ref, landingLatDeg, landingLonDeg, true, rPlanet);


    %% 6. Divert Plot
    if divertEnabled
        % Plot 3D trajectories from 1.5x divert altitude to landing
        figure('Name', 'Divert Trajectories - 3D'); 
        hold on; grid on;
        
        altThreshold = 1.5 * altDivert; % m (dimensional)
        altThresholdND = altThreshold / refVals.L_ref; % Non-dimensional
        
        colors = lines(size(problemParams.divertPoints, 1));
        legendEntries = {};
        
        for idx = 1:size(problemParams.divertPoints, 1)
            traj = divertData{idx}.stateTraj(:,1:3);
            
            % Calculate altitude for each point
            altitudes = vecnorm(traj, 2, 2) - rPlanetND;
            
            % Find indices below threshold altitude
            belowThreshold = altitudes <= altThresholdND;
            
            if any(belowThreshold)
                trajFiltered = traj(belowThreshold, :);
                trajENU = zeros(size(trajFiltered));
                for i = 1:size(trajFiltered, 1)
                    rMCMF = trajFiltered(i,:)' * refVals.L_ref;
                    rENU = MCMF2ENU(rMCMF, landingLatDeg, landingLonDeg, true, rPlanet);
                    trajENU(i, :) = rENU';
                end
                if idx == 1
                    plot3(trajENU(:,1), trajENU(:,2), trajENU(:,3), ...
                     'LineStyle','--', 'LineWidth', 2, 'Color', colors(idx,:));
                else
                    plot3(trajENU(:,1), trajENU(:,2), trajENU(:,3), ...
                    'LineWidth', 2, 'Color', colors(idx,:));
                end
                legendEntries{end+1} = sprintf('Point %d: E=%.0fm, N=%.0fm', ...
                    idx, problemParams.divertPoints(idx,1), problemParams.divertPoints(idx,2));
            end
        end
        
        xlabel('East (m)'); ylabel('North (m)'); zlabel('Up (m)');
        title(sprintf('Divert Trajectories (Below %.0f m altitude)', altThreshold));
        subtitle("Positive North is Flight Origin");
        view([51 41]);
        legend(legendEntries, 'Location', 'bestoutside');
        
        % Plot 2D final landing positions
        figure('Name', 'Divert Landing Positions - 2D'); 
        hold on; grid on;
        
        for idx = 1:size(problemParams.divertPoints, 1)
            traj = divertData{idx}.stateTraj(:,1:3) * refVals.L_ref; % Convert to dimensional (MCMF)
            % Convert final position from MCMF to ENU
            finalPosMCMF = traj(end, 1:3)';
            finalPosENU = MCMF2ENU(finalPosMCMF, landingLatDeg, landingLonDeg, true, rPlanet);
            if idx==1
                plot(finalPosENU(1), finalPosENU(2), 'o', 'MarkerSize', 10, ...
                'MarkerFaceColor', colors(idx,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            else
                plot(finalPosENU(1), finalPosENU(2), 'o', 'MarkerSize', 10, ...
                'MarkerFaceColor', colors(idx,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'HandleVisibility','off');
            end

        end
        
        % Plot target divert points for reference
        for idx = 1:size(problemParams.divertPoints, 1)
            targetPoint = problemParams.divertPoints(idx, 1:2);
            plot(targetPoint(1), targetPoint(2), 'x', 'MarkerSize', 15, ...
                'Color', colors(idx,:), 'LineWidth', 3,'DisplayName',num2str(idx));
        end
        
        xlabel('East (m)'); ylabel('North (m)');
        title('Final Landing Positions vs Target Divert Points');
        legend('Actual Landing', 'Target Points', 'Location', 'best');
        axis equal; grid on;
        
        % Print divert info
        baseCaseMass = divertData{1}.stateTraj(end,7) * refVals.M_ref;
        baseCaseTOF = divertData{1}.tTraj(end) * refVals.T_ref;
        fprintf('\n=== Divert Landing Errors ===\n');
        for idx = 1:size(problemParams.divertPoints, 1)
            traj = divertData{idx}.stateTraj(:,1:3) * refVals.L_ref; % MCMF coordinates
            % Convert from MCMF to ENU for error calculation
            finalPosMCMF = traj(end, 1:3)';
            finalPosENU = MCMF2ENU(finalPosMCMF, landingLatDeg, landingLonDeg, true, rPlanet);
            targetPos = problemParams.divertPoints(idx, 1:2)';
            error = norm(finalPosENU(1:2) - targetPos);
            massDelta =  baseCaseMass - divertData{idx}.stateTraj(end,7)*refVals.M_ref;
            TOFDelta = divertData{idx}.tTraj(end)*refVals.T_ref - baseCaseTOF;
            fprintf('Point %d (E=%.0fm, N=%.0fm): Error = %.4f m, Increased Fuel Cost: %.1f kg, Increased TOF: %.1f s\n', ...
                idx, targetPos(1), targetPos(2), error, massDelta, TOFDelta);
        end
        fprintf('=============================\n');

        % Plot Throttle Profiles
        figure('Name', 'Divert Throttle Profiles');
        hold on; grid on;
        title(sprintf('Throttle Profiles (Post-Divert, <%.0fm)', altDivert));
        ylim([ 0 105]);
        xlabel("Time (s)"); ylabel("Throttle (%)");

        for idx = 1:size(problemParams.divertPoints, 1)
            tTraj = divertData{idx}.tTraj;
            mass = divertData{idx}.stateTraj(:,7);
            aT = divertData{idx}.aTSim;
            r = divertData{idx}.stateTraj(:,1:3);
            tDim = tTraj * refVals.T_ref;

            rMagND = vecnorm(r, 2, 2);
            altDim = (rMagND - rPlanetND) * refVals.L_ref;
            mask = altDim <= problemParams.altDivert;

            if any(mask)
                aTMag = vecnorm(aT, 2, 1)';
                throttle = (mass(mask) .* aTMag(mask)) ./ (nonDimParams.maxThrustND) * 100;

                if idx == 1
                    plot(tDim(mask), throttle, 'k-', 'LineWidth', 2, 'DisplayName', 'Original Trajectory');
                else
                    plot(tDim(mask), throttle, 'LineWidth', 1.5, 'DisplayName', sprintf("Pt %d", idx), 'Color', colors(idx, :));
                end
            end
        end
        legend('Location','bestoutside');
        

        % Plot Complete Thrust Profiles
        figure('Name', 'Full Flight Thrust Accel Profiles');
        hold on; grid on;
        xlabel('Time (s)'); ylabel('Thrust Accel (m/s^2)');
        title('Thrust Acceleration Profiles, Full-Flight');

        for idx = 1:size(problemParams.divertPoints, 1)
            tTraj = divertData{idx}.tTraj;
            aT = divertData{idx}.aTSim;
            tDim = tTraj * refVals.T_ref;
            aTMag = vecnorm(aT, 2, 1)';
            if idx == 1
                plot(tDim, aTMag * refVals.A_ref, 'k-', 'LineWidth', 2, 'DisplayName', 'Original Trajectory');
            else
                plot(tDim, aTMag * refVals.A_ref, 'LineWidth', 1.5, 'DisplayName', sprintf("Pt %d", idx), 'Color', colors(idx,:));
            end
        end
        legend('Location','bestoutside');

        % Plot Gamma1 vs Gamma 2 History
        figure('Name','Gamma 1 vs Gamma 2 History');
        hold on; grid on;
        xlabel('Gamma 1'); ylabel('Gamma 2');
        title('Parameter Evolution During Divert (Gamma 1 vs Gamma 2)');
    
        for idx = 1:size(problemParams.divertPoints, 1)
            currHist = divertData{idx}.optHistory;
            
            if ~isempty(currHist)
                g1 = currHist(:,2);
                g2 = currHist(:,3);
    
                % Store 'idx' in UserData so the tooltip can access it later
                if idx == 1
                    plot(g1, g2, 'ko', 'LineWidth', 1.5, 'MarkerSize', 4, ...
                         'MarkerFaceColor', 'k', 'DisplayName', 'Original Trajectory', ...
                         'UserData', idx);
                else
                    plot(g1, g2, 'o', 'LineWidth', 1.5, 'MarkerSize', 4, ...
                         'Color', colors(idx,:), 'MarkerFaceColor', colors(idx,:), ...
                         'DisplayName', sprintf('Pt %d', idx), ...
                         'UserData', idx);
                end
            end
        end
        legend('Location','bestoutside');
    
        % Enable Data Cursor Mode with Custom Tooltip
        dcm = datacursormode(gcf);
        dcm.Enable = 'on';
        dcm.UpdateFcn = @(obj, event_obj) { ...
            ['Gamma 1: ', num2str(event_obj.Position(1), '%.4f')], ...
            ['Gamma 2: ', num2str(event_obj.Position(2), '%.4f')], ...
            ['Pt: ', num2str(get(event_obj.Target, 'UserData'))] ...
        };
end