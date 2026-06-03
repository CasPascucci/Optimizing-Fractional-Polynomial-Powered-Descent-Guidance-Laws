function plotting(tTraj, stateTraj, optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, aTList, refVals, problemParams, nonDimParams, optimParams, flag_thrustGotLimited, optHistory, ICstates, betaParam)


    %% Preprocess
    runTime = datestr(now, 'HH:MM:SS');
    
    descString = string.empty;
    
    descString(end+1) = sprintf("Beta: %.2f", betaParam);
    descString(end+1) = sprintf("Throttle Limit");
    if optimParams.updateOpt
        descString(end+1) = "ReOpt";

    end
    if optimParams.pointingEnabled
        descString(end+1) = "Pointing";
    end
    if optimParams.glideSlopeEnabled
        descString(end+1) = "GS";
    end
    if flag_thrustGotLimited
        descString(end+1) = "Saturated";
    end
    descString = join(descString, "+");
    legendTag = descString;% + " @ " + runTime;


    % Define Color for each run
    fig1Handle = findobj(0, 'Name', 'Optim Throttle');
    existingRuns = 0;
    if ~isempty(fig1Handle)
        % Count data entries on Fig 1
        fig1Entries = findobj(fig1Handle, 'Type', 'line');
        for i = 1:length(fig1Entries)
            if contains(fig1Entries(i).DisplayName, 'B')
                existingRuns = existingRuns + 1;
            end
        end
    end
    availableColors = colororder;
    colorIdx = mod(existingRuns, size(availableColors,1)) + 1;
    runColor = availableColors(colorIdx, :);
    % Check for simulation data availability
    hasSimData = ~isempty(tTraj) && ~isempty(stateTraj) && ~isempty(aTList);
    
    % Extract Parameters
    tgo0 = optParams(3);
    L_ref = refVals.L_ref;  
    T_ref = refVals.T_ref;  
    A_ref = refVals.A_ref;  
    V_ref = refVals.V_ref;  
    M_ref = refVals.M_ref;  
    rMoon = problemParams.rMoon; 
    [~, ~, U0] = enuBasis(deg2rad(problemParams.landingLatDeg),deg2rad(problemParams.landingLonDeg));
    isp = nonDimParams.ispND;
    

    % Process Simulation Data
    if hasSimData
        rState = stateTraj(:,1:3);
        rDim = rState * L_ref;
        vState = stateTraj(:,4:6);
        vDim = vState * V_ref;
        mState = stateTraj(:,7);
        mDim = mState * M_ref;
        %tgoState = tgo0 - tTraj;
        %tgoDim = tgoState * T_ref;

        aTDim = aTList * A_ref;
        aTDimNorm = vecnorm(aTDim,2,1);

        deltaR = MCMF2ENU(rDim',problemParams.landingLatDeg,problemParams.landingLonDeg,true,rMoon);
        East  = deltaR(1,:)';
        North = deltaR(2,:)';
        Up    = deltaR(3,:)';

        alt_m = vecnorm(rDim,2,2) - rMoon;

        vDimENU = MCMF2ENU(vDim',problemParams.landingLatDeg,problemParams.landingLonDeg,false,rMoon);
        vE = vDimENU(1,:)';
        vN = vDimENU(2,:)';
        vU = vDimENU(3,:)';

        aTDimENU = MCMF2ENU(aTDim,problemParams.landingLatDeg,problemParams.landingLonDeg,false,rMoon);
        aE = aTDimENU(1,:)';
        aN = aTDimENU(2,:)';
        aU = aTDimENU(3,:)';
        aT_ENU = [aE, aN, aU];
        aT_norm_ENU = vecnorm(aT_ENU,2,2);
    end

    maxThrustDim = problemParams.maxThrustDim;
    minThrustDim = problemParams.minThrustDim;

    %% Optimization Figures
    nodeCount = optimParams.nodeCount;
    aTOptim = aTOptim';
    rdOptim = rdOptim';
    vdOptim = vdOptim';
    mOptim = mOptim';
    
    rdOptimTOPO = MCMF2ENU(rdOptim',problemParams.landingLatDeg,problemParams.landingLonDeg,true,rMoon/L_ref);
    rdOptimTOPODim = rdOptimTOPO*L_ref;

    tgospanOpt = linspace(0,tgo0,nodeCount);
    tspanOpt = tgo0 - tgospanOpt;
    aTNormOpt = vecnorm(aTOptim,2,2);

    EOpt = rdOptimTOPODim(1, :);
    NOpt = rdOptimTOPODim(2, :);
    UOpt = rdOptimTOPODim(3, :);
    
    % Figure 1: Optimization Throttle Profile
    figure(1); hold on; grid on;
    set(gcf, 'Name', 'Optim Throttle')
    thrustDim = aTNormOpt .* mOptim *(refVals.M_ref*refVals.A_ref);
    
    plot(tspanOpt(1:end)*T_ref, thrustDim(1:end)/problemParams.maxThrustDim, '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName', legendTag);
    
    if isempty(findobj(gca, 'DisplayName', 'Max Throttle'))
        yline(1.0, 'k--', 'LineWidth', 1, 'DisplayName', 'Max Throttle');
        yline(problemParams.minThrustDim/problemParams.maxThrustDim, 'k--', 'LineWidth', 1, 'DisplayName', 'Min Throttle');
    end

    xlabel('Time s'); ylabel('Throttle Fraction');
    title('Optimization Throttle Profile');
    ylim([0, 1.2]);
    legend('Location','best');
    set(gca, 'FontSize', 20);

    % Figure 2: Optimization Range vs Altitude
    rdOptimDim = rdOptim*L_ref;
    rhatOpt = rdOptimDim ./ vecnorm(rdOptimDim,2,2);
    centralOpt = acos(rhatOpt*U0);
    arcLengthOpt = rMoon * centralOpt;
    alt_opt = vecnorm(rdOptimDim,2,2) - rMoon;

    figure(2); hold on; grid on;
    set(gcf, 'Name', 'Optimization Ground Range vs Alt');

    plot(arcLengthOpt/1000, alt_opt/1000, '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName', legendTag);
    xlabel('Range km'); ylabel('Up km');
    title('Optimization Flight Path (Range vs Alt)');
    legend('Location','best');
    set(gca, 'FontSize', 20);
    
    % Figure 3: Optimization 3D Trajectory
    UPlotMax = 3000; % Alt filter threshold
    t_nodes = linspace(0, tgo0, nodeCount);
    t_highres = linspace(0, tgo0, 2000); 
    
    ETraj_HR = spline(t_nodes, EOpt, t_highres);
    NTraj_HR = spline(t_nodes, NOpt, t_highres);
    UTraj_HR = spline(t_nodes, UOpt, t_highres);
    
    idx_plot = UTraj_HR <= UPlotMax; 
    if ~any(idx_plot); idx_plot = true(size(UTraj_HR)); end 

    figure(3); hold on; grid on; axis equal;
    set(gcf, 'Name', 'Optimization 3D');
    set(gca, 'FontSize', 20);
    
    %Generate Glideslope Cone Visualization
    if isempty(findobj(gca, 'Type', 'Surface'))
        theta_fun = @(u) min(89.9, 45 + 45 * max(0, min(1, (u - 250) ./ 250))); 
        Umin = 0; 
        Umax = min(500, max(UTraj_HR(idx_plot), [], 'omitnan')); 
        if isempty(Umax) || isnan(Umax); Umax = 500; end
        U_samp = linspace(Umin, Umax, 360);
        theta_deg = theta_fun(U_samp);
        cosBound = cosd(theta_deg);
        cosBound = max(cosBound, 1e-3);
        R_samp = U_samp .* sqrt(1 - cosBound.^2) ./ cosBound; 
        Rcap = 3000;
        R_samp = min(R_samp, Rcap);
        azimuth = linspace(0, 2*pi, 360);
        [azGrid, altGrid] = meshgrid(azimuth, U_samp);
        [~, rGrid] = meshgrid(azimuth, R_samp);
        eastGrid  = rGrid .* cos(azGrid);
        northGrid = rGrid .* sin(azGrid);

        cone1 = surf(eastGrid/1000, northGrid/1000, altGrid/1000,altGrid/1000,'EdgeAlpha',0.15,'FaceAlpha',0.4,'MeshStyle','row','LineWidth',0.8, 'HandleVisibility','off');
        % Plateau Ring
        UPlat = Umax; RPlat = Rcap;
        azPlat = linspace(0, 2*pi, 360);
        eastPlat = RPlat * cos(azPlat);
        northPlat = RPlat * sin(azPlat);
        UPlat = UPlat * ones(size(azPlat));
        cone2 = plot3(eastPlat/1000, northPlat/1000, UPlat/1000, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    end
    
    % Main Trajectory
    plot3(ETraj_HR(idx_plot)/1000, NTraj_HR(idx_plot)/1000, UTraj_HR(idx_plot)/1000, '.-', 'Color', runColor, 'LineWidth', 2, 'MarkerSize',5, 'DisplayName', legendTag);
    
   
    xlabel('East (km)'); ylabel('North (km)'); zlabel('Up (km)');
    title('Optimization 3D Trajectory (Final 2km)');
    subtitle('With Glideslope Cone');
    view(80, 15); camproj orthographic;
    zlim([0,2]); xlim([-3,3]); ylim([-3,3]);
    axis square;
    legend('Location','north');
    set(gca, 'FontSize', 20);

    % Figure 4: Optimization Velocity Profile
    vdOptimDim = vdOptim * V_ref;
    vdOptimTOPO = MCMF2ENU(vdOptimDim',problemParams.landingLatDeg,problemParams.landingLonDeg,false,rMoon);
    figure(4); hold on; grid on;
    set(gcf, 'Name', 'Opt Vel');
    set(gca, 'FontSize', 20);
    
    plot(tspanOpt*T_ref, vecnorm(vdOptimTOPO, 2, 1), '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName', legendTag);
    % plot(tspanOpt*T_ref, vdOptimTOPO(1,:), 'b-', 'LineWidth', 1.0, 'DisplayName', [legendTag ' East']);
    % plot(tspanOpt*T_ref, vdOptimTOPO(2,:), 'b--', 'LineWidth', 1.0, 'DisplayName', [legendTag ' North']);
    % plot(tspanOpt*T_ref, vdOptimTOPO(3,:), 'b:', 'LineWidth', 1.0, 'DisplayName', [legendTag ' Up']);
    

    legend('Location', 'best');
    xlabel('Time s'); ylabel('Velocity m/s'); title('Optimization Velocity Profile');
    grid on;
    xlim([0, ceil(max(tspanOpt*T_ref)/100)*100]);

    % Figure 5: Optimization TOPO Acceleration
    aTOptimDim = aTOptim * A_ref;
    aTOptimTOPO = MCMF2ENU(aTOptimDim',problemParams.landingLatDeg,problemParams.landingLonDeg,false,rMoon);
    figure(5); hold on; grid on;
    set(gcf, 'Name', 'Opt Accel');
    set(gca, 'FontSize', 20);
    
    plot(tspanOpt(1:end)*T_ref, vecnorm(aTOptimTOPO(:,1:end),2,1), '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName',legendTag);
    % plot(tspanOpt(2:end)*T_ref, aTOptimTOPO(1,2:end), 'b-', 'LineWidth', 1.0, 'DisplayName', [labelMainOpt ' East']);
    % plot(tspanOpt(2:end)*T_ref, aTOptimTOPO(2,2:end), 'b--', 'LineWidth', 1.0, 'DisplayName', [labelMainOpt ' North']);
    % plot(tspanOpt(2:end)*T_ref, aTOptimTOPO(3,2:end), 'b:', 'LineWidth', 1.0, 'DisplayName', [labelMainOpt ' Up']);

    if isempty(findobj(gca, 'DisplayName', 'Desired Final Acceleration, a_f^*'))
        afStarDim = norm(nonDimParams.afStarND) * A_ref;
        yline(afStarDim, 'k--', 'LineWidth', 1, 'DisplayName', 'Desired Final Acceleration, a_f^*');
    end
    legend('Location', 'best');
    xlabel('Time s'); ylabel('Acceleration m/s^2'); title('Optimization Acceleration Profile');
    grid on;

    %% Simulation Figures
    if hasSimData
        % Figure 6: 3D Trajectory
        idx2km = find(alt_m <= 2000, 1, 'first');
        if isempty(idx2km)
            NtotalSim = numel(alt_m);
            idx2km = max(1, NtotalSim - floor(0.20*NtotalSim) + 1);
        end
        idx_sim = idx2km:numel(alt_m);

        ETrajSim = East(idx_sim);
        NTrajSim = North(idx_sim);
        UTrajSim = Up(idx_sim);

        figure(6); hold on; grid on; axis equal;
        set(gcf, 'Name', 'Simulation 3D');
        set(gca, 'FontSize', 20);

        if isempty(findobj(gca, 'Type', 'Surface'))
            theta_fun = @(u) min(89.9, 45 + 45 * max(0, min(1, (u - 250) ./ 250))); 
            Umin = 0; 
            Umax = min(500, max(UTraj_HR(idx_plot), [], 'omitnan')); 
            if isempty(Umax) || isnan(Umax); Umax = 500; end
            U_samp = linspace(Umin, Umax, 360);
            theta_deg = theta_fun(U_samp);
            cosBound = cosd(theta_deg);
            cosBound = max(cosBound, 1e-3);
            R_samp = U_samp .* sqrt(1 - cosBound.^2) ./ cosBound; 
            Rcap = 3000;
            R_samp = min(R_samp, Rcap);
            azimuth = linspace(0, 2*pi, 360);
            [azGrid, altGrid] = meshgrid(azimuth, U_samp);
            [~, rGrid] = meshgrid(azimuth, R_samp);
            eastGrid  = rGrid .* cos(azGrid);
            northGrid = rGrid .* sin(azGrid);
            
            surf(eastGrid/1000, northGrid/1000, altGrid/1000,altGrid/1000,'EdgeAlpha',0.15,'FaceAlpha',0.4,'MeshStyle','row','LineWidth',0.8, 'HandleVisibility','off');
            % Plateau Ring
            UPlat = Umax; RPlat = Rcap;
            azPlat = linspace(0, 2*pi, 360);
            eastPlat = RPlat * cos(azPlat);
            northPlat = RPlat * sin(azPlat);
            UPlat = UPlat * ones(size(azPlat));
            plot3(eastPlat/1000, northPlat/1000, UPlat/1000, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
        end

        plot3(ETrajSim/1000, NTrajSim/1000, UTrajSim/1000, '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName', legendTag);


        xlabel('East (km)'); ylabel('North (km)'); zlabel('Up (km)');
        title('Simulation 3D Trajectory (Final 2km)');
        view(80, 15); camproj orthographic;
        zlim([0, 2]); xlim([-3,3]); ylim([-3,3]); axis square;
        legend('Location','bestoutside');

        % Figure 7: Range vs Altitude
        rhat      = rDim ./ vecnorm(rDim,2,2);
        central   = acos(rhat*U0);
        arcLength = rMoon * central;
        
        figure(7); hold on; grid on;
        set(gcf, 'Name', 'Simulation Ground Range vs Alt');
        set(gca, 'FontSize', 20);
        plot(arcLength/1000, alt_m/1000, '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName', legendTag);
        if norm([alt_m(end), arcLength(end)]) > 1
            fprintf("Landing site error > 1m: %.2fm\n", norm([alt_m(end), arcLength(end)]));
        end

        xlabel('Range (km)'); ylabel('Up (km)');
        title('Simulation Flight Path (Range vs Alt)'); grid on;
        legend('Location','best');

        % Figure 8: TOPO Acceleration
        figure(8); hold on; grid on;
        set(gcf, 'Name', 'Sim Accel');
        set(gca, 'FontSize', 20);

        plot(tTraj*T_ref, aT_norm_ENU, '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName', legendTag);
        % plot(tTraj*T_ref, aE, 'b-', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' East']);
        % plot(tTraj*T_ref, aN, 'b--', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' North']);
        % plot(tTraj*T_ref, aU, 'b:', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' Up']);
        
        if isempty(findobj(gca, 'DisplayName', 'Desired Final Acceleration, a_f^*'))
            afStarDim = norm(nonDimParams.afStarND) * A_ref;
            yline(afStarDim, 'k--', 'LineWidth', 1, 'DisplayName', 'Desired Final Acceleration, a_f^*');
        end
        legend('Location', 'best');
        xlabel('Time s'); ylabel('Accel m/s^2'); title('Simulation Acceleration Profile');
        if flag_thrustGotLimited
            subtitle("Thrust is being Throttled");
        end

        % Figure 9: Velocity Profile
        figure(9); hold on; grid on;
        set(gcf, 'Name', 'Sim Velocity');
        set(gca, 'FontSize', 20);

        plot(tTraj*T_ref, vecnorm(vDim, 2, 2), '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName', legendTag);
        % plot(tTraj*T_ref, vE, 'b-', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' East']);
        % plot(tTraj*T_ref, vN, 'b--', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' North']);
        % plot(tTraj*T_ref, vU, 'b:', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' Up']);

        legend('Location', 'best');
        xlabel('Time s'); ylabel('Velocity m/s'); title('Simulation Velocity Profile');
        grid on;
        
        % Figure 10: Throttle
        figure(10);
        set(gcf, 'Name', 'Sim Throttle');
        set(gca, 'FontSize', 20);
        hold on; grid on;
        thrustDim = aTDimNorm' .* mDim;
        plot(tTraj*T_ref, thrustDim/maxThrustDim, '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName', legendTag);
        if isempty(findobj(gca, 'DisplayName', 'Max Throttle'))
            yline(1.0, 'k-', 'LineWidth', 1, 'DisplayName', 'Max Throttle');
            yline(minThrustDim/maxThrustDim, 'k-', 'LineWidth', 1, 'DisplayName', 'Min Throttle');
        end
        xlabel('Time s'); ylabel('Throttle Fraction'); title('Simulation Throttle Profile');
        ylim([0, 1.2]);
        legend('Location','best');
        
        % Figure 11: Mass
        figure(11);
        set(gcf, 'Name', 'Sim Mass');
        set(gca, 'FontSize', 20);
        hold on; grid on;
        plot(tTraj*T_ref, mDim, '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName', legendTag);
        xlabel('Time s'); ylabel('Mass kg'); title('Vehicle Mass vs Time');
        legend('Location','best');
        propUsedC = max(0, mDim(1) - mDim(end));
        subtitle(sprintf('Propellant used = %.1f kg', propUsedC));

        % Figure 12: Time vs Altitude
        figure(12); hold on; grid on;
        set(gcf, 'Name', 'Sim Alt');
        set(gca, 'FontSize', 20);
        plot(tTraj*T_ref, alt_m/1000, '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName', legendTag);

        xlabel('Time s'); ylabel('Altitude m'); title('Time vs Altitude (Dimensional)');
        legend('Location','best');
        grid on;

        % Figure 13: Pointing Angle Analysis
    
        % Use optimization trajectory instead of sim
        aE_opt = aTOptimTOPO(1,:)';
        aN_opt = aTOptimTOPO(2,:)';
        aU_opt = aTOptimTOPO(3,:)';
        aMag_opt = vecnorm([aE_opt aN_opt aU_opt], 2, 2);
        thrustU_opt = [aE_opt aN_opt aU_opt] ./ aMag_opt;
        dotUp_opt = max(-1, min(1, thrustU_opt(:,3)));
        phiOpt = acosd(dotUp_opt);


        figure(13);
        set(gcf, 'Name', 'Pointing Analysis');
        set(gca, 'FontSize', 20);
        % subplot(2,1,1); hold on; grid on;
        % set(gca, 'FontSize', 20);
        plot(tspanOpt*T_ref, phiOpt, '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName', legendTag);
        hold on; grid on;
        if optimParams.pointingEnabled
            phi0_deg   = optimParams.minPointing;
            tgoOpt = tgospanOpt' * T_ref;
            phiA_opt = 0.5 * optimParams.maxTiltAccel .* (tgoOpt.^2);
            ThetaOpt = min(180, phi0_deg + phiA_opt);
            plot(tspanOpt*T_ref, ThetaOpt, '--', 'Color', runColor, 'LineWidth', 2, 'DisplayName','\Theta Limit');

            xlabel('Time s'); ylabel('Degrees');
            title('Pointing Angle vs Limit');
            legend('Location','best');

        else
            xlabel('Time s'); ylabel('Degrees');
            title('Pointing Angle vs Limit');
            legend('Location','bestoutside');
        end
        % Figure 14: Pointing Angle Analysis (Simulation)
        aMag_sim = vecnorm(aT_ENU, 2, 2);
        thrustU_sim = aT_ENU ./ aMag_sim;
        dotUp_sim = max(-1, min(1, thrustU_sim(:,3)));
        phiSim = acosd(dotUp_sim);

        figure(14);
        set(gcf, 'Name', 'Pointing Analysis (Simulation)');
        set(gca, 'FontSize', 20);
        plot(tTraj*T_ref, phiSim, '-', 'Color', runColor, 'LineWidth', 2, 'DisplayName', legendTag);
        hold on; grid on;
        if optimParams.pointingEnabled
            phi0_deg  = optimParams.minPointing;
            tgoSim    = (tgo0 - tTraj) * T_ref;
            phiA_sim  = 0.5 * optimParams.maxTiltAccel .* (tgoSim.^2);
            ThetaSim  = min(180, phi0_deg + phiA_sim);
            plot(tTraj*T_ref, ThetaSim, '--', 'Color', runColor, 'LineWidth', 2, 'DisplayName', '\Theta Limit');
        end
        xlabel('Time s'); ylabel('Degrees');
        title('Pointing Angle vs Limit (Simulation)');
        legend('Location','best');

    else
        fprintf('\n=== NOTE: Simulation plots skipped (runSimulation = false) ===\n');
    end

    %% 4. Parameter History (Re-Optimization Only)
    if ~isempty(optHistory) && exist('ICstates', 'var') && ~isempty(ICstates)
        figure(15);
        set(gcf, 'Name', 'ReOpt Parameter History');
        
        t_elapsedND = table2array(optHistory(:,1));
        x_opt = t_elapsedND .* refVals.T_ref;
        gamma_hist = table2array(optHistory(:,2));
        gamma2_hist = table2array(optHistory(:,3));
        tgoSolved_ND = table2array(optHistory(:,5));
        tgoSolved = tgoSolved_ND .* refVals.T_ref;
        
        ICstates_array = table2array(ICstates);
        nPoints = size(optHistory, 1);
        c1_norm = zeros(nPoints, 1);
        c2_norm = zeros(nPoints, 1);
        
        rfStar = nonDimParams.rfStarND;
        vfStar = nonDimParams.vfStarND;
        afStar = nonDimParams.afStarND;
        gConst = nonDimParams.gConst;
        
        for i = 1:nPoints
            r_i = ICstates_array(i, 1:3)';
            v_i = ICstates_array(i, 4:6)';
            [c1, c2] = calculateCoeffs(r_i, v_i, tgoSolved_ND(i), gamma_hist(i), gamma2_hist(i), afStar, rfStar, vfStar, gConst);
            c1_norm(i) = norm(c1);
            c2_norm(i) = norm(c2);
        end
        
        % Continue curve through static time at end of flight
        if hasSimData
            tFinal = tTraj(end) * refVals.T_ref;
            x_opt = [x_opt; tFinal];
            gamma_hist = [gamma_hist; gamma_hist(end)];
            gamma2_hist = [gamma2_hist; gamma2_hist(end)];
            tgoSolved = [tgoSolved; 0]; % Tgo goes to 0
            c1_norm = [c1_norm; c1_norm(end)];
            c2_norm = [c2_norm; c2_norm(end)];
        end

        tgoInitial = tgoSolved(1);
        tgoStaticAtUpdate = tgoInitial - x_opt;
        tgoDiff = tgoStaticAtUpdate - tgoSolved;
        
        subplot(2,1,1); hold on; grid on; plot(x_opt, gamma_hist, 'LineWidth', 2, 'Color', runColor, 'DisplayName', legendTag); title('\gamma_1'); grid on; xlabel('Time s'); legend('Location','northwest');
        subplot(2,1,2); hold on; grid on; plot(x_opt, gamma2_hist, 'LineWidth', 2, 'Color', runColor, 'DisplayName', legendTag); title('\gamma_2'); grid on; xlabel('Time s');
        % subplot(3,2,3); hold on; grid on; plot(x_opt, tgoSolved, 'LineWidth', 2, 'Color', runColor, 'DisplayName', legendTag); title('t_{go} Solved'); grid on; xlabel('Time s');
        % subplot(3,2,4); hold on; grid on; plot(x_opt, c1_norm, 'LineWidth', 2, 'Color', runColor, 'DisplayName', legendTag); title('||c_1||'); grid on; xlabel('Time s');
        % subplot(3,2,5); hold on; grid on; plot(x_opt, c2_norm, 'LineWidth', 2, 'Color', runColor, 'DisplayName', legendTag); title('||c_2||'); grid on; xlabel('Time s');
        % subplot(3,2,6); hold on; grid on; stairs(x_opt, tgoDiff, 'LineWidth', 2, 'Color', runColor, 'DisplayName', legendTag); title('\Delta t_{go} (Static - ReOpt)'); ylabel('s'); xlabel('Time s'); grid on;
        
        sgtitle('Optimization Parameters Through Reoptimization');
        set(gca, 'FontSize', 20);
    end

end