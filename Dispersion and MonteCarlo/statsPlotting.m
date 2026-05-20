function statsPlotting(Results, beta, optimizationParams)
exitFlags = [Results.exitflag];
validFlags = exitFlags == 1 | exitFlags == 2;
Results = Results(validFlags);

ok   = [Results.exit_ok]' == 1;
gamma= [Results.gamma]';
gamma2=[Results.gamma2]';
kr   = [Results.kr]';
tgo  = [Results.tgo]';        % seconds
fuel_opt = [Results.fuel_opt]';
fuel_sim = [Results.fuel_sim]';
coeff1   = [Results.coeff1]';
coeff2   = [Results.coeff2]';
coeff3   = [Results.coeff3]';
coeff4   = [Results.coeff4]';

nOK = sum(ok); N = numel(Results);
fprintf('Plotted: %d / %d  (%.1f%%), converged cases\n', nOK, N, 100*nOK/max(N,1));

vars = {'gamma','gamma2','tgo','fuel_opt','fuel_sim','coeff1','coeff2','coeff3','coeff4'};
Xall = [gamma,gamma2,tgo,fuel_opt,fuel_sim,coeff1,coeff2,coeff3,coeff4];

condTags = {};
    if isfield(optimizationParams,'glideSlopeEnabled') && optimizationParams.glideSlopeEnabled
        condTags{end+1} = 'GS';
    end
    if isfield(optimizationParams,'pointingEnabled') && optimizationParams.pointingEnabled
        condTags{end+1} = 'PT';
    end
    if isempty(condTags)
        condTags = {'Thrust Only'};
    end
    
    % Create a unique timestamp for figure names so they don't overwrite
    runID = datestr(now, 'HH:MM:SS');
    mainTitle = sprintf('Dispersion Study (Beta: %.2f | Constraints: %s)', beta, strjoin(condTags, ', '));
%% Histograms
Xs = Xall(ok,:);
figure('Name',[strjoin(condTags, ', '), ' Histograms (success only) - ' runID]);
t = tiledlayout(3,3,'Padding','compact','TileSpacing','compact');
for i = 1:numel(vars)
    nexttile; histogram(Xs(:,i), 50); grid on;
    xlabel(vars{i},"Interpreter","none"); ylabel('count'); title(['Hist ' vars{i}],'Interpreter','none');
end
sgtitle(t, {mainTitle; 'Parameter Histograms (Success Only)'});
%% Boxcharts
figure('Name',[strjoin(condTags, ', '), ' Boxcharts by exit_ok - ' runID]);
t = tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

for i = 1:numel(vars)
    v = [Results.(vars{i})]';
    ax = nexttile;
    boxchart(ax, v); grid(ax,'on');
    title(ax, ['Box ' vars{i}],"Interpreter","none"); ylabel(ax, vars{i});
    xticklabels("");
end
sgtitle(t, {mainTitle; 'Parameter Boxcharts by Exit Status'});
%% Final Position Dispersion
% Calculate and plot the final position dispersion
finalPositions = [Results.final_error]';
figure('Name',[strjoin(condTags, ', '), ' Final Position Dispersion - ', runID]);
h = plot3(finalPositions(:,1),finalPositions(:,2),finalPositions(:,3),'x');
hold on;
xlabel('East Position'); ylabel('North Position'); zlabel('Up Position');
title({mainTitle; 'Final Position Dispersion for Converged Cases'});

idx = [Results.k]';
h.UserData = idx;
dataTip = h.DataTipTemplate;
dataTip.DataTipRows(end+1) = dataTipTextRow('Index', idx);
grid on;
end