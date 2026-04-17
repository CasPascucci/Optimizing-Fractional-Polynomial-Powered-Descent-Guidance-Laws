folder = 'Thesis Figures/OptVsSim/GS';

pairs = {
    'Optim Throttle.png',                'Sim Throttle.png',                  'CombinedVert Throttle.png';
};

for i = 1:size(pairs, 1)
    top    = imread(fullfile(folder, pairs{i,1}));
    bottom = imread(fullfile(folder, pairs{i,2}));

    % Match widths if they differ (width is the 2nd dimension)
    wT = size(top, 2);  
    wB = size(bottom, 2);
    
    if wT ~= wB
        w = max(wT, wB);
        % imresize with [NaN, w] scales the height automatically to maintain aspect ratio
        top    = imresize(top,    [NaN, w]);
        bottom = imresize(bottom, [NaN, w]);
    end

    % Stack vertically using a semicolon
    combined = [top; bottom];
    
    imwrite(combined, fullfile(folder, pairs{i,3}));
    fprintf('Saved %s\n', pairs{i,3});
end