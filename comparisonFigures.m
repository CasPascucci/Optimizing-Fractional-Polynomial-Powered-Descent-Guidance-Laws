folder = 'Thesis Figures/OptVsSim/';

pairs = {
    'Optimization 3D.png',                'Simulation 3D.png',                  'Combined 3D.png';
};

for i = 1:size(pairs, 1)
    left  = imread(fullfile(folder, pairs{i,1}));
    right = imread(fullfile(folder, pairs{i,2}));

    % Match heights if they differ
    hL = size(left,1);  hR = size(right,1);
    if hL ~= hR
        h = max(hL, hR);
        left  = imresize(left,  [h, NaN]);
        right = imresize(right, [h, NaN]);
    end

    combined = [left, right];
    imwrite(combined, fullfile(folder, pairs{i,3}));
    fprintf('Saved %s\n', pairs{i,3});
end