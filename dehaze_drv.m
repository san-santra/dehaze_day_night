%
% Copyright (C) 2019  Sanchayan Santra
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% You should have received a copy of the GNU Lesser General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% script to call the dehazing function

path_prefix = 'images';
out_folder = 'output';

name = 'Oberfallenberg';
files = dir(fullfile(path_prefix, ['*', name, '*']));
files = {files.name}';

im = imread(fullfile(path_prefix, files{1}));

% display the debuggin figures
display_fig = false;
write_im = true;

% default parameters (all possible)
window_size = 8;
ransac_thresh = 0.02;

inlier_frac = 0.4;
modality_thresh = 0.06;
origin_test_thresh = 0.0005;

dc_thr_frac = 0.1;
t_step = 3;
hough_thr_frac = 0.3;

err_sigma = 1.0/30;
intersect_angle_thresh = 15;
intersect_thresh = 0.05;
shading_variability_thresh = 0.006;

min_i_diff_sq = 0.0001;
lambda = 1;
numsample = 5;
alpha = 1/5000;
beta = 0.000001;

gamma = 0.8;

[enh_im, air_comp, air_est, rem, inlier_fit_step] ...
    = dehaze_general(im, display_fig, window_size, ransac_thresh, ...
    inlier_frac, modality_thresh, origin_test_thresh, dc_thr_frac, ...
    t_step, hough_thr_frac, err_sigma, intersect_angle_thresh, ...
    intersect_thresh, shading_variability_thresh, min_i_diff_sq, ...
    lambda, numsample, alpha, beta);

if write_im
    imwrite(enh_im.^gamma, fullfile(out_folder, [name, '_fixed.png']));
    imwrite(air_comp, fullfile(out_folder, [name, '_airlight.png']));
    imwrite(air_est, fullfile(out_folder, [name, '_airlight_est.png']));
    imwrite(rem, fullfile(out_folder, [name, '_rem.png']));
    imwrite(inlier_fit_step, fullfile(out_folder, [name, '_inlier_fit_step.png']));
end
