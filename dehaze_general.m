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

function [enh_im, air_comp, air_est, rem, inlier_fit_step] ...
    = dehaze_general(in_im, display, window_size, ransac_thresh, ...
    inlier_frac, modality_thresh, origin_test_thresh, dc_thr_frac, ...
    t_step, hough_thr_frac, err_sigma, intersect_angle_thresh, ...
    intersect_thresh, shading_variability_thresh, min_i_diff_sq, ...
    lambda, numsample, alpha, beta)

    im = im2double(in_im);

    % fit line and get the patches and related data
    [color_line_D, color_line_P, plane_normal, patch_dc, ...
        inlier_mask, xrow, ycol] = fit_line_patches(im, window_size, ...
        ransac_thresh, inlier_frac, modality_thresh, origin_test_thresh, ...
        shading_variability_thresh);

    % compute A_dir_mat
    % display starts at figure 1
    [patch_labels, A_dir_mat] = compute_A_dir_mat_hough( ...
        plane_normal, patch_dc, dc_thr_frac, t_step, hough_thr_frac, display);
    
    % calculate shift
    [air_est_col, est_conf_col] = compute_line_shift(im, window_size, ...
        xrow, ycol, color_line_D, color_line_P, inlier_mask, A_dir_mat, ...
        patch_labels, err_sigma, intersect_angle_thresh, intersect_thresh, ...
        shading_variability_thresh);
    
    % aggregate
    [nrow, ncol, ~] = size(im);
    [labels_mat, air_est, est_conf, inlier_fit_step] = ...
        aggregate_from_patches(patch_labels, nrow, ncol, window_size, ...
        xrow, ycol, air_est_col, est_conf_col, inlier_mask);
    
    if(display)
        figure(10);
        imagesc(labels_mat);
    end
    
    % interpolate
    % fill up \hat{A}
    numlabel = size(A_dir_mat, 2);
    [~, influence] = interp_A_laplacian(labels_mat, ...
        numlabel, im, lambda, min_i_diff_sq);
    normalized_infl = influence./repmat(sum(influence, 2), [1, numlabel]);

    % a_est interpolation
    a_interp = interp_a_est_laplacian(im, air_est, est_conf, ...
        numsample, alpha, beta, min_i_diff_sq);
    
    a_interp(a_interp < 0) = 0;
    
    im_A_col = zeros(3, nrow*ncol);
    % weighted sum of the inlfluences
    for ii=1:3
        res = bsxfun(@times, normalized_infl, A_dir_mat(ii, :));
        res = sum(res, 2);
        im_A_col(ii, :) = res';
    end
    
    im_A_col_r = reshape(im_A_col, 3, nrow, ncol);
    im_A_mat = permute(im_A_col_r, [2 3 1]);
    air_comp = repmat(a_interp, [1 1 3]).*im_A_mat;
    rem = im - 0.90*air_comp;
    
    rem(rem < 0) = 0;
    rem(rem > 1) = 1;
    
    % contrast enhancement
    l = rgb2gray(air_comp);
    enh_im = rem./repmat(1 - l, [1 1 3]);
 
end
