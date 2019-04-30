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

function [air_est_col, est_conf_col] = compute_line_shift(im, patch_size, ...
    xrow, ycol, color_line_D, color_line_P, inlier_mask, A_dir_mat, ...
    labels, err_sigma, intersect_angle_thresh, intersect_thresh, ...
    shading_variability_tresh)

    % A_dir_mat ~> 3 x num_A_dir

    numpatch = size(color_line_D, 2);
    numlight = size(A_dir_mat, 2);
    air_est_col = zeros(patch_size^2, numpatch);
    est_conf_col = zeros(patch_size^2, numpatch);
    
    for ii=1:numpatch
        D = color_line_D(:, ii);
        P = color_line_P(:, ii);
        
        if(all(P) == 0)
            continue;
        end
        
        if(numlight ~= 1)
            if(labels(ii) == 0)
                continue;
            end
            A_dir = A_dir_mat(:, labels(ii));
        else
            A_dir = A_dir_mat(:, 1);
        end
        
        % compute intersection (shift)
        a_dot_d = A_dir'*D;
        d_dot_p = D'*P;
        a_dot_p = A_dir'*P;
        
        r = (1/(1-a_dot_d^2))*([1, a_dot_d; a_dot_d, 1]*[-d_dot_p; a_dot_p]);
        
        l = r(1);
        s = r(2);
        
        obj_fn_vec = (l*D + P - s*A_dir);
        obj_fn_val = obj_fn_vec'*obj_fn_vec;
        
        % confidence of the estimate
        num = A_dir - D*a_dot_d;
        conf = 0.5*err_sigma^2*(num'*num)/(1 - (a_dot_d*a_dot_d))^2;
        
        % some more tests to validate
        % intersection angle > intersect_angle_thresh
        cos_angle = D'*A_dir;  % both are unit vectors
        intersect_angle_test = (cos_angle <= cosd(intersect_angle_thresh));
        
        % range_test
        row = xrow(ii);
        col = ycol(ii);
        patch = im(row:row+patch_size -1, col:col+patch_size -1, :);
        patch_data = reshape(permute(patch, [3 1 2]), [3, patch_size^2]);
        inlier_data = patch_data(:, inlier_mask(:, ii));
        
        patch_i_norm = sqrt(sum(inlier_data.*inlier_data, 1));
        min_i_norm = min(patch_i_norm);
        range_test = s > 0 && (s < 1.05*min_i_norm);
        
        % intersection error test
        intersect_test = obj_fn_val < intersect_thresh;
        
        test_val = intersect_angle_test && intersect_test && range_test; ...
        
        if(test_val)
            air_est_col(:, ii) = inlier_mask(:, ii).*s;
            est_conf_col(:, ii) = conf.*inlier_mask(:, ii);
            
            % should the label be discarded also ?
            % labels(ii) = 0;
        end
    end
end
