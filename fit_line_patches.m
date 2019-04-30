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


function [color_line_D, color_line_P, plane_normal, patch_dc, ...
    inlier_mask, xrow, ycol] = fit_line_patches(im, window_size, ...
    ransac_thresh, inlier_frac, modality_thresh, origin_test_thresh, ...
    shading_variability_thresh)

    [nrow, ncol, ch] = size(im);
    
    if(ch ~= 3)
        error('The method works for RGB images only');
    end
    
    addpath('RANSAC_code'); % this is same as the one in dehaze_varlight
    
    % break into patches
    patch_size = window_size;
    numoverlap = 2;
    
    overlap = floor(patch_size/numoverlap);
    x = 1:patch_size - overlap:nrow - patch_size+1;
    y = 1:patch_size - overlap:ncol - patch_size+1;
    [X, Y] = meshgrid(x, y);
    
    xrow = X(:);
    ycol = Y(:);
    
    plane_normal = zeros(3, length(xrow));
    color_line_D = zeros(3, length(xrow));
    color_line_P = zeros(3, length(xrow));
    patch_dc = zeros(1, length(xrow));
    inlier_mask = false(patch_size^2, length(xrow));
    
    % for all patches estimate line
    for ii=1:length(xrow)
        % Fit line considering the RGB values as points 
        % also compute the plane containing this line and the origin
        
        row = xrow(ii);
        col = ycol(ii);
        
        patch = im(row:row+patch_size -1, col:col+patch_size -1, :);
        patch_data_vec = reshape(permute(patch, [3 1 2]), [3, patch_size^2]);
        
        % RANSAC
        [~, L, inliers] = ransacfitline(patch_data_vec, ransac_thresh);
        
        % L is empty if RANSAC is unable to fit the model
        model_fit_test = ~isempty(L);
        
        if(model_fit_test)
            % line parameters and the normal
            if(L(:, 1) > L(:, 2))
                D = L(:, 1) - L(:, 2);
                P = L(:, 2);
            else
                D = L(:, 2) - L(:, 1);
                P = L(:, 1);
            end
            
            normal = mycrossp(L(:, 1), L(:, 2));
            % unit vector
            normal = normal./sqrt(normal'*normal);
            D = D./sqrt(D'*D);
            
            % validate the fitted line
            % # inlier > inlier_frac
            num_inlier_test = length(inliers) >= inlier_frac*size(patch_data_vec, 2);

            % the remaining test are done for inliers only
            inlier_data = patch_data_vec(:, inliers);
            
            % D positive ( all component )
            dir_pos_test = all(D >= 0);
            
            % unimodality
	    % by threshodling gradient
            % option 1 - edge() of patch, if it contains an edge pixel among
            % the inliners
            edge_patch = edge(rgb2gray(patch), 'sobel', modality_thresh);
            modality_test = ~any(edge_patch(inliers));
            
            
            % line through origin ?
            % changing to perpenticular distance from origin
            rho = - sum(D.*P, 1);
            d_sq = (rho.*D + P).^2;
            % one sum should be here
            
            not_through_origin_test = d_sq >= origin_test_thresh;
            % as all components of D and P are positive this will do
            
            % shading variablity
            I_m_P = bsxfun(@minus, inlier_data, P);
            mul_D = bsxfun(@times, I_m_P, D);
            ImP_dot_D = sum(mul_D, 1);
         
            variability = std(ImP_dot_D);
            variability_test = variability > shading_variability_thresh;
            
            test_val = num_inlier_test && dir_pos_test && modality_test ...
                && variability_test;
            
            patch_dc(1, ii) = min(patch(:));
            
            if(test_val)
                if(not_through_origin_test)
                    plane_normal(:, ii) = normal;
                end
                color_line_D(:, ii) = D;
                color_line_P(:, ii) = P;
                
                inlier_mask(inliers, ii) = 1;
            end
%         else % not required as the values are already 0 for others
        end
    end
    
    rmpath('RANSAC_code');
    
    assert(~any(isnan(plane_normal(:))));
    % now the return the values
    % color_line_D
    % color_line_P
    % plane_normal
    % inlier_mask
    % patch_dc
    
end
