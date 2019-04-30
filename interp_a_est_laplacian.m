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

function [interp_val] = interp_a_est_laplacian(im, est_vals, conf, ...
    numsample, alpha, beta, min_i_diff_sq)
    
    max_edge_dist = 1;
    i_diff_thresh = 0.01;
    
    [nrow, ncol, ~] = size(im);
    numnode = nrow * ncol;
    ind = 1:numnode;
    ind_mat = reshape(ind, nrow, ncol);
    
    reshaped_im = reshape(permute(im, [3 1 2]), 3, nrow*ncol);
    
    addpath('RANSAC_code');
    % for myrandomsample()
    
    wt = sparse(numnode, numnode);
    % length 2^i down edge
    for ii=0:max_edge_dist
        down_neigh = ind_mat + 2^ii;
        down_neigh_excl = down_neigh(1:end -2^ii, :);
        ind_mat_excl = ind_mat(1:end -2^ii, :);

        neigh_i_diff = reshaped_im(:, ind_mat_excl(:)) - reshaped_im(:, down_neigh_excl(:));
        i_d_norm_sq = sum(neigh_i_diff.*neigh_i_diff, 1);
        down_wt = 1./max(i_d_norm_sq, min_i_diff_sq);

        down_edge = sparse(ind_mat_excl(:), down_neigh_excl(:), down_wt, numnode, numnode);

        % also include up edges
        wt = wt + down_edge + down_edge';
    end
    
    % length 2^i right edge
    for ii=0:max_edge_dist
        right_neigh = ind_mat + (2^ii)*nrow;
        right_neigh_excl = right_neigh(:, 1:end -(2^ii));
        ind_mat_excl = ind_mat(:, 1:end -(2^ii));

        neigh_i_diff = reshaped_im(:, ind_mat_excl(:)) - reshaped_im(:, right_neigh_excl(:));
        i_d_norm_sq = sum(neigh_i_diff.*neigh_i_diff, 1);
        right_wt = 1./max(i_d_norm_sq, min_i_diff_sq);

        right_edge = sparse(ind_mat_excl(:), right_neigh_excl(:), right_wt, numnode, numnode);

        % also include left edges
        wt = wt + right_edge + right_edge';
    end
    
    % Long range connections
    % edges are added between similar intensity pixels by sampling within a
    % big window
    window_size = floor(min(nrow, ncol)*0.15);
    
    patch_size = window_size;
    overlap = floor(patch_size*0.95);
    x = 1:max(1, patch_size - overlap):nrow - patch_size + 1;
    y = 1:max(1, patch_size - overlap):ncol - patch_size + 1;
    [X, Y] = meshgrid(x, y);
    
    xrow = X(:);
    ycol = Y(:);
    
    i_idx = ones(1, length(xrow)); 
    j_idx = ones(1, length(xrow));
    s_val = zeros(1, length(xrow));
    % Causing problem with 0's as i and j value. So used 1. But this won't
    % be a problem as if some s is 0 them the correspoding i and j is
    % ignored
    for ii=1:length(xrow)
        row = xrow(ii);
        col = ycol(ii);
        
        ind_patch = ind_mat(row:row+patch_size-1, col:col+patch_size-1);
        
%         sampled_ind = randsample(ind_patch(:), numsample);
        sampled_ind = ind_patch(myrandomsample(patch_size^2, numsample));
        
        % edges are added between the center pixel and one of the pixels
        % sampled from the window
        x = ind_patch(floor(numel(ind_patch)/2));
        
        % remove x and Nx from sampled_ind
        sampled_ind = sampled_ind(sampled_ind ~= x);
        sampled_ind = sampled_ind(sampled_ind ~= x+1 & sampled_ind ~= x-1 & ...
                                sampled_ind ~= x+nrow & sampled_ind ~= x-row);
        % as currently +/-1 and +/-nrow of x
        
        % if no sample left continue
        if(isempty(sampled_ind))
            continue;
        end
        
        i_diff = bsxfun(@minus, reshaped_im(:, x), reshaped_im(:, sampled_ind));
        i_d_norm_sq = sum(i_diff.*i_diff, 1);
        
        sel_mask = i_d_norm_sq < i_diff_thresh;
        
        % if no sample is selected continue
        if(sum(sel_mask) == 0)
            continue;
        end
        
        sel_idx = sampled_ind(sel_mask);
        sel_norm_sq = i_d_norm_sq(sel_mask);
        
        % In call to sparse(i,j,s,m,n,nzmax)
        % Any elements of s that are zero are ignored, along with 
        % the corresponding values of i and j.
        i_idx(ii) = x;
        j_idx(ii) = sel_idx(1);
        s_val(ii) = 1./max(sel_norm_sq(1), min_i_diff_sq);
    end
    
    wt_l = sparse(i_idx, j_idx, s_val, numnode, numnode);
    
    wt = wt + wt_l + wt_l';
    
    fprintf('# long range edges:%d\n', nnz(wt_l));
    
    % Graph laplacian    
    degree = sum(wt, 2);
    D = spdiags(degree, 0, numnode, numnode);
    
    laplacian_mat = D - wt;
    
    sigma = spdiags(conf(:), 0, numnode, numnode);
    
 
    % Little modification
    intensity = sqrt(sum(im.*im, 3));
    b = 1./intensity;
    b(intensity == 0) = 0;

    interp_val_col = (sigma + alpha*laplacian_mat)\(sigma*est_vals(:) - beta*(b(:)));
    interp_val = reshape(interp_val_col, nrow, ncol);
    
    rmpath('RANSAC_code');
end
