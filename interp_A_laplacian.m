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

function [A_labels, influence] = interp_A_laplacian(label_est_mat, ...
    numlabel, im, lambda, min_i_diff_sq)

    max_edge_dist = 1;

    % create graph
    [nrow, ncol, ~] = size(im);
    numnode = nrow*ncol;
    ind = 1:numnode;
    ind_mat = reshape(ind, nrow, ncol);
    
    reshaped_im = reshape(permute(im, [3 1 2]), 3, numnode);
    
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
    
    % degree
    degree = sum(wt, 2);
    D = spdiags(degree, 0, numnode, numnode);
    
    nz_label = (label_est_mat ~= 0);
    Y = sparse(ind_mat(nz_label), label_est_mat(nz_label), 1, numnode, numlabel);
    
    laplacian_mat = D - wt;
    
    L = (lambda*laplacian_mat + speye(numnode, numnode))\Y;
    
    [~, out_label] = max(L, [], 2);
    A_labels = reshape(out_label, nrow, ncol);
    
    influence = L;
end
