function [labels_mat, air_est, est_conf, inlier_fit_step] = ...
    aggregate_from_patches(labels, nrow, ncol, patch_size, xrow, ycol, ...
    air_est_col, est_conf_col, inlier_mask)

    labels_mat = zeros(nrow, ncol);
    air_est = zeros(nrow, ncol);
    est_conf = zeros(nrow, ncol);
    inlier_fit_step = zeros(nrow, ncol);
    
    for ii=1:length(xrow)
        row = xrow(ii);
        col = ycol(ii);
        

        air_est_patch = air_est(row:row+patch_size-1, col:col+patch_size-1);
        air_est_patch_c = reshape(air_est_col(:, ii), patch_size, patch_size);

        max_val = max(air_est_patch, air_est_patch_c);
        update_mask = (max_val == air_est_patch_c); % these indices need update

        air_est_patch(update_mask) = air_est_patch_c(update_mask);

        air_est(row:row+patch_size-1, col:col+patch_size-1) = air_est_patch;

        % corresponding confidence 
        est_conf_patch = est_conf(row:row+patch_size-1, col:col+patch_size-1);
        est_conf_col_r = reshape(est_conf_col(:, ii), patch_size, patch_size);

        est_conf_patch(update_mask) = est_conf_col_r(update_mask);
        est_conf(row:row+patch_size-1, col:col+patch_size-1) = est_conf_patch;

        % corresponding label
        labels_mat_patch = labels_mat(row:row+patch_size-1, col:col+patch_size-1);
        labels_mat_npatch = labels(ii)*ones(patch_size);

        labels_mat_patch(update_mask) = labels_mat_npatch(update_mask);

        labels_mat(row:row+patch_size -1, col:col+patch_size -1) =  labels_mat_patch;
        
        % inlier pixels from fit step
        inlier_fit_step(row:row+patch_size-1, col:col+patch_size-1) = ...
            max(inlier_fit_step(row:row+patch_size-1, col:col+patch_size-1), ...
            reshape(inlier_mask(:, ii), patch_size, patch_size));
    end
end
