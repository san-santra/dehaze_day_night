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

function [labels, atmlight_dir_mat] = compute_A_dir_mat_hough( ...
    plane_normal, patch_dc, dc_thr_frac, t_step, hough_thr_frac, display)

    labels = zeros(1, size(plane_normal, 2));

    patch_idx = 1:size(plane_normal, 2);
    
    nz_col = (plane_normal(1, :) ~= 0) & (plane_normal(2, :) ~= 0) & ...
        (plane_normal(3, :) ~= 0);
    
    nz_normal = plane_normal(:, nz_col);
    nz_normal_patch_dc = patch_dc(:, nz_col);
    nz_patch_idx = patch_idx(nz_col);
    
    % remove based on patch_dc value
    sel_mask = nz_normal_patch_dc > dc_thr_frac*max(nz_normal_patch_dc);
    
    sel_normals = nz_normal(:, sel_mask);
    sel_patch_idx = nz_patch_idx(sel_mask);
    
    % apply hough transform to get A_dirs
    % plane eqn -> x_i*cos(t)*sin(p) + y_i*sin(t)*sin(p) + z_i*cos(p) = 0
    % we have (x_i, y_i, z_i) from the normals and will vary t
    % t and p is know to vary between 0 and 90
    t = 0:t_step:90; % taking degrees, easier to index
    hough_accum = zeros(length(t), length(t)); 
    
    for ii=1:size(sel_normals, 2)
        for t_idx = 1:length(t)
            t_val = t(t_idx);
            
            tan_p = -sel_normals(3, ii)/(sel_normals(1, ii)*cosd(t_val) ...
                + sel_normals(2, ii)*sind(t_val));
            p = atand(tan_p);
            
            p_idx = floor((round(p) - min(t))/t_step) +1;
            if(p>=0 && p<=90)
                hough_accum(t_idx, p_idx) = hough_accum(t_idx, p_idx) +1;
            end
        end
    end
    
    peak = max(hough_accum(:));
    
    thr_out  = hough_accum > hough_thr_frac*peak;
    sel_idx = imregionalmax(hough_accum.*thr_out);
    [theta_idx, phi_idx] = find(sel_idx);
    
    % get the value
    theta = t(theta_idx);
    phi = t(phi_idx);
    
    x = cosd(theta).*sind(phi);
    y = sind(theta).*sind(phi);
    z = cosd(phi);
    
    atmlight_dir_mat = cat(1, x, y, z);
    
    % recompute A_dir_mat from the assigned normals
    done = false;
    while ~done
        done = true;
        
        prev_atmlight_dir_mat = atmlight_dir_mat;
        normal_adir_idx = assign_normals(sel_normals, atmlight_dir_mat);
        labels(sel_patch_idx) = normal_adir_idx;
        
        for ii=1:size(atmlight_dir_mat, 2)
            Aii_normal = plane_normal(:, labels ==ii);

            if(~isempty(Aii_normal))
                [vec, val] = eig(cov(Aii_normal'))

                [~, idx] = sort(diag(val), 'ascend');
                vec = vec(:, idx);

                a_dir = vec(:, 1);

                if(all(a_dir < 0))
                    a_dir = -a_dir;
                end
            end

    %         assert(all(a_dir > 0));
            if(isempty(Aii_normal) || any(a_dir <= 0))
                % bad A_dir. Discard and assign these normals to other a_dir
                % can just remove this a_dir and assign all the normals to some
                % a_dir or just assign the these removed normals. If we
                % assign just the removed normals then need to think about
                % label shifting to account for the deleted label
                
                % discard modification made to atmlight_dir_mat
                % just update atmlight_dir_mat. the assignment will be done
                % when the while loop repeats.
                keep_idx = true(1, size(atmlight_dir_mat, 2));
                keep_idx(ii) = false;
                atmlight_dir_mat = prev_atmlight_dir_mat(:, keep_idx);

                done = false;
                break;
            end
            atmlight_dir_mat(:, ii) = a_dir;
        end
    end

    atmlight_dir_mat
% One component negative in a_dir means fitting problem. What should be done ?
% may change the fitting function. But in this case the negative component
% may just go to zero. Better reassign the normals. 
    
    if(display)
        figure(1);
        imagesc(hough_accum);
        xlabel('Hough space accumulator');
        figure(3);
        imagesc(sel_idx);
        xlabel('Detected peaks in Hough space');
        figure(2);
        clf(figure(2));
        hold on;
        for ii=1:size(atmlight_dir_mat, 2)
            color_vec = rand(1, 3);

            x = atmlight_dir_mat(1, ii);
            y = atmlight_dir_mat(2, ii);
            z = atmlight_dir_mat(3, ii);

            sel_normal = plane_normal(:, labels == ii);
            scatter3(sel_normal(1, :), sel_normal(2, :), sel_normal(3, :), 'MarkerEdgeColor', color_vec);
            line([0, x], [0, y], [0, z], 'color', color_vec);
        end
        xlabel('Normals and corresponding A_{dir}')
        hold off;
    end
end
