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


function [labels] = assign_normals(normals, atmlight_dir_mat)
    
    normal_light_error = pdist2(normals', atmlight_dir_mat', 'cosine');
    % pdist2 returns (1 - cos(\theta))
    normal_light_error = abs(normal_light_error - 1);
    
    [~, labels] = min(normal_light_error, [], 2);
end
