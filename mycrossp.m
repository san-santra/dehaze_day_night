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

function [cross_p] = mycrossp(v1, v2)
% cross product with out much checks

    cross_p = [v1(2)*v2(3) - v1(3)*v2(2);
                v1(3)*v2(1) - v1(1)*v2(3);
                v1(1)*v2(2) - v1(2)*v2(1)];

end
