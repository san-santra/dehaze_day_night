% RANSACFITLINE - fits line to 3D array of points using RANSAC
%
% Usage  [L, inliers] = ransacfitline(XYZ, t, feedback)
%
% This function uses the RANSAC algorithm to robustly fit a line
% to a set of 3D data points.
%
% Arguments:
%          XYZ - 3xNpts array of xyz coordinates to fit line to.
%          t   - The distance threshold between data point and the line
%                used to decide whether a point is an inlier or not.
%          feedback - Optional flag 0 or 1 to turn on RANSAC feedback
%                     information.
%
% Returns:.
%           V - Line obtained by a simple fitting on the points that
%               are considered inliers.  The line goes through the
%               calculated mean of the inlier points, and is parallel to
%               the principal eigenvector.  The line is scaled by the
%               square root of the largest eigenvalue.
%               This line is a n*2 matrix.  The first column is the
%               beginning point, the second column is the end point of the
%               line.
%           L - The two points in the data set that were found to
%               define a line having the most number of inliers.
%               The two columns of L defining the two points.
%           inliers - The indices of the points that were considered
%                     inliers to the fitted line.
%
% See also:  RANSAC, FITPLANE, RANSACFITPLANE

% Copyright (c) 2003-2006 Peter Kovesi and Felix Duvallet (CMU)
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% Aug  2006 - created ransacfitline from ransacfitplane
%             author: Felix Duvallet

function [V, L, inliers] = ransacfitline(XYZ, t, feedback)
    
    if nargin == 2
	feedback = 0;
    end
    
    [rows, npts] = size(XYZ);
    
    if rows ~=3
        error('data is not 3D');
    end
    
    if npts < 2
        error('too few points to fit line');
    end
    
    s = 2;  % Minimum No of points needed to fit a line.
        
    fittingfn = @defineline;
    distfn    = @lineptdist;
    degenfn   = @isdegenerate;

    [L, inliers] = ransac(XYZ, fittingfn, distfn, degenfn, s, t, feedback);
    
    % % % ransac can return empty when solution is not found
    if(isempty(L))
        V = [];
        return;
    end
    % % % So returns empty in case no solution is found
    
    % Find the line going through the mean, parallel to the major
    % eigenvector
%     V = fitline3d(XYZ(:, inliers));
    V = []; %not used in our implementation

    
%------------------------------------------------------------------------
% Function to define a line given 2 data points as required by
% RANSAC.

function L = defineline(X)
    L = X;
    
%------------------------------------------------------------------------
% Function to calculate distances between a line and an array of points.
% The line is defined by a 3x2 matrix, L.  The two columns of L defining
% two points that are the endpoints of the line.
%
% A line can be defined with two points as:
%        lambda*p1 + (1-lambda)*p2
% Then, the distance between the line and another point (p3) is:
%        norm( lambda*p1 + (1-lambda)*p2 - p3 )
% where
%                  (p2-p1).(p2-p3)
%        lambda =  ---------------
%                  (p1-p2).(p1-p2)
%
% lambda can be found by taking the derivative of:
%      (lambda*p1 + (1-lambda)*p2 - p3)*(lambda*p1 + (1-lambda)*p2 - p3)
% with respect to lambda and setting it equal to zero

function [inliers, L] = lineptdist(L, X, t)

    p1 = L(:,1);
    p2 = L(:,2);
    
%     npts = length(X);

    p1_m_p2 = p1 - p2;
%     deno = ((p1 - p2)'*(p1 - p2));
    deno = (p1_m_p2'*p1_m_p2);
    
%     p2_m_p3_mat2 = bsxfun(@minus, p2, X);
    p2_m_p3_mat = p2(:, ones(size(X, 2), 1)) - X;
    
    D = p2 - p1;
%     tmp2 = bsxfun(@times, D, p2_m_p3_mat);
    tmp = D(:, ones(size(p2_m_p3_mat, 2), 1)).*p2_m_p3_mat; 
    
    dot_p = sum(tmp, 1);
    lambda_mat = dot_p ./ deno;
%     lambda_mat = repmat(lambda_mat, size(p1, 1), 1);
    lambda_mat = lambda_mat(ones(size(p1, 1), 1), :);

%     tmp2 = bsxfun(@times, lambda_mat, p1) + bsxfun(@times, 1-lambda_mat, p2) - X;
    tmp = lambda_mat.*p1(:, ones(size(lambda_mat, 2), 1)) + ...
            (1-lambda_mat).*p2(:, ones(size(lambda_mat, 2), 1)) - X;

    
    d_alt = sqrt(sum(tmp.*tmp, 1));

% %     d = zeros(npts, 1);
% %     for i = 1:npts
% %         p3 = X(:,i);
% %       
% % %         lambda = dot((p2 - p1), (p2-p3)) / dot( (p1-p2), (p1-p2) );
% %         lambda = ((p2 - p1)'*(p2 - p3)) / deno;
% %         
% %         d(i) = norm(lambda*p1 + (1-lambda)*p2 - p3);
% %     end

% % %   To test the new implementation without for loop
% %     if(mean(abs(d - d_alt')) > 0.0001)
% %         p2;
% %     end
% %     assert(mean(abs(d - d_alt')) < 0.0001);
    
    d = d_alt;
    inliers = find(abs(d) < t);
    
%------------------------------------------------------------------------
% Function to determine whether a set of 2 points are in a degenerate
% configuration for fitting a line as required by RANSAC.
% In this case two points are degenerate if they are the same point
% or if they are exceedingly close together.

function r = isdegenerate(X)
    %find the norm of the difference of the two points
    % this will be 0 iff the two points are the same (the norm of their
    % difference is zero)
    r = norm(X(:,1) - X(:,2)) < eps;
