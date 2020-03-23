function [repr, errs] = reproject(P, p3D, p2D)
%REPROJECTION compute reprojected poits
% Parameters:
% - P: 3x4 projection matrix
% - p3D: 4xN or 3xN matrix whose cols are 3D points in homog coords or not
% - p2D: optional parameter, if given returns also reprojection errors
% considering p2D corresponding by order to p3D
% Returns
% - repr: reprojected points, 3xN
% - errs: reprojection errors

if size(p3D, 1) == 3
    p3D = to_homogeneus(p3D);
end

repr = P * p3D;

if nargin > 2
    repr = from_homogeneous(repr);
    if size(p2D, 1) == 3
        p2D = from_homogeneous(p2D);
    end
    errs = vecnorm(repr - p2D);
end
end

