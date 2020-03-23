function [R, t, inl, R_guess, t_guess] = localize_3D2D(pts3D, pts2D, K)
%LOCALIZE_3D2D outliers rejection and refinement on 3D2D correspondences
% Parameters
% - pts3D: world points, 3xN or 4xN (homog coords) matrix
% - pts2D: image points, 2xN or 3xN (homog coords) matrix
% - K: calibration matrix
% Returns:
% - R, t: refined extrinsics
% - inl: logical array of inliers indices
% - R_guess, t_guess: extrinsics computed in outliers rejection


global REPROJ_MINIMIZ_METHOD;

if size(pts3D, 1) == 4
    pts3D = from_homogeneous(pts3D);
end
if size(pts2D, 1) == 3
    pts2D = from_homogeneous(pts2D);
end

if size(pts3D, 2) < 4
    R = nan(3);
    t = nan(3, 1);
    inl = zeros(1, size(pts3D, 2));
    return;
end

[R_guess, t_guess, inl] = p3p_ransac(pts2D', pts3D', K);

if sum(inl) < 0
    R = nan(3);
    t = nan(3, 1);
    inl = zeros(1, size(pts3D, 2));
    return;
end

if strcmp(REPROJ_MINIMIZ_METHOD, 'LM')
    [R, t] = LM_reproj(pts3D(:, inl), pts2D(:, inl), K, R_guess, t_guess);
elseif strcmp(REPROJ_MINIMIZ_METHOD, 'LHM')
    pts2D_norm = K \ to_homogeneus(pts2D);
    pts2D_norm = pts2D_norm ./ pts2D_norm(3, :);
    [R, t] = LHM(pts3D(:, inl), pts2D_norm(:, inl), R_guess);
else 
    error('Specify extrinsics refinement method in parameters.m as LM or LHM');
end

