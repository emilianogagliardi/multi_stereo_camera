function [R, t, inl, R_guess, t_guess, errs] = localize_ms_coupled_3D2D(pts3D, pts2D, cam_idx, ...
    multistereo)
% LOCALIZE_COUPLED_3D2D outliers rejection and position refinement 3D-2D coupled solution
% Parameters
% - pts3D: 3xN matrix whose cols are 3D points in world frame
% - pts2D: 3xN matrix whose cols are 2D points (image points in hom
% coords), matched with pts3D. Left cameras perceptions.
% - cam_idx: N elements array whose ii element is the index of the stereo 
% pair perceiving the ii pts2D
% - multistereo: multistereo object
% Returns:
% - R, t extrinsics of principal camera in world frame
% - inl_idx: inliers matches, logical array

n = size(pts2D, 2);
if n < 4
    R = nan(3);
    t = nan(3, 1);
    inl = [];
    R_guess = nan(3);
    t_guess = nan(3, 1);
    errs = nan;
    return
end

% prepare data for gp3p ransac (opengv data style)
pts2D_norm = pts2D;
ts = zeros(3, n);
for ii = 1:multistereo.stereoPairsNumber()
    idx = cam_idx == ii;
    n_cam_points = sum(idx);
    if n_cam_points > 0
        [R, t] = multistereo.Tcam2principal(ii);
        K = multistereo.getKs(ii);
        pts2D_norm(:, idx) = K \ pts2D_norm(:, idx);
        pts2D_norm(:, idx) = pts2D_norm(:, idx) ./ vecnorm(pts2D_norm(:, idx));
        pts2D_norm(:, idx) = R * pts2D_norm(:, idx);
        ts(:, idx) = repmat(t, 1, n_cam_points);
    end
end

% outliers rejection
[R_guess, t_guess, inl] = gp3p_ransac(pts2D_norm', pts3D', ts');

if sum(inl) < 3
    R = nan(3);
    t = nan(3, 1);
    inl = [];
    R_guess = nan(3);
    t_guess = nan(3, 1);
    errs = nan;
    return
end

% minimize reprojection error over all inliers
inl2D = pts2D(:, inl);
inl3D = pts3D(:, inl);
inl_idx = cam_idx(inl);
[R, t] = LM_reproj_multicam(inl3D, inl2D, inl_idx, multistereo, R_guess, t_guess);


if nargout > 5
    p3D = R * inl3D + t;
    errs = [];
    for ii = 1:multistereo.stereoPairsNumber()
        P = multistereo.projMatrix(ii);
        p3D_ = p3D(:, inl_idx == ii);
        p2D_ = inl2D(:, inl_idx == ii);
        [~, e] = reproject(P, p3D_, p2D_);
        errs = [e, errs];
    end
end

end
