function [R, t, inl] = localize_ms_decoupled_3D2D(pts3D, pts2D, cam_idx, ...
    multistereo, min_inliers, depths)
% LOCALIZE_DECOUPLED_3D2D outliers rejection and position refinement for each camera, then join solutions
% Parameters
% - pts3D: 3xN matrix whose cols are 3D points in world frame
% - pts2D: 3xN matrix whose cols are 2D points (image points in hom
% coords), matched with pts3D. Left cameras perceptions.
% - cam_idx: N elements array whose ii element is the index of the stereo 
% pair perceiving the ii pts2D
% - multistereo: multistereo object
% - min_inliers: min inliers per single camera
% - depths: depth of p3D 
% Returns:
% - R, t extrinsics of principal camera in world frame
% - inl_idx: inliers matches, logical array

n = size(pts2D, 2);
n_cams = multistereo.stereoPairsNumber();

if nargin < 6 
    depths = ones(1, n);
end

% localize each camera
inl = false(1, n);
inl_n = zeros(1, n_cams);
Rs = zeros(3, 3, n_cams); 
ts = zeros(3, n_cams);
avg_depths = zeros(1, n_cams);
for ii = 1:n_cams
    idx = cam_idx == ii;
    n_cam_points = sum(idx);
    if n_cam_points > min_inliers
        % localize camera
        Kl = multistereo.getKs(ii);
        [Rii, tii, inlii] = localize_3D2D(pts3D(:, idx), pts2D(:, idx), Kl); 
        inl_n(ii) = sum(inlii);
        if inl_n(ii) > min_inliers
            % bring result in principal camera frame
            [Rc2p, tc2p] = multistereo.Tcam2principal(ii);
            Rs(:, :, ii) = Rc2p * Rii;
            ts(:, ii) = Rc2p * tii + tc2p;
            % store inliers for current camera
            inl(idx) = inlii;
            % store average depth
            depths_ = depths(idx);
            avg_depths(ii) = mean(depths_(inlii));
        end
    end
end


valid = inl_n > min_inliers;
if sum(valid) == 0
    R = nan(3);
    t = nan(3);
    inl = [];
    return;
end

Rs = Rs(:, :, valid);
ts = ts(:, valid);
inl_n = inl_n(valid);
avg_depths = avg_depths(valid);
w = inl_n ./ sqrt(avg_depths);
%w = 1 ./ sqrt(avg_depths);
%w = ones(1, sum(valid));
w = w / sum(w);

R = Rmean_chordal(Rs, w);
t = sum(w.*ts, 2);
end
