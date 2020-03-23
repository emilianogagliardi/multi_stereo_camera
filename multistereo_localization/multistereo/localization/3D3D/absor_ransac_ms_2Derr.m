function [R, t, inl_idx] = absor_ransac_ms_2Derr(p3Dw, p3Dc, p2D, cam_idx, multistereo, w1, w2, minimize_reproj)
%MULTI_ABSOR_RANSAC_2DERR absor ransac with reprojection error evaluation for inliers selection
% Parameters:
% - p3Dw: Nx3 matrix whose rows are world points
% - p3Dc: Nx3 matrix whose rows are camera points
% - p2D: Nx2 matrix whose rows are keypoints
% - cam_idx: Nx1 array whose i-th value is the indices of the camera
% that perceived the i-th keypoint. Indices goes from 1 to n-cameras
% - multistereo: MultiStereo object
% - minimize_reproj: bool, if true use LM reprojection error minimization
% over inliers instead of point clouds alignment. In that case weights are
% not considered
% Returns inliers indices and estimation of principal camera extrinsics in
% world frame

global ABSOR_RANSAC_SAMPLE_NUMBER;
global ABSOR_RANSAC_MAX_ITERATIONS;
global ABSOR_RANSAC_CONFIDENCE;
global ABSOR_RANSAC_MAX_REPROJ_ERR;

n_points = size(p3Dw, 1);
    
if nargin == 5 || isempty(w1)
    w1 = ones(n_points, 1);
    w2 = ones(n_points, 1);
else
    w1 = w1  / sum(w1);
    w2 = w2  / sum(w2);
end

if nargin < 8 
    minimize_reproj = false;
end

% ransac parameters
params.sampleSize = ABSOR_RANSAC_SAMPLE_NUMBER;
params.recomputeModelFromInliers = false;
params.defaultModel.R = nan(3);
params.defaultModel.t = nan(3, 1);
params.maxNumTrials = ABSOR_RANSAC_MAX_ITERATIONS;
params.confidence = ABSOR_RANSAC_CONFIDENCE;
params.maxDistance = ABSOR_RANSAC_MAX_REPROJ_ERR;

% RANSAC function handles
funcs.fitFunc = @solveCameraPose;
funcs.evalFunc = @evalCameraPose;
funcs.checkFunc = @check;
points = pack(p3Dw, p3Dc, p2D, cam_idx);
%[isFound, pose, inl_idx] = vision.internal.ransac.msac(...
%     points, params, funcs, multistereo);
[isFound, pose, inl_idx] = ransac(points, params, funcs, multistereo);
if isFound 
    % compute the model from inliers
    if ~minimize_reproj
        [R, t] = absor(p3Dw(inl_idx, :)', p3Dc(inl_idx, :)', w1(inl_idx)', w2(inl_idx)');
    else
        [R, t] =  LM_reproj_multicam(p3Dw(inl_idx, :)', p2D(inl_idx, :)', cam_idx(inl_idx), multistereo, pose.R, pose.t);
    end
else
    R = nan(3);
    t = nan(3, 1);
    inl_idx = false(size(p3Dw, 1));
end
end

%--------------------------------------------------------------------------
function pose = solveCameraPose(points, varargin) % absolute orientation
[p3Dw, p3Dc, ~, ~] = unpack(points);
[R, t] = absor(p3Dw', p3Dc');
pose.R = R;
pose.t = t;
end

%--------------------------------------------------------------------------
function errs = evalCameraPose(pose, points, varargin) % reprojection error
multistereo = varargin{1};
n = size(points, 1);
errs = zeros(n, 1);
R = pose.R;
t = pose.t;
[p3Dw, ~, p2D, cam_idx] = unpack(points);
for ii = 1:n
    idx = cam_idx(ii);
    K = multistereo.getKs(idx);
    [Rc, tc] = multistereo.Tprincipal2cam(idx);
    repr = (K * [Rc * R, Rc * t + tc] * [p3Dw(ii, :), 1]')';
    repr = repr(1:2) ./ repr(3);
    errs(ii) = norm(repr - p2D(ii, :))^2;
end
end

%--------------------------------------------------------------------------
function o = pack(p3Dw, p3Dc, p2D, cam_idx)
o = [p3Dw, p3Dc, p2D, cam_idx];
end

%--------------------------------------------------------------------------
function [p3Dw, p3Dc, p2D, cam_idx] = unpack(input)
p3Dw = input(:, 1:3);
p3Dc = input(:, 4:6);
p2D = input(:, 7:8);
cam_idx = input(:, 9);
end

%--------------------------------------------------------------------------
function r = check(pose, varargin)
r = all(all(~isnan(pose.R))) && all(~isnan(pose.t));
end


