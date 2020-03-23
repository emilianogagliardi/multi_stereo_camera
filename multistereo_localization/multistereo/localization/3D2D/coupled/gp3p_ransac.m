function [R, t, inl_idx] = gp3p_ransac(p2D, p3D, ts)
% GP3P_RANSAC ransac noncentral p3p
% Parameters:
% - p3D: Nx3 matrix whose rows are world points
% - p2D: Nx3 matrix whose rows are normalized viewing rays 
% Rc'*(Kc\im_pt) where Rc is the position of camera c in multicamera body
%frame, Kc its calibration matrix, and im_pt the homogeneaous image point
% detected by cam c corresponding to a 3D point
% - ts: Nx3 matrix whose row ii is the translation of the camera that
% perceived the point in position ii. If Pc = Kc*[Rc, tc], ts contains
% -Rc'tc

% parameters
global GP3P_RANSAC_MAX_ITERATIONS;
global GP3P_RANSAC_CONFIDENCE;
global GP3P_RANSAC_MAX_REPROJ_ERROR;

% ransac parameters
params.sampleSize = 4;
params.recomputeModelFromInliers = false;
params.defaultModel.R = nan(3);
params.defaultModel.t = nan(3, 1);
params.maxNumTrials = GP3P_RANSAC_MAX_ITERATIONS;
params.confidence = GP3P_RANSAC_CONFIDENCE;
params.maxDistance = GP3P_RANSAC_MAX_REPROJ_ERROR;


% RANSAC function handles
funcs.fitFunc = @solveCameraPose;
funcs.evalFunc = @evalCameraPose;
funcs.checkFunc = @check;
points = pack(p3D, p2D, ts);

%[isFound, pose, inl_idx] = vision.internal.ransac.msac(...
%    points, params, funcs);
[isFound, pose, inl_idx] = ransac(points, params, funcs);
if isFound 
        R = pose.R;
        t = pose.t;
else
    R = nan(3);
    t = nan(3, 1);
    inl_idx = false(size(p2D, 1));
end
end

%--------------------------------------------------------------------------
function pose = solveCameraPose(points)
[p3D, p2D, ts] = unpack(points);
X = opengv('gp3p', double(p3D'), double([p2D';ts']));
if all(all(all(~isnan(X))))
    % use the 4th pair to select the best solution
    p3D_ = p3D(4, :)';
    p2D_ = p2D(4, :)';
    t_ = ts(4, :)';
    n_solutions = size(X, 3);
    errors = zeros(1, n_solutions);
    for ii = 1:n_solutions
        R = X(1:3, 1:3, ii)';
        t = - R * X(1:3, 4, ii);
        proj = R * p3D_ + t - t_;
        proj = proj / norm(proj);
        errors(ii) = 1 - p2D_' * proj;
    end

    [~, idx] = min(errors);
    pose.R = X(1:3, 1:3, idx)';
    pose.t = -pose.R * X(1:3, 4, idx);
else
    pose.R = nan(3);
    pose.t = nan(1, 3);
end
end

%--------------------------------------------------------------------------
function errs = evalCameraPose(pose, points, varargin)

[p3D, p2D, ts] = unpack(points);

R = pose.R;
t = pose.t;
proj = R * p3D' + t - ts';
norms = sqrt(sum(proj.*proj));
proj = proj./repmat(norms, 3, 1);
errs = 1 - dot(proj, p2D');
end

%--------------------------------------------------------------------------
function o = pack(p3D, p2D, ts)
o = [p3D, p2D, ts];
end

%--------------------------------------------------------------------------
function [p3D, p2D, ts] = unpack(input)
p3D = input(:, 1:3);
p2D = input(:, 4:6);
ts = input(:, 7:9);
end

%--------------------------------------------------------------------------
function r = check(pose, varargin)
r = all(all(~isnan(pose.R))) && all(~isnan(pose.t));
end


