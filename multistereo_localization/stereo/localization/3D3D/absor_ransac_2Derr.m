function [R, t, inl_idx] = absor_ransac_2Derr(pc1, pc2, kpts, K, w1, w2)
%ABSOR_RANSAC_2Drejection Compute the transformation that aligns the two point cloud using absolute orientation solution quaternion inside RANSAC
% matches are considered to be inliers if the distance of the points
% in pc2 projected in camera generating pc1 from matched points in
% kpts is less than a threshold in pixel
% Parameters
% - pc1, pc2: two points clouds Nx3 whose rows are 3d points, the
% matching is given by the order of the pointclouds
% - kpts: Nx2 keypoints locations in the left camera that generated
% pc2 (matched with pc2 points, the match is given by the order)
% - K: calibration matrix of the left camera 
% - w1: estimation of the precision of points in pc1, wieghts error
% funciton
% - w2: the same for point 2. Resulting weight are w1 .* w2
% Returns
% - the rototraslation that applied to pc1 minimize the weighted square
% distances between matched points in pc1 (transformed) and pc2
% (considering only inliers)
    
    global ABSOR_RANSAC_SAMPLE_NUMBER
    global ABSOR_RANSAC_MAX_ITERATIONS;
    global ABSOR_RANSAC_CONFIDENCE;
    global ABSOR_RANSAC_MAX_REPROJ_ERR;
    
    n_points = size(pc1, 1);
    
    if nargin == 4
        w1 = ones(n_points, 1);
        w2 = ones(n_points, 1);
    else
        w1 = w1  * n_points / sum(w1);
        w2 = w2  * n_points / sum(w2);
    end
    
    % ransac parameters
    params.sampleSize = ABSOR_RANSAC_SAMPLE_NUMBER;
    params.recomputeModelFromInliers = false;
    params.defaultModel.R = nan(3);
    params.defaultModel.t = nan(1, 3);
    params.maxNumTrials = ABSOR_RANSAC_MAX_ITERATIONS;
    params.confidence = ABSOR_RANSAC_CONFIDENCE;
    params.maxDistance = ABSOR_RANSAC_MAX_REPROJ_ERR;

    % RANSAC function handles
    funcs.fitFunc = @solveCameraPose;
    funcs.evalFunc = @evalCameraPose;
    funcs.checkFunc = @check;
    
    points = pack(pc1, pc2, kpts);
    [isFound, ~, inl_idx] = vision.internal.ransac.msac(...
        points, params, funcs, K);
    
    if isFound 
        % compute model from inliers
        [R, t] = absor(pc1(inl_idx, :)', pc2(inl_idx, :)', w1(inl_idx)', w2(inl_idx)');
    else
        R = nan(3);
        t = nan(3, 1);
        inl_idx = false(size(pc1, 1));
    end
end

%--------------------------------------------------------------------------
function pose = solveCameraPose(input, varargin)

    [pc1, pc2] = unpack(input);

    [R, t] = absor(pc1', pc2');

    pose.R = R;
    pose.t = t;
end

%--------------------------------------------------------------------------
function dis = evalCameraPose(pose, points, varargin)

    [pc1, ~, kpts] = unpack(points);
    
    K = varargin{1};
    R = pose.R;
    t = pose.t;
    
    P = K * [R, t];
    projections = [pc1, ones(size(pc1, 1), 1)] * P';
    projections = projections(:, 1:2) ./ projections(:, 3);
    
    diffs = kpts - projections;
    dis = sum(diffs.^2, 2);
end


%--------------------------------------------------------------------------
function points = pack(pc1, pc2, kpts)
    points = [pc1, pc2, kpts];
end

%--------------------------------------------------------------------------
function [pc1, pc2, kpts] = unpack(input)
    pc1 = input(:, 1:3);
    pc2 = input(:, 4:6);
    kpts = input(:, 7:8);
end

%--------------------------------------------------------------------------
function r = check(pose, varargin)
    r = ~isempty(pose) && ~isempty(pose.R) && ~isempty(pose.t);
end