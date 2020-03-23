function [R, t, inl_idx] = absor_ransac(pc1, pc2, w1, w2)
%ABSOR_RANSAC Compute the transformation that aligns the two point cloud using absolute orientation solution quaternion inside RANSAC
% matches are considered to be inliers if the distance of 3D points in
% a pair is less than a threshold
% Parameters
% - pc1, pc2: two points clouds Nx3 whose rows are 3d points, the
% matching is given by the order of the pointclouds
% - w1: estimation of the precision of points in pc1, wieghts error
% funciton
% - w2: the same for point 2. Resulting weight are product of w1 and w2
% Returns
% - the rototraslation that applied to pc1 minimize the weighted square
% distances between matched points in pc1 (transformed) and pc2
% (considering only inliers)
    
    global ABSOR_RANSAC_SAMPLE_NUMBER;
    global ABSOR_RANSAC_MAX_ITERATIONS;
    global ABSOR_RANSAC_CONFIDENCE;
    global ABSOR_RANSAC_MAX_3D_DISTANCE;
    
    n_points = size(pc1, 1);
    
    if nargin == 2
        w1 = ones(n_points, 1);
        w2 = ones(n_points, 1);
    else
        w1 = w1  / sum(w1);
        w2 = w2  / sum(w2);
    end
    
    % ransac parameters
    params.sampleSize = ABSOR_RANSAC_SAMPLE_NUMBER;
    params.recomputeModelFromInliers = false;
    params.defaultModel.R = nan(3);
    params.defaultModel.t = nan(1, 3);
    params.maxNumTrials = ABSOR_RANSAC_MAX_ITERATIONS;
    params.confidence = ABSOR_RANSAC_CONFIDENCE;
    params.maxDistance = ABSOR_RANSAC_MAX_3D_DISTANCE;

    % RANSAC function handles
    funcs.fitFunc = @solveCameraPose;
    funcs.evalFunc = @evalCameraPose;
    funcs.checkFunc = @check;
    
    points = pack(pc1, pc2);
    [isFound, ~, inl_idx] = vision.internal.ransac.msac(...
        points, params, funcs);
    if isFound 
        % compute the model from inliers
        [R, t] = absor(pc1(inl_idx, :)', pc2(inl_idx, :)', w1(inl_idx)', w2(inl_idx)');
    else
        R = nan(3);
        t = nan(3, 1);
        inl_idx = false(size(pc1, 1));
    end
end

%--------------------------------------------------------------------------
function pose = solveCameraPose(input)

    [pc1, pc2] = unpack(input);
    
    [R, t] = absor(pc1', pc2');
    pose.R = R;
    pose.t = t;
end

%--------------------------------------------------------------------------
function dis = evalCameraPose(pose, points)

    [pc1, pc2] = unpack(points);
    
    R = pose.R;
    t = pose.t;
    
    aligned_pc1 = (R * pc1' + t)';
    diffs = (aligned_pc1 - pc2);
    dis = sum(diffs.^2, 2);
end


%--------------------------------------------------------------------------
function points = pack(pc1, pc2)
    points = [pc1, pc2];
end

%--------------------------------------------------------------------------
function [pc1, pc2] = unpack(input)
    pc1 = input(:, 1:3);
    pc2 = input(:, 4:6);
end

%--------------------------------------------------------------------------
function r = check(pose, varargin)
    r = ~isempty(pose) && ~isempty(pose.R) && ~isempty(pose.t);
end