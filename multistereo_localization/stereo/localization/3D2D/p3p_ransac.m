function [R, t, inl_idx] = p3p_ransac(kpts, pc, K)
%P3P_RANSAC Run p3p inside ransac (this is the same as estimateWorldCameraPose
% Parameters:
% - pc: point cloud Nx3 matrix whose rows are points
% - kpts: image points, Nx2, paired with points in pc (left camera)
% - K: calibration matrix
% returns R, t estimation  of camera extrinsics P = K*[R, t] and
% inliers indices in logical array

    % parameters
    global P3P_RANSAC_MAX_ITERATIONS;
    global P3P_RANSAC_CONFIDENCE;
    global P3P_RANSAC_MAX_REPROJ_ERROR;
    
    % ransac parameters
    params.sampleSize = 4;
    params.recomputeModelFromInliers = false;
    params.defaultModel.R = nan(3);
    params.defaultModel.t = nan(1, 3);
    params.maxNumTrials = P3P_RANSAC_MAX_ITERATIONS;
    params.confidence = P3P_RANSAC_CONFIDENCE;
    params.maxDistance = P3P_RANSAC_MAX_REPROJ_ERROR;
    
    
    % RANSAC function handles
    funcs.fitFunc = @solveCameraPose;
    funcs.evalFunc = @evalCameraPose;
    funcs.checkFunc = @check;
    points = pack(pc, kpts);
    [isFound, pose, inl_idx] = vision.internal.ransac.msac(...
        points, params, funcs, K);
    
    if isFound 
            R = pose.R';
            t = pose.t';
    else
        R = nan(3);
        t = nan(3, 1);
        inl_idx = false(1, size(pc, 1));
    end
end

%--------------------------------------------------------------------------
function pose = solveCameraPose(points, varargin)
% SOLVECAMERAPOSE: use solveP3P then use the 4th point to select the best solution returned by solveP3P
    
    K = varargin{1};
    [pc, kpts] = unpack(points);

    [Rs, ts] = vision.internal.calibration.solveP3P(...
        kpts(1:3, :), pc(1:3, :), K');
    if ~isempty(Rs)
        % use the 4th pair to select the best solution
        p3D = [pc(4, :), 1];
        p2D = kpts(4, :);
        n_solutions = size(ts, 1);
        errors = zeros(1, n_solutions);
        for ii = 1:n_solutions
            P = [Rs(:,:,ii); ts(ii,:)] * K';
            projectedPoint = p3D * P;
            projectedPoint = projectedPoint(1:2) ./ projectedPoint(3);
            d = p2D - projectedPoint;
            errors(ii) = d * d';
        end

        [~, idx] = min(errors);
        pose.R = Rs(:, :, idx);
        pose.t = ts(idx, :);
    else
        pose.R = nan(3);
        pose.t = nan(1, 3);
    end
end

%--------------------------------------------------------------------------
function dis = evalCameraPose(pose, points, varargin)

    [pc, kpts] = unpack(points);
    
    K = varargin{1};
    R = pose.R;
    t = pose.t;
    %pc_cam = pc * R + t;
    %depths = pc_cam(:, 3);
    %not_valid = depths < 0;
    
    P = [R; t] * K';
    projections = [pc, ones(size(pc, 1), 1)] * P;
    projections = projections(:, 1:2) ./ projections(:, 3);
    
    diffs = kpts - projections;
    dis = sum(diffs.^2, 2);
    %dis(not_valid) = inf(1, sum(not_valid));
end

%--------------------------------------------------------------------------
function points = pack(pc, kpts)
    points = [pc, kpts];
end

%--------------------------------------------------------------------------
function [pc, kpts] = unpack(input)
    pc = input(:, 1:3);
    kpts = input(:, 4:5);
end

%--------------------------------------------------------------------------
function r = check(pose, varargin)
    r = ~isempty(pose) && ~isempty(pose.R) && ~isempty(pose.t);
end
