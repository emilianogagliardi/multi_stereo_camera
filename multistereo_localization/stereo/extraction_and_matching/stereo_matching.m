function [kLm, kRm, dLm, dRm] = stereo_matching(dL, dR, kL, kR)
%STEREO_MATCHING
% Parameters
% - dL, dR: NxX and MxX matrices of descriptors where X is the descriptor
% dimensionality
% - kL, kR: 3xN and 3xM matrices whose columns are keypoints
% Returns
% - kLm, kRm: 3xP matrices of keypoitns, matched points on the same row
% - dLm, dRm: PxX matrices of descriptors, corresponding to kLm and kRm

global STEREO_DESC_DISTANCE;
global STEREO_MAX_RATIO;
global USE_EPIPOLAR_SEARCH;

if ~USE_EPIPOLAR_SEARCH
    
    m = matchFeatures(dL, dR, 'MatchThreshold', STEREO_DESC_DISTANCE, ...
        'MaxRatio', STEREO_MAX_RATIO, 'Method', 'Approximate');
    kLm = kL(:, m(:, 1));
    kRm = kR(:, m(:, 2));
    if nargout > 2
        dLm = dL(m(:, 1), :);
        if nargout > 3
            dRm = dR(m(:, 2), :);
        end
    end
    
else
    d = @(x, y) norm(x-y);
    [kLm, kRm, dLm, dRm] = epipolar_search2(dL, dR, from_homogeneous(kL)', ...
        from_homogeneous(kR)', d);
    kLm = to_homogeneus(kLm');
    kRm = to_homogeneus(kRm');
end
end

