function [kL, dL, kR, dR] = detect_keypoints(imL, imR)
    % DETECT_KEYPOINTS
    % Parameters
    % - imL: left camera image
    % - imR: right camera image
    % Returns
    % - kL: 3xN matrix whose columns are homogeneous coordinates points of
    % left images
    % - dL: NxX descriptors where dL(ii, :) is associated to kL(:, ii)
    % - kR, dR: the same as for left camera, but for right camera.
    %
    % You can avoid to pass right images and obtain only kL, dL, camidxL
    
    global SURF_THRESH;
    global SURF_OCTAVES;
    global SURF_SCALES;
    global SELECT_UNIFORM;
    global KEYPOINTS_NUMBER;
    
    % compute left features
    kL = detectSURFFeatures(imL, ...
        'MetricThreshold', SURF_THRESH, ...
        'NumOctaves', SURF_OCTAVES, ...
        'NumScaleLevels', SURF_SCALES);
    if SELECT_UNIFORM
        kL = selectUniform(kL, KEYPOINTS_NUMBER, size(imL));
    end
    [dL, kL] = extractFeatures(imL, kL);
    kL = kL.Location';
    kL = [kL; ones(1, size(kL, 2))]; % homogeneous coords
    
    if nargin > 1 % compute right features
        kR = detectSURFFeatures(imR, ...
            'MetricThreshold', SURF_THRESH, ...
            'NumOctaves', SURF_OCTAVES, ...
            'NumScaleLevels', SURF_SCALES);
        [dR, kR] = extractFeatures(imR, kR);
        kR = kR.Location';
        kR = [kR; ones(1, size(kR, 2))]; % homogeneous coords
    end
end

