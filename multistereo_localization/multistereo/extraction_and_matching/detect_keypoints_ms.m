function [kL, dL, camidxL, kR, dR, camidxR] = detect_keypoints_ms(imagesL, imagesR)
    % DETECT_KEYPOINTS_MS
    % Parameters
    % - imagesL: concatenation over third dimension of left camera images
    % - imagesR: concatenation over third dimension of right camera images
    % Returns
    % - kL: 3xN matrix whose columns are homogeneous coordinates points of
    % left images
    % - dL: NxX descriptors where dL(ii, :) is associated to kL(:, ii)
    % - camidxL: N elements vector, camidxL(ii) is the index of the image
    % in which kL(:, ii) has been found. 
    % The image is imagesL(:, :, camidxL(ii))
    % - kR, dR, camidxR: the same as for left camera, but for right camera.
    %
    % You can avoid to pass right images and obtain only kL, dL, camidxL
    
    global SURF_THRESH;
    global SURF_OCTAVES;
    global SURF_SCALES;
    global SELECT_UNIFORM;
    global KEYPOINTS_NUMBER;
    
    n = size(imagesL, 3);
    
    % compute left features
    for ii = n:-1:1
        k = detectSURFFeatures(imagesL(:, :, ii), ...
            'MetricThreshold', SURF_THRESH, ...
            'NumOctaves', SURF_OCTAVES, ...
            'NumScaleLevels', SURF_SCALES);
        if SELECT_UNIFORM
            k = selectUniform(k, KEYPOINTS_NUMBER, size(imagesL(:, :, ii)));
        end
        [d, k] = extractFeatures(imagesL(:, :, ii), k);
        features(ii).desc = d;
        features(ii).kpts = k.Location;
        features(ii).camidx = ones(size(k, 1), 1) * ii;
    end
    kL = cat(1, features.kpts)';
    kL = [kL; ones(1, size(kL, 2))]; % homogeneous coords
    dL = cat(1, features.desc);
    camidxL = cat(1, features.camidx);
    
    if nargin > 1 % compute right features
        for ii = n:-1:1
            k = detectSURFFeatures(imagesR(:, :, ii), ...
                'MetricThreshold', SURF_THRESH, ...
                'NumOctaves', SURF_OCTAVES, ...
                'NumScaleLevels', SURF_SCALES);
            [d, k] = extractFeatures(imagesR(:, :, ii), k);
            features(ii).desc = d;
            features(ii).kpts = k.Location ;
            features(ii).camidx = ones(size(k, 1), 1) * ii;
        end
        kR = cat(1, features.kpts)';
        kR = [kR; ones(1, size(kR, 2))]; % homogeneous coords
        dR = cat(1, features.desc);
        camidxR = cat(1, features.camidx);
    end
end

