function [kpts, desc] = detect_keypoints(im, is_pattern)
global AFFINE_INVARIANT_FEATURES;
global SURF_METRIC_THRESHOLD;
global NUM_OCTAVES_PATTERN;
global NUM_SCALE_LEVELS_PATTERN;
if nargin == 1
    is_pattern = false;
end
if AFFINE_INVARIANT_FEATURES
    [kpts, desc] = affineDetect(im, SURF_METRIC_THRESHOLD, true);
else
    if is_pattern
        kpts = detectSURFFeatures(im, 'MetricThreshold', SURF_METRIC_THRESHOLD, 'NumScaleLevels', NUM_SCALE_LEVELS_PATTERN, 'NumOctaves', NUM_OCTAVES_PATTERN);
    else
        kpts = detectSURFFeatures(im, 'MetricThreshold', SURF_METRIC_THRESHOLD);
    end
    [desc, kpts] = extractFeatures(im, kpts);
    kpts = kpts.Location;
end
end

