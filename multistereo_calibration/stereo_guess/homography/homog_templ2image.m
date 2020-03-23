function [H, inl_templ, inl_image, inl_templ_scaled] = homog_templ2image(desc_templ, ...
    desc_image, pts_templ, pts_image, template, image)
%HOMOG_TEMPL2IMAGE compute homography from template to image
    
global DEBUG_STEREO_GUESS;
global HOMOG_MAX_DISTANCE;
global HOMOG_TRIALS;
global MAX_RATIO;
global MATCH_THRESHOLD;

matches = matchFeatures(desc_templ, desc_image, 'MatchThreshold', ...
    MATCH_THRESHOLD, 'MaxRatio', MAX_RATIO, 'method', 'approximate', ...
    'unique', true);

if size(matches, 1) < 4
    H = nan;
    inl_templ = [];
    inl_image = [];
    inl_templ_scaled = [];
    return;
end

templ_matched = pts_templ(matches(:, 1), :);
templ_matched_scaled = scale_template(templ_matched, template);
im_matched = pts_image(matches(:, 2), :);

if DEBUG_STEREO_GUESS
    f = figure();
    showMatchedFeatures(template, image, templ_matched, ...
        im_matched, 'montage');
    pause
    close(f);
end

[H, inl_mask] = ...
    estimateGeometricTransform2(templ_matched_scaled, ...
    im_matched, 'projective', 'MaxDistance', ...
    HOMOG_MAX_DISTANCE, 'MaxNumTrials', HOMOG_TRIALS);

inl_templ = templ_matched(inl_mask, :);
inl_templ_scaled = templ_matched_scaled(inl_mask, :);
inl_image = im_matched(inl_mask, :);

if DEBUG_STEREO_GUESS
    f = figure;
    showMatchedFeatures(template, image, inl_templ, inl_image, 'montage');
    title('inliers matches');
    pause
    close(f);
end

H = H.T';
end

