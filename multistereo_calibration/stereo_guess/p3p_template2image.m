function [R, t, inl_templ, inl_image, inl_templ_scaled] = p3p_template2image(desc_templ, ...
    desc_image, pts_templ, pts_image, K, template, image)

global DEBUG_STEREO_GUESS;
global MAX_RATIO;
global MATCH_THRESHOLD;

matches = matchFeatures(desc_templ, desc_image, 'MatchThreshold', ...
    MATCH_THRESHOLD, 'MaxRatio', MAX_RATIO, 'method', 'approximate', ...
    'unique', true);

if size(matches, 1) < 4
    R = nan;
    t = nan;
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

[R, t, inl_mask] = p3p_ransac(im_matched, [templ_matched_scaled, zeros(size(templ_matched_scaled, 1), 1)], K);

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
end

