function [R, t, ptsL, ptsR, templ_3DL, templ_3DR, err_l, err_r] = stereo_guess3(...
    stereo, pts_template, desc_template, stereo_im)
% STEREO_GUESS compute a position guess of the left camera wrt a given planar template
% Parameters
% - stereo: matlab stereo object
% - pts_template: template keypoints
% - desc_template: template keypoints descriptors
% - stereo_im: concatenation over first dimension of left and right images
% of the template taken with a stereo camera
% Returns:
% - R, t: rototraslation guess of the camera with respect to the template
% - ptsL, ptsR: inlier images points
% - pts_templ3DL, pts_templ3DL: inlier template 3D points
% - err_l, err_r: mean reprojection errors


global DEBUG_STEREO_GUESS;
global TEMPLATE;
global MIN_INLIERS_SINGLE_IMAGE;

imL = stereo_im(:, 1:size(stereo_im, 2)/2);
imR = stereo_im(:, size(stereo_im, 2)/2+1:end);

imL = undistortImage(imL, stereo.CameraParameters1);
imR = undistortImage(imR, stereo.CameraParameters2);

[Kl, Kr, Rs, ts] = convert_stereo_params(stereo);


[kptsL, descL] = detect_keypoints(imL);
[kptsR, descR] = detect_keypoints(imR);


if isempty(kptsL) || isempty(kptsR)
    R = nan;
    t = nan;
    ptsL = [];
    ptsR = [];
    templ_3DL = [];
    templ_3DR = [];
    err_l = Inf;
    err_r = Inf;
    return;
end

[~, ~, ~, ptsL, templL_scaled] = p3p_template2image(desc_template, descL, pts_template, kptsL, Kl, TEMPLATE, imL);
[~, ~, ~, ptsR, templR_scaled] = p3p_template2image(desc_template, descR, pts_template, kptsR, Kr, TEMPLATE, imR);

if isempty(ptsL) || isempty(ptsR)
    R = nan;
    t = nan;
    ptsL = [];
    ptsR = [];
    templ_3DL = [];
    templ_3DR = [];
    err_l = Inf;
    err_r = Inf;
    return;
end

templ_3DL = [templL_scaled, zeros(size(templL_scaled, 1), 1)]';
templ_3DR = [templR_scaled, zeros(size(templR_scaled, 1), 1)]';


if size(ptsL, 1) > size(ptsR, 1) % use left for refinement
    ptsL_norm = Kl \ [ptsL, ones(size(ptsL, 1), 1)]';
    ptsL_norm = ptsL_norm./ptsL_norm(3, :);
    try 
        T = MLPnP(templ_3DL, ptsL_norm);
    catch
        R = nan;
        t = nan;
        ptsL = [];
        ptsR = [];
        templ_3DL = [];
        templ_3DR = [];
        err_l = Inf;
        err_r = Inf;
        return;
    end
    R_ = T(1:3, 1:3);
    t_ = T(1:3, 4);
else % use right for refinement
    ptsR_norm = Kr \ [ptsR, ones(size(ptsR, 1), 1)]';
    ptsR_norm = ptsR_norm ./ ptsR_norm(3, :);
    try 
        T = MLPnP(templ_3DR, ptsR_norm);
    catch
        R = nan;
        t = nan;
        ptsL = [];
        ptsR = [];
        templ_3DL = [];
        templ_3DR = [];
        err_l = Inf;
        err_r = Inf;
        return;
    end
    T = [Rs, ts; zeros(1, 3), 1] \ [T; zeros(1, 3), 1];
    R_ = T(1:3, 1:3);
    t_ = T(1:3, 4);
end

if size(ptsL, 1) >= MIN_INLIERS_SINGLE_IMAGE && size(ptsR, 1) >= MIN_INLIERS_SINGLE_IMAGE
    % jointly refine
    [el, jl] = ej_wrapper_reprojection(Kl, eye(3), zeros(3, 1), ...
        [templL_scaled, zeros(size(templL_scaled, 1), 1)], ptsL);
    [er, jr] = ej_wrapper_reprojection(Kr, Rs, ts, ...
        [templR_scaled, zeros(size(templR_scaled, 1), 1)], ptsR);
    e = @(p) [el(p); er(p)];
    j = @(p) [jl(p); jr(p)];
    w_ = rodrigues(R_);
    p = [w_', t_'];
    o = @(x) lsqnonlin_helper(x, e, j);
    % solve
    options=optimset('Display','off',...
                     'DerivativeCheck','off',...
                     'Jacobian','on',...
                     'MaxIter', 100, ...
                     'MaxFunEvals', inf, ...
                     'TolX', 1e-5, ...
                     'TolFun', 1e-5, ...
                     'Algorithm', 'levenberg-marquardt'); 

    sol = lsqnonlin(o, double(p), [], [], options);
    R = rodrigues(sol(1:3));
    t = sol(4:6)';
else
    R = R_;
    t = t_;
end
    

% compute the reprojection error
Pl_ = Kl * [R_, t_];
repr_left_ = Pl_ * [templ_3DL; ones(1, size(templ_3DL, 2))];
repr_left_ = repr_left_(1:2, :) ./ repr_left_(3, :);
Pl = Kl * [R, t];
repr_left = Pl * [templ_3DL; ones(1, size(templ_3DL, 2))];
repr_left = repr_left(1:2, :) ./ repr_left(3, :);

Pr = Kr * [Rs, ts];
repr_right_ = Pr * [R_ * templ_3DR + t_; ones(1, size(templ_3DR, 2))];
repr_right_ = repr_right_(1:2, :) ./ repr_right_(3, :);
repr_right = Pr * [R * templ_3DR + t; ones(1, size(templ_3DR, 2))];
repr_right = repr_right(1:2, :) ./ repr_right(3, :);

err_l = mean(vecnorm(repr_left' - ptsL, 2, 2));
err_r = mean(vecnorm(repr_right' - ptsR, 2, 2));

if DEBUG_STEREO_GUESS
    
    fl = figure;
    imshow(imL);
    hold on
    plot(ptsL(:, 1), ptsL(:, 2), 'g+', 'MarkerSize', 1);
    plot(repr_left_(1, :), repr_left_(2, :), 'b+', 'MarkerSize', 1);
    plot(repr_left(1, :), repr_left(2, :), 'r+', 'MarkerSize', 1);
    title('Reprojection on left image');
%     
%     disp(['Reprojection error left homography localization: ' num2str(...
%         mean(vecnorm(repr_left_' - ptsL, 2, 2)))]);
%     disp(['Reprojection error left LHM: ' num2str(...
%         mean(vecnorm(repr_left' - ptsL, 2, 2)))]);
    
    pause
    
    fr = figure;
    imshow(imR);
    hold on
    plot(ptsR(:, 1), ptsR(:, 2), 'g+', 'MarkerSize', 1);
    plot(repr_right_(1, :), repr_right_(2, :), 'b+', 'MarkerSize', 1);
    plot(repr_right(1, :), repr_right(2, :), 'r+', 'MarkerSize', 1);
    title('Reprojection on right image');
%     
%     disp(['Reprojection error right homography localization: ' num2str(...
%         mean(vecnorm(repr_right_' - ptsR, 2, 2)))]);
%     disp(['Reprojection error right LHM: ' num2str(...
%         mean(vecnorm(repr_right' - ptsR, 2, 2)))]);
    
    pause
    close(fl, fr);
end
end

function [e, j] = lsqnonlin_helper(p, err_fun, jacob_fun)
    if nargout < 2
        e = double(err_fun([p, 1]));
    else
        e = double(err_fun([p, 1]));
        j = double(jacob_fun([p, 1]));
    end
end