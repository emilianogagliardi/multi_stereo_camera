% try LM for extrinsic refinement by reprojection error minimization
clear
close all
clc

% read image and camera parameters
parameters
load('guidance_stereos.mat');
stereo = cam0;
[Kl, Kr, Rs, ts] = convert_stereo_params(stereo);
template = imread('pattern.png');
template = imresize(template, 0.05);
[kpts_templ, desc_templ] = affineDetect(template, 1000, true);
stereo_im = imread('/home/emiliano/Desktop/guidance/multicam_joined/0-0.png');
imL = stereo_im(:, 1:size(stereo_im, 2)/2);
imR = stereo_im(:, size(stereo_im, 2)/2+1:end);
imL = undistortImage(imL, stereo.CameraParameters1);
imR = undistortImage(imR, stereo.CameraParameters2);

% compute an homography on the left image for position guess
[kptsL, descL] = affineDetect(imL, 1000, true);
[kptsR, descR] = affineDetect(imR, 1000, true);
[Hl, ~, ptsL, pts3DL] = homog_templ2image(desc_templ, descL, kpts_templ, kptsL, template, imL);
[~, ~, ptsR, pts3DR] = homog_templ2image(desc_templ, descR, kpts_templ, kptsR, template, imR);

% compute a pose guess from the homography
[R_guess, t_guess] = localize_planar_from_H(Kl, Hl);
theta = acos((trace(R_guess)-1)/2);
w_guess = (theta/(2*sin(theta)))*[R_guess(3,2)-R_guess(2,3); R_guess(1,3)-R_guess(3,1); R_guess(2,1)-R_guess(1,2)];


%%
[e, ~, j] = err_and_jacobian_wrapper(cam0, ptsL, ptsR, pts3DL, pts3DR);

% put the pattern pose to the identity
e = @(p) e([[0, 0, 0.0001, 0, 0, 0], p]);
j = @(p) j([[0, 0, 0.0001, 0, 0, 0], p]);

global FUN;
global JACOBIAN;
FUN = e;
JACOBIAN = j;

options=optimset('Display','iter',...
                 'DerivativeCheck','off',...
                 'Jacobian','on',...
                 'MaxIter',2000, 'MaxFunEvals', 20000, 'TolX', 1e-6, 'TolFun', 1e-6, 'Algorithm', 'levenberg-marquardt'); 


sol = lsqnonlin(@(x) objective_function(x), [w_guess', t_guess'], [], [], options);
    
R = rodrigues(double(sol(1:3)));
t = double(sol(4:6))';

% left camera visualization
projected = Kl * [R, t] * [pts3DL, zeros(size(pts3DL, 1), 1), ones(size(pts3DL, 1), 1)]';
projected = projected(1:2, :) ./ projected(3, :);
projected_guess = Kl * [R_guess, t_guess] * [pts3DL, zeros(size(pts3DL, 1), 1), ones(size(pts3DL, 1), 1)]';
projected_guess = projected_guess(1:2, :) ./ projected_guess(3, :);

figure, imshow(imL);
hold on
plot(ptsL(:, 1), ptsL(:, 2), 'g+', 'MarkerSize', 1);
plot(projected(1, :), projected(2, :), 'r+', 'MarkerSize', 1);
plot(projected_guess(1, :), projected_guess(2, :), 'b+', 'MarkerSize', 1);

err_guess_left = mean(vecnorm(ptsL - projected_guess', 2, 2))
err_left = mean(vecnorm(ptsL - projected', 2, 2))

% right camera visualization
projected = Kr * [Rs, ts] * [R, t; zeros(1, 3), 1] * [pts3DR, zeros(size(pts3DR, 1), 1), ones(size(pts3DR, 1), 1)]';
projected = projected(1:2, :) ./ projected(3, :);
projected_guess = Kr * [Rs, ts] * [R_guess, t_guess; zeros(1, 3), 1] * [pts3DR, zeros(size(pts3DR, 1), 1), ones(size(pts3DR, 1), 1)]';
projected_guess = projected_guess(1:2, :) ./ projected_guess(3, :);

figure, imshow(imR);
hold on
plot(ptsR(:, 1), ptsR(:, 2), 'g+', 'MarkerSize', 1);
plot(projected(1, :), projected(2, :), 'r+', 'MarkerSize', 1);
plot(projected_guess(1, :), projected_guess(2, :), 'b+', 'MarkerSize', 1);

err_guess = mean(vecnorm(ptsR - projected_guess', 2, 2))
err = mean(vecnorm(ptsR - projected', 2, 2))


