clear
close all
%clc

% load calibration and create multistereo object
load('examples/data/calibration/stereos.mat'); stereos = {cam1, cam2, cam3, cam4}; % no bottom camera
load('examples/data/calibration/multicams.mat'); multicalib = sol; clear sol;
multistereo = MultiStereo(stereos, multicalib);

% read images
base_im1L = imread('examples/data/pose4/im1L.png');
base_im2L = imread('examples/data/pose4/im2L.png');
base_im3L = imread('examples/data/pose4/im3L.png');
base_im4L = imread('examples/data/pose4/im4L.png');
base_imagesL = cat(3, base_im1L, base_im2L, base_im3L, base_im4L);
base_im1R = imread('examples/data/pose4/im1R.png');
base_im2R = imread('examples/data/pose4/im2R.png');
base_im3R = imread('examples/data/pose4/im3R.png');
base_im4R = imread('examples/data/pose4/im4R.png');
base_imagesR = cat(3, base_im1R, base_im2R, base_im3R, base_im4R);

% load configuration parameters
parameters;

% detect keypoints in all images
[kL, dL, camidxL, kR, dR, camidxR] = detect_keypoints_ms(base_imagesL, base_imagesR);
% perform 3D reconstruction, get p3D, associated desc, camera index, and
% left camera image point coordinates
[p3D, p3D_desc, p3D_idx, p3D_kL, depths] = reconstruct_ms(kL, kR, dL, dR, camidxL, camidxR, multistereo);
% visualize reconstruction
f = visualize_ms_reconstruction(p3D, p3D_idx, multistereo);
pause;
close(f);

% extract keypoints in second view 
im1L = imread('examples/data/pose3/im1L.png');
im2L = imread('examples/data/pose3/im2L.png');
im3L = imread('examples/data/pose3/im3L.png');
im4L = imread('examples/data/pose3/im4L.png');
imagesL = cat(3, im1L, im2L, im3L, im4L);
im1R = imread('examples/data/pose3/im1R.png');
im2R = imread('examples/data/pose3/im2R.png');
im3R = imread('examples/data/pose3/im3R.png');
im4R = imread('examples/data/pose3/im4R.png');
imagesR = cat(3, im1R, im2R, im3R, im4R);
[kL2, dL2, camidxL2, kR2, dR2, camidxR2] = detect_keypoints_ms(imagesL, imagesR);
% perform reconstruction in second view
[p3D2, p3D_desc2, p3D_idx2, p3D_kL2, depths2] = reconstruct_ms(kL2, kR2, dL2, dR2, camidxL2, camidxR2, multistereo);

% match features
m = matchFeatures(p3D_desc, p3D_desc2, 'method', 'approximate', 'MatchThreshold', 1.5, 'MaxRatio', 0.8);

% visualize all matches
p3D_kL_m = p3D_kL(:, m(:, 1));
p3D_idx_m = p3D_idx(m(:, 1));
p3D_kL2_m = p3D_kL2(:, m(:, 2));
p3D_idx2_m = p3D_idx2(m(:, 2));
fs = visualize_ms_matches(base_imagesL, imagesL, p3D_kL_m, p3D_kL2_m, p3D_idx_m, p3D_idx2_m);
pause;
close(fs);

% localize the multicamera in second view with respec fist view
p3D_m = p3D(:, m(:, 1));
p3D2_m = p3D2(:, m(:, 2));
depths_m = depths(m(:, 1));
depths2_m = depths2(m(:, 2));
%[R, t, inl] = absor_ransac(p3D_m', p3D2_m');%, 1./real(sqrt(depths_m))', 1./real(sqrt(depths2_m))');

% localize using reprojection error in outliers rejection 
[R_, t_, inl_] = absor_ransac_ms_2Derr(p3D_m', p3D2_m', from_homogeneous(p3D_kL2_m)', p3D_idx2_m, multistereo);%, 1./real(sqrt(depths_m))', 1./real(sqrt(depths2_m))');

% visualize inliers
%fs = visualize_ms_matches(base_imagesL, imagesL, p3D_kL_m(:, inl), p3D_kL2_m(:, inl), p3D_idx_m(inl), p3D_idx2_m(inl));
%pause, close(fs);

% visualize inliers with reprojection error in outliers rejection
fs = visualize_ms_matches(base_imagesL, imagesL, p3D_kL_m(:, inl_), p3D_kL2_m(:, inl_), p3D_idx_m(inl_), p3D_idx2_m(inl_));
pause, close(fs);

% visualize aligned point clouds
f = figure();
p3D_ = R_ * p3D + t_;
plot3(p3D_(1, :), p3D_(2, :), p3D_(3, :), 'ro', 'MarkerSize', 1);
hold on
plot3(p3D2(1, :), p3D2(2, :), p3D2(3, :), 'bo', 'MarkerSize', 1);
axis equal
