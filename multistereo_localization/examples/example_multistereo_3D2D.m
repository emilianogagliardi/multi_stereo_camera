clear
close all
%clc

% load calibration and create multistereo object
load('examples/data/calibration/stereos.mat'); stereos = {cam1, cam2, cam3, cam4}; % no bottom camera
load('examples/data/calibration/multicams3.mat'); %multicalib = sol; clear sol;
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

% extract keypoitns in second view (only left cameras)
im1L = imread('examples/data/pose3/im1L.png');
im2L = imread('examples/data/pose3/im2L.png');
im3L = imread('examples/data/pose3/im3L.png');
im4L = imread('examples/data/pose3/im4L.png');
imagesL = cat(3, im1L, im2L, im3L, im4L);
[p2D, p2D_desc, p2D_idx] = detect_keypoints_ms(imagesL);

% match features
m = matchFeatures(p3D_desc, p2D_desc, 'method', 'approximate', 'MatchThreshold', 1.5, 'MaxRatio', 0.8);

% visualize all matches
p3D_kL_m = p3D_kL(:, m(:, 1));
p3D_idx_m = p3D_idx(m(:, 1));
p2D_m = p2D(:, m(:, 2));
p2D_idx_m = p2D_idx(m(:, 2));
fs = visualize_ms_matches(base_imagesL, imagesL, p3D_kL_m, p2D_m, p3D_idx_m, p2D_idx_m);
pause;
close(fs);

% localize the multicamera in second view with respec first view
p3D_m = p3D(:, m(:, 1));
[R, t, inl, R_guess, t_guess] = localize_ms_coupled_3D2D(p3D_m, p2D_m, p2D_idx_m, multistereo);

% localize using a decoupled approach instead
depths_m = depths(m(:, 1));
[R_, t_, inl_] = localize_ms_decoupled_3D2D(p3D_m, p2D_m, p2D_idx_m, multistereo, 5, depths_m);

% visualize inliers
fs = visualize_ms_matches(base_imagesL, imagesL, p3D_kL_m(:, inl), p2D_m(:, inl), p3D_idx_m(inl), p2D_idx_m(inl));
pause, close(fs);

% visualize inliers for the decoupled solution
fs = visualize_ms_matches(base_imagesL, imagesL, p3D_kL_m(:, inl_), p2D_m(:, inl_), p3D_idx_m(inl_), p2D_idx_m(inl_));
pause, close(fs);

% visualize reprojection error
p3D_m_i = p3D_m(:, inl);
p2D_m_i = p2D_m(:, inl);
p2D_idx_m_i = p2D_idx_m(inl);
visualize_ms_reprojection(imagesL, p3D_m_i, p2D_m_i, p2D_idx_m_i, R, t, multistereo, R_guess, t_guess);

% visualize reprojection error for the decoupled solution
fprintf('\n');
p3D_m_i = p3D_m(:, inl_);
p2D_m_i = p2D_m(:, inl_);
p2D_idx_m_i = p2D_idx_m(inl_);
visualize_ms_reprojection(imagesL, p3D_m_i, p2D_m_i, p2D_idx_m_i, R_, t_, multistereo);

