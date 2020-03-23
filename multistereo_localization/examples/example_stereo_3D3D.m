clear
close all
clc

% load calibration 
load('examples/data/calibration/stereos.mat'); stereo = cam4;
[Kl, Kr, Rs, ts] = convert_stereo_params(stereo);
% read images
base_imL = imread('examples/data/pose1/im4L.png');
base_imR = imread('examples/data/pose1/im4R.png');

% load configuration parameters
parameters;

% detect keypoints in left and right images
[kL, dL, kR, dR] = detect_keypoints(base_imL, base_imR);
% perform sparse reconstruction
[p3D, p3D_desc, kL] = reconstruct(kL, kR, dL, dR, stereo);

% read second pose images
imL = imread('examples/data/pose2/im4L.png');
imR = imread('examples/data/pose2/im4R.png');
% detect keypoints in left and right images
[kL2, dL2, kR2, dR2] = detect_keypoints(imL, imR);
% perform sparse reconstruction
[p3D2, p3D_desc2, kL2] = reconstruct(kL2, kR2, dL2, dR2, stereo);

% match the features
m = matchFeatures(p3D_desc, p3D_desc2, 'method', 'approximate', 'MatchThreshold', 1.2, 'MaxRatio', 0.8);
kLm = kL(:, m(:, 1));
kL2m = kL2(:, m(:, 2));
f = figure();
showMatchedFeatures(base_imL, imL, from_homogeneous(kLm)', from_homogeneous(kL2m)', 'montage');
pause;
close(f);

% localize camera (absor_ransac performs outliers rejection and global
% optimization)
p3Dm = p3D(:, m(:, 1));
p3Dm2 = p3D2(:, m(:, 2));
[R, t, inl] =  absor_ransac(p3Dm', p3Dm2');

% localize camera with absor_ransac_2Derr, that uses absor for solving
% camera pose and reprojection error for outliers rejection
[R_, t_, inl_] = absor_ransac_2Derr(p3Dm', p3Dm2', from_homogeneous(kL2m)', Kl);

% visualize inliers matches
kLm_i = kLm(:, inl);
kL2m_i = kL2m(:, inl);
f = figure();
showMatchedFeatures(base_imL, imL, from_homogeneous(kLm_i)', from_homogeneous(kL2m_i)', 'montage');
title('inliers with object space error');

kLm_i_ = kLm(:, inl_);
kL2m_i_ = kL2m(:, inl_);
f_ = figure();
showMatchedFeatures(base_imL, imL, from_homogeneous(kLm_i_)', from_homogeneous(kL2m_i_)', 'montage');