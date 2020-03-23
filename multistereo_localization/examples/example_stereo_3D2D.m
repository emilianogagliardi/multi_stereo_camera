clear
close all


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

tic
[kLm, kRm, dLm, dRm] = stereo_matching(dL, dR, kL, kR);
toc
figure;
showMatchedFeatures(base_imL, base_imR, from_homogeneous(kLm)', from_homogeneous(kRm)');
pause
tic
m = matchFeatures(dL, dR,  'MatchThreshold', 1, 'MaxRatio', 1);
toc
figure;
showMatchedFeatures(base_imL, base_imR, from_homogeneous(kL(:, m(:, 1)))', from_homogeneous(kR(:, m(:, 2)))');


%%
% perform sparse reconstruction
[p3D, p3D_desc, p3D_kL] = reconstruct(kL, kR, dL, dR, stereo);

% visualize reconstruction
f = figure();
plot3(p3D(1, :), p3D(2, :), p3D(3, :), 'ro', 'MarkerSize', 1);
pause, close(f);

% read second pose image left
imL = imread('examples/data/pose2/im3L.png');

% extract keypoints 
[p2D, p2D_desc] = detect_keypoints(imL);

% match the features
m = matchFeatures(p3D_desc, p2D_desc, 'method', 'approximate', 'MatchThreshold', 1.2, 'MaxRatio', 0.8);
p3Dm = p3D(:, m(:, 1));
p3D_kLm = p3D_kL(:, m(:, 1));
p2Dm = p2D(:, m(:, 2));
f = figure();
showMatchedFeatures(base_imL, imL, from_homogeneous(p3D_kLm)', from_homogeneous(p2Dm)', 'montage');
pause;
close(f);

% localize camera
[R, t, inl, R_guess, t_guess] =  localize_3D2D(p3Dm, p2Dm, Kl);

% visualize inliers matches
p3D_kLm_i = p3D_kLm(:, inl);
p2Dm_i = p2Dm(:, inl);
f = figure();
showMatchedFeatures(base_imL, imL, from_homogeneous(p3D_kLm_i)', from_homogeneous(p2Dm_i)', 'montage');

% visualize reprojection errors
p3Dm_i = p3Dm(:, inl);
visualize_reprojection(imL, p3Dm_i, p2Dm_i, Kl, R, t, R_guess, t_guess);


