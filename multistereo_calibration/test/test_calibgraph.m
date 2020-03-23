%% compute the guess
clear
close all
clc

% change the value of images folder in calibration_parameters to try
% different sets of images
calibration_parameters


load('guidance_stereos.mat');


cams = {cam1, cam2, cam3, cam4};


g = CalibrationGraph(cams);

%%

[g, sol1] = g.Calibrate();

%%
ptsL00 = g.G.Edges.edgeData(1).pts2DL;
ptsR00 = g.G.Edges.edgeData(1).pts2DR;
pts3DL00 = g.G.Edges.edgeData(1).pts3DL;
pts3DR00 = g.G.Edges.edgeData(1).pts3DR;

ptsL01 = g.G.Edges.edgeData(2).pts2DL;
ptsR01 = g.G.Edges.edgeData(2).pts2DR;
pts3DL01 = g.G.Edges.edgeData(2).pts3DL;
pts3DR01 = g.G.Edges.edgeData(2).pts3DR;

ptsL10 = g.G.Edges.edgeData(3).pts2DL;
ptsR10 = g.G.Edges.edgeData(3).pts2DR;
pts3DL10 = g.G.Edges.edgeData(3).pts3DL;
pts3DR10 = g.G.Edges.edgeData(3).pts3DR;

ptsL11 = g.G.Edges.edgeData(4).pts2DL;
ptsR11 = g.G.Edges.edgeData(4).pts2DR;
pts3DL11 = g.G.Edges.edgeData(4).pts3DL;
pts3DR11 = g.G.Edges.edgeData(4).pts3DR;

cam1_guess = [g.G.Nodes.Tguess(2).R, g.G.Nodes.Tguess(2).t; 0, 0, 0, 1];
pat0_guess = [g.G.Nodes.Tguess(3).R, g.G.Nodes.Tguess(3).t; 0, 0, 0, 1];
pat1_guess = [g.G.Nodes.Tguess(4).R, g.G.Nodes.Tguess(4).t; 0, 0, 0, 1];

% compute the reprojection errors

% cam0
[Kl, Kr, Rs, ts] = convert_stereo_params(cam0);
% left
% guess pattern0
[proj_left00_guess, err_left00_guess] = reprojection(Kl, pat0_guess(1:3, 1:3), pat0_guess(1:3, 4), pts3DL00, ptsL00);
% guess pattern1
[proj_left01_guess, err_left01_guess] = reprojection(Kl, pat1_guess(1:3, 1:3), pat1_guess(1:3, 4), pts3DL01, ptsL01);
% optimized pattern0
sol_pat0 = sol(7:12);
pat0 = SO3R3_2_T(sol_pat0);
[proj_left00, err_left00] = reprojection(Kl, pat0(1:3, 1:3), pat0(1:3, 4), pts3DL00, ptsL00);
% optimized pattern1
sol_pat1 = sol(13:18);
pat1 = SO3R3_2_T(sol_pat1);
[proj_left01, err_left01] = reprojection(Kl, pat1(1:3, 1:3), pat1(1:3, 4), pts3DL01, ptsL01);
% right
Tr_guess = [Rs, ts; zeros(1, 3), 1] * pat0_guess;
[proj_right00_guess, err_right00_guess] = reprojection(Kr, Tr_guess(1:3, 1:3), Tr_guess(1:3, 4), pts3DR00, ptsR00);
Tr_guess = [Rs, ts; zeros(1, 3), 1] * pat1_guess;
[proj_right01_guess, err_right01_guess] = reprojection(Kr, Tr_guess(1:3, 1:3), Tr_guess(1:3, 4), pts3DR01, ptsR01);
Tr = [Rs, ts; zeros(1, 3), 1] * pat0;
[proj_right00, err_right00] = reprojection(Kr, Tr(1:3, 1:3), Tr(1:3, 4), pts3DR00, ptsR00);
Tr = [Rs, ts; zeros(1, 3), 1] * pat1;
[proj_right01, err_right01] = reprojection(Kr, Tr(1:3, 1:3), Tr(1:3, 4), pts3DR01, ptsR01);

% cam1
[Kl, Kr, Rs, ts] = convert_stereo_params(cam1);
% left
T_guess = cam1_guess \ pat0_guess;
[proj_left10_guess, err_left10_guess] = reprojection(Kl, T_guess(1:3, 1:3), T_guess(1:3, 4), pts3DL10, ptsL10);
T_guess = cam1_guess \ pat1_guess;
[proj_left11_guess, err_left11_guess] = reprojection(Kl, T_guess(1:3, 1:3), T_guess(1:3, 4), pts3DL11, ptsL11);
Tc = SO3R3_2_T(sol(1:6));
T = Tc * pat0;
[proj_left10, err_left10] = reprojection(Kl, T(1:3, 1:3), T(1:3, 4), pts3DL10, ptsL10);
T = Tc * pat1;
[proj_left11, err_left11] = reprojection(Kl, T(1:3, 1:3), T(1:3, 4), pts3DL11, ptsL11);
% right
Tr_guess = [Rs, ts; zeros(1, 3), 1] * (cam1_guess \ pat0_guess);
[proj_right10_guess, err_right10_guess] = reprojection(Kr, Tr_guess(1:3, 1:3), Tr_guess(1:3, 4), pts3DR10, ptsR10);
Tr_guess = [Rs, ts; zeros(1, 3), 1] * (cam1_guess \ pat1_guess);
[proj_right11_guess, err_right11_guess] = reprojection(Kr, Tr_guess(1:3, 1:3), Tr_guess(1:3, 4), pts3DR11, ptsR11);
Tr = [Rs, ts; zeros(1, 3), 1] * Tc * pat0;
[proj_right10, err_right10] = reprojection(Kr, Tr(1:3, 1:3), Tr(1:3, 4), pts3DR10, ptsR10);
Tr = [Rs, ts; zeros(1, 3), 1] * Tc * pat1;
[proj_right11, err_right11] = reprojection(Kr, Tr(1:3, 1:3), Tr(1:3, 4), pts3DR11, ptsR11);


disp(['Cam0 left image pattern 0: ', num2str(mean(err_left00_guess)) ' -> ' num2str(mean(err_left00))]);
disp(['Cam0 right image pattern 0: ', num2str(mean(err_right00_guess)) ' -> ' num2str(mean(err_right00))]);
disp(['Cam0 left image pattern 1: ', num2str(mean(err_left01_guess)) ' -> ' num2str(mean(err_left01))]);
disp(['Cam0 right image pattern 1: ', num2str(mean(err_right01_guess)) ' -> ' num2str(mean(err_right01))]);
disp(['Cam1 left image pattern 0: ', num2str(mean(err_left10_guess)) ' -> ' num2str(mean(err_left10))]);
disp(['Cam1 right image pattern 0: ', num2str(mean(err_right10_guess)) ' -> ' num2str(mean(err_right10))]);
disp(['Cam1 left image pattern 1: ', num2str(mean(err_left11_guess)) ' -> ' num2str(mean(err_left11))]);
disp(['Cam1 right image pattern 1: ', num2str(mean(err_right11_guess)) ' -> ' num2str(mean(err_right11))]);

