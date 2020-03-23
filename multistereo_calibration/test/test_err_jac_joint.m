%% compute the guess
clear
close all
clc

% change the value of images folder in calibration_parameters
calibration_parameters

load('guidance_stereos.mat');
cams = {cam0, cam1};

g = CalibrationGraph(cams);

%% refine

% obtain error function for image taken by cam0
ptsL0 = g.G.Edges.edgeData(1).pts2DL;
ptsR0 = g.G.Edges.edgeData(1).pts2DR;
pts3DL0 = g.G.Edges.edgeData(1).pts3DL;
pts3DR0 = g.G.Edges.edgeData(1).pts3DR;
[e0, jp0] = err_and_jacobian_wrapper_cam0(cam0, ptsL0, ptsR0, pts3DL0, pts3DR0);
wp_guess = rodrigues(g.G.Nodes.Tguess(3).R);
tp_guess = g.G.Nodes.Tguess(3).t;

ptsL1 = g.G.Edges.edgeData(2).pts2DL;
ptsR1 = g.G.Edges.edgeData(2).pts2DR;
pts3DL1 = g.G.Edges.edgeData(2).pts3DL;
pts3DR1 = g.G.Edges.edgeData(2).pts3DR;
[e1, jc1, jp1] = err_and_jacobian_wrapper(cam1, ptsL1, ptsR1, pts3DL1, pts3DR1);
wc1_guess = rodrigues(g.G.Nodes.Tguess(2).R');
tc1_guess = -g.G.Nodes.Tguess(2).R' * g.G.Nodes.Tguess(2).t;

% two poses to be optimized: pattern and cam1
% a variable vector has the form
% wcx wcy wcz tcx tcy tcz wpx wpy wpz tcx tcy tcz

global FUN;
FUN = @(p) [e0(p(7:12)); e1(p)];

global JACOBIAN;
JACOBIAN = @(p) [zeros(2*(size(pts3DL0, 1)+size(pts3DR0, 1)), 6) , jp0(p(7:12));
                 jc1(p), jp1(p)];

guess = [wc1_guess', tc1_guess', wp_guess', tp_guess'];

options=optimset('Display','iter',...
                 'DerivativeCheck','off',...
                 'Jacobian','on',...
                 'MaxIter',2000, 'MaxFunEvals', 20000, 'TolX', 1e-6, 'TolFun', 1e-6, 'Algorithm', 'levenberg-marquardt'); 


sol = lsqnonlin(@(x) objective_function(x), guess, [], [], options);

%% compute the errors
Rc1inv = rodrigues(sol(1:3));
tc1inv = sol(4:6)';
Rp = rodrigues(sol(7:9));
tp = sol(10:12)';


% cam0
[Kl, Kr, Rs, ts] = convert_stereo_params(cam0);
% guess
Rp_guess = g.G.Nodes.Tguess(3).R;
tp_guess = g.G.Nodes.Tguess(3).t;
[proj_left0_guess, err_left0_guess] = reprojection(Kl, Rp_guess, tp_guess, pts3DL0, ptsL0);
Tr_guess = [Rs, ts; 0, 0, 0, 1] * [Rp_guess, tp_guess; 0, 0, 0, 1];
[proj_right0_guess, err_right0_guess] = reprojection(Kr, Tr_guess(1:3, 1:3), Tr_guess(1:3, 4), pts3DR0, ptsR0);
% optimized
[proj_left0, err_left0] = reprojection(Kl, Rp, tp, pts3DL0, ptsL0);
Tr = [Rs, ts; 0, 0, 0, 1] * [Rp, tp; 0, 0, 0, 1];
[proj_right0, err_right0] = reprojection(Kr, Tr(1:3, 1:3), Tr(1:3, 4), pts3DR0, ptsR0);

% cam1
[Kl, Kr, Rs, ts] = convert_stereo_params(cam1);
% guess
Rc_guess = g.G.Nodes.Tguess(2).R;
tc_guess = g.G.Nodes.Tguess(2).t;
Tl_guess = [Rc_guess, tc_guess; 0, 0, 0, 1] \ [Rp_guess, tp_guess; 0, 0, 0, 1];
[proj_left1_guess, err_left1_guess] = reprojection(Kl, Tl_guess(1:3, 1:3), Tl_guess(1:3, 4), pts3DL1, ptsL1);
Tr_guess = [Rs, ts; 0, 0, 0, 1] * Tl_guess;
[proj_right1_guess, err_right1_guess] = reprojection(Kr, Tr_guess(1:3, 1:3), Tr_guess(1:3, 4), pts3DR1, ptsR1);
% optimized
Tl = [Rc1inv, tc1inv; 0, 0, 0, 1] * [Rp, tp; 0, 0, 0, 1];
[proj_left1, err_left1] = reprojection(Kl, Tl(1:3, 1:3), Tl(1:3, 4), pts3DL1, ptsL1);
Tr = [Rs, ts; 0, 0, 0, 1] * Tl;
[proj_right1, err_right1] = reprojection(Kr, Tr(1:3, 1:3), Tr(1:3, 4), pts3DR1, ptsR1);

disp(['Cam0 error guess left image: ', num2str(mean(err_left0_guess))]);
disp(['Cam0 error optimized left image: ', num2str(mean(err_left0))]);
disp(['Cam0 error guess right image: ', num2str(mean(err_right0_guess))]);
disp(['Cam0 error optimized right image: ', num2str(mean(err_right0))]);
disp(['Cam1 error guess left image: ', num2str(mean(err_left1_guess))]);
disp(['Cam1 error optimized left image: ', num2str(mean(err_left1))]);
disp(['Cam1 error guess right image: ', num2str(mean(err_right1_guess))]);
disp(['Cam1 error optimized right image: ', num2str(mean(err_right1))]);

disp(['Total error guess: ', num2str(mean([err_left0_guess, err_left1_guess, err_right0_guess, err_right1_guess]))]);
disp(['Total error optimizrd: ', num2str(mean([err_left0, err_left1, err_right0, err_right1]))]);
