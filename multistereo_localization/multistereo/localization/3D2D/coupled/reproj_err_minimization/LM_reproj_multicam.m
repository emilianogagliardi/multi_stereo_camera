function [R, t] = LM_reproj_multicam(p3D, p2D, cam_idx, multistereo, R_guess, t_guess)
% LM_REPROJ_MULTICAM multicamera iterative extrinsic refinement
% Parameters
% - p3D: 3xN world points
% - p2D: 3xN image points homog coords
% - cam_idx: N elements vector in ii contains the index of the stereo pair
% that perceived pts2D(:, ii)
% - mutlistereo: MultiStereo object
% - R_guess, t_guess: extrinsic parameters guess of principal camera
% Returns refined extrinsics of the principal camera

global LM_TOLX;
global LM_TOLFUN;
global LM_MAX_ITERATIONS;

% obtain error and jacobian functions
e = @(p) [];
j = @(p) [];

for ii = multistereo.stereoPairsNumber():-1:1
    idx = cam_idx == ii;
    if sum(idx) > 0
        p3Dii = p3D(:, idx);
        p2Dii = p2D(:, idx);
        K = multistereo.getKs(ii);
        [R, t] = multistereo.Tprincipal2cam(ii);
        [e_ii, j_ii] = ej_wrapper_reprojection(K, R, t, p3Dii', p2Dii');
        e = @(p) [e(p); e_ii(p)];
        j = @(p) [j(p); j_ii(p)];
    end
end

w_guess = rodrigues(R_guess);

o = @(x) objective(x, e, j);
    
% solve
options=optimset('Display','off',...
                 'DerivativeCheck','off',...
                 'Jacobian','on',...
                 'MaxIter', LM_MAX_ITERATIONS, ...
                 'MaxFunEvals', Inf, ...
                 'TolX', LM_TOLX, ...
                 'TolFun', LM_TOLFUN, ...
                 'Algorithm', 'levenberg-marquardt'); 

sol = lsqnonlin(o, double([w_guess', t_guess']), [], [], options);

R = rodrigues(double(sol(1:3)));
t = double(sol(4:6))';
end


function [e, j] = objective(p, err_fun, jacob_fun)
    if nargout < 2
        e = double(err_fun(p));
    else
        e = double(err_fun(p));
        j = double(jacob_fun(p));
    end
end


