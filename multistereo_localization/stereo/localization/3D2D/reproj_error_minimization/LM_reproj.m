function [R, t] = LM_reproj(pts3D, pts2D, K,...
    R_guess, t_guess)
%LM_reproj iterative extrinsics refinement
% Parameters:
% - K: calibration matrix
% - pts3D: world points, 4xN homog coords or 3xN matrix
% - pts2D: image points 3xN homog coords or 2xN matrix
% - R_guess, t_guess: extrinsic guess, starting solution for iterative
% refinement
% Returns:
% - R, t: refined extrinsics

global LM_TOLX;
global LM_TOLFUN;
global LM_MAX_ITERATIONS;

if size(pts2D, 1) == 3
    pts2D = from_homogeneous(pts2D);
end
if size(pts3D, 1) == 4
    pts3D = from_homogeneous(pts3D);
end

w_guess = rodrigues(R_guess);

% use multicamera wrapper considering only one camera for computing
% error and jacobian functions
[e, j] = ej_wrapper_reprojection(K, eye(3), zeros(3, 1), pts3D', pts2D');

o = @(x) objective(x, e, j);

% solve
options=optimset('Display','off',...
                 'DerivativeCheck','off',...
                 'Jacobian','on',...
                 'MaxIter', LM_MAX_ITERATIONS, ...
                 'MaxFunEvals', inf, ...
                 'TolX', LM_TOLX, ...
                 'TolFun', LM_TOLFUN, ...
                 'Algorithm', 'levenberg-marquardt'); 

sol = lsqnonlin(o, double([w_guess', t_guess']), [], [], options);

R = rodrigues(sol(1:3));
t = sol(4:6)';
end

function [e, j] = objective(p, err_fun, jacob_fun)
    if nargout < 2
        e = double(err_fun(p));
    else
        e = double(err_fun(p));
        j = double(jacob_fun(p));
    end
end
