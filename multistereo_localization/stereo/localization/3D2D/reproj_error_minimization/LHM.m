function [R, t, it, new_err] = LHM(P, V, R)
% LHM iterative reprojection error refinement
% Parameters
% - P: 3xN matrix containing object points
% - V: normalized image coordinates, 3xN matrix [u; v; 1]
% - R: rotation guess
% Returns
% - R, t: refined extrinsics
% - it: number of iterations
% - new_err: object space collinearity error


    global LHM_EPSILON; % error threshold to stop
    global LHM_TOLERANCE; % stop if abs(old_error-new_err)/old_err < LHM_TOLERANCE
    global LHM_MAX_ITERATIONS;
    
    % compute projection matrices
    F = F_(V);

    % compute first element in t(R) expression
    t_factor = inv(eye(3) - mean(F, 3))/size(P, 2);
    
    % compute t from R guess
    t = t_(P, F, R, t_factor);
    
    % compute projection on viewing rays of R*p+t
    Q = Q_(P, F, R, t);
    
    % error of guess rotation
    old_err = mean(vecnorm(R*P+t - Q, 2, 1).^2);
    
    % one iteration
    [R, t] = kernel(P, Q, F, t_factor);
    Q = Q_(P, F, R, t);
    new_err = mean(vecnorm(R*P+t - Q, 2, 1).^2);

    % iterate
    it = 1;
    while new_err > LHM_EPSILON && abs(old_err - new_err) / old_err > LHM_TOLERANCE
        old_err = new_err;
        [R, t] = kernel(P, Q, F, t_factor);
        Q = Q_(P, F, R, t);
        new_err = mean(vecnorm(R*P+t - Q, 2, 1).^2);
        it = it + 1;
        if it > LHM_MAX_ITERATIONS 
            break
        end
    end

end

function F = F_(V)
    n = size(V, 2);
    F = zeros(3, 3, n);
    for ii = 1:n
        V_i = V(:, ii);
        F(:,:,ii) = (V_i*V_i')/(V_i'*V_i);
    end
end

function Q = Q_(P, F, R, t)
    n = size(P, 2);
    Q = zeros(3, n);
    for ii = 1:n
        Q(:, ii) = F(:, :, ii) * (R * P(:, ii) + t);
    end
end

function t = t_(P, F, R, t_factor)
    n = size(P, 2);
    sum_ = zeros(3, 1);
    for ii = 1:n
        sum_ = sum_ + (F(:,:,ii)-eye(3)) * R * P(:, ii);
    end
    t = t_factor * sum_;
end

function [R, t] = kernel(P, Q, F, t_factor)
    R = absor(P, Q);
    t = t_(P, F, R, t_factor);
end
