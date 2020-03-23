function [R, t, error, weighted_error] = absor( ...
    pts1, pts2, weights1, weights2)
% ABSOLUTE-ORIENTATION: implementation of "Closed form solution of absolute orientation using unit quaternions" from Berthold K. P. Horn
    % Here no scale factor is considered.
    % weighted centroids are computed then the product of the weight of the
    % pair is kept for computing rotation
    % Parameters:
    % - pts1, pts2: 3xN matrices, where columns are 3D vectors
    % representing the same points in different coordinates system.
    % pts1(:, ii) needs to represent the same point as pts2(:, ii)
    % - weights: 1xN vector, weighting the error function for each pair
    % Returns:
    % The rototralsation that applied to pts1 minimizes:
    % sum_ii{ weights(ii) * norm(R * pts1(:, ii) + t - pts2(:, ii))^2 }
    
    if nargin == 2
        weights1 = ones(1, size(pts1, 2));
        weights2 = ones(1, size(pts1, 2));
    end
    
    % normalize weights
    weights1 = weights1 / sum(weights1);
    weights2 = weights2 / sum(weights2);
    
    % weighted centroids
    c1 = mean(pts1, 2);
    c2 = mean(pts2, 2);
    
    % ROTATION
    pts1_ = pts1 - c1;
    pts2_ = pts2 - c2;
    
    weights = weights1 .* weights2;
    weights = weights/sum(weights);
    
    M = pts1_ * pts2_';
    % compute the matrix N
    diagonal = diag( [M(1, 1) + M(2, 2) + M(3, 3), ...
        M(1, 1) - M(2, 2) - M(3, 3), ...
        -M(1, 1) + M(2, 2) - M(3, 3), ...
        -M(1, 1) - M(2, 2) + M(3, 3)]);
    up_left = ...
        [0, M(2, 3) - M(3, 2), M(3, 1) - M(1, 3), M(1, 2) - M(2, 1);
         0,                 0, M(1, 2) + M(2, 1), M(1, 3) + M(3, 1);
         0,                 0,                 0, M(2, 3) + M(3, 2);
         0,                 0,                 0,                 0];
    N = diagonal + up_left + up_left';
    % eigenvector corresponding to the most positive eigenvalue of N
	[V,D]=eig(N);
    [~,emax]=max(real(diag(D))); 
    emax=emax(1);
    q = V(:,emax); % Gets eigenvector corresponding to maximum eigenvalue
    q = real(q);   % Get rid of imaginary part caused by numerical error
    [~, ii]=max(abs(q)); 
    sgn=sign(q(ii(1)));
    q=q*sgn; % Sign ambiguity
    R = quat2rotm(q');
    
    % TRANSLATION
    t = c2 - R * c1;
    
    % compute the error
    error_vecs = R * pts1 + t - pts2;
    error_norms = vecnorm(error_vecs, 2, 1);
    error_norms_square = error_norms .^ 2;
    error = mean(error_norms_square);
    
    weighted_error = mean(weights1 .* error_norms_square);
end

