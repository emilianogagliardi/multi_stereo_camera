function [projected, err] = reprojection(K, R, t, p3D, p2D)
if size(p3D, 2) == 3
    p3D = p3D';
    p2D = p2D';
end
projected = K * [R, t] * [p3D; ones(1, size(p3D, 2))];
projected = projected(1:2, :) ./ projected(3, :);
err = vecnorm(p2D-projected, 2, 1);
projected = projected';
end

