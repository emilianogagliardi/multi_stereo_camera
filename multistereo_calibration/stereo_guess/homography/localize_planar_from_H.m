function [R, t] = localize_planar_from_H(K, H)
%LOCALIZE_PLANAR_FROM_H compute camera pose from homography
ijo = K \ H;
% normalize the i and j vectors
n = mean([norm(ijo(:, 1), 2), norm(ijo(:, 2), 2)]);
i = ijo(:, 1) / norm(ijo(:, 1), 2);
j = ijo(:, 2) / norm(ijo(:, 2), 2);
k = cross(i, j) / norm(cross(i, j), 2);
t = ijo(:, 3)/n;
R = [i, j, k];
% due to numerical error R might not be a rotation matrix, approximate it
[U, ~, V] = svd(R);
R = (U * V');
end

