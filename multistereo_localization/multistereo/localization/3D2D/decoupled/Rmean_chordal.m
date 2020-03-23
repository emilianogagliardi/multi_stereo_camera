function R = Rmean_chordal(Rs, weights)
%RMEAN_CHORDAL 

if nargin == 1
    weights = ones(1, size(Rs, 3));
end

weights = weights ./ sum(weights);

Ce = zeros(3);
n = size(Rs, 3);
for ii = 1:n
    if all(all(~isnan(Rs(:, :, ii))))
        Ce = Ce + weights(ii) * Rs(:, :, ii);
    end
end
[U, ~, V] = svd(Ce);
UVt = U * V';
if det(UVt) >= 0
    R = UVt;
else
    R = U * diag([1, 1, -1]) * V';
end
end

