function R = Rmean_geodesic(Rs, weights, epsilon)
%RMEAN_GEODESIC compute rotations mean according to SO(3) geodesic distance
R = Rs(:, :, 1);
n = size(Rs, 3);
while 1
    r = zeros(3, 1);
    for ii = 1:n
        r = r + weights(ii) * rodrigues(R' * Rs(:, :, ii));
    end
    r = 1/n * r;
    if norm(r) < epsilon
        return
    end
    R = R * rodrigues(r);
end
end