function f = visualize_ms_reconstruction(p3D, cam_idx, multistereo)
f = figure();
hold on;
for ii = 1:multistereo.stereoPairsNumber()
    idx = cam_idx == ii;
    p = p3D(:, idx);
    plot3(p(1, :), p(2, :), p(3, :), 'o', 'Color', color_map(ii), 'MarkerSize', 1);
end
axis equal
axis tight
xlabel('X');
ylabel('Y');
zlabel('Z');
end

