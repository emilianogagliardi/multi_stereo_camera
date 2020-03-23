function [valid, R, t, err, im_out] = checkerboard_gt(im, cam, square_size, checkerboard_size, show_check_points)
% CHECKERBOARD_GT localize wrt checkerboard and remove the checkerboard from the image 
% valid = true if checkerboard detected, R, t are extrinsics in
% checkerboard frame, err is reprojection error, im_out is the input image
% in which the checkerboard has been covered


[imagePoints,checkSize] = detectCheckerboardPoints(im);
if checkSize ~= checkerboard_size
    valid = false;
    R = nan(3); t = nan(3, 1); err = nan(); im_out = im;
    return;
end
worldPoints = generateCheckerboardPoints(checkerboard_size, square_size);
if size(imagePoints, 1) ~= size(worldPoints, 1)
    valid = false;
    R = nan(3); t = nan(3, 1); err = nan(); im_out = im;
    return;
end
[R, t] = extrinsics(imagePoints,worldPoints,cam);
R = R';
t = t';
K = cam.IntrinsicMatrix';
P = K * [R, t];
reproj = P * [worldPoints, zeros(size(worldPoints, 1), 1), ones(size(worldPoints, 1), 1)]';
reproj = reproj(1:2, :) ./ reproj(3, :);
err = mean(vecnorm(imagePoints - reproj', 2, 2));
valid = true;

% remove the checkerboard from the image
pts = [-2*square_size, -2*square_size, 0; 
    -2*square_size, checkerboard_size(1) * square_size, 0;
    checkerboard_size(2) * square_size, checkerboard_size(1) * square_size, 0;
    checkerboard_size(2) * square_size, -2*square_size, 0];
im_pts = P * [pts'; ones(1, 4)];
im_pts = im_pts(1:2, :) ./ im_pts(3, :);
im_out = insertShape(im,'FilledPolygon',[im_pts(1, 1), im_pts(2, 1), im_pts(1, 2), im_pts(2, 2), im_pts(1, 3), im_pts(2, 3), im_pts(1, 4), im_pts(2, 4)], 'SmoothEdges', true, 'Color', [0, 0, 0], 'Opacity', 1);
im_out = rgb2gray(im_out);

if show_check_points
    f = figure();
    imshow(im);
    hold on
    plot(imagePoints(:, 1), imagePoints(:, 2), 'go');
    plot(reproj(1, :), reproj(2, :), 'r+');
    pause;
    close(f);
end
end





