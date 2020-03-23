function f = visualize_reprojection(im, p3D, p2D, K, R, t, R_guess, t_guess)
%VISUALIZE_REPROJECTION 
% R_guess and t_guess are optional

P = K * [R, t];
[repr, err] = reprojection(P, p3D, p2D);
figure, imshow(im);

if size(p2D, 1) == 3
    p2D = from_homogeneous(p2D);
end
hold on
plot(p2D(1, :), p2D(2, :), 'g+');
plot(repr(1, :), repr(2, :), 'r+');

if nargin > 6
    P_guess = K * [R_guess, t_guess];
    [repr_guess, err_guess] = reproject(P_guess, p3D, p2D);
    plot(repr_guess(1, :), repr_guess(2, :), 'b+');
    disp(['Guess mean reprojection error ' num2str(mean(err_guess))]);
    title('GREEN: detected, BLUE: reprojected guess, RED: reprojected optimized');
else
    title('GREEN: detected, RED: reprojected');
end
disp(['Mean reprojection error ' num2str(mean(err))]);


end

