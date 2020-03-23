function visualize_ms_reprojection(ims, p3D_in, p2D, cam_idx, R, t, multistereo, R_guess, t_guess)
%VISUALIZE_MS_REPROJECTION
fs = [];

if size(p2D, 2) == 3
    p2D = from_homogeneous(p2D);
end

p3D = R * p3D_in + t;

if nargin > 7
    p3D_guess = R_guess * p3D_in + t_guess;
end

for ii = 1:size(ims, 3)
    idx = cam_idx == ii;
    p2Dii = p2D(:, idx);
    p3Dii = p3D(:, idx);
    P = multistereo.projMatrix(ii);
    [reproj, errs] = reproject(P, p3Dii, p2Dii);
    fs = [fs, figure()];
    imshow(ims(:, :, ii));
    hold on
    plot(p2Dii(1, :), p2Dii(2, :), 'g+');
    plot(reproj(1, :), reproj(2, :), 'r+');
    if nargin > 7
        p3Dii_guess = p3D_guess(:, idx);
        [reproj_, errs_] = reproject(P, p3Dii_guess, p2Dii);
        plot(reproj_(1, :), reproj_(2, :), 'b+');
        disp(['Image ', num2str(ii), ' guess mean reprojection error ' num2str(mean(errs_))]);
        title('GREEN: detected, BLUE: reprojected guess, RED: reprojected optimized');
    else
        title('GREEN: detected, RED: reprojected');
    end
    disp(['Image ', num2str(ii), ' mean reprojection error ' num2str(mean(errs))]);
    pause
end
close(fs)
end

