function [found, Rout, tout, err, images_out] = multicamera_checkerboard_gt(images, cams, multistereo, square_size, checkerboard_size, show_check_points)
% find the checkerboard in the images and give the position of the
% multicamera with respect to it in the principal camera frame
% cams is an array of cameraParameters, images a NxMxNum_cams (concatenated
% images, multistereo is a MultiStereo object

images_out = [];
found = false;
for ii = 1:multistereo.stereoPairsNumber
    [valid, R, t, err, im_out] = checkerboard_gt(images(:, :, ii), cams{ii}.CameraParameters1, square_size, checkerboard_size, show_check_points);
    if valid 
        found = true;
        T_ = [R, t; zeros(1, 3), 1];
        [Rcam, tcam] = multistereo.Tcam2principal(ii);
        Tcam = [Rcam, tcam; zeros(1, 3), 1];
        T = Tcam * T_;
        Rout = T(1:3, 1:3); tout = T(1:3, 4);
    end
    images_out = cat(3, images_out, im_out);
end

if ~found
    Rout = nan(3);
    tout = nan(3, 1);
    err = nan;
    images_out = images;
end

end

