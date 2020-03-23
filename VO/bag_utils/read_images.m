function [imagesL, imagesR, ts] = read_images(topicsL, idx, cams, topicsR)
% read stereo images from bag, images are assumed to have all the same
% size, topicsR is optional 
n = length(cams);
for ii = n:-1:1
    msg = readMessages(topicsL(ii), idx);
    im = readImage(msg{1});
    im = undistortImage(im, cams{ii}.CameraParameters1);
    imagesL(:, :, ii) = im;
end
if nargout > 1
   for ii = n:-1:1
        msg = readMessages(topicsR(ii), idx);
        im = readImage(msg{1});
        im = undistortImage(im, cams{ii}.CameraParameters1);
        imagesR(:, :, ii) = im;
    end 
end
ts = msg{1}.Header.Stamp;
end

