function [p3D, desc, kLm] = reconstruct(kL, kR, dL, dR, stereo)
%RECONSTRUCT
% Parameters
% - kL, kR: 3xN and 3xM matrices keypoints extracted in left and right
% images
% - dLs, dRs: corresponding descriptors NxX and MxX where X is descriptor
% dimensionality
% - stereo: stereoParameters object (matlab stereo camera calibration)
% Returns:
% - p3D: reconstructed point in left camera frame 4xN (homog coords)
% - desc: descriptors from left cameras of corresponding reconstructed
% points
% - kLm: keypoints on left cameras used for reconstruction, can be used
% to visualize matching with other views

global REPROJ_ERROR_OUTLIERS;

% matching
[kLm, kRm, desc] = stereo_matching(dL, dR, kL, kR);

% perform reconstruction by linear triangulation (SVD) of matched features
[Pl, Pr] = stereo_Ps(stereo);

if ~REPROJ_ERROR_OUTLIERS
    p3D = triangulate(from_homogeneous(kLm)', from_homogeneous(kRm)', Pl', Pr')';
else % reprojection error trimming
    [p3D, errs] = triangulate(from_homogeneous(kLm)', from_homogeneous(kRm)', Pl', Pr');
    valid_idx = errs <= mean(errs);
    p3D = p3D(valid_idx, :)';
    desc = desc(valid_idx, :);
    kLm = kLm(:, valid_idx);
end
end

