function [p3D, desc, cam_idx, kL_out, depth] = reconstruct_ms(kLs, kRs, dLs, dRs, iLs, iRs, multistereo)
%RECONSTRUCT_MS 
% Parameters
% - kLs, kRs: 3xN and 3xM matrices keypoints extracted in left and right
% images
% - dLs, dRs: corresponding descriptors NxX and MxX where X is descriptor
% dimensionality
% - iLs, iRs: indices of the stereo pair perceiving the point
% - multistereo: MultiStereo object
% Returns:
% - p3D: reconstructed point in principal camera frame 4xN (homog coords)
% - desc: descriptors from left cameras of corresponding reconstructed
% points
% - cam_idx: indices of the camera used for reconstructing corresponding
% point, N elements vector
% - kL_out: keypoints on left cameras used for reconstruction, can be used
% to visualize matching with other views

global REPROJ_ERROR_OUTLIERS;

depth_required = nargout > 4;

% iterate over the cameras and perform reconstruction
for ii = multistereo.stereoPairsNumber():-1:1
    
    % get keypoints and descriptors of current camera
    curr_camL = iLs == ii;
    curr_camR = iRs == ii;
    kL = kLs(:, curr_camL);
    kR = kRs(:, curr_camR);
    dL = dLs(curr_camL, :);
    dR = dRs(curr_camR, :);
    
    % matching
    [kLm, kRm, dLm] = stereo_matching(dL, dR, kL, kR);
    
    % perform reconstruction by linear triangulation of matched features
    if depth_required
        [Rs, ts] = multistereo.stereoLeft2Right(ii);
        [Kl, Kr] = multistereo.getKs(ii);
        Pl = [Kl, zeros(3, 1)];
        Pr = Kr * [Rs, ts];
    else
        [Pl, Pr] = multistereo.projMatrix(ii);
    end
    
    if ~REPROJ_ERROR_OUTLIERS
        p3Dii = triangulate(from_homogeneous(kLm)', from_homogeneous(kRm)', Pl', Pr')';
        
    else % reprojection error trimming
        [p3Dii, errs] = triangulate(from_homogeneous(kLm)', from_homogeneous(kRm)', Pl', Pr');
        valid_idx = errs <= mean(errs);
        p3Dii = p3Dii(valid_idx, :)';
        dLm = dLm(valid_idx, :);
        kLm = kLm(:, valid_idx);
    end
    
    if depth_required
        out(ii).depth = p3Dii(3, :);
        % points have been computed in stereo frame, bring to rig frame
        [Rc, tc] = multistereo.Tcam2principal(ii);
        p3Dii = Rc * p3Dii + tc;
    end
    out(ii).p3D = p3Dii;
    out(ii).dL = dLm;
    out(ii).kL = kLm;
    out(ii).cam_idx = ii * ones(size(p3Dii, 2), 1);
end

p3D = cat(2, out.p3D);
desc = cat(1, out.dL);
cam_idx = cat(1, out.cam_idx);
if nargout > 3
    kL_out = cat(2, out.kL);
end
if depth_required
    depth = cat(2, out.depth);
end
end

