function m = track_keypoints_ms(p3D, d3D, p2D, d2D, p2D_cam_idx, ...
    multistereo, R_guess, t_guess, im_size, norm_fun, max_dist)

global TRACK_SEARCH_SIZE;

n_cameras = multistereo.stereoPairsNumber();

m = zeros(min(size(p3D, 2), size(p2D, 2)), 2);
match_idx = 0;

if size(p2D, 1) == 3
    p2D = from_homogeneous(p2D);
end

if size(p3D, 1) == 4
    p3D = from_homogeneous(p3D);
end
p3D = R_guess * p3D + t_guess;

% build a tensor with [im_height, im_width, n_camera] size of all zeros,
% then for each p2D, if round(p2D(:, ii)) = [col, row] and p2D comes from camera
% jj p2D_tensor(row, col, jj) = ii
p2D_tensor = zeros([im_size, n_cameras]);
for ii = 1:size(p2D, 2)
    col = ceil(p2D(1, ii));
    row = ceil(p2D(2, ii));
    cam = p2D_cam_idx(ii);
    p2D_tensor(row, col, cam) = ii;
end


% iterate over 2D points and project them on all cameras, then look for
% matching in all images in a squared region
for p3D_idx = 1:size(p3D, 2)
    curr_p3D = p3D(:, p3D_idx);
    curr_p3D_desc = d3D(p3D_idx, :);
    
    for cam_idx = 1:n_cameras
        % project the point on the current camera according to the guessed
        % pose
        K = multistereo.getKs(cam_idx);
        [R_, t_] = multistereo.Tprincipal2cam(cam_idx);
        curr_p3D_ = R_ * curr_p3D + t_;
        
        if curr_p3D_(3) <= 0
            continue;
        end
        
        guess = from_homogeneous(reproject([K, zeros(3, 1)], curr_p3D_));
        
        % compute the corners of a patch around the guessed point
        top_left_row = ceil(guess(2) - TRACK_SEARCH_SIZE /2);
        top_left_col = ceil(guess(1) - TRACK_SEARCH_SIZE /2);
        bottom_right_row = ceil(guess(2) + TRACK_SEARCH_SIZE /2);
        bottom_right_col = ceil(guess(1) + TRACK_SEARCH_SIZE /2);
        
        % check if it is in field of view
        if top_left_row < 1 || top_left_col < 1 || ...
            bottom_right_row > im_size(1) || bottom_right_col > im_size(2)
            continue;
        end
        
        % retrieve all keypoints in the predicted region
        candidates_idx = p2D_tensor(top_left_row:bottom_right_row, ...
            top_left_col:bottom_right_col, cam_idx);
        candidates_idx = candidates_idx(:);
        % if candidate idx is 0 no keypoint falls in the corresponding bin
        candidates_idx = candidates_idx(candidates_idx ~= 0); 
        
%         % esc if no candidates
        if isempty(candidates_idx) 
            continue;
        end
        
        best_dist = max_dist;
        best_idx = 0;
        for ii = 1:length(candidates_idx)
            curr_p2D_idx = candidates_idx(ii);
            curr_p2D_desc = d2D(curr_p2D_idx, :);
            dist = norm_fun(curr_p3D_desc, curr_p2D_desc);
            if dist < best_dist
                best_dist = dist;
                best_idx = curr_p2D_idx;
            end
        end
        
        % if something found add match
        if best_idx ~= 0
            m(match_idx + 1, :) = [p3D_idx, best_idx];
            match_idx = match_idx + 1;
        end
    end
end

% if some match found
if match_idx > 0
    m = m(1:match_idx, :);
else
    m = [];
end