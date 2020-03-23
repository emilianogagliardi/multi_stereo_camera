clear
close all
clc

VISUALIZATION = true;

% load calibration and create multistereo object
load('examples/data/calibration/stereos.mat'); 
stereos = {cam1, cam2, cam3, cam4}; % no bottom camera
load('examples/data/calibration/multicams.mat');
multicalib = sol; clear sol;
multistereo = MultiStereo(stereos, multicalib);

% load localization parameters
parameters;

% load the bag containing odometry data
bag = rosbag(['/home/emiliano/Development/bags/guidance' ...
    '/handeye_odom/1/odom/odom1.bag']);
topic1L = select(bag, 'Topic', '/guidance/1/left');
topic2L = select(bag, 'Topic', '/guidance/2/left');
topic3L = select(bag, 'Topic', '/guidance/3/left');
topic4L = select(bag, 'Topic', '/guidance/4/left');
topicsL = [topic1L, topic2L, topic3L, topic4L];
topic1R = select(bag, 'Topic', '/guidance/1/right');
topic2R = select(bag, 'Topic', '/guidance/2/right');
topic3R = select(bag, 'Topic', '/guidance/3/right');
topic4R = select(bag, 'Topic', '/guidance/4/right');
topicsR = [topic1R, topic2R, topic3R, topic4R];

first_frame = 1;

% .......................FIRST FRAME---------------------------------------
% perform reconstruction in the first frame
[imsL_old, imsR_old, image_time] = read_images(topicsL, first_frame, stereos, topicsR);
[kL, dL, camidxL, kR, dR, camidxR] = detect_keypoints_ms(imsL_old, imsR_old);
[p3D_old, desc_old, cam_idx_old, kL_old] = reconstruct_ms(kL, kR, dL, ...
    dR, camidxL, camidxR, multistereo);

% assing a unique id to each 3D point
landmarks_ID = 1:size(p3D_old, 2);
% starting point for the creation of new landmark indices
first_new_index = size(p3D_old) + 1;

% add the 3D points to the map, with their corresponding perceptions
untracked = struct('p3D', p3D_old, 'kpts', from_homogeneous(kL_old), 'cam_idx', cam_idx_old, ...
    'landmarks_ID', landmarks_ID);
% no tracked points at first frame
tracked = struct('kpts', [], 'cam_idx', [], 'landmarks_id', []);
map = MapObject(multistereo);
map = map.AddFrame(tracked, untracked, eye(3), zeros(3, 1));


images_time = image_time;

N = topic1L.NumMessages; % loop over the bag
for timestamp = first_frame+1:N-1
    [imsL, imsR, image_time] = read_images(topicsL, timestamp, stereos, topicsR);
    %[~, Rcheck, tcheck, ~, imsL] = multicamera_checkerboard_gt(imsL, cams, multistereo, 108, [7, 9], VISUALIZATION);
    %T_check = cat(3, T_check, [Rcheck, tcheck; zeros(1, 3), 1]);
    images_time = [images_time, image_time];
    
    %--------------------EXTRACTION AND RECONSTRUCTION---------------------
    % extract keypoints and perform reconstruction, for simplicity only
    % left image points surviving to stereo matching are kept for subsequent
    % stepts, and right image points are discarded
    [kL, dL, camidxL, kR, dR, camidxR] = detect_keypoints_ms(imsL, imsR);
    [p3D, desc, cam_idx, kL] = reconstruct_ms(kL, kR, dL, dR, camidxL, ...
        camidxR, multistereo);
    
    %----------------------------LOCALIZATION------------------------------
    % localize the camera wrt previous reconstruction (3D2D)
    % first time no constant motion model can be applied
    %if timestamp == first_frame + 1
        m = matchFeatures(desc_old, desc, 'method', 'approximate', ...
            'MatchThreshold', 1, 'MaxRatio', 0.6);
    %else % apply constant motion model for matching
    %    m = track_keypoints_ms(p3D_old, desc_old, kL, desc, ...
    %        cam_idx, multistereo, R, t, [size(imsL, 1), size(imsL, 2)], ...
    %        @(x, y) norm(x-y), 1);
%          m = matchFeatures(desc_old, desc, 'method', 'exhaustive', ...
%              'MatchThreshold', 1.2, 'MaxRatio', 1);
    %end
    p3D_old_m = p3D_old(:, m(:, 1));
    kL_m = kL(:, m(:, 2));
    cam_idx_m = cam_idx(m(:, 2));
    [R, t, inliers] = localize_ms_coupled_3D2D(p3D_old_m, kL_m, ...
        cam_idx_m, multistereo);
    
    %------------------------UPDATE THE MAP--------------------------------
    % inliers image points need to be added as perception of an already
    % existing landmark (3D point) in the map. Outliers or non matched image
    % points and correspondind reconstruction are instead new points,
    % (new landmarks with corresponding perceptions) that need a new unique 
    % id to be created
    
    % outliers and non matched
    non_matched_idx = setdiff(1:size(kL, 2), m(:, 2));
    matched_outliers_idx = m(~inliers, 2)';
    % indices of the current p3D, kL, descriptors that have not been
    % matched with the previous perception, or that have been matched but
    % resulted as outliers in the localization step
    not_tracked_idx = [non_matched_idx, matched_outliers_idx];
    non_tracked_p3D = p3D(:, not_tracked_idx);
    non_tracked_kL = kL(:, not_tracked_idx);
    non_tracked_cam_idx = cam_idx(not_tracked_idx);
    non_tracked_desc = desc(not_tracked_idx, :);
    new_landmarks_ID = ...
        first_new_index : first_new_index + size(non_tracked_p3D, 2)-1;
    first_new_index = first_new_index + size(non_tracked_p3D, 2);
    % create the untracked points structure
    untracked = struct('p3D', non_tracked_p3D, 'kpts', from_homogeneous(non_tracked_kL), ...
        'cam_idx', non_tracked_cam_idx, 'landmarks_ID', new_landmarks_ID);
    
    % inliers
    % indices of the current p3D, kL, descriptors that resulted as inliers
    % in the localization step
    tracked_idx = m(inliers, 2);
    tracked_p3D = p3D(:, tracked_idx);
    tracked_kL = kL(:, tracked_idx);
    tracked_cam_idx = cam_idx(tracked_idx);
    tracked_desc = desc(tracked_idx, :);
    % id of the tracked landmarks
    tracked_landmarks_ID = landmarks_ID(m(inliers, 1));
    % create the tracked points structure
    tracked = struct('kpts', from_homogeneous(tracked_kL), 'cam_idx', tracked_cam_idx, ...
        'landmarks_ID', tracked_landmarks_ID);
    
    map = map.AddFrame(tracked, untracked, R, t);
    map = map.Optimize(3);
    
    %------------------------------VISUALIZATION---------------------------
    if VISUALIZATION
        % visualize matches
        kL_old_m = kL_old(:, m(:, 1));
        cam_idx_old_m = cam_idx_old(m(:, 1));
        % show all 
%         visualize_ms_matches(imsL, imsL_old, kL_m, ...
%             kL_old_m, cam_idx_m, cam_idx_old_m, 'falsecolor');
        % show inliers
        visualize_ms_matches(imsL, imsL_old, kL_m(:, inliers), ...
            kL_old_m(:, inliers), cam_idx_m(inliers), ...
            cam_idx_old_m(inliers), 'falsecolor');
        pause 
        close all
        % visualize reprojections
        visualize_ms_reprojection(imsL, p3D_old_m(:, inliers), ...
            kL_m(:, inliers), cam_idx_m(inliers), R, t, multistereo);
    end
    
    %---------------UPDATE PERCEPTION AT t-1 INSTANT-----------------------
    % current reconstruction and relative image point and descriptors
    % become old perception in the next loop iteration
    imsL_old = imsL;
    p3D_old = [tracked_p3D, non_tracked_p3D];
    kL_old = [tracked_kL, non_tracked_kL];
    cam_idx_old = [tracked_cam_idx; non_tracked_cam_idx];
    desc_old = [tracked_desc; non_tracked_desc];
    landmarks_ID = [tracked_landmarks_ID, new_landmarks_ID];
end

%% 
% get optitrack poses and timestamps
optis_time = [];
opti_poses = zeros(length(images_time), 7); % tx ty tz qw qx qy qz
for ii = 1:length(images_time)
    msg = nearest_msg(bag, images_time(ii), '/Robot_1/pose', 0.1);
    optis_time = [optis_time, msg.Header.Stamp];
    t_opti = [msg.Pose.Position.X, msg.Pose.Position.Y, ...
        msg.Pose.Position.Z];
    q_opti = [msg.Pose.Orientation.X, ...
              msg.Pose.Orientation.Y, ...
              msg.Pose.Orientation.Z, ...
              msg.Pose.Orientation.W];
    opti_poses(ii, :) = [t_opti*1000, q_opti];
end

% write images poses in tx ty tz qx qy qz format
[Rs, ts] = map.GetPoses();
image_poses = zeros(length(images_time), 7); % tx ty tz qw qx qy qz
for ii = 1:size(Rs, 3)
    q = rotm2quat(Rs(:, :, ii));
    q = [q(2:4), q(1)];
    ts(:, ii);
    image_poses(ii, :) = [ts(:, ii)', q];
end

% write trajectories on file as 
% timestamp tx ty tz qx qy qz qw
f_camera = fopen('stamped_traj_estimate.txt', 'wt');
f_ground = fopen('stamped_groundtruth.txt', 'wt');
for ii = 1:length(images_time)
  t_camera = rostime2sec(images_time(ii));
  t_opti = rostime2sec(optis_time(ii));
  fprintf(f_camera, '%.20d %.20d %.20d %.20d %.20d %.20d %.20d %.20d\n', [t_camera, image_poses(ii, :)]);
  fprintf(f_ground, '%.20d %.20d %.20d %.20d %.20d %.20d %.20d %.20d\n', [t_opti, opti_poses(ii, :)]);
end
fclose(f_camera);
fclose(f_ground);


