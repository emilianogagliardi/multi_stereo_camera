clear
close all
clc

VISUALIZATION = false;

% load calibration and create multistereo object
load('examples/data/calibration/stereos.mat'); 
stereos = {cam1, cam2, cam3, cam4}; % no bottom camera
load('examples/data/calibration/multicams5.mat');
multicalib = sol; clear sol;
multistereo = MultiStereo(stereos, multicalib);

% load localization parameters
parameters;

e_3D2D = [];
e_dec3D2D = [];
e_3D3D = [];
e_hybrid = [];
e_single = [];

for ii = 1:10
    bag1 = rosbag('/home/emiliano/Development/bags/guidance/evaluation/eval_checkerboard7/pose1.bag');
    topic1L = select(bag1, 'Topic', '/guidance/1/left');
    topic2L = select(bag1, 'Topic', '/guidance/2/left');
    topic3L = select(bag1, 'Topic', '/guidance/3/left');
    topic4L = select(bag1, 'Topic', '/guidance/4/left');
    topicsL1 = [topic1L, topic2L, topic3L, topic4L];
    topic1R = select(bag1, 'Topic', '/guidance/1/right');
    topic2R = select(bag1, 'Topic', '/guidance/2/right');
    topic3R = select(bag1, 'Topic', '/guidance/3/right');
    topic4R = select(bag1, 'Topic', '/guidance/4/right');
    topicsR1 = [topic1R, topic2R, topic3R, topic4R];
    bag2 = rosbag(['/home/emiliano/Development/bags/guidance/evaluation/eval_checkerboard7/pose', num2str(ii+1), '.bag']);
    topic1L = select(bag2, 'Topic', '/guidance/1/left');
    topic2L = select(bag2, 'Topic', '/guidance/2/left');
    topic3L = select(bag2, 'Topic', '/guidance/3/left');
    topic4L = select(bag2, 'Topic', '/guidance/4/left');
    topicsL2 = [topic1L, topic2L, topic3L, topic4L];
    topic1R = select(bag2, 'Topic', '/guidance/1/right');
    topic2R = select(bag2, 'Topic', '/guidance/2/right');
    topic3R = select(bag2, 'Topic', '/guidance/3/right');
    topic4R = select(bag2, 'Topic', '/guidance/4/right');
    topicsR2 = [topic1R, topic2R, topic3R, topic4R];
    
    % perform reconstruction in bag 1
    [imsL1, imsR1] = read_images(topicsL1, 1, stereos, topicsR1);
    [~, R_gt1, t_gt1, ~, imsL1] = multicamera_checkerboard_gt(imsL1, stereos, multistereo, 108, [7, 9], VISUALIZATION);
    T_gt1 = [R_gt1, t_gt1; zeros(1, 3), 1];
    [kL1, dL1, camidxL1, kR1, dR1, camidxR1] = detect_keypoints_ms(imsL1, imsR1);
    [p3D1, dL1, cam_idx1, kL1, depths1] = reconstruct_ms(kL1, kR1, dL1, dR1, camidxL1, camidxR1, multistereo);
    
    
    e_angles1 = [];
    e_ts1 = [];
    inl_num1 = [];
    e_angles2 = [];
    e_ts2 = [];
    inl_num2 = [];
    e_angles3 = [];
    e_ts3 = [];
    inl_num3 = [];
    e_angles4 = [];
    e_ts4 = [];
    inl_num4 = [];
    e_angles5 = [];
    e_ts5 = [];
    inl_num5 = [];
    
    for jj = 1:10
        % localization in bag 2 wrt bag 1 reconstruction
        imsL2 = read_images(topicsL2, 1, stereos);
        imsR2 = read_images(topicsR2, 1, stereos);
        [~, R_gt, t_gt, ~, imsL2] = multicamera_checkerboard_gt(imsL2, stereos, multistereo, 108, [7, 9], VISUALIZATION);
        T_gt = [R_gt, t_gt; zeros(1, 3), 1] / T_gt1;
        R_gt = T_gt(1:3, 1:3);
        t_gt = T_gt(1:3, 4);
        [kL2, dL2, camidxL2, kR2, dR2, camidxR2] = detect_keypoints_ms(imsL2, imsR2);
        
        [p3D2, dL2, camidxL2, kL2] = reconstruct_ms(kL2, kR2, dL2, dR2, camidxL2, camidxR2, multistereo);
        
        m = matchFeatures(dL1, dL2, 'method', 'exhaustive', 'MatchThreshold', 1.5, 'MaxRatio', 0.8);
        p3D1_m = p3D1(:, m(:, 1));
        kL1_m = kL1(:, m(:, 1));
        camidxL1_m = cam_idx1(m(:, 1));
        kL2_m = kL2(:, m(:, 2));
        camidxL2_m = camidxL2(m(:, 2));
        
        [R1, t1, inliers1, ~, ~, errs] = localize_ms_coupled_3D2D(p3D1_m, kL2_m, camidxL2_m, multistereo);
        [R2, t2, inliers2] = localize_ms_decoupled_3D2D(p3D1_m, kL2_m, camidxL2_m, multistereo, 6); 
        p3D2_m = p3D2(:, m(:, 2));
        [R3, t3, inliers3] = absor_ransac(p3D1_m', p3D2_m');
        [R4, t4, inliers4] = absor_ransac_ms_2Derr(p3D1_m', p3D2_m', from_homogeneous(kL2_m)', camidxL2_m, multistereo, [], [], true);
        
        R = R4; t = t4; inliers = inliers4;
        
        
        e = inf;
        for cam = 1:4
            kL2_m_cam = kL2_m(:, camidxL2_m == cam);
            p3D1_m_cam = p3D1_m(:, camidxL2_m == cam);
            K = multistereo.getKs(cam);
            [R_, t_, inliers_] = localize_3D2D(p3D1_m_cam, kL2_m_cam, K);
            [R__, t__] = multistereo.Tcam2principal(cam);
            R_ = R__ * R_;
            t_ = R__ * t_ + t__;
            thiserr = norm(t_-t_gt);
            if thiserr < e
                e = thiserr;
                R5 = R_;
                t5 = t_;
                inliers5 = inliers_;
            end
        end
        
        angle_gt = rad2deg(rotm2eul(R_gt))';

        e_ts1 = [e_ts1, norm(t1-t_gt)];
        e_angles1 = [e_angles1, rad2deg(rotm2eul(R1))' - angle_gt];
        inl_num1 = [inl_num1, sum(inliers1)];
        e_ts2 = [e_ts2, norm(t2-t_gt)];
        e_angles2 = [e_angles2, rad2deg(rotm2eul(R2))' - angle_gt];
        inl_num2 = [inl_num2, sum(inliers2)];
        e_ts3 = [e_ts3, norm(t3-t_gt)];
        e_angles3 = [e_angles3, rad2deg(rotm2eul(R3))' - angle_gt];
        inl_num3 = [inl_num3, sum(inliers3)];
        e_ts4 = [e_ts4, norm(t4-t_gt)];
        e_angles4 = [e_angles4, rad2deg(rotm2eul(R4))' - angle_gt];
        inl_num4 = [inl_num4, sum(inliers4)];
        e_ts5 = [e_ts5, norm(t5-t_gt)];
        e_angles5 = [e_angles5, rad2deg(rotm2eul(R5))' - angle_gt];
        inl_num5 = [inl_num5, sum(inliers5)];
        
        if VISUALIZATION
            visualize_ms_matches(imsL1, imsL2, kL1_m, kL2_m, ...
                camidxL1_m, camidxL2_m);
            
            pause
            close all
            
            % visualize matches
            kL1_inl = kL1_m(:, inliers);
            camidxL1_inl = camidxL1_m(inliers);
            kL2_inl = kL2_m(:, inliers);
            camidxL2_inl = camidxL2_m(inliers);
            visualize_ms_matches(imsL1, imsL2, kL1_inl, kL2_inl, ...
                camidxL1_inl, camidxL2_inl);
            
            pause
            close all
            
            % visualize reprojections
            p3D1_inl = p3D1_m(:, inliers);
            visualize_ms_reprojection(imsL2, p3D1_inl, ...
                kL2_inl, camidxL2_inl, R, t, multistereo);
        end
    end
    
    e_3D2D =    [e_3D2D,    [mean(e_ts1); mean(abs(e_angles1), 2); mean(inl_num1)]];
    e_dec3D2D = [e_dec3D2D, [mean(e_ts2); mean(abs(e_angles2), 2); mean(inl_num2)]];
    e_3D3D =    [e_3D3D,    [mean(e_ts3); mean(abs(e_angles3), 2); mean(inl_num3)]];
    e_hybrid =  [e_hybrid,  [mean(e_ts4); mean(abs(e_angles4), 2); mean(inl_num4)]];
    e_single =  [e_single,  [mean(e_ts5); mean(abs(e_angles5), 2); mean(inl_num5)]];
end

%% plot trasl err
close all
bar([e_3D2D(1, :); e_hybrid(1, :); e_3D3D(1, :); e_dec3D2D(1, :); e_single(1, :)]');
xticklabels({'pose1', 'pose2', 'pose3', 'pose4', 'pose5', 'pose6', 'pose7', 'pose8', 'pose9'});
legend({'3D2D', 'Hybrid', '3D3D', 'dec3D2D', 'Best Stereo'}, 'Location', 'northwest');
ylabel('Translation Error (mm)');

%% plot angles
figure();
bar([e_3D2D(2, :); e_hybrid(2, :); e_3D3D(2, :); e_dec3D2D(2, :); e_single(2, :)]');
xticklabels({'pose1', 'pose2', 'pose3', 'pose4', 'pose5', 'pose6', 'pose7', 'pose8', 'pose9'});
legend({'3D2D', 'Hybrid', '3D3D', 'dec3D2D', 'Best Stereo'}, 'Location', 'northwest');
ylabel('Roll Error (deg)');
%%
figure();
bar([e_3D2D(3, :); e_hybrid(3, :); e_3D3D(3, :); e_dec3D2D(3, :); e_single(3, :)]');
xticklabels({'pose1', 'pose2', 'pose3', 'pose4', 'pose5', 'pose6', 'pose7', 'pose8', 'pose9'});
legend({'3D2D', 'Hybrid', '3D3D', 'dec3D2D', 'Best Stereo'}, 'Location', 'northwest');
ylabel('Pitch Error (deg)');
%%
figure();
bar([e_3D2D(4, :); e_hybrid(4, :); e_3D3D(4, :); e_dec3D2D(4, :); e_single(4, :)]');
xticklabels({'pose1', 'pose2', 'pose3', 'pose4', 'pose5', 'pose6', 'pose7', 'pose8', 'pose9'});
legend({'3D2D', 'Hybrid', '3D3D', 'dec3D2D', 'Best Stereo'}, 'Location', 'northwest');
ylabel('Yaw Error (deg)');

%% plot inliers
figure();
bar([e_3D2D(5, :); e_hybrid(5, :); e_3D3D(5, :); e_dec3D2D(5, :); e_single(5, :)]');
xticklabels({'pose1', 'pose2', 'pose3', 'pose4', 'pose5', 'pose6', 'pose7', 'pose8', 'pose9'});
legend({'3D2D', 'Hybrid', '3D3D', 'dec3D2D', 'Best Stereo'}, 'Location', 'northeast');
ylabel('Inliers Number');




% e_ts4 =
% 
%    97.8950   24.1192   55.2281   97.8951   33.1361
% 
% e_angles4
% 
% e_angles4 =
% 
%    -0.7272   -0.4890   -0.7063   -0.7272   -0.5367
%     0.1261    0.3166    0.3943    0.1261    0.3241
%     0.7853    0.3987    0.6273    0.7853    0.4452
