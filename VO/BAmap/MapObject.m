classdef MapObject
    %MAP_OBJECT A set of Landmarks and camera poses that can be optimized with bundle adjustment
    
    properties
        Landmarks           % containers.map object of Landmarks
        Extrinsics          % 4x4xN, homog extrinsics at each timestamp
        Timestamp           % ts of the last inserted frame
        Multistereo         % MultiStereo object
        PointsNumber
    end
    
    methods
        function obj = MapObject(multistereo)
            %MAP Construct an instance of this class
            obj.Landmarks = containers.Map('KeyType', 'int32', ...
                'ValueType', 'any');
            obj.Timestamp = 0;
            obj.Extrinsics = [];
            obj.Multistereo = multistereo;
            obj.PointsNumber = 0;
        end
        
        function obj = AddFrame(obj, tracked, untracked, R, t)
            %ADD_FRAME add a new perception with new pose to the map
            % Parameters
            % - tracked: struct with fields kpts, cam_idx, landmarks_ID.
            %   kpts is an array of image points, cam_idx is an array of
            %   indices that in position ii contains the indices of the
            %   camera perceiving kpts(:, ii), landmarks_ID is an array of
            %   IDs of the corresponding 3D point, already present in the
            %   map
            % - untracked: struct with fields p3D, kpts, cam_idx,
            %   landmarks_ID.
            %   p3D are new 3D points to be added to the map, kpts and
            %   cam_idx the corresponding image points and perceiving
            %   camera indices, landmark_ID are new ID for the new
            %   landmarks
            % - R, t: extrinsics at time t with respect time t-1
            obj.Timestamp
            obj.Timestamp = obj.Timestamp + 1;
            
            if obj.Timestamp == 1 % first frame
                obj.Extrinsics(:, :, 1) = eye(4);
            else
                obj.Extrinsics(:, :, obj.Timestamp) = ...
                    [R, t; zeros(1, 3), 1] * obj.Extrinsics(:, :, obj.Timestamp - 1);
            end
            
            % add a new perception for the tracked landmarks
            for ii = 1:size(tracked.kpts, 2)
                ID = tracked.landmarks_ID(ii);
                kpt = tracked.kpts(:, ii);
                cam_idx = tracked.cam_idx(ii);
                obj.Landmarks(ID) = obj.Landmarks(ID).AddPerception(kpt, ...
                    obj.Timestamp, cam_idx);
            end
            
            % add new landmarks with their first perception
            for ii = 1:size(untracked.kpts, 2)
                ID = untracked.landmarks_ID(ii);
                % TODO assert ID not present in map keys
                kpt = untracked.kpts(:, ii);
                cam_idx = untracked.cam_idx(ii);
                % bring 3D points from current camera frame to world frame
                R_ = obj.Extrinsics(1:3, 1:3, obj.Timestamp)';
                t_ = -R_ * obj.Extrinsics(1:3, 4, obj.Timestamp);
                p3D = R_ * untracked.p3D(:, ii) + t_;
                obj.Landmarks(ID) = Landmark(p3D, kpt, obj.Timestamp, cam_idx);
                obj.PointsNumber = obj.PointsNumber + 1;
            end
        end
        
        function obj = Optimize(obj, windowSize)
            % OPTIMIZE: perform bundle adjustment over the last window_size
            % frames
            
            % do not optimize if not enough time has passed
            % first pose must not be optimized
            if obj.Timestamp <= windowSize 
                return;
            end
            
            % retrieve optimized poses
            optimizedPoses = zeros(6*windowSize, 1);
            for ii = 1:windowSize
                R = obj.Extrinsics(1:3, 1:3, obj.Timestamp - windowSize + ii);
                t = obj.Extrinsics(1:3, 4, obj.Timestamp - windowSize + ii);
                w = rodrigues(R);
                optimizedPoses(6*(ii-1)+1: 6*ii) = [w; t];
            end
            % select all the landmarks with some perceptions inside the
            % time window [obj.Timestamp, obj.Timestamp-window_size]
            % landmarks coordinates are stacked in a 1x3N vector where
            % N is the number of landmarks. To retrieve coordinates after
            % optimization we also build an array of landmarks ID in the
            % same order as the one in the vector
            optimizedIDs = [];
            optimizedCoords = [];
            perceptions = {};
            percNumber = 0;
            pointsNumber = 0;
            for ID_ = obj.Landmarks.keys()
                ID = ID_{1};
                landmark = obj.Landmarks(ID);
                percs = landmark.Perceptions;
                validPercs = [];
                % retrieve percs inside window time
                for ii = 1:length(percs)
                    if percs(ii).Ts > obj.Timestamp - windowSize
                        validPercs = [validPercs, percs(ii)];
                    end
                end
                if length(validPercs) < 2 % do not optimize this point
                    continue;
                end
                % enough perceptions for this landmark in this time window,
                % add to optimized landmarks 
                percNumber = percNumber + length(validPercs);
                pointsNumber = pointsNumber + 1;
                perceptions{end+1} = validPercs;
                optimizedCoords = [optimizedCoords; landmark.P3D];
                optimizedIDs = [optimizedIDs, ID_{1}];
            end

            guess = double([optimizedPoses; optimizedCoords]);
            objfun = @(x) obj.objectiveFunction(x, windowSize, perceptions, percNumber, pointsNumber);
            
            options=optimset('Display','iter',...
                 'Jacobian','off',...
                 'MaxIter',10, 'MaxFunEvals', 100000, 'TolX', 1e-4, ...
                 'TolFun', 1e-3, 'Algorithm', 'levenberg-marquardt'); 

            sol = lsqnonlin(objfun, guess, [], [], options);
            
            % retrieve optimized poses and points from solution and store
            % them again in map object
            poses = sol(1:6*windowSize);
            for ii = 1:windowSize
                pose = poses(6*(ii-1)+1:6*(ii-1)+6);
                R = rodrigues(pose(1:3));
                t = pose(4:6);
                T = [R, t; zeros(1, 3), 1];
                obj.Extrinsics(:, :, obj.Timestamp-windowSize+ii) = T;
            end
            points = sol(6*windowSize+1:end);
            points = reshape(points, 3, pointsNumber);
            for ii = 1:size(points, 2)
                p = obj.Landmarks(optimizedIDs(ii));
                p.P3D = points(:, ii);
                obj.Landmarks(optimizedIDs(ii)) = p;
            end
        end
        
        function [e, j] = objectiveFunction(obj, optVar, windowSize, perceptions, percsNumber, pointNumber)
            
            % retrieve poses from optimized variable
            poses = reshape(optVar(1:6*windowSize), 6, windowSize);
            % retrieve points from optimized variable
            points = reshape(optVar(6*windowSize+1:end), 3, pointNumber);
            
            e = nan(2*percsNumber, 1);
            
            percIdx = 1;
            
            if nargout == 1 % compute only error
                % for each 3D point
                for p3Didx = 1:pointNumber
                    p3D = points(:, p3Didx);
                    p3Dperceptions = perceptions{p3Didx}; 
                    % compute reprojection error for each perception of the
                    % current 3D point
                    for p3DpercIdx = 1:length(p3Dperceptions)
                        perc = p3Dperceptions(p3DpercIdx);
                        p2D = perc.P2D; % image point
                        camIdx = perc.CamIdx; % perceiving camera index
                        ts = perc.Ts;% time instant at which has perceived
                        % obtain pose index in poses matrix from perception
                        % timestamp
                        poseIdx = ts -(obj.Timestamp - windowSize);
                        pose = poses(:, poseIdx); % extrinsics guess in ts
                        % x, y reprojection error
                        e_ = ej_wrapper_BA(obj.Multistereo, camIdx, pose(1:3), ...
                            pose(4:6), p3D, p2D);
                        e(2*(percIdx-1)+1:2*percIdx) = e_;
                        percIdx = percIdx + 1;
                    end
                end
            else % compute error and jacobian
                j = zeros(2*percsNumber, length(optVar));
                % for each 3D point
                for p3Didx = 1:pointNumber
                    p3D = points(:, p3Didx);
                    p3Dperceptions = perceptions{p3Didx}; 
                    % compute reprojection error for each perception of the
                    % current 3D point
                    for p3DpercIdx = 1:length(p3Dperceptions)
                        perc = p3Dperceptions(p3DpercIdx);
                        p2D = perc.P2D; % image point
                        camIdx = perc.CamIdx; % perceiving camera index
                        ts = perc.Ts; % time instant at which has perceived
                        % obtain pose index in poses matrix from perception
                        % timestamp
                        poseIdx = ts -(obj.Timestamp - windowSize);
                        pose = poses(:, poseIdx); % extrinsics guess in ts
                        % x, y reprojection error
                        [e_, jpose, jpoint] = ej_wrapper_BA(obj.Multistereo, ...
                            camIdx, pose(1:3), pose(4:6), p3D, p2D);
                        e(2*(percIdx-1)+1:2*percIdx) = e_;
                        % add jacobian of pose
                        j(2*(percIdx-1)+1:2*percIdx, 6*(poseIdx-1)+1:6*poseIdx) = jpose;
                        % add jacobian of point
                        j(2*(percIdx-1)+1:2*percIdx, 6*windowSize + (3*(p3Didx-1)+1:3*p3Didx)) = jpoint;
                        percIdx = percIdx + 1;
                    end
                end
            end
        end
        
        function  obj = OptimizePoseOnly(obj, windowSize)
            % OPTIMIZE: perform bundle adjustment over the last window_size
            % frames without landmarks optimization
            
            % do not optimize if not enough time has passed
            % first pose must not be optimized
            if obj.Timestamp <= windowSize 
                return;
            end
            
            % retrieve optimized poses
            optimizedPoses = zeros(6*windowSize, 1);
            for ii = 1:windowSize
                R = obj.Extrinsics(1:3, 1:3, obj.Timestamp - windowSize + ii);
                t = obj.Extrinsics(1:3, 4, obj.Timestamp - windowSize + ii);
                w = rodrigues(R);
                if all(w==0) % apply perturbation
                    w = rand(3, 1)*eps;
                end
                optimizedPoses(6*(ii-1)+1: 6*ii) = [w; t];
            end
            % select all the landmarks with some perceptions inside the
            % time window [obj.Timestamp, obj.Timestamp-window_size]
            % landmarks coordinates are stacked in a 1x3N vector where
            % N is the number of landmarks
            landmarksCoords = [];
            perceptions = {};
            percNumber = 0;
            pointsNumber = 0;
            for ID_ = obj.Landmarks.keys()
                ID = ID_{1};
                landmark = obj.Landmarks(ID);
                percs = landmark.Perceptions;
                validPercs = [];
                % retrieve percs inside window time
                for ii = 1:length(percs)
                    if percs(ii).Ts > obj.Timestamp - windowSize
                        validPercs = [validPercs, percs(ii)];
                    end
                end
                if length(validPercs) < 2 % do not optimize this point
                    continue;
                end
                % enough perceptions for this landmark in this time window,
                % add to optimized landmarks 
                percNumber = percNumber + length(validPercs);
                pointsNumber = pointsNumber + 1;
                perceptions{end+1} = validPercs;
                landmarksCoords = [landmarksCoords, landmark.P3D];
            end
            
            guess = double(optimizedPoses);
            objfun = @(x) obj.objectiveFunctionPoseOnly(x, windowSize, landmarksCoords, perceptions, percNumber, pointsNumber);
            
            options=optimset('Display','iter',...
                 'Jacobian','on',...
                 'MaxIter',100, 'MaxFunEvals', 100000, 'TolX', 1e-5, ...
                 'TolFun', 1e-5, 'Algorithm', 'levenberg-marquardt'); 

            sol = lsqnonlin(objfun, guess, [], [], options);
            
            % retrieve optimized poses and points from solution and store
            % them again in map object
            poses = sol(1:6*windowSize);
            for ii = 1:windowSize
                pose = poses(6*(ii-1)+1:6*(ii-1)+6);
                R = rodrigues(pose(1:3));
                t = pose(4:6);
                T = [R, t; zeros(1, 3), 1];
                obj.Extrinsics(:, :, obj.Timestamp-windowSize+ii) = T;
            end
        end
        
        function [e, j] = objectiveFunctionPoseOnly(obj, optVar, windowSize, ...
                points, perceptions, percsNumber, pointNumber)
            % retrieve poses from optimized variable
            poses = reshape(optVar(1:6*windowSize), 6, windowSize);
            
            e = nan(2*percsNumber, 1);
            
            percIdx = 1;
            
            if nargout == 1 % compute only error
                % for each 3D point
                for p3Didx = 1:pointNumber
                    p3D = points(:, p3Didx);
                    p3Dperceptions = perceptions{p3Didx}; 
                    % compute reprojection error for each perception of the
                    % current 3D point
                    for p3DpercIdx = 1:length(p3Dperceptions)
                        perc = p3Dperceptions(p3DpercIdx);
                        p2D = perc.P2D; % image point
                        camIdx = perc.CamIdx; % perceiving camera index
                        ts = perc.Ts;% time instant at which has perceived
                        % obtain pose index in poses matrix from perception
                        % timestamp
                        poseIdx = ts -(obj.Timestamp - windowSize);
                        pose = poses(:, poseIdx); % extrinsics guess in ts
                        % x, y reprojection error
                        e_ = ej_wrapper_BA(obj.Multistereo, camIdx, pose(1:3), ...
                            pose(4:6), p3D, p2D);
                        e(2*(percIdx-1)+1:2*percIdx) = e_;
                        percIdx = percIdx + 1;
                    end
                end
            else % compute error and jacobian
                j = zeros(2*percsNumber, length(optVar));
                % for each 3D point
                for p3Didx = 1:pointNumber
                    p3D = points(:, p3Didx);
                    p3Dperceptions = perceptions{p3Didx}; 
                    % compute reprojection error for each perception of the
                    % current 3D point
                    for p3DpercIdx = 1:length(p3Dperceptions)
                        perc = p3Dperceptions(p3DpercIdx);
                        p2D = perc.P2D; % image point
                        camIdx = perc.CamIdx; % perceiving camera index
                        ts = perc.Ts; % time instant at which has perceived
                        % obtain pose index in poses matrix from perception
                        % timestamp
                        poseIdx = ts -(obj.Timestamp - windowSize);
                        pose = poses(:, poseIdx); % extrinsics guess in ts
                        % x, y reprojection error
                        [e_, jpose] = ej_wrapper_BA(obj.Multistereo, ...
                            camIdx, pose(1:3), pose(4:6), p3D, p2D);
                        e(2*(percIdx-1)+1:2*percIdx) = e_;
                        % add jacobian of pose
                        j(2*(percIdx-1)+1:2*percIdx, 6*(poseIdx-1)+1:6*poseIdx) = jpose;
                        percIdx = percIdx + 1;
                    end
                end
            end
        end
        
        function errors = computeReprojectionErrors(obj)
            errors = [];
            for ii = 1:length(obj.Landmarks)
                p3D = obj.Landmarks(ii).P3D;
                percs = obj.Landmarks(ii).Perceptions;
                for jj = 1:length(percs)
                    if length(percs) == 1
                        perc = percs;
                    else
                        perc = percs(jj);
                    end
                    p2D = perc.P2D;
                    ts = perc.Ts;
                    camIdx = perc.CamIdx;
                    R = obj.Extrinsics(1:3, 1:3, ts);
                    t = obj.Extrinsics(1:3, 4, ts);
                    K = obj.Multistereo.getKs(camIdx);
                    [Rc, tc] = obj.Multistereo.Tprincipal2cam(camIdx);
                    proj = K * [Rc * R, Rc * t + tc] * [p3D;1];
                    proj = proj(1:2)./proj(3);
                    err = norm(proj-p2D);
                    errors = [errors, err];
                    if err > 8
                        ii
                        jj
                        err
                    end
                end
            end
        end
        
        function [Rs, ts] = GetPoses(obj)
            Rs = zeros(3, 3, size(obj.Extrinsics, 3));
            ts = zeros(3, size(obj.Extrinsics, 3));
            for ii =1:size(obj.Extrinsics, 3)
                Rs(:, :, ii) = obj.Extrinsics(1:3, 1:3, ii)';
                t = -obj.Extrinsics(1:3, 1:3, ii)' * obj.Extrinsics(1:3, 4, ii);
                ts(:, ii) = t;
            end
        end
    end
    
end

