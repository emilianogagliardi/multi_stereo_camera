 classdef CalibrationGraph
    %CALIBRATION_GRAPH build and optimize a pose graph for calibrating a system of stereo cameras with non overlapping field of view
    % Each stereo camera is assumed to be alredy calibrated with the matlab
    % stereo calibration app.
    % Take pictures of a known pattern, such that at least two stereo
    % camera are seeing it simultaneusly. Only a small part of the pattern
    % could be in the FOV of a single stereo. 
    
    
    properties %(Access = private)
        G
        CamsMap
        PointsNumber
        ErrsFun
        JacobPatternPosesFun
        JacobCamPosesFun
        JacobScaleFun
        VarIndicesCam
        VarIndicesPat
        CanCalibrate
        IsCalibrated
        Sol
    end
    
    %-------------------------PUBLIC-METHODS------------------------------
    methods
        function obj = CalibrationGraph(cams)
            %CALIBRATIONobj.GRAPH Construct an instance of this class
            % Parameters:
            % - cams: cell array containing stero cameras structs (ouput of
            % matlab calibration), in the order as their id in filenames
            % (see calibration_parameters.m)
            %
            % Builds a graph with cams and pattern poses nodes, a cam is 
            % connected to a pattern pose if the cam is seeing it in such 
            % pose. Each edge contains a rototraslation guess with respect 
            % cam0. 
            % obj.Graph nodes have fields 
            % Name, Type, Tguess.
            % example1: 'Cam0',     'c',    struct('R', ..., 't', ...)
            % example2: 'Pattern3', 'p',    struct('R', ..., 't', ...)
            
            calibration_parameters;

            obj.IsCalibrated = false;
            
            obj.G = graph();
            obj.PointsNumber = 0;
            
            obj = obj.AddCamsNodes(cams);
            
            obj = obj.AddPatternsNodes(cams);
            
            if obj.AllCamerasConnected()
                obj.CanCalibrate = true;
                obj = obj.ComputeGuesses();
                obj = obj.BuildErrorAndJacobian();
                fprintf('\nAll cams are connected, can calibrate\n');
            else
                obj.CanCalibrate = false;
                fprintf('\nCan not calibrate, some cameras are not connected\n');
            end
        end
        
        function [obj, multicalib] = Calibrate(obj)
            % obtain the guess as SO(3)xR, note that camera poses are
            % inverted
            if ~obj.CanCalibrate
                error ('Can not calibrate, some cameras are not connected');
            end
            
            if obj.IsCalibrated
                guess = obj.Sol;
            else
                guess = zeros(1, 6*(size(obj.G.Nodes, 1) - 1) + 1);
                guess(end) = 1; % scale guess is one
                for node_idx = 2:size(obj.G.Nodes, 1) % convert node pose guess in rodrigues
                    R = obj.G.Nodes.Tguess(node_idx).R;
                    t = obj.G.Nodes.Tguess(node_idx).t;
                    type = obj.G.Nodes.Type{node_idx};
                    var_idx = obj.NodeIdxToVarIdx(node_idx);
                    if type == 'c'
                        t = R'*t;
                        R = R';
                    end
                    r = rodrigues(R);
                    guess(var_idx) = [r', t'];
                end
            end
            
            options=optimset('Display','iter',...
                 'Jacobian','on',...
                 'MaxIter',2000, 'MaxFunEvals', 20000, 'TolX', 1e-8, 'TolFun', 1e-8, 'Algorithm', 'levenberg-marquardt'); 

            obj.Sol = lsqnonlin(@(x) obj.ObjectiveFunction(x), guess, [], [], options);
            
            obj.IsCalibrated = true;
            
            multicalib = obj.BuildOutput(obj.Sol);
            
            obj = obj.ComputeReprojectionErrors();
        end

        
        function [e, j] = ObjectiveFunction(obj, poses)
            % poses is a 1x(6*node_number) vector, where poses(6*i:6*(i+1))
            % is the pose of node i expressed as [wx, wy, wz, tx, ty, tz]
            % The variable order is the same as the nodes of the graph,
            % thus derivatives (for each element of the sum of square) need
            % to have the same order
            % the last element of poses is a scale factor to be applyed to
            % the 3D points
            
            % err in position i of ErrsFun corresponds to edge i in the
            % graph, and is funciton of the poses o the two connected nodes
            % (a camera and a pattern) node j and node k. 
            
            e = zeros(obj.PointsNumber*2, 1); % for each point x and y error
            
            global OPTIMIZE_SCALE;
            
            if OPTIMIZE_SCALE
                s = poses(end);
            else
                s = 1;
            end
            
            if nargout == 1 % compute only the error
                idx = 1;
                for ii = 1:size(obj.ErrsFun, 1)
                    err_fun = obj.ErrsFun{ii};
                    % obtain the arguments from the optimized var poses
                    pat_pose_arg = poses(obj.VarIndicesPat(ii, :));
                    % cam pose is variable idx are nan if the current
                    % camera node is cam0
                    if all(isnan(obj.VarIndicesCam(ii, :))) % cam0 edge, no cam variable
                        curr_e = err_fun([pat_pose_arg, s]);
                    else
                        cam_pose_arg = poses(obj.VarIndicesCam(ii, :));
                        curr_e = err_fun([cam_pose_arg, pat_pose_arg, s]);
                    end
                    n = size(curr_e, 1);
                    e(idx:idx+n-1) = curr_e;
                    idx = idx + n;
                end
            else % compute error and jacobian
                j = zeros(obj.PointsNumber*2, 6 * (size(obj.G.Nodes, 1)-1) + 1); % 6 var for each pose, 1 var for scale
                idx = 1;
                for ii = 1:size(obj.ErrsFun, 1)
                    err_fun = obj.ErrsFun{ii};
                    j_cam = obj.JacobCamPosesFun{ii};
                    j_pat = obj.JacobPatternPosesFun{ii};
                    j_scale = obj.JacobScaleFun{ii};
                    % obtain the arguments from the optimized var poses
                    pat_pose_arg = poses(obj.VarIndicesPat(ii, :));
                    % cam pose is variable idx are nan if the current
                    % camera node is cam0
                    if all(isnan(obj.VarIndicesCam(ii, :))) % cam0 edge, no cam variable
                        curr_e = err_fun([pat_pose_arg, s]);
                        curr_j_pat = j_pat([pat_pose_arg, s]);
                        if OPTIMIZE_SCALE
                            curr_j_scale = j_scale([pat_pose_arg, s]);
                        else
                            curr_j_scale = zeros(size(curr_e, 1), 1);
                        end
                        curr_j = [zeros(size(curr_e, 1), 6 * (size(obj.G.Nodes, 1)-1)), curr_j_scale];
                        curr_j(:, obj.VarIndicesPat(ii, :)) = curr_j_pat;
                    else
                        cam_pose_arg = poses(obj.VarIndicesCam(ii, :));
                        curr_e = err_fun([cam_pose_arg, pat_pose_arg, s]);
                        curr_j_pat = j_pat([cam_pose_arg, pat_pose_arg, s]);
                        curr_j_cam = j_cam([cam_pose_arg, pat_pose_arg, s]);
                        if OPTIMIZE_SCALE
                            curr_j_scale = j_scale([cam_pose_arg, pat_pose_arg, s]);
                        else
                            curr_j_scale = zeros(size(curr_e, 1), 1);
                        end
                        curr_j = [zeros(size(curr_e, 1), 6 * (size(obj.G.Nodes, 1)-1)), curr_j_scale];
                        curr_j(:, obj.VarIndicesPat(ii, :)) = curr_j_pat;
                        curr_j(:, obj.VarIndicesCam(ii, :)) = curr_j_cam;
                    end
                    n = size(curr_e, 1);
                    e(idx:idx+n-1) = curr_e;
                    j(idx:idx+n-1, :) = curr_j;
                    idx = idx + n;
                end
            end
        end
        
        function obj = RemoveBadPoints(obj, max_reproj)
            % after the object is calibrated, remove all the points with
            % reprojection error greather than max_reproj
            if ~obj.IsCalibrated
                error('To remove outliers first calibrate');
            end
            
            n_points_old = 0;
            n_points_new = 0;
            
            to_remove = [];
            
            for edge_idx = 1:size(obj.G.Edges, 1)
                
                nL = size(obj.G.Edges.edgeData(edge_idx).pts2DL, 1);
                nR = size(obj.G.Edges.edgeData(edge_idx).pts2DR, 1);
                
                valid_left = obj.G.Edges.edgeData(edge_idx).errorsL < max_reproj;
                valid_errL = obj.G.Edges.edgeData(edge_idx).errorsL(valid_left);
                valid_pts2DL = obj.G.Edges.edgeData(edge_idx).pts2DL(valid_left, :);
                valid_pts3DL = obj.G.Edges.edgeData(edge_idx).pts3DL(valid_left, :);
                obj.G.Edges.edgeData(edge_idx).errorsL = valid_errL;
                obj.G.Edges.edgeData(edge_idx).pts2DL = valid_pts2DL;
                obj.G.Edges.edgeData(edge_idx).pts3DL = valid_pts3DL;
                
                valid_right = obj.G.Edges.edgeData(edge_idx).errorsR < max_reproj;
                valid_errR = obj.G.Edges.edgeData(edge_idx).errorsR(valid_right);
                valid_pts2DR = obj.G.Edges.edgeData(edge_idx).pts2DR(valid_right, :);
                valid_pts3DR = obj.G.Edges.edgeData(edge_idx).pts3DR(valid_right, :);
                obj.G.Edges.edgeData(edge_idx).errorsR = valid_errR;
                obj.G.Edges.edgeData(edge_idx).pts2DR = valid_pts2DR;
                obj.G.Edges.edgeData(edge_idx).pts3DR = valid_pts3DR;
                
                disp(['Edge ', num2str(edge_idx), ':']);
                fprintf(['\tLeft:  kept ', num2str(sum(valid_left)), ' of ', num2str(nL), '\n']);
                fprintf(['\tRight: kept ', num2str(sum(valid_right)), ' of ', num2str(nR), '\n']);
                
                n_survived_points = sum(valid_left) + sum(valid_right);
                if n_survived_points  == 0
                    disp('Removed');
                    to_remove = [to_remove, edge_idx];
                end
                
                n_points_old = n_points_old + nL + nR;
                n_points_new = n_points_new + n_survived_points;

            end
            
            for ii = 1:numel(to_remove)
                obj.G = rmedge(obj.G, to_remove(ii));
                to_remove = to_remove - 1;
            end
            
            assert(n_points_old == obj.PointsNumber);
            obj.PointsNumber = n_points_new;
            
            % check if can calibrate
             if obj.AllCamerasConnected()
                obj.CanCalibrate = true;
                obj = BuildErrorAndJacobian(obj);
                fprintf('\nAll cams are connected, can calibrate\n');
            else
                obj.CanCalibrate = false;
                fprintf('\nCan not calibrate, some cameras are not connected\n');
            end
        end
        
    end % END OF PUBLIC METHODS
    
    %------------------------PRIVATE-METHODS-------------------------------
    methods%(Access = private)
        
        function obj = AddCamsNodes(obj, cams)
            % add camera nodes to the graph
            n_cams = numel(cams);
            % from 'Cam0' to 'CamN' where N is n_cams
            names = cellstr([repmat('Cam', n_cams, 1), int2str([0:n_cams-1]')])';
            % cams nodes have type field equal to 'c'
            types = cellstr(repmat('c', n_cams, 1));
            % pose guesses are unkown until the full graph is built. Set
            % Cam0 pose as the identity and leave the other nan
            guesses = [struct('R', eye(3), 't', zeros(3, 1)); ...
                       repmat(struct('R', nan(3), 't', nan(3, 1)), n_cams-1, 1)];
            obj.G = addnode(obj.G, names);
            obj.G.Nodes.Type = types;
            obj.G.Nodes.Tguess = guesses;            

            % build a Map containing camera node name and stereo
            % calibration data
            obj.CamsMap = containers.Map;
            for ii = 1:numel(names)
                obj.CamsMap(names{ii}) = cams{ii};
            end
        end
        
        function obj = AddPatternsNodes(obj, cams)
            % loop over the images (images directory defined in
            % calibration_parameters.m) and compute a position guess for
            % cameras with respect to pattern poses, store pattern nodes
            % and edges
            
            global IMAGES_FOLDER;
            global TEMPLATE;
            global MIN_INLIERS_SINGLE_IMAGE;
            global MAX_REPR_ERR_GUESS;
            
            % template keypoints needed for computing guesses
            [kpts_templ, desc_templ] = detect_keypoints(TEMPLATE, true);
            
            % get images filenames
            d = dir([IMAGES_FOLDER, '/*.png']);
            names = char(d.name);

            % obtain all the timestamps
            timestamps = [];
            for ii = 1:size(names, 1)
                [~, timestamp] = split_image_name(names(ii, :));
                timestamps = [timestamps, timestamp];
            end

            % loop over the timestamps, and add corresponding pattern pose node to the
            % graph, connected to each cam that is seeing it
            for ts = 0:max(timestamps)

                % find all images at timestamp ii
                images_d = d(timestamps == ts);

                % if the pattern is in the FOV of at least 2 cams
                if numel(images_d) >= 2 
                    disp(['*********** Timestamp ', num2str(ts), ' ***********']);
                    
                    % compute the guess over the images containing this
                    % pattern pose, and check if the result is ok. If the
                    % check is passed for at least two cameras then add the
                    % pattern pose node and connect it to the cameras
                    valid = [];
                    valid_camnodes = [];
                    valid_names = [];
                    for jj = 1:numel(images_d)
                        
                        [cam_idx, ts_cam] = split_image_name(images_d(jj).name);
                        assert(ts_cam == ts);

                        % compute a rototraslation guess
                        cam_node_idx = findnode(obj.G, ['Cam' num2str(cam_idx)]);
                        stereo_params = cams{cam_node_idx};
                        stereo_im = imread([images_d(jj).folder, '/', images_d(jj).name]);
                        
                        [R, t, pts2DL, pts2DR, pts3DL, pts3DR, repr_l, repr_r] = stereo_guess(stereo_params, kpts_templ, desc_templ, stereo_im);
                        fprintf(['Image ', images_d(jj).name, ':\n \t', ...
                            'left image inliers:  ', num2str(size(pts2DL, 1)), '\n \t', ...
                            'right image inliers: ', num2str(size(pts2DR, 1)), '\n \t', ...
                            'left image reprojection error:  ', num2str(repr_l), '\n \t', ...
                            'right image reprojection error: ', num2str(repr_r), '\n']);
                        
                        % check the number of inliers and the guess reprojection error
                        if (size(pts2DL, 1) >= MIN_INLIERS_SINGLE_IMAGE && repr_l < MAX_REPR_ERR_GUESS) || ...
                                (size(pts2DR, 1) >= MIN_INLIERS_SINGLE_IMAGE && repr_r < MAX_REPR_ERR_GUESS) % if one of the image in the stereo pair is valid
                            edgeData = struct('R', R, 't', t, 'pts2DL', [], 'pts2DR', [], 'pts3DL', []', 'pts3DR', []', ...
                                'errorsL', [], 'errorsR', []);
                            if size(pts2DL, 1) >= MIN_INLIERS_SINGLE_IMAGE && repr_l < MAX_REPR_ERR_GUESS % left is valid
                                edgeData.pts2DL = pts2DL;
                                edgeData.pts3DL = pts3DL';
                                edgeData.errorsL = nan(size(pts2DL, 1), 1);
                                disp('Left is valid');
                            else
                                disp('Left discarded');
                            end
                            if size(pts2DR, 1) >= MIN_INLIERS_SINGLE_IMAGE && repr_r < MAX_REPR_ERR_GUESS % right is valid
                                edgeData.pts2DR = pts2DR;
                                edgeData.pts3DR = pts3DR';
                                edgeData.errorsR = nan(size(pts2DR, 1), 1);
                                disp('Right is valid');
                            else
                                disp('Right discarded');
                            end
                            valid = [valid, edgeData];
                            valid_camnodes = [valid_camnodes, cam_node_idx];
                            valid_names = [valid_names, images_d(jj).name, ', '];
                        end
                    end
                    
                    
                    % if succesfull guess on at least two images add node
                    % and connect
                    if numel(valid) >= 2
                        fprintf(['Adding pattern pose ', num2str(ts), ' with images ', valid_names(1:end-2), '\n\n']);
                        obj.G = addnode(obj.G, table({['Pattern', num2str(ts)]}, 'p', struct('R', nan(3), 't', nan(3, 1)),...
                        'VariableNames', {'Name', 'Type', 'Tguess'}));
                        pattern_node_idx = findnode(obj.G, ['Pattern', num2str(ts)]);
                        for kk = 1:numel(valid)
                            edgeData = valid(kk);
                            cam_node_idx = valid_camnodes(kk);
                            obj.G = addedge(obj.G, cam_node_idx, pattern_node_idx, table(edgeData));
                            obj.PointsNumber = obj.PointsNumber + size(edgeData.pts2DL, 1) + size(edgeData.pts2DR, 1);
                        end
                    else
                        fprintf(['Pattern pose ', num2str(ts), ' discarded: not enough valid images to build a path Cam-Pattern-Cam\n\n']);
                    end
                end
            end
        end
        
        function obj = ComputeGuesses(obj)
            % compute for each node a rototraslation guess with respect to 
            % cam0, by composing the transformations on the path in the 
            % spanning tree of obj.G. Note that each edge connect a camera 
            % node and a pattern node, and the transformations always are 
            % from pattern to camera, thus some need to be inverted
            
            % obtain a tree
            root = findnode(obj.G, 'Cam0');
            T = minspantree(obj.G, 'Root', root);
            figure, plot(T);

            
            % this solution does a lot of useless computation, does not exploit the
            % tree structure and rototraslation composition but....who cares
            for g_node_idx = 1:size(obj.G.Nodes, 1)

                % obtain the nodes that connect the current node to the root 
                path = shortestpath(T, 1, g_node_idx);

                % loop over the nodes on the path, and compose rototraslations
                Trasf = eye(4);
                for jj = 1:length(path)-1
                    node1_idx = path(jj);
                    node2_idx = path(jj+1);
                    edge = findedge(T, node1_idx, node2_idx);
                    R = T.Edges.edgeData(edge).R;
                    t = T.Edges.edgeData(edge).t;
                    currT = [R, t; zeros(1, 3), 1];
                    if T.Nodes.Type{node2_idx} == 'c' % rototr on edges from pattern to cam need to be inverted
                        Trasf = Trasf / currT;
                    else 
                        Trasf = Trasf * currT;
                    end
                end

                obj.G.Nodes.Tguess(g_node_idx).R = Trasf(1:3, 1:3);
                obj.G.Nodes.Tguess(g_node_idx).t = Trasf(1:3, 4);
            end
        end
        
        function obj = BuildErrorAndJacobian(obj)
            % traverse the graph and build an array of error functions and
            % jacobian functions. 
            % See err_and_jacobian_wrapper.m and
            % err_and_jacobian_wrapper_cam0.m for description of the
            % functions
            % while building cell arrays of functions also fills the 
            % VarIndicesCam and VarIndicesPat matrices, in such a way that
            % the input of ErrsFun{ii} is 
            % [poses([VarIndicesCam(ii, :), VarIndicesPat(ii, :)]), s]
            % for images taken from a camera different from cam0, and
            % [poses(VarIndicesPat(ii, :)), s] for images taken by cam0
            n_edges = size(obj.G.Edges, 1);
            
            obj.ErrsFun = cell(n_edges, 1);
            obj.JacobCamPosesFun = cell(n_edges, 1);
            obj.JacobPatternPosesFun = cell(n_edges, 1);
            obj.JacobScaleFun = cell(n_edges, 1);
            
            obj.VarIndicesCam = nan(n_edges, 6);
            obj.VarIndicesPat = nan(n_edges, 6);
            
            for edge_idx = 1:n_edges
                edge_data = obj.G.Edges.edgeData(edge_idx);
                % get the nodes connected
                [node1_idx, node2_idx] = obj.ConnectedNodes(edge_idx);
                % find which one is the camera node
                if obj.G.Nodes.Type{node1_idx} == 'c'
                    cam_node_idx = node1_idx;
                    pat_node_idx = node2_idx;
                    assert(obj.G.Nodes.Type{node2_idx} == 'p');
                else
                    cam_node_idx = node2_idx;
                    pat_node_idx = node1_idx;
                    assert(obj.G.Nodes.Type{node2_idx} == 'c');
                    assert(obj.G.Nodes.Type{node1_idx} == 'p');
                end
                
                % check if the camera is Cam0
                if all(obj.G.Nodes.Name{cam_node_idx} == 'Cam0')
                    [e, jp, js] = err_and_jacobian_wrapper_cam0(...
                        obj.CamsMap(obj.G.Nodes.Name{cam_node_idx}), ...
                        edge_data.pts2DL, edge_data.pts2DR, ...
                        edge_data.pts3DL, edge_data.pts3DR);
                    
                    obj.ErrsFun{edge_idx} = e;
                    obj.JacobCamPosesFun{edge_idx} = nan; % cam0 pose is not a variable, no jacob defined
                    obj.JacobPatternPosesFun{edge_idx} = jp;
                    obj.JacobScaleFun{edge_idx} = js;
                    
                    obj.VarIndicesPat(edge_idx, :) = obj.NodeIdxToVarIdx(pat_node_idx);
                else
                    [e, jc, jp, js] = err_and_jacobian_wrapper(...
                        obj.CamsMap(obj.G.Nodes.Name{cam_node_idx}), ...
                        edge_data.pts2DL, edge_data.pts2DR, ...
                        edge_data.pts3DL, edge_data.pts3DR);
                    
                    obj.ErrsFun{edge_idx} = e;
                    obj.JacobCamPosesFun{edge_idx} = jc;
                    obj.JacobPatternPosesFun{edge_idx} = jp;
                    obj.JacobScaleFun{edge_idx} = js;
                    
                    obj.VarIndicesCam(edge_idx, :) = obj.NodeIdxToVarIdx(cam_node_idx);
                    obj.VarIndicesPat(edge_idx, :) = obj.NodeIdxToVarIdx(pat_node_idx);
                end
            end
        end
        
        function indices = NodeIdxToVarIdx(obj, node_idx)
            % get in input the graph node idx and returns the array of
            % indidces of the corresponding pose variable in the big
            % optimized variable vector
            
            % the variable vector has 6 elements for each node in the graph
            % except cam0, and the pose variable are concatenated in the
            % order of graph nodes indices
            % node2 -> 1:6
            % node3 -> 7:12
            % node4 -> 13:18
            % ....
            start = (node_idx - 2) * 6 + 1;
            indices = start : start+5;
        end
        
        function node_idx = VarIdxToNodeIdx(obj, indices)
            start = indices(1);
            node_idx = ((start - 1) / 6) + 2;
        end
        
        function [node1_idx, node2_idx] = ConnectedNodes(obj, edge_idx)
            node1_idx = findnode(obj.G, obj.G.Edges.EndNodes{edge_idx, 1});
            node2_idx = findnode(obj.G, obj.G.Edges.EndNodes{edge_idx, 2});
        end
        
        function test = AllCamerasConnected(obj)
            test = true;
            n_cams = obj.CamsMap.Count;
            for ii = 1:n_cams
                if isempty(outedges(obj.G, ii))
                    test = false;
                    break;
                end
            end
        end
        
        function poses = BuildOutput(obj, sol)
            n_cams = obj.CamsMap.Count;
            poses = zeros(4, 4, n_cams);
            poses(:, :, 1) = eye(4); % cam0
            for ii = 2:n_cams
                poses(:, :, ii) = inv(SO3R3_2_T(sol(obj.NodeIdxToVarIdx(ii))));
            end
        end
        
        function [obj, camera_idx, image_idx, all_errs] = ComputeReprojectionErrors(obj)
            
            disp('Optimization result:');
            sol = obj.Sol;
            s = sol(end);
            
            image_idx = [];
            camera_idx = [];
            all_errs = [];
            
            for edge_idx = 1:size(obj.G.Edges, 1)
                % get the nodes connected
                [node1_idx, node2_idx] = obj.ConnectedNodes(edge_idx);
                % find which one is the camera node
                if obj.G.Nodes.Type{node1_idx} == 'c'
                    cam_node_idx = node1_idx;
                    pat_node_idx = node2_idx;
                    assert(obj.G.Nodes.Type{node2_idx} == 'p');
                else
                    cam_node_idx = node2_idx;
                    pat_node_idx = node1_idx;
                    assert(obj.G.Nodes.Type{node2_idx} == 'c');
                    assert(obj.G.Nodes.Type{node1_idx} == 'p');
                end
                
                [Kl, Kr, Rs, ts] =  convert_stereo_params(obj.CamsMap(obj.G.Nodes.Name{cam_node_idx}));
                pts2DL = obj.G.Edges.edgeData(edge_idx).pts2DL;
                pts2DR = obj.G.Edges.edgeData(edge_idx).pts2DR;
                pts3DL = obj.G.Edges.edgeData(edge_idx).pts3DL;
                pts3DR = obj.G.Edges.edgeData(edge_idx).pts3DR;
                p = sol(obj.NodeIdxToVarIdx(pat_node_idx));
                Tp = SO3R3_2_T(p); % pattern pose
                
                % check if the camera is Cam0
                if all(obj.G.Nodes.Name{cam_node_idx} == 'Cam0')
                    Tl = Tp;
                else
                    c = sol(obj.NodeIdxToVarIdx(cam_node_idx));
                    Tc = SO3R3_2_T(c); % camera inverse pose
                    Tl = Tc * Tp;
                end
                
                disp(['Edge ' num2str(edge_idx), ':']);
                if size(pts2DL, 1) > 0
                    Rl = Tl(1:3, 1:3); tl = Tl(1:3, 4);
                    [~, errsL] = reprojection(Kl, Rl, tl, s * pts3DL, pts2DL);
                    obj.G.Edges.edgeData(edge_idx).errorsL = errsL;
                    all_errs = [all_errs, errsL];
                    camera_idx =  [camera_idx, repmat([cam_node_idx], 1, length(errsL))];
                    image_idx =  [image_idx, repmat([edge_idx], 1, length(errsL))];
                    fprintf('\tLeft camera:\n')
                    fprintf(['\t\tNum points:', num2str(length(errsL)), '\n']);
                    fprintf(['\t\tMean error:', num2str(mean(errsL)), '\n']);
                    fprintf(['\t\tMax error:', num2str(max(errsL)), '\n']);
                end
                
                if size(pts2DR, 1) > 0
                    Tr = [Rs, ts; zeros(1, 3), 1] * Tl;
                    Rr = Tr(1:3, 1:3); tr = Tr(1:3, 4);
                    [~, errsR] = reprojection(Kr, Rr, tr, s * pts3DR, pts2DR);
                    obj.G.Edges.edgeData(edge_idx).errorsR = errsR;
                    all_errs = [all_errs, errsR];
                    camera_idx =  [camera_idx, repmat([cam_node_idx], 1, length(errsR))];
                    image_idx =  [image_idx, repmat([edge_idx], 1, length(errsR))];
                    fprintf('\tRight camera:\n')
                    fprintf(['\t\tNum points:', num2str(length(errsR)), '\n']);
                    fprintf(['\t\tMean error:', num2str(mean(errsR)), '\n']);
                    fprintf(['\t\tMax error:', num2str(max(errsR)), '\n']);
                end
                
                disp('');
            end
        end
    end 
end

