function [kps,descrs] = affineDetect_gpu(gi,SURF_THRESH,maxPoints)

    % SAMPLING PARAMETERS
    nt = 5;
    
    % creating params vector
    params = zeros(47, 2);
    
    ts = 2.^((0.5).*[0:nt]);
    for ii = 1:length(ts)
        t = ts(ii);
        ths = 0:72.0 / t:180;
        for th = ths
            params(ii,:) = [t,th];
        end
    end

    descrs = zeros(maxPoints, 64);
    kps = zeros(maxPoints, 2);
    valid_kpts_idx = 1;
    for idx=1:length(params)

        t = params(idx,1)
        th = params(idx,2)

        [imgS,Ai,mask] = affineSkew(t,th,gi);
        md = bwdist(~mask);


        sfs = detectSURFFeatures(imgS,'MetricThreshold',SURF_THRESH);
        [sds,~] = extractFeatures(imgS, sfs);

        % the detected keypoints have to be taken back to original
        % image coordinates
        sfs_pts = zeros(maxPoints, 2);
        fs = zeros(maxPoints, 64);
        
        valid_sfs_idx = 1;
        for i=1:length(sfs)
            e_pts = sfs(i).Location; 
            s = sfs(i).Scale;

            % check if the point is inside the boundaries of the
            % mask 6sqrt(2) scale (see original SIFT paper)
            if(md(floor(e_pts(2)),floor(e_pts(1))) >= 6*sqrt(2)*s)
                e_pts(1,3) = 1; %adding homogeneous
                o_pts = Ai * e_pts';
                sfs_pts(valid_sfs_idx,:) = round(o_pts(1:2));
                fs(valid_sfs_idx,:) = sds(i,:);
                valid_sfs_idx = valid_sfs_idx + 1;
            end
        end
        
        sfs_pts = sfs_pts(1:valid_sfs_idx-1, :);
        fs = fs(1:valid_sfs_idx-1, :);

        if valid_kpts_idx + valid_sfs_idx - 1 < maxPoints
            kps(valid_kpts_idx+1:valid_kpts_idx+valid_sfs_idx-1, :) = sfs_pts;
            descrs(valid_kpts_idx+1:valid_kpts_idx+valid_sfs_idx-1, :) = fs;
            valid_kpts_idx = valid_kpts_idx + valid_sfs_idx;
        else
            break;
        end
    end
    kps = kpts(1:valid_kpts_idx, :);
    descrs = descrs(1:valid_kpts_idx, :);
end             
