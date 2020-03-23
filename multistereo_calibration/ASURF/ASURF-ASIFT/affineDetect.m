function [kps,descrs] = affineDetect(gi,SURF_THRESH,parallel)
    
    %SURF_THRESH = 500;
    
    % SAMPLING PARAMETERS
    nt = 5;
    
    %creating params vector
    params = [];
    
    ts = 2.^((0.5).*[0:nt]);
    for t = ts
        ths = 0:72.0 / t:180;
        for th = ths
            params(end+1,:) = [t,th];
        end
    end
    
    kps = [];
    descrs = [];
    
    if(parallel == true)
        p = gcp();
        for i=1:length(params)
            param = params(i,:);
            f(i) = parfeval(p,@computeFeatures,2,gi,param,SURF_THRESH); 
        end       

        % Collect the results as they become available.
        for i=1:length(params)
          % fetchNext blocks until next results are available.
          [~,pts,ds] = fetchNext(f);
          kps = vertcat(kps,pts);
          descrs = vertcat(descrs,ds);
        end
        
    else
        
        for idx=1:length(params)
            
            % params(idx,:)
            t = params(idx,1)
            th = params(idx,2)
            
            [imgS,Ai,mask] = affineSkew(t,th,gi);
            md = bwdist(~mask);
            
            %imgSG = rgb2gray(imgS);

            sfs = detectSURFFeatures(imgS,'MetricThreshold',SURF_THRESH);
            [sds,~] = extractFeatures(imgS,sfs);

            %[fs1,~] = extractFeatures(imgS(:,:,1),sfs);
            %[fs2,~] = extractFeatures(imgS(:,:,2),sfs);
            %[fs3,~] = extractFeatures(imgS(:,:,3),sfs);
            %fs = [fs1,fs2,fs3];

            % the detected keypoints have to be taken back to original
            % image coordinates
            sfs_pts = [];
            fs = [];

            for i=1:length(sfs)
                e_pts = sfs(i).Location; 
                s = sfs(i).Scale;

                % check if the point is inside the boundaries of the
                % mask 6sqrt(2) scale (see original SIFT paper)
                if(md(floor(e_pts(2)),floor(e_pts(1))) >= 6*sqrt(2)*s)
                    e_pts(1,3) = 1; %adding homogeneous
                    o_pts = Ai * e_pts';
                    sfs_pts(end+1,:) = round(o_pts(1:2));
                    fs(end+1,:) = sds(i,:);
                end
            end
            
            kps = vertcat(kps,sfs_pts);
            descrs = vertcat(descrs,fs);
                
            %{   
            elseif(strcmp(method,'SIFT') || (strcmp(method,'PHOW')))
                
                if(strcmp(method,'SIFT'))
                    % detect SIFT features
                    imgS = rgb2gray(imgS);
                    [sfs,fs] = vl_sift(single(imgS)); %
                else
                    [sfs,fs]=vl_phow(im2single(imgS),'Color','rgb','Sizes', [4 8]) ;
                end
                
                % the detected keypoints have to be taken back to original
                % image coordinates
                sfs_pts = zeros(length(sfs),2);
                for i=1:length(sfs)
                    e_pts = [sfs(1,i),sfs(2,i)]; 
                    e_pts(1,3) = 1; %adding homogeneous
                    o_pts = Ai * e_pts';
                    sfs_pts(i,:) = round(o_pts(1:2));
                end
            end
            %}
        end
    end             
end

function [pts,ds] = computeFeatures(imgGray,param,SURF_THRESH)
    t = param(1);
    th = param(2);
    [imgS,Ai,mask] = affineSkew(t,th,imgGray);
    md = bwdist(~mask);
      
    sfs = detectSURFFeatures(imgS,'MetricThreshold',SURF_THRESH);
    %sfs = detectHarrisFeatures(imgS);
    [sds,~] = extractFeatures(imgS,sfs,'Method','SURF');

    % the detected keypoints have to be taken back to original
    % image coordinates
    pts = [];
    ds = [];

    for i=1:length(sfs)
        e_pts = sfs(i).Location; 
        s = sfs(i).Scale;
        %s = 1.6;

        % check if the point is inside the boundaries of the
        % mask 6sqrt(2) scale (see original SIFT paper)
        if(md(floor(e_pts(2)),floor(e_pts(1))) >= 6*sqrt(2)*s)
            e_pts(1,3) = 1; %adding homogeneous
            o_pts = Ai * e_pts';
            pts(end+1,:) = o_pts(1:2);
            ds(end+1,:) = sds(i,:);
        end
    end
    
end