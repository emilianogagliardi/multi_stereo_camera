function [img,Ai,mask] = affineSkew(t,th,img)
    [h,w,~] = size(img);
    
    % creating a mask sized as original image
    mask = zeros(h,w);
    mask(:) = 255;
    
    % incrementally creating the transformation matrix
    A = [1,0,0;0,1,0];
    Ai = [1,0,0;0,1,0;0,0,1];
    
    if(th ~= 0)
        th = pi*th/180;
        s = sin(th);
        c = cos(th);
        A = [c,-s;  s, c; 0,0];
        corners = [0, 0; w, 0; w, h; 0, h];
        tcorners = floor( corners * A' );
        
        pgon = polyshape(tcorners(:,1)',tcorners(:,2)');
        [xlim,ylim] = boundingbox(pgon);
        
        x = min(xlim);
        y = min(ylim);
        w = max(xlim)-min(xlim);
        h = max(ylim)-min(ylim);
        
        % stack the third column
        A(:,3) = [-x,-y,1];
        tform = affine2d(A');
        img = imwarp(img,tform,'OutputView',imref2d([h,w]));
    end
    
    if (t ~= 1.0)
        s = 0.8*sqrt(t*t-1);
        
        img = imgaussfilt(img,[s,0.01]);
        
        img = imresize(img,[h,w * 1.0/t],'Method','Nearest');
        A(1,:) = A(1,:) ./ t;
    end
    
    if (th ~= 0.0 || t ~= 1.0)
        [h, w,~] = size(img);
        tform = affine2d(A');
        mask = imwarp(mask,tform,'OutputView',imref2d([h,w]));
        Ai = invert(tform);
        Ai = Ai.T';
    end  
    
end