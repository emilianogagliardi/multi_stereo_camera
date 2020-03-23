function [e, jp, js] = err_and_jacobian_wrapper_cam0(stereo, pts2DL, pts2DR, pts3DL, pts3DR)
%ERR_AND_JACOBIAN_WRAPPER compute error and jacobian
    
% same as err_and_jacobian_wrapper, but the returned functions takes in
% input only the pattern pose and the scale

    [Kl, Kr, Rs, ts] = convert_stereo_params(stereo);
    ts1 = ts(1); ts2 = ts(2); ts3 = ts(3);
    rs1_1 = Rs(1, 1); rs1_2 = Rs(1, 2); rs1_3 = Rs(1, 3);
    rs2_1 = Rs(2, 1); rs2_2 = Rs(2, 2); rs2_3 = Rs(2, 3);
    rs3_1 = Rs(3, 1); rs3_2 = Rs(3, 2); rs3_3 = Rs(3, 3);

    elx = @(p) ...
        err_left_x_cam0(pts3DL(:, 1), pts3DL(:, 2), Kl(1, 1), p(4), p(6), p(1), p(2), p(3), Kl(1, 3), p(7), pts2DL(:, 1));
    
    ely = @(p) ...
        err_left_y_cam0(pts3DL(:, 1), pts3DL(:, 2), Kl(2, 2), p(5), p(6), p(1), p(2), p(3), Kl(2, 3), p(7), pts2DL(:, 2));
    
    erx = @(p) ...
        err_right_x_cam0(pts3DR(:, 1), pts3DR(:, 2), Kr(1, 1), p(4), p(5), p(6), p(1), p(2), p(3), Kr(1, 3), rs1_1, rs1_2, rs1_3, rs3_1, rs3_2, rs3_3, p(7), ts1, ts3, pts2DR(:, 1));
    
    ery = @(p) ...
        err_right_y_cam0(pts3DR(:, 1), pts3DR(:, 2), Kr(2, 2), p(4), p(5), p(6), p(1), p(2), p(3), Kr(2, 3), rs2_1, rs2_2, rs2_3, rs3_1, rs3_2, rs3_3, p(7), ts2, ts3, pts2DR(:, 2));

    
    if size(pts2DL, 1) > 0 && size(pts2DR, 1) > 0
        e = @(x) [elx(x); ely(x); erx(x); ery(x)];
    elseif size(pts2DL, 1) > 0
        e = @(x) [elx(x); ely(x)];
    else
        e = @(x) [erx(x); ery(x)];
    end
    

    
    jlpx = @(X,Y,ptx,ptz,pwx,pwy,pwz,s) ...
        jacob_left_patpose_x_cam0(X,Y,Kl(1, 1),ptx,ptz,pwx,pwy,pwz,Kl(1, 3),s);
    jlpx = @(p) cell2mat(arrayfun(jlpx, pts3DL(:, 1), pts3DL(:, 2), ...
        repmat(p(4), [size(pts3DL, 1), 1]), repmat(p(6), [size(pts3DL, 1), 1]), repmat(p(1), [size(pts3DL, 1), 1]), repmat(p(2), [size(pts3DL, 1), 1]), repmat(p(3), [size(pts3DL, 1), 1]), repmat(p(7), [size(pts3DL, 1), 1]), 'UniformOutput', false));

    jlpy = @(X,Y,pty,ptz,pwx,pwy,pwz,s) ...
        jacob_left_patpose_y_cam0(X,Y,Kl(2, 2),pty,ptz,pwx,pwy,pwz,Kl(2, 3),s);
    jlpy = @(p) cell2mat(arrayfun(jlpy, pts3DL(:, 1), pts3DL(:, 2), ...
        repmat(p(5), [size(pts3DL, 1), 1]), repmat(p(6), [size(pts3DL, 1), 1]), repmat(p(1), [size(pts3DL, 1), 1]), repmat(p(2), [size(pts3DL, 1), 1]), repmat(p(3), [size(pts3DL, 1), 1]), repmat(p(7), [size(pts3DL, 1), 1]), 'UniformOutput', false));
     
    jrpx = @(X,Y,ptx,pty,ptz,pwx,pwy,pwz,s) ...
        jacob_right_patpose_x_cam0(X,Y,Kr(1, 1),ptx,pty,ptz,pwx,pwy,pwz,Kr(1, 3),rs1_1,rs1_2,rs1_3,rs3_1,rs3_2,rs3_3,s,ts1,ts3);
    jrpx = @(p) cell2mat(arrayfun(jrpx, pts3DR(:, 1), pts3DR(:, 2), ...
        repmat(p(4), [size(pts3DR, 1), 1]), repmat(p(5), [size(pts3DR, 1), 1]), repmat(p(6), [size(pts3DR, 1), 1]), repmat(p(1), [size(pts3DR, 1), 1]), repmat(p(2), [size(pts3DR, 1), 1]), repmat(p(3), [size(pts3DR, 1), 1]), repmat(p(7), [size(pts3DR, 1), 1]), 'UniformOutput', false));

    jrpy = @(X,Y,ptx,pty,ptz,pwx,pwy,pwz,s) ...
        jacob_right_patpose_y_cam0(X,Y,Kr(2, 2),ptx,pty,ptz,pwx,pwy,pwz,Kr(2, 3),rs2_1,rs2_2,rs2_3,rs3_1,rs3_2,rs3_3,s,ts2,ts3);
    jrpy = @(p) cell2mat(arrayfun(jrpy, pts3DR(:, 1), pts3DR(:, 2), ...
        repmat(p(4), [size(pts3DR, 1), 1]), repmat(p(5), [size(pts3DR, 1), 1]), repmat(p(6), [size(pts3DR, 1), 1]), repmat(p(1), [size(pts3DR, 1), 1]), repmat(p(2), [size(pts3DR, 1), 1]), repmat(p(3), [size(pts3DR, 1), 1]), repmat(p(7), [size(pts3DR, 1), 1]), 'UniformOutput', false));

    if size(pts2DL, 1) > 0 && size(pts2DR, 1) > 0
        jp = @(p) [jlpx(p); jlpy(p); jrpx(p); jrpy(p)];
    elseif size(pts2DL, 1) > 0
        jp = @(p) [jlpx(p); jlpy(p)];
    else
        jp = @(p) [jrpx(p); jrpy(p)];
    end
    
    jlsx = @(p) ...
        jacob_left_s_x_cam0(pts3DL(:, 1),pts3DL(:, 2),Kl(1, 1),p(4),p(6),p(1),p(2),p(3),Kl(1, 3),p(7));
    
    jlsy = @(p) ...
        jacob_left_s_y_cam0(pts3DL(:, 1),pts3DL(:, 2),Kl(2, 2),p(5),p(6),p(1),p(2),p(3),Kl(2, 3),p(7));
    
    jrsx = @(p) ...
        jacob_right_s_x_cam0(pts3DR(:, 1),pts3DR(:, 2),Kr(1, 1),p(4),p(5),p(6),p(1),p(2),p(3),Kr(1, 3),rs1_1,rs1_2,rs1_3,rs3_1,rs3_2,rs3_3,p(7),ts1,ts3);
    
    jrsy = @(p) ...
        jacob_right_s_y_cam0(pts3DR(:, 1),pts3DR(:, 2),Kr(2, 2),p(4),p(5),p(6),p(1),p(2),p(3),Kr(2, 3),rs2_1,rs2_2,rs2_3,rs3_1,rs3_2,rs3_3,p(7),ts2,ts3);
 
    if size(pts2DL, 1) > 0 && size(pts2DR, 1) > 0
        js = @(p) [jlsx(p); jlsy(p); jrsx(p); jrsy(p)];
    elseif size(pts2DL, 1) > 0
        js = @(p) [jlsx(p); jlsy(p)];
    else
        js = @(p) [jrsx(p); jrsy(p)];
    end
end

