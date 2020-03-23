function [e,j] = ej_wrapper_reprojection(K, R, t, p3D, p2D)
% return error function and jacobian function for reprojection error minimization in a calibrated multicamera, for one of its camera
% The returned function can be passed to lsqnonlin (or all the functions
% together, one for each camera with it corresponding point)
% Parameters:
% - R and t bring point from principal camera frame to the considered camera frame
% - p3D, p2D: 3D-2D pairs perceived by the considered camera

% err_fun_x_ @(Ri1_1,Ri1_2,Ri1_3,Ri3_1,Ri3_2,Ri3_3,X,Y,Z,fx,px,ti1,ti3,tx,ty,tz,wx,wy,wz,x)
% err_fun_y_  @(Ri2_1,Ri2_2,Ri2_3,Ri3_1,Ri3_2,Ri3_3,X,Y,Z,fy,py,ti2,ti3,tx,ty,tz,wx,wy,wz,y)
% jacob_fun_x_ @(Ri1_1,Ri1_2,Ri1_3,Ri3_1,Ri3_2,Ri3_3,X,Y,Z,fx,px,ti1,ti3,tx,ty,tz,wx,wy,wz)
% jacob_fun_y_ @(Ri2_1,Ri2_2,Ri2_3,Ri3_1,Ri3_2,Ri3_3,X,Y,Z,fy,py,ti2,ti3,tx,ty,tz,wx,wy,wz)

ti1 = t(1); ti2 = t(2); ti3 = t(3);
Ri1_1 = R(1, 1); Ri1_2 = R(1, 2); Ri1_3 = R(1, 3);
Ri2_1 = R(2, 1); Ri2_2 = R(2, 2); Ri2_3 = R(2, 3);
Ri3_1 = R(3, 1); Ri3_2 = R(3, 2); Ri3_3 = R(3, 3);

fx = K(1, 1);
px = K(1, 3);
fy = K(2, 2);
py = K(2, 3);

ex = @(p) ...
    err_fun_x_(Ri1_1,Ri1_2,Ri1_3,Ri3_1,Ri3_2,Ri3_3,p3D(:, 1),p3D(:, 2),p3D(:, 3),fx,px,ti1,ti3,p(4),p(5),p(6),p(1),p(2),p(3),p2D(:, 1));
ey = @(p) ...
    err_fun_y_(Ri2_1,Ri2_2,Ri2_3,Ri3_1,Ri3_2,Ri3_3,p3D(:, 1),p3D(:, 2),p3D(:, 3),fy,py,ti2,ti3,p(4),p(5),p(6),p(1),p(2),p(3),p2D(:, 2));

e = @(p) [ex(p); ey(p)];

jx = @(X,Y,Z,tx,ty,tz,wx,wy,wz) ...
    jacob_fun_x_(Ri1_1,Ri1_2,Ri1_3,Ri3_1,Ri3_2,Ri3_3,X,Y,Z,fx,px,ti1,ti3,tx,ty,tz,wx,wy,wz);
jx = @(p) cell2mat(arrayfun(jx, p3D(:, 1), p3D(:, 2), p3D(:, 3), ...
        repmat(p(4), [size(p3D, 1), 1]), repmat(p(5), [size(p3D, 1), 1]), repmat(p(6), [size(p3D, 1), 1]), repmat(p(1), [size(p3D, 1), 1]), repmat(p(2), [size(p3D, 1), 1]), repmat(p(3), [size(p3D, 1), 1]), 'UniformOutput', false));

jy = @(X,Y,Z,tx,ty,tz,wx,wy,wz) ...
    jacob_fun_y_(Ri2_1,Ri2_2,Ri2_3,Ri3_1,Ri3_2,Ri3_3,X,Y,Z,fy,py,ti2,ti3,tx,ty,tz,wx,wy,wz);
jy = @(p) cell2mat(arrayfun(jy, p3D(:, 1), p3D(:, 2), p3D(:, 3), ...
        repmat(p(4), [size(p3D, 1), 1]), repmat(p(5), [size(p3D, 1), 1]), repmat(p(6), [size(p3D, 1), 1]), repmat(p(1), [size(p3D, 1), 1]), repmat(p(2), [size(p3D, 1), 1]), repmat(p(3), [size(p3D, 1), 1]), 'UniformOutput', false));
    
j = @(p) [jx(p); jy(p)];
end

