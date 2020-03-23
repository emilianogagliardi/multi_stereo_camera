% this scripts generates matlab functions that computes the reprojection
% error and the jacobian as a function of:
% - the position of the multicamera, defined as the position of the
% principal camera, expressed in SO(3)xR^3
% - the intrinsic parameter of the camera
% - the position of the considered camera with respect to the principal
% camera, known from multicamera calibration
% - the 3D point
% - the 2D point

% All the quantities, except the position of the multicamera will be known.
% Thus you can partially apply this function with the known quanitites to
% obtain an error function with respect to only the multicamera pose

clear
close all
clc

% optimized variable
wx = sym('wx'); % optimized rotation variable, SO(3)
wy = sym('wy');
wz = sym('wz');
tx = sym('tx'); % optimized translation variable R^3
ty = sym('ty');
tz = sym('tz');

% intrinsic
% pose
Ri = sym('Ri', [3, 3]);
ti = sym('ti', [3, 1]);
% calibration matrix
fx = sym('fx');
fy = sym('fy');
px = sym('px');
py = sym('py');

% point
X = sym('X');
Y = sym('Y');
Z = sym('Z');
x = sym('x');
y = sym('y');


 K = [fx, 0, px;
     0, fy, py;
     0, 0, 1];

% rodrigues
theta = sqrt(wx^2+wy^2+wz^2);
omega = [0 -wz wy;
         wz 0 -wx;
        -wy wx 0;];
R = eye(3) + (sin(theta)/theta)*omega + ((1-cos(theta))/theta^2)*(omega*omega);
t = [tx;ty;tz];

% projection
proj = K*[Ri * R,  Ri * t + ti] * [X; Y; Z; 1];
proj = proj(1:2, :) ./proj(3, :);
diff = [x; y] - proj;

err_sym_x = diff(1); 
err_sym_y = diff(2);
jacob_sym_x = jacobian(err_sym_x,[wx wy wz tx ty tz]);
jacob_sym_y = jacobian(err_sym_y,[wx wy wz tx ty tz]);

 % err_fun_x_ @(Ri1_1,Ri1_2,Ri1_3,Ri3_1,Ri3_2,Ri3_3,X,Y,Z,fx,px,ti1,ti3,tx,ty,tz,wx,wy,wz,x)
err_fun_x_ = matlabFunction(err_sym_x, 'File', 'multistereo_localization/multistereo/localization/3D2D/coupled/reproj_err_minimization/errs_and_jacobs/err_fun_x_');

% err_fun_y_  @(Ri2_1,Ri2_2,Ri2_3,Ri3_1,Ri3_2,Ri3_3,X,Y,Z,fy,py,ti2,ti3,tx,ty,tz,wx,wy,wz,y)
err_fun_y_ = matlabFunction(err_sym_y, 'File', 'multistereo_localization/multistereo/localization/3D2D/coupled/reproj_err_minimization/errs_and_jacobs/err_fun_y_');

% jacob_fun_x_ @(Ri1_1,Ri1_2,Ri1_3,Ri3_1,Ri3_2,Ri3_3,X,Y,Z,fx,px,ti1,ti3,tx,ty,tz,wx,wy,wz)
jacob_fun_x_ = matlabFunction(jacob_sym_x, 'File', 'multistereo_localization/multistereo/localization/3D2D/coupled/reproj_err_minimization/errs_and_jacobs/jacob_fun_x_'); 

 % jacob_fun_y_ @(Ri2_1,Ri2_2,Ri2_3,Ri3_1,Ri3_2,Ri3_3,X,Y,Z,fy,py,ti2,ti3,tx,ty,tz,wx,wy,wz)
jacob_fun_y_ = matlabFunction(jacob_sym_y, 'File', 'multistereo_localization/multistereo/localization/3D2D/coupled/reproj_err_minimization/errs_and_jacobs/jacob_fun_y_');



