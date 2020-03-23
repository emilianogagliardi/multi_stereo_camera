clear 
close all
clc

wx = sym('wx'); % optimized rotation variable, SO(3)
wy = sym('wy');
wz = sym('wz');
tx = sym('tx'); % optimized translation variable R^3
ty = sym('ty');
tz = sym('tz');
fx = sym('fx'); % intrinsics
fy = sym('fy');
px = sym('px');
py = sym('py');
X = sym('X');
Y = sym('Y');
Z = sym('Z');
x = sym('x');
y = sym('y');
Rc = sym('R', 3);
tc = sym('t', [3, 1]);

K = [fx 0  px;
     0  fy py;
     0  0  1];

theta = sqrt(wx^2 + wy^2 + wz^2);
omega=[0   -wz wy;
       wz  0   -wx;
       -wy wx  0];
R = eye(3)+(sin(theta)/theta)*omega+((1-cos(theta))/theta^2)*(omega*omega);
t = [tx; ty; tz];

proj = K * [Rc * R, Rc * t + tc] * [X; Y; Z; 1];
proj = proj(1:2) ./ proj(3);

ex = x - proj(1);
ey = y - proj(2);

jx_pose = jacobian(ex, [wx, wy, wz, tx, ty, tz]);
jy_pose = jacobian(ey, [wx, wy, wz, tx, ty, tz]);
jx_point = jacobian(ex, [X, Y, Z]);
jy_point = jacobian(ey, [X, Y, Z]);

matlabFunction(ex, 'File', '~/Development/matlab/multi_stereo_camera/VO/BAmap/errs_and_jacobs/ex_BA.m');
matlabFunction(ey, 'File', '~/Development/matlab/multi_stereo_camera/VO/BAmap/errs_and_jacobs/ey_BA.m');
matlabFunction(jx_pose, 'File', '~/Development/matlab/multi_stereo_camera/VO/BAmap/errs_and_jacobs/jx_pose_BA.m');
matlabFunction(jy_pose, 'File', '~/Development/matlab/multi_stereo_camera/VO/BAmap/errs_and_jacobs/jy_pose_BA.m');
matlabFunction(jx_point, 'File', '~/Development/matlab/multi_stereo_camera/VO/BAmap/errs_and_jacobs/jx_point_BA.m');
matlabFunction(jy_point, 'File', '~/Development/matlab/multi_stereo_camera/VO/BAmap/errs_and_jacobs/jy_point_BA.m');