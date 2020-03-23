function [e, j_pose, j_point] = ej_wrapper_BA(multistereo, cam_idx, w, t, p3D, p2D)

[Rc, tc] = multistereo.Tprincipal2cam(cam_idx);
R1_1 = Rc(1, 1); R1_2 = Rc(1, 2); R1_3 = Rc(1, 3);
R2_1 = Rc(2, 1); R2_2 = Rc(2, 2); R2_3 = Rc(2, 3);
R3_1 = Rc(3, 1); R3_2 = Rc(3, 2); R3_3 = Rc(3, 3);
t1 = tc(1); t2 = tc(2); t3=tc(3);
K = multistereo.getKs(cam_idx);
fx = K(1, 1); fy = K(2, 2); px = K(1, 3); py = K(2, 3);
wx = w(1); wy = w(2); wz = w(3);
tx = t(1); ty = t(2); tz = t(3);
X = p3D(1); Y = p3D(2); Z = p3D(3);
x = p2D(1); y = p2D(2);

ex = ex_BA(R1_1,R1_2,R1_3,R3_1,R3_2,R3_3,X,Y,Z,fx,px,t1,t3,tx,ty,tz,wx,wy,wz,x);
ey = ey_BA(R2_1,R2_2,R2_3,R3_1,R3_2,R3_3,X,Y,Z,fy,py,t2,t3,tx,ty,tz,wx,wy,wz,y);
e = double([ex; ey]);

if nargout > 1
    jx_pose = jx_pose_BA(R1_1,R1_2,R1_3,R3_1,R3_2,R3_3,X,Y,Z,fx,px,t1,t3,tx,ty,tz,wx,wy,wz);
    jy_pose = jy_pose_BA(R2_1,R2_2,R2_3,R3_1,R3_2,R3_3,X,Y,Z,fy,py,t2,t3,tx,ty,tz,wx,wy,wz);
    jx_point = jx_point_BA(R1_1,R1_2,R1_3,R3_1,R3_2,R3_3,X,Y,Z,fx,px,t1,t3,tx,ty,tz,wx,wy,wz);
    jy_point = jy_point_BA(R2_1,R2_2,R2_3,R3_1,R3_2,R3_3,X,Y,Z,fy,py,t2,t3,tx,ty,tz,wx,wy,wz);
    j_pose = double([jx_pose; jy_pose]);
    j_point = double([jx_point; jy_point]);
end

end

