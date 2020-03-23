function out1 = err_right_y_cam0(X,Y,fy,ptx,pty,ptz,pwx,pwy,pwz,py,rs2_1,rs2_2,rs2_3,rs3_1,rs3_2,rs3_3,s,ts2,ts3,y)
%ERR_RIGHT_Y_CAM0
%    OUT1 = ERR_RIGHT_Y_CAM0(X,Y,FY,PTX,PTY,PTZ,PWX,PWY,PWZ,PY,RS2_1,RS2_2,RS2_3,RS3_1,RS3_2,RS3_3,S,TS2,TS3,Y)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    28-Nov-2018 10:39:37

t2 = pwx.^2;
t3 = pwy.^2;
t4 = pwz.^2;
t5 = t2+t3+t4;
t6 = sqrt(t5);
t7 = sin(t6);
t8 = 1.0./sqrt(t5);
t9 = cos(t6);
t10 = t9-1.0;
t11 = 1.0./t5;
t12 = pwz.*t7.*t8;
t13 = 1.0./s;
t14 = ptx.*rs3_1;
t15 = pty.*rs3_2;
t16 = ptz.*rs3_3;
t17 = t14+t15+t16+ts3;
t18 = pwx.*pwy.*t10.*t11;
t19 = pwy.*t7.*t8;
t20 = pwx.*pwz.*t10.*t11;
t21 = t19+t20;
t22 = t3+t4;
t23 = t10.*t11.*t22;
t24 = t23+1.0;
t25 = t12-t18;
t26 = rs3_1.*t24;
t27 = t12+t18;
t28 = pwx.*t7.*t8;
t33 = pwy.*pwz.*t10.*t11;
t29 = t28-t33;
t30 = t2+t4;
t31 = t10.*t11.*t30;
t32 = t31+1.0;
t34 = rs3_3.*t29;
t35 = rs3_2.*t32;
t36 = t34+t35-rs3_1.*t27;
t37 = rs3_2.*t25;
t38 = t26+t37-rs3_3.*t21;
out1 = y-(X.*(py.*t38+fy.*(-rs2_3.*t21+rs2_1.*t24+rs2_2.*t25))+Y.*(py.*t36+fy.*(-rs2_1.*t27+rs2_3.*t29+rs2_2.*t32))+t13.*(py.*t17+fy.*(ts2+ptx.*rs2_1+pty.*rs2_2+ptz.*rs2_3)))./(X.*t38+Y.*t36+t13.*t17);