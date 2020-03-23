function out1 = jacob_left_s_y_cam0(X,Y,fy,pty,ptz,pwx,pwy,pwz,py,s)
%JACOB_LEFT_S_Y_CAM0
%    OUT1 = JACOB_LEFT_S_Y_CAM0(X,Y,FY,PTY,PTZ,PWX,PWY,PWZ,PY,S)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    28-Nov-2018 10:39:58

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
t12 = 1.0./s.^2;
t13 = pwy.*t7.*t8;
t14 = pwx.*pwz.*t10.*t11;
t15 = t13+t14;
t16 = pwx.*t7.*t8;
t25 = pwy.*pwz.*t10.*t11;
t17 = t16-t25;
t18 = Y.*t17;
t19 = 1.0./s;
t20 = ptz.*t19;
t21 = t18+t20-X.*t15;
t22 = fy.*pty;
t23 = ptz.*py;
t24 = t22+t23;
out1 = (t12.*t24)./t21-ptz.*t12.*1.0./t21.^2.*(Y.*(py.*t17+fy.*(t10.*t11.*(t2+t4)+1.0))+t19.*t24+X.*(fy.*(pwz.*t7.*t8-pwx.*pwy.*t10.*t11)-py.*t15));
