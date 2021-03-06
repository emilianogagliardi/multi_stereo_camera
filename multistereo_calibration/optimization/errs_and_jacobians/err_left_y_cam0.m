function out1 = err_left_y_cam0(X,Y,fy,pty,ptz,pwx,pwy,pwz,py,s,y)
%ERR_LEFT_Y_CAM0
%    OUT1 = ERR_LEFT_Y_CAM0(X,Y,FY,PTY,PTZ,PWX,PWY,PWZ,PY,S,Y)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    28-Nov-2018 10:39:36

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
t12 = 1.0./s;
t13 = pwy.*t7.*t8;
t14 = pwx.*pwz.*t10.*t11;
t15 = t13+t14;
t16 = pwx.*t7.*t8;
t17 = t16-pwy.*pwz.*t10.*t11;
out1 = y-(Y.*(py.*t17+fy.*(t10.*t11.*(t2+t4)+1.0))+X.*(fy.*(pwz.*t7.*t8-pwx.*pwy.*t10.*t11)-py.*t15)+t12.*(fy.*pty+ptz.*py))./(-X.*t15+Y.*t17+ptz.*t12);
