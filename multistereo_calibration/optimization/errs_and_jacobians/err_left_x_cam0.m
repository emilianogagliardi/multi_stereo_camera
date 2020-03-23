function out1 = err_left_x_cam0(X,Y,fx,ptx,ptz,pwx,pwy,pwz,px,s,x)
%ERR_LEFT_X_CAM0
%    OUT1 = ERR_LEFT_X_CAM0(X,Y,FX,PTX,PTZ,PWX,PWY,PWZ,PX,S,X)

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
t13 = pwx.*t7.*t8;
t14 = t13-pwy.*pwz.*t10.*t11;
t15 = pwy.*t7.*t8;
t16 = pwx.*pwz.*t10.*t11;
t17 = t15+t16;
out1 = x+(X.*(px.*t17-fx.*(t10.*t11.*(t3+t4)+1.0))+Y.*(fx.*(pwz.*t7.*t8+pwx.*pwy.*t10.*t11)-px.*t14)-t12.*(fx.*ptx+ptz.*px))./(-X.*t17+Y.*t14+ptz.*t12);
