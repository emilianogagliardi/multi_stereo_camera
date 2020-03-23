function out1 = jacob_left_s_x(X,Y,ctx,ctz,cwx,cwy,cwz,fx,ptx,pty,ptz,pwx,pwy,pwz,px,s)
%JACOB_LEFT_S_X
%    OUT1 = JACOB_LEFT_S_X(X,Y,CTX,CTZ,CWX,CWY,CWZ,FX,PTX,PTY,PTZ,PWX,PWY,PWZ,PX,S)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    28-Nov-2018 10:05:14

t2 = cwy.^2;
t3 = cwz.^2;
t4 = cwx.^2;
t5 = t2+t3+t4;
t6 = sqrt(t5);
t7 = cos(t6);
t8 = t7-1.0;
t9 = 1.0./t5;
t10 = sin(t6);
t11 = 1.0./sqrt(t5);
t12 = cwy.*t10.*t11;
t13 = t2+t4;
t14 = t8.*t9.*t13;
t15 = t14+1.0;
t16 = ptz.*t15;
t17 = cwx.*cwz.*t8.*t9;
t18 = t12+t17;
t19 = cwx.*t10.*t11;
t33 = cwy.*cwz.*t8.*t9;
t20 = t19-t33;
t21 = pty.*t20;
t37 = ptx.*t18;
t22 = ctz+t16+t21-t37;
t23 = pwx.^2;
t24 = pwy.^2;
t25 = pwz.^2;
t26 = t23+t24+t25;
t27 = sqrt(t26);
t28 = sin(t27);
t29 = 1.0./sqrt(t26);
t30 = cos(t27);
t31 = t30-1.0;
t32 = 1.0./t26;
t34 = pwz.*t28.*t29;
t35 = 1.0./s.^2;
t36 = 1.0./s;
t38 = t22.*t36;
t39 = pwy.*t28.*t29;
t40 = pwx.*pwz.*t31.*t32;
t41 = t39+t40;
t42 = t15.*t41;
t43 = pwx.*pwy.*t31.*t32;
t44 = t24+t25;
t45 = t31.*t32.*t44;
t46 = t45+1.0;
t47 = t18.*t46;
t48 = pwx.*t28.*t29;
t67 = pwy.*pwz.*t31.*t32;
t49 = t48-t67;
t50 = t15.*t49;
t51 = t34+t43;
t52 = t18.*t51;
t53 = t23+t25;
t54 = t31.*t32.*t53;
t55 = t54+1.0;
t56 = t20.*t55;
t57 = t50+t52+t56;
t58 = Y.*t57;
t59 = t34-t43;
t60 = cwz.*t10.*t11;
t61 = cwx.*cwy.*t8.*t9;
t62 = t60+t61;
t63 = t2+t3;
t64 = t8.*t9.*t63;
t65 = t64+1.0;
t66 = t42+t47-t20.*t59;
t68 = t12-t17;
t69 = ptx.*t65;
t70 = px.*t22;
out1 = (t35.*(t70+fx.*(ctx+t69-pty.*t62+ptz.*(t12-cwx.*cwz.*t8.*t9))))./(t38+t58-X.*(t42+t47-t20.*(t34-pwx.*pwy.*t31.*t32)))-t22.*t35.*(t36.*(t70+fx.*(ctx+t69-pty.*t62+ptz.*t68))-X.*(px.*t66+fx.*(t41.*t68-t46.*t65+t59.*t62))+Y.*(px.*t57-fx.*(t51.*t65-t49.*t68+t55.*t62))).*1.0./(t38+t58-X.*t66).^2;