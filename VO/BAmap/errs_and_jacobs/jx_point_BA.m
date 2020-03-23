function jx_point = jx_point_BA(R1_1,R1_2,R1_3,R3_1,R3_2,R3_3,X,Y,Z,fx,px,t1,t3,tx,ty,tz,wx,wy,wz)
%JX_POINT_BA
%    JX_POINT = JX_POINT_BA(R1_1,R1_2,R1_3,R3_1,R3_2,R3_3,X,Y,Z,FX,PX,T1,T3,TX,TY,TZ,WX,WY,WZ)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    27-Feb-2019 12:26:07

t4 = wx.^2;
t5 = wy.^2;
t6 = wz.^2;
t7 = t4+t5+t6;
t8 = sqrt(t7);
t9 = sin(t8);
t10 = 1.0./sqrt(t7);
t11 = cos(t8);
t12 = t11-1.0;
t13 = 1.0./t7;
t14 = t9.*t10.*wz;
t22 = t12.*t13.*wx.*wy;
t15 = t14-t22;
t16 = t9.*t10.*wy;
t17 = t12.*t13.*wx.*wz;
t18 = t16+t17;
t19 = t5+t6;
t20 = t12.*t13.*t19;
t21 = t20+1.0;
t23 = R3_2.*t15;
t24 = R3_1.*t21;
t27 = R3_3.*t18;
t25 = t23+t24-t27;
t26 = t9.*t10.*wx;
t28 = R3_1.*tx;
t29 = R3_2.*ty;
t30 = R3_3.*tz;
t31 = R1_2.*t15;
t32 = R1_1.*t21;
t73 = R1_3.*t18;
t33 = t31+t32-t73;
t34 = fx.*t33;
t35 = px.*t25;
t36 = t34+t35;
t37 = t14+t22;
t38 = t12.*t13.*wy.*wz;
t39 = t4+t6;
t40 = t12.*t13.*t39;
t41 = t40+1.0;
t42 = t26-t38;
t43 = R3_2.*t41;
t44 = t16-t17;
t45 = t26+t38;
t46 = t4+t5;
t47 = t12.*t13.*t46;
t48 = t47+1.0;
t49 = R3_1.*t44;
t50 = R3_3.*t48;
t55 = R3_2.*t45;
t51 = t49+t50-t55;
t52 = X.*t25;
t53 = R3_3.*t42;
t61 = R3_1.*t37;
t54 = t43+t53-t61;
t56 = Z.*t51;
t57 = R1_3.*t42;
t58 = R1_2.*t41;
t75 = R1_1.*t37;
t59 = t57+t58-t75;
t60 = fx.*t59;
t62 = px.*t54;
t63 = t60+t62;
t64 = Y.*t54;
t65 = t3+t28+t29+t30+t52+t56+t64;
t66 = R1_1.*tx;
t67 = R1_2.*ty;
t68 = R1_3.*tz;
t69 = t1+t66+t67+t68;
t70 = fx.*t69;
t71 = t3+t28+t29+t30;
t72 = px.*t71;
t74 = X.*t36;
t76 = Y.*t63;
t77 = R1_1.*t44;
t78 = R1_3.*t48;
t86 = R1_2.*t45;
t79 = t77+t78-t86;
t80 = fx.*t79;
t81 = px.*t51;
t82 = t80+t81;
t83 = Z.*t82;
t84 = t70+t72+t74+t76+t83;
t85 = 1.0./t65.^2;
t87 = 1.0./t65;
jx_point = [-t36./(t3+t28+t29+t30+t52+t56+Y.*(t43-R3_1.*t37+R3_3.*(t26-t12.*t13.*wy.*wz)))+t25.*t84.*t85,-t63.*t87+t54.*t84.*t85,-t82.*t87+t51.*t84.*t85];
