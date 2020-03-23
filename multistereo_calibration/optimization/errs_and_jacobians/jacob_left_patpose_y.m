function out1 = jacob_left_patpose_y(X,Y,cty,ctz,cwx,cwy,cwz,fy,ptx,pty,ptz,pwx,pwy,pwz,py,s)
%JACOB_LEFT_PATPOSE_Y
%    OUT1 = JACOB_LEFT_PATPOSE_Y(X,Y,CTY,CTZ,CWX,CWY,CWZ,FY,PTX,PTY,PTZ,PWX,PWY,PWZ,PY,S)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    28-Nov-2018 10:02:15

t2 = cwx.^2;
t3 = cwz.^2;
t4 = cwy.^2;
t5 = t2+t3+t4;
t6 = pwx.^2;
t7 = pwy.^2;
t8 = pwz.^2;
t9 = t6+t7+t8;
t10 = sqrt(t9);
t11 = cos(t10);
t12 = 1.0./t9;
t13 = sin(t10);
t14 = 1.0./t9.^(3.0./2.0);
t15 = t11-1.0;
t16 = sqrt(t5);
t17 = cos(t16);
t18 = t17-1.0;
t19 = 1.0./t5;
t20 = 1.0./t9.^2;
t21 = sin(t16);
t22 = 1.0./sqrt(t5);
t23 = t7+t8;
t24 = pwx.*pwy.*t13.*t14;
t25 = pwz.*t6.*t13.*t14;
t26 = pwz.*t6.*t15.*t20.*2.0;
t69 = pwz.*t12.*t15;
t70 = pwx.*pwy.*t11.*t12;
t27 = t24+t25+t26-t69-t70;
t28 = cwx.*t21.*t22;
t29 = cwy.*cwz.*t18.*t19;
t30 = pwx.*pwz.*t11.*t12;
t31 = pwy.*t6.*t13.*t14;
t32 = pwy.*t6.*t15.*t20.*2.0;
t39 = pwy.*t12.*t15;
t40 = pwx.*pwz.*t13.*t14;
t33 = t30+t31+t32-t39-t40;
t34 = pwx.*t13.*t14.*t23;
t35 = pwx.*t15.*t20.*t23.*2.0;
t36 = t34+t35;
t37 = cwz.*t21.*t22;
t99 = cwx.*cwy.*t18.*t19;
t38 = t37-t99;
t41 = t2+t3;
t42 = t18.*t19.*t41;
t43 = t42+1.0;
t44 = t6+t8;
t45 = t28+t29;
t46 = cwy.*t21.*t22;
t47 = cwx.*cwz.*t18.*t19;
t48 = t46+t47;
t49 = -t30+t31+t32-t39+t40;
t50 = t2+t4;
t51 = t18.*t19.*t50;
t52 = t51+1.0;
t53 = 1.0./sqrt(t9);
t54 = t13.*t53;
t55 = t6.*t11.*t12;
t56 = pwx.*pwy.*pwz.*t13.*t14;
t57 = pwx.*pwy.*pwz.*t15.*t20.*2.0;
t65 = t6.*t13.*t14;
t58 = t54+t55+t56+t57-t65;
t59 = t28-t29;
t60 = pwx.*t13.*t14.*t44;
t61 = pwx.*t15.*t20.*t44.*2.0;
t66 = pwx.*t12.*t15.*2.0;
t62 = t60+t61-t66;
t63 = pwz.*t13.*t53;
t64 = t48.*t49;
t67 = t59.*t62;
t68 = t64+t67-t52.*t58;
t71 = t27.*t52;
t72 = t33.*t59;
t73 = t36.*t48;
t74 = t71+t72+t73;
t75 = 1.0./s;
t76 = ptz.*t52;
t77 = pty.*t59;
t102 = ptx.*t48;
t78 = ctz+t76+t77-t102;
t79 = t75.*t78;
t80 = pwy.*t13.*t53;
t81 = pwx.*pwz.*t12.*t15;
t82 = t80+t81;
t83 = t52.*t82;
t84 = pwx.*pwy.*t12.*t15;
t85 = t12.*t15.*t23;
t86 = t85+1.0;
t87 = t48.*t86;
t88 = pwx.*t13.*t53;
t101 = pwy.*pwz.*t12.*t15;
t89 = t88-t101;
t90 = t52.*t89;
t91 = t63+t84;
t92 = t48.*t91;
t93 = t12.*t15.*t44;
t94 = t93+1.0;
t95 = t59.*t94;
t96 = t90+t92+t95;
t97 = Y.*t96;
t98 = t63-t84;
t121 = t59.*t98;
t100 = t83+t87-t121;
t103 = pwz.*t7.*t13.*t14;
t104 = pwz.*t7.*t15.*t20.*2.0;
t105 = -t24-t69+t70+t103+t104;
t106 = pwy.*pwz.*t13.*t14;
t107 = pwx.*t7.*t13.*t14;
t108 = pwx.*t7.*t15.*t20.*2.0;
t113 = pwx.*t12.*t15;
t114 = pwy.*pwz.*t11.*t12;
t109 = t106+t107+t108-t113-t114;
t110 = pwy.*t13.*t14.*t44;
t111 = pwy.*t15.*t20.*t44.*2.0;
t112 = t110+t111;
t115 = -t106+t107+t108-t113+t114;
t116 = t7.*t13.*t14;
t124 = t7.*t11.*t12;
t117 = -t54+t56+t57+t116-t124;
t118 = pwy.*t13.*t14.*t23;
t119 = pwy.*t15.*t20.*t23.*2.0;
t126 = pwy.*t12.*t15.*2.0;
t120 = t118+t119-t126;
t132 = X.*t100;
t122 = t79+t97-t132;
t123 = t59.*t115;
t125 = t52.*t117;
t127 = t48.*t120;
t128 = t123+t125+t127;
t129 = t48.*t109;
t130 = t59.*t112;
t131 = t129+t130-t52.*t105;
t133 = 1.0./t122.^2;
t134 = t43.*t98;
t135 = t45.*t82;
t136 = t38.*t86;
t137 = t134+t135+t136;
t138 = fy.*t137;
t179 = py.*t100;
t139 = t138-t179;
t140 = X.*t139;
t141 = t38.*t91;
t142 = t45.*t89;
t180 = t43.*t94;
t143 = t141+t142-t180;
t144 = fy.*t143;
t181 = py.*t96;
t145 = t144-t181;
t146 = pty.*t43;
t147 = ptx.*t38;
t183 = ptz.*t45;
t148 = cty+t146+t147-t183;
t149 = fy.*t148;
t150 = py.*t78;
t151 = t149+t150;
t152 = t75.*t151;
t182 = Y.*t145;
t153 = t140+t152-t182;
t154 = pwx.*t8.*t13.*t14;
t155 = pwx.*t8.*t15.*t20.*2.0;
t156 = t106-t113-t114+t154+t155;
t157 = t8.*t11.*t12;
t163 = t8.*t13.*t14;
t158 = t54+t56+t57+t157-t163;
t159 = pwz.*t13.*t14.*t23;
t160 = pwz.*t15.*t20.*t23.*2.0;
t162 = pwz.*t12.*t15.*2.0;
t161 = t159+t160-t162;
t164 = pwy.*t8.*t13.*t14;
t165 = pwy.*t8.*t15.*t20.*2.0;
t166 = t30-t39-t40+t164+t165;
t167 = -t54+t56+t57-t157+t163;
t168 = pwz.*t13.*t14.*t44;
t169 = pwz.*t15.*t20.*t44.*2.0;
t170 = -t162+t168+t169;
t171 = 1.0./t122;
t172 = t52.*t156;
t173 = t59.*t158;
t174 = t48.*t161;
t175 = t172+t173+t174;
t176 = t48.*t167;
t177 = t59.*t170;
t178 = t176+t177-t52.*t166;
out1 = [-(X.*(py.*t74-fy.*(t27.*t45+t36.*t38-t33.*t43))-Y.*(py.*t68+fy.*(-t38.*t49+t45.*t58+t43.*t62)))./(t79+t97-X.*(t83+t87-t59.*(t63-pwx.*pwy.*t12.*t15)))+t133.*t153.*(X.*t74-Y.*t68),-t171.*(X.*(py.*t128-fy.*(t38.*t120-t43.*t115+t45.*t117))-Y.*(py.*t131+fy.*(-t38.*t109+t45.*t105+t43.*t112)))+t133.*t153.*(X.*t128-Y.*t131),-t171.*(X.*(py.*t175-fy.*(t38.*t161-t43.*t158+t45.*t156))-Y.*(py.*t178+fy.*(-t38.*t167+t45.*t166+t43.*t170)))+t133.*t153.*(X.*t175-Y.*t178),-t75.*t171.*(fy.*t38-py.*t48)-t48.*t75.*t133.*t153,-t75.*t171.*(fy.*t43+py.*t59)+t59.*t75.*t133.*t153,t75.*t171.*(fy.*t45-py.*t52)+t52.*t75.*t133.*t153];
